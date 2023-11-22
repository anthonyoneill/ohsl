use core::ops::{Index, IndexMut};
pub use crate::vector::{Vector, Vec64};

// Sparse matrices and sparse linear algebra - https://crates.io/crates/sparse21 with modifications

use std::usize::MAX;
use std::cmp::{min, max};

use more_asserts::assert_gt;

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct Eindex(pub usize);

// Triplets are a type alias for tuples of (row, col, val).
type Triplet = (usize, usize, f64);

#[derive(Debug, Copy, Clone)]
pub enum Axis { ROWS = 0, COLS }

use Axis::*;

impl Axis {
    fn other(&self) -> Axis {
        match self {
            Axis::ROWS => Axis::COLS,
            Axis::COLS => Axis::ROWS,
        }
    }
}

struct AxisPair<T> {
    rows: T,
    cols: T,
}

impl<T> Index<Axis> for AxisPair<T> {
    type Output = T;

    fn index(&self, ax: Axis) -> &Self::Output {
        match ax {
            Axis::ROWS => &self.rows,
            Axis::COLS => &self.cols,
        }
    }
}

impl<T> IndexMut<Axis> for AxisPair<T> {
    fn index_mut(&mut self, ax: Axis) -> &mut Self::Output {
        match ax {
            Axis::ROWS => &mut self.rows,
            Axis::COLS => &mut self.cols,
        }
    }
}

#[derive(PartialEq, Debug, Copy, Clone)]
pub enum MatrixState { CREATED = 0, FACTORING, FACTORED }

#[derive(Debug, Clone)]
pub struct Element {
    index: Eindex,
    row: usize,
    col: usize,
    val: f64,
    fillin: bool,
    //orig: (usize, usize, f64),
    next_in_row: Option<Eindex>,
    next_in_col: Option<Eindex>,
}

impl PartialEq for Element {
    fn eq(&self, other: &Self) -> bool {
        return self.row == other.row &&
            self.col == other.col &&
            self.val == other.val;
    }
}

impl Element {
    /// Create a new element
    pub fn new(index: Eindex, row: usize, col: usize, val: f64, fillin: bool) -> Element {
        Element {
            index,
            row,
            col,
            val,
            fillin,
            //orig: (row, col, val),
            next_in_row: None,
            next_in_col: None,
        }
    }
    /// Return the element index 
    pub fn index(&self) -> Eindex {
        self.index
    }
    /// Return the element row 
    pub fn row(&self) -> usize {
        self.row
    }
    /// Return the element column 
    pub fn col(&self) -> usize {
        self.col
    }
    /// Return the value stored in the element
    pub fn val(&self) -> f64 {
        self.val
    }
    fn loc(&self, ax: Axis) -> usize {
        match ax {
            Axis::ROWS => self.row,
            Axis::COLS => self.col,
        }
    }
    fn set_loc(&mut self, ax: Axis, to: usize) {
        match ax {
            Axis::ROWS => self.row = to,
            Axis::COLS => self.col = to,
        }
    }
    fn next(&self, ax: Axis) -> Option<Eindex> {
        match ax {
            Axis::ROWS => self.next_in_row,
            Axis::COLS => self.next_in_col,
        }
    }
    fn set_next(&mut self, ax: Axis, e: Option<Eindex>) {
        match ax {
            Axis::ROWS => self.next_in_row = e,
            Axis::COLS => self.next_in_col = e,
        }
    }
}

struct AxisMapping {
    e2i: Vec<usize>,
    i2e: Vec<usize>,
    history: Vec<(usize, usize)>,
}

impl AxisMapping {
    fn new(size: usize) -> AxisMapping {
        AxisMapping {
            e2i: (0..size).collect(),
            i2e: (0..size).collect(),
            history: vec![],
        }
    }
    fn swap_int(&mut self, x: usize, y: usize) {
        // Swap internal indices x and y
        let tmp = self.i2e[x];
        self.i2e[x] = self.i2e[y];
        self.i2e[y] = tmp;
        self.e2i[self.i2e[x]] = x;
        self.e2i[self.i2e[y]] = y;
        self.history.push((x, y));
    }
}

pub struct AxisData {
    hdrs: Vec<Option<Eindex>>,
    qtys: Vec<usize>,
    markowitz: Vec<usize>,
    mapping: Option<AxisMapping>,
}

impl AxisData {
    fn new() -> AxisData {
        AxisData {
            hdrs: vec![],
            qtys: vec![],
            markowitz: vec![],
            mapping: None,
        }
    }
    fn grow(&mut self, to: usize) {
        if to <= self.hdrs.len() { return; }
        let by = to - self.hdrs.len();
        for _ in 0..by {
            self.hdrs.push(None);
            self.qtys.push(0);
            self.markowitz.push(0);
        }
    }
    fn setup_factoring(&mut self) {
        self.markowitz.copy_from_slice(&self.qtys);
        self.mapping = Some(AxisMapping::new(self.hdrs.len()));
    }
    fn swap(&mut self, x: usize, y: usize) {
        self.hdrs.swap(x, y);
        self.qtys.swap(x, y);
        self.markowitz.swap(x, y);
        if let Some(m) = &mut self.mapping {
            m.swap_int(x, y);
        }
    }
}

type SpResult<T> = Result<T, &'static str>;

/// Sparse matrix
pub struct Sparse {
    // Sparse.elements is the owner of all `Element`s.
    // Everything else gets referenced via `Eindex`es.
    state: MatrixState,
    elements: Vec<Element>,
    axes: AxisPair<AxisData>,
    diag: Vec<Option<Eindex>>,
    fillins: Vec<Eindex>,
}

impl Sparse {
    /// Create a new, initially empty Sparse matrix
    pub fn new() -> Sparse {
        Sparse {
            state: MatrixState::CREATED,
            axes: AxisPair {                //TODO do we need axes why not just rows and cols ? 
                rows: AxisData::new(),
                cols: AxisData::new(),
            },
            diag: vec![],
            elements: vec![],
            fillins: vec![],
        }
    }
    /// Return the state of the Sparse matrix
    pub fn state(&self) -> MatrixState {
        self.state
    }
    /// Return the diagonal elements 
    pub fn diag(&self) -> Vec<Option<Eindex>> {
        self.diag.clone()
    }
    /// Return the vector of elements 
    pub fn elements(&self) -> Vec<Element> {
        self.elements.clone()
    }
    /// Create a new Sparse matrix from a vector of (row, col, val) `entries`.
    pub fn from_triplets(entries: Vec<Triplet>) -> Sparse {
        let mut m = Sparse::new();
        for e in entries.iter() {
            m.add_element(e.0, e.1, e.2);
        }
        return m;
    }
    /// Create an n*n identity Sparse matrix 
    pub fn identity(n: usize) -> Sparse {
        let mut m = Sparse::new();
        for k in 0..n { m.add_element(k, k, 1.0); }
        return m;
    }
    /// Add an element at location `(row, col)` with value `val`. 
    pub fn add_element(&mut self, row: usize, col: usize, val: f64) {
        self._add_element(row, col, val, false);
    }
    /// Add elements correspoding to each triplet `(row, col, val)`
    /// Rows and columns are `usize`, and `vals` are `f64`.
    pub fn add_elements(&mut self, elements: Vec<Triplet>) {
        for e in elements.iter() {
            self.add_element(e.0, e.1, e.2);
        }
    }

    fn insert(&mut self, e: &mut Element) {
        let mut expanded = false;
        if e.row + 1 > self.num_rows() {
            self.axes[Axis::ROWS].grow(e.row + 1);
            expanded = true;
        }
        if e.col + 1 > self.num_cols() {
            self.axes[Axis::COLS].grow(e.col + 1);
            expanded = true;
        }
        if expanded {
            let new_diag_len = std::cmp::min(self.num_rows(), self.num_cols());
            for _ in 0..new_diag_len - self.diag.len() {
                self.diag.push(None);
            }
        }
        // Insert along each Axis
        self.insert_axis(Axis::COLS, e);
        self.insert_axis(Axis::ROWS, e);
        // Update row & col qtys
        self.axes[ROWS].qtys[e.row] += 1;
        self.axes[COLS].qtys[e.col] += 1;
        if self.state == MatrixState::FACTORING {
            self.axes[ROWS].markowitz[e.row] += 1;
            self.axes[COLS].markowitz[e.col] += 1;
        }
        // Update our special arrays
        if e.row == e.col { self.diag[e.row] = Some(e.index); }
        if e.fillin { self.fillins.push(e.index); }
    }
    fn insert_axis(&mut self, ax: Axis, e: &mut Element) {
        // Insert Element `e` along Axis `ax`

        let head_ptr = self.axes[ax].hdrs[e.loc(ax)];
        let head_idx = match head_ptr {
            Some(h) => h,
            None => {
                // Adding first element in this row/col
                return self.set_hdr(ax, e.loc(ax), Some(e.index));
            }
        };
        let off_ax = ax.other();
        if self[head_idx].loc(off_ax) > e.loc(off_ax) {
            // `e` is the new first element
            e.set_next(ax, head_ptr);
            return self.set_hdr(ax, e.loc(ax), Some(e.index));
        }

        // `e` comes after at least one Element.  Search for its position.
        let mut prev = head_idx;
        while let Some(next) = self[prev].next(ax) {
            if self[next].loc(off_ax) >= e.loc(off_ax) { break; }
            prev = next;
        }
        // And splice it in-between `prev` and `nxt`
        e.set_next(ax, self[prev].next(ax));
        self[prev].set_next(ax, Some(e.index));
    }
    fn add_fillin(&mut self, row: usize, col: usize) -> Eindex {
        return self._add_element(row, col, 0.0, true);
    }
    fn _add_element(&mut self, row: usize, col: usize, val: f64, fillin: bool) -> Eindex {
        // Element creation & insertion, used by `add_fillin` and the public `add_element`.
        let index = Eindex(self.elements.len());
        let mut e = Element::new(index.clone(), row, col, val, fillin);
        self.insert(&mut e);
        self.elements.push(e);
        return index;
    }
    /// Returns the Element-value at `(row, col)` if present, or None if not.
    pub fn get(&self, row: usize, col: usize) -> Option<f64> {
        if row >= self.num_rows() { return None; }
        if col >= self.num_cols() { return None; }

        if row == col { // On diagonal; easy access
            return match self.diag[row] {
                None => None,
                Some(d) => Some(self[d].val),
            };
        }
        // Off-diagonal. Search across `row`. 
        let mut ep = self.hdr(ROWS, row);
        while let Some(ei) = ep {
            let e = &self[ei];
            if e.col == col { return Some(e.val); } 
            else if e.col > col { return None; }
            ep = e.next_in_row;
        }
        return None;
    }
    /// Make major state transitions
    fn set_state(&mut self, state: MatrixState) -> Result<(), &'static str> {
        match state {
            MatrixState::CREATED => Err("Matrix State Error"),
            MatrixState::FACTORING => {
                if self.state == MatrixState::FACTORING { return Ok(()); }
                if self.state == MatrixState::FACTORED { return Err("Already Factored"); }

                self.axes[Axis::ROWS].setup_factoring();
                self.axes[Axis::COLS].setup_factoring();

                self.state = state;
                return Ok(());
            }
            MatrixState::FACTORED => {
                if self.state == MatrixState::FACTORING {
                    self.state = state;
                    return Ok(());
                } else { return Err("Matrix State Error"); }
            }
        }
    }
    fn move_element(&mut self, ax: Axis, idx: Eindex, to: usize) {
        let loc = self[idx].loc(ax);
        if loc == to { return; }
        let off_ax = ax.other();
        let y = self[idx].loc(off_ax);

        if loc < to {
            let br = match self.before_loc(off_ax, y, to, Some(idx)) {
                Some(ei) => ei,
                None => panic!("ERROR"),
            };
            if br != idx {
                let be = self.prev(off_ax, idx, None);
                let nxt = self[idx].next(off_ax);
                match be {
                    None => self.set_hdr(off_ax, y, nxt),
                    Some(be) => self[be].set_next(off_ax, nxt),
                };
                let brn = self[br].next(off_ax);
                self[idx].set_next(off_ax, brn);
                self[br].set_next(off_ax, Some(idx));
            }
        } else {
            let br = self.before_loc(off_ax, y, to, None);
            let be = self.prev(off_ax, idx, None);

            if br != be { // We (may) need some pointer updates
                if let Some(ei) = be {
                    let nxt = self[idx].next(off_ax);
                    self[ei].set_next(off_ax, nxt);
                }
                match br {
                    None => { // New first in row/col
                        let first = self.hdr(off_ax, y);
                        self[idx].set_next(off_ax, first);
                        self.axes[off_ax].hdrs[y] = Some(idx);
                    }
                    Some(br) => {
                        if br != idx { // Splice `idx` in after `br`
                            let nxt = self[br].next(off_ax);
                            self[idx].set_next(off_ax, nxt);
                            self[br].set_next(off_ax, Some(idx));
                        }
                    }
                };
            }
        }

        // Update the moved-Element's location
        self[idx].set_loc(ax, to);

        if loc == y { // If idx was on our diagonal, remove it
            self.diag[loc] = None;
        } else if to == y { // Or if it's now on the diagonal, add it
            self.diag[to] = Some(idx);
        }
    }
    fn exchange_elements(&mut self, ax: Axis, ix: Eindex, iy: Eindex) {
        // Swap two elements `ax` indices.
        // Elements must be in the same off-axis vector,
        // and the first argument `ex` must be the lower-indexed off-axis.
        // E.g. exchange_elements(Axis.rows, ex, ey) exchanges the rows of ex and ey.

        let off_ax = ax.other();
        let off_loc = self[ix].loc(off_ax);

        let bx = self.prev(off_ax, ix, None);
        let by = match self.prev(off_ax, iy, Some(ix)) {
            Some(e) => e,
            None => panic!("ERROR!"),
        };

        let locx = self[ix].loc(ax);
        let locy = self[iy].loc(ax);
        self[iy].set_loc(ax, locx);
        self[ix].set_loc(ax, locy);

        match bx {
            None => {
                // If `ex` is the *first* entry in the column, replace it to our header-list
                self.set_hdr(off_ax, off_loc, Some(iy));
            }
            Some(bxe) => {
                // Otherwise patch ey into bx
                self[bxe].set_next(off_ax, Some(iy));
            }
        }

        if by == ix { // `ex` and `ey` are adjacent
            let tmp = self[iy].next(off_ax);
            self[iy].set_next(off_ax, Some(ix));
            self[ix].set_next(off_ax, tmp);
        } else { // Elements in-between `ex` and `ey`.  Update the last one.
            let xnxt = self[ix].next(off_ax);
            let ynxt = self[iy].next(off_ax);
            self[iy].set_next(off_ax, xnxt);
            self[ix].set_next(off_ax, ynxt);
            self[by].set_next(off_ax, Some(ix));
        }

        // Update our diagonal array, if necessary
        if locx == off_loc {
            self.diag[off_loc] = Some(iy);
        } else if locy == off_loc {
            self.diag[off_loc] = Some(ix);
        }
    }
    fn prev(&self, ax: Axis, idx: Eindex, hint: Option<Eindex>) -> Option<Eindex> {
        // Find the element previous to `idx` along axis `ax`. 
        // If provided, `hint` *must* be before `idx`, or search will fail. 
        let prev: Option<Eindex> = match hint {
            Some(_) => hint,
            None => self.hdr(ax, self[idx].loc(ax)),
        };
        let mut pi: Eindex = match prev {
            None => { return None; }
            Some(pi) if pi == idx => { return None; }
            Some(pi) => pi,
        };
        while let Some(nxt) = self[pi].next(ax) {
            if nxt == idx { break; }
            pi = nxt;
        }
        return Some(pi);
    }
    fn before_loc(&self, ax: Axis, loc: usize, before: usize, hint: Option<Eindex>) -> Option<Eindex> {
        let prev: Option<Eindex> = match hint {
            Some(_) => hint,
            None => self.hdr(ax, loc),
        };
        let off_ax = ax.other();
        let mut pi: Eindex = match prev {
            None => { return None; }
            Some(pi) if self[pi].loc(off_ax) >= before => { return None; }
            Some(pi) => pi,
        };
        while let Some(nxt) = self[pi].next(ax) {
            if self[nxt].loc(off_ax) >= before { break; }
            pi = nxt;
        }
        return Some(pi);
    }
    fn swap(&mut self, ax: Axis, a: usize, b: usize) {
        if a == b { return; }
        let x = min(a, b);
        let y = max(a, b);

        let hdrs = &self.axes[ax].hdrs;
        let mut ix = hdrs[x];
        let mut iy = hdrs[y];
        let off_ax = ax.other();

        loop {
            match (ix, iy) {
                (Some(ex), Some(ey)) => {
                    let ox = self[ex].loc(off_ax);
                    let oy = self[ey].loc(off_ax);
                    if ox < oy {
                        self.move_element(ax, ex, y);
                        ix = self[ex].next(ax);
                    } else if oy < ox {
                        self.move_element(ax, ey, x);
                        iy = self[ey].next(ax);
                    } else {
                        self.exchange_elements(ax, ex, ey);
                        ix = self[ex].next(ax);
                        iy = self[ey].next(ax);
                    }
                }
                (None, Some(ey)) => {
                    self.move_element(ax, ey, x);
                    iy = self[ey].next(ax);
                }
                (Some(ex), None) => {
                    self.move_element(ax, ex, y);
                    ix = self[ex].next(ax);
                }
                (None, None) => { break; }
            }
        }
        // Swap all the relevant pointers & counters
        self.axes[ax].swap(x, y);
    }
    /// Updates self to S = L + U - I.
    /// Diagonal entries are those of U;
    /// L has diagonal entries equal to one.
    fn lu_factorize(&mut self) -> SpResult<()> { 
        assert_gt!( self.diag.len(), 0 );
        for k in 0..self.axes[ROWS].hdrs.len() {
            if self.hdr(ROWS, k).is_none() { return Err("Singular Matrix"); }
        }
        for k in 0..self.axes[COLS].hdrs.len() {
            if self.hdr(COLS, k).is_none() { return Err("Singular Matrix"); }
        }
        self.set_state(MatrixState::FACTORING)?;

        for n in 0..self.diag.len() - 1 {
            let pivot = match self.search_for_pivot(n) {
                None => return Err("Pivot Search Fail"),
                Some(p) => p,
            };
            self.swap(ROWS, self[pivot].row, n);
            self.swap(COLS, self[pivot].col, n);
            self.row_col_elim(pivot, n)?;
        }
        self.set_state(MatrixState::FACTORED)?;
        return Ok(());
    }
    fn search_for_pivot(&self, n: usize) -> Option<Eindex> {
        let mut ei = self.markowitz_search_diagonal(n);
        if let Some(_) = ei { return ei; }
        ei = self.markowitz_search_submatrix(n);
        if let Some(_) = ei { return ei; }
        return self.find_max(n);
    }

    fn max_after(&self, ax: Axis, after: Eindex) -> Eindex {
        let mut best = after;
        let mut best_val = self[after].val.abs();
        let mut e = self[after].next(ax);

        while let Some(ei) = e {
            let val = self[ei].val.abs();
            if val > best_val {
                best = ei;
                best_val = val;
            }
            e = self[ei].next(ax);
        }
        return best;
    }

    fn markowitz_product(&self, ei: Eindex) -> usize {
        let e = &self[ei];
        let mr = self[Axis::ROWS].markowitz[e.row];
        let mc = self[Axis::COLS].markowitz[e.col];
        assert_gt!( mr, 0 );
        assert_gt!( mc, 0 );
        return (mr - 1) * (mc - 1);
    }

    fn markowitz_search_diagonal(&self, n: usize) -> Option<Eindex> {
        let rel_threshold = 1e-3;
        let abs_threshold = 0.0;
        let ties_mult = 5;

        let mut best_elem = None;
        let mut best_mark = MAX; // Actually use usize::MAX!
        let mut best_ratio = 0.0;
        let mut num_ties = 0;

        for k in n..self.diag.len() {
            let d = match self.diag[k] {
                None => { continue; }
                Some(d) => d,
            };

            // Check whether this element meets our threshold criteria
            let max_in_col = self.max_after(COLS, d);
            let threshold = rel_threshold * self[max_in_col].val.abs() + abs_threshold;
            if self[d].val.abs() < threshold { continue; }

            // If so, compute and compare its Markowitz product to our best
            let mark = self.markowitz_product(d);
            if mark < best_mark {
                num_ties = 0;
                best_elem = self.diag[k];
                best_mark = mark;
                best_ratio = (self[d].val / self[max_in_col].val).abs();
            } else if mark == best_mark {
                num_ties += 1;
                let ratio = (self[d].val / self[max_in_col].val).abs();
                if ratio > best_ratio {
                    best_elem = self.diag[k];
                    best_mark = mark;
                    best_ratio = ratio;
                }
                if num_ties >= best_mark * ties_mult { return best_elem; }
            }
        }
        return best_elem;
    }

    fn markowitz_search_submatrix(&self, n: usize) -> Option<Eindex> {
        //let rel_threshold = 1e-3;
        //let abs_threshold = 0.0;
        //let ties_mult = 5;

        let mut best_elem = None;
        let mut best_mark = MAX; // Actually use usize::MAX!
        let mut best_ratio = 0.0;
        //let mut num_ties = 0;

        for _k in n..self.axes[COLS].hdrs.len() {
            let mut e = self.hdr(COLS, n);
            // Advance to a row ≥ n
            while let Some(ei) = e {
                if self[ei].row >= n { break; }
                e = self[ei].next_in_col;
            }
            let ei = match e {
                None => { continue; }
                Some(d) => d,
            };

            // Check whether this element meets our threshold criteria
            let max_in_col = self.max_after(COLS, ei);
            //let threshold = rel_threshold * self[max_in_col].val.abs() + abs_threshold;

            while let Some(ei) = e {
                // If so, compute and compare its Markowitz product to our best
                let mark = self.markowitz_product(ei);
                if mark < best_mark {
                    //num_ties = 0;
                    best_elem = e;
                    best_mark = mark;
                    best_ratio = (self[ei].val / self[max_in_col].val).abs();
                } else if mark == best_mark {
                    //num_ties += 1;
                    let ratio = (self[ei].val / self[max_in_col].val).abs();
                    if ratio > best_ratio {
                        best_elem = e;
                        best_mark = mark;
                        best_ratio = ratio;
                    }
//                    // FIXME: do we want tie-counting in here?
//                    if num_ties >= best_mark * ties_mult { return best_elem; }
                }
                e = self[ei].next_in_col;
            }
        }
        return best_elem;
    }
    /// Find the max (abs value) element in sub-matrix of indices ≥ `n`.
    /// Returns `None` if no elements present. 
    fn find_max(&self, n: usize) -> Option<Eindex> {
        let mut max_elem = None;
        let mut max_val = 0.0;

        // Search each column ≥ n
        for k in n..self.axes[COLS].hdrs.len() {
            let mut ep = self.hdr(COLS, k);

            // Advance to a row ≥ n
            while let Some(ei) = ep {
                if self[ei].row >= n { break; }
                ep = self[ei].next_in_col;
            }
            // And search over remaining elements
            while let Some(ei) = ep {
                let val = self[ei].val.abs();
                if val > max_val {
                    max_elem = ep;
                    max_val = val;
                }
                ep = self[ei].next_in_col;
            }
        }
        return max_elem;
    }

    fn row_col_elim(&mut self, pivot: Eindex, n: usize) -> SpResult<()> {
        let de = match self.diag[n] {
            Some(de) => de,
            None => return Err("Singular Matrix"),
        };
        assert_eq!( de, pivot );
        let pivot_val = self[pivot].val;
        assert_ne!( pivot_val, 0.0 );

        // Divide elements in the pivot column by the pivot-value
        let mut plower = self[pivot].next_in_col;
        while let Some(ple) = plower {
            self[ple].val /= pivot_val;
            plower = self[ple].next_in_col;
        }

        let mut pupper = self[pivot].next_in_row;
        while let Some(pue) = pupper {
            let pupper_col = self[pue].col;
            plower = self[pivot].next_in_col;
            let mut psub = self[pue].next_in_col;
            while let Some(ple) = plower {
                // Walk `psub` down to the lower pointer
                while let Some(pse) = psub {
                    if self[pse].row >= self[ple].row { break; }
                    psub = self[pse].next_in_col;
                }
                let pse = match psub {
                    None => self.add_fillin(self[ple].row, pupper_col),
                    Some(pse) if self[pse].row > self[ple].row => {
                        self.add_fillin(self[ple].row, pupper_col)
                    }
                    Some(pse) => pse,
                };

                // Update the `psub` element value
                self[pse].val -= self[pue].val * self[ple].val;
                psub = self[pse].next_in_col;
                plower = self[ple].next_in_col;
            }
            self.axes[COLS].markowitz[pupper_col] -= 1;
            pupper = self[pue].next_in_row;
        }
        // Update remaining Markowitz counts
        self.axes[ROWS].markowitz[n] -= 1;
        self.axes[COLS].markowitz[n] -= 1;
        plower = self[pivot].next_in_col;
        while let Some(ple) = plower {
            let plower_row = self[ple].row;
            self.axes[ROWS].markowitz[plower_row] -= 1;
            plower = self[ple].next_in_col;
        }
        return Ok(());
    }
    /// Solve the system `Ax=b`, where: 
    /// * `A` is `self` 
    /// * `b` is argument `rhs`
    /// * `x` is the return value. 
    /// 
    /// Returns a `Result` containing the `Vec64` representing `x` if successful. 
    /// Returns an `Err` if unsuccessful. 
    /// 
    /// Performs LU factorization, forward and backward substitution. 
    pub fn solve(&mut self, rhs: Vec64 ) -> SpResult<Vec64> {
        let result = self._solve( rhs.vec );
        let vec = match result {
            Ok(vec) => vec,
            Err(_error) => return Err("Sparse solve method failed."),
        };
        let ans = Vec64::create( vec );
        Ok(ans)
    }  
    pub fn _solve(&mut self, rhs: Vec<f64>) -> SpResult<Vec<f64>> {
        if self.state == MatrixState::CREATED { self.lu_factorize()?; }
        assert_eq!( self.state, MatrixState::FACTORED );

        // Unwind any row-swaps
        let mut c: Vec<f64> = vec![0.0; rhs.len()];
        let row_mapping = self.axes[ROWS].mapping.as_ref().unwrap();
        for k in 0..c.len() {
            c[row_mapping.e2i[k]] = rhs[k];
        }

        // Forward substitution: Lc=b
        for k in 0..self.diag.len() {
            // Walk down each column, update c
            if c[k] == 0.0 { continue; } // No updates to make on this iteration

            let di = match self.diag[k] {
                Some(di) => di,
                None => return Err("Singular Matrix"),
            };
            let mut e = self[di].next_in_col;
            while let Some(ei) = e {
                c[self[ei].row] -= c[k] * self[ei].val;
                e = self[ei].next_in_col;
            }
        }

        // Backward substitution: Ux=c
        for k in (0..self.diag.len()).rev() {
            // Walk each row, update c
            let di = match self.diag[k] {
                Some(di) => di,
                None => return Err("Singular Matrix"),
            };
            let mut ep = self[di].next_in_row;
            while let Some(ei) = ep {
                c[k] -= c[self[ei].col] * self[ei].val;
                ep = self[ei].next_in_row;
            }
            c[k] /= self[di].val;
        }

        // Unwind any column-swaps
        let mut soln: Vec<f64> = vec![0.0; c.len()];
        let col_mapping = self.axes[COLS].mapping.as_ref().unwrap();
        for k in 0..c.len() {
            soln[k] = c[col_mapping.e2i[k]];
        }
        return Ok(soln);
    }
    pub fn hdr(&self, ax: Axis, loc: usize) -> Option<Eindex> { self.axes[ax].hdrs[loc] }
    fn set_hdr(&mut self, ax: Axis, loc: usize, ei: Option<Eindex>) { self.axes[ax].hdrs[loc] = ei; }
    pub fn num_rows(&self) -> usize { self.axes[ROWS].hdrs.len() }
    pub fn num_cols(&self) -> usize { self.axes[COLS].hdrs.len() }
    pub fn size(&self) -> (usize, usize) { (self.num_rows(), self.num_cols()) }
}

impl Index<Eindex> for Sparse {
    type Output = Element;
    fn index(&self, index: Eindex) -> &Self::Output { &self.elements[index.0] }
}

impl IndexMut<Eindex> for Sparse {
    fn index_mut(&mut self, index: Eindex) -> &mut Self::Output { &mut self.elements[index.0] }
}

impl Index<Axis> for Sparse {
    type Output = AxisData;
    fn index(&self, ax: Axis) -> &Self::Output { &self.axes[ax] }
}

impl IndexMut<Axis> for Sparse {
    fn index_mut(&mut self, ax: Axis) -> &mut Self::Output { &mut self.axes[ax] }
}