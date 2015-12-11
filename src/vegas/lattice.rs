//! Useful functions and data structures to build lattices

use util::super_mod;


#[derive(Debug)]
/// Represents a site, the cell represents where the site is
/// in the `Lattice` and the atom represents wich atom it is
/// within the unitcell
pub struct Site {
    cell: (i64, i64, i64),
    atom: u32,
}


impl Site {
    pub fn atom(&self) -> u32 {
        self.atom
    }
}


/// Represents a lattice, it requires the periodicity of the lattice,
/// its shape, its number of atoms per site and the vertices as a list
/// of vertex descriptors
pub struct Lattice {
    pbc: (bool, bool, bool),
    shape: (u32, u32, u32),
    natoms: u32,
    vertices: Vec<Vertex>,
}


impl Lattice {

    /// Returns a site in the first image of a lattice according the lattice's
    /// periodicity, returns none if the site is outside the lattice
    pub fn inside(&self, site: &Site) -> Option<Site> {

        if !(site.atom < self.natoms) {
            return None
        }

        let (mut x, mut y, mut z) = site.cell;
        let (sx, sy, sz) = (
            self.shape.0 as i64,
            self.shape.1 as i64,
            self.shape.2 as i64,
            );

        if !self.pbc.0  && (x < 0 || sx <= x) {
            return None
        } else {
            x = super_mod(x, sx);
        }

        if !self.pbc.1  && (y < 0 || sy <= y) {
            return None
        } else {
            y = super_mod(y, sy);
        }

        if !self.pbc.2  && (z < 0 || sz <= z) {
            return None
        } else {
            z = super_mod(z, sz);
        }

        Some(Site { cell: (x, y, z), atom: site.atom} )
    }

    pub fn sites(&self) -> SiteIterator {
        SiteIterator::new(self)
    }

    pub fn index(&self, site: &Site) -> Option<usize> {
        self.inside(site).map(|site| {
            let atom = site.atom as usize;
            let natm = self.natoms as usize;
            let (sx, sy) = (self.shape.0 as usize, self.shape.1 as usize);
            let (cx, cy, cz) = (
                site.cell.0 as usize,
                site.cell.1 as usize,
                site.cell.2 as usize,
                );
            natm * (sx * (sy * cz + cy) + cx) + atom
        })
    }

    pub fn tgts(&self, site: &Site) -> Option<Vec<(Site, Option<f64>)>> {
        self.inside(site).map(|site| {
            let mut tgts: Vec<(Site, Option<f64>)> = vec![];
            for vx in &self.vertices {
                match vx.tgt_for(&site) {
                    None => continue,
                    Some(tgt) => {
                        match self.inside(&tgt) {
                            None => continue,
                            Some(tgt) => {
                                tgts.push((tgt, vx.exch));
                            },
                        }
                    },
                }
            }
            tgts
        })
    }

    pub fn nsites(&self) -> usize {
        self.natoms as usize *
            self.shape.0 as usize *
            self.shape.1 as usize *
            self.shape.2 as usize
    }

    pub fn map_for_sites<T>(&self, op: fn(&Site) -> T) -> Vec<T> {
        let mut data: Vec<T> = vec![];
        for site in self.sites() {
            data.push(op(&site));
        }
        data
    }

}


/// A little builder for lattices
pub struct LatticeBuilder {
    pbc: (bool, bool, bool),
    shape: (u32, u32, u32),
    natoms: u32,
    vertices: Vec<Vertex>,
}


/// A little builder for lattices
impl LatticeBuilder {
    pub fn new() -> LatticeBuilder {
        LatticeBuilder {
            pbc: (true, true, true),
            shape: (10u32, 10u32, 10u32),
            natoms: 1u32,
            vertices: Vec::new(),
        }
    }

    pub fn pbc(mut self, pbc: (bool, bool, bool)) -> LatticeBuilder {
        self.pbc = pbc;
        self
    }

    pub fn shape(mut self, shape: (u32, u32, u32)) -> LatticeBuilder {
        self.shape = shape;
        self
    }

    pub fn natoms(mut self, natoms: u32) -> LatticeBuilder {
        self.natoms = natoms;
        self
    }

    pub fn vertices(mut self, vertices: Vec<Vertex>) -> LatticeBuilder {
        self.vertices = vertices;
        self
    }

    pub fn finalize(self) -> Lattice {
        Lattice {
            pbc: self.pbc,
            shape: self.shape,
            natoms: self.natoms,
            vertices: self.vertices,
        }
    }
}


/// Iterates over the cells of a lattice
struct CellIterator {
    cur: u32,
    max: (u32, u32, u32),
}


impl CellIterator {
    pub fn new(lattice: &Lattice) -> CellIterator {
        CellIterator { cur: 0, max: lattice.shape }
    }
}


impl Iterator for CellIterator {

    type Item = (u32, u32, u32);

    fn next(&mut self) -> Option<(u32, u32, u32)> {
        if self.cur == self.max.0 * self.max.1 * self.max.2 {
            return None;
        }
        let x =  self.cur % self.max.0;
        let y = (self.cur / self.max.0) % self.max.1;
        let z =  self.cur / self.max.0  / self.max.1 ;
        self.cur += 1;
        Some((x, y, z))
    }

}


/// Iterates over the sites of a lattice
pub struct SiteIterator {
    cell_it: CellIterator,
    cur_cell: Option<<CellIterator as Iterator>::Item>,
    cur_at: u32,
    max_at: u32,
}


impl SiteIterator {

    pub fn new(lattice: &Lattice) -> SiteIterator {
        let mut iter = CellIterator::new(lattice);
        SiteIterator {
            cur_cell: iter.next(),
            cell_it: iter,
            cur_at: 0,
            max_at: lattice.natoms,
        }
    }

}


impl Iterator for SiteIterator {

    type Item = Site;

    fn next(&mut self) -> Option<Site> {
        if self.max_at == 0 {
            return None
        }
        if self.cur_at == self.max_at {
            self.cur_at = 0;
            self.cur_cell = self.cell_it.next();
        }
        let at = self.cur_at;
        self.cur_at = self.cur_at + 1;
        match self.cur_cell {
            None => None,
            Some((x, y, z)) => Some(Site {
                cell: (x as i64, y as i64, z as i64),
                atom: at,
            }),
        }
    }

}


/// Represents a vertex descriptor, for a vertex that can go beyond the
/// unit cell of a lattice.
pub struct Vertex {
    src: u32,
    tgt: u32,
    delta: (i64, i64, i64),
    exch: Option<f64>,
}


impl Vertex {

    fn tgt_for(&self, site: &Site) -> Option<Site> {
        if site.atom != self.src {
            return None
        }
        Some(Site {
            cell: (site.cell.0 + self.delta.0,
                   site.cell.1 + self.delta.1,
                   site.cell.2 + self.delta.2),
                   atom: self.tgt,
        })
    }

    pub fn list_for_cubic() -> Vec<Vertex> {
        vec![
            Vertex { src: 0, tgt: 0, delta: (1, 0, 0), exch: None },
            Vertex { src: 0, tgt: 0, delta: (0, 1, 0), exch: None },
            Vertex { src: 0, tgt: 0, delta: (0, 0, 1), exch: None },
            Vertex { src: 0, tgt: 0, delta: (-1, 0, 0), exch: None },
            Vertex { src: 0, tgt: 0, delta: (0, -1, 0), exch: None },
            Vertex { src: 0, tgt: 0, delta: (0, 0, -1), exch: None },
        ]
    }

    pub fn list_for_hcp() -> Vec<Vertex> {
        vec![
            // Zero in plane
            Vertex { src: 0, tgt: 0, delta: (1, 0, 0), exch: None },
            Vertex { src: 0, tgt: 0, delta: (0, 1, 0), exch: None },
            Vertex { src: 0, tgt: 0, delta: (1, 1, 0), exch: None },
            // Zero in plane backwards
            Vertex { src: 0, tgt: 0, delta: (-1,  0, 0), exch: None },
            Vertex { src: 0, tgt: 0, delta: ( 0, -1, 0), exch: None },
            Vertex { src: 0, tgt: 0, delta: (-1, -1, 0), exch: None },
            // One in plane
            Vertex { src: 1, tgt: 1, delta: (1, 0, 0), exch: None },
            Vertex { src: 1, tgt: 1, delta: (0, 1, 0), exch: None },
            Vertex { src: 1, tgt: 1, delta: (1, 1, 0), exch: None },
            // One in plane backwards
            Vertex { src: 1, tgt: 1, delta: (-1,  0, 0), exch: None },
            Vertex { src: 1, tgt: 1, delta: ( 0, -1, 0), exch: None },
            Vertex { src: 1, tgt: 1, delta: (-1, -1, 0), exch: None },
            // Zero with one
            Vertex { src: 0, tgt: 1, delta: ( 0,  0, 0), exch: None },
            Vertex { src: 0, tgt: 1, delta: (-1,  0, 0), exch: None },
            Vertex { src: 0, tgt: 1, delta: (-1, -1, 0), exch: None },
            // Zero with one downwards
            Vertex { src: 0, tgt: 1, delta: ( 0,  0, -1), exch: None },
            Vertex { src: 0, tgt: 1, delta: (-1,  0, -1), exch: None },
            Vertex { src: 0, tgt: 1, delta: (-1, -1, -1), exch: None },
            // One with zero
            Vertex { src: 1, tgt: 0, delta: ( 0,  0, 0), exch: None },
            Vertex { src: 1, tgt: 0, delta: ( 1,  0, 0), exch: None },
            Vertex { src: 1, tgt: 0, delta: ( 1,  1, 0), exch: None },
            // Zero with one upwards
            Vertex { src: 1, tgt: 0, delta: ( 0,  0,  1), exch: None },
            Vertex { src: 1, tgt: 0, delta: ( 1,  0,  1), exch: None },
            Vertex { src: 1, tgt: 0, delta: ( 1,  1,  1), exch: None },
        ]
    }

    pub fn list_for_honeycomb() -> Vec<Vertex> {
        vec![
            Vertex { src: 0, tgt: 1, delta: (0, 0, 0), exch: None },
            Vertex { src: 1, tgt: 0, delta: (0, 0, 0), exch: None },
            Vertex { src: 0, tgt: 1, delta: (1, 0, 0), exch: None },
            Vertex { src: 1, tgt: 0, delta: (0, 1, 0), exch: None },
            Vertex { src: 0, tgt: 1, delta: (-1, 0, 0), exch: None },
            Vertex { src: 1, tgt: 0, delta: (0, -1, 0), exch: None },

            // For the 3D lulz
            Vertex { src: 0, tgt: 0, delta: (0, 0, 1), exch: None },
            Vertex { src: 1, tgt: 1, delta: (0, 0, 1), exch: None },
            Vertex { src: 0, tgt: 0, delta: (0, 0, -1), exch: None },
            Vertex { src: 1, tgt: 1, delta: (0, 0, -1), exch: None },
        ]
    }

    pub fn list_for_manganite() -> Vec<Vertex> {
        vec![
            Vertex { src: 0,  tgt: 2,  delta: (-1, 0,  0,  ), exch: Some(4.65), },
            Vertex { src: 0,  tgt: 6,  delta: (0,  -1, 0,  ), exch: Some(4.65), },
            Vertex { src: 0,  tgt: 18, delta: (0,  0,  -1, ), exch: Some(4.65), },
            Vertex { src: 0,  tgt: 1,  delta: (0,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 0,  tgt: 3,  delta: (0,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 0,  tgt: 9,  delta: (0,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 1,  tgt: 7,  delta: (0,  -1, 0,  ), exch: Some(1.35), },
            Vertex { src: 1,  tgt: 19, delta: (0,  0,  -1, ), exch: Some(1.35), },
            Vertex { src: 1,  tgt: 0,  delta: (0,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 1,  tgt: 2,  delta: (0,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 1,  tgt: 4,  delta: (0,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 1,  tgt: 10, delta: (0,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 2,  tgt: 8,  delta: (0,  -1, 0,  ), exch: Some(7.77), },
            Vertex { src: 2,  tgt: 20, delta: (0,  0,  -1, ), exch: Some(7.77), },
            Vertex { src: 2,  tgt: 1,  delta: (0,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 2,  tgt: 5,  delta: (0,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 2,  tgt: 11, delta: (0,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 2,  tgt: 0,  delta: (1,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 3,  tgt: 5,  delta: (-1, 0,  0,  ), exch: Some(1.35), },
            Vertex { src: 3,  tgt: 21, delta: (0,  0,  -1, ), exch: Some(1.35), },
            Vertex { src: 3,  tgt: 0,  delta: (0,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 3,  tgt: 4,  delta: (0,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 3,  tgt: 6,  delta: (0,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 3,  tgt: 12, delta: (0,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 4,  tgt: 22, delta: (0,  0,  -1, ), exch: Some(7.77), },
            Vertex { src: 4,  tgt: 1,  delta: (0,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 4,  tgt: 3,  delta: (0,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 4,  tgt: 5,  delta: (0,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 4,  tgt: 7,  delta: (0,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 4,  tgt: 13, delta: (0,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 5,  tgt: 23, delta: (0,  0,  -1, ), exch: Some(4.65), },
            Vertex { src: 5,  tgt: 2,  delta: (0,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 5,  tgt: 4,  delta: (0,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 5,  tgt: 8,  delta: (0,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 5,  tgt: 14, delta: (0,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 5,  tgt: 3,  delta: (1,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 6,  tgt: 8,  delta: (-1, 0,  0,  ), exch: Some(7.77), },
            Vertex { src: 6,  tgt: 24, delta: (0,  0,  -1, ), exch: Some(7.77), },
            Vertex { src: 6,  tgt: 3,  delta: (0,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 6,  tgt: 7,  delta: (0,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 6,  tgt: 15, delta: (0,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 6,  tgt: 0,  delta: (0,  1,  0,  ), exch: Some(4.65), },
            Vertex { src: 7,  tgt: 25, delta: (0,  0,  -1, ), exch: Some(4.65), },
            Vertex { src: 7,  tgt: 4,  delta: (0,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 7,  tgt: 6,  delta: (0,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 7,  tgt: 8,  delta: (0,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 7,  tgt: 16, delta: (0,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 7,  tgt: 1,  delta: (0,  1,  0,  ), exch: Some(1.35), },
            Vertex { src: 8,  tgt: 26, delta: (0,  0,  -1, ), exch: Some(1.35), },
            Vertex { src: 8,  tgt: 5,  delta: (0,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 8,  tgt: 7,  delta: (0,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 8,  tgt: 17, delta: (0,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 8,  tgt: 2,  delta: (0,  1,  0,  ), exch: Some(7.77), },
            Vertex { src: 8,  tgt: 6,  delta: (1,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 9,  tgt: 11, delta: (-1, 0,  0,  ), exch: Some(1.35), },
            Vertex { src: 9,  tgt: 15, delta: (0,  -1, 0,  ), exch: Some(1.35), },
            Vertex { src: 9,  tgt: 0,  delta: (0,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 9,  tgt: 10, delta: (0,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 9,  tgt: 12, delta: (0,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 9,  tgt: 18, delta: (0,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 10, tgt: 16, delta: (0,  -1, 0,  ), exch: Some(7.77), },
            Vertex { src: 10, tgt: 1,  delta: (0,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 10, tgt: 9,  delta: (0,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 10, tgt: 11, delta: (0,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 10, tgt: 13, delta: (0,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 10, tgt: 19, delta: (0,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 11, tgt: 17, delta: (0,  -1, 0,  ), exch: Some(4.65), },
            Vertex { src: 11, tgt: 2,  delta: (0,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 11, tgt: 10, delta: (0,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 11, tgt: 14, delta: (0,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 11, tgt: 20, delta: (0,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 11, tgt: 9,  delta: (1,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 12, tgt: 14, delta: (-1, 0,  0,  ), exch: Some(7.77), },
            Vertex { src: 12, tgt: 3,  delta: (0,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 12, tgt: 9,  delta: (0,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 12, tgt: 13, delta: (0,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 12, tgt: 15, delta: (0,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 12, tgt: 21, delta: (0,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 13, tgt: 4,  delta: (0,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 13, tgt: 10, delta: (0,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 13, tgt: 12, delta: (0,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 13, tgt: 14, delta: (0,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 13, tgt: 16, delta: (0,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 13, tgt: 22, delta: (0,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 14, tgt: 5,  delta: (0,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 14, tgt: 11, delta: (0,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 14, tgt: 13, delta: (0,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 14, tgt: 17, delta: (0,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 14, tgt: 23, delta: (0,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 14, tgt: 12, delta: (1,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 15, tgt: 17, delta: (-1, 0,  0,  ), exch: Some(4.65), },
            Vertex { src: 15, tgt: 6,  delta: (0,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 15, tgt: 12, delta: (0,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 15, tgt: 16, delta: (0,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 15, tgt: 24, delta: (0,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 15, tgt: 9,  delta: (0,  1,  0,  ), exch: Some(1.35), },
            Vertex { src: 16, tgt: 7,  delta: (0,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 16, tgt: 13, delta: (0,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 16, tgt: 15, delta: (0,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 16, tgt: 17, delta: (0,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 16, tgt: 25, delta: (0,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 16, tgt: 10, delta: (0,  1,  0,  ), exch: Some(7.77), },
            Vertex { src: 17, tgt: 8,  delta: (0,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 17, tgt: 14, delta: (0,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 17, tgt: 16, delta: (0,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 17, tgt: 26, delta: (0,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 17, tgt: 11, delta: (0,  1,  0,  ), exch: Some(4.65), },
            Vertex { src: 17, tgt: 15, delta: (1,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 18, tgt: 20, delta: (-1, 0,  0,  ), exch: Some(7.77), },
            Vertex { src: 18, tgt: 24, delta: (0,  -1, 0,  ), exch: Some(7.77), },
            Vertex { src: 18, tgt: 9,  delta: (0,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 18, tgt: 19, delta: (0,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 18, tgt: 21, delta: (0,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 18, tgt: 0,  delta: (0,  0,  1,  ), exch: Some(4.65), },
            Vertex { src: 19, tgt: 25, delta: (0,  -1, 0,  ), exch: Some(4.65), },
            Vertex { src: 19, tgt: 10, delta: (0,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 19, tgt: 18, delta: (0,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 19, tgt: 20, delta: (0,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 19, tgt: 22, delta: (0,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 19, tgt: 1,  delta: (0,  0,  1,  ), exch: Some(1.35), },
            Vertex { src: 20, tgt: 26, delta: (0,  -1, 0,  ), exch: Some(1.35), },
            Vertex { src: 20, tgt: 11, delta: (0,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 20, tgt: 19, delta: (0,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 20, tgt: 23, delta: (0,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 20, tgt: 2,  delta: (0,  0,  1,  ), exch: Some(7.77), },
            Vertex { src: 20, tgt: 18, delta: (1,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 21, tgt: 23, delta: (-1, 0,  0,  ), exch: Some(4.65), },
            Vertex { src: 21, tgt: 12, delta: (0,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 21, tgt: 18, delta: (0,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 21, tgt: 22, delta: (0,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 21, tgt: 24, delta: (0,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 21, tgt: 3,  delta: (0,  0,  1,  ), exch: Some(1.35), },
            Vertex { src: 22, tgt: 13, delta: (0,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 22, tgt: 19, delta: (0,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 22, tgt: 21, delta: (0,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 22, tgt: 23, delta: (0,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 22, tgt: 25, delta: (0,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 22, tgt: 4,  delta: (0,  0,  1,  ), exch: Some(7.77), },
            Vertex { src: 23, tgt: 14, delta: (0,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 23, tgt: 20, delta: (0,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 23, tgt: 22, delta: (0,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 23, tgt: 26, delta: (0,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 23, tgt: 5,  delta: (0,  0,  1,  ), exch: Some(4.65), },
            Vertex { src: 23, tgt: 21, delta: (1,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 24, tgt: 26, delta: (-1, 0,  0,  ), exch: Some(1.35), },
            Vertex { src: 24, tgt: 15, delta: (0,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 24, tgt: 21, delta: (0,  0,  0,  ), exch: Some(1.35), },
            Vertex { src: 24, tgt: 25, delta: (0,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 24, tgt: 6,  delta: (0,  0,  1,  ), exch: Some(7.77), },
            Vertex { src: 24, tgt: 18, delta: (0,  1,  0,  ), exch: Some(7.77), },
            Vertex { src: 25, tgt: 16, delta: (0,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 25, tgt: 22, delta: (0,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 25, tgt: 24, delta: (0,  0,  0,  ), exch: Some(7.77), },
            Vertex { src: 25, tgt: 26, delta: (0,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 25, tgt: 7,  delta: (0,  0,  1,  ), exch: Some(4.65), },
            Vertex { src: 25, tgt: 19, delta: (0,  1,  0,  ), exch: Some(4.65), },
            Vertex { src: 26, tgt: 17, delta: (0,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 26, tgt: 23, delta: (0,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 26, tgt: 25, delta: (0,  0,  0,  ), exch: Some(4.65), },
            Vertex { src: 26, tgt: 8,  delta: (0,  0,  1,  ), exch: Some(1.35), },
            Vertex { src: 26, tgt: 20, delta: (0,  1,  0,  ), exch: Some(1.35), },
            Vertex { src: 26, tgt: 24, delta: (1,  0,  0,  ), exch: Some(1.35), },
        ]
    }

    pub fn list_for_magnetite() -> Vec<Vertex> {
        vec![
            Vertex { src: 0,   tgt: 19,  delta: (-1,  -1,  -1, ), exch: Some(0.11), },
            Vertex { src: 0,   tgt: 4,   delta: (-1,  -1,  0 , ), exch: Some(2.92), },
            Vertex { src: 0,   tgt: 5,   delta: (-1,  -1,  0 , ), exch: Some(2.92), },
            Vertex { src: 0,   tgt: 11,  delta: (-1,  -1,  0 , ), exch: Some(2.92), },
            Vertex { src: 0,   tgt: 17,  delta: (-1,  0 ,  -1, ), exch: Some(2.92), },
            Vertex { src: 0,   tgt: 22,  delta: (-1,  0 ,  -1, ), exch: Some(2.92), },
            Vertex { src: 0,   tgt: 23,  delta: (-1,  0 ,  -1, ), exch: Some(2.92), },
            Vertex { src: 0,   tgt: 7,   delta: (-1,  0 ,  0 , ), exch: Some(0.11), },
            Vertex { src: 0,   tgt: 14,  delta: (0,   -1,  -1, ), exch: Some(2.92), },
            Vertex { src: 0,   tgt: 20,  delta: (0,   -1,  -1, ), exch: Some(2.92), },
            Vertex { src: 0,   tgt: 21,  delta: (0,   -1,  -1, ), exch: Some(2.92), },
            Vertex { src: 0,   tgt: 6,   delta: (0,   -1,  0 , ), exch: Some(0.11), },
            Vertex { src: 0,   tgt: 18,  delta: (0,   0 ,  -1, ), exch: Some(0.11), },
            Vertex { src: 0,   tgt: 2,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 0,   tgt: 3,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 0,   tgt: 8,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 1,   tgt: 15,  delta: (0,   0 ,  -1, ), exch: Some(2.92), },
            Vertex { src: 1,   tgt: 16,  delta: (0,   0 ,  -1, ), exch: Some(2.92), },
            Vertex { src: 1,   tgt: 18,  delta: (0,   0 ,  -1, ), exch: Some(0.11), },
            Vertex { src: 1,   tgt: 19,  delta: (0,   0 ,  -1, ), exch: Some(0.11), },
            Vertex { src: 1,   tgt: 20,  delta: (0,   0 ,  -1, ), exch: Some(2.92), },
            Vertex { src: 1,   tgt: 21,  delta: (0,   0 ,  -1, ), exch: Some(2.92), },
            Vertex { src: 1,   tgt: 22,  delta: (0,   0 ,  -1, ), exch: Some(2.92), },
            Vertex { src: 1,   tgt: 23,  delta: (0,   0 ,  -1, ), exch: Some(2.92), },
            Vertex { src: 1,   tgt: 2,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 1,   tgt: 3,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 1,   tgt: 4,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 1,   tgt: 5,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 1,   tgt: 6,   delta: (0,   0 ,  0 , ), exch: Some(0.11), },
            Vertex { src: 1,   tgt: 7,   delta: (0,   0 ,  0 , ), exch: Some(0.11), },
            Vertex { src: 1,   tgt: 9,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 1,   tgt: 10,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 2,   tgt: 23,  delta: (-1,  0 ,  -1, ), exch: Some(-0.63), },
            Vertex { src: 2,   tgt: 5,   delta: (-1,  0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 2,   tgt: 7,   delta: (-1,  0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 2,   tgt: 18,  delta: (0,   0 ,  -1, ), exch: Some(2.92), },
            Vertex { src: 2,   tgt: 20,  delta: (0,   0 ,  -1, ), exch: Some(-0.63), },
            Vertex { src: 2,   tgt: 0,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 2,   tgt: 1,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 2,   tgt: 3,   delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 2,   tgt: 6,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 2,   tgt: 8,   delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 2,   tgt: 9,   delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 2,   tgt: 12,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 3,   tgt: 21,  delta: (0,   -1,  -1, ), exch: Some(-0.63), },
            Vertex { src: 3,   tgt: 4,   delta: (0,   -1,  0 , ), exch: Some(-0.63), },
            Vertex { src: 3,   tgt: 6,   delta: (0,   -1,  0 , ), exch: Some(2.92), },
            Vertex { src: 3,   tgt: 18,  delta: (0,   0 ,  -1, ), exch: Some(2.92), },
            Vertex { src: 3,   tgt: 22,  delta: (0,   0 ,  -1, ), exch: Some(-0.63), },
            Vertex { src: 3,   tgt: 0,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 3,   tgt: 1,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 3,   tgt: 2,   delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 3,   tgt: 7,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 3,   tgt: 8,   delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 3,   tgt: 9,   delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 3,   tgt: 13,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 4,   tgt: 19,  delta: (0,   0 ,  -1, ), exch: Some(2.92), },
            Vertex { src: 4,   tgt: 21,  delta: (0,   0 ,  -1, ), exch: Some(-0.63), },
            Vertex { src: 4,   tgt: 1,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 4,   tgt: 5,   delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 4,   tgt: 6,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 4,   tgt: 10,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 4,   tgt: 11,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 4,   tgt: 22,  delta: (0,   1 ,  -1, ), exch: Some(-0.63), },
            Vertex { src: 4,   tgt: 3,   delta: (0,   1 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 4,   tgt: 7,   delta: (0,   1 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 4,   tgt: 13,  delta: (0,   1 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 4,   tgt: 0,   delta: (1,   1 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 5,   tgt: 19,  delta: (0,   0 ,  -1, ), exch: Some(2.92), },
            Vertex { src: 5,   tgt: 23,  delta: (0,   0 ,  -1, ), exch: Some(-0.63), },
            Vertex { src: 5,   tgt: 1,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 5,   tgt: 4,   delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 5,   tgt: 7,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 5,   tgt: 10,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 5,   tgt: 11,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 5,   tgt: 20,  delta: (1,   0 ,  -1, ), exch: Some(-0.63), },
            Vertex { src: 5,   tgt: 2,   delta: (1,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 5,   tgt: 6,   delta: (1,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 5,   tgt: 12,  delta: (1,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 5,   tgt: 0,   delta: (1,   1 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 6,   tgt: 5,   delta: (-1,  0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 6,   tgt: 11,  delta: (-1,  0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 6,   tgt: 20,  delta: (0,   0 ,  -1, ), exch: Some(2.92), },
            Vertex { src: 6,   tgt: 21,  delta: (0,   0 ,  -1, ), exch: Some(2.92), },
            Vertex { src: 6,   tgt: 1,   delta: (0,   0 ,  0 , ), exch: Some(0.11), },
            Vertex { src: 6,   tgt: 2,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 6,   tgt: 4,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 6,   tgt: 9,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 6,   tgt: 10,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 6,   tgt: 12,  delta: (0,   0 ,  0 , ), exch: Some(0.11), },
            Vertex { src: 6,   tgt: 14,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 6,   tgt: 15,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 6,   tgt: 0,   delta: (0,   1 ,  0 , ), exch: Some(0.11), },
            Vertex { src: 6,   tgt: 3,   delta: (0,   1 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 6,   tgt: 8,   delta: (0,   1 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 6,   tgt: 13,  delta: (0,   1 ,  0 , ), exch: Some(0.11), },
            Vertex { src: 7,   tgt: 4,   delta: (0,   -1,  0 , ), exch: Some(2.92), },
            Vertex { src: 7,   tgt: 11,  delta: (0,   -1,  0 , ), exch: Some(2.92), },
            Vertex { src: 7,   tgt: 22,  delta: (0,   0 ,  -1, ), exch: Some(2.92), },
            Vertex { src: 7,   tgt: 23,  delta: (0,   0 ,  -1, ), exch: Some(2.92), },
            Vertex { src: 7,   tgt: 1,   delta: (0,   0 ,  0 , ), exch: Some(0.11), },
            Vertex { src: 7,   tgt: 3,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 7,   tgt: 5,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 7,   tgt: 9,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 7,   tgt: 10,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 7,   tgt: 13,  delta: (0,   0 ,  0 , ), exch: Some(0.11), },
            Vertex { src: 7,   tgt: 16,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 7,   tgt: 17,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 7,   tgt: 0,   delta: (1,   0 ,  0 , ), exch: Some(0.11), },
            Vertex { src: 7,   tgt: 2,   delta: (1,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 7,   tgt: 8,   delta: (1,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 7,   tgt: 12,  delta: (1,   0 ,  0 , ), exch: Some(0.11), },
            Vertex { src: 8,   tgt: 11,  delta: (-1,  -1,  0 , ), exch: Some(-0.63), },
            Vertex { src: 8,   tgt: 7,   delta: (-1,  0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 8,   tgt: 17,  delta: (-1,  0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 8,   tgt: 6,   delta: (0,   -1,  0 , ), exch: Some(2.92), },
            Vertex { src: 8,   tgt: 14,  delta: (0,   -1,  0 , ), exch: Some(-0.63), },
            Vertex { src: 8,   tgt: 0,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 8,   tgt: 2,   delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 8,   tgt: 3,   delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 8,   tgt: 9,   delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 8,   tgt: 12,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 8,   tgt: 13,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 8,   tgt: 18,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 9,   tgt: 1,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 9,   tgt: 2,   delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 9,   tgt: 3,   delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 9,   tgt: 6,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 9,   tgt: 7,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 9,   tgt: 8,   delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 9,   tgt: 10,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 9,   tgt: 12,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 9,   tgt: 13,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 9,   tgt: 15,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 9,   tgt: 16,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 9,   tgt: 18,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 10,  tgt: 1,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 10,  tgt: 4,   delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 10,  tgt: 5,   delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 10,  tgt: 6,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 10,  tgt: 7,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 10,  tgt: 9,   delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 10,  tgt: 11,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 10,  tgt: 15,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 10,  tgt: 16,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 10,  tgt: 19,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 10,  tgt: 13,  delta: (0,   1 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 10,  tgt: 12,  delta: (1,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 11,  tgt: 4,   delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 11,  tgt: 5,   delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 11,  tgt: 10,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 11,  tgt: 19,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 11,  tgt: 7,   delta: (0,   1 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 11,  tgt: 13,  delta: (0,   1 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 11,  tgt: 17,  delta: (0,   1 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 11,  tgt: 6,   delta: (1,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 11,  tgt: 12,  delta: (1,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 11,  tgt: 14,  delta: (1,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 11,  tgt: 0,   delta: (1,   1 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 11,  tgt: 8,   delta: (1,   1 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 12,  tgt: 5,   delta: (-1,  0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 12,  tgt: 7,   delta: (-1,  0 ,  0 , ), exch: Some(0.11), },
            Vertex { src: 12,  tgt: 10,  delta: (-1,  0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 12,  tgt: 11,  delta: (-1,  0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 12,  tgt: 16,  delta: (-1,  0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 12,  tgt: 17,  delta: (-1,  0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 12,  tgt: 19,  delta: (-1,  0 ,  0 , ), exch: Some(0.11), },
            Vertex { src: 12,  tgt: 23,  delta: (-1,  0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 12,  tgt: 2,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 12,  tgt: 6,   delta: (0,   0 ,  0 , ), exch: Some(0.11), },
            Vertex { src: 12,  tgt: 8,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 12,  tgt: 9,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 12,  tgt: 14,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 12,  tgt: 15,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 12,  tgt: 18,  delta: (0,   0 ,  0 , ), exch: Some(0.11), },
            Vertex { src: 12,  tgt: 20,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 13,  tgt: 4,   delta: (0,   -1,  0 , ), exch: Some(2.92), },
            Vertex { src: 13,  tgt: 6,   delta: (0,   -1,  0 , ), exch: Some(0.11), },
            Vertex { src: 13,  tgt: 10,  delta: (0,   -1,  0 , ), exch: Some(2.92), },
            Vertex { src: 13,  tgt: 11,  delta: (0,   -1,  0 , ), exch: Some(2.92), },
            Vertex { src: 13,  tgt: 14,  delta: (0,   -1,  0 , ), exch: Some(2.92), },
            Vertex { src: 13,  tgt: 15,  delta: (0,   -1,  0 , ), exch: Some(2.92), },
            Vertex { src: 13,  tgt: 19,  delta: (0,   -1,  0 , ), exch: Some(0.11), },
            Vertex { src: 13,  tgt: 21,  delta: (0,   -1,  0 , ), exch: Some(2.92), },
            Vertex { src: 13,  tgt: 3,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 13,  tgt: 7,   delta: (0,   0 ,  0 , ), exch: Some(0.11), },
            Vertex { src: 13,  tgt: 8,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 13,  tgt: 9,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 13,  tgt: 16,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 13,  tgt: 17,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 13,  tgt: 18,  delta: (0,   0 ,  0 , ), exch: Some(0.11), },
            Vertex { src: 13,  tgt: 22,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 14,  tgt: 11,  delta: (-1,  0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 14,  tgt: 19,  delta: (-1,  0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 14,  tgt: 17,  delta: (-1,  1 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 14,  tgt: 6,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 14,  tgt: 12,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 14,  tgt: 15,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 14,  tgt: 20,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 14,  tgt: 21,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 14,  tgt: 8,   delta: (0,   1 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 14,  tgt: 13,  delta: (0,   1 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 14,  tgt: 18,  delta: (0,   1 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 14,  tgt: 0,   delta: (0,   1 ,  1 , ), exch: Some(2.92), },
            Vertex { src: 15,  tgt: 6,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 15,  tgt: 9,   delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 15,  tgt: 10,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 15,  tgt: 12,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 15,  tgt: 14,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 15,  tgt: 16,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 15,  tgt: 18,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 15,  tgt: 19,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 15,  tgt: 20,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 15,  tgt: 21,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 15,  tgt: 1,   delta: (0,   0 ,  1 , ), exch: Some(2.92), },
            Vertex { src: 15,  tgt: 13,  delta: (0,   1 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 16,  tgt: 7,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 16,  tgt: 9,   delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 16,  tgt: 10,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 16,  tgt: 13,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 16,  tgt: 15,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 16,  tgt: 17,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 16,  tgt: 18,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 16,  tgt: 19,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 16,  tgt: 22,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 16,  tgt: 23,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 16,  tgt: 1,   delta: (0,   0 ,  1 , ), exch: Some(2.92), },
            Vertex { src: 16,  tgt: 12,  delta: (1,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 17,  tgt: 11,  delta: (0,   -1,  0 , ), exch: Some(-0.63), },
            Vertex { src: 17,  tgt: 19,  delta: (0,   -1,  0 , ), exch: Some(2.92), },
            Vertex { src: 17,  tgt: 7,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 17,  tgt: 13,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 17,  tgt: 16,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 17,  tgt: 22,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 17,  tgt: 23,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 17,  tgt: 14,  delta: (1,   -1,  0 , ), exch: Some(-0.63), },
            Vertex { src: 17,  tgt: 8,   delta: (1,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 17,  tgt: 12,  delta: (1,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 17,  tgt: 18,  delta: (1,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 17,  tgt: 0,   delta: (1,   0 ,  1 , ), exch: Some(2.92), },
            Vertex { src: 18,  tgt: 17,  delta: (-1,  0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 18,  tgt: 23,  delta: (-1,  0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 18,  tgt: 14,  delta: (0,   -1,  0 , ), exch: Some(2.92), },
            Vertex { src: 18,  tgt: 21,  delta: (0,   -1,  0 , ), exch: Some(2.92), },
            Vertex { src: 18,  tgt: 8,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 18,  tgt: 9,   delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 18,  tgt: 12,  delta: (0,   0 ,  0 , ), exch: Some(0.11), },
            Vertex { src: 18,  tgt: 13,  delta: (0,   0 ,  0 , ), exch: Some(0.11), },
            Vertex { src: 18,  tgt: 15,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 18,  tgt: 16,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 18,  tgt: 20,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 18,  tgt: 22,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 18,  tgt: 0,   delta: (0,   0 ,  1 , ), exch: Some(0.11), },
            Vertex { src: 18,  tgt: 1,   delta: (0,   0 ,  1 , ), exch: Some(0.11), },
            Vertex { src: 18,  tgt: 2,   delta: (0,   0 ,  1 , ), exch: Some(2.92), },
            Vertex { src: 18,  tgt: 3,   delta: (0,   0 ,  1 , ), exch: Some(2.92), },
            Vertex { src: 19,  tgt: 10,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 19,  tgt: 11,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 19,  tgt: 15,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 19,  tgt: 16,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 19,  tgt: 21,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 19,  tgt: 23,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 19,  tgt: 1,   delta: (0,   0 ,  1 , ), exch: Some(0.11), },
            Vertex { src: 19,  tgt: 4,   delta: (0,   0 ,  1 , ), exch: Some(2.92), },
            Vertex { src: 19,  tgt: 5,   delta: (0,   0 ,  1 , ), exch: Some(2.92), },
            Vertex { src: 19,  tgt: 13,  delta: (0,   1 ,  0 , ), exch: Some(0.11), },
            Vertex { src: 19,  tgt: 17,  delta: (0,   1 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 19,  tgt: 22,  delta: (0,   1 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 19,  tgt: 12,  delta: (1,   0 ,  0 , ), exch: Some(0.11), },
            Vertex { src: 19,  tgt: 14,  delta: (1,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 19,  tgt: 20,  delta: (1,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 19,  tgt: 0,   delta: (1,   1 ,  1 , ), exch: Some(0.11), },
            Vertex { src: 20,  tgt: 19,  delta: (-1,  0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 20,  tgt: 23,  delta: (-1,  0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 20,  tgt: 5,   delta: (-1,  0 ,  1 , ), exch: Some(-0.63), },
            Vertex { src: 20,  tgt: 12,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 20,  tgt: 14,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 20,  tgt: 15,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 20,  tgt: 18,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 20,  tgt: 21,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 20,  tgt: 1,   delta: (0,   0 ,  1 , ), exch: Some(2.92), },
            Vertex { src: 20,  tgt: 2,   delta: (0,   0 ,  1 , ), exch: Some(-0.63), },
            Vertex { src: 20,  tgt: 6,   delta: (0,   0 ,  1 , ), exch: Some(2.92), },
            Vertex { src: 20,  tgt: 0,   delta: (0,   1 ,  1 , ), exch: Some(2.92), },
            Vertex { src: 21,  tgt: 14,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 21,  tgt: 15,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 21,  tgt: 19,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 21,  tgt: 20,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 21,  tgt: 1,   delta: (0,   0 ,  1 , ), exch: Some(2.92), },
            Vertex { src: 21,  tgt: 4,   delta: (0,   0 ,  1 , ), exch: Some(-0.63), },
            Vertex { src: 21,  tgt: 6,   delta: (0,   0 ,  1 , ), exch: Some(2.92), },
            Vertex { src: 21,  tgt: 13,  delta: (0,   1 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 21,  tgt: 18,  delta: (0,   1 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 21,  tgt: 22,  delta: (0,   1 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 21,  tgt: 0,   delta: (0,   1 ,  1 , ), exch: Some(2.92), },
            Vertex { src: 21,  tgt: 3,   delta: (0,   1 ,  1 , ), exch: Some(-0.63), },
            Vertex { src: 22,  tgt: 19,  delta: (0,   -1,  0 , ), exch: Some(2.92), },
            Vertex { src: 22,  tgt: 21,  delta: (0,   -1,  0 , ), exch: Some(-0.63), },
            Vertex { src: 22,  tgt: 4,   delta: (0,   -1,  1 , ), exch: Some(-0.63), },
            Vertex { src: 22,  tgt: 13,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 22,  tgt: 16,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 22,  tgt: 17,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 22,  tgt: 18,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 22,  tgt: 23,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 22,  tgt: 1,   delta: (0,   0 ,  1 , ), exch: Some(2.92), },
            Vertex { src: 22,  tgt: 3,   delta: (0,   0 ,  1 , ), exch: Some(-0.63), },
            Vertex { src: 22,  tgt: 7,   delta: (0,   0 ,  1 , ), exch: Some(2.92), },
            Vertex { src: 22,  tgt: 0,   delta: (1,   0 ,  1 , ), exch: Some(2.92), },
            Vertex { src: 23,  tgt: 16,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 23,  tgt: 17,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 23,  tgt: 19,  delta: (0,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 23,  tgt: 22,  delta: (0,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 23,  tgt: 1,   delta: (0,   0 ,  1 , ), exch: Some(2.92), },
            Vertex { src: 23,  tgt: 5,   delta: (0,   0 ,  1 , ), exch: Some(-0.63), },
            Vertex { src: 23,  tgt: 7,   delta: (0,   0 ,  1 , ), exch: Some(2.92), },
            Vertex { src: 23,  tgt: 12,  delta: (1,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 23,  tgt: 18,  delta: (1,   0 ,  0 , ), exch: Some(2.92), },
            Vertex { src: 23,  tgt: 20,  delta: (1,   0 ,  0 , ), exch: Some(-0.63), },
            Vertex { src: 23,  tgt: 0,   delta: (1,   0 ,  1 , ), exch: Some(2.92), },
            Vertex { src: 23,  tgt: 2,   delta: (1,   0 ,  1 , ), exch: Some(-0.63), },
            ]
    }
}


pub struct Adjacency {
    lims: Vec<usize>,
    nbhs: Vec<usize>,
    exch: Vec<f64>,
}


impl Adjacency {
    pub fn new(lattice: &Lattice) -> Adjacency
    {
        let mut lims = vec![0];
        let mut nbhs = vec![];
        let mut exch = vec![];
        for site in lattice.sites() {
            let targets = lattice.tgts(&site).unwrap();
            let mut pnbhs: Vec<usize> = vec![];
            let mut pexch: Vec<f64> = vec![];
            for (site, exchange) in targets {
                pnbhs.push(lattice.index(&site).unwrap());
                pexch.push(exchange.unwrap_or(0.0));
            }
            let last = lims.last().unwrap().clone();
            lims.push(last + pnbhs.len());
            nbhs.extend(pnbhs.iter());
            exch.extend(pexch.iter());
        }
        Adjacency { lims: lims, nbhs: nbhs, exch: exch, }
    }

    pub fn nbhs_of<'a>(&'a self, item: usize) -> Option<&'a[usize]> {
        if item >= self.lims.len() - 1 {
            return None
        }
        let low = self.lims[item] as usize;
        let hi = self.lims[item + 1] as usize;
        Some(&self.nbhs[low..hi])
    }

    pub fn exch_of<'a>(&'a self, item: usize) -> Option<&'a[f64]> {
        if item >= self.lims.len() - 1 {
            return None
        }
        let low = self.lims[item] as usize;
        let hi = self.lims[item + 1] as usize;
        Some(&self.exch[low..hi])
    }
}


struct Locator {
    a1: (f64, f64, f64),
    a2: (f64, f64, f64),
    a3: (f64, f64, f64),
    basis: Vec<(f64, f64, f64)>,
}


impl Locator {

    fn locate(&self, site: &Site) -> Option<(f64, f64, f64)> {
        let at = site.atom as usize;
        if at >= self.basis.len() {
            return None
        }
        let (bx, by, bz) = self.basis[at];
        let (cx, cy, cz) = (
            site.cell.0 as f64,
            site.cell.1 as f64,
            site.cell.2 as f64,
            );
        let pos = (
            cx * self.a1.0 + cy * self.a2.0 + cz * self.a3.0 + bx,
            cx * self.a1.1 + cy * self.a2.1 + cz * self.a3.1 + by,
            cx * self.a1.2 + cy * self.a2.2 + cz * self.a3.2 + bz
            );
        Some(pos)
    }

    pub fn for_cubic(a: f64) -> Locator {
        Locator {
            a1: (a, 0.0f64, 0.0f64),
            a2: (0.0f64, a, 0.0f64),
            a3: (0.0f64, 0.0f64, a),
            basis: vec![(0.0f64, 0.0f64, 0.0f64)],
        }
    }

}


// Tests

#[test]
fn testing_the_inside() {
    let latt = LatticeBuilder::new()
        .pbc((true, true, false))
        .vertices(Vertex::list_for_cubic())
        .finalize();
    assert!(latt.inside(&Site { cell: (10, 10, 9), atom: 0 }).is_some());
    assert!(latt.inside(&Site { cell: (10, -1, 9), atom: 0 }).is_some());
    assert!(latt.inside(&Site { cell: (10, 10, 10), atom: 0 }).is_none());
    assert!(latt.inside(&Site { cell: (10, 10, 9), atom: 2 }).is_none());
}

#[test]
fn test_big_magnetite() {
    let latt = LatticeBuilder::new()
        .pbc((true, true, true))
        .shape((4, 4, 4))
        .natoms(27)
        .vertices(Vertex::list_for_manganite())
        .finalize();
    let adj = Adjacency::new(&latt);
    assert_eq!(1729, adj.lims.len());
    assert_eq!(10368, adj.nbhs.len());
}
