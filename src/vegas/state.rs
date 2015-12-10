//! Describes states for spin systems, for now only the Heisenberg-like
//! state is implemented.

use std::ops::Mul;

extern crate rand;

use rand::distributions::normal::StandardNormal;


pub trait StateConstructors {
    fn up(size: usize) -> Self;
    fn rand(size: usize) -> Self;
    fn rand_with_norms(size: usize, norms: &Vec<f64>) -> Self;
}


pub trait SpinConstructors {
    fn up() -> Self;
    fn rand() -> Self;
}


#[derive(Copy, Clone)]
pub struct Spin {
    x: f64,
    y: f64,
    z: f64,
}


impl Spin {

    pub fn norm(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    pub fn normalized(&self) -> Spin {
        let norm = self.norm();
        Spin {
            x: self.x / norm,
            y: self.y / norm,
            z: self.z / norm,
        }
    }

    pub fn with_norm(&self, norm: f64) -> Spin {
        let temp = self.normalized();
        Spin {
            x: temp.x / norm,
            y: temp.y / norm,
            z: temp.z / norm,
        }
    }

    pub fn dot(&self, other: &Spin) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    pub fn z(&self) -> f64 {
        self.z
    }

}


impl SpinConstructors for Spin {
    fn up() -> Spin {
        Spin { x: 0.0f64, y: 0.0f64, z: 1.0f64,  }
    }

    fn rand() -> Spin {
        let StandardNormal(x) = rand::random();
        let StandardNormal(y) = rand::random();
        let StandardNormal(z) = rand::random();
        let norm = (x * x + y * y + z * z).sqrt();
        Spin { x: x / norm, y: y / norm, z: z / norm, }
    }

}


impl Mul for Spin {

    type Output = f64;

    fn mul(self, other: Spin) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

}


pub type State = Vec<Spin>;


impl StateConstructors for State {

    fn up(size: usize) -> State {
        vec![Spin::up(); size]
    }

    fn rand(size: usize) -> State {
        (0..size).map(|_| { Spin::rand() }).collect()
    }

    fn rand_with_norms(size: usize, norms: &Vec<f64>) -> State {
        (0..size).map(|i| { Spin::rand().with_norm(norms[i]) }).collect::<State>()
    }

}

pub trait CommonObservables {
    fn mag(&self) -> (f64, f64, f64);
    fn mag_len(&self) -> f64 {
        let (x, y, z) = self.mag();
        (x * x + y * y + z * z).sqrt()
    }
}

impl CommonObservables for State {
    fn mag(&self) -> (f64, f64, f64) {
        let (mut x, mut y, mut z) = (0.0f64, 0.0f64, 0.0f64);
        for item in self.iter() {
            x += item.x;
            y += item.y;
            z += item.z;
        }
        (x, y, z)
    }

}
