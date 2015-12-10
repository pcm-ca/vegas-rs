extern crate rand;
use rand::distributions::{IndependentSample, Range};

use state::Spin;
use state::State;
use state::SpinConstructors;
use state::StateConstructors;
use energy::EnergyComponent;


pub trait Integrator<T: EnergyComponent> {
    fn step(&self, energy: &T, state: &State) -> State;
}


pub struct MetropolisIntegrator {
    temp: f64,
}


impl MetropolisIntegrator {
    pub fn new(temp: f64) -> MetropolisIntegrator {
        MetropolisIntegrator { temp: temp }
    }

    pub fn temp(&self) -> f64 {
        self.temp
    }

    pub fn heat(&mut self, delta: f64) {
        self.temp += delta;
    }

    pub fn cool(&mut self, delta: f64) {
        self.heat( - delta);
    }
}


impl<T: EnergyComponent> Integrator<T> for MetropolisIntegrator {
    fn step(&self, energy: &T, state: &State) -> State {
        let mut new_state = (*state).clone();
        let sites = Range::new(0, new_state.len());
        let mut rng = rand::thread_rng();
        for _ in 0..new_state.len() {
            let site = sites.ind_sample(&mut rng);
            let old_energy = energy.energy(&new_state, site);
            new_state[site] = Spin::rand().with_norm(state[site].norm());
            let new_energy = energy.energy(&new_state, site);
            let delta = new_energy - old_energy;
            if delta < 0.0 {
                continue
            }
            if rand::random::<f64>() < (- delta / self.temp).exp() {
                continue
            }
            new_state[site] = state[site];
        }
        new_state
    }
}
