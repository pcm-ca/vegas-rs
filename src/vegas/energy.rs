//! Exposes types and traits to create energy components for Monte Carlo simulations
//! they can be agnostic as well and can be used in deterministic simulators as well


use lattice::Adjacency;
use state::State;


pub trait EnergyComponent {
    fn energy(&self, state: &State, index: usize) -> f64;
    fn total_energy(&self, state: &State) -> f64 {
        let mut total = 0f64;
        for i in 0..state.len() {
            total += self.energy(state, i);
        }
        total
    }
}


pub struct ExchangeComponent {
    adjacency: Adjacency,
    exchange: f64,
}


impl ExchangeComponent {
    pub fn new(adjacency: Adjacency, exchange: f64) -> ExchangeComponent {
        ExchangeComponent { adjacency: adjacency, exchange: exchange, }
    }
}


impl EnergyComponent for ExchangeComponent {

    fn energy(&self, state: &State, index: usize) -> f64 {
        let mut ene = 0f64;
        let s = state[index];
        for nb in self.adjacency.nbhs_of(index).unwrap() {
            ene -= s * state[*nb];
        }
        ene * self.exchange
    }

    fn total_energy(&self, state: &State) -> f64 {
        let mut total = 0f64;
        for i in 0..state.len() {
            total += self.energy(state, i);
        }
        0.5 * total
    }
}


pub struct ComplexExchangeComponent {
    adjacency: Adjacency,
}


impl ComplexExchangeComponent {
    pub fn new(adjacency: Adjacency) -> ComplexExchangeComponent {
        ComplexExchangeComponent { adjacency: adjacency }
    }
}


impl EnergyComponent for ComplexExchangeComponent {

    fn energy(&self, state: &State, index: usize) -> f64 {
        let mut ene = 0f64;
        let s = state[index];
        let nbhs = self.adjacency.nbhs_of(index).unwrap().iter();
        let exch = self.adjacency.exch_of(index).unwrap().iter();
        for (nb, exc) in nbhs.zip(exch) {
            ene -= exc * (s * state[*nb]);
        }
        ene
    }

    fn total_energy(&self, state: &State) -> f64 {
        let mut total = 0f64;
        for i in 0..state.len() {
            total += self.energy(state, i);
        }
        0.5 * total
    }
}


pub struct ComposedEnergy<T1, T2> where T1: EnergyComponent, T2: EnergyComponent {
    comp1: T1,
    comp2: T2,
}


impl<T1, T2> ComposedEnergy<T1, T2> where T1: EnergyComponent, T2: EnergyComponent {
    pub fn new(a: T1, b: T2) -> ComposedEnergy<T1, T2> {
        ComposedEnergy { comp1: a, comp2: b }
    }
}


impl<T1, T2> EnergyComponent for ComposedEnergy<T1, T2> where T1: EnergyComponent, T2: EnergyComponent {
    fn energy(&self, state: &State, index: usize) -> f64 {
        self.comp1.energy(state, index) + self.comp2.energy(state, index)
    }

    fn total_energy(&self, state: &State) -> f64 {
        self.comp1.total_energy(state) + self.comp2.total_energy(state)
    }
}
