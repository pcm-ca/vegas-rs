#[macro_use] extern crate vegas;

use vegas::lattice::Adjacency;
use vegas::lattice::LatticeBuilder;
use vegas::lattice::Vertex;
use vegas::state::State;
use vegas::state::StateConstructors;
use vegas::state::CommonObservables;
use vegas::energy::EnergyComponent;
use vegas::energy::ComplexExchangeComponent;
use vegas::energy::ZAxisAnisotropy;
use vegas::integrator::Integrator;
use vegas::integrator::MetropolisIntegrator;


pub fn main() {

    let latt = LatticeBuilder::new()
        .pbc((true, true, true))
        .shape((4, 4, 4))
        .vertices(Vertex::list_for_magnetite())
        .natoms(24)
        .finalize();

    let _manganite_spins: Vec<f64> = latt.sites().map(|site| {
        match site.atom() {
            0  => 2.0,
            1  => 1.5,
            2  => 2.0,
            3  => 1.5,
            4  => 2.0,
            5  => 2.0,
            6  => 2.0,
            7  => 2.0,
            8  => 1.5,
            9  => 1.5,
            10 => 2.0,
            11 => 2.0,
            12 => 2.0,
            13 => 2.0,
            14 => 1.5,
            15 => 2.0,
            16 => 1.5,
            17 => 2.0,
            18 => 2.0,
            19 => 2.0,
            20 => 1.5,
            21 => 2.0,
            22 => 1.5,
            23 => 2.0,
            24 => 1.5,
            25 => 2.0,
            26 => 2.0,
            _ => panic!("oh margoth"),
        }
    }).collect();

    let spins: Vec<f64> = latt.sites().map(|site| {
        match site.atom() {
            0  => 5.0,
            1  => 5.0,
            2  => 6.0,
            3  => 6.0,
            4  => 6.0,
            5  => 6.0,
            6  => 5.0,
            7  => 5.0,
            8  => 6.0,
            9  => 6.0,
            10 => 6.0,
            11 => 6.0,
            12 => 5.0,
            13 => 5.0,
            14 => 6.0,
            15 => 6.0,
            16 => 6.0,
            17 => 6.0,
            18 => 5.0,
            19 => 5.0,
            20 => 6.0,
            21 => 6.0,
            22 => 6.0,
            23 => 6.0,
            _ => panic!("oh margoth"),
        }
    }).collect();

    let mut state = State::rand_with_norms(latt.nsites(), &spins);

    let hamiltonian = ComplexExchangeComponent::new(Adjacency::new(&latt));

    let mut integrator = MetropolisIntegrator::new(250.0);

    loop {
        for _ in 0..1000 {
            state = integrator.step(&hamiltonian, &state);
            println!("{} {}", hamiltonian.total_energy(&state), state.mag_len());
        }
        if integrator.temp() < 0.1 {
            break
        }
        integrator.cool(10.0);
        println!("");
    }

}
