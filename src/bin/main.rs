#[macro_use] extern crate vegas;

use vegas::lattice::Adjacency;
use vegas::lattice::LatticeBuilder;
use vegas::lattice::Vertex;
use vegas::state::State;
use vegas::state::StateConstructors;
use vegas::state::CommonObservables;
use vegas::energy::EnergyComponent;
use vegas::energy::ComplexExchangeComponent;
use vegas::energy::ZAxisAnisotropyComponent;
use vegas::integrator::Integrator;
use vegas::integrator::MetropolisIntegrator;


pub fn main() {

    let latt = LatticeBuilder::new()
        .pbc((true, true, true))
        .shape((5, 5, 5))
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
            0  => 2.5,
            1  => 2.5,
            2  => 2.0,
            3  => 2.5,
            4  => 2.5,
            5  => 2.0,
            6  => 2.5,
            7  => 2.5,
            8  => 2.5,
            9  => 2.5,
            10 => 2.5,
            11 => 2.5,
            12 => 2.5,
            13 => 2.5,
            14 => 2.0,
            15 => 2.0,
            16 => 2.5,
            17 => 2.0,
            18 => 2.5,
            19 => 2.5,
            20 => 2.5,
            21 => 2.0,
            22 => 2.0,
            23 => 2.0,
            _ => panic!("oh margoth"),
        }
    }).collect();

//    let _cobalt_spins: Vec<f64> = latt.sites().map(|site| {
//        match site.atom() {
//            0  => 1.0,
//            1  => 1.0,
//            _ => panic!("oh margoth"),
//        }
//    }).collect();


    let hamiltonian = hamiltonian!(
        ComplexExchangeComponent::new(Adjacency::new(&latt))
        // ZAxisAnisotropyComponent::new(1.0)
        );

    let mut integrator = MetropolisIntegrator::new(100.0);

    let mut state = State::rand_with_norms(latt.nsites(), &spins);

    println!("# number of sites: {}", state.len());

    loop {
        for _ in 0..1000 {
            state = integrator.step(&hamiltonian, &state);
            println!("{} {}", hamiltonian.total_energy(&state), state.mag_len());
        }
        if integrator.temp() < 10.0 {
            break
        }
        integrator.cool(2.5);
        println!("");
    }

}
