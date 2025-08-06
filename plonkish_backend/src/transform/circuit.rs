use halo2_curves::bn256::Fr;
use halo2_proofs::plonk::ConstraintSystem;
use zkwasm_halo2::{
    arithmetic::MultiMillerLoop,
    plonk::{
        convert_constraint_system_fr, from_scalar, Circuit as ZkCircuit,
        ConstraintSystem as ZkConstraintSystem,
    },
};

use crate::backend::WitnessEncoding;

#[derive(Debug)]
pub struct ZKWASMCircuit<'a, E: MultiMillerLoop, C: ZkCircuit<E::Scalar>> {
    pub circuit: &'a C,
    pub config: C::Config,
    pub cs: ConstraintSystem<Fr>,
    pub k: u32,
    pub instances: Vec<Vec<Fr>>,
    pub instances_scalar: Vec<Vec<E::Scalar>>,
    pub row_mapping: Vec<usize>,
}

pub fn get_zkwasm_circuit<D: WitnessEncoding, E: MultiMillerLoop, T>(
    k: u32,
    circuit: &'_ [T],
    instances: Vec<E::Scalar>,
) -> ZKWASMCircuit<'_, E, T>
where
    T: ZkCircuit<E::Scalar>,
{
    let circuit = &circuit[0];
    let mut cs = ZkConstraintSystem::default();
    let config = T::configure(&mut cs);

    let cs: ConstraintSystem<Fr> = convert_constraint_system_fr::<E>(cs);
    let row_mapping: Vec<usize> = (0..(1 << k)).collect();
    // Convert Gate Constraints.
    ZKWASMCircuit {
        circuit,
        config,
        cs,
        k,
        instances: vec![instances
            .clone()
            .into_iter()
            .map(|scalar| from_scalar::<E>(&scalar))
            .collect::<Vec<Fr>>()],
        instances_scalar: vec![instances],
        // row_mapping: D::row_mapping(k as usize),
        row_mapping,
    }
}
