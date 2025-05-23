// Add these functions to plonkish_backend/src/backend/hyperplonk/util.rs

use crate::{
    backend::{PlonkishCircuit, PlonkishCircuitInfo},
    util::{
        arithmetic::{Field, PrimeField},
        expression::{Expression, Query, Rotation},
    },
};
use std::array;

/// Create a PlonkishCircuitInfo for Anemoi hash with TurboPlonk constraints
/// Based on the paper "An efficient verifiable state for zk-EVM and beyond from the Anemoi hash function"
pub fn anemoi_hash_circuit_info<F: PrimeField>(
    num_vars: usize,
    num_rounds: usize,
) -> PlonkishCircuitInfo<F> {
    // We need 5 wires as described in the paper: w1, w2, w3, w4, wo
    // And selector polynomials for the Anemoi constraints
    
    // Selector polynomials:
    // q1, q2, q3, q4: linear combination selectors
    // qm1, qm2: multiplication selectors
    // qc: constant selector
    // qo: output selector
    // qprk1, qprk2, qprk3, qprk4: processed round key selectors
    let num_selectors = 14;
    
    // Create expressions for wire polynomials
    let [w1, w2, w3, w4, wo] = &array::from_fn(|i| {
        Query::new(i, Rotation::cur())
    }).map(Expression::Polynomial);
    
    // Create expressions for selector polynomials
    let [q1, q2, q3, q4, qo, qm1, qm2, qc, qecc, qb, qprk1, qprk2, qprk3, qprk4] = 
        &array::from_fn(|i| Query::new(i + 5, Rotation::cur())).map(Expression::Polynomial);
    
    // Create expressions for next rotation (for Anemoi constraints)
    let [w1_next, w2_next, w3_next, wo_next] = &[
        Query::new(0, Rotation::next()),
        Query::new(1, Rotation::next()),
        Query::new(2, Rotation::next()),
        Query::new(4, Rotation::next()),
    ].map(Expression::Polynomial);
    
    // Base constraint (from vanilla plonk)
    let base_constraint = q1 * w1 + q2 * w2 + q3 * w3 + q4 * w4
        + qm1 * w1 * w2 + qm2 * w3 * w4 + qc
        - qo * wo;
    
    // Boolean constraints
    let bool_constraints = vec![
        qb * w2 * (w2 - Expression::one()),
        qb * w3 * (w3 - Expression::one()),
        qb * w4 * (w4 - Expression::one()),
    ];
    
    // Anemoi constraints (from Section 5.3 of the paper)
    // These implement the closed butterfly S-box verification
    let g = Expression::<F>::Constant(F::GENERATOR);
    let g_inv = Expression::<F>::Constant(F::GENERATOR.invert().unwrap());
    
    // Helper expressions for Anemoi round
    let c_prime_1 = w1 + w4 + g.clone() * (w2 + w3) + qprk3;
    let c_prime_2 = g.clone() * (w1 + w4) + (g.clone() * g.clone() + Expression::one()) * (w2 + w3) + qprk4;
    
    // First Anemoi equation
    let anemoi_1 = qprk3.clone() * (
        (c_prime_1.clone() - w3_next).pow(5)
        + g.clone() * c_prime_1.clone().pow(2)
        - (Expression::Constant(F::from(2u64)) * w1 + w4 + g.clone() * (Expression::Constant(F::from(2u64)) * w2 + w3) + qprk1)
    );
    
    // Second Anemoi equation  
    let anemoi_2 = qprk3.clone() * (
        (c_prime_2.clone() - wo).pow(5)
        + g.clone() * c_prime_2.clone().pow(2)
        - (g.clone() * (Expression::Constant(F::from(2u64)) * w1 + w4) 
           + (g.clone() * g.clone() + Expression::one()) * (Expression::Constant(F::from(2u64)) * w2 + w3) + qprk2)
    );
    
    // Third Anemoi equation
    let anemoi_3 = qprk3.clone() * (
        (c_prime_1 - w3_next.clone()).pow(5)
        + g.clone() * w3_next.pow(2)
        + g_inv.clone()
        - w1_next
    );
    
    // Fourth Anemoi equation
    let anemoi_4 = qprk3.clone() * (
        (c_prime_2 - wo.clone()).pow(5)
        + g * wo.pow(2)
        + g_inv
        - w2_next
    );
    
    // Collect all constraints
    let mut constraints = vec![base_constraint];
    constraints.extend(bool_constraints);
    constraints.extend(vec![anemoi_1, anemoi_2, anemoi_3, anemoi_4]);
    
    // Create preprocessed polynomials (selectors)
    let preprocess_polys = vec![vec![F::ZERO; 1 << num_vars]; num_selectors];
    
    PlonkishCircuitInfo {
        k: num_vars,
        num_instances: vec![0],
        preprocess_polys,
        num_witness_polys: vec![5], // w1, w2, w3, w4, wo
        num_challenges: vec![0],
        constraints,
        lookups: Vec::new(),
        permutations: Vec::new(), // Will be filled based on circuit structure
        max_degree: Some(7), // Due to pow(5) in Anemoi constraints
    }
}

/// Generate a random Anemoi hash circuit for testing
pub fn rand_anemoi_hash_circuit<F: PrimeField, R: Rotatable + From<usize>>(
    num_vars: usize,
    num_rounds: usize,
    mut preprocess_rng: impl RngCore,
    mut witness_rng: impl RngCore,
) -> (PlonkishCircuitInfo<F>, impl PlonkishCircuit<F>) {
    let size = 1 << num_vars;
    let gates_per_round = 16; // From the paper: 16 gates per Jive compression
    
    // Initialize polynomials
    let mut w1_values = vec![F::ZERO; size];
    let mut w2_values = vec![F::ZERO; size];
    let mut w3_values = vec![F::ZERO; size];
    let mut w4_values = vec![F::ZERO; size];
    let mut wo_values = vec![F::ZERO; size];
    
    // Selector polynomials
    let mut q1 = vec![F::ZERO; size];
    let mut q2 = vec![F::ZERO; size];
    let mut q3 = vec![F::ZERO; size];
    let mut q4 = vec![F::ZERO; size];
    let mut qo = vec![F::ZERO; size];
    let mut qm1 = vec![F::ZERO; size];
    let mut qm2 = vec![F::ZERO; size];
    let mut qc = vec![F::ZERO; size];
    let mut qecc = vec![F::ZERO; size];
    let mut qb = vec![F::ZERO; size];
    let mut qprk1 = vec![F::ZERO; size];
    let mut qprk2 = vec![F::ZERO; size];
    let mut qprk3 = vec![F::ZERO; size];
    let mut qprk4 = vec![F::ZERO; size];
    
    // Generate round constants (from π digits as mentioned in the paper)
    let round_constants = generate_anemoi_round_constants::<F>(num_rounds);
    
    let mut permutation = Permutation::default();
    
    // Set up gates for each round
    for round in 0..num_rounds {
        let gate_offset = round * gates_per_round;
        
        // First gate of the round: input setup
        if round == 0 {
            // Initial input values
            w1_values[gate_offset] = F::random(&mut witness_rng);
            w2_values[gate_offset] = F::random(&mut witness_rng);
            w3_values[gate_offset] = F::random(&mut witness_rng);
            w4_values[gate_offset] = F::ZERO; // Salt/padding constant
        }
        
        // Set up Anemoi round constraints
        qprk1[gate_offset] = round_constants[round].0;
        qprk2[gate_offset] = round_constants[round].1;
        qprk3[gate_offset] = round_constants[round].2;
        qprk4[gate_offset] = round_constants[round].3;
        
        // Connect output of this round to input of next round via copy constraints
        if round < num_rounds - 1 {
            let next_gate = (round + 1) * gates_per_round;
            permutation.copy((0, gate_offset + 1), (0, next_gate)); // w1
            permutation.copy((1, gate_offset + 1), (1, next_gate)); // w2
            permutation.copy((2, gate_offset + 1), (2, next_gate)); // w3
            permutation.copy((4, gate_offset), (3, next_gate));     // wo -> w4
        }
    }
    
    // Create circuit info
    let circuit_info = PlonkishCircuitInfo {
        k: num_vars,
        num_instances: vec![0],
        preprocess_polys: vec![
            q1, q2, q3, q4, qo, qm1, qm2, qc, qecc, qb,
            qprk1, qprk2, qprk3, qprk4
        ],
        num_witness_polys: vec![5],
        num_challenges: vec![0],
        constraints: anemoi_hash_circuit_info::<F>(num_vars, num_rounds).constraints,
        lookups: Vec::new(),
        permutations: permutation.into_cycles(),
        max_degree: Some(7),
    };
    
    let witness = vec![w1_values, w2_values, w3_values, w4_values, wo_values];
    
    (
        circuit_info,
        MockCircuit::new(vec![], witness),
    )
}

/// Generate Anemoi round constants from π digits
fn generate_anemoi_round_constants<F: PrimeField>(num_rounds: usize) -> Vec<(F, F, F, F)> {
    // In practice, these would be derived from π digits as described in the paper
    // This is a placeholder implementation
    (0..num_rounds)
        .map(|i| {
            let base = F::from((i * 4) as u64);
            (
                base,
                base + F::ONE,
                base + F::from(2u64),
                base + F::from(3u64),
            )
        })
        .collect()
}

/// Create a circuit that computes Anemoi Jive compression (3-to-1)
pub fn anemoi_jive_compression_circuit<F: PrimeField>(
    num_vars: usize,
) -> PlonkishCircuitInfo<F> {
    const NUM_ROUNDS: usize = 14; // From the paper
    
    // Create base Anemoi circuit
    let mut circuit_info = anemoi_hash_circuit_info::<F>(num_vars, NUM_ROUNDS);
    
    // Add constraint for final sum computation (Jive mode)
    // The output is sum of input and output of permutation
    let output_sum_gate = NUM_ROUNDS * 16; // After all round gates
    
    // Additional constraints for Jive compression can be added here
    
    circuit_info
}

/// Helper trait extension for Expression to add pow method
trait ExpressionExt<F: Field> {
    fn pow(self, n: u32) -> Expression<F>;
}

impl<F: Field> ExpressionExt<F> for Expression<F> {
    fn pow(self, n: u32) -> Expression<F> {
        match n {
            0 => Expression::one(),
            1 => self,
            _ => {
                let half = self.clone().pow(n / 2);
                let half_squared = half.clone() * half;
                if n % 2 == 0 {
                    half_squared
                } else {
                    half_squared * self
                }
            }
        }
    }
}