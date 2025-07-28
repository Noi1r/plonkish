use crate::{
    backend::{
        hyperplonk::util::common::*,
        mock::MockCircuit,
        PlonkishCircuit, PlonkishCircuitInfo,
    },
    util::{
        arithmetic::PrimeField,
        expression::{rotate::Rotatable, Expression, Query, Rotation},
    },
};

use rand::RngCore;
use std::array;

// Import the ExpressionExt trait for pow method
use super::common::ExpressionExt;

/// Simple index test circuit: verifies w_next = w_current + 1
pub fn rand_index_test_circuit<F: PrimeField, R: Rotatable + From<usize>>(
    num_vars: usize,
    _preprocess_rng: impl RngCore,
    mut witness_rng: impl RngCore,
) -> (PlonkishCircuitInfo<F>, impl PlonkishCircuit<F>) {
    let size = 1 << num_vars;
    let usable_indices = R::from(num_vars).usable_indices();

    // Initialize polynomials
    let mut q_seq = vec![F::ZERO; size]; // Selector polynomial
    let mut w_values = vec![F::ZERO; size]; // Witness polynomial

    let mut permutation = Permutation::default();
    permutation.copy((1, usable_indices[0]), (1, usable_indices[0])); // Self-reference

    // Set up incrementing sequence on usable_indices
    for i in 0..std::cmp::min(usable_indices.len().saturating_sub(1), 8) {
        let current_idx = usable_indices[i];
        let next_idx = usable_indices[i + 1];

        // Set values: 0, 1, 2, 3, ...
        w_values[current_idx] = F::from(i as u64);
        w_values[next_idx] = F::from((i + 1) as u64);

        // Activate constraint
        if i != 0 {
            q_seq[current_idx] = F::ONE;
        }
    }

    // Constraint: q_seq * (w_next - w_current - 1) = 0
    let constraint: Expression<F> = Expression::Polynomial::<F>(Query::new(0, Rotation::cur()))
        * (Expression::Polynomial::<F>(Query::new(1, Rotation::next()))
            - Expression::Polynomial::<F>(Query::new(1, Rotation::prev()))
            - Expression::one()
            - Expression::one());

    let circuit_info = PlonkishCircuitInfo {
        k: num_vars,
        num_instances: vec![],
        preprocess_polys: vec![q_seq],
        num_witness_polys: vec![1],
        num_challenges: vec![0],
        constraints: vec![constraint],
        lookups: Vec::new(),
        permutations: permutation.into_cycles(),
        max_degree: Some(2),
    };

    (circuit_info, MockCircuit::new(vec![], vec![w_values]))
}

/// Build a complete Jive CRH circuit for Merkle tree level
pub fn build_jive_crh_circuit<F: PrimeField>(
    left_hash: F,
    middle_hash: F,
    right_hash: F,
    level_padding: F,
    constants: &AnemoiConstants<F>,
    round_keys: &AnemoiRoundKeys<F>,
    usable_indices: &[usize],
    gate_counter: &mut usize,
    w1_values: &mut [F],
    w2_values: &mut [F],
    w3_values: &mut [F],
    w4_values: &mut [F],
    wo_values: &mut [F],
    qprk1: &mut [F],
    qprk2: &mut [F],
    qprk3: &mut [F],
    qprk4: &mut [F],
    q1: &mut [F],
    q2: &mut [F],
    q3: &mut [F],
    q4: &mut [F],
    qo: &mut [F],
    permutation: &mut Permutation,
) -> F {
    const NUM_ROUNDS: usize = 14;

    // Input gate for Jive CRH
    let jive_input_idx = usable_indices[*gate_counter];
    *gate_counter += 1;

    let mut jive_x = [left_hash, middle_hash];
    let mut jive_y = [right_hash, level_padding];
    let sum_before_perm = jive_x[0] + jive_x[1] + jive_y[0] + jive_y[1];

    w1_values[jive_input_idx] = jive_x[0];
    w2_values[jive_input_idx] = jive_x[1];
    w3_values[jive_input_idx] = jive_y[0];
    w4_values[jive_input_idx] = jive_y[1];
    wo_values[jive_input_idx] = sum_before_perm;

    // Enforce constraints for input gate
    q1[jive_input_idx] = F::ONE;
    q2[jive_input_idx] = F::ONE;
    q3[jive_input_idx] = F::ONE;
    q4[jive_input_idx] = F::ONE;
    qo[jive_input_idx] = F::ONE;

    // Create self-referencing permutations for this Jive instance
    for poly_idx in 14..19 {
        permutation.copy((poly_idx, jive_input_idx), (poly_idx, jive_input_idx));
    }

    // 14 Anemoi rounds
    for round in 0..NUM_ROUNDS {
        let round_gate_idx = usable_indices[*gate_counter];
        *gate_counter += 1;

        // Set preprocessed round keys
        qprk1[round_gate_idx] = round_keys.preprocessed_round_keys_x[round][0];
        qprk2[round_gate_idx] = round_keys.preprocessed_round_keys_x[round][1];
        qprk3[round_gate_idx] = round_keys.preprocessed_round_keys_y[round][0];
        qprk4[round_gate_idx] = round_keys.preprocessed_round_keys_y[round][1];

        // Store pre-round state
        w1_values[round_gate_idx] = jive_x[0];
        w2_values[round_gate_idx] = jive_x[1];
        w3_values[round_gate_idx] = jive_y[0];
        w4_values[round_gate_idx] = jive_y[1];

        // Apply Anemoi round transformation using helper function
        apply_anemoi_round(&mut jive_x, &mut jive_y, round, constants, round_keys);
    }

    // Next round values gate for rotation::next() constraints
    let next_round_idx = usable_indices[*gate_counter];
    *gate_counter += 1;

    w1_values[next_round_idx] = jive_x[0];
    w2_values[next_round_idx] = jive_x[1];
    w3_values[next_round_idx] = jive_y[0];
    w4_values[next_round_idx] = jive_y[1];

    // Final state constraint
    q1[next_round_idx] = F::from(3u64) * constants.generator + F::from(3u64);
    q2[next_round_idx] = F::from(3u64) * (constants.generator * constants.generator + constants.generator + F::ONE);
    q3[next_round_idx] = F::from(2u64) * (constants.generator * constants.generator + constants.generator + F::ONE);
    q4[next_round_idx] = F::from(2u64) * constants.generator + F::from(2u64);
    qo[next_round_idx] = F::ONE;

    // Apply final transformations using helper function
    apply_final_transformations(&mut jive_x, &mut jive_y, constants);

    let sum_after_perm = jive_x[0] + jive_x[1] + jive_y[0] + jive_y[1];
    wo_values[next_round_idx] = sum_after_perm;

    // Final Jive output gate
    let jive_output_idx = usable_indices[*gate_counter];
    *gate_counter += 1;

    let current_hash = sum_before_perm + sum_after_perm; // Jive output

    w1_values[jive_output_idx] = sum_before_perm;
    w2_values[jive_output_idx] = sum_after_perm;
    w3_values[jive_output_idx] = F::ZERO;
    w4_values[jive_output_idx] = F::ZERO;
    wo_values[jive_output_idx] = current_hash;

    // Constraint: wo = w1 + w2
    q1[jive_output_idx] = F::ONE;
    q2[jive_output_idx] = F::ONE;
    qo[jive_output_idx] = F::ONE;

    permutation.copy((18, next_round_idx), (15, jive_output_idx));
    permutation.copy((14, jive_output_idx), (18, usable_indices[*gate_counter - 17]));

    current_hash
}

pub fn rand_jive_crh_circuit<F: PrimeField, R: Rotatable + From<usize>>(
    input1: F,
    input2: F,
    input3: F,
    padding_constant: F, // For domain separation
    _preprocess_rng: impl RngCore,
    _witness_rng: impl RngCore,
) -> (PlonkishCircuitInfo<F>, impl PlonkishCircuit<F>) {
    const NUM_ROUNDS: usize = 14;
    const TOTAL_GATES: usize = NUM_ROUNDS + 3;

    let num_vars = (TOTAL_GATES as f64).log2().ceil() as usize;
    let size = 1 << num_vars;
    let usable_indices = R::from(num_vars).usable_indices();

    if usable_indices.len() < TOTAL_GATES {
        panic!(
            "Not enough usable indices: {} < {}",
            usable_indices.len(),
            TOTAL_GATES
        );
    }

    // Initialize polynomials using helper functions
    let (mut w1_values, mut w2_values, mut w3_values, mut w4_values, mut wo_values) = 
        init_witness_polynomials(size);
    let [mut q1, mut q2, mut q3, mut q4, mut qo, mut qm1, mut qm2, mut qc, mut qecc, mut qb, mut qprk1, mut qprk2, mut qprk3, mut qprk4] = 
        init_selector_polynomials(size);

    let mut permutation = Permutation::default();

    // Use helper structs for constants and keys
    let constants = AnemoiConstants::new();
    let round_keys = AnemoiRoundKeys::new();

    // Set initial state: x[1]=input1, x[2]=input2, y[1]=input3, y[2]=padding_constant
    let mut current_x = [input1, input2];
    let mut current_y = [input3, padding_constant];
    let sum_before_perm = current_x[0] + current_x[1] + current_y[0] + current_y[1];

    // Input gate (Gate 0) - enforce salt constraint
    let input_idx = usable_indices[0];
    w1_values[input_idx] = current_x[0];
    w2_values[input_idx] = current_x[1];
    w3_values[input_idx] = current_y[0];
    w4_values[input_idx] = current_y[1];
    wo_values[input_idx] = sum_before_perm; // Not used in input gate

    // Enforce w4 = padding_constant
    q1[input_idx] = F::ONE;
    q2[input_idx] = F::ONE;
    q3[input_idx] = F::ONE;
    q4[input_idx] = F::ONE;
    qo[input_idx] = F::ONE;
    // qc[input_idx] = -padding_constant;

    for poly_idx in 14..19 {
        // w1, w2, w3, w4, wo
        permutation.copy((poly_idx, input_idx), (poly_idx, input_idx));
    }

    // Anemoi rounds (Gates 1-14)
    for round in 0..NUM_ROUNDS {
        let gate_idx = usable_indices[round + 1];

        // Set preprocessed round keys
        qprk1[gate_idx] = round_keys.preprocessed_round_keys_x[round][0];
        qprk2[gate_idx] = round_keys.preprocessed_round_keys_x[round][1];
        qprk3[gate_idx] = round_keys.preprocessed_round_keys_y[round][0];
        qprk4[gate_idx] = round_keys.preprocessed_round_keys_y[round][1];

        // Store pre-round state
        w1_values[gate_idx] = current_x[0];
        w2_values[gate_idx] = current_x[1];
        w3_values[gate_idx] = current_y[0];
        w4_values[gate_idx] = current_y[1];

        // Apply Anemoi round using helper function
        apply_anemoi_round(&mut current_x, &mut current_y, round, &constants, &round_keys);
    }
    // Set up next round values for the last round to enable Rotation::next() constraints
    let final_idx = usable_indices[NUM_ROUNDS + 1];
    w1_values[final_idx] = current_x[0];
    w2_values[final_idx] = current_x[1];
    w3_values[final_idx] = current_y[0];
    w4_values[final_idx] = current_y[1];

    // final state constraint
    // TODO: For better security, this should be split into multiple constraints.
    // Currently implemented as PoC with combined constraint
    q1[final_idx] = F::from(3u64) * constants.generator + F::from(3u64);
    q2[final_idx] = F::from(3u64) * (constants.generator * constants.generator + constants.generator + F::ONE);
    q3[final_idx] = F::from(2u64) * (constants.generator * constants.generator + constants.generator + F::ONE);
    q4[final_idx] = F::from(2u64) * constants.generator + F::from(2u64);
    qo[final_idx] = F::ONE;

    // Apply final transformations using helper function
    apply_final_transformations(&mut current_x, &mut current_y, &constants);

    let sum_after_perm = current_x[0] + current_x[1] + current_y[0] + current_y[1];
    wo_values[final_idx] = sum_after_perm;

    // Final output gate (Gate NUM_ROUNDS+2) - Jive output
    let output_gate_idx = usable_indices[NUM_ROUNDS + 2];
    let jive_output = sum_before_perm + sum_after_perm;

    w1_values[output_gate_idx] = sum_before_perm;
    w2_values[output_gate_idx] = sum_after_perm;
    w3_values[output_gate_idx] = F::ZERO;
    w4_values[output_gate_idx] = F::ZERO;
    wo_values[output_gate_idx] = jive_output;
    permutation.copy((18, final_idx), (15, output_gate_idx));
    permutation.copy(
        (14, output_gate_idx),
        (18, usable_indices[NUM_ROUNDS + 2 - 16]),
    );

    // Constraint: wo = w1 + w2 (final Jive output)
    q1[output_gate_idx] = F::ONE;
    q2[output_gate_idx] = F::ONE;
    qo[output_gate_idx] = F::ONE;

    // Set up copy constraints for round transitions
    for round in 1..=NUM_ROUNDS {
        let current_idx = usable_indices[round];
        let next_idx = usable_indices[round + 1];

        wo_values[current_idx] = w4_values[next_idx];
        permutation.copy((18, current_idx), (17, next_idx)); // wo -> w4
    }

    let circuit_info = PlonkishCircuitInfo {
        k: num_vars,
        num_instances: vec![],
        preprocess_polys: vec![
            q1, q2, q3, q4, qo, qm1, qm2, qc, qecc, qb, qprk1, qprk2, qprk3, qprk4,
        ],
        num_witness_polys: vec![5],
        num_challenges: vec![0],
        constraints: anemoi_hash_circuit_info::<F>(num_vars).constraints,
        lookups: Vec::new(),
        permutations: permutation.into_cycles(),
        max_degree: Some(7),
    };

    let witness = vec![w1_values, w2_values, w3_values, w4_values, wo_values];

    (circuit_info, MockCircuit::new(vec![], witness))
}

/// Create a PlonkishCircuitInfo for Anemoi hash with TurboPlonk constraints
/// Based on the paper "An efficient verifiable state for zk-EVM and beyond from the Anemoi hash function"
pub fn anemoi_hash_circuit_info_original<F: PrimeField>(
    num_vars: usize,
) -> PlonkishCircuitInfo<F> {
    //14 selector polynomials for the Anemoi constraints
    // Create expressions for selector polynomials
    let [q1, q2, q3, q4, qo, qm1, qm2, qc, qecc, qb, qprk1, qprk2, qprk3, qprk4] =
        &array::from_fn(|i| Query::new(i, Rotation::cur())).map(Expression::Polynomial);

    // We need 5 wires as described in the paper: w1, w2, w3, w4, wo
    // Create expressions for wire polynomials
    let [w1, w2, w3, w4, wo] =
        &array::from_fn(|i| Query::new(i + 14, Rotation::cur())).map(Expression::Polynomial);

    // Create expressions for next rotation (for Anemoi constraints)
    let [w1_next, w2_next, w3_next, w4_next] = &[
        Query::new(14, Rotation::next()),
        Query::new(15, Rotation::next()),
        Query::new(16, Rotation::next()),
        Query::new(17, Rotation::next()),
    ]
    .map(Expression::Polynomial);

    // Base constraint (TurboPlonk gate)
    let base_constraint = q1 * w1
        + q2 * w2
        + q3 * w3
        + q4 * w4
        + qm1 * w1 * w2
        + qm2 * w3 * w4
        + qc
        + qecc * w1 * w2 * w3 * w4 * wo
        - qo * wo;

    // Boolean constraints
    let bool_constraints = vec![
        qb * w2 * (w2 - Expression::one()),
        qb * w3 * (w3 - Expression::one()),
        qb * w4 * (w4 - Expression::one()),
    ];

    // Anemoi constraints (from Section 6 of the paper)
    let g = Expression::<F>::Constant(F::from(5u64)); // Generator
    let g_inv = Expression::<F>::Constant(
        F::from_str_vartime(
            "8755297148735710088898562298102910035419345760166413737479281674630323398247",
        )
        .unwrap(),
    ); // Delta

    // Helper expressions for Anemoi round
    let c_prime_1 = w1 + w4 + g.clone() * (w2 + w3) + qprk3;
    let c_prime_2 =
        g.clone() * (w1 + w4) + (g.clone() * g.clone() + Expression::one()) * (w2 + w3) + qprk4;

    // First Anemoi equation
    let anemoi_1 = qprk3.clone()
        * ((c_prime_1.clone() - w3_next).pow(5) + g.clone() * c_prime_1.clone().pow(2)
            - (Expression::Constant(F::from(2u64)) * w1
                + w4
                + g.clone() * (Expression::Constant(F::from(2u64)) * w2 + w3)
                + qprk1));

    // Second Anemoi equation
    let anemoi_2 = qprk3.clone()
        * ((c_prime_2.clone() - w4_next).pow(5) + g.clone() * c_prime_2.clone().pow(2)
            - (g.clone() * (Expression::Constant(F::from(2u64)) * w1 + w4)
                + (g.clone() * g.clone() + Expression::one())
                    * (Expression::Constant(F::from(2u64)) * w2 + w3)
                + qprk2));

    // Third Anemoi equation
    let anemoi_3 = qprk3.clone()
        * ((c_prime_1 - w3_next.clone()).pow(5)
            + g.clone() * w3_next.clone().pow(2)
            + g_inv.clone()
            - w1_next);

    // Fourth Anemoi equation
    let anemoi_4 = qprk3.clone()
        * ((c_prime_2 - w4_next.clone()).pow(5)
            + g.clone() * w4_next.clone().pow(2)
            + g_inv.clone()
            - w2_next);

    // Collect all constraints
    let mut constraints = vec![base_constraint];
    constraints.extend(bool_constraints);
    constraints.extend(vec![anemoi_1, anemoi_2, anemoi_3, anemoi_4]);

    // Create preprocessed polynomials (selectors) - 14 selectors
    let preprocess_polys = vec![vec![F::ZERO; 1 << num_vars]; 14];

    PlonkishCircuitInfo {
        k: num_vars,
        num_instances: vec![0],
        preprocess_polys,
        num_witness_polys: vec![5], // w1, w2, w3, w4, wo
        num_challenges: vec![0],
        constraints,
        lookups: Vec::new(),
        permutations: Vec::new(),
        max_degree: Some(7), // Due to pow(5) in Anemoi constraints
    }
}

/// Create a PlonkishCircuitInfo for Anemoi hash with TurboPlonk constraints
/// Based on the paper "An efficient verifiable state for zk-EVM and beyond from the Anemoi hash function"
pub fn anemoi_hash_circuit_info<F: PrimeField>(
    num_vars: usize,
) -> PlonkishCircuitInfo<F> {
    //14 selector polynomials for the Anemoi constraints
    // Create expressions for selector polynomials
    let [q1, q2, q3, q4, qo, qm1, qm2, qc, qecc, qb, qprk1, qprk2, qprk3, qprk4] =
        &array::from_fn(|i| Query::new(i, Rotation::cur())).map(Expression::Polynomial);

    // We need 5 wires as described in the paper: w1, w2, w3, w4, wo
    // Create expressions for wire polynomials
    let [w1, w2, w3, w4, wo] =
        &array::from_fn(|i| Query::new(i + 14, Rotation::cur())).map(Expression::Polynomial);

    // Create expressions for next rotation (for Anemoi constraints)
    let [w1_next, w2_next, w3_next] = &[
        Query::new(14, Rotation::next()),
        Query::new(15, Rotation::next()),
        Query::new(16, Rotation::next()),
    ]
    .map(Expression::Polynomial);

    // Base constraint (TurboPlonk gate)
    let base_constraint = q1 * w1
        + q2 * w2
        + q3 * w3
        + q4 * w4
        + qm1 * w1 * w2
        + qm2 * w3 * w4
        + qc
        + qecc * w1 * w2 * w3 * w4 * wo
        - qo * wo;

    // Boolean constraints
    let bool_constraints = vec![
        qb * w2 * (w2 - Expression::one()),
        qb * w3 * (w3 - Expression::one()),
        qb * w4 * (w4 - Expression::one()),
    ];

    // Anemoi constraints (from Section 6 of the paper)
    let g = Expression::<F>::Constant(F::from(5u64)); // Generator
    let g_inv = Expression::<F>::Constant(
        F::from_str_vartime(
            "8755297148735710088898562298102910035419345760166413737479281674630323398247",
        )
        .unwrap(),
    ); // Delta

    // Helper expressions for Anemoi round
    let c_prime_1 = w1 + w4 + g.clone() * (w2 + w3) + qprk3;
    let c_prime_2 =
        g.clone() * (w1 + w4) + (g.clone() * g.clone() + Expression::one()) * (w2 + w3) + qprk4;

    // First Anemoi equation
    let anemoi_1 = qprk3.clone()
        * ((c_prime_1.clone() - w3_next).pow(5) + g.clone() * c_prime_1.clone().pow(2)
            - (Expression::Constant(F::from(2u64)) * w1
                + w4
                + g.clone() * (Expression::Constant(F::from(2u64)) * w2 + w3)
                + qprk1));

    // Second Anemoi equation
    let anemoi_2 = qprk3.clone()
        * ((c_prime_2.clone() - wo).pow(5) + g.clone() * c_prime_2.clone().pow(2)
            - (g.clone() * (Expression::Constant(F::from(2u64)) * w1 + w4)
                + (g.clone() * g.clone() + Expression::one())
                    * (Expression::Constant(F::from(2u64)) * w2 + w3)
                + qprk2));

    // Third Anemoi equation
    let anemoi_3 = qprk3.clone()
        * ((c_prime_1 - w3_next.clone()).pow(5)
            + g.clone() * w3_next.clone().pow(2)
            + g_inv.clone()
            - w1_next);

    // Fourth Anemoi equation
    let anemoi_4 = qprk3.clone()
        * ((c_prime_2 - wo.clone()).pow(5) + g.clone() * wo.clone().pow(2) + g_inv.clone()
            - w2_next);

    // Collect all constraints
    let mut constraints = vec![base_constraint];
    constraints.extend(bool_constraints);
    constraints.extend(vec![anemoi_1, anemoi_2, anemoi_3, anemoi_4]);

    // Create preprocessed polynomials (selectors) - 14 selectors
    let preprocess_polys = vec![vec![F::ZERO; 1 << num_vars]; 14];

    PlonkishCircuitInfo {
        k: num_vars,
        num_instances: vec![0],
        preprocess_polys,
        num_witness_polys: vec![5], // w1, w2, w3, w4, wo
        num_challenges: vec![0],
        constraints,
        lookups: Vec::new(),
        permutations: Vec::new(),
        max_degree: Some(7), // Due to pow(5) in Anemoi constraints
    }
}

/// Generate a complete Anemoi hash circuit implementing the paper's specifications
/// Automatically determines the required num_vars based on gate count
pub fn rand_anemoi_hash_circuit_with_flatten<F: PrimeField, R: Rotatable + From<usize>>(
    _preprocess_rng: impl RngCore,
    mut witness_rng: impl RngCore,
) -> (PlonkishCircuitInfo<F>, impl PlonkishCircuit<F>) {
    const NUM_ROUNDS: usize = 14; // From AnemoiJive256 specification
    const GATES_NEEDED: usize = NUM_ROUNDS + 4; // Input + rounds + output + padding

    // Calculate minimum num_vars needed
    let num_vars = (GATES_NEEDED as f64).log2().ceil() as usize;
    let size = 1 << num_vars;

    // Get usable indices
    let usable_indices = R::from(num_vars).usable_indices();

    if usable_indices.len() < GATES_NEEDED {
        panic!(
            "Not enough usable indices: {} < {}",
            usable_indices.len(),
            GATES_NEEDED
        );
    }

    // Initialize polynomials
    let (mut w1_values, mut w2_values, mut w3_values, mut w4_values, mut wo_values) = 
        init_witness_polynomials(size);
    let [mut q1, mut q2, mut q3, mut q4, mut qo, mut qm1, mut qm2, mut qc, mut qecc, mut qb, mut qprk1, mut qprk2, mut qprk3, mut qprk4] = 
        init_selector_polynomials(size);

    let mut permutation = Permutation::default();

    // Use real Anemoi constants
    let constants = AnemoiConstants::new();
    let round_keys = AnemoiRoundKeys::new();

    // Create initial state (2 field elements each for x and y)
    let mut current_x = [F::from(1u64), F::from(1u64)];
    let mut current_y = [F::from(3u64), F::from(423u64)];

    // Use usable_indices to set initial state
    let init_idx = usable_indices[0];
    w1_values[init_idx] = current_x[0];
    w2_values[init_idx] = current_x[1];
    w3_values[init_idx] = current_y[0];
    w4_values[init_idx] = current_y[1];
    wo_values[init_idx] = current_x[0] + current_x[1] + current_y[0] + current_y[1];

    // Simple constraint for input gate
    q1[init_idx] = F::ONE;
    q2[init_idx] = F::ONE;
    q3[init_idx] = F::ONE;
    q4[init_idx] = F::ONE;
    qo[init_idx] = F::ONE;

    // Create self-referencing permutations for witness polynomials
    for poly_idx in 14..19 {
        // w1, w2, w3, w4, wo
        permutation.copy((poly_idx, init_idx), (poly_idx, init_idx));
    }

    // Anemoi rounds - simulate the actual Anemoi permutation
    for round in 0..NUM_ROUNDS {
        let current_gate_idx = usable_indices[round + 1];

        // Ensure next index exists (for Rotation::next())
        if round + 2 >= usable_indices.len() {
            panic!("Not enough usable indices for round {}", round);
        }

        // Set preprocessed round key selectors (use preprocessed values for qprk)
        qprk1[current_gate_idx] = round_keys.preprocessed_round_keys_x[round][0];
        qprk2[current_gate_idx] = round_keys.preprocessed_round_keys_x[round][1];
        qprk3[current_gate_idx] = round_keys.preprocessed_round_keys_y[round][0];
        qprk4[current_gate_idx] = round_keys.preprocessed_round_keys_y[round][1];

        // Store pre-round values
        w1_values[current_gate_idx] = current_x[0];
        w2_values[current_gate_idx] = current_x[1];
        w3_values[current_gate_idx] = current_y[0];
        w4_values[current_gate_idx] = current_y[1];

        // Apply one round of Anemoi permutation
        apply_anemoi_round(&mut current_x, &mut current_y, round, &constants, &round_keys);
    }

    // Set up next round values for the last round to enable Rotation::next() constraints
    let final_idx = usable_indices[NUM_ROUNDS + 1];
    w1_values[final_idx] = current_x[0];
    w2_values[final_idx] = current_x[1];
    w3_values[final_idx] = current_y[0];
    w4_values[final_idx] = current_y[1];

    // final state constraint
    q1[final_idx] = F::from(3u64) * constants.generator + F::from(3u64);
    q2[final_idx] = F::from(3u64) * (constants.generator * constants.generator + constants.generator + F::ONE);
    q3[final_idx] = F::from(2u64) * (constants.generator * constants.generator + constants.generator + F::ONE);
    q4[final_idx] = F::from(2u64) * constants.generator + F::from(2u64);
    qo[final_idx] = F::ONE;

    // Apply final transformations using helper function
    apply_final_transformations(&mut current_x, &mut current_y, &constants);

    let sum_after_perm = current_x[0] + current_x[1] + current_y[0] + current_y[1];
    wo_values[final_idx] = sum_after_perm;

    // Final output gate (Gate NUM_ROUNDS+2) - Jive output
    let output_gate_idx = usable_indices[NUM_ROUNDS + 2];
    let jive_output = w1_values[final_idx] + w2_values[final_idx] + w3_values[final_idx] + w4_values[final_idx] + sum_after_perm;

    w1_values[output_gate_idx] = w1_values[final_idx] + w2_values[final_idx] + w3_values[final_idx] + w4_values[final_idx];
    w2_values[output_gate_idx] = sum_after_perm;
    w3_values[output_gate_idx] = F::ZERO;
    w4_values[output_gate_idx] = F::ZERO;
    wo_values[output_gate_idx] = jive_output;
    permutation.copy((18, final_idx), (15, output_gate_idx));
    permutation.copy(
        (14, output_gate_idx),
        (18, usable_indices[NUM_ROUNDS + 2 - 16]),
    );

    // Constraint: wo = w1 + w2 (final Jive output)
    q1[output_gate_idx] = F::ONE;
    q2[output_gate_idx] = F::ONE;
    qo[output_gate_idx] = F::ONE;

    // Set up copy constraints using usable indices
    for round in 1..=NUM_ROUNDS {
        let current_idx = usable_indices[round];
        let next_idx = usable_indices[round + 1];

        // wo -> w4 connection (Section 5.3 optimization)
        wo_values[current_idx] = w4_values[next_idx];
        permutation.copy((18, current_idx), (17, next_idx)); // wo -> w4 (next round)
    }

    // Create circuit info
    let circuit_info = PlonkishCircuitInfo {
        k: num_vars,
        num_instances: vec![], // No instances
        preprocess_polys: vec![
            q1, q2, q3, q4, qo, qm1, qm2, qc, qecc, qb, qprk1, qprk2, qprk3, qprk4,
        ],
        num_witness_polys: vec![5],
        num_challenges: vec![0],
        constraints: anemoi_hash_circuit_info::<F>(num_vars).constraints,
        lookups: Vec::new(),
        permutations: permutation.into_cycles(),
        max_degree: Some(7),
    };

    let witness = vec![w1_values, w2_values, w3_values, w4_values, wo_values];

    (circuit_info, MockCircuit::new(vec![], witness))
}

/// Generic Anemoi Jive hash implementation using PrimeField
/// This extracts the logic from merkle_membership_proof_circuit to make it reusable
pub fn anemoi_jive_hash<F: PrimeField>(x: [F; 2], y: [F; 2]) -> F {
    const NUM_ROUNDS: usize = 14;

    // Use helper structs for constants and keys
    let constants = AnemoiConstants::new();
    let round_keys = AnemoiRoundKeys::new();

    // Set initial state
    let mut current_x = x;
    let mut current_y = y;
    let sum_before_perm = current_x[0] + current_x[1] + current_y[0] + current_y[1];

    // Anemoi rounds - apply transformations using helper function
    for round in 0..NUM_ROUNDS {
        apply_anemoi_round(&mut current_x, &mut current_y, round, &constants, &round_keys);
    }

    // Apply final transformations using helper function
    apply_final_transformations(&mut current_x, &mut current_y, &constants);

    let sum_after_perm = current_x[0] + current_x[1] + current_y[0] + current_y[1];

    // Jive output: sum of input and output of permutation
    sum_before_perm + sum_after_perm
}