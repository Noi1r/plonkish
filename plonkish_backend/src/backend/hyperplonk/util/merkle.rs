use crate::{
    backend::{
        hyperplonk::util::{
            anemoi::{anemoi_hash_circuit_info_original, anemoi_jive_hash, build_jive_crh_circuit},
            common::*,
        },
        mock::MockCircuit,
        PlonkishCircuit, PlonkishCircuitInfo,
    },
    util::{
        arithmetic::PrimeField,
        expression::rotate::Rotatable,
    },
};

use rand::RngCore;

/// Helper structures for Merkle tree proofs
#[derive(Debug, Clone)]
pub struct MerkleProofNode<F: PrimeField> {
    pub position: MerklePosition,
    pub current_hash: F, // Current node hash
    pub sibling1: F,     // Left or first sibling
    pub sibling2: F,     // Right or second sibling
}

#[derive(Debug, Clone)]
pub enum MerklePosition {
    Left,   // Current node is left child
    Middle, // Current node is middle child
    Right,  // Current node is right child
}

/// Position verification gates for Merkle proof
fn build_tree_location_verification_gates<F: PrimeField>(
    proof_node: &MerkleProofNode<F>,
    usable_indices: &[usize],
    gate_counter: &mut usize,
    w1_values: &mut [F],
    w2_values: &mut [F],
    w3_values: &mut [F],
    w4_values: &mut [F],
    wo_values: &mut [F],
    q2: &mut [F],
    q3: &mut [F],
    q4: &mut [F],
    qc: &mut [F],
    qb: &mut [F],
    qm1: &mut [F],
    qm2: &mut [F],
    qo: &mut [F],
) -> (F, F, F, F) {
    // Gate 1: Boolean constraint for position indicators
    let pos_gate_idx = usable_indices[*gate_counter];
    *gate_counter += 1;

    let (is_left, is_middle, is_right) = match proof_node.position {
        MerklePosition::Left => (F::ONE, F::ZERO, F::ZERO),
        MerklePosition::Middle => (F::ZERO, F::ONE, F::ZERO),
        MerklePosition::Right => (F::ZERO, F::ZERO, F::ONE),
    };

    w1_values[pos_gate_idx] = F::ZERO;
    w2_values[pos_gate_idx] = is_left;
    w3_values[pos_gate_idx] = is_middle;
    w4_values[pos_gate_idx] = is_right;
    wo_values[pos_gate_idx] = F::ZERO;

    // Constraint: is_left + is_middle + is_right = 1
    q2[pos_gate_idx] = F::ONE;
    q3[pos_gate_idx] = F::ONE;
    q4[pos_gate_idx] = F::ONE;
    qc[pos_gate_idx] = -F::ONE;
    qb[pos_gate_idx] = F::ONE;

    // Determine the three children based on position
    let (left_hash, middle_hash, right_hash) = match proof_node.position {
        MerklePosition::Left => (
            proof_node.current_hash,
            proof_node.sibling1,
            proof_node.sibling2,
        ),
        MerklePosition::Middle => (
            proof_node.sibling1,
            proof_node.current_hash,
            proof_node.sibling2,
        ),
        MerklePosition::Right => (
            proof_node.sibling1,
            proof_node.sibling2,
            proof_node.current_hash,
        ),
    };

    // Gate 2: First position selection
    let sel1_gate_idx = usable_indices[*gate_counter];
    *gate_counter += 1;

    w1_values[sel1_gate_idx] = left_hash;
    w2_values[sel1_gate_idx] = is_left;
    w3_values[sel1_gate_idx] = middle_hash;
    w4_values[sel1_gate_idx] = is_middle;
    wo_values[sel1_gate_idx] = left_hash * is_left + middle_hash * is_middle;

    qm1[sel1_gate_idx] = F::ONE;
    qm2[sel1_gate_idx] = F::ONE;
    qo[sel1_gate_idx] = F::ONE;

    // Gate 3: Second position selection
    let sel2_gate_idx = usable_indices[*gate_counter];
    *gate_counter += 1;

    let partial_sum = wo_values[sel1_gate_idx];
    w1_values[sel2_gate_idx] = right_hash;
    w2_values[sel2_gate_idx] = is_right;
    w3_values[sel2_gate_idx] = partial_sum;
    w4_values[sel2_gate_idx] = F::ZERO;
    wo_values[sel2_gate_idx] = proof_node.current_hash;

    qm1[sel2_gate_idx] = F::ONE;
    q3[sel2_gate_idx] = F::ONE;
    qo[sel2_gate_idx] = F::ONE;

    (left_hash, middle_hash, right_hash, proof_node.current_hash)
}

/// Complete Merkle tree membership proof circuit with full Jive CRH implementation
/// Verifies that a leaf is part of a 3-ary Merkle tree with given root
pub fn merkle_membership_proof_circuit<F: PrimeField, R: Rotatable + From<usize>>(
    leaf_value: F,
    merkle_path: Vec<MerkleProofNode<F>>,
    expected_root: F,
    _preprocess_rng: impl RngCore,
    _witness_rng: impl RngCore,
) -> (PlonkishCircuitInfo<F>, impl PlonkishCircuit<F>) {
    const NUM_ROUNDS: usize = 14;
    const GATES_PER_JIVE: usize = NUM_ROUNDS + 4;
    const GATES_PER_LEVEL: usize = GATES_PER_JIVE + 10;

    let tree_depth = merkle_path.len();
    let padding_constants = generate_padding_constants(tree_depth);

    // Calculate total gates needed
    let total_gates = tree_depth * GATES_PER_LEVEL + 10;
    let num_vars = (total_gates as f64).log2().ceil() as usize;
    let size = 1 << num_vars;
    let usable_indices = R::from(num_vars).usable_indices();

    if usable_indices.len() < total_gates {
        panic!(
            "Not enough usable indices: {} < {}",
            usable_indices.len(),
            total_gates
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

    let mut gate_counter = 0;
    let mut current_hash = leaf_value;

    // Process each level of the Merkle path
    for (level, proof_node) in merkle_path.iter().enumerate() {
        let level_padding = padding_constants[level];

        // Position verification using helper function
        let (left_hash, middle_hash, right_hash, _) = build_tree_location_verification_gates(
            proof_node,
            &usable_indices,
            &mut gate_counter,
            &mut w1_values,
            &mut w2_values,
            &mut w3_values,
            &mut w4_values,
            &mut wo_values,
            &mut q2,
            &mut q3,
            &mut q4,
            &mut qc,
            &mut qb,
            &mut qm1,
            &mut qm2,
            &mut qo,
        );

        // Jive CRH implementation using helper function
        current_hash = build_jive_crh_circuit(
            left_hash,
            middle_hash,
            right_hash,
            level_padding,
            &constants,
            &round_keys,
            &usable_indices,
            &mut gate_counter,
            &mut w1_values,
            &mut w2_values,
            &mut w3_values,
            &mut w4_values,
            &mut wo_values,
            &mut qprk1,
            &mut qprk2,
            &mut qprk3,
            &mut qprk4,
            &mut q1,
            &mut q2,
            &mut q3,
            &mut q4,
            &mut qo,
            &mut permutation,
        );
    }

    // Final root verification gate
    let root_gate_idx = usable_indices[gate_counter];
    w1_values[root_gate_idx] = current_hash;
    w2_values[root_gate_idx] = expected_root;
    w3_values[root_gate_idx] = F::ZERO;
    w4_values[root_gate_idx] = F::ZERO;
    wo_values[root_gate_idx] = F::ZERO;

    // Constraint: current_hash - expected_root = 0
    q1[root_gate_idx] = F::ONE;
    q2[root_gate_idx] = -F::ONE;

    let circuit_info = PlonkishCircuitInfo {
        k: num_vars,
        num_instances: vec![],
        preprocess_polys: vec![
            q1, q2, q3, q4, qo, qm1, qm2, qc, qecc, qb, qprk1, qprk2, qprk3, qprk4,
        ],
        num_witness_polys: vec![5],
        num_challenges: vec![0],
        constraints: anemoi_hash_circuit_info_original::<F>(num_vars).constraints,
        lookups: Vec::new(),
        permutations: permutation.into_cycles(),
        max_degree: Some(7),
    };

    let witness = vec![w1_values, w2_values, w3_values, w4_values, wo_values];

    (circuit_info, MockCircuit::new(vec![], witness))
}

/// Generate random Merkle membership proof test data using generic PrimeField
/// Returns (leaf_value, merkle_path, expected_root)
pub fn generate_random_merkle_proof<F: PrimeField>(
    depth: usize,
    mut rng: impl RngCore,
) -> (F, Vec<MerkleProofNode<F>>, F) {
    // Generate padding constants for domain separation
    let padding_constants = generate_padding_constants(depth);

    // Generate random leaf value
    let leaf_value = F::random(&mut rng);

    // Generate random path (position at each level from leaf to root)
    let path_positions: Vec<MerklePosition> = (0..depth)
        .map(|_| match rng.next_u32() % 3 {
            0 => MerklePosition::Left,
            1 => MerklePosition::Middle,
            _ => MerklePosition::Right,
        })
        .collect();

    // Build merkle path from leaf to root
    let mut current_hash = leaf_value;
    let mut merkle_path = Vec::new();

    for (level, position) in path_positions.iter().enumerate() {
        // Generate random siblings
        let sibling1 = F::random(&mut rng);
        let sibling2 = F::random(&mut rng);

        // Determine the three children based on current node's position
        let (left_hash, middle_hash, right_hash) = match position {
            MerklePosition::Left => (current_hash, sibling1, sibling2),
            MerklePosition::Middle => (sibling1, current_hash, sibling2),
            MerklePosition::Right => (sibling1, sibling2, current_hash),
        };

        // Use level-specific padding constant for domain separation
        let level_padding = padding_constants[level];

        // Compute parent hash using Anemoi Jive CRH
        // This follows the 3-ary Merkle tree construction in Section 7.1
        let parent_hash = {
            let x = [left_hash, middle_hash]; // First 2 children
            let y = [right_hash, level_padding]; // Third child + domain separation
            anemoi_jive_hash(x, y)
        };

        // Add to merkle path (this represents the current level)
        merkle_path.push(MerkleProofNode {
            position: position.clone(),
            current_hash, // Hash of current node (before going up to parent)
            sibling1,     // First sibling
            sibling2,     // Second sibling
        });

        // Move up to parent for next iteration
        current_hash = parent_hash;
    }

    let expected_root = current_hash;

    (leaf_value, merkle_path, expected_root)
}

/// Generate random Merkle membership proof circuit with proper generic implementation
pub fn rand_merkle_membership_proof_circuit<F: PrimeField, R: Rotatable + From<usize>>(
    depth: usize,
    preprocess_rng: impl RngCore,
    mut witness_rng: impl RngCore,
) -> (PlonkishCircuitInfo<F>, impl PlonkishCircuit<F>) {
    // Generate test data using generic implementation
    let (leaf_value, merkle_path, expected_root) =
        generate_random_merkle_proof(depth, &mut witness_rng);

    merkle_membership_proof_circuit::<F, R>(
        leaf_value,
        merkle_path,
        expected_root,
        preprocess_rng,
        witness_rng,
    )
}