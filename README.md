# Plonkish: HyperPlonk with Zeromorph Shift Optimization

[![Rust](https://img.shields.io/badge/rust-1.70+-blue.svg)](https://www.rust-lang.org)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Build Status](https://img.shields.io/badge/build-passing-brightgreen.svg)](#)

**Plonkish** is a high-performance implementation of the HyperPlonk zero-knowledge proof system with Zeromorph shift optimization. This repository focuses on zero-knowledge proof systems, specifically optimizing the shift operations in HyperPlonk using the Zeromorph approach to reduce exponential costs to just one additional commitment.

## 🚀 Quick Start

```bash
# Clone the repository
git clone https://github.com/Noi1r/plonkish.git
cd plonkish

# Build the project
cargo build --release

# Run tests
cargo test --release

# Run benchmarks
cargo bench --bench proof_system
```

## 📋 Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Documentation](#documentation)
- [Architecture](#architecture)
- [Benchmarks](#benchmarks)
- [Contributing](#contributing)
- [License](#license)

## ✨ Features

### Core Features
- **HyperPlonk Implementation**: Complete implementation of the HyperPlonk proving system
- **Zeromorph PCS**: Advanced polynomial commitment scheme with shift support
- **Shift Optimization**: Reduces exponential shift costs to a single additional commitment
- **Anemoi Hash Integration**: Efficient hash function implementation for Merkle trees
- **Multiple Circuit Types**: Support for various zero-knowledge circuits



## 🛠 Installation

### Prerequisites

- **Rust**: Version 1.70 or higher
- **System Requirements**: 
  - Memory: 8GB RAM minimum (16GB recommended)
  - Storage: 2GB free space
  - OS: Linux, macOS, or Windows

### Install Rust

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source ~/.cargo/env
```

### Build from Source

```bash
# Clone the repository
git clone https://github.com/Noi1r/plonkish.git
cd plonkish

# Build in release mode (recommended)
cargo build --release

# Or build in development mode
cargo build
```

## 🎯 Usage

### Basic Example

```rust
use plonkish_backend::{
    backend::{PlonkishBackend, hyperplonk::HyperPlonk},
    pcs::multilinear::Zeromorph,
    util::transcript::Keccak256Transcript,
};
use halo2_curves::bn256::Bn256;

// Initialize the proving system
type ProofSystem = HyperPlonk<Zeromorph<UnivariateKzg<Bn256>>>;

// Setup phase
let circuit_info = create_circuit_info();
let param = ProofSystem::setup(&circuit_info, &mut rng)?;
let (prover_param, verifier_param) = ProofSystem::preprocess(&param, &circuit_info)?;

// Proving phase
let circuit = create_circuit();
let mut transcript = Keccak256Transcript::new();
ProofSystem::prove(&prover_param, &circuit, &mut transcript, &mut rng)?;

// Verification phase
let mut transcript = Keccak256Transcript::new();
ProofSystem::verify(&verifier_param, &instances, &mut transcript, &mut rng)?;
```

### Circuit Examples

#### Merkle Tree Membership Proof

```rust
use plonkish_backend::backend::hyperplonk::util::merkle::*;

// Create a Merkle membership proof circuit
let (circuit_info, circuit) = merkle_membership_proof_circuit::<Fr, Lexical>(
    num_vars,
    &proof_nodes,
    &positions,
    &leaf_value,
    &root_hash,
    rng
);
```

#### Anemoi Hash Circuit

```rust
use plonkish_backend::backend::hyperplonk::util::anemoi::*;

// Create an Anemoi hash circuit
let (circuit_info, circuit) = rand_anemoi_hash_circuit_with_flatten::<Fr, Lexical>(
    preprocess_rng,
    witness_rng
);
```

## Commands

### Testing
```bash
# Run specific test with output
cargo test --release --package plonkish_backend --lib -- backend::hyperplonk::test::merkle_membership_proof_zeromorph_kzg --exact --show-output 

# Run all tests
cargo test --release
```

### Benchmarking
```bash
# Basic benchmark for proof systems
cargo bench --bench proof_system -- \
    --system hyperplonk \
    --system halo2 \
    --system espresso_hyperplonk \
    --circuit vanilla_plonk \
    --k 20..24

# Benchmark with timing breakdown (requires gnuplot)
cargo bench --bench proof_system --features timer -- \
    --system hyperplonk \
    --circuit vanilla_plonk \
    --k 20..24 \
  | cargo run plotter -- -

# PCS benchmarks
cargo bench --bench pcs

# Zero check benchmarks  
cargo bench --bench zero_check
```

### Build Profiles
```bash
# Release build
cargo build --release

# Flamegraph profiling build
cargo build --profile flamegraph
```

## Architecture

### Core Components

1. **Backend Layer** (`plonkish_backend/src/backend/`)
   - `PlonkishBackend` trait: Core interface for proof systems
   - `HyperPlonk`: Main implementation with shift support via `prove_with_shift()` and `verify_with_shift()`
   - Preprocessor, Prover, and Verifier modules

2. **Polynomial Commitment Schemes** (`plonkish_backend/src/pcs/`)
   - `PolynomialCommitmentScheme` trait with shift extensions (`open_shift`, `verify_shifted_evaluation`)
   - **Zeromorph**: Primary PCS implementation supporting shifts
   - Univariate KZG: Underlying commitment scheme

3. **Frontend Layer** (`plonkish_backend/src/frontend/`)
   - Halo2 frontend integration (feature-gated)
   - Circuit abstraction via `PlonkishCircuit` trait

4. **Polynomial Representations** (`plonkish_backend/src/poly/`)
   - `MultilinearPolynomial`: Main polynomial type
   - `UnivariatePolynomial`: For KZG operations
   - Rotation evaluation support in `poly/multilinear.rs`

5. **Utilities** (`plonkish_backend/src/util/`)
   - Arithmetic operations (FFT, MSM)
   - Expression evaluation and rotation handling
   - Transcript management
   - Parallel processing utilities

### Key Shift Implementation

The repository's main contribution is the shift optimization:
- **Problem**: Original HyperPlonk shift requires exponential number of opens
- **Solution**: Zeromorph-based shift reduces cost to one additional commitment
- **Implementation**: 
  - `open_shift()` and `verify_shifted_evaluation()` in Zeromorph PCS
  - `prove_with_shift()` and `verify_with_shift()` in HyperPlonk backend
  - `Evaluation_for_shift` struct for batch operations

### Circuit Information Structure

`PlonkishCircuitInfo` contains:
- `k`: Circuit size (2^k)
- `num_instances`: Instance polynomial counts
- `preprocess_polys`: Preprocessed polynomials
- `num_witness_polys`: Witness polynomials per phase
- `constraints`: Circuit constraints as expressions
- `lookups`: Vector lookup arguments
- `permutations`: Permutation arguments

## Development Notes

- **Workspace Structure**: Two main packages - `plonkish_backend` (core) and `benchmark` (performance testing)
- **Feature Flags**: 
  - `timer`: Enables timing instrumentation
  - `parallel`: Enables rayon parallelization
  - `frontend-halo2`: Enables Halo2 circuit frontend
  - `benchmark`: Required for benchmark builds
- **Testing**: Most critical tests are in `backend::hyperplonk::test` module
- **Dependencies**: Uses custom forks of halo2 and other ZK libraries for compatibility

## Limitations

- Only Zeromorph PCS is currently supported due to shift requirements
- Other PCS implementations are disabled until shift support is added
- Requires specific versions of halo2 and related libraries from custom forks

## Acknowledgements

- Types for plonkish circuit structure are ported from https://github.com/zcash/halo2.
- Most part of [HyperPlonk](https://eprint.iacr.org/2022/1355.pdf) and multilinear KZG PCS implementation are ported from https://github.com/EspressoSystems/hyperplonk with reorganization and extension to support `halo2` constraint system.

