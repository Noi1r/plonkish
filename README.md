# HyperPlonk with Zeromorph


This is a fork of [plonkish](https://github.com/han0110/plonkish) implementing HyperPlonk with Zeromorph shift optimization. The repository focuses on zero-knowledge proof systems, specifically optimizing the shift operations in HyperPlonk using the Zeromorph approach to reduce exponential costs to just one additional commitment.

Key implementations:
- **Zeromorph PCS**: `plonkish_backend/src/pcs/multilinear/zeromorph.rs` - Main polynomial commitment scheme with shift support
- **HyperPlonk Backend**: `plonkish_backend/src/backend/hyperplonk/` - Complete HyperPlonk implementation
- **Anemoi Hash**: `plonkish_backend/src/anemoi_hash/` and `plonkish_backend/src/backend/hyperplonk/util.rs` - Anemoi hash and jive CRH circuits

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
- Most part of [Brakedown](https://eprint.iacr.org/2021/1043.pdf) specification and multilinear PCS implementation are ported from https://github.com/conroi/lcpc.
