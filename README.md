# HyperPlonk with Zeromorph

This repository is a fork of [plonkish](https://github.com/han0110/plonkish) with Zeromorph shift implementation.

In the original plonkish, the shift is implemented by using the method in [HyperPlonk](https://eprint.iacr.org/2022/1355.pdf), which will result in exponential cost(requires an exponential number of opens), with the shift approach in [Zeromorph](https://eprint.iacr.org/2023/917.pdf) the cost is reduced to just one more commitment(plonkish_backend/src/pcs/multilinear/zeromorph.rs).

We also implement the [anemoi hash and jive crh](https://eprint.iacr.org/2023/1043.pdf) circuits in the plonkish backend(plonkish_backend/src/backend/hyperplonk/util.rs).

[TODO] because of the shift approach requirement, the other pcs are not supported yet. Only the zeromorph pcs is supported.

## How to run

```sh
cargo test --release --package plonkish_backend --lib -- backend::hyperplonk::test::merkle_membership_proof_zeromorph_kzg --exact --show-output
```

## Benchmark

### Proof systems

On different proof systems with KZG polynomial commitment scheme.

```sh
Usage: cargo bench --bench proof_system -- [OPTIONS]

Options:
  --system <SYSTEM>    Proof system(s) to run. [possible values: hyperplonk, halo2, espresso_hyperplonk]
  --circuit <CIRCUIT>  Circuit to run. [possible values: vanilla_plonk, aggregation]
  --k <K>              (Range of) log number of rows.
```

For example to compare different proof systems on vanilla PLONK, run:

```sh
cargo bench --bench proof_system -- \
    --system hyperplonk \
    --system halo2 \
    --system espresso_hyperplonk \
    --circuit vanilla_plonk \
    --k 20..24
```

Then the proving time (in millisecond) will be written to `target/bench/{hyperplonk,halo2,espresso_hyperplonk}` respectively.

To further see cost breakdown of proving time without witness collecting time, run the same bench commanad with an extra cargo flag `--features timer`, then pipe the output to plotter `cargo run plotter -- -`, and the result will be rendered in `target/bench`. For example:

```sh
cargo bench --bench proof_system --features timer -- ... \
  | cargo run plotter -- -
```

Note that `plotter` requires `gnuplot` installed already.

## Acknowledgements

- Types for plonkish circuit structure are ported from https://github.com/zcash/halo2.
- Most part of [HyperPlonk](https://eprint.iacr.org/2022/1355.pdf) and multilinear KZG PCS implementation are ported from https://github.com/EspressoSystems/hyperplonk with reorganization and extension to support `halo2` constraint system.
- Most part of [Brakedown](https://eprint.iacr.org/2021/1043.pdf) specification and multilinear PCS implementation are ported from https://github.com/conroi/lcpc.
