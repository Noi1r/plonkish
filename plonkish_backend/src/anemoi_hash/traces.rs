use crate::util::arithmetic::PrimeField;
use std::fmt::{Debug, Formatter, Result};

/// The structure for the trace of the Anemoi-Jive sponge hash function.
#[derive(Clone)]
pub struct AnemoiVLHTrace<F: PrimeField, const N: usize, const NUM_ROUNDS: usize> {
    /// The input sequence.
    pub input: Vec<F>,
    /// The state before each permutation.
    pub before_permutation: Vec<([F; N], [F; N])>,
    /// The intermediate values for each permutation.
    pub intermediate_values_before_constant_additions:
        Vec<([[F; N]; NUM_ROUNDS], [[F; N]; NUM_ROUNDS])>,
    /// The state after each permutation.
    pub after_permutation: Vec<([F; N], [F; N])>,
    /// The output.
    pub output: F,
}

impl<F: PrimeField, const N: usize, const NUM_ROUNDS: usize> Default
    for AnemoiVLHTrace<F, N, NUM_ROUNDS>
{
    fn default() -> Self {
        Self {
            input: vec![],
            before_permutation: vec![],
            intermediate_values_before_constant_additions: vec![],
            after_permutation: vec![],
            output: F::default(),
        }
    }
}

impl<F: PrimeField, const N: usize, const NUM_ROUNDS: usize> Debug
    for AnemoiVLHTrace<F, N, NUM_ROUNDS>
{
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        f.write_str("input:\n")?;
        for (i, elem) in self.input.iter().enumerate() {
            f.write_fmt(format_args!("\r x[{}] = {:?}\n", i, elem))?;
        }

        let chunk_len = if self.input.len() % (2 * N - 1) == 0 {
            self.input.len() / (2 * N - 1)
        } else {
            self.input.len() / (2 * N - 1) + 1
        };

        for i in 0..chunk_len {
            f.write_fmt(format_args!("before permutation: {}\n", i))?;

            for (i, elem) in self.before_permutation[i].0.iter().enumerate() {
                f.write_fmt(format_args!("\r\r x[{}] = {:?}\n", i, elem))?;
            }

            for (i, elem) in self.before_permutation[i].1.iter().enumerate() {
                f.write_fmt(format_args!("\r \r y[{}] = {:?}\n", i, elem))?;
            }

            for r in 0..NUM_ROUNDS {
                f.write_fmt(format_args!("round {}: intermediate permutation\n", r))?;

                for (i, elem) in self.intermediate_values_before_constant_additions[i].0[r]
                    .iter()
                    .enumerate()
                {
                    f.write_fmt(format_args!("\r\r x[{}] = {:?}\n", i, elem))?;
                }

                for (i, elem) in self.intermediate_values_before_constant_additions[i].1[r]
                    .iter()
                    .enumerate()
                {
                    f.write_fmt(format_args!("\r\r y[{}] = {:?}\n", i, elem))?;
                }
            }

            f.write_fmt(format_args!("after permutation: {}\n", i))?;

            for (i, elem) in self.after_permutation[i].0.iter().enumerate() {
                f.write_fmt(format_args!("\r\r x[{}] = {:?}\n", i, elem))?;
            }

            for (i, elem) in self.after_permutation[i].1.iter().enumerate() {
                f.write_fmt(format_args!("\r \r y[{}] = {:?}\n", i, elem))?;
            }
        }

        f.write_fmt(format_args!("output = {:?}\n", self.output))
    }
}

/// The structure for the trace of the Anemio-Jive CRH.
#[derive(Clone, PartialEq)]
pub struct JiveTrace<F: PrimeField, const N: usize, const NUM_ROUNDS: usize> {
    /// The first half of the input.
    pub input_x: [F; N],
    /// The second half of the input.
    pub input_y: [F; N],
    /// The first half of the intermediate values in the rounds.
    pub intermediate_x_before_constant_additions: [[F; N]; NUM_ROUNDS],
    /// The second half of the intermediate values in the rounds.
    pub intermediate_y_before_constant_additions: [[F; N]; NUM_ROUNDS],
    /// The first half of the final output (after the linear layer).
    pub final_x: [F; N],
    /// The second half of the final output (after the linear layer).
    pub final_y: [F; N],
    /// The output of the Jive CRH.
    pub output: F,
}

impl<F: PrimeField, const N: usize, const NUM_ROUNDS: usize> Default
    for JiveTrace<F, N, NUM_ROUNDS>
{
    fn default() -> Self {
        Self {
            input_x: [F::default(); N],
            input_y: [F::default(); N],
            intermediate_x_before_constant_additions: [[F::default(); N]; NUM_ROUNDS],
            intermediate_y_before_constant_additions: [[F::default(); N]; NUM_ROUNDS],
            final_x: [F::default(); N],
            final_y: [F::default(); N],
            output: F::default(),
        }
    }
}

impl<F: PrimeField, const N: usize, const NUM_ROUNDS: usize> Debug for JiveTrace<F, N, NUM_ROUNDS> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        f.write_str("input_x:\n")?;
        for (i, elem) in self.input_x.iter().enumerate() {
            f.write_fmt(format_args!("\r x[{}] = {:?}\n", i, elem))?;
        }

        f.write_str("input_y:\n")?;
        for (i, elem) in self.input_y.iter().enumerate() {
            f.write_fmt(format_args!("\r y[{}] = {:?}\n", i, elem))?;
        }

        for r in 0..NUM_ROUNDS {
            f.write_fmt(format_args!("round {}:\n", r))?;

            for (i, elem) in self.intermediate_x_before_constant_additions[r]
                .iter()
                .enumerate()
            {
                f.write_fmt(format_args!("\r x[{}] = {:?}\n", i, elem))?;
            }

            for (i, elem) in self.intermediate_y_before_constant_additions[r]
                .iter()
                .enumerate()
            {
                f.write_fmt(format_args!("\r y[{}] = {:?}\n", i, elem))?;
            }
        }

        f.write_str("final_x:\n")?;
        for (i, elem) in self.final_x.iter().enumerate() {
            f.write_fmt(format_args!("\r x[{}] = {:?}\n", i, elem))?;
        }

        f.write_str("final_y:\n")?;
        for (i, elem) in self.final_y.iter().enumerate() {
            f.write_fmt(format_args!("\r y[{}] = {:?}\n", i, elem))?;
        }

        f.write_fmt(format_args!("output: {:?}\n", self.output))
    }
}

/// The structure for the trace of the Anemoi-Jive stream cipher.
#[derive(Clone)]
pub struct AnemoiStreamCipherTrace<F: PrimeField, const N: usize, const NUM_ROUNDS: usize> {
    /// The input sequence.
    pub input: Vec<F>,
    /// The state before each permutation.
    pub before_permutation: Vec<([F; N], [F; N])>,
    /// The intermediate values for each permutation.
    pub intermediate_values_before_constant_additions:
        Vec<([[F; N]; NUM_ROUNDS], [[F; N]; NUM_ROUNDS])>,
    /// The state after each permutation.
    pub after_permutation: Vec<([F; N], [F; N])>,
    /// The output.
    pub output: Vec<F>,
}

impl<F: PrimeField, const N: usize, const NUM_ROUNDS: usize> Default
    for AnemoiStreamCipherTrace<F, N, NUM_ROUNDS>
{
    fn default() -> Self {
        Self {
            input: vec![],
            before_permutation: vec![],
            intermediate_values_before_constant_additions: vec![],
            after_permutation: vec![],
            output: vec![],
        }
    }
}

impl<F: PrimeField, const N: usize, const NUM_ROUNDS: usize> Debug
    for AnemoiStreamCipherTrace<F, N, NUM_ROUNDS>
{
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        f.write_str("input:\n")?;
        for (i, elem) in self.input.iter().enumerate() {
            f.write_fmt(format_args!("\r x[{}] = {:?}\n", i, elem))?;
        }
        let chunk_len = if self.input.len() % (2 * N - 1) == 0 {
            self.input.len() / (2 * N - 1)
        } else {
            self.input.len() / (2 * N - 1) + 1
        };

        for i in 0..chunk_len {
            f.write_fmt(format_args!("before permutation: {}\n", i))?;

            for (i, elem) in self.before_permutation[i].0.iter().enumerate() {
                f.write_fmt(format_args!("\r\r x[{}] = {:?}\n", i, elem))?;
            }

            for (i, elem) in self.before_permutation[i].1.iter().enumerate() {
                f.write_fmt(format_args!("\r \r y[{}] = {:?}\n", i, elem))?;
            }

            for r in 0..NUM_ROUNDS {
                f.write_fmt(format_args!("round {}: intermediate permutation\n", r))?;

                for (i, elem) in self.intermediate_values_before_constant_additions[i].0[r]
                    .iter()
                    .enumerate()
                {
                    f.write_fmt(format_args!("\r\r x[{}] = {:?}\n", i, elem))?;
                }

                for (i, elem) in self.intermediate_values_before_constant_additions[i].1[r]
                    .iter()
                    .enumerate()
                {
                    f.write_fmt(format_args!("\r\r y[{}] = {:?}\n", i, elem))?;
                }
            }

            f.write_fmt(format_args!("after permutation: {}\n", i))?;

            for (i, elem) in self.after_permutation[i].0.iter().enumerate() {
                f.write_fmt(format_args!("\r\r x[{}] = {:?}\n", i, elem))?;
            }

            for (i, elem) in self.after_permutation[i].1.iter().enumerate() {
                f.write_fmt(format_args!("\r \r y[{}] = {:?}\n", i, elem))?;
            }
        }

        f.write_str("output:\n")?;
        for (i, elem) in self.output.iter().enumerate() {
            f.write_fmt(format_args!("\r x[{}] = {:?}\n", i, elem))?;
        }

        Ok(())
    }
}
