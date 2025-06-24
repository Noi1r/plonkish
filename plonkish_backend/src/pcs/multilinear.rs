use crate::{
    poly::multilinear::MultilinearPolynomial,
    util::{arithmetic::Field, end_timer, izip, parallel::parallelize, start_timer, Itertools},
    Error,
};

mod brakedown;
mod gemini;
mod hyrax;
mod ipa;
mod kzg;
mod zeromorph;

pub use brakedown::{
    MultilinearBrakedown, MultilinearBrakedownCommitment, MultilinearBrakedownParam,
};
pub use gemini::Gemini;
pub use hyrax::{MultilinearHyrax, MultilinearHyraxCommitment, MultilinearHyraxParam};
pub use ipa::{MultilinearIpa, MultilinearIpaCommitment, MultilinearIpaParam};
pub use kzg::{
    MultilinearKzg, MultilinearKzgCommitment, MultilinearKzgParam, MultilinearKzgProverParam,
    MultilinearKzgVerifierParam,
};
pub use zeromorph::{Zeromorph, ZeromorphKzgProverParam, ZeromorphKzgVerifierParam};

fn validate_input<'a, F: Field>(
    function: &str,
    param_num_vars: usize,
    polys: impl IntoIterator<Item = &'a MultilinearPolynomial<F>>,
    points: impl IntoIterator<Item = &'a Vec<F>>,
) -> Result<(), Error> {
    let polys = polys.into_iter().collect_vec();
    let points = points.into_iter().collect_vec();
    for poly in polys.iter() {
        if param_num_vars < poly.num_vars() {
            return Err(err_too_many_variates(
                function,
                param_num_vars,
                poly.num_vars(),
            ));
        }
    }
    let input_num_vars = polys
        .iter()
        .map(|poly| poly.num_vars())
        .chain(points.iter().map(|point| point.len()))
        .next()
        .expect("To have at least 1 poly or point");
    for point in points.into_iter() {
        if point.len() != input_num_vars {
            return Err(Error::InvalidPcsParam(format!(
                "Invalid point (expect point to have {input_num_vars} variates but got {})",
                point.len()
            )));
        }
    }
    Ok(())
}

fn err_too_many_variates(function: &str, upto: usize, got: usize) -> Error {
    Error::InvalidPcsParam(if function == "trim" {
        format!(
            "Too many variates to {function} (param supports variates up to {upto} but got {got})"
        )
    } else {
        format!(
            "Too many variates of poly to {function} (param supports variates up to {upto} but got {got})"
        )
    })
}

fn quotients<F: Field, T>(
    poly: &MultilinearPolynomial<F>,
    point: &[F],
    f: impl Fn(usize, Vec<F>) -> T,
) -> (Vec<T>, F) {
    assert_eq!(poly.num_vars(), point.len());

    let mut remainder = poly.evals().to_vec();
    let mut quotients = point
        .iter()
        .zip(0..poly.num_vars())
        .rev()
        .map(|(x_i, num_vars)| {
            let timer = start_timer(|| "quotients");
            // 这里的修改会改变remainder的值
            let (remaimder_lo, remainder_hi) = remainder.split_at_mut(1 << num_vars);
            let mut quotient = vec![F::ZERO; remaimder_lo.len()];

            parallelize(&mut quotient, |(quotient, start)| {
                izip!(quotient, &remaimder_lo[start..], &remainder_hi[start..])
                    .for_each(|(q, r_lo, r_hi)| *q = *r_hi - r_lo);
            });
            parallelize(remaimder_lo, |(remaimder_lo, start)| {
                izip!(remaimder_lo, &remainder_hi[start..])
                    .for_each(|(r_lo, r_hi)| *r_lo += (*r_hi - r_lo as &_) * x_i);
            });

            remainder.truncate(1 << num_vars);
            end_timer(timer);

            f(num_vars, quotient)
        })
        .collect_vec();
    quotients.reverse();

    (quotients, remainder[0])
}

mod additive {
    use crate::{
        pcs::{
            multilinear::validate_input, Additive, Evaluation, Evaluation_for_shift, Point,
            PolynomialCommitmentScheme,
        },
        piop::sum_check::{
            classic::{ClassicSumCheck, CoefficientsProver},
            eq_xy_eval, SumCheck as _, VirtualPolynomial,
        },
        poly::multilinear::{rotation_eval, MultilinearPolynomial},
        util::{
            arithmetic::{fe_to_bytes, inner_product, PrimeField},
            end_timer,
            expression::{Expression, Query, Rotation},
            start_timer,
            transcript::{TranscriptRead, TranscriptWrite},
            Itertools,
        },
        Error,
    };
    use std::{borrow::Cow, collections::{BTreeMap, HashMap}, ops::Deref, ptr::addr_of};

    type SumCheck<F> = ClassicSumCheck<CoefficientsProver<F>>;

    pub fn batch_open<F, Pcs>(
        pp: &Pcs::ProverParam,
        num_vars: usize,
        polys: Vec<&Pcs::Polynomial>,
        comms: Vec<&Pcs::Commitment>,
        points: &[Point<F, Pcs::Polynomial>],
        evals: &[Evaluation<F>],
        transcript: &mut impl TranscriptWrite<Pcs::CommitmentChunk, F>,
    ) -> Result<(), Error>
    where
        F: PrimeField,
        Pcs: PolynomialCommitmentScheme<F, Polynomial = MultilinearPolynomial<F>>,
        Pcs::Commitment: Additive<F>,
    {
        validate_input("batch open", num_vars, polys.clone(), points)?;

        if cfg!(feature = "sanity-check") {
            assert_eq!(
                points
                    .iter()
                    .map(|point| point.iter().map(fe_to_bytes::<F>).collect_vec())
                    .unique()
                    .count(),
                points.len()
            );
            for eval in evals {
                let (poly, point) = (&polys[eval.poly()], &points[eval.point()]);
                assert_eq!(poly.evaluate(point), *eval.value());
            }
        }

        let ell = evals.len().next_power_of_two().ilog2() as usize;
        let t = transcript.squeeze_challenges(ell);

        let timer = start_timer(|| "merged_polys");
        let eq_xt = MultilinearPolynomial::eq_xy(&t);
        let merged_polys = evals.iter().zip(eq_xt.evals().iter()).fold(
            vec![(F::ONE, Cow::<MultilinearPolynomial<_>>::default()); points.len()],
            |mut merged_polys, (eval, eq_xt_i)| {
                if merged_polys[eval.point()].1.is_empty() {
                    merged_polys[eval.point()] = (*eq_xt_i, Cow::Borrowed(polys[eval.poly()]));
                } else {
                    let coeff = merged_polys[eval.point()].0;
                    if coeff != F::ONE {
                        merged_polys[eval.point()].0 = F::ONE;
                        *merged_polys[eval.point()].1.to_mut() *= &coeff;
                    }
                    *merged_polys[eval.point()].1.to_mut() += (eq_xt_i, polys[eval.poly()]);
                }
                merged_polys
            },
        );
        end_timer(timer);

        let unique_merged_polys = merged_polys
            .iter()
            .unique_by(|(_, poly)| addr_of!(*poly.deref()))
            .collect_vec();
        let unique_merged_poly_indices = unique_merged_polys
            .iter()
            .enumerate()
            .map(|(idx, (_, poly))| (addr_of!(*poly.deref()), idx))
            .collect::<HashMap<_, _>>();
        let expression = merged_polys
            .iter()
            .enumerate()
            .map(|(idx, (scalar, poly))| {
                let poly = unique_merged_poly_indices[&addr_of!(*poly.deref())];
                Expression::<F>::eq_xy(idx)
                    * Expression::Polynomial(Query::new(poly, Rotation::cur()))
                    * scalar
            })
            .sum();

        let virtual_poly = VirtualPolynomial::new(
            &expression,
            unique_merged_polys.iter().map(|(_, poly)| poly.deref()),
            &[],
            points,
        );
        let tilde_gs_sum =
            inner_product(evals.iter().map(Evaluation::value), &eq_xt[..evals.len()]);
        let (g_prime_eval, challenges, _) =
            SumCheck::prove(&(), num_vars, virtual_poly, tilde_gs_sum, transcript)?;

        let timer = start_timer(|| "g_prime");
        let eq_xy_evals = points
            .iter()
            .map(|point| eq_xy_eval(&challenges, point))
            .collect_vec();
        let g_prime = merged_polys
            .into_iter()
            .zip(eq_xy_evals.iter())
            .map(|((scalar, poly), eq_xy_eval)| (scalar * eq_xy_eval, poly.into_owned()))
            .sum::<MultilinearPolynomial<_>>();
        end_timer(timer);

        let g_prime_comm = if cfg!(feature = "sanity-check") {
            let scalars = evals
                .iter()
                .zip(eq_xt.evals())
                .map(|(eval, eq_xt_i)| eq_xy_evals[eval.point()] * eq_xt_i)
                .collect_vec();
            let bases = evals.iter().map(|eval| comms[eval.poly()]);
            Pcs::Commitment::msm(&scalars, bases)
        } else {
            Pcs::Commitment::default()
        };

        Pcs::open(
            pp,
            &g_prime,
            &g_prime_comm,
            &challenges,
            &g_prime_eval,
            transcript,
        )
    }

    pub fn batch_open_for_shift<F, Pcs>(
        pp: &Pcs::ProverParam,
        num_vars: usize,
        polys: Vec<&Pcs::Polynomial>,
        comms: Vec<&Pcs::Commitment>,
        points: &[Point<F, Pcs::Polynomial>],
        evals: &[Evaluation_for_shift<F>],
        transcript: &mut impl TranscriptWrite<Pcs::CommitmentChunk, F>,
    ) -> Result<(), Error>
    where
        F: PrimeField,
        Pcs: PolynomialCommitmentScheme<F, Polynomial = MultilinearPolynomial<F>>,
        Pcs::Commitment: Additive<F>,
    {
        // 验证多项式和点的维度
        validate_input("batch open", num_vars, polys.clone(), points)?;

        if cfg!(feature = "sanity-check") {
            assert_eq!(
                points
                    .iter()
                    .map(|point| point.iter().map(fe_to_bytes::<F>).collect_vec())
                    .unique()
                    .count(),
                points.len()
            );
        }

        // --- Part 1: Handle Rotation::cur() evaluations ---
        let ell = evals.len().next_power_of_two().ilog2() as usize;
        let t = transcript.squeeze_challenges(ell);

        let timer = start_timer(|| "merged_polys (cur)");
        let eq_xt = MultilinearPolynomial::eq_xy(&t);

        let evals_cur = evals
            .iter()
            .filter(|eval| eval.rotation() == Rotation::cur())
            .collect_vec();

        // Check if there are any current rotation evaluations
        if !evals_cur.is_empty() {
            let merged_polys_cur = evals_cur.iter().zip(eq_xt.evals().iter()).fold(
                (F::ONE, Cow::<MultilinearPolynomial<_>>::default()),
                |mut merged_polys, (eval, eq_xt_i)| {
                    let poly_ref = polys[eval.poly()];
                    if merged_polys.1.is_empty() {
                        merged_polys = (*eq_xt_i, Cow::Borrowed(poly_ref));
                    } else {
                        let coeff = merged_polys.0;
                        if coeff != F::ONE {
                            merged_polys.0 = F::ONE;
                            *merged_polys.1.to_mut() *= &coeff;
                        }
                        // Ensure the polynomial being added has the correct number of variables
                        assert_eq!(
                            merged_polys.1.num_vars(),
                            poly_ref.num_vars(),
                            "Mismatched num_vars in merging"
                        );
                        *merged_polys.1.to_mut() += (*eq_xt_i, poly_ref);
                    }
                    merged_polys
                },
            );
            end_timer(timer);

            let expression_cur = Expression::<F>::eq_xy(0)
                * Expression::Polynomial(Query::new(0, Rotation::cur()))
                * merged_polys_cur.0;

            // Ensure merged polynomial has correct num_vars if it wasn't empty
            if !merged_polys_cur.1.is_empty() {
                assert_eq!(
                    merged_polys_cur.1.num_vars(),
                    num_vars,
                    "Merged polynomial has incorrect num_vars"
                );
            }

            let virtual_poly_cur = VirtualPolynomial::new(
                &expression_cur,
                if merged_polys_cur.1.is_empty() {
                    vec![] // Handle case where there are no Rotation::cur evaluations gracefully
                } else {
                    vec![merged_polys_cur.1.deref()]
                },
                &[],
                points,
            );

            let tilde_gs_sum_cur = inner_product(
                evals_cur.iter().map(|eval| eval.value()),
                &eq_xt[..evals_cur.len()],
            );

            let commitment_cur = if cfg!(feature = "sanity-check") {
                let scalars = evals_cur
                    .iter()
                    .zip(eq_xt.evals())
                    .map(|(eval, eq_xt_i)| *eq_xt_i) // Dereference eq_xt_i
                    .collect_vec();
                let bases = evals_cur.iter().map(|eval| comms[eval.poly()]);
                Pcs::Commitment::msm(&scalars, bases)
            } else {
                Pcs::Commitment::default()
            };

            Pcs::open(
                pp,
                &merged_polys_cur.1,
                &commitment_cur,
                &points[0],
                &tilde_gs_sum_cur,
                transcript,
            )?;

            // println!("open_for_shift_cur done");
        } // End if !evals_cur.is_empty()

        // --- Part 2: Handle Rotated evaluations ---
        let eval_rotations = evals
        .iter()
        .filter(|eval| eval.rotation() != Rotation::cur())
        .collect_vec();

        if eval_rotations.is_empty() {
            return Ok(()); // No rotated evaluations to process
        }

        // Group evaluations by rotation
        let mut evals_by_rotation: BTreeMap<Rotation, Vec<&Evaluation_for_shift<F>>> = BTreeMap::new();
        for eval in eval_rotations {
            evals_by_rotation.entry(eval.rotation()).or_default().push(eval);
        }

            // Process each rotation group
        for (rotation, rotation_evals) in evals_by_rotation {
                // Extract necessary info for this rotation group
                let polys_rotated: Vec<&MultilinearPolynomial<F>> = rotation_evals.iter().map(|eval| polys[eval.poly()]).collect();
                let comms_rotated: Vec<&Pcs::Commitment> = rotation_evals.iter().map(|eval| comms[eval.poly()]).collect();
                // Get references to the values
                let values_rotated: Vec<&F> = rotation_evals.iter().map(|eval| eval.value()).collect();

                if polys_rotated.is_empty() {
                    continue; // Skip if no polynomials for this rotation
                }

                // Squeeze challenges *only* for combining polynomials within this group
                let num_rotated = rotation_evals.len();

                let ell_rotated = if num_rotated == 1 {
                    2
                } else {
                    num_rotated.next_power_of_two().ilog2() as usize
                };

                let challenges_rotated_combine = transcript.squeeze_challenges(ell_rotated);
                let eq_xt_rotated = MultilinearPolynomial::eq_xy(&challenges_rotated_combine);
                
                // Combine polynomials for the current rotation
                let timer = start_timer(|| format!("merged_polys ({:?})", rotation));
        
                let merged_poly_rotated_cow = rotation_evals.iter().zip(eq_xt_rotated.evals().iter()).fold(
                    // Initialize with scalar 1 and default (empty) polynomial
                    (F::ONE, Cow::<MultilinearPolynomial<_>>::default()),
                    |mut merged, (eval, eq_xt_i)| {
                        let poly_ref = polys[eval.poly()];
                        if merged.1.is_empty() {
                            // First polynomial, borrow it with the coefficient
                            merged = (*eq_xt_i, Cow::Borrowed(poly_ref));
                        } else {
                            // Subsequent polynomials, ensure we have a mutable owned version
                            let coeff = merged.0;
                            if coeff != F::ONE {
                                // Apply previous scalar before adding the new poly
                                *merged.1.to_mut() *= &coeff;
                                merged.0 = F::ONE; // Reset scalar after applying
                            }
                            // Add the new polynomial scaled by its eq_xt_i coefficient
                            assert_eq!(merged.1.num_vars(), poly_ref.num_vars(), "Mismatched num_vars in merging rotated");
                            *merged.1.to_mut() += (*eq_xt_i, poly_ref);
                        }
                        merged
                    },
                );
                end_timer(timer);


                // Handle the final scalar if it wasn't applied in the loop
                let (merged_scalar, merged_poly_cow) = merged_poly_rotated_cow;
                let mut merged_poly_owned = merged_poly_cow.into_owned(); // Get owned version
                if merged_scalar != F::ONE {
                    merged_poly_owned *= &merged_scalar; // Apply final scalar
                }

                // Compute the combined evaluation value for the merged polynomial
                // Pass references using eval.value() and copy eq_xt values
                let merged_value = inner_product(
                    values_rotated.iter().copied(), // Dereference to get &F
                    eq_xt_rotated[..num_rotated].iter(), // Pass iterator of F
                );
                // Apply the overall scalar to the combined value as well
                // let final_merged_value = merged_value * merged_scalar;

                // Calculate the commitment to the merged polynomial for this rotation using MSM
                let merged_comm_rotated = if cfg!(feature = "sanity-check") {
                    // Calculate scalars for MSM: eq_xt_i * overall_scalar
                    let scalars = eq_xt_rotated.evals()[..num_rotated]
                        .iter()
                        .map(|eq_val| *eq_val * merged_scalar)
                        .collect_vec();
                    Pcs::Commitment::msm(&scalars, comms_rotated) // Ensure signature matches Additive trait
                } else {
                    // Should not happen if polys_rotated is not empty, but provide default
                    Pcs::Commitment::default()
                };

                // --- Apply Zeromorph Logic (Adapted) ---
                // Defer challenge squeezing (y, z) and a_0 calculation to the PCS function.
                // Pass the transcript mutable reference.

                Pcs::prove_shifted_evaluation(
                    pp,
                    &merged_poly_owned, // The combined polynomial f (already scaled)
                    &merged_comm_rotated, // Commitment to f
                    &points[0], // Target evaluation point u (Vec<F>)
                    &merged_value, // Target evaluation value v = f_shifted(u) (scaled)
                    &rotation, // The specific rotation being proven
                    transcript, // Pass the transcript for internal challenge squeezing
                )?;
                /* --- End Placeholder --- */
            } // End loop over rotations

        Ok(())
    }

    pub fn batch_verify<F, Pcs>(
        vp: &Pcs::VerifierParam,
        num_vars: usize,
        comms: Vec<&Pcs::Commitment>,
        points: &[Point<F, Pcs::Polynomial>],
        evals: &[Evaluation<F>],
        transcript: &mut impl TranscriptRead<Pcs::CommitmentChunk, F>,
    ) -> Result<(), Error>
    where
        F: PrimeField,
        Pcs: PolynomialCommitmentScheme<F, Polynomial = MultilinearPolynomial<F>>,
        Pcs::Commitment: Additive<F>,
    {
        validate_input("batch verify", num_vars, [], points)?;

        let ell = evals.len().next_power_of_two().ilog2() as usize;
        let t = transcript.squeeze_challenges(ell);

        let eq_xt = MultilinearPolynomial::eq_xy(&t);
        let tilde_gs_sum =
            inner_product(evals.iter().map(Evaluation::value), &eq_xt[..evals.len()]);
        let (g_prime_eval, challenges) =
            SumCheck::verify(&(), num_vars, 2, tilde_gs_sum, transcript)?;

        let eq_xy_evals = points
            .iter()
            .map(|point| eq_xy_eval(&challenges, point))
            .collect_vec();
        let g_prime_comm = {
            let scalars = evals
                .iter()
                .zip(eq_xt.evals())
                .map(|(eval, eq_xt_i)| eq_xy_evals[eval.point()] * eq_xt_i)
                .collect_vec();
            let bases = evals.iter().map(|eval| comms[eval.poly()]);
            Pcs::Commitment::msm(&scalars, bases)
        };
        Pcs::verify(vp, &g_prime_comm, &challenges, &g_prime_eval, transcript)
    }

    pub fn batch_verify_for_shift<F, Pcs>(
        vp: &Pcs::VerifierParam,
        num_vars: usize,
        comms: Vec<&Pcs::Commitment>,
        points: &[Point<F, Pcs::Polynomial>],
        evals: &[Evaluation_for_shift<F>],
        transcript: &mut impl TranscriptRead<Pcs::CommitmentChunk, F>,
    ) -> Result<(), Error>
    where
        F: PrimeField,
        Pcs: PolynomialCommitmentScheme<F, Polynomial = MultilinearPolynomial<F>>,
        Pcs::Commitment: Additive<F>,
    {
        validate_input("batch verify", num_vars, [], points)?;

        let ell = evals.len().next_power_of_two().ilog2() as usize;
        let t = transcript.squeeze_challenges(ell);

        let eq_xt = MultilinearPolynomial::eq_xy(&t);

        let evals_cur = evals
            .iter()
            .filter(|eval| eval.rotation() == Rotation::cur())
            .collect_vec();

        let tilde_gs_sum = inner_product(
            evals_cur.iter().map(|eval| eval.value()),
            &eq_xt[..evals_cur.len()],
        );

        let commitment_cur = {
            let scalars = evals_cur
                .iter()
                .zip(eq_xt.evals())
                .map(|(eval, eq_xt_i)| *eq_xt_i) // Dereference eq_xt_i
                .collect_vec();
            let bases = evals_cur.iter().map(|eval| comms[eval.poly()]);
            Pcs::Commitment::msm(&scalars, bases)
        };

        Pcs::verify(vp, &commitment_cur, &points[0], &tilde_gs_sum, transcript).is_ok();


        // --- Part 2: 验证旋转的求值 ---
        let eval_rotations = evals
        .iter()
        .filter(|eval| eval.rotation() != Rotation::cur())
        .collect_vec();

       if !eval_rotations.is_empty() {
           // 按 rotation 分组
           let mut evals_by_rotation : BTreeMap<Rotation, Vec<&Evaluation_for_shift<F>>> = BTreeMap::new();
           for eval in eval_rotations {
               evals_by_rotation.entry(eval.rotation()).or_default().push(eval);
           }

           // 逐个处理 rotation 分组
           for (rotation, rotation_evals) in evals_by_rotation {
               let num_rotated = rotation_evals.len();

               let ell_rotated = if num_rotated == 1 {
                    2
                } else {
                    num_rotated.next_power_of_two().ilog2() as usize
                };

               let challenges_rotated_combine = transcript.squeeze_challenges(ell_rotated);
               let eq_xt_rotated = MultilinearPolynomial::eq_xy(&challenges_rotated_combine);

               // 计算该分组的合并求值结果
               let merged_value = inner_product(
                   rotation_evals.iter().map(|eval| eval.value()), // 输入 &F
                   eq_xt_rotated.evals()[..num_rotated].iter(), // 需要 F
               );

               // 计算该分组的合并承诺
               let merged_comm_rotated = {
                   let scalars = eq_xt_rotated.evals()[..num_rotated].to_vec();
                   let bases = rotation_evals.iter().map(|eval| comms[eval.poly()]);
                   Pcs::Commitment::msm(&scalars, bases)
               };

               // 调用特定于 Pcs 的移位验证函数
               Pcs::verify_shifted_evaluation(
                    vp,
                    &merged_comm_rotated, // 对合并后的 f 的承诺
                    &points[0],    // 求值点 u
                    &merged_value,       // 合并后的声称值 v = f_d(u)
                    &rotation,           // 当前的移位信息
                    transcript,          // 包含证明数据的 transcript
                )?;
           }
       }

       Ok(())
   }
}
