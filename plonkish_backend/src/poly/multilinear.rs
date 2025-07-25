use crate::{
    pcs::Additive,
    poly::Polynomial,
    util::{
        arithmetic::{div_ceil, usize_from_bits_le, Field},
        chain,
        expression::{rotate::BinaryField, Rotation},
        impl_index, izip_eq,
        parallel::{num_threads, parallelize, parallelize_iter},
        BitIndex, Deserialize, Itertools, Serialize,
    },
};
use num_integer::Integer;
use rand::RngCore;
use std::{
    borrow::{Borrow, Cow},
    iter::{self, Sum},
    mem,
    ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign},
};

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct MultilinearPolynomial<F> {
    evals: Vec<F>,
    num_vars: usize,
}

impl<F> Default for MultilinearPolynomial<F> {
    fn default() -> Self {
        MultilinearPolynomial::zero()
    }
}

impl<F: Field> Additive<F> for MultilinearPolynomial<F> {
    fn msm<'a, 'b>(
        scalars: impl IntoIterator<Item = &'a F>,
        bases: impl IntoIterator<Item = &'b Self>,
    ) -> Self
    where
        Self: 'b,
    {
        izip_eq!(scalars, bases).sum()
    }
}

impl<F> MultilinearPolynomial<F> {
    pub fn new(evals: Vec<F>) -> Self {
        let num_vars = if evals.is_empty() {
            0
        } else {
            let num_vars = evals.len().ilog2() as usize;
            assert_eq!(evals.len(), 1 << num_vars);
            num_vars
        };

        Self { evals, num_vars }
    }

    pub const fn zero() -> Self {
        Self {
            evals: Vec::new(),
            num_vars: 0,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.evals.is_empty()
    }

    pub fn num_vars(&self) -> usize {
        self.num_vars
    }

    pub fn evals(&self) -> &[F] {
        &self.evals
    }

    pub fn into_evals(self) -> Vec<F> {
        self.evals
    }

    pub fn iter(&self) -> impl Iterator<Item = &F> {
        self.evals.iter()
    }
}

impl<F: Field> Polynomial<F> for MultilinearPolynomial<F> {
    type Point = Vec<F>;

    fn coeffs(&self) -> &[F] {
        self.evals.as_slice()
    }

    fn evaluate(&self, point: &Self::Point) -> F {
        MultilinearPolynomial::evaluate(self, point)
    }

    #[cfg(any(test, feature = "benchmark"))]
    fn rand(n: usize, rng: impl rand::RngCore) -> Self {
        Self::new(crate::util::test::rand_vec(n, rng))
    }

    #[cfg(any(test, feature = "benchmark"))]
    fn rand_point(k: usize, rng: impl rand::RngCore) -> Self::Point {
        crate::util::test::rand_vec(k, rng)
    }

    #[cfg(any(test, feature = "benchmark"))]
    fn squeeze_point(
        k: usize,
        transcript: &mut impl crate::util::transcript::FieldTranscript<F>,
    ) -> Self::Point {
        iter::repeat_with(|| transcript.squeeze_challenge())
            .take(k)
            .collect()
    }
}

impl<F: Field> MultilinearPolynomial<F> {
    pub fn eq_xy(y: &[F]) -> Self {
        if y.is_empty() {
            return Self::zero();
        }

        let expand_serial = |next_evals: &mut [F], evals: &[F], y_i: &F| {
            for (next_evals, eval) in next_evals.chunks_mut(2).zip(evals.iter()) {
                next_evals[1] = *eval * y_i;
                next_evals[0] = *eval - &next_evals[1];
            }
        };

        let mut evals = vec![F::ONE];
        for y_i in y.iter().rev() {
            let mut next_evals = vec![F::ZERO; 2 * evals.len()];
            if evals.len() < 32 {
                expand_serial(&mut next_evals, &evals, y_i);
            } else {
                let mut chunk_size = div_ceil(evals.len(), num_threads());
                if chunk_size.is_odd() {
                    chunk_size += 1;
                }
                parallelize_iter(
                    next_evals
                        .chunks_mut(chunk_size)
                        .zip(evals.chunks(chunk_size >> 1)),
                    |(next_evals, evals)| expand_serial(next_evals, evals, y_i),
                );
            }
            evals = next_evals;
        }

        Self {
            evals,
            num_vars: y.len(),
        }
    }

    pub fn rand(num_vars: usize, mut rng: impl RngCore) -> Self {
        Self::new(
            iter::repeat_with(|| F::random(&mut rng))
                .take(1 << num_vars)
                .collect(),
        )
    }

    pub fn evaluate(&self, x: &[F]) -> F {
        assert_eq!(x.len(), self.num_vars);
        evaluate(&self.evals, x)
    }

    pub fn fix_last_vars(&self, x: &[F]) -> Self {
        assert!(!x.is_empty() && x.len() <= self.num_vars);

        let mut output = self.evals.clone();
        for (x_i, len) in x.iter().rev().zip((1..).map(|i| 1 << (self.num_vars - i))) {
            let (lo, hi) = output.split_at_mut(len);
            let chunk_size = div_ceil(len, num_threads());
            parallelize_iter(
                lo.chunks_mut(chunk_size).zip(hi.chunks(chunk_size)),
                |(lo, hi)| {
                    lo.iter_mut()
                        .zip(hi.iter())
                        .for_each(|(lo, hi)| *lo += (*hi - lo as &_) * x_i);
                },
            );
        }

        output.truncate(1 << (self.num_vars - x.len()));
        Self::new(output)
    }

    pub fn fix_var(&self, x_i: &F) -> Self {
        let mut output = Vec::with_capacity(1 << (self.num_vars - 1));
        merge_into(&mut output, self.evals(), x_i, 1, 0);
        Self::new(output)
    }

    pub fn fix_var_in_place(&mut self, x_i: &F, buf: &mut Self) {
        merge_into(&mut buf.evals, self.evals(), x_i, 1, 0);
        buf.num_vars = self.num_vars - 1;
        mem::swap(self, buf);
    }

    pub fn evaluate_for_rotation(&self, x: &[F], rotation: Rotation) -> Vec<F> {
        assert_eq!(x.len(), self.num_vars);
        if rotation == Rotation::cur() {
            return vec![self.evaluate(x)];
        }

        let signed_d = rotation.0;
        let abs_d = rotation.distance();
        let n_evals = 1 << self.num_vars();

        // 1. 旋转求值列表
        let mut rotated_evals = self.evals().to_vec();
        if signed_d > 0 {
            rotated_evals.rotate_left(abs_d);
        } else {
            // signed_d < 0
            rotated_evals.rotate_right(abs_d);
        }

        // 2. 构造旋转后的多项式
        let poly_rotated = MultilinearPolynomial::new(rotated_evals);

        // 3. 对旋转后的多项式在原始点 point 求值
        vec![poly_rotated.evaluate(x)]
    }

    // pub fn evaluate_for_rotation(&self, x: &[F], rotation: Rotation) -> Vec<F> {
    //     assert_eq!(x.len(), self.num_vars);
    //     if rotation == Rotation::cur() {
    //         return vec![self.evaluate(x)];
    //     }

    //     let distance = rotation.distance();
    //     let num_x = self.num_vars - distance;
    //     let mut evals = vec![F::ZERO; 1 << distance];
    //     let chunk_size = div_ceil(evals.len(), num_threads());
    //     if rotation < Rotation::cur() {
    //         let x = &x[distance..];
    //         let flipped_x = x.iter().map(flip).collect_vec();
    //         let pattern = rotation_eval_point_pattern::<false>(self.num_vars, distance);
    //         let offset_mask = (1 << self.num_vars) - (1 << num_x);
    //         parallelize_iter(
    //             evals.chunks_mut(chunk_size).zip(pattern.chunks(chunk_size)),
    //             |(evals, pattern)| {
    //                 let mut buf = Vec::with_capacity(1 << (num_x - 1));
    //                 let mut last_buf = Some(Vec::with_capacity(1 << (num_x - 1)));
    //                 for (eval, pat) in evals.iter_mut().zip(pattern.iter()) {
    //                     let offset = pat & offset_mask;
    //                     let mut evals = Cow::Borrowed(&self[offset..offset + (1 << num_x)]);
    //                     for (idx, (x_i, flipped_x_i)) in x.iter().zip(flipped_x.iter()).enumerate()
    //                     {
    //                         let x_i = if pat.nth_bit(idx) { flipped_x_i } else { x_i };
    //                         merge_into(&mut buf, &evals, x_i, 1, 0);
    //                         if let Cow::Owned(_) = evals {
    //                             mem::swap(evals.to_mut(), &mut buf);
    //                         } else {
    //                             evals = mem::replace(&mut buf, last_buf.take().unwrap()).into();
    //                         }
    //                     }
    //                     *eval = evals[0];
    //                     last_buf = Some(evals.into_owned());
    //                 }
    //             },
    //         );
    //     } else {
    //         let x = &x[..num_x];
    //         let flipped_x = x.iter().map(flip).collect_vec();
    //         let pattern = rotation_eval_point_pattern::<true>(self.num_vars, distance);
    //         let skip_mask = (1 << distance) - 1;
    //         parallelize_iter(
    //             evals.chunks_mut(chunk_size).zip(pattern.chunks(chunk_size)),
    //             |(evals, pattern)| {
    //                 let mut buf = Vec::with_capacity(1 << (num_x - 1));
    //                 let mut last_buf = Some(Vec::with_capacity(1 << (num_x - 1)));
    //                 for (eval, pat) in evals.iter_mut().zip(pattern.iter()) {
    //                     let mut evals = Cow::Borrowed(self.evals());
    //                     let skip = pat & skip_mask;
    //                     let x_0 = if pat.nth_bit(distance) {
    //                         &flipped_x[0]
    //                     } else {
    //                         &x[0]
    //                     };
    //                     merge_into(&mut buf, &evals, x_0, distance + 1, skip);
    //                     evals = mem::replace(&mut buf, last_buf.take().unwrap()).into();

    //                     for ((x_i, flipped_x_i), idx) in
    //                         x.iter().zip(flipped_x.iter()).zip(distance..).skip(1)
    //                     {
    //                         let x_i = if pat.nth_bit(idx) { flipped_x_i } else { x_i };
    //                         merge_in_place(&mut evals, x_i, 1, 0, &mut buf);
    //                     }
    //                     *eval = evals[0];
    //                     last_buf = Some(evals.into_owned());
    //                 }
    //             },
    //         );
    //     }
    //     evals
    // }
}

impl<F: Field, P: Borrow<MultilinearPolynomial<F>>> Add<P> for &MultilinearPolynomial<F> {
    type Output = MultilinearPolynomial<F>;

    fn add(self, rhs: P) -> MultilinearPolynomial<F> {
        let mut output = self.clone();
        output += rhs;
        output
    }
}

impl<F: Field, P: Borrow<MultilinearPolynomial<F>>> AddAssign<P> for MultilinearPolynomial<F> {
    fn add_assign(&mut self, rhs: P) {
        let rhs = rhs.borrow();
        match (self.is_empty(), rhs.is_empty()) {
            (_, true) => {}
            (true, false) => *self = rhs.clone(),
            (false, false) => {
                assert_eq!(self.num_vars, rhs.num_vars);

                parallelize(&mut self.evals, |(lhs, start)| {
                    for (lhs, rhs) in lhs.iter_mut().zip(rhs[start..].iter()) {
                        *lhs += rhs;
                    }
                });
            }
        }
    }
}

impl<F: Field, BF: Borrow<F>, P: Borrow<MultilinearPolynomial<F>>> AddAssign<(BF, P)>
    for MultilinearPolynomial<F>
{
    fn add_assign(&mut self, (scalar, rhs): (BF, P)) {
        let (scalar, rhs) = (scalar.borrow(), rhs.borrow());
        match (self.is_empty(), rhs.is_empty() | (scalar == &F::ZERO)) {
            (_, true) => {}
            (true, false) => {
                *self = rhs.clone();
                *self *= scalar;
            }
            (false, false) => {
                assert_eq!(self.num_vars, rhs.num_vars);

                if scalar == &F::ONE {
                    *self += rhs;
                } else if scalar == &-F::ONE {
                    *self -= rhs;
                } else {
                    parallelize(&mut self.evals, |(lhs, start)| {
                        for (lhs, rhs) in lhs.iter_mut().zip(rhs[start..].iter()) {
                            *lhs += &(*scalar * rhs);
                        }
                    });
                }
            }
        }
    }
}

impl<F: Field, P: Borrow<MultilinearPolynomial<F>>> Sub<P> for &MultilinearPolynomial<F> {
    type Output = MultilinearPolynomial<F>;

    fn sub(self, rhs: P) -> MultilinearPolynomial<F> {
        let mut output = self.clone();
        output -= rhs;
        output
    }
}

impl<F: Field, P: Borrow<MultilinearPolynomial<F>>> SubAssign<P> for MultilinearPolynomial<F> {
    fn sub_assign(&mut self, rhs: P) {
        let rhs = rhs.borrow();
        match (self.is_empty(), rhs.is_empty()) {
            (_, true) => {}
            (true, false) => {
                *self = rhs.clone();
                *self *= &-F::ONE;
            }
            (false, false) => {
                assert_eq!(self.num_vars, rhs.num_vars);

                parallelize(&mut self.evals, |(lhs, start)| {
                    for (lhs, rhs) in lhs.iter_mut().zip(rhs[start..].iter()) {
                        *lhs -= rhs;
                    }
                });
            }
        }
    }
}

impl<F: Field, BF: Borrow<F>, P: Borrow<MultilinearPolynomial<F>>> SubAssign<(BF, P)>
    for MultilinearPolynomial<F>
{
    fn sub_assign(&mut self, (scalar, rhs): (BF, P)) {
        *self += (-*scalar.borrow(), rhs);
    }
}

impl<F: Field, BF: Borrow<F>> Mul<BF> for &MultilinearPolynomial<F> {
    type Output = MultilinearPolynomial<F>;

    fn mul(self, rhs: BF) -> MultilinearPolynomial<F> {
        let mut output = self.clone();
        output *= rhs;
        output
    }
}

impl<F: Field, BF: Borrow<F>> MulAssign<BF> for MultilinearPolynomial<F> {
    fn mul_assign(&mut self, rhs: BF) {
        let rhs = rhs.borrow();
        if rhs == &F::ZERO {
            self.evals = vec![F::ZERO; self.evals.len()]
        } else if rhs == &-F::ONE {
            parallelize(&mut self.evals, |(evals, _)| {
                for eval in evals.iter_mut() {
                    *eval = -*eval;
                }
            });
        } else if rhs != &F::ONE {
            parallelize(&mut self.evals, |(lhs, _)| {
                for lhs in lhs.iter_mut() {
                    *lhs *= rhs;
                }
            });
        }
    }
}

impl<F: Field, P: Borrow<MultilinearPolynomial<F>>> Sum<P> for MultilinearPolynomial<F> {
    fn sum<I: Iterator<Item = P>>(mut iter: I) -> MultilinearPolynomial<F> {
        let init = match (iter.next(), iter.next()) {
            (Some(lhs), Some(rhs)) => lhs.borrow() + rhs.borrow(),
            (Some(lhs), None) => return lhs.borrow().clone(),
            _ => return Self::zero(),
        };
        iter.fold(init, |mut acc, poly| {
            acc += poly.borrow();
            acc
        })
    }
}

impl<F: Field, BF: Borrow<F>, P: Borrow<MultilinearPolynomial<F>>> Sum<(BF, P)>
    for MultilinearPolynomial<F>
{
    fn sum<I: Iterator<Item = (BF, P)>>(mut iter: I) -> MultilinearPolynomial<F> {
        let init = match iter.next() {
            Some((scalar, poly)) => {
                let mut poly = poly.borrow().clone();
                poly *= scalar.borrow();
                poly
            }
            _ => return Self::zero(),
        };
        iter.fold(init, |mut acc, (scalar, poly)| {
            acc += (scalar.borrow(), poly.borrow());
            acc
        })
    }
}

impl_index!(MultilinearPolynomial, evals);

pub(crate) fn evaluate<F: Field>(evals: &[F], x: &[F]) -> F {
    assert_eq!(1 << x.len(), evals.len());

    let mut evals = Cow::Borrowed(evals);
    let mut bits = Vec::new();
    let mut buf = Vec::with_capacity(evals.len() >> 1);
    for x_i in x.iter() {
        if x_i == &F::ZERO || x_i == &F::ONE {
            bits.push(x_i == &F::ONE);
            continue;
        }

        let distance = bits.len() + 1;
        let skip = usize_from_bits_le(&bits);
        merge_in_place(&mut evals, x_i, distance, skip, &mut buf);
        bits.clear();
    }

    evals[usize_from_bits_le(&bits)]
}

pub fn rotation_eval<F: Field>(x: &[F], rotation: Rotation, evals_for_rotation: &[F]) -> F {
    // if rotation == Rotation::cur() {
    // assert!(evals_for_rotation.len() == 1);
    return evals_for_rotation[0];
    // }

    // let num_vars = x.len();
    // let distance = rotation.distance();
    // assert!(evals_for_rotation.len() == 1 << distance);
    // assert!(distance <= num_vars);

    // let (pattern, nths, x) = if rotation < Rotation::cur() {
    //     (
    //         rotation_eval_coeff_pattern::<false>(num_vars, distance),
    //         (1..=distance).rev().collect_vec(),
    //         x[0..distance].iter().rev().collect_vec(),
    //     )
    // } else {
    //     (
    //         rotation_eval_coeff_pattern::<true>(num_vars, distance),
    //         (num_vars - 1..).take(distance).collect(),
    //         x[num_vars - distance..].iter().collect(),
    //     )
    // };
    // x.into_iter().zip(nths).enumerate().fold(
    //     Cow::Borrowed(evals_for_rotation),
    //     |evals, (idx, (x_i, nth))| {
    //         pattern
    //             .iter()
    //             .step_by(1 << idx)
    //             .map(|pat| pat.nth_bit(nth))
    //             .zip(zip_self!(evals.iter()))
    //             .map(|(bit, (eval_0, eval_1))| {
    //                 if bit {
    //                     (*eval_0 - eval_1) * x_i + eval_1
    //                 } else {
    //                     (*eval_1 - eval_0) * x_i + eval_0
    //                 }
    //             })
    //             .collect_vec()
    //             .into()
    //     },
    // )[0]
}

pub fn rotation_eval_points<F: Field>(x: &[F], rotation: Rotation) -> Vec<Vec<F>> {
    if rotation == Rotation::cur() {
        return vec![x.to_vec()];
    }

    let distance = rotation.distance();
    let num_x = x.len() - distance;
    if rotation < Rotation::cur() {
        let pattern = rotation_eval_point_pattern::<false>(x.len(), distance);
        let x = &x[distance..];
        let flipped_x = x.iter().map(flip).collect_vec();
        pattern
            .iter()
            .map(|pat| {
                chain![
                    (0..num_x).map(|idx| {
                        if pat.nth_bit(idx) {
                            flipped_x[idx]
                        } else {
                            x[idx]
                        }
                    }),
                    (0..distance).map(|idx| bit_to_field(pat.nth_bit(idx + num_x)))
                ]
                .collect_vec()
            })
            .collect()
    } else {
        let pattern = rotation_eval_point_pattern::<true>(x.len(), distance);
        let x = &x[..num_x];
        let flipped_x = x.iter().map(flip).collect_vec();
        pattern
            .iter()
            .map(|pat| {
                chain![
                    (0..distance).map(|idx| bit_to_field(pat.nth_bit(idx))),
                    (0..num_x).map(|idx| {
                        if pat.nth_bit(idx + distance) {
                            flipped_x[idx]
                        } else {
                            x[idx]
                        }
                    })
                ]
                .collect_vec()
            })
            .collect()
    }
}

pub(crate) fn rotation_eval_point_pattern<const NEXT: bool>(
    num_vars: usize,
    distance: usize,
) -> Vec<usize> {
    let bf = BinaryField::new(num_vars);
    let remainder = if NEXT { bf.primitive() } else { bf.x_inv() };
    let mut pattern = vec![0; 1 << distance];
    for depth in 0..distance {
        for (e, o) in zip_self!(0..pattern.len(), 1 << (distance - depth)) {
            let rotated = if NEXT {
                pattern[e] << 1
            } else {
                pattern[e] >> 1
            };
            pattern[o] = rotated ^ remainder;
            pattern[e] = rotated;
        }
    }
    pattern
}

pub(crate) fn rotation_eval_coeff_pattern<const NEXT: bool>(
    num_vars: usize,
    distance: usize,
) -> Vec<usize> {
    let bf = BinaryField::new(num_vars);
    let remainder = if NEXT {
        bf.primitive() - (1 << num_vars)
    } else {
        bf.x_inv() << distance
    };
    let mut pattern = vec![0; 1 << (distance - 1)];
    for depth in 0..distance - 1 {
        for (e, o) in zip_self!(0..pattern.len(), 1 << (distance - depth - 1)) {
            let rotated = if NEXT {
                pattern[e] << 1
            } else {
                pattern[e] >> 1
            };
            pattern[o] = rotated ^ remainder;
            pattern[e] = rotated;
        }
    }
    pattern
}

fn flip<F: Field>(x: &F) -> F {
    F::ONE - x
}

fn bit_to_field<F: Field>(bit: bool) -> F {
    if bit {
        F::ONE
    } else {
        F::ZERO
    }
}

fn merge_in_place<F: Field>(
    evals: &mut Cow<[F]>,
    x_i: &F,
    distance: usize,
    skip: usize,
    buf: &mut Vec<F>,
) {
    merge_into(buf, evals, x_i, distance, skip);
    if let Cow::Owned(_) = evals {
        mem::swap(evals.to_mut(), buf);
    } else {
        *evals = mem::replace(buf, Vec::with_capacity(buf.len() >> 1)).into();
    }
}

pub(crate) fn merge_into<F: Field>(
    target: &mut Vec<F>,
    evals: &[F],
    x_i: &F,
    distance: usize,
    skip: usize,
) {
    assert!(target.capacity() >= evals.len() >> distance);
    target.resize(evals.len() >> distance, F::ZERO);

    let step = 1 << distance;
    parallelize(target, |(target, start)| {
        let start = (start << distance) + skip;
        for (target, (eval_0, eval_1)) in
            target.iter_mut().zip(zip_self!(evals.iter(), step, start))
        {
            *target = (*eval_1 - eval_0) * x_i + eval_0;
        }
    });
}

macro_rules! zip_self {
    (@ $iter:expr, $step:expr, $skip:expr) => {
        $iter.skip($skip).step_by($step).zip($iter.skip($skip + ($step >> 1)).step_by($step))
    };
    ($iter:expr) => {
        zip_self!(@ $iter, 2, 0)
    };
    ($iter:expr, $step:expr) => {
        zip_self!(@ $iter, $step, 0)
    };
    ($iter:expr, $step:expr, $skip:expr) => {
        zip_self!(@ $iter, $step, $skip)
    };
}

pub(crate) use zip_self;

#[cfg(test)]
mod test {
    use crate::{
        poly::multilinear::{rotation_eval, zip_self, MultilinearPolynomial},
        util::{
            arithmetic::Field,
            expression::{
                rotate::{BinaryField, Rotatable},
                Rotation,
            },
            test::rand_vec,
            Itertools,
        },
    };
    use halo2_curves::bn256::Fr;
    use rand::{rngs::OsRng, RngCore};
    use std::iter;

    fn fix_vars<F: Field>(evals: &[F], x: &[F]) -> Vec<F> {
        x.iter().fold(evals.to_vec(), |evals, x_i| {
            zip_self!(evals.iter())
                .map(|(eval_0, eval_1)| (*eval_1 - eval_0) * x_i + eval_0)
                .collect_vec()
        })
    }

    #[test]
    fn fix_var() {
        let rand_x_i = || match OsRng.next_u32() % 3 {
            0 => Fr::ZERO,
            1 => Fr::ONE,
            2 => Fr::random(OsRng),
            _ => unreachable!(),
        };
        for num_vars in 0..16 {
            for _ in 0..10 {
                let poly = MultilinearPolynomial::rand(num_vars, OsRng);
                let x = iter::repeat_with(rand_x_i).take(num_vars).collect_vec();
                let eval = fix_vars(poly.evals(), &x)[0];
                assert_eq!(poly.evaluate(&x), eval);
                assert_eq!(x.iter().fold(poly, |poly, x_i| poly.fix_var(x_i))[0], eval);
            }
        }
    }

    #[test]
    fn evaluate_for_rotation() {
        let mut rng = OsRng;
        for num_vars in 1..16 {
            let bf = BinaryField::new(num_vars);
            let rotate = |f: &Vec<Fr>| {
                (0..1 << num_vars)
                    .map(|idx| f[bf.rotate(idx, Rotation::next())])
                    .collect_vec()
            };
            let f = rand_vec(1 << num_vars, &mut rng);
            let fs = iter::successors(Some(f), |f| Some(rotate(f)))
                .map(MultilinearPolynomial::new)
                .take(num_vars)
                .collect_vec();
            let x = rand_vec::<Fr>(num_vars, &mut rng);

            for rotation in -(num_vars as i32) + 1..num_vars as i32 {
                let rotation = Rotation(rotation);
                let (f, f_rotated) = if rotation < Rotation::cur() {
                    (fs.last().unwrap(), &fs[fs.len() - rotation.distance() - 1])
                } else {
                    (fs.first().unwrap(), &fs[rotation.distance()])
                };
                assert_eq!(
                    rotation_eval(&x, rotation, &f.evaluate_for_rotation(&x, rotation)),
                    f_rotated.evaluate(&x),
                );
            }
        }
    }
}
