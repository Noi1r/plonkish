use crate::{
    pcs::{
        multilinear::{additive, quotients},
        univariate::{
            err_too_large_deree, UnivariateKzg, UnivariateKzgProverParam,
            UnivariateKzgVerifierParam,
        },
        Evaluation, Evaluation_for_shift, Point, PolynomialCommitmentScheme,
    },
    poly::{multilinear::{rotation_eval, MultilinearPolynomial}, univariate::UnivariatePolynomial},
    util::{
        arithmetic::{
            powers, squares, variable_base_msm, BatchInvert, Curve, Field, MultiMillerLoop,
        },
        chain, izip,
        parallel::parallelize,
        transcript::{TranscriptRead, TranscriptWrite},
        Deserialize, DeserializeOwned, Itertools, Serialize,
    },
    Error,
};
use rand::RngCore;
use std::{marker::PhantomData, process::exit};

#[derive(Clone, Debug)]
pub struct Zeromorph<Pcs>(PhantomData<Pcs>);

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(bound(
    serialize = "M::G1Affine: Serialize",
    deserialize = "M::G1Affine: DeserializeOwned",
))]
pub struct ZeromorphKzgProverParam<M: MultiMillerLoop> {
    commit_pp: UnivariateKzgProverParam<M>,
    open_pp: UnivariateKzgProverParam<M>,
}

impl<M: MultiMillerLoop> ZeromorphKzgProverParam<M> {
    pub fn degree(&self) -> usize {
        self.commit_pp.degree()
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(bound(
    serialize = "M::G1Affine: Serialize, M::G2Affine: Serialize",
    deserialize = "M::G1Affine: DeserializeOwned, M::G2Affine: DeserializeOwned",
))]
pub struct ZeromorphKzgVerifierParam<M: MultiMillerLoop> {
    vp: UnivariateKzgVerifierParam<M>,
    s_offset_g2: M::G2Affine,
}

impl<M: MultiMillerLoop> ZeromorphKzgVerifierParam<M> {
    pub fn g1(&self) -> M::G1Affine {
        self.vp.g1()
    }

    pub fn g2(&self) -> M::G2Affine {
        self.vp.g2()
    }

    pub fn s_g2(&self) -> M::G2Affine {
        self.vp.s_g2()
    }
}

impl<M> PolynomialCommitmentScheme<M::Scalar> for Zeromorph<UnivariateKzg<M>>
where
    M: MultiMillerLoop,
    M::Scalar: Serialize + DeserializeOwned,
    M::G1Affine: Serialize + DeserializeOwned,
    M::G2Affine: Serialize + DeserializeOwned,
{
    type Param = <UnivariateKzg<M> as PolynomialCommitmentScheme<M::Scalar>>::Param;
    type ProverParam = ZeromorphKzgProverParam<M>;
    type VerifierParam = ZeromorphKzgVerifierParam<M>;
    // 多项式的形式还是MultilinearPolynomial
    type Polynomial = MultilinearPolynomial<M::Scalar>;
    type Commitment = <UnivariateKzg<M> as PolynomialCommitmentScheme<M::Scalar>>::Commitment;
    type CommitmentChunk =
        <UnivariateKzg<M> as PolynomialCommitmentScheme<M::Scalar>>::CommitmentChunk;

    fn setup(poly_size: usize, batch_size: usize, rng: impl RngCore) -> Result<Self::Param, Error> {
        assert!(poly_size.is_power_of_two());
        UnivariateKzg::<M>::setup((poly_size+1).next_power_of_two(), batch_size, rng)
    }

    fn trim(
        param: &Self::Param,
        poly_size: usize,
        batch_size: usize,
    ) -> Result<(Self::ProverParam, Self::VerifierParam), Error> {
        assert!(poly_size.is_power_of_two());

        let (commit_pp, vp) = UnivariateKzg::<M>::trim(param, (poly_size+1).next_power_of_two(), batch_size)?;
        let offset = param.monomial_g1().len() - (poly_size+1).next_power_of_two();
        let open_pp = {
            let monomial_g1 = param.monomial_g1()[offset..].to_vec();
            UnivariateKzgProverParam::new((poly_size.ilog2()+1) as usize, monomial_g1, Vec::new())
        };
        let s_offset_g2 = param.powers_of_s_g2()[offset];

        Ok((
            ZeromorphKzgProverParam { commit_pp, open_pp },
            ZeromorphKzgVerifierParam { vp, s_offset_g2 },
        ))
    }

    fn commit(pp: &Self::ProverParam, poly: &Self::Polynomial) -> Result<Self::Commitment, Error> {
        if pp.degree() + 1 < poly.evals().len() {
            let got = poly.evals().len() - 1;
            return Err(err_too_large_deree("commit", pp.degree(), got));
        }

        Ok(UnivariateKzg::commit_monomial(&pp.commit_pp, poly.evals()))
    }

    fn batch_commit<'a>(
        pp: &Self::ProverParam,
        polys: impl IntoIterator<Item = &'a Self::Polynomial>,
    ) -> Result<Vec<Self::Commitment>, Error> {
        polys
            .into_iter()
            .map(|poly| Self::commit(pp, poly))
            .collect()
    }

    fn open(
        pp: &Self::ProverParam,
        poly: &Self::Polynomial,
        comm: &Self::Commitment,
        point: &Point<M::Scalar, Self::Polynomial>,
        eval: &M::Scalar,
        transcript: &mut impl TranscriptWrite<Self::CommitmentChunk, M::Scalar>,
    ) -> Result<(), Error> {
        // 多变量多项式返回 num_vars >=2
        let num_vars = poly.num_vars();
        // 多项式degree 需要小于等于commit_pp的degree
        if pp.degree() + 1 < poly.evals().len() {
            let got = poly.evals().len() - 1;
            return Err(err_too_large_deree("open", pp.degree(), got));
        }

        if cfg!(feature = "sanity-check") {
            assert_eq!(Self::commit(pp, poly).unwrap().0, comm.0);
            assert_eq!(poly.evaluate(point), *eval);
            //println!("===============================");
        }

        // point.len() == poly.num_vars()
        let (quotients, remainder) =
            quotients(poly, point, |_, q| UnivariatePolynomial::monomial(q));
        UnivariateKzg::batch_commit_and_write(&pp.commit_pp, &quotients, transcript)?;

        if cfg!(feature = "sanity-check") {
            assert_eq!(&remainder, eval);
        }

        let y = transcript.squeeze_challenge();

        let q_hat = {
            let mut q_hat = vec![M::Scalar::ZERO; 1 << num_vars];
            for (idx, (power_of_y, q)) in izip!(powers(y), &quotients).enumerate() {
                let offset = (1 << num_vars) - (1 << idx);
                parallelize(&mut q_hat[offset..], |(q_hat, start)| {
                    izip!(q_hat, q.iter().skip(start))
                        .for_each(|(q_hat, q)| *q_hat += power_of_y * q)
                });
            }
            UnivariatePolynomial::monomial(q_hat)
        };
        // println!("UnivariateKzg::commit_and_write(&pp.commit_pp, &q_hat, transcript)?;");
        UnivariateKzg::commit_and_write(&pp.commit_pp, &q_hat, transcript)?;

        let x = transcript.squeeze_challenge();
        let z = transcript.squeeze_challenge();

        let (eval_scalar, q_scalars) = eval_and_quotient_scalars(y, x, z, point);

        let mut f = UnivariatePolynomial::monomial(poly.evals().to_vec());
        f *= z;
        f += &q_hat;
        f[0] += eval_scalar * eval;
        izip!(&quotients, &q_scalars).for_each(|(q, scalar)| f += (scalar, q));

        // UnivariateKzg::commit_and_write(&pp.open_pp, &f, transcript)?;
        
        //println!("f: {:?}", f.coeffs());
        //println!("f[0]: {:?}", f[0]);
        //println!("f.evaluate(): {:?}", f.evaluate(&x));
        //println!("=======================================================");

        let comm = if cfg!(feature = "sanity-check") {
            assert_eq!(f.evaluate(&x), M::Scalar::ZERO);
            UnivariateKzg::commit_monomial(&pp.open_pp, f.coeffs())
        } else {
            Default::default()
        };

        UnivariateKzg::<M>::open(&pp.open_pp, &f, &comm, &x, &M::Scalar::ZERO, transcript)
    }

    fn batch_open<'a>(
        pp: &Self::ProverParam,
        polys: impl IntoIterator<Item = &'a Self::Polynomial>,
        comms: impl IntoIterator<Item = &'a Self::Commitment>,
        points: &[Point<M::Scalar, Self::Polynomial>],
        evals: &[Evaluation<M::Scalar>],
        transcript: &mut impl TranscriptWrite<Self::CommitmentChunk, M::Scalar>,
    ) -> Result<(), Error>
    where
        Self::Commitment: 'a,
    {
        let polys = polys.into_iter().collect_vec();
        let comms = comms.into_iter().collect_vec();
        let num_vars = points.first().map(|point| point.len()).unwrap_or_default();
        additive::batch_open::<_, Self>(pp, num_vars, polys, comms, points, evals, transcript)
    }

    fn batch_open_for_shift<'a>(
        pp: &Self::ProverParam,
        polys: impl IntoIterator<Item = &'a Self::Polynomial>,
        comms: impl IntoIterator<Item = &'a Self::Commitment>,
        points: &[Point<M::Scalar, Self::Polynomial>],
        evals: &[Evaluation_for_shift<M::Scalar>],
        transcript: &mut impl TranscriptWrite<Self::CommitmentChunk, M::Scalar>,
    ) -> Result<(), Error>
    where
        Self::Commitment: 'a,
    {
        let polys = polys.into_iter().collect_vec();
        let comms = comms.into_iter().collect_vec();
        let num_vars = points.first().map(|point| point.len()).unwrap_or_default();
        additive::batch_open_for_shift::<_, Self>(
            pp, num_vars, polys, comms, points, evals, transcript,
        )
    }

    fn read_commitments(
        vp: &Self::VerifierParam,
        num_polys: usize,
        transcript: &mut impl TranscriptRead<Self::CommitmentChunk, M::Scalar>,
    ) -> Result<Vec<Self::Commitment>, Error> {
        UnivariateKzg::read_commitments(&vp.vp, num_polys, transcript)
    }

    fn prove_shifted_evaluation(
        pp: &Self::ProverParam,                       // ZeromorphKzgProverParam<M>
        poly: &Self::Polynomial, // MultilinearPolynomial<M::Scalar> (merged and scaled)
        comm: &Self::Commitment, // Commitment<M::Scalar, UnivariateKzg<M>> (to the original merged poly)
        point: &Point<M::Scalar, Self::Polynomial>, // Vec<M::Scalar> (the point 'u')
        value: &M::Scalar,       // Claimed value v = f_shifted(u)
        rotation: &crate::util::expression::Rotation, // Use the provided Rotation struct
        transcript: &mut impl TranscriptWrite<Self::CommitmentChunk, M::Scalar>,
    ) -> Result<(), Error> {
        // Use the i32 value to check the sign, and distance() for magnitude
        let signed_d = rotation.0;
        let abs_d = rotation.distance(); // Get magnitude (usize)

        let num_vars = poly.num_vars();
        let n_evals = 1 << num_vars; // 2^num_vars

        // --- Pre-checks ---

        if cfg!(feature = "sanity-check") {
            assert_eq!(rotation_eval(point, *rotation, &poly.evaluate_for_rotation(point, *rotation)), *value);
        }

        if signed_d == 0 {
            return Err(Error::InternalError(
                "Rotation distance is zero in prove_shifted_evaluation".to_string(),
            ));
        }
        if pp.commit_pp.degree() + 1 < n_evals {
            let got = n_evals - 1;
            return Err(err_too_large_deree(
                "prove_shifted_evaluation (poly degree)",
                pp.commit_pp.degree(),
                got,
            ));
        }
        if abs_d >= n_evals {
            return Err(Error::InvalidInput(format!(
                "Rotation distance {} magnitude must be smaller than polynomial size {}",
                signed_d, n_evals
            )));
        }

        // --- Shift the polynomial f to get f_d ---
        let mut poly_d_evals = poly.evals().to_vec();
        if signed_d > 0 {
            poly_d_evals.rotate_left(abs_d);
        } else {
            poly_d_evals.rotate_right(abs_d);
        }
        let poly_d = MultilinearPolynomial::new(poly_d_evals.clone());

        // --- Compute and Commit Multilinear Quotients q_{d,k} ---
        let (quotients_d, remainder_d): (
            Vec<crate::poly::univariate::UnivariatePolynomial<M::Scalar>>,
            M::Scalar,
        ) = quotients(&poly_d, point, |_, q| UnivariatePolynomial::monomial(q));

        UnivariateKzg::batch_commit_and_write(&pp.commit_pp, &quotients_d, transcript)?;

        println!("remainder_d: {:?}", remainder_d);
        println!("value: {:?}", value);

        if cfg!(feature = "sanity-check") {
            let value_d = poly_d.evaluate(point);
            assert_eq!(
                &remainder_d, &value_d,
                "Shifted polynomial evaluation mismatch inside prover (f_d(u) != v)"
            );
        }

        // --- Combine Quotients with challenge y ---
        let y = transcript.squeeze_challenge();
        let q_d_hat = {
            let mut q_d_hat_coeffs = vec![M::Scalar::ZERO; n_evals];
            for (idx, (power_of_y, q_k)) in izip!(powers(y), &quotients_d).enumerate() {
                let offset = n_evals - (1 << idx);
                let q_k_coeffs = q_k.coeffs();
                parallelize(&mut q_d_hat_coeffs[offset..], |(q_hat_chunk, start)| {
                    izip!(q_hat_chunk, q_k_coeffs.iter().skip(start))
                        .for_each(|(q_hat_val, q_k_coeff)| *q_hat_val += power_of_y * q_k_coeff)
                });
            }
            UnivariatePolynomial::monomial(q_d_hat_coeffs)
        };

        crate::pcs::univariate::UnivariateKzg::<M>::commit_and_write(
            &pp.commit_pp,
            &q_d_hat,
            transcript,
        )?;

        // --- Compute Scalars and Construct Final Check Polynomial F_d ---
        let x = transcript.squeeze_challenge();
        let z = transcript.squeeze_challenge();

        let (eval_scalar, q_scalars) = eval_and_quotient_scalars(y, x, z, point);

         // 原始多项式 f 的单变量表示 (U_n(f))
         let f_uni = crate::poly::univariate::UnivariatePolynomial::monomial(poly.evals().to_vec());

         // 构造多项式 X^N
         let x_n = {
             let mut coeffs = vec![M::Scalar::ZERO; n_evals + 1];
             if let Some(coeff) = coeffs.get_mut(n_evals) {
                  *coeff = M::Scalar::ONE;
             } else if n_evals == 0 { // 处理 num_vars=0 的边界情况
                  coeffs = vec![M::Scalar::ONE]; // X^0 = 1
             }
              crate::poly::univariate::UnivariatePolynomial::monomial(coeffs)
         };

 
         // 构造检查多项式中涉及 v 和 q_{d,k} 的部分
         // term3_inner = v * eval_scalar + sum(q_{d,k} * q_scalar_k)
         let value_d = poly_d.evaluate(point);
         let mut term3_inner =
              // 使用 monomial 构造常数项
             crate::poly::univariate::UnivariatePolynomial::monomial(vec![value_d * eval_scalar]);
         izip!(&quotients_d, &q_scalars).for_each(|(q, scalar)| {
              term3_inner += &(q * *scalar)
         });
 
 
         // 根据移位方向计算 F_d
         let F_d = if signed_d > 0 {
             // 左移 d = abs_d
             // P_Ad: 使用 f 的前 d 个系数构造
             let p_ad =
                 crate::poly::univariate::UnivariatePolynomial::monomial(poly.evals()[..abs_d].to_vec());
             // X^d
             let x_d = {
                 let mut coeffs = vec![M::Scalar::ZERO; abs_d + 1];
                  if let Some(coeff) = coeffs.get_mut(abs_d) { *coeff = M::Scalar::ONE; }
                  else if abs_d == 0 { coeffs = vec![M::Scalar::ONE]; } // X^0 = 1
                 crate::poly::univariate::UnivariatePolynomial::monomial(coeffs)
             };
 
             // Term1 = f_uni - P_Ad + P_Ad * X^N
             let mut term1 = f_uni; // 克隆 f_uni
             term1 -= &p_ad; 
             term1 += &(&p_ad * &x_n); 
 
             // F_d = z * Term1 + X^d * q_d_hat + X^d * term3_inner
             let mut F_d = &term1 * z;
             F_d += &(&x_d * &q_d_hat);      
             F_d += &(&x_d * &term3_inner); 
             F_d
         } else {
             // 右移 d' = abs_d = -signed_d
             // P_Ad: 使用 f 的后 d' 个系数构造
             let p_ad = crate::poly::univariate::UnivariatePolynomial::monomial(
                 poly.evals()[n_evals - abs_d..].to_vec(),
             );
             // X^{d'}
             let x_d_prime = {
                  let mut coeffs = vec![M::Scalar::ZERO; abs_d + 1];
                  if let Some(coeff) = coeffs.get_mut(abs_d) { *coeff = M::Scalar::ONE; }
                  else if abs_d == 0 { coeffs = vec![M::Scalar::ONE]; } // X^0 = 1
                 crate::poly::univariate::UnivariatePolynomial::monomial(coeffs)
             };
 
             // Term1 = X^d' * f_uni - X^d' * P_Ad + P_Ad * X^N 
             let mut term1 = &x_d_prime * &f_uni;    
             term1 -= &(&x_d_prime * &p_ad); 
             term1 += &(&p_ad * &x_n);        
 
             // F_d = z * Term1 + q_d_hat + term3_inner 
             let mut F_d = &term1 * z;
             F_d += &q_d_hat;           
             F_d += &term3_inner;       
             F_d
         };
         
        // let comm = UnivariateKzg::commit_and_write(&pp.commit_pp, &F_d, transcript)?;

        let comm = if cfg!(feature = "sanity-check") {
            assert_eq!(F_d.evaluate(&x), M::Scalar::ZERO);
            UnivariateKzg::commit_monomial(&pp.open_pp, F_d.coeffs())
        } else {
            Default::default()
        };

        UnivariateKzg::<M>::open(
            &pp.open_pp,
            &F_d,
            &comm,
            &x,
            &M::Scalar::ZERO, // 声称的求值结果是零
            transcript,
        )
    } // End of prove_shifted_evaluation

    fn verify(
        vp: &Self::VerifierParam,
        comm: &Self::Commitment,
        point: &Point<M::Scalar, Self::Polynomial>,
        eval: &M::Scalar,
        transcript: &mut impl TranscriptRead<Self::CommitmentChunk, M::Scalar>,
    ) -> Result<(), Error> {
        let num_vars = point.len();

        let q_comms = transcript.read_commitments(num_vars)?;

        let y = transcript.squeeze_challenge();

        let q_hat_comm = transcript.read_commitment()?;

        let x = transcript.squeeze_challenge();
        let z = transcript.squeeze_challenge();

        let (eval_scalar, q_scalars) = eval_and_quotient_scalars(y, x, z, point);

        let scalars = chain![[M::Scalar::ONE, z, eval_scalar * eval], q_scalars].collect_vec();
        let bases = chain![[q_hat_comm, comm.0, vp.g1()], q_comms].collect_vec();
        let c: M::G1Affine = variable_base_msm(&scalars, &bases).into();

        // let c = transcript.read_commitment()?;

        let pi = transcript.read_commitment()?;

        M::pairings_product_is_identity(&[
            (&c, &(-vp.s_offset_g2).into()),
            (&pi, &(vp.s_g2() - (vp.g2() * x).into()).to_affine().into()),
        ])
        .then_some(())
        .ok_or_else(|| Error::InvalidPcsOpen("Invalid Zeromorph KZG open".to_string()))
    }



    fn verify_shifted_evaluation(
        vp: &Self::VerifierParam,                       // ZeromorphKzgVerifierParam<M>
        comm: &Self::Commitment, // 对原始合并多项式 f 的承诺 C_f
        point: &Point<M::Scalar, Self::Polynomial>, // 求值点 u
        value: &M::Scalar,       // 声称的求值 v = f_d(u)
        rotation: &crate::util::expression::Rotation, 
        transcript: &mut impl TranscriptRead<Self::CommitmentChunk, M::Scalar>,
    ) -> Result<(), Error>
    {
        let num_vars = point.len();

        // 1. 从 transcript 读取商多项式的承诺 C_{q_{d,k}}
        let q_comms_d = crate::pcs::univariate::UnivariateKzg::<M>::read_commitments(&vp.vp, num_vars, transcript)?;

        // 2. Squeeze challenge y, 读取合并商多项式的承诺 C_{q_d_hat}
        let y = transcript.squeeze_challenge();
        let q_d_hat_comm = crate::pcs::univariate::UnivariateKzg::<M>::read_commitment(&vp.vp, transcript)?;

        // 3. Squeeze challenges x, z, 计算验证所需的标量
        let x = transcript.squeeze_challenge();
        let z = transcript.squeeze_challenge();
        let (eval_scalar, q_scalars) = eval_and_quotient_scalars(y, x, z, point); 

        // println!("eval_scalar: {:?}", eval_scalar);
        // println!("q_scalars: {:?}", q_scalars);

        //    C 对应于 f_check = q_d_hat + z*f + eval_scalar*v*1 + sum(q_scalar_k * q_{d,k})
        let scalars = chain![
            [M::Scalar::ONE, z, eval_scalar * value],
            q_scalars
        ].collect_vec();
        let bases = chain![
            [q_d_hat_comm.0, comm.0, vp.g1()], // Use .0 for comms, vp.g1() is already G1Affine
            q_comms_d.iter().map(|c| c.0)      // Use .0 for quotient comms
        ].collect_vec();
        
        // 使用 MSM 计算重构的承诺 C
        let reconstructed_commitment_c: M::G1Affine = variable_base_msm(&scalars, &bases).into();


        // 5. 从 transcript 读取最终的单变量 KZG 打开证明 pi_d (由 prover 的 UnivariateKzg::open 生成)
        let pi_d = transcript.read_commitment()?;

        // 6. 执行最终的配对检查 (验证 C 在点 x 打开为 0，使用证明 pi_d)
        // if !M::pairings_product_is_identity(&[
        //     (&reconstructed_commitment_c, &(-vp.s_offset_g2).into()),
        //     (&pi_d, &(vp.s_g2() - (vp.g2() * x).into()).to_affine().into()),
        // ]) {
        //     return Err(Error::InvalidPcsOpen(format!(
        //         "Invalid Zeromorph KZG shifted open for rotation {}", rotation.0
        //     )));
        // }
        // Ok(())

        M::pairings_product_is_identity(&[
            (&reconstructed_commitment_c, &(-vp.s_offset_g2).into()),
            (&pi_d, &(vp.s_g2() - (vp.g2() * x).into()).to_affine().into()),
        ])
        .then_some(())
        .ok_or_else(|| {
            Error::InvalidPcsOpen(format!(
                "Invalid Zeromorph KZG shifted open for rotation {}", rotation.0))} // 使用 rotation.0 获取带符号距
        )
        // Ok(())
    } 


    fn batch_verify<'a>(
        vp: &Self::VerifierParam,
        comms: impl IntoIterator<Item = &'a Self::Commitment>,
        points: &[Point<M::Scalar, Self::Polynomial>],
        evals: &[Evaluation<M::Scalar>],
        transcript: &mut impl TranscriptRead<Self::CommitmentChunk, M::Scalar>,
    ) -> Result<(), Error> {
        let num_vars = points.first().map(|point| point.len()).unwrap_or_default();
        let comms = comms.into_iter().collect_vec();
        additive::batch_verify::<_, Self>(vp, num_vars, comms, points, evals, transcript)
    }

    fn batch_verify_for_shift<'a>(
        vp: &Self::VerifierParam,
        comms: impl IntoIterator<Item = &'a Self::Commitment>,
        points: &[Point<M::Scalar, Self::Polynomial>],
        evals: &[Evaluation_for_shift<M::Scalar>],
        transcript: &mut impl TranscriptRead<Self::CommitmentChunk, M::Scalar>,
    ) -> Result<(), Error> {
        let num_vars = points.first().map(|point| point.len()).unwrap_or_default();
        let comms = comms.into_iter().collect_vec();
        additive::batch_verify_for_shift::<_, Self>(vp, num_vars, comms, points, evals, transcript)
    }
}

fn eval_and_quotient_scalars<F: Field>(y: F, x: F, z: F, u: &[F]) -> (F, Vec<F>) {
    let num_vars = u.len();

    let squares_of_x = squares(x).take(num_vars + 1).collect_vec();
    //offsets＿of＿x ：计算一个偏移量序列 $O=\left[O_0, \ldots, O_{n-1}\right]$ ，其中 $O_k=\prod_{j=k+1}^{n-1} s_j=$ $\prod_{j=k+1}^{n-1} x^{2^j}$ 。这个序列与多线性求值和转换有关。
    let offsets_of_x = {
        let mut offsets_of_x = squares_of_x
            .iter()
            .rev()
            .skip(1)
            .scan(F::ONE, |state, power_of_x| {
                *state *= power_of_x;
                Some(*state)
            })
            .collect_vec();
        offsets_of_x.reverse();
        offsets_of_x
    };
    // vs ：计算另一个序列 $V=\left[V_0, \ldots, V_n\right]$ ，其中 $V_k=\left(s_n-1\right) /\left(s_k-1\right)=\left(x^{2^n}-\right.$ 1）$/\left(x^{2^k}-1\right)$ 。可以证明 $V_k=\sum_{j=0}^{2^{n-k}-1}\left(x^{2^k}\right)^j=\Phi_{n-k}\left(x^{2^k}\right)$（参见论文 中 $\Phi$ 的定义）。
    let vs = {
        let v_numer = squares_of_x[num_vars] - F::ONE;
        let mut v_denoms = squares_of_x
            .iter()
            .map(|square_of_x| *square_of_x - F::ONE)
            .collect_vec();
        v_denoms.batch_invert();
        v_denoms
            .iter()
            .map(|v_denom| v_numer * v_denom)
            .collect_vec()
    };


    let q_scalars = izip!(powers(y), offsets_of_x, squares_of_x, &vs, &vs[1..], u)
        .map(|(power_of_y, offset_of_x, square_of_x, v_i, v_j, u_i)| {
            -(power_of_y * offset_of_x + z * (square_of_x * v_j - *u_i * v_i))
        })
        .collect_vec();

    (-vs[0] * z, q_scalars)
}


#[cfg(test)]
mod test {
    use crate::{
        pcs::{
            multilinear::zeromorph::Zeromorph,
            test::{run_batch_commit_open_verify, run_commit_open_verify},
            univariate::UnivariateKzg,
        },
        util::transcript::Keccak256Transcript,
    };
    use halo2_curves::bn256::Bn256;

    type Pcs = Zeromorph<UnivariateKzg<Bn256>>;

    #[test]
    fn commit_open_verify() {
        run_commit_open_verify::<_, Pcs, Keccak256Transcript<_>>();
    }

    #[test]
    fn batch_commit_open_verify() {
        run_batch_commit_open_verify::<_, Pcs, Keccak256Transcript<_>>();
    }
}
