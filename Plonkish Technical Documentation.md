# Plonkish Technical Documentation

## 1. Introduction

Plonkish is a high-performance, Rust-based implementation of a zero-knowledge proof system based on the **HyperPlonk** proving system. Its core innovations include the integration of the **Zeromorph** polynomial commitment scheme (PCS) and a novel optimization for `shift` operations. The project aims to address the performance bottlenecks found in traditional Plonk-like systems when dealing with large-scale circuits, particularly the overhead associated with Fast Fourier Transforms (FFTs) and proving zero-knowledge properties.

By adapting the Plonk protocol from the traditional univariate polynomial setting to multilinear polynomials over the boolean hypercube, Plonkish achieves several key advantages:

* **Linear-Time Prover**: It eliminates the need for large FFTs, reducing the prover's runtime complexity from $O(N \log N)$ to $O(N)$, where $N$ is the circuit size.
* **Efficient High-Degree Custom Gates**: It supports custom gates of a much higher degree than standard Plonk without a significant impact on prover time. This allows for more efficient representations of complex circuits, such as cryptographic hash functions.
* **Optimized Shift Operations**: Leveraging the Zeromorph PCS, Plonkish dramatically reduces the cost of proving polynomial shift relations (e.g., `f(x) = g(x_shifted)`). This is critical for efficiently implementing lookup arguments and state machines.

The codebase provides a complete backend for HyperPlonk, and also integrates the cutting-edge **Anemoi** hash function to build efficient Merkle trees, enhancing its utility for applications requiring verifiable state, such as zk-EVMs.

## 2. System Architecture

Plonkish is designed with a modular architecture, decoupling the different layers and functionalities of the proving system to enhance extensibility and maintainability. The core architecture is divided into the following components:

* **Backend Layer** (`plonkish_backend/src/backend/`):
    * This is the heart of the proving system. It defines the `PlonkishBackend` trait, which provides a unified interface for different proof systems.
    * The primary implementation is the `HyperPlonk` module, containing the complete logic for the preprocessor, prover, and verifier.
    * This layer explicitly supports efficient shift proofs through methods like `prove_with_shift()` and `verify_with_shift()`.

* **Polynomial Commitment Schemes (PCS)** (`plonkish_backend/src/pcs/`):
    * Defines a generic `PolynomialCommitmentScheme` trait, extended with interfaces for shift operations like `open_shift` and `verify_shifted_evaluation`.
    * **Zeromorph** (`pcs/multilinear/zeromorph.rs`) is the flagship feature of this layer. It serves as the primary PCS for multilinear polynomials, providing efficient commitment and proof mechanisms with native support for optimized shift proofs.
    * It is built upon a univariate KZG commitment scheme (`pcs/univariate/kzg.rs`).

* **Frontend Layer** (`plonkish_backend/src/frontend/`):
    * This layer is responsible for translating specific computational problems (circuits) into a format that the Plonkish backend can process.
    * Currently, it supports integration with Zcash's **Halo2** circuit frontend via the `frontend-halo2` feature flag, allowing developers to define circuits using the familiar Halo2 API.

* **Polynomial Representations** (`plonkish_backend/src/poly/`):
    * Defines the core data structures used throughout the system, primarily `MultilinearPolynomial` and `UnivariatePolynomial`.
    * The multilinear polynomial module includes support for evaluating rotated polynomials, which is key to handling wire relations in HyperPlonk.

* **Utilities** (`plonkish_backend/src/util/`):
    * Provides a collection of foundational tools, including:
        * **Arithmetic**: FFTs, multi-scalar multiplications (MSMs), etc.
        * **Expression Engine**: A system for defining and evaluating polynomial constraints.
        * **Transcript Management**: An implementation of the Fiat-Shamir transform to convert interactive protocols into non-interactive ones.
        * **Parallelization**: Leverages the `rayon` library to accelerate compute-intensive tasks like proof generation.

## 3. Core Concepts & Technology

### 3.1 HyperPlonk and the Sum-Check Protocol

HyperPlonk is an evolution of the Plonk protocol that shifts the underlying mathematical foundation from univariate polynomials to multilinear polynomials over the boolean hypercube. The primary benefit of this transition is the replacement of the costly FFT with the **Sum-check protocol**. The Sum-check protocol is a highly efficient interactive protocol for verifying the sum of a multilinear polynomial's evaluations over the entire boolean hypercube. Its prover has a linear time complexity, offering a significant performance improvement over the quasi-linear complexity of the FFT.

#### The Limitations of Plonk

In the standard Plonk protocol, the execution trace (witness) and constraints (such as gate constraints and copy constraints) are encoded as univariate polynomials. To prove that these constraints hold over the entire evaluation domain, the protocol must perform numerous polynomial operations, the most expensive of which are the FFT and its inverse (IFFT). As the circuit size $N$ increases, the $O(N \log N)$ complexity of the FFT becomes a significant performance bottleneck for proof generation.

#### HyperPlonk's Paradigm Shift: Embracing Multilinear Polynomials

HyperPlonk takes a different approach. Instead of viewing the $N$ values of a circuit as evaluations over a multiplicative subgroup, it treats them as evaluations over the **boolean hypercube** $\{0, 1\}^k$ (where $N = 2^k$). Consequently, the circuit's witnesses and selectors are represented as multilinear polynomials in $k$ variables.

A multilinear polynomial $p(X_1, ..., X_k)$ is special because the degree of each variable $X_i$ is at most 1. This structure aligns perfectly with the binary nature of the boolean hypercube.

All circuit constraints—whether custom gate constraints or copy constraints between wires—can ultimately be combined into a single, large multilinear polynomial identity. This identity must hold (typically, evaluate to zero) for all $2^k$ points on the boolean hypercube. For example, a total constraint polynomial $C(X_1, ..., X_k)$ must satisfy:
$C(x) = 0$ for all $x \in \{0, 1\}^k$.

The core challenge for the prover then becomes: **How can one efficiently prove to a verifier that a given multilinear polynomial $C$ is identically zero over the entire boolean hypercube?**

This is precisely where the Sum-check protocol comes into play.

#### The Sum-Check Protocol: An Efficient Tool for Summation Verification

The Sum-check protocol is a concise yet powerful interactive protocol. It allows a Prover to convince a Verifier that the sum of a multilinear polynomial $g(X_1, ..., X_k)$ over all points in the boolean hypercube equals a publicly claimed value $S$.

> **Goal**: Prove that $\sum_{x \in \{0, 1\}^k} g(x) = S$

The protocol proceeds as follows:

1.  **Initialization**: The Prover claims to the Verifier that $\sum_{x_1,...,x_k \in \{0,1\}} g(x_1, ..., x_k) = S$.

2.  **Interactive Rounds (k rounds in total)**:
    * **Round 1**:
        * **Prover**: Fixes the first variable $X_1$ and sums $g$ over all possible values of the other variables $X_2, ..., X_k \in \{0, 1\}$. This results in a univariate polynomial in $X_1$, denoted $p_1(X_1) = \sum_{x_2,...,x_k \in \{0,1\}} g(X_1, x_2, ..., x_k)$. The degree of $p_1$ is low (at most the degree of $g$). The Prover sends $p_1(X_1)$ to the Verifier.
        * **Verifier**: Upon receiving $p_1$, it performs a consistency check: it verifies if $p_1(0) + p_1(1)$ equals the claimed sum $S$ from the previous step. If the check fails, the proof is rejected. If it succeeds, the Verifier generates a random challenge $r_1$ and sends it to the Prover.

    * **Round 2**:
        * **Prover**: "Binds" the received challenge $r_1$ to the first variable of $g$, creating a new polynomial with fewer variables: $g_1(X_2, ..., X_k) = g(r_1, X_2, ..., X_k)$. The goal for this round is to prove that $\sum_{x_2,...,x_k \in \{0,1\}} g_1(x_2, ..., x_k)$ equals $p_1(r_1)$. The Prover repeats the process, calculating and sending $p_2(X_2)$.
        * **Verifier**: Performs the same consistency check, $p_2(0) + p_2(1) == p_1(r_1)$, and then generates a new random challenge $r_2$.

    * **... and so on ...**

3.  **Final Check**: After $k$ rounds of interaction, all variables $X_1, ..., X_k$ have been replaced by random challenges $r_1, ..., r_k$. The Verifier only needs to perform one final check: it directly computes the value of $g(r_1, ..., r_k)$ and verifies that it equals the evaluation of the polynomial sent by the Prover in the final round at point $r_k$, i.e., $g(r_1, ..., r_k) == p_k(r_k)$.

**The magic of the Sum-check protocol lies in:**

* **Dimensionality Reduction**: It transforms a complex problem involving $2^k$ summations into a single evaluation of the polynomial $g$ at one random point $(r_1, ..., r_k)$ after $k$ rounds of interaction.
* **Efficiency**: In each round, the Prover's main task is to compute $p_i$. Therefore, the Prover's total complexity is $O(N)$, where $N = 2^k$. 

#### How HyperPlonk Utilizes the Sum-Check Protocol

Returning to HyperPlonk, the objective is to prove that the constraint polynomial $C(x)$ is zero for all $x \in \{0, 1\}^k$. This is equivalent to proving that $\forall y,\sum_{x \in \{0, 1\}^k} C(x)  eq(x,y) = 0$, where $eq(x,y)=\prod_i^k(x_iy_i + (1-x_i)(1-y_i))$ . By the Schwartz-Zippel lemma, this equality holds with high probability only if $C(x)$ is identically zero.

Therefore, HyperPlonk bundles all circuit constraints into $C(x)$ and then leverages the Sum-check protocol to prove to the verifier that the sum of $C(x) eq(x,y)$ is indeed zero polynomial. This completely replaces the quotient polynomial and associated FFTs used in Plonk for handling permutation arguments and gate constraints.

In summary, by shifting its computational foundation from univariate polynomials and FFTs to multilinear polynomials and the Sum-check protocol, HyperPlonk successfully reduces the prover's computational complexity from $O(N \log N)$ to $O(N)$. This paves the way for building larger and more complex zero-knowledge applications.



### 3.2 Polynomial Commitment Schemes (PCS) and Zeromorph

A Polynomial Commitment Scheme (PCS) is a core component of all modern ZK-SNARKs. It allows a prover to "commit" to a polynomial, producing a short digest. Later, the prover can prove the polynomial's evaluation at any point without revealing the entire polynomial.

**Zeromorph** introduces an elegant and generic construction that **reduces the problem of committing to and evaluating multilinear polynomials to the well-understood domain of univariate polynomials**. This approach is not only conceptually clean but also offers exponential performance improvements for achieving the zero-knowledge property.

#### The Core Idea: The Multilinear-to-Univariate Isomorphism

The foundation of Zeromorph is a linear isomorphism, denoted $\mathcal{U}_n$, between the vector space of $n$-variate multilinear polynomials and the space of univariate polynomials with degree less than $2^n$.

> **The $\mathcal{U}_n$ Map**: $\mathbb{F}[X_0, ..., X_{n-1}]^{\le1} \rightarrow \mathbb{F}[X]^{<2^n}$

Intuitively, this map takes the $2^n$ evaluations of a multilinear polynomial over the boolean hypercube $\{0, 1\}^n$ and treats them as the $2^n$ coefficients of its corresponding univariate polynomial.

For example, for a 2-variable multilinear polynomial $f(X_0, X_1)$ with evaluations $a_{00}, a_{10}, a_{01}, a_{11}$ at points $(0,0), (1,0), (0,1), (1,1)$ respectively, the $U_2$ map produces the following univariate polynomial:
$\hat{f}(X) = a_{00} + a_{10}X + a_{01}X^2 + a_{11}X^3$

This $\hat{f}(X)$ is referred to as the **"univariatisation"** of $f$.

#### The Zeromorph Commitment Mechanism

Leveraging this map, the commitment process in Zeromorph becomes straightforward:

1.  **Commit**: To commit to an $n$-variate multilinear polynomial $f$, the prover first computes its corresponding univariate polynomial $\hat{f} = \mathcal{U}_n(f)$.
2.  Then, it uses any **additively homomorphic univariate PCS** (like the hiding version of KZG used in Plonkish) to commit to $\hat{f}$, generating the commitment $C$.

This design brilliantly converts a complex multivariate problem into a univariate one, for which mature and highly efficient tools already exist.

#### The Zeromorph Evaluation Proof Protocol

The true power of Zeromorph is revealed in its evaluation proof protocol. To prove that a multilinear polynomial $f$ evaluates to $v$ at a point $u = (u_0, ..., u_{n-1})$, i.e., $f(u) = v$, the protocol ingeniously combines two core theoretical concepts:

1.  **Euclidean Division for Multilinear Polynomials**:
    An $n$-variate multilinear polynomial $f$ satisfies $f(u) = v$ if and only if there exist $n$ quotient polynomials $q_0, ..., q_{n-1}$ (where $q_k$ is a $k$-variate multilinear polynomial) such that the following identity holds:
    $f - v = \sum_k (X_k - u_k) * q_k$

2.  **Linearity of the Isomorphism**:
    Because $\mathcal{U}_n$ is a linear map, the multilinear polynomial identity above holds if and only if its "univariatised" version also holds:
    $\mathcal{U}_n(f) - \mathcal{U}_n(v) = \sum_k [\mathcal{U}_n(X_k * q_k) - u_k * \mathcal{U}_n(q_k)]$

The core of the Zeromorph protocol is to **prove that this second, univariate polynomial identity holds**.

**Protocol Overview:**

1.  **Prover's Preparation**:
    * The Prover computes all quotient polynomials $q_k$.
    * For each $q_k$, it computes its univariatisation, $\mathcal{U}_k(q_k)$.
    * It then uses the underlying univariate PCS (e.g., KZG) to commit to each $\mathcal{U}_k(q_k)$, sending these commitments $[\mathcal{U}_k(q_k)]$ to the Verifier.

2.  **Verifier's Challenge**:
    * Upon receiving the commitments to the quotient polynomials, the Verifier must first run a **degree-check protocol** to ensure that each committed $\mathcal{U}_k(q_k)$ indeed has a degree no greater than $2^k - 1$. This is a crucial step for the protocol's soundness.
    * After the degree checks pass, the Verifier sends a random challenge point $z$ to the Prover.

3.  **Final Single-Point Evaluation**:
    * Using the homomorphic properties of the commitments $[\mathcal{U}_n(f)]$ and $[\mathcal{U}_k(q_k)]$, the Prover and Verifier can jointly compute a commitment to the "target polynomial" derived from the main univariate identity evaluated at $z$.
    * The Prover's final task is to provide a **standard, univariate evaluation proof** showing that this combined target polynomial evaluates to 0 at the point $z$.

Through this sequence of steps, a complex multilinear evaluation problem is elegantly reduced to committing to several smaller univariate polynomials, performing degree checks, and finally proving a single evaluation on a combined univariate polynomial.

This construction drastically lowers the overhead for achieving zero-knowledge. Traditional methods require a large, random "masking polynomial" to hide the witness, which is costly. In contrast, the zero-knowledge version of Zeromorph incurs extra costs primarily from committing to the $n$ small quotient polynomials, reducing the overhead from an exponential $O(2^n)$ to a near-constant $O(n)$. This represents a monumental leap in performance.



### 3.3 The Shift Operation Optimization

In many advanced zero-knowledge applications, such as state machines (zk-EVMs) and lookup arguments, a critical requirement is to prove the correctness of transitions between different states. In the context of multilinear polynomials, this often boils down to proving an evaluation on a "shifted" version of a polynomial.

Specifically, if we represent a state vector of size $N = 2^n$ with an $n$-variate multilinear polynomial $f$, its "shift," denoted as $f_{\leftarrow}$, represents the evolution of the state. The evaluation vector of $f_{\leftarrow}$ is $(a_1, ..., a_{N-1}, a_0)$. The prover needs to convince the verifier that for a given point $u$, the shifted polynomial $f_{\leftarrow}$ evaluates to $v$, i.e., $f_{\leftarrow}(u) = v$.

#### The Inefficiency of Traditional Methods

Directly proving an evaluation on a shifted polynomial like $f_{\leftarrow}$ is highly inefficient. A naive approach would be for the prover to commit to a new polynomial $f_{\leftarrow}$, but this introduces significant computational and communication overhead. And it also needs to prove the relation between $f_{\leftarrow}$ and $f$. The original methods proposed in earlier works, like HyperPlonk, while functional, incur an exponential cost in a long-distance shift setting, which is prohibitive.

#### Zeromorph's Breakthrough: A Univariate Identity

Zeromorph ingeniously bypasses the need to introduce a new commitment for the shifted polynomial. The core insight comes from a simple but powerful **univariate polynomial identity** that connects the original polynomial $f$ with its shifted version $f_{\leftarrow}$.

Recall that Zeromorph's foundation is the mapping $\mathcal{U}_n$ which converts an $n$-variate multilinear polynomial $f$ into a univariate polynomial $\mathcal{U}_n(f)$ of degree less than $2^n$. The Zeromorph paper establishes the following crucial relationship:

$X \cdot \mathcal{U}_n(f_{\leftarrow}) = \mathcal{U}_n(f) - a_0 + a_0 X^N$

Here, $a_0 = f(0, ..., 0)$ is the evaluation of the original polynomial at the zero point. This identity is the cornerstone of Zeromorph's efficient shift proof.

#### Transforming a "Shift Proof" into a "Standard Proof"

Using this identity, a proof of evaluation on a **shifted polynomial $f_{\leftarrow}$** can be transformed into a proof involving only the **original polynomial commitment $f$** and simple algebraic manipulations.

The protocol proceeds as follows:

1. **Commit $a_0$**: The prover first needs to commit $a_0$. 

2. **Reframe the Proof Goal**: With $a_0$ committed, the task of proving $f_{\leftarrow}(u) = v$ is equivalent to proving the identity $f_{\leftarrow} - v = \sum_k (X_k - u_k)q_k$. By applying the $\mathcal{U}_n$ map and combining it with the core univariate identity, the final goal becomes verifying the following univariate polynomial identity:

   $\mathcal{U}_n(f) - a_0 + a_0X^N - X \cdot \mathcal{U}_n(v) = X \cdot \sum_k ( ... ) \mathcal{U}_n(q_k)^{<2^k}$

3. **Run the Standard Protocol**: Now, the prover and verifier can execute a protocol almost identical to the standard Zeromorph evaluation protocol, centered around this new univariate identity. The prover commits to the quotient polynomials $q_k$, the verifier issues a random challenge $z$, and the prover provides a final univariate evaluation proof.

Compared to the exponential cost of previous methods, Zeromorph reduces the cost of a shift proof to be nearly identical to that of a single standard evaluation proof. This optimization is crucial for the performance of Plonkish in applications that involve numerous state transitions, such as zk-EVMs.



### 3.4 The Anemoi Hash Function: A Custom Gate for SNARK-Optimized Hashing

In complex applications requiring "verifiable state," such as a zk-EVM, an efficient hash function is the cornerstone for building Merkle trees to authenticate memory access and state transitions. However, standard cryptographic hash functions like SHA-256 or Keccak are prohibitively expensive to implement in arithmetic circuits, often requiring tens of thousands of gates, which becomes a major bottleneck for proof generation. To overcome this, the field of cryptography has developed "SNARK-friendly" hash functions, whose algebraic structures are meticulously designed for optimal representation within zero-knowledge circuits. **Anemoi** is a state-of-the-art example in this domain.

A core advantage of Plonkish is its ability to support high-degree "custom gates." This allows us to "hard-code" the complex logic of the Anemoi hash algorithm directly into a highly-optimized gate constraint, thereby enabling the verification of hash computations with minimal circuit overhead.

#### The Anemoi Permutation

Anemoi follows the classic **Substitution-Permutation Network (SPN)** design. For a state consisting of $2\ell$ field elements, each round $r$ of the computation involves the following steps:

- **Constant addition**: For input $(\vec{x}, \vec{y}) \in\left(\mathbb{F}^{\ell}, \mathbb{F}^{\ell}\right)$, the $r$-th round has some round-specific constants $\vec{c}_r \in \mathbb{F}^{\ell}$ and $\vec{d}_r \in \mathbb{F}^{\ell}$. It outputs $\left(\vec{x}+\vec{c}_r, \vec{y}+\vec{d}_r\right)$.
- **MDS matrix**: The MDS matrix in the Anemoi permutation is a fixed matrix $M$ of size $\mathbb{F}^{\ell \times \ell}$. For input $(\vec{x}, \vec{y}) \in\left(\mathbb{F}^{\ell}, \mathbb{F}^{\ell}\right)$, we have $\vec{u}=M \cdot \vec{x}$ and $\vec{v}=M \cdot \vec{y}_\omega$ where $\vec{y}_\omega$ a shifted version of $\vec{y}$ (in the Amenoi permutation, moving the first element to the end). The Anemoi permutation uses a shift instead of having another MDS matrix for $\vec{y}$.
- **Pseudo-Hadamard transform**: For input $(\vec{x}, \vec{y}) \in\left(\mathbb{F}^{\ell}, \mathbb{F}^{\ell}\right)$, an additional step is made after the MDS matrix to mix $\vec{x}$ and $\vec{y}$. The mix needs to be invertible. In Anemoi, this is done by having $\vec{v}:=\vec{y}+\vec{x}$ and $\vec{u}:=\vec{y}+2 \vec{x}$. It outputs $(\vec{u}, \vec{v})$.
- **S-box**: For input $(\vec{x}, \vec{y}) \in\left(\mathbb{F}^{\ell}, \mathbb{F}^{\ell}\right)$, and we let the S -box be $S(x, y) \rightarrow\left(x^{\prime}, y^{\prime}\right)$, then the output $\vec{u}$ and $\vec{v}$ can be computed by letting $(u[i], v[i])=S(x[i], y[i])$ for $i=1, \ldots, \ell$.

The key innovation of Anemoi is its **Flystel S-box**, which is based on an "extended Feistel" construction. While the forward computation of the S-box (used offline to generate the hash) can involve complex exponentiations (e.g., $x^{1/\alpha}$), its inverse relation (used for in-circuit verification) can be expressed as a **low-degree polynomial**. This allows the prover to efficiently compute the hash offline while only needing to satisfy a low-degree constraint within the circuit.

#### Anemoi Constraints for Plonkish

To integrate Anemoi into a system like Plonkish, we must translate its round permutation into a set of polynomial constraints. As we will use it in the Merkle tree, we set the state contain 4 elements. One round of the Anemoi permutation can be encoded using **4 core constraints** over **5 witness columns** $(w_1, w_2, w_3, w_4, w_o)$ to store the intermediate states.

Assume that row $i$ of the circuit represents one round of the Anemoi permutation. The witness columns at this row store the input state and some intermediate values, while the witness columns at the next row (row $i+1$) store the output state of the current round.

**1. State Mapping and Processed Round Keys (PRK):**

First, the 4-element state is mapped to 4 of the witness columns. For a given row $X$ in the circuit (representing one round of computation), the witness columns store the input state, while the next row (accessed via a shift, $X\omega$) stores the output state.

* Input state $(x[1], x[2], y[1], y[2])$ corresponds to $(w_1(X), w_2(X), w_3(X), w_4(X))$.
* Output state $(x'[1], x'[2], y'[1], y'[2])$ corresponds to $(w_1(X\omega), w_2(X\omega), w_3(X\omega), w_4(X\omega))$.

The linear layer of the Anemoi permutation (which includes constant addition, MDS matrix multiplication, and a pseudo-Hadamard transform) can be "pre-processed" and their effects inlined directly into the constraints. The results of this pre-computation are called "Processed Round Keys" (PRKs). In the Plonkish circuit, these are encoded as **preprocessed polynomials (selectors)**, denoted $q_{prk1}(X)$, $q_{prk2}(X)$, $q_{prk3}(X)$, and $q_{prk4}(X)$.

**2. The Core Constraint Equations:**

The correctness of the Anemoi permutation is enforced by the following four constraint equations. These equations leverage the low-degree verification property of the Flystel S-box to relate the input state, output state, and the PRKs. The PRK polynomials also act as **selectors**: they are non-zero only on the rows where an Anemoi computation is being performed, thereby "activating" the constraints.

* **Constraint 1**:
    $q_{prk3}(X) \cdot [ ( (w_1(X) + w_4(X) + g \cdot (w_2(X) + w_3(X)) + q_{prk3}(X)) - w_3(X\omega) )^5 + g \cdot (w_1(X) + w_4(X) + g \cdot (w_2(X) + w_3(X)) + q_{prk3}(X))^2 - (2w_1(X) + w_4(X) + g \cdot (2w_2(X) + w_3(X)) + q_{prk1}(X)) ] = 0$

* **Constraint 2**:
    $q_{prk3}(X) \cdot [ ( (g \cdot (w_1(X) + w_4(X)) + (g^2 + 1) \cdot (w_2(X) + w_3(X)) + q_{prk4}(X)) - w_4(X\omega) )^5 + g \cdot (g \cdot (w_1(X) + w_4(X)) + (g^2 + 1) \cdot (w_2(X) + w_3(X)) + q_{prk4}(X))^2 - (g \cdot (2w_1(X) + w_4(X)) + (g^2 + 1) \cdot (2w_2(X) + w_3(X)) + q_{prk2}(X)) ] = 0$

* **Constraint 3**:
    $q_{prk3}(X) \cdot [ ( (w_1(X) + w_4(X) + g \cdot (w_2(X) + w_3(X)) + q_{prk3}(X)) - w_3(X\omega) )^5 + g \cdot (w_3(X\omega))^2 + g^{-1} - w_1(X\omega) ] = 0$

* **Constraint 4**:
    $q_{prk3}(X) \cdot [ ( (g \cdot (w_1(X) + w_4(X)) + (g^2 + 1) \cdot (w_2(X) + w_3(X)) + q_{prk4}(X)) - w_4(X\omega) )^5 + g \cdot (w_4(X\omega))^2 + g^{-1} - w_2(X\omega) ] = 0$

Here, $g$ is a generator of the field. With these four constraints, we can verify one round of the Anemoi permutation with extreme efficiency inside the Plonkish circuit.



#### Jive Mode and Merkle Trees

To build an efficient **Collision-Resistant Hash (CRH)**, the Anemoi paper introduces the **Jive** mode of operation. Jive utilizes the Anemoi permutation function, denoted $P$, to construct a hash function, with the core idea being a linear summation after the permutation is applied.

<img src="https://anemoi-hash.github.io/Figures/mode.png" alt="Jive Mode" style="zoom:50%;" />

For instance, for a 2-to-1 hash that maps two elements $(x, y)$ to a single hash output $h$, the Jive mode computation is:
1.  Compute the permutation: $(u, v) = P(x, y)$
2.  Compute the final hash: $h = u + v + x + y$

This construction allows Plonkish to verify a **3-to-1 Merkle tree compression** (where a parent node is the hash of three children) with remarkable efficiency. Plonkish implement Anemoi-Jive as a custom gate, verifying one 3-to-1 compression requires **only about 20 gates**, which averages to about 1 gate per round. This represents a monumental performance leap compared to traditional hash functions, which can require hundreds or thousands of gates, and it dramatically reduces the cost of verifying state storage and transitions in applications like zk-EVMs.

## 4. Implementation Details & Usage

### 4.1 Language and Ecosystem

The project is written entirely in **Rust** and is deeply integrated with the **`arkworks`** ecosystem, a comprehensive suite of Rust libraries for zero-knowledge cryptography. `arkworks` provides the foundational finite field arithmetic, elliptic curve operations, and polynomial manipulation routines that underpin Plonkish's efficient implementation.

### 4.2 Main Interfaces and Data Structures

* **`PlonkishBackend` Trait**: Defines the core lifecycle of the proof system, including methods like `setup`, `preprocess`, `prove`, `verify`, and the specialized `prove_with_shift` and `verify_with_shift`.
* **`PlonkishCircuit` Trait**: The frontend abstraction for a circuit, responsible for providing circuit information (`circuit_info`) and synthesizing the witness (`synthesize`).
* **`PlonkishCircuitInfo` Struct**: A descriptor containing all metadata for a circuit, including:
    * `k`: The circuit size (as a power of two, 2^k).
    * `num_instances`: The number of public input polynomials.
    * `preprocess_polys`: Preprocessed polynomials (e.g., selector polynomials).
    * `num_witness_polys`: The number of witness polynomials.
    * `constraints`: The circuit constraints, represented as `Expression<F>`.
    * `lookups`: Parameters for lookup arguments.
    * `permutations`: Parameters for permutation arguments (copy constraints).

## 5. Conclusion

By integrating the cutting-edge technologies of HyperPlonk, Zeromorph, and Anemoi, Plonkish stands as a high-performance and extensible zero-knowledge proving system. Its novel optimization for shift operations and its support for efficient hash functions give it a distinct advantage in applications that handle complex state transitions and large datasets.



## Reference

[1] Chen B, Bünz B, Boneh D, et al. Hyperplonk: Plonk with linear-time prover and high-degree custom gates[C]//Annual International Conference on the Theory and Applications of Cryptographic Techniques. Cham: Springer Nature Switzerland, 2023: 499-530.

[2] Kohrita T, Towa P. Zeromorph: Zero-Knowledge Multilinear-Evaluation Proofs from Homomorphic Univariate Commitments: T. Kohrita, P. Towa[J]. Journal of Cryptology, 2024, 37(4): 38.

[3] Liu J, Patil H, Peddireddy A S, et al. An efficient verifiable state for zk-EVM and beyond from the Anemoi hash function[J]. Cryptology ePrint Archive, 2022.
