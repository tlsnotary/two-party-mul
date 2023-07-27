//! Implements the two-party multiplication protocol (Protocol 1, page 12) from
//! <https://eprint.iacr.org/2019/523>

use mpz_share_conversion_core::Field;
use rand::{rngs::ThreadRng, thread_rng, Rng};
use std::marker::PhantomData;

mod dot;
use dot::DotProduct;

// Some constants defined in the paper
//
// Bits needed for to represent elements of the field
const KAPPA: usize = 256;

// Security parameter
const S: usize = 80;

// Redundancy factor
const ZETA: usize = KAPPA + 2 * S;

// Batch size
const L: usize = 1;

// Number of OTs
const ETA: usize = ZETA * L;

pub struct M2A<T: Field> {
    _field: PhantomData<T>,
    rng: ThreadRng,
}

impl<T: Field> M2A<T> {
    fn gadget(&mut self) -> [T; ZETA] {
        std::array::from_fn(|_| T::rand(&mut self.rng))
    }

    fn alpha(&mut self, a_tilde: [T; L], a_head: [T; L]) -> [[T; ETA]; 2] {
        let mut alpha = [[T::zero(); ETA]; 2];

        for k in 0..L {
            alpha[0][k * ZETA..(k + 1) * ZETA].copy_from_slice(&[a_tilde[k]; ZETA]);
            alpha[1][k * ZETA..(k + 1) * ZETA].copy_from_slice(&[a_head[k]; ZETA]);
        }

        alpha
    }

    fn a_tilde(&mut self) -> [T; L] {
        std::array::from_fn(|_| T::rand(&mut self.rng))
    }

    fn a_head(&mut self) -> [T; L] {
        std::array::from_fn(|_| T::rand(&mut self.rng))
    }

    fn beta(&mut self) -> [T; ETA] {
        let mut beta = [T::zero(); ETA];
        beta.iter_mut().for_each(|el| {
            if self.rng.gen() {
                *el = T::one()
            }
        });
        beta
    }

    fn b_tilde(&self, gadget: [T; ZETA], beta: [T; ETA]) -> [T; L] {
        let mut b_tilde = [T::zero(); L];
        for k in 0..L {
            let mut beta_part = [T::zero(); ZETA];
            beta_part.copy_from_slice(&beta[k * ZETA..(k + 1) * ZETA]);
            b_tilde[k] = gadget.dot(beta_part);
        }
        b_tilde
    }

    fn omega_a(&mut self) -> [[T; ETA]; 2] {
        let first = std::array::from_fn(|_| T::rand(&mut self.rng));
        let second = std::array::from_fn(|_| T::rand(&mut self.rng));

        [first, second]
    }

    fn chi_head(&self) -> [T; L] {
        std::array::from_fn(|_| T::rand(&mut thread_rng()))
    }

    fn chi_tilde(&self) -> [T; L] {
        std::array::from_fn(|_| T::rand(&mut thread_rng()))
    }

    fn r(
        &self,
        chi_tilde: [T; L],
        chi_head: [T; L],
        z_tilde_a: [T; ETA],
        z_head_a: [T; ETA],
    ) -> [T; ZETA] {
        let mut r = [T::zero(); ZETA];

        for j in 0..ZETA {
            r[j] = (0..L).fold(T::zero(), |acc, i| {
                acc + chi_tilde[i] * z_tilde_a[i * ZETA + j] + chi_head[i] * z_head_a[i * ZETA + j]
            });
        }

        r
    }

    fn u(&self, chi_tilde: [T; L], chi_head: [T; L], a_tilde: [T; L], a_head: [T; L]) -> [T; L] {
        let mut u = [T::zero(); L];

        u.iter_mut()
            .enumerate()
            .for_each(|(k, el)| *el = chi_tilde[k] * a_tilde[k] + chi_head[k] * a_head[k]);

        u
    }

    fn check(
        &self,
        r: [T; ZETA],
        u: [T; L],
        beta: [T; ETA],
        chi_tilde: [T; L],
        chi_head: [T; L],
        z_tilde_b: [T; ETA],
        z_head_b: [T; ETA],
    ) -> bool {
        for j in 0..ZETA {
            let lhs = r[j]
                + (0..L).fold(T::zero(), |acc, i| {
                    acc + chi_tilde[i] * z_tilde_b[i * ZETA + j]
                        + chi_head[i] * z_head_b[i * ZETA + j]
                });

            let rhs = (0..L).fold(T::zero(), |acc, i| acc + beta[i * ZETA + j] * u[i]);

            if lhs != rhs {
                return false;
            }
        }

        true
    }

    fn gamma(&self, input: [T; L], tilde: [T; L]) -> [T; L] {
        let mut gamma = [T::zero(); L];

        gamma
            .iter_mut()
            .enumerate()
            .for_each(|(k, el)| *el = input[k] + -tilde[k]);

        gamma
    }

    fn output(&self, input: [T; L], gamma: [T; L], gadget: [T; ZETA], z_tilde: [T; ETA]) -> [T; L] {
        let mut z = [T::zero(); L];

        for (i, el) in z.iter_mut().enumerate() {
            *el = input[i] * gamma[i]
                + (0..ZETA).fold(T::zero(), |acc, j| acc + gadget[j] * z_tilde[i * ZETA + j]);
        }
        z
    }
}

impl<T: Field> Default for M2A<T> {
    fn default() -> Self {
        let rng = thread_rng();

        Self {
            _field: PhantomData,
            rng,
        }
    }
}

fn func_cote<T: Field>(
    alpha: [[T; ETA]; 2],
    beta: [T; ETA],
    omega_a: [[T; ETA]; 2],
) -> ([[T; ETA]; 2], [[T; ETA]; 2]) {
    let mut omega_b: [[T; ETA]; 2] = [[T::zero(); ETA]; 2];

    for k in 0..ETA {
        omega_b[0][k] = alpha[0][k] * beta[k] + -omega_a[0][k];
        omega_b[1][k] = alpha[1][k] * beta[k] + -omega_a[1][k];
    }
    (omega_a, omega_b)
}

pub fn protocol<T: Field>(
    a: [T; L],
    b: [T; L],
    malicious_omega_a: Option<[[T; ETA]; 2]>,
) -> ([T; L], [T; L]) {
    //Setup
    let mut m2a = M2A::<T>::default();
    let gadget = m2a.gadget();

    // Step 1
    let beta = m2a.beta();
    let b_tilde = m2a.b_tilde(gadget, beta);

    // Step 2
    let a_tilde = m2a.a_tilde();
    let a_head = m2a.a_head();
    let alpha = m2a.alpha(a_tilde, a_head);

    // Step 3
    // If Alice is malicious, she can choose omega_a in the implementation of F_COTE
    let (omega_a, omega_b) = if let Some(omega_a) = malicious_omega_a {
        func_cote(alpha, beta, omega_a)
    } else {
        // If she is not malicious she supplies random input
        func_cote(alpha, beta, m2a.omega_a())
    };

    let [z_tilde_a, z_head_a] = [omega_a[0], omega_a[1]];
    let [z_tilde_b, z_head_b] = [omega_b[0], omega_b[1]];

    // Step 4
    let (chi_tilde, chi_head) = (m2a.chi_tilde(), m2a.chi_head());

    // Step 5
    let r = m2a.r(chi_tilde, chi_head, z_tilde_a, z_head_a);
    let u = m2a.u(chi_tilde, chi_head, a_tilde, a_head);

    // Step 6
    let check = m2a.check(r, u, beta, chi_tilde, chi_head, z_tilde_b, z_head_b);
    if !check {
        panic!("Consistency check failed!");
    }

    // Step 7
    let gamma_a = m2a.gamma(a, a_tilde);
    let gamma_b = m2a.gamma(b, b_tilde);

    // Step 8
    let z_a = m2a.output(a, gamma_b, gadget, z_tilde_a);
    let z_b = m2a.output(b_tilde, gamma_a, gadget, z_tilde_b);

    (z_a, z_b)
}

#[cfg(test)]
mod tests {
    use mpz_share_conversion_core::fields::{p256::P256, UniformRand};

    use super::*;

    #[test]
    fn test_two_party_mul() {
        // Get rng
        let mut rng = thread_rng();

        // Random input
        let a = [P256::rand(&mut rng)];
        let b = [P256::rand(&mut rng)];

        // Execute protocol
        // Alice is honest
        let (z_a, z_b) = protocol(a, b, None);

        // Check result
        for k in 0..L {
            assert_eq!(z_a[k] + z_b[k], a[k] * b[k]);
        }
    }

    #[test]
    fn test_alice_malicious() {
        // Get rng
        let mut rng = thread_rng();

        // Alice is malicious and supplies zero as input. This allows her to infer Bob's output z_b
        let a = [P256::zero()];
        let b = [P256::rand(&mut rng)];

        // Execute protocol
        // No need for Alice to choose omega_a
        let (z_a, z_b) = protocol(a, b, None);

        // Check result
        for k in 0..L {
            //Protocol works
            assert_eq!(z_a[k] + z_b[k], a[k] * b[k]);

            //Check Bob's output. Indeed z_b = -z_a
            assert_ne!(z_b[k], -z_a[k]);
        }
    }
}
