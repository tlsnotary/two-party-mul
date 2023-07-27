//! Implements the two-party multiplication protocol (Protocol 1, page 12) from
//! <https://eprint.iacr.org/2019/523>

use mpz_share_conversion_core::Field;
use rand::{rngs::ThreadRng, thread_rng, Rng};
use std::marker::PhantomData;

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

    // Note that the return dimension is [[[T; 2]; ZETA]; L] which we simplify a little bit
    fn alpha(&mut self, a_tilde: [T; L], a_head: [T; L]) -> [[T; 2]; ETA] {
        let mut alpha = [[T::zero(); 2]; ETA];

        for k in 0..L {
            for i in 0..ZETA {
                alpha[k * ZETA + i] = [a_tilde[k], a_head[k]];
            }
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

    fn omega_a(&mut self) -> [[T; 2]; ETA] {
        std::array::from_fn(|_| [T::rand(&mut self.rng), T::rand(&mut self.rng)])
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

    fn bob_check(
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

    //TODO: Is there a mistake in the paper with b_tilde ?
    fn output(&self, input: [T; L], gamma: [T; L], gadget: [T; ZETA], z_tilde: [T; ETA]) -> [T; L] {
        let mut z = [T::zero(); L];

        for (i, el) in z.iter_mut().enumerate() {
            *el = input[i] * gamma[i]
                + (0..ZETA).fold(T::zero(), |acc, j| acc + gadget[i] * z_tilde[i * ZETA + j]);
        }
        z
    }
}

impl<T: Field> Default for M2A<T> {
    fn default() -> Self {
        let mut rng = thread_rng();

        Self {
            _field: PhantomData,
            rng,
        }
    }
}

pub fn func_cote<T: Field>(
    alpha: [[T; 2]; ETA],
    beta: [T; ETA],
    omega_a: [[T; 2]; ETA],
) -> ([[T; 2]; ETA], [[T; 2]; ETA]) {
    let mut omega_b: [[T; 2]; ETA] = [[T::zero(); 2]; ETA];

    for k in 0..ETA {
        omega_b[k] = [
            alpha[k][0] * beta[k] + -omega_a[k][0],
            alpha[k][1] * beta[k] + -omega_a[k][1],
        ];
    }
    (omega_a, omega_b)
}

trait DotProduct: Copy {
    type Output;

    fn dot(self, other: Self) -> Self::Output;
}

impl<T: Field, const N: usize> DotProduct for [T; N] {
    type Output = T;

    fn dot(self, other: Self) -> Self::Output {
        self.into_iter()
            .zip(other)
            .fold(T::zero(), |acc, (a, b)| acc + a * b)
    }
}

#[cfg(test)]
mod tests {
    use mpz_share_conversion_core::fields::p256::P256;

    use super::*;

    #[test]
    fn test_two_party_mul() {
        let mut m2a = M2A::<P256>::default();
        let gadget = m2a.gadget();

        // Step 1
        let beta = m2a.beta();
        let b_tilde = m2a.b_tilde(gadget, beta);

        // Step 2
        let a_tilde = m2a.a_tilde();
        let a_head = m2a.a_head();
        let alpha = m2a.alpha(a_tilde, a_head);

        // Step 3
        let (omega_a, omega_b) = func_cote(alpha, beta, m2a.omega_a());
        let [z_tilde_a, z_head_a] = [omega_a[0], omega_a[1]];
        let [z_tilde_b, z_head_b] = [omega_b[0], omega_b[1]];

        // Step 4
        let (chi_tilde, chi_head) = (m2a.chi_tilde(), m2a.chi_head());

        // Step 5
        let r = m2a.r(chi_tilde, chi_head, z_tilde_a, z_head_a);
        let u = m2a.u(chi_tilde, chi_head, a_tilde, a_head);
    }
}
