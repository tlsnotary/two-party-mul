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
    gadget: [T; ZETA],
}

impl<T: Field> M2A<T> {
    // Note that the return dimension is [[[T; 2]; ZETA]; L] which we simplify a little bit
    fn alpha(&mut self) -> [[T; 2]; ETA] {
        let mut alpha = [[T::zero(); 2]; ETA];

        let a_tilde = self.a_tilde();
        let a_head = self.a_head();

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

    fn b_tilde(gadget: [T; ZETA], beta: [T; ETA]) -> [T; L] {
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

    fn chi_head() -> [T; L] {
        std::array::from_fn(|_| T::rand(&mut thread_rng()))
    }

    fn chi_tilde() -> [T; L] {
        std::array::from_fn(|_| T::rand(&mut thread_rng()))
    }

    fn r(
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

    fn u(chi_tilde: [T; L], chi_head: [T; L], a_tilde: [T; L], a_head: [T; L]) -> [T; L] {
        let mut u = [T::zero(); L];

        u.iter_mut()
            .enumerate()
            .for_each(|(k, el)| *el = chi_tilde[k] * a_tilde[k] + chi_head[k] * a_head[k]);

        u
    }

    fn bob_check(
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

    fn gamma(input: [T; L], tilde: [T; L]) -> [T; L] {
        let mut gamma = [T::zero(); L];

        gamma
            .iter_mut()
            .enumerate()
            .for_each(|(k, el)| *el = input[k] + -tilde[k]);

        gamma
    }

    //TODO: Is there a mistake in the paper with b_tilde ?
    fn output(input: [T; L], gamma: [T; L], gadget: [T; ZETA], z_tilde: [T; ETA]) -> [T; L] {
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
        let gadget = std::array::from_fn(|_| T::rand(&mut rng));

        Self {
            _field: PhantomData,
            rng,
            gadget,
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
