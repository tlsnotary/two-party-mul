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

        for k in 0..L {
            let a_tilde = self.a_tilde();
            let a_head = self.a_head();

            for i in 0..ZETA {
                alpha[k * ZETA + i] = [a_tilde[i], a_head[i]];
            }
        }

        alpha
    }

    fn a_tilde(&mut self) -> [T; ZETA] {
        std::array::from_fn(|_| T::rand(&mut self.rng))
    }

    fn a_head(&mut self) -> [T; ZETA] {
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
        todo!();
    }

    fn u(chi_tilde: [T; L], chi_head: [T; L], a_tilde: [T; ZETA], a_head: [T; ZETA]) -> [T; L] {
        todo!();
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
