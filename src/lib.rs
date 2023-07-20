//! Implements the two-party multiplication protocol from here: <https://eprint.iacr.org/2019/523>

use mpz_ot::{ObliviousReceive, ObliviousSend};
use mpz_share_conversion_core::Field;
use rand::{rngs::ThreadRng, thread_rng, CryptoRng, RngCore};
use std::marker::PhantomData;

pub struct M2A<T: Field, U, V: CryptoRng + RngCore = ThreadRng> {
    _field: PhantomData<T>,
    _role: PhantomData<U>,
    rng: V,
}

impl<T: Field, U: ObliviousSend<T>, V: CryptoRng + RngCore> M2A<T, U, V> {}

impl<T: Field, U: ObliviousReceive<bool, T>, V: CryptoRng + RngCore> M2A<T, U, V> {
    pub fn beta(&mut self) -> T {
        todo!()
    }

    pub fn b_tilde() -> T {
        todo!()
    }
}

impl<T: Field, U, V: CryptoRng + RngCore> M2A<T, U, V> {
    // Bits needed for to represent elements of the field
    const KAPPA: u32 = T::BIT_SIZE;

    // Security parameter
    const S: u32 = 80;

    // Redundancy factor
    const ZETA: u32 = Self::KAPPA + 2 * Self::S;

    // Batch size
    const L: u32 = 1;

    // Number of OTs
    const ETA: u32 = Self::ZETA * Self::L;

    fn gadget(&self) -> T {
        todo!()
    }
}

impl<T: Field, U> Default for M2A<T, U, ThreadRng> {
    fn default() -> Self {
        let rng = thread_rng();
        Self {
            _field: PhantomData,
            _role: PhantomData,
            rng,
        }
    }
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
