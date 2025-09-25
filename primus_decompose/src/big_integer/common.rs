use std::iter::FusedIterator;

use integer::{BigIntegerOps, UnsignedInteger, izip};

pub struct BigUintScalar<T: UnsignedInteger>(pub Vec<T>);

/// An iterator over scalars.
pub struct BigUintScalarIter<T: UnsignedInteger> {
    scalar: Vec<T>,
    length: usize,
    log_basis: u32,
}

impl<T: UnsignedInteger> BigUintScalarIter<T> {
    /// Creates a new [`ScalarIter<T>`].
    #[inline]
    pub fn new(scalar: Vec<T>, length: usize, log_basis: u32) -> Self {
        Self {
            scalar,
            length,
            log_basis,
        }
    }
}

impl<T: UnsignedInteger> Iterator for BigUintScalarIter<T> {
    type Item = BigUintScalar<T>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if self.length == 0 {
            None
        } else {
            let next = self.scalar.clone();
            self.length -= 1;
            if self.length != 0 {
                self.scalar.slice_left_shift_assign(self.log_basis);
            }
            Some(BigUintScalar(next))
        }
    }
}

impl<T: UnsignedInteger> FusedIterator for BigUintScalarIter<T> {}
