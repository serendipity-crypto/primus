use std::fmt::Debug;

use integer::UnsignedInteger;

mod error;

pub mod lazy_ops;
pub mod ops;

pub use error::ReduceError;

use lazy_ops::*;
use ops::*;

/// Trait for types that represent a modulus.
pub trait Modulus: Copy {
    type ValueT;

    /// Returns the modulus value.
    fn value(self) -> Option<Self::ValueT>;

    /// Returns the modulus value without checking.
    fn value_unchecked(self) -> Self::ValueT;

    /// Returns the value of the modulus minus one.
    fn minus_one(self) -> Self::ValueT;
}

/// An trait indicate the modulus can perform operation like a ring.
pub trait RingAdapter<T>:
    Sized
    + Debug
    + Clone
    + Copy
    + Send
    + Sync
    + Modulus<ValueT = T>
    + Reduce<T, Output = T>
    + ReduceAssign<T>
    + ReduceAdd<T, Output = T>
    + ReduceAddAssign<T>
    + ReduceSub<T, Output = T>
    + ReduceSubAssign<T>
    + ReduceDouble<T, Output = T>
    + ReduceDoubleAssign<T>
    + ReduceNeg<T, Output = T>
    + ReduceNegAssign<T>
    + ReduceMul<T, Output = T>
    + ReduceMulAssign<T>
    + ReduceMulAdd<T, Output = T>
    + ReduceMulAddAssign<T>
    + ReduceSquare<T, Output = T>
    + ReduceSquareAssign<T>
    + ReduceExp<T>
    + ReduceExpPowOf2<T>
    + ReduceDotProduct<T>
    + TryReduceInv<T, Output = T>
{
}

impl<T: UnsignedInteger, M> RingAdapter<T> for M where
    M: Sized
        + Debug
        + Clone
        + Copy
        + Send
        + Sync
        + Modulus<ValueT = T>
        + Reduce<T, Output = T>
        + ReduceAssign<T>
        + ReduceAdd<T, Output = T>
        + ReduceAddAssign<T>
        + ReduceSub<T, Output = T>
        + ReduceSubAssign<T>
        + ReduceDouble<T, Output = T>
        + ReduceDoubleAssign<T>
        + ReduceNeg<T, Output = T>
        + ReduceNegAssign<T>
        + ReduceMul<T, Output = T>
        + ReduceMulAssign<T>
        + ReduceMulAdd<T, Output = T>
        + ReduceMulAddAssign<T>
        + ReduceSquare<T, Output = T>
        + ReduceSquareAssign<T>
        + ReduceExp<T>
        + ReduceExpPowOf2<T>
        + ReduceDotProduct<T>
        + TryReduceInv<T, Output = T>
{
}

/// An trait indicate the modulus can perform operation like a field.
pub trait FieldAdapter<T>:
    RingAdapter<T>
    + LazyReduce<T, Output = T>
    + LazyReduceAssign<T>
    + LazyReduceMul<T, Output = T>
    + LazyReduceMulAssign<T>
    + LazyReduceMulAdd<T, Output = T>
    + LazyReduceMulAddAssign<T>
    + for<'a> LazyReduce<&'a [T], Output = T>
    + for<'a> Reduce<&'a [T], Output = T>
    + ReduceInv<T, Output = T>
    + ReduceInvAssign<T>
    + ReduceDiv<T, Output = T>
    + ReduceDivAssign<T>
{
}

impl<T: UnsignedInteger, M> FieldAdapter<T> for M where
    M: RingAdapter<T>
        + LazyReduce<T, Output = T>
        + LazyReduceAssign<T>
        + LazyReduceMul<T, Output = T>
        + LazyReduceMulAssign<T>
        + LazyReduceMulAdd<T, Output = T>
        + LazyReduceMulAddAssign<T>
        + for<'a> LazyReduce<&'a [T], Output = T>
        + for<'a> Reduce<&'a [T], Output = T>
        + ReduceInv<T, Output = T>
        + ReduceInvAssign<T>
        + ReduceDiv<T, Output = T>
        + ReduceDivAssign<T>
{
}
