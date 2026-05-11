use std::fmt::Debug;

use primus_integer::UnsignedInteger;

use super::*;

/// A marker trait indicating the modulus can perform ring operations
/// (add, sub, neg, mul, square, exp, dot-product).
pub trait RingContext<T>:
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
{
}

impl<T: UnsignedInteger, M> RingContext<T> for M where
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
{
}

/// A marker trait indicating the modulus can perform field operations
/// (ring + lazy reduce, multiplicative inverse, division).
pub trait FieldContext<T>:
    RingContext<T>
    + LazyReduce<T, Output = T>
    + LazyReduceAssign<T>
    + LazyReduceMul<T, Output = T>
    + LazyReduceMulAssign<T>
    + LazyReduceMulAdd<T, Output = T>
    + LazyReduceMulAddAssign<T>
    + for<'a> LazyReduce<&'a [T], Output = T>
    + for<'a> Reduce<&'a [T], Output = T>
    + TryReduceInv<T, Output = T>
    + ReduceInv<T, Output = T>
    + ReduceInvAssign<T>
    + ReduceDiv<T, Output = T>
    + ReduceDivAssign<T>
{
}

impl<T: UnsignedInteger, M> FieldContext<T> for M where
    M: RingContext<T>
        + LazyReduce<T, Output = T>
        + LazyReduceAssign<T>
        + LazyReduceMul<T, Output = T>
        + LazyReduceMulAssign<T>
        + LazyReduceMulAdd<T, Output = T>
        + LazyReduceMulAddAssign<T>
        + for<'a> LazyReduce<&'a [T], Output = T>
        + for<'a> Reduce<&'a [T], Output = T>
        + TryReduceInv<T, Output = T>
        + ReduceInv<T, Output = T>
        + ReduceInvAssign<T>
        + ReduceDiv<T, Output = T>
        + ReduceDivAssign<T>
{
}
