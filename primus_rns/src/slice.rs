use integer::UnsignedInteger;

pub trait BigIntOps {
    type ValueT;
    fn add_value_inplace(&mut self, value: Self::ValueT) -> bool;
    fn sub_value_inplace(&mut self, value: Self::ValueT) -> bool;
    fn mul_value_inplace(&mut self, value: Self::ValueT) -> Self::ValueT;
}

impl<T: UnsignedInteger> BigIntOps for [T] {
    type ValueT = T;

    fn add_value_inplace(&mut self, value: Self::ValueT) -> bool {
        let mut carry;
        match self {
            [first, other @ ..] => {
                (*first, carry) = first.overflowing_add(value);
                for v in other.iter_mut() {
                    if !carry {
                        return false;
                    }
                    (*v, carry) = (*v).overflowing_add(T::ONE);
                }
                carry
            }
            _ => unreachable!(),
        }
    }

    fn sub_value_inplace(&mut self, value: Self::ValueT) -> bool {
        let mut borrow;
        match self {
            [first, other @ ..] => {
                (*first, borrow) = first.overflowing_sub(value);
                for v in other.iter_mut() {
                    if !borrow {
                        return false;
                    }
                    (*v, borrow) = (*v).overflowing_sub(T::ONE);
                }
                borrow
            }
            _ => unreachable!(),
        }
    }

    fn mul_value_inplace(&mut self, value: Self::ValueT) -> Self::ValueT {
        if value.is_zero() {
            self.fill(T::ZERO);
            return T::ZERO;
        }

        let mut carry = T::ZERO;
        for ele in self.iter_mut() {
            let (low, high) = value.carrying_mul(*ele, carry);
            *ele = low;
            carry = high;
        }

        carry
    }
}
