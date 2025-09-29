use super::CarryingAdd;

macro_rules! impl_uint_carrying_add {
    ($($T:ty),*) => {
        $(
            impl CarryingAdd for $T {
                type CarryT = bool;

                #[inline]
                fn carrying_add(self, rhs: Self, carry: Self::CarryT) -> (Self, Self::CarryT) {
                    let (a, b) = self.overflowing_add(rhs);
                    let (c, d) = a.overflowing_add(carry as Self);
                    (c, b || d)
                }
            }
        )*
    };
}

impl_uint_carrying_add! {u8, u16, u32, u64, u128, usize}
