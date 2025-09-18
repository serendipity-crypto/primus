use proc_macro2::TokenStream;
use quote::quote;
use syn::Ident;

pub(crate) fn impl_reduce_ops(name: &Ident, modulus: &TokenStream, ty: &syn::Type) -> TokenStream {
    quote! {
        impl ::barrett::reduce::ops::Reduce<#ty> for #name {
            type Output = #ty;

            /// Calculates `value (mod modulus)`.
            #[inline(always)]
            fn reduce(self, value: #ty) -> Self::Output {
                use ::barrett::reduce::ops::ReduceOnce;
                use ::barrett::reduce::lazy_ops::LazyReduce;
                ::barrett::uint_modulus::UintModulus(#modulus).reduce_once(self.lazy_reduce(value))
            }
        }

        impl ::barrett::reduce::ops::Reduce<[#ty; 2]> for #name {
            type Output = #ty;

            /// Calculates `value (mod modulus)`.
            #[inline(always)]
            fn reduce(self, value: [#ty; 2]) -> Self::Output {
                use ::barrett::reduce::ops::ReduceOnce;
                use ::barrett::reduce::lazy_ops::LazyReduce;
                ::barrett::uint_modulus::UintModulus(#modulus).reduce_once(self.lazy_reduce(value))
            }
        }

        impl ::barrett::reduce::ops::Reduce<(#ty, #ty)> for #name {
            type Output = #ty;

            /// Calculates `value (mod modulus)`.
            #[inline(always)]
            fn reduce(self, value: (#ty, #ty)) -> Self::Output {
                use ::barrett::reduce::ops::ReduceOnce;
                use ::barrett::reduce::lazy_ops::LazyReduce;
                ::barrett::uint_modulus::UintModulus(#modulus).reduce_once(self.lazy_reduce(value))
            }
        }

        impl ::barrett::reduce::ops::Reduce<&[#ty]> for #name {
            type Output = #ty;

            /// Calculates `value (mod modulus)` when value's length > 0.
            #[inline(always)]
            fn reduce(self, value: &[#ty]) -> Self::Output {
                use ::barrett::reduce::ops::ReduceOnce;
                use ::barrett::reduce::lazy_ops::LazyReduce;
                ::barrett::uint_modulus::UintModulus(#modulus).reduce_once(self.lazy_reduce(value))
            }
        }

        impl ::barrett::reduce::ops::ReduceAssign<#ty> for #name {
            /// Calculates `value (mod modulus)`.
            #[inline]
            fn reduce_assign(self, value: &mut #ty) {
                use ::barrett::reduce::ops::Reduce;
                *value = self.reduce(*value);
            }
        }

        impl ::barrett::reduce::ops::ReduceOnce<#ty> for #name {
            type Output = #ty;

            #[inline(always)]
            fn reduce_once(self, value: #ty) -> Self::Output {
                use ::barrett::reduce::ops::ReduceOnce;
                ::barrett::uint_modulus::UintModulus(#modulus).reduce_once(value)
            }
        }

        impl ::barrett::reduce::ops::ReduceOnceAssign<#ty> for #name {
            #[inline(always)]
            fn reduce_once_assign(self, value: &mut #ty) {
                ::barrett::uint_modulus::UintModulus(#modulus).reduce_once_assign(value);
            }
        }

        impl ::barrett::reduce::ops::ReduceAdd<#ty> for #name {
            type Output = #ty;

            #[inline(always)]
            fn reduce_add(self, a: #ty, b: #ty) -> Self::Output {
                ::barrett::uint_modulus::UintModulus(#modulus).reduce_add(a, b)
            }
        }

        impl ::barrett::reduce::ops::ReduceAddAssign<#ty> for #name {
            #[inline(always)]
            fn reduce_add_assign(self, a: &mut #ty, b: #ty) {
                ::barrett::uint_modulus::UintModulus(#modulus).reduce_add_assign(a, b);
            }
        }

        impl ::barrett::reduce::ops::ReduceDouble<#ty> for #name {
            type Output = #ty;

            #[inline(always)]
            fn reduce_double(self, value: #ty) -> Self::Output {
                ::barrett::uint_modulus::UintModulus(#modulus).reduce_double(value)
            }
        }

        impl ::barrett::reduce::ops::ReduceDoubleAssign<#ty> for #name {
            #[inline(always)]
            fn reduce_double_assign(self, value: &mut #ty) {
                ::barrett::uint_modulus::UintModulus(#modulus).reduce_double_assign(value);
            }
        }

        impl ::barrett::reduce::ops::ReduceSub<#ty> for #name {
            type Output = #ty;

            #[inline(always)]
            fn reduce_sub(self, a: #ty, b: #ty) -> Self::Output {
                ::barrett::uint_modulus::UintModulus(#modulus).reduce_sub(a, b)
            }
        }

        impl ::barrett::reduce::ops::ReduceSubAssign<#ty> for #name {
            #[inline(always)]
            fn reduce_sub_assign(self, a: &mut #ty, b: #ty) {
                ::barrett::uint_modulus::UintModulus(#modulus).reduce_sub_assign(a, b);
            }
        }

        impl ::barrett::reduce::ops::ReduceNeg<#ty> for #name {
            type Output = #ty;

            #[inline(always)]
            fn reduce_neg(self, value: #ty) -> Self::Output {
                ::barrett::uint_modulus::UintModulus(#modulus).reduce_neg(value)
            }
        }

        impl ::barrett::reduce::ops::ReduceNegAssign<#ty> for #name {
            #[inline(always)]
            fn reduce_neg_assign(self, value: &mut #ty) {
                ::barrett::uint_modulus::UintModulus(#modulus).reduce_neg_assign(value);
            }
        }

        impl ::barrett::reduce::ops::ReduceMul<#ty> for #name {
            type Output = #ty;

            #[inline]
            fn reduce_mul(self, a: #ty, b: #ty) -> Self::Output {
                use ::barrett::reduce::ops::Reduce;
                use ::barrett::integer::WideningMul;
                self.reduce(WideningMul::widening_mul(a, b))
            }
        }

        impl ::barrett::reduce::ops::ReduceMulAssign<#ty> for #name {
            #[inline]
            fn reduce_mul_assign(self, a: &mut #ty, b: #ty) {
                use ::barrett::reduce::ops::Reduce;
                use ::barrett::integer::WideningMul;
                *a = self.reduce(WideningMul::widening_mul(*a, b));
            }
        }

        impl ::barrett::reduce::ops::ReduceSquare<#ty> for #name {
            type Output = #ty;

            #[inline]
            fn reduce_square(self, value: #ty) -> Self::Output {
                use ::barrett::reduce::ops::Reduce;
                use ::barrett::integer::WideningMul;
                self.reduce(WideningMul::widening_mul(value, value))
            }
        }

        impl ::barrett::reduce::ops::ReduceSquareAssign<#ty> for #name {
            #[inline]
            fn reduce_square_assign(self, value: &mut #ty) {
                use ::barrett::reduce::ops::Reduce;
                use ::barrett::integer::WideningMul;
                *value = self.reduce(WideningMul::widening_mul(*value, *value));
            }
        }

        impl ::barrett::reduce::ops::ReduceMulAdd<#ty> for #name {
            type Output = #ty;

            #[inline]
            fn reduce_mul_add(self, a: #ty, b: #ty, c: #ty) -> Self::Output {
                use ::barrett::reduce::ops::Reduce;
                use ::barrett::integer::CarryingMul;
                self.reduce(a.carrying_mul(b, c))
            }
        }

        impl ::barrett::reduce::ops::ReduceMulAddAssign<#ty> for #name {
            #[inline]
            fn reduce_mul_add_assign(self, a: &mut #ty, b: #ty, c: #ty) {
                use ::barrett::reduce::ops::Reduce;
                use ::barrett::integer::CarryingMul;
                *a = self.reduce(a.carrying_mul(b, c));
            }
        }

        impl ::barrett::reduce::ops::ReduceInv<#ty> for #name {
            type Output = #ty;

            #[inline(always)]
            fn reduce_inv(self, value: #ty) -> Self::Output {
                ::barrett::uint_modulus::UintModulus(#modulus).reduce_inv(value)
            }
        }

        impl ::barrett::reduce::ops::ReduceInvAssign<#ty> for #name {
            #[inline(always)]
            fn reduce_inv_assign(self, value: &mut #ty) {
                ::barrett::uint_modulus::UintModulus(#modulus).reduce_inv_assign(value);
            }
        }

        impl ::barrett::reduce::ops::ReduceDiv<#ty> for #name {
            type Output = #ty;

            #[inline]
            fn reduce_div(self, a: #ty, b: #ty) -> Self::Output {
                use ::barrett::reduce::ops::{ReduceMul, ReduceInv};
                self.reduce_mul(a, self.reduce_inv(b))
            }
        }

        impl ::barrett::reduce::ops::ReduceDivAssign<#ty> for #name {
            #[inline]
            fn reduce_div_assign(self, a: &mut #ty, b: #ty) {
                use ::barrett::reduce::ops::{ReduceMulAssign, ReduceInv};
                self.reduce_mul_assign(a, self.reduce_inv(b));
            }
        }

        impl ::barrett::reduce::ops::ReduceExp<#ty> for #name{
            #[inline]
            fn reduce_exp<E: ::barrett::integer::UnsignedInteger>(self, base: #ty, mut exp: E) -> #ty {
                use ::barrett::reduce::ops::{ReduceSquareAssign, ReduceMulAssign};
                if exp.is_zero() {
                    return 1;
                }

                if base == 0 {
                    return 0;
                }

                debug_assert!(base < #modulus);

                let mut power: #ty = base;

                let exp_trailing_zeros = exp.trailing_zeros();
                if exp_trailing_zeros > 0 {
                    for _ in 0..exp_trailing_zeros {
                        self.reduce_square_assign(&mut power);
                    }
                    exp >>= exp_trailing_zeros;
                }

                if exp.is_one() {
                    return power;
                }

                let mut intermediate: #ty = power;
                for _ in 1..(E::BITS - exp.leading_zeros()) {
                    exp >>= 1;
                    self.reduce_square_assign(&mut power);
                    if !(exp & E::ONE).is_zero() {
                        self.reduce_mul_assign(&mut intermediate, power);
                    }
                }
                intermediate
            }
        }

        impl ::barrett::reduce::ops::ReduceExpPowOf2<#ty> for #name {
            #[inline]
            fn reduce_exp_power_of_2(self, base: #ty, exp_log: u32) -> #ty {
                use ::barrett::reduce::ops::ReduceSquareAssign;
                if base == 0 {
                    return 0;
                }

                let mut power = base;

                for _ in 0..exp_log {
                    self.reduce_square_assign(&mut power);
                }

                power
            }
        }

        impl ::barrett::reduce::ops::ReduceDotProduct<#ty> for #name {
            #[inline]
            fn reduce_dot_product(self, a: impl AsRef<[#ty]>, b: impl AsRef<[#ty]>) -> #ty {
                use ::barrett::reduce::ops::*;
                /// `c += a * b`
                fn multiply_add(c: &mut [#ty; 2], a: #ty, b: #ty) {
                    use ::barrett::integer::WideningMul;
                    let (lw, hw) = a.widening_mul(b);
                    let carry;
                    (c[0], carry) = c[0].overflowing_add(lw);
                    (c[1], _) = c[1].carrying_add(hw, carry);
                }

                let a = a.as_ref();
                let b = b.as_ref();

                debug_assert_eq!(a.len(), b.len());

                let mut a_iter = a.chunks_exact(16);
                let mut b_iter = b.chunks_exact(16);

                let inter = (&mut a_iter)
                    .zip(&mut b_iter)
                    .map(|(a_s, b_s)| {
                        let mut c: [#ty; 2] = [0, 0];
                        for (&a, &b) in a_s.iter().zip(b_s) {
                            multiply_add(&mut c, a, b);
                        }
                        self.reduce(c)
                    })
                    .fold(0, |acc: #ty, b| self.reduce_add(acc, b));

                let mut c: [#ty; 2] = [0, 0];
                a_iter
                    .remainder()
                    .iter()
                    .zip(b_iter.remainder())
                    .for_each(|(&a, &b)| {
                        multiply_add(&mut c, a, b);
                    });
                self.reduce_add(self.reduce(c), inter)
            }

            #[inline]
            fn reduce_dot_product_iter(
                self,
                a: impl IntoIterator<Item = #ty>,
                b: impl IntoIterator<Item = #ty>,
            ) -> #ty {
                use ::barrett::reduce::ops::*;
                /// `c += a * b`
                fn multiply_add(c: &mut [#ty; 2], a: #ty, b: #ty) {
                    use ::barrett::integer::WideningMul;
                    let (lw, hw) = a.widening_mul(b);
                    let carry;
                    (c[0], carry) = c[0].overflowing_add(lw);
                    (c[1], _) = c[1].carrying_add(hw, carry);
                }

                let mut a_iter = a.into_iter();
                let mut b_iter = b.into_iter();

                let mut a_temp_array = [0; 16];
                let mut b_temp_array = [0; 16];

                let mut i = 0;
                let mut result = 0;

                while let (Some(a_next), Some(b_next)) = (a_iter.next(), b_iter.next()) {
                    if i < 16 {
                        a_temp_array[i] = a_next;
                        b_temp_array[i] = b_next;
                        i += 1;
                    } else {
                        let mut c: [#ty; 2] = [0, 0];
                        for (&a, b) in a_temp_array.iter().zip(b_temp_array) {
                            multiply_add(&mut c, a, b);
                        }
                        self.reduce_add_assign(&mut result, self.reduce(c));

                        a_temp_array.fill(0);
                        b_temp_array.fill(0);
                        a_temp_array[0] = a_next;
                        b_temp_array[0] = b_next;
                        i = 1;
                    }
                }

                let mut c: [#ty; 2] = [0, 0];
                for (&a, &b) in a_temp_array[..i].iter().zip(b_temp_array[..i].iter()) {
                    multiply_add(&mut c, a, b);
                }
                self.reduce_add_assign(&mut result, self.reduce(c));

                result
            }
        }
    }
}
