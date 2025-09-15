use proc_macro2::TokenStream;
use quote::quote;
use syn::Ident;

pub(crate) fn impl_lazy_reduce_ops(
    name: &Ident,
    modulus: &TokenStream,
    ty: &syn::Type,
    ratio: &[TokenStream; 2],
) -> TokenStream {
    let [r0, r1] = ratio;
    quote! {
        impl ::reduce::lazy_ops::LazyReduce<#ty> for #name {
            type Output = #ty;

            /// Calculates `value (mod 2*modulus)`.
            #[inline]
            fn lazy_reduce(self, value: #ty) -> #ty {
                use ::integer::{CarryingMul, WideningMul};
                // Step 1.
                //              ratio[1]  ratio[0]
                //         *               value
                //   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                //            +-------------------+
                //            |  tmp1   |         |    <-- value * ratio[0]
                //            +-------------------+
                //   +------------------+
                //   |      tmp2        |              <-- value * ratio[1]
                //   +------------------+
                //   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                //   +--------+
                //   |   q₃   |
                //   +--------+
                let tmp = value.widening_mul_hw(#r0); // tmp1
                let q = value.carrying_mul_hw(#r1, tmp); // q₃

                // Step 2.
                value.wrapping_sub(q.wrapping_mul(#modulus)) // r = r₁ - r₂
            }
        }

        impl ::reduce::lazy_ops::LazyReduce<[#ty; 2]> for #name {
            type Output = #ty;

            /// Calculates `value (mod 2*modulus)`.
            #[inline]
            fn lazy_reduce(self, value: [#ty; 2]) -> Self::Output {
                use ::integer::{CarryingMul, WideningMul};
                // Step 1.
                //                        ratio[1]  ratio[0]
                //                   *    value[1]  value[0]
                //   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                //                      +-------------------+
                //                      |         a         |    <-- value[0] * ratio[0]
                //                      +-------------------+
                //             +------------------+
                //             |        b         |              <-- value[0] * ratio[1]
                //             +------------------+
                //             +------------------+
                //             |        c         |              <-- value[1] * ratio[0]
                //             +------------------+
                //   +------------------+
                //   |        d         |                        <-- value[1] * ratio[1]
                //   +------------------+
                //   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                //             +--------+
                //             |   q₃   |
                //             +--------+
                let ah = value[0].widening_mul_hw(#r0);

                let b = CarryingMul::carrying_mul(value[0], #r1, ah);
                let c = WideningMul::widening_mul(value[1], #r0);

                let d = value[1].wrapping_mul(#r1);

                let bch = b.1 + c.1 + b.0.overflowing_add(c.0).1 as #ty;

                let q = d.wrapping_add(bch);

                // Step 2.
                value[0].wrapping_sub(q.wrapping_mul(#modulus))
            }
        }

        impl ::reduce::lazy_ops::LazyReduce<(#ty, #ty)> for #name {
            type Output = #ty;

            /// Calculates `value (mod 2*modulus)`.
            #[inline]
            fn lazy_reduce(self, value: (#ty, #ty)) -> Self::Output {
                use ::integer::{CarryingMul, WideningMul};
                // Step 1.
                //                        ratio[1]  ratio[0]
                //                   *    value.1   value.0
                //   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                //                      +-------------------+
                //                      |         a         |    <-- value.0 * ratio[0]
                //                      +-------------------+
                //             +------------------+
                //             |        b         |              <-- value.0 * ratio[1]
                //             +------------------+
                //             +------------------+
                //             |        c         |              <-- value.1 * ratio[0]
                //             +------------------+
                //   +------------------+
                //   |        d         |                        <-- value.1 * ratio[1]
                //   +------------------+
                //   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                //             +--------+
                //             |   q₃   |
                //             +--------+
                let ah = value.0.widening_mul_hw(#r0);

                let b = CarryingMul::carrying_mul(value.0, #r1, ah);
                let c = WideningMul::widening_mul(value.1, #r0);

                let d = value.1.wrapping_mul(#r1);

                let bch = b.1 + c.1 + b.0.overflowing_add(c.0).1 as #ty;

                let q = d.wrapping_add(bch);

                // Step 2.
                value.0.wrapping_sub(q.wrapping_mul(#modulus))
            }
        }

        impl ::reduce::lazy_ops::LazyReduce<&[#ty]> for #name {
            type Output = #ty;

            /// Calculates `value (mod 2*modulus)` when value's length > 0.
            #[inline]
            fn lazy_reduce(self, value: &[#ty]) -> Self::Output {
                match value {
                    &[] => unreachable!(),
                    &[v] => {
                        if v < #modulus << 1u32 {
                            v
                        } else {
                            self.lazy_reduce(v)
                        }
                    }
                    [other @ .., last] => other
                        .iter()
                        .rfold(*last, |acc, &x| self.lazy_reduce([x, acc])),
                }
            }
        }

        impl ::reduce::lazy_ops::LazyReduceAssign<#ty> for #name {
            /// Calculates `value (mod 2*modulus)`.
            #[inline]
            fn lazy_reduce_assign(self, value: &mut #ty) {
                use ::reduce::lazy_ops::LazyReduce;
                *value = self.lazy_reduce(*value);
            }
        }
    }
}

pub(crate) fn impl_reduce_add_ops(
    name: &Ident,
    modulus: &TokenStream,
    ty: &syn::Type,
) -> TokenStream {
    quote! {
        impl ::reduce::ops::ReduceAdd<#ty> for #name {
            type Output = #ty;

            #[inline(always)]
            fn reduce_add(self, a: #ty, b: #ty) -> Self::Output {
                ::uint_modulus::UintModulus(#modulus).reduce_add(a, b)
            }
        }
    }
}
