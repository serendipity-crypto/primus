use proc_macro2::TokenStream;
use quote::quote;
use syn::Ident;

pub(crate) fn impl_lazy_reduce_ops(
    name: &Ident,
    modulus: &TokenStream,
    ty: &syn::Path,
    ratio: &[TokenStream; 2],
) -> TokenStream {
    let [r0, r1] = ratio;
    quote! {
        impl ::primus_modulus::reduce::lazy_ops::LazyReduce<#ty> for #name {
            type Output = #ty;

            /// Calculates `value (mod 2*modulus)`.
            #[inline]
            fn lazy_reduce(self, value: #ty) -> #ty {
                use ::primus_modulus::integer::{CarryingMul, WideningMul};
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

        impl ::primus_modulus::reduce::lazy_ops::LazyReduce<[#ty; 2]> for #name {
            type Output = #ty;

            /// Calculates `value (mod 2*modulus)`.
            #[inline]
            fn lazy_reduce(self, value: [#ty; 2]) -> Self::Output {
                use ::primus_modulus::integer::{CarryingMul, WideningMul};
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

        impl ::primus_modulus::reduce::lazy_ops::LazyReduce<(#ty, #ty)> for #name {
            type Output = #ty;

            /// Calculates `value (mod 2*modulus)`.
            #[inline]
            fn lazy_reduce(self, value: (#ty, #ty)) -> Self::Output {
                use ::primus_modulus::integer::{CarryingMul, WideningMul};
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

        impl ::primus_modulus::reduce::lazy_ops::LazyReduce<&[#ty]> for #name {
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

        impl ::primus_modulus::reduce::lazy_ops::LazyReduceAssign<#ty> for #name {
            /// Calculates `value (mod 2*modulus)`.
            #[inline]
            fn lazy_reduce_assign(self, value: &mut #ty) {
                use ::primus_modulus::reduce::lazy_ops::LazyReduce;
                *value = self.lazy_reduce(*value);
            }
        }

        impl ::primus_modulus::reduce::lazy_ops::LazyReduceMul<#ty> for #name {
            type Output = #ty;

            #[inline]
            fn lazy_reduce_mul(self, a: #ty, b: #ty) -> Self::Output {
                use ::primus_modulus::reduce::lazy_ops::LazyReduce;
                use ::primus_modulus::integer::WideningMul;
                self.lazy_reduce(WideningMul::widening_mul(a, b))
            }
        }

        impl ::primus_modulus::reduce::lazy_ops::LazyReduceMulAssign<#ty> for #name {
            #[inline]
            fn lazy_reduce_mul_assign(self, a: &mut #ty, b: #ty) {
                use ::primus_modulus::reduce::lazy_ops::LazyReduce;
                use ::primus_modulus::integer::WideningMul;
                *a = self.lazy_reduce(WideningMul::widening_mul(*a, b));
            }
        }

        impl ::primus_modulus::reduce::lazy_ops::LazyReduceMulAdd<#ty> for #name {
            type Output = #ty;

            #[inline]
            fn lazy_reduce_mul_add(self, a: #ty, b: #ty, c: #ty) -> Self::Output {
                use ::primus_modulus::reduce::lazy_ops::LazyReduce;
                use ::primus_modulus::integer::CarryingMul;
                self.lazy_reduce(CarryingMul::carrying_mul(a, b, c))
            }
        }

        impl ::primus_modulus::reduce::lazy_ops::LazyReduceMulAddAssign<#ty> for #name {
            #[inline]
            fn lazy_reduce_mul_add_assign(self, a: &mut #ty, b: #ty, c: #ty) {
                use ::primus_modulus::reduce::lazy_ops::LazyReduce;
                use ::primus_modulus::integer::CarryingMul;
                *a = self.lazy_reduce(CarryingMul::carrying_mul(*a, b, c));
            }
        }
    }
}
