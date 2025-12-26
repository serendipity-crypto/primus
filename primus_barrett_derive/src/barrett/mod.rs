use proc_macro2::TokenStream;
use quote::ToTokens;
use syn::Result;

use crate::{BarrettModulusInput, Modulus};

mod ratio;

mod basic;

mod lazy_ops;
mod ops;

#[inline]
pub(super) fn derive(input: &BarrettModulusInput) -> Result<TokenStream> {
    let modulus = Modulus::from_syn(&input.value, &input.ty)?;
    modulus.check_leading_zeros(&input.value)?;

    Ok(impl_barrett(input, modulus))
}

fn impl_barrett(input: &BarrettModulusInput, modulus: Modulus) -> TokenStream {
    let vis = &input.vis;

    let name = &input.ident;

    let ty = &input.ty;

    let ratio = match modulus {
        Modulus::U8(m) => ratio::gen_ratio_u8(m).map(|v| v.into_token_stream()),
        Modulus::U16(m) => ratio::gen_ratio_u16(m).map(|v| v.into_token_stream()),
        Modulus::U32(m) => ratio::gen_ratio_u32(m).map(|v| v.into_token_stream()),
        Modulus::U64(m) => ratio::gen_ratio_u64(m).map(|v| v.into_token_stream()),
    };

    let modulus = modulus.into_token_stream();

    let impl_basic = basic::basic(vis, name, &modulus, ty, &ratio);

    let lazy_reduce_ops = lazy_ops::impl_lazy_reduce_ops(name, &modulus, ty, &ratio);

    let reduce_ops = ops::impl_reduce_ops(name, &modulus, ty);

    quote::quote! {
        #impl_basic

        #lazy_reduce_ops

        #reduce_ops
    }
}
