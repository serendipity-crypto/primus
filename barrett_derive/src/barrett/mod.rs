use proc_macro2::TokenStream;
use quote::ToTokens;
use syn::{DeriveInput, Result};

use crate::{BarrettModulusDeriveData, Modulus};

mod basic;
mod ops;
mod ratio;

#[inline]
pub(super) fn derive(input: &DeriveInput) -> Result<TokenStream> {
    let barrett_modulus_derive_data = BarrettModulusDeriveData::from_syn(input)?;

    Ok(impl_barrett(barrett_modulus_derive_data))
}

fn impl_barrett(data: BarrettModulusDeriveData) -> TokenStream {
    let name = &data.ident;

    let modulus = data.modulus.into_token_stream();

    let ratio = match data.modulus {
        Modulus::U8(m) => ratio::gen_ratio_u8(m).map(|v| v.into_token_stream()),
        Modulus::U16(m) => ratio::gen_ratio_u16(m).map(|v| v.into_token_stream()),
        Modulus::U32(m) => ratio::gen_ratio_u32(m).map(|v| v.into_token_stream()),
        Modulus::U64(m) => ratio::gen_ratio_u64(m).map(|v| v.into_token_stream()),
    };

    let ty = data.ty;

    let impl_basic = basic::basic(name, &modulus);

    let impl_lazy_reduce_ops = ops::impl_lazy_reduce_ops(name, &modulus, &ty, &ratio);

    let impl_add_ops = ops::impl_reduce_add_ops(name, &modulus, &ty);

    quote::quote! {
        #impl_basic

        #impl_lazy_reduce_ops

        #impl_add_ops
    }
}
