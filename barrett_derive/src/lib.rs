use proc_macro::TokenStream;
use syn::{DeriveInput, parse_macro_input};

mod ast;
mod barrett;
mod modulus;

pub(crate) use ast::BarrettModulusDeriveData;
pub(crate) use modulus::Modulus;

#[proc_macro_derive(Barrett, attributes(modulus, value_type))]
pub fn derive_barrett(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);

    barrett::derive(&input)
        .unwrap_or_else(|err| err.to_compile_error())
        .into()
}
