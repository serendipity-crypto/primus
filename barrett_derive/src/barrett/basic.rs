use proc_macro2::TokenStream;
use quote::quote;
use syn::Ident;

pub(crate) fn basic(name: &Ident, modulus: &TokenStream) -> TokenStream {
    quote! {
        impl ::std::clone::Clone for #name {
            #[inline]
            fn clone(&self) -> Self {
                *self
            }
        }

        impl ::std::marker::Copy for #name {}

        impl ::std::fmt::Debug for #name {
            #[inline]
            fn fmt(&self, f: &mut ::std::fmt::Formatter<'_>) -> ::std::fmt::Result {
                write!(f, "{}", #modulus)
            }
        }

        impl ::std::hash::Hash for #name {
            #[inline]
            fn hash<H: ::std::hash::Hasher>(&self, state: &mut H){
                #modulus.hash(state)
            }
        }
    }
}
