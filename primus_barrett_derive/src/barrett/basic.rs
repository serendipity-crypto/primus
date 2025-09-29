use proc_macro2::TokenStream;
use quote::quote;
use syn::Ident;

pub(crate) fn basic(
    name: &Ident,
    modulus: &TokenStream,
    ty: &syn::Type,
    ratio: &[TokenStream; 2],
) -> TokenStream {
    let [r0, r1] = ratio;
    quote! {
        impl #name {
            /// Retures the modulus value.
            pub const fn value() -> #ty {
                #modulus
            }

            /// Returns the ratio of this [`BarrettModulus<T>`].
            #[inline]
            pub const fn ratio() -> [#ty; 2] {
                [#r0, #r1]
            }
        }

        impl ::primus_barrett_modulus::reduce::Modulus for #name {
            type ValueT = #ty;

            #[inline(always)]
            fn value(self) -> Option<Self::ValueT> {
                Some(#modulus)
            }

            #[inline(always)]
            fn value_unchecked(self) -> Self::ValueT {
                #modulus
            }

            #[inline(always)]
            fn minus_one(self) -> Self::ValueT {
                #modulus - 1
            }
        }

        impl ::std::marker::Copy for #name {}

        impl ::std::clone::Clone for #name {
            #[inline]
            fn clone(&self) -> Self {
                *self
            }
        }

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
