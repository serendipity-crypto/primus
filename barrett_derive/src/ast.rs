use syn::{DeriveInput, Meta, Token, punctuated::Punctuated};

use crate::Modulus;

pub(crate) struct BarrettModulusDeriveData {
    pub(crate) ident: syn::Ident,
    pub(crate) ty: syn::Type,
    pub(crate) modulus: Modulus,
}

impl BarrettModulusDeriveData {
    pub(crate) fn from_syn(input: &DeriveInput) -> syn::Result<Self> {
        let mut ty = None;
        let mut modulus_lit = None;

        for attr in &input.attrs {
            if attr.path().is_ident("modulus") {
                let nested =
                    attr.parse_args_with(Punctuated::<Meta, Token![,]>::parse_terminated)?;
                for meta in nested {
                    if meta.path().is_ident("value")
                        && let Meta::NameValue(nv) = &meta
                        && let syn::Expr::Lit(expr) = &nv.value
                        && let syn::Lit::Int(lit_int) = &expr.lit
                    {
                        modulus_lit = Some(lit_int.clone());
                    } else {
                        ty = Some(meta.path().clone());
                    }
                }
            }
        }

        if let (Some(modulus_lit), Some(ty)) = (modulus_lit, ty) {
            let modulus = Modulus::from_syn(&modulus_lit, &ty)?;
            modulus.check_leading_zeros(&modulus_lit)?;
            let ty = modulus.ty();
            return Ok(Self {
                ident: input.ident.clone(),
                modulus,
                ty,
            });
        }

        Err(syn::Error::new_spanned(
            input,
            "modulus should be in the form of `#[modulus(u8|u16|u32|u64, value = <value>)]`",
        ))
    }
}
