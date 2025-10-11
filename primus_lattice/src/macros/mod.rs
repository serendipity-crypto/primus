// macro_rules! foo {
//     (@$trait_:ident [$($args:ident,)*] where [$($preds:tt)+]) => {
//         foo! {
//             @as_item
//             impl<$($args),*> $trait_<$($args),*>
//                 where $($args: ::std::any::Any + 'static,)*
//                       $($preds)*
//             {
//                 #[allow(non_camel_case_types, dead_code)]
//                 pub fn bar<__foo_T: $trait_<$($args),*>>(&self) {}
//             }
//         }
//     };
//     (@as_item $i:item) => { $i };

//     (
//         $trait_:ident < $($args:ident),* $(,)* >
//         where $($preds:tt)+
//     ) => {
//         foo! { @$trait_ [$($args,)*] where [$($preds)*] }
//     };
// }

macro_rules! impl_bytes_conversion {
    ($cipher:ident < $s:ident, $t:ident >) => {
        impl<$s, $t> $cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + DataOwned,
            $t: UnsignedInteger,
        {
            /// Creates from bytes `data`.
            #[inline]
            pub fn from_bytes(data: &[u8]) -> Self {
                let converted_data: &[$t] = bytemuck::cast_slice(data);

                Self {
                    data: ArrayBase::from_slice(converted_data),
                }
            }
        }

        impl<$s, $t> $cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + DataMut,
            $t: UnsignedInteger,
        {
            /// Copy from bytes `data`.
            #[inline]
            pub fn from_bytes_assign(&mut self, data: &[u8]) {
                let converted_data: &[$t] = bytemuck::cast_slice(data);

                self.data.copy_from_slice(converted_data);
            }
        }

        impl<$s, $t> $cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + Data,
            $t: UnsignedInteger,
        {
            /// Converts `self` into bytes.
            #[inline]
            pub fn to_bytes(&self) -> Vec<u8> {
                let converted_data: &[u8] = bytemuck::cast_slice(self.data.as_ref());

                converted_data.to_vec()
            }

            /// Converts `self` into bytes, stored in `data`.
            #[inline]
            pub fn to_bytes_inplace(&self, data: &mut [u8]) {
                let converted_data: &[u8] = bytemuck::cast_slice(self.data.as_ref());

                data.copy_from_slice(converted_data);
            }

            /// Returns the bytes count.
            #[inline]
            pub fn bytes_count(&self) -> usize {
                self.data.len() * <$t>::BYTES_COUNT
            }
        }
    };
}

// macro_rules! impl_temp {
//     ($cipher:ident < $s:ident, $t:ident >) => {
//         impl<$s, $t> $cipher<$s, $t>
//         where
//             $s: RawData<Elem = $t> + DataOwned,
//             $t: UnsignedInteger,
//         {
//         }

//         impl<$s, $t> $cipher<$s, $t>
//         where
//             $s: RawData<Elem = $t> + DataMut,
//             $t: UnsignedInteger,
//         {
//         }

//         impl<$s, $t> $cipher<$s, $t>
//         where
//             $s: RawData<Elem = $t> + Data,
//             $t: UnsignedInteger,
//         {
//         }
//     };
// }

macro_rules! impl_basic_operation_single_modulus {
    ($cipher:ident < $s:ident, $t:ident >) => {
        impl<$s, $t> $cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + DataOwned,
            $t: UnsignedInteger,
        {
            /// Perform element-wise modular addition `self + rhs`.
            #[inline]
            pub fn add_element_wise<M, A>(mut self, rhs: &$cipher<A>, modulus: M) -> Self
            where
                M: FieldContext<$t>,
                A: RawData<Elem = $t> + Data,
            {
                self.data.add_assign(&rhs.data, modulus);
                self
            }

            /// Perform element-wise modular subtraction `self - rhs`.
            #[inline]
            pub fn sub_element_wise<M, A>(mut self, rhs: &$cipher<A>, modulus: M) -> Self
            where
                M: FieldContext<$t>,
                A: RawData<Elem = $t> + Data,
            {
                self.data.sub_assign(&rhs.data, modulus);
                self
            }
        }

        impl<$s, $t> $cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + DataMut,
            $t: UnsignedInteger,
        {
            /// Performs an element-wise modular addition assignment `self += rhs`.
            #[inline]
            pub fn add_element_wise_assign<M, A>(&mut self, rhs: &$cipher<A>, modulus: M)
            where
                M: FieldContext<$t>,
                A: RawData<Elem = $t> + Data,
            {
                self.data.add_assign(&rhs.data, modulus);
            }

            /// Performs an element-wise modular subtraction assignment `self -= rhs`
            #[inline]
            pub fn sub_element_wise_assign<M, A>(&mut self, rhs: &$cipher<A>, modulus: M)
            where
                M: FieldContext<$t>,
                A: RawData<Elem = $t> + Data,
            {
                self.data.sub_assign(&rhs.data, modulus);
            }
        }

        impl<$s, $t> $cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + Data,
            $t: UnsignedInteger,
        {
            /// Performs in-place element-wise modular addition:`result = self + rhs`,
            #[inline]
            pub fn add_inplace<M, A>(&self, rhs: &Self, result: &mut $cipher<A>, modulus: M)
            where
                M: FieldContext<$t>,
                A: RawData<Elem = $t> + DataMut,
            {
                self.data.add_inplace(&rhs.data, &mut result.data, modulus)
            }

            /// Performs in-place element-wise modular addition:`result = self - rhs`,
            #[inline]
            pub fn sub_inplace<M, A>(&self, rhs: &Self, result: &mut $cipher<A>, modulus: M)
            where
                M: FieldContext<$t>,
                A: RawData<Elem = $t> + DataMut,
            {
                self.data.sub_inplace(&rhs.data, &mut result.data, modulus)
            }
        }
    };
}

macro_rules! impl_basic_operation_multiple_modulus {
    ($cipher:ident < $s:ident, $t:ident >) => {
        impl<$s, $t> $cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + DataOwned,
            $t: UnsignedInteger,
        {
            /// Perform element-wise modular addition `self + rhs`.
            #[inline]
            pub fn add_element_wise<M, A>(
                mut self,
                rhs: &$cipher<A>,
                cipher_single_modulus_len: usize,
                moduli: &[M],
            ) -> Self
            where
                M: FieldContext<$t>,
                A: RawData<Elem = $t> + Data,
            {
                self.add_element_wise_assign(rhs, cipher_single_modulus_len, moduli);
                self
            }

            /// Perform element-wise modular subtraction `self - rhs`.
            #[inline]
            pub fn sub_element_wise<M, A>(
                mut self,
                rhs: &$cipher<A>,
                cipher_single_modulus_len: usize,
                moduli: &[M],
            ) -> Self
            where
                M: FieldContext<$t>,
                A: RawData<Elem = $t> + Data,
            {
                self.sub_element_wise_assign(rhs, cipher_single_modulus_len, moduli);
                self
            }
        }

        impl<$s, $t> $cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + DataMut,
            $t: UnsignedInteger,
        {
            /// Performs an element-wise modular addition assignment `self += rhs`.
            #[inline]
            pub fn add_element_wise_assign<M, A>(
                &mut self,
                rhs: &$cipher<A>,
                cipher_single_modulus_len: usize,
                moduli: &[M],
            ) where
                M: FieldContext<T>,
                A: RawData<Elem = T> + Data,
            {
                izip!(
                    self.data.chunks_exact_mut(cipher_single_modulus_len),
                    rhs.data.chunks_exact(cipher_single_modulus_len),
                    moduli
                )
                .for_each(|(x, y, m)| {
                    ArrayBase(x).add_assign(&ArrayBase(y), *m);
                });
            }

            /// Performs an element-wise modular subtraction assignment `self -= rhs`.
            #[inline]
            pub fn sub_element_wise_assign<M, A>(
                &mut self,
                rhs: &$cipher<A>,
                cipher_single_modulus_len: usize,
                moduli: &[M],
            ) where
                M: FieldContext<T>,
                A: RawData<Elem = T> + Data,
            {
                izip!(
                    self.data.chunks_exact_mut(cipher_single_modulus_len),
                    rhs.data.chunks_exact(cipher_single_modulus_len),
                    moduli
                )
                .for_each(|(x, y, m)| {
                    ArrayBase(x).add_assign(&ArrayBase(y), *m);
                });
            }
        }

        impl<$s, $t> $cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + Data,
            $t: UnsignedInteger,
        {
            /// Performs element-wise modular addition `result = self + rhs`.
            #[inline]
            pub fn add_inplace<M, A, B>(
                &self,
                rhs: &$cipher<A>,
                result: &mut $cipher<B>,
                cipher_single_modulus_len: usize,
                moduli: &[M],
            ) where
                M: FieldContext<T>,
                A: RawData<Elem = T> + Data,
                B: RawData<Elem = T> + DataMut,
            {
                izip!(
                    self.data.chunks_exact(cipher_single_modulus_len),
                    rhs.data.chunks_exact(cipher_single_modulus_len),
                    result.data.chunks_exact_mut(cipher_single_modulus_len),
                    moduli
                )
                .for_each(|(x, y, z, m)| {
                    ArrayBase(x).add_inplace(&ArrayBase(y), &mut ArrayBase(z), *m);
                });
            }

            /// Performs element-wise modular subtraction `result = self - rhs`.
            #[inline]
            pub fn sub_inplace<M, A, B>(
                &self,
                rhs: &$cipher<A>,
                result: &mut $cipher<B>,
                cipher_single_modulus_len: usize,
                moduli: &[M],
            ) where
                M: FieldContext<T>,
                A: RawData<Elem = T> + Data,
                B: RawData<Elem = T> + DataMut,
            {
                izip!(
                    self.data.chunks_exact(cipher_single_modulus_len),
                    rhs.data.chunks_exact(cipher_single_modulus_len),
                    result.data.chunks_exact_mut(cipher_single_modulus_len),
                    moduli
                )
                .for_each(|(x, y, z, m)| {
                    ArrayBase(x).sub_inplace(&ArrayBase(y), &mut ArrayBase(z), *m);
                });
            }
        }
    };
}

macro_rules! impl_ntt {
    ($cipher:ident < $s:ident, $t:ident >,$ntt_cipher:ident) => {
        impl<$s, $t> $cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + DataOwned,
            $t: UnsignedInteger,
        {
            /// ntt transform
            #[inline]
            pub fn into_ntt_form<Table>(mut self, ntt_table: &Table) -> $ntt_cipher<S>
            where
                Table: NttTable<ValueT = $t> + Ntt,
            {
                let poly_length = ntt_table.poly_length();
                self.data.chunks_exact_mut(poly_length).for_each(|poly| {
                    ntt_table.transform_slice(poly);
                });

                $ntt_cipher::new(self.data)
            }
        }

        impl<$s, $t> $cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + Data,
            $t: UnsignedInteger,
        {
            /// ntt transform
            #[inline]
            pub fn to_ntt_form_inplace<Table, A>(
                &self,
                result: &mut $ntt_cipher<A>,
                ntt_table: &Table,
            ) where
                A: RawData<Elem = $t> + DataMut,
                Table: NttTable<ValueT = $t> + Ntt,
            {
                let poly_length = ntt_table.poly_length();

                result.data.copy_from_slice(self.data.as_ref());

                result.data.chunks_exact_mut(poly_length).for_each(|poly| {
                    ntt_table.transform_slice(poly);
                });
            }
        }
    };
}

macro_rules! impl_intt {
    ($ntt_cipher:ident < $s:ident, $t:ident >,$cipher:ident) => {
        impl<$s, $t> $ntt_cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + DataOwned,
            $t: UnsignedInteger,
        {
            /// ntt inverse transform
            #[inline]
            pub fn into_coeff_form<Table>(mut self, ntt_table: &Table) -> $cipher<S>
            where
                Table: NttTable<ValueT = $t> + Ntt,
            {
                let poly_length = ntt_table.poly_length();

                self.data.chunks_exact_mut(poly_length).for_each(|poly| {
                    ntt_table.inverse_transform_slice(poly);
                });

                $cipher::new(self.data)
            }
        }

        impl<$s, $t> $ntt_cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + Data,
            $t: UnsignedInteger,
        {
            /// ntt inverse transform
            #[inline]
            pub fn to_coeff_form_inplace<Table, A>(
                &self,
                result: &mut $cipher<A>,
                ntt_table: &Table,
            ) where
                A: RawData<Elem = $t> + DataMut,
                Table: NttTable<ValueT = $t> + Ntt,
            {
                result.data.copy_from_slice(self.data.as_ref());

                let poly_length = ntt_table.poly_length();

                result
                    .data
                    .chunks_exact_mut(poly_length)
                    .for_each(|values| {
                        ntt_table.inverse_transform_slice(values);
                    });
            }
        }
    };
}

macro_rules! impl_crt_ntt {
    ($cipher:ident < $s:ident, $t:ident >,$ntt_cipher:ident) => {
        impl<$s, $t> $cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + DataOwned,
            $t: UnsignedInteger,
        {
            /// ntt transform
            #[inline]
            pub fn into_ntt_form<Table>(
                self,
                table: &Table,
                cipher_single_modulus_len: usize,
            ) -> $ntt_cipher<$s>
            where
                Table: DcrtTable<ValueT = $t> + Dcrt,
            {
                let poly_length = table.poly_length();

                let Self { mut data } = self;

                data.chunks_exact_mut(cipher_single_modulus_len)
                    .zip(table.iter())
                    .for_each(|(cipher, ntt_table)| {
                        cipher.chunks_exact_mut(poly_length).for_each(|poly| {
                            ntt_table.transform_slice(poly);
                        });
                    });

                $ntt_cipher::new(data)
            }
        }

        impl<$s, $t> $cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + Data,
            $t: UnsignedInteger,
        {
            /// ntt transform
            #[inline]
            pub fn to_ntt_form_inplace<Table, A>(
                &self,
                result: &mut $ntt_cipher<A>,
                table: &Table,
                cipher_single_modulus_len: usize,
            ) where
                Table: DcrtTable<ValueT = $t> + Dcrt,
                A: RawData<Elem = $t> + DataMut,
            {
                result.data.copy_from_slice(self.data.as_ref());

                let poly_length = table.poly_length();

                result
                    .data
                    .chunks_exact_mut(cipher_single_modulus_len)
                    .zip(table.iter())
                    .for_each(|(cipher, ntt_table)| {
                        cipher.chunks_exact_mut(poly_length).for_each(|a| {
                            ntt_table.transform_slice(a);
                        });
                    });
            }
        }
    };
}

macro_rules! impl_crt_intt {
    ($ntt_cipher:ident < $s:ident, $t:ident >,$cipher:ident) => {
        impl<$s, $t> $ntt_cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + DataOwned,
            $t: UnsignedInteger,
        {
            /// ntt inverse transform
            #[inline]
            pub fn into_coeff_form<Table>(
                self,
                table: &Table,
                cipher_single_modulus_len: usize,
            ) -> $cipher<$s>
            where
                Table: DcrtTable<ValueT = $t> + Dcrt,
            {
                let poly_length = table.poly_length();

                let Self { mut data } = self;

                data.chunks_exact_mut(cipher_single_modulus_len)
                    .zip(table.iter())
                    .for_each(|(cipher, ntt_table)| {
                        cipher.chunks_exact_mut(poly_length).for_each(|poly| {
                            ntt_table.inverse_transform_slice(poly);
                        });
                    });

                $cipher::new(data)
            }
        }

        impl<$s, $t> $ntt_cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + Data,
            $t: UnsignedInteger,
        {
            /// ntt inverse transform
            #[inline]
            pub fn to_coeff_form_inplace<Table, A>(
                &self,
                result: &mut $cipher<A>,
                table: &Table,
                cipher_single_modulus_len: usize,
            ) where
                Table: DcrtTable<ValueT = $t> + Dcrt,
                A: RawData<Elem = $t> + DataMut,
            {
                result.data.copy_from_slice(self.data.as_ref());

                let poly_length = table.poly_length();

                result
                    .data
                    .chunks_exact_mut(cipher_single_modulus_len)
                    .zip(table.iter())
                    .for_each(|(cipher, ntt_table)| {
                        cipher.chunks_exact_mut(poly_length).for_each(|a| {
                            ntt_table.inverse_transform_slice(a);
                        });
                    });
            }
        }
    };
}

macro_rules! impl_zero {
    ($cipher:ident < $s:ident, $t:ident >) => {
        impl<$s, $t> $cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + DataOwned,
            $t: UnsignedInteger,
        {
            /// Creates with all entries equal to zero.
            #[inline]
            pub fn zero(cipher_len: usize) -> Self {
                Self {
                    data: ArrayBase::from_vec(vec![T::ZERO; cipher_len]),
                }
            }
        }

        impl<$s, $t> $cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + DataMut,
            $t: UnsignedInteger,
        {
            /// Set all entries equal to zero.
            #[inline]
            pub fn set_zero(&mut self) {
                self.data.set_zero();
            }
        }
    };
}
