macro_rules! impl_common {
    ($cipher:ident < $s:ident, $t:ident >) => {
        impl<$s, $t> $cipher<$s, $t>
        where
            $s: RawData<Elem = $t>,
            $t: UnsignedInteger,
        {
            #[doc = concat!(r" Creates a new [`",stringify!($cipher),"<",stringify!($s),", ",stringify!($t),">`].")]
            #[inline]
            pub fn new(data: ArrayBase<S>) -> Self {
                Self { data }
            }
        }

        impl<$s, $t> AsRef<[$t]> for $cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + Data,
            $t: UnsignedInteger,
        {
            #[inline]
            fn as_ref(&self)->&[$t]{
                self.data.as_ref()
            }

        }

        impl<$s, $t> AsMut<[$t]> for $cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + DataMut,
            $t: UnsignedInteger,
        {
            #[inline]
            fn as_mut(&mut self)->&mut [$t]{
                self.data.as_mut()
            }

        }
    };
}

macro_rules! impl_bytes_conversion {
    ($cipher:ident < $s:ident, $t:ident >) => {
        impl<$s, $t> $cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + DataOwned,
            $t: UnsignedInteger,
        {
            #[doc = concat!(r" Creates a new [`",stringify!($cipher),"<",stringify!($s),", ",stringify!($t),">`] from bytes `data`.")]
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

macro_rules! impl_zero {
    ($cipher:ident < $s:ident, $t:ident >) => {
        impl<$s, $t> $cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + DataOwned,
            $t: UnsignedInteger,
        {
            #[doc = concat!(r" Creates a new [`",stringify!($cipher),"<",stringify!($s),", ",stringify!($t),">`] with all values or coefficients equal to zero.")]
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
            /// Set all values or coefficients equal to zero.
            #[inline]
            pub fn set_zero(&mut self) {
                self.data.set_zero();
            }
        }
    };
}

macro_rules! impl_basic_operation_single_modulus {
    ($cipher:ident < $s:ident, $t:ident >) => {
        impl<$s, $t> $cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + DataMut,
            $t: UnsignedInteger,
        {
            /// Perform element-wise modular addition `self + rhs`.
            #[inline]
            pub fn add_element_wise<M, A>(mut self, rhs: &$cipher<A>, modulus: M) -> Self
            where
                M: FieldContext<$t>,
                A: RawData<Elem = $t> + Data,
            {
                self.data.add_element_wise_assign(&rhs.data, modulus);
                self
            }

            /// Perform element-wise modular subtraction `self - rhs`.
            #[inline]
            pub fn sub_element_wise<M, A>(mut self, rhs: &$cipher<A>, modulus: M) -> Self
            where
                M: FieldContext<$t>,
                A: RawData<Elem = $t> + Data,
            {
                self.data.sub_element_wise_assign(&rhs.data, modulus);
                self
            }

            /// Performs an element-wise modular addition assignment `self += rhs`.
            #[inline]
            pub fn add_element_wise_assign<M, A>(&mut self, rhs: &$cipher<A>, modulus: M)
            where
                M: FieldContext<$t>,
                A: RawData<Elem = $t> + Data,
            {
                self.data.add_element_wise_assign(&rhs.data, modulus);
            }

            /// Performs an element-wise modular subtraction assignment `self -= rhs`
            #[inline]
            pub fn sub_element_wise_assign<M, A>(&mut self, rhs: &$cipher<A>, modulus: M)
            where
                M: FieldContext<$t>,
                A: RawData<Elem = $t> + Data,
            {
                self.data.sub_element_wise_assign(&rhs.data, modulus);
            }
        }

        impl<$s, $t> $cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + Data,
            $t: UnsignedInteger,
        {
            /// Performs in-place element-wise modular addition:`result = self + rhs`,
            #[inline]
            pub fn add_element_wise_inplace<M, A, B>(
                &self,
                rhs: &$cipher<A>,
                result: &mut $cipher<B>,
                modulus: M,
            ) where
                M: FieldContext<$t>,
                A: RawData<Elem = $t> + Data,
                B: RawData<Elem = $t> + DataMut,
            {
                self.data
                    .add_element_wise_inplace(&rhs.data, &mut result.data, modulus)
            }

            /// Performs in-place element-wise modular addition:`result = self - rhs`,
            #[inline]
            pub fn sub_element_wise_inplace<M, A, B>(
                &self,
                rhs: &$cipher<A>,
                result: &mut $cipher<B>,
                modulus: M,
            ) where
                M: FieldContext<$t>,
                A: RawData<Elem = $t> + Data,
                B: RawData<Elem = $t> + DataMut,
            {
                self.data
                    .sub_element_wise_inplace(&rhs.data, &mut result.data, modulus)
            }
        }
    };
}

macro_rules! impl_basic_operation_multiple_modulus {
    ($cipher:ident < $s:ident, $t:ident >) => {
        impl<$s, $t> $cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + DataMut,
            $t: UnsignedInteger,
        {
            // #[inline]
            // pub fn iter_crt_poly_mut(
            //     &mut self,
            //     crt_poly_length: usize,
            // ) -> std::slice::ChunksExactMut<'_, T> {
            //     self.data.chunks_exact_mut(crt_poly_length)
            // }

            /// Perform element-wise modular addition `self + rhs`.
            #[inline]
            pub fn add_element_wise<M, A>(
                mut self,
                rhs: &$cipher<A>,
                poly_length: usize,
                crt_poly_length: usize,
                moduli: &[M],
            ) -> Self
            where
                M: FieldContext<$t>,
                A: RawData<Elem = $t> + Data,
            {
                self.add_element_wise_assign(rhs, poly_length, crt_poly_length, moduli);
                self
            }

            /// Perform element-wise modular subtraction `self - rhs`.
            #[inline]
            pub fn sub_element_wise<M, A>(
                mut self,
                rhs: &$cipher<A>,
                poly_length: usize,
                crt_poly_length: usize,
                moduli: &[M],
            ) -> Self
            where
                M: FieldContext<$t>,
                A: RawData<Elem = $t> + Data,
            {
                self.sub_element_wise_assign(rhs, poly_length, crt_poly_length, moduli);
                self
            }

            /// Performs an element-wise modular addition assignment `self += rhs`.
            #[inline]
            pub fn add_element_wise_assign<M, A>(
                &mut self,
                rhs: &$cipher<A>,
                poly_length: usize,
                crt_poly_length: usize,
                moduli: &[M],
            ) where
                M: FieldContext<T>,
                A: RawData<Elem = T> + Data,
            {
                izip!(
                    self.data.chunks_exact_mut(crt_poly_length),
                    rhs.data.chunks_exact(crt_poly_length),
                )
                .for_each(|(x, y)| {
                    izip!(
                        x.chunks_exact_mut(poly_length),
                        y.chunks_exact(poly_length),
                        moduli
                    )
                    .for_each(|(a, b, &modulus)| {
                        ArrayBase(a).add_element_wise_assign(&ArrayBase(b), modulus);
                    });
                });
            }

            /// Performs an element-wise modular subtraction assignment `self -= rhs`.
            #[inline]
            pub fn sub_element_wise_assign<M, A>(
                &mut self,
                rhs: &$cipher<A>,
                poly_length: usize,
                crt_poly_length: usize,
                moduli: &[M],
            ) where
                M: FieldContext<T>,
                A: RawData<Elem = T> + Data,
            {
                izip!(
                    self.data.chunks_exact_mut(crt_poly_length),
                    rhs.data.chunks_exact(crt_poly_length),
                )
                .for_each(|(x, y)| {
                    izip!(
                        x.chunks_exact_mut(poly_length),
                        y.chunks_exact(poly_length),
                        moduli
                    )
                    .for_each(|(a, b, &modulus)| {
                        ArrayBase(a).sub_element_wise_assign(&ArrayBase(b), modulus);
                    });
                });
            }
        }

        impl<$s, $t> $cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + Data,
            $t: UnsignedInteger,
        {
            // #[inline]
            // pub fn iter_each_modulus(
            //     &self,
            //     single_modulus_len: usize,
            // ) -> std::slice::ChunksExact<'_, T> {
            //     self.data.chunks_exact(single_modulus_len)
            // }

            /// Performs element-wise modular addition `result = self + rhs`.
            #[inline]
            pub fn add_element_wise_inplace<M, A, B>(
                &self,
                rhs: &$cipher<A>,
                result: &mut $cipher<B>,
                poly_length: usize,
                crt_poly_length: usize,
                moduli: &[M],
            ) where
                M: FieldContext<T>,
                A: RawData<Elem = T> + Data,
                B: RawData<Elem = T> + DataMut,
            {
                izip!(
                    self.data.chunks_exact(crt_poly_length),
                    rhs.data.chunks_exact(crt_poly_length),
                    result.data.chunks_exact_mut(crt_poly_length),
                )
                .for_each(|(x, y, z)| {
                    izip!(
                        x.chunks_exact(poly_length),
                        y.chunks_exact(poly_length),
                        z.chunks_exact_mut(poly_length),
                        moduli
                    )
                    .for_each(|(a, b, c, &modulus)| {
                        ArrayBase(a).add_element_wise_inplace(
                            &ArrayBase(b),
                            &mut ArrayBase(c),
                            modulus,
                        );
                    });
                });
            }

            /// Performs element-wise modular subtraction `result = self - rhs`.
            #[inline]
            pub fn sub_element_wise_inplace<M, A, B>(
                &self,
                rhs: &$cipher<A>,
                result: &mut $cipher<B>,
                poly_length: usize,
                crt_poly_length: usize,
                moduli: &[M],
            ) where
                M: FieldContext<T>,
                A: RawData<Elem = T> + Data,
                B: RawData<Elem = T> + DataMut,
            {
                izip!(
                    self.data.chunks_exact(crt_poly_length),
                    rhs.data.chunks_exact(crt_poly_length),
                    result.data.chunks_exact_mut(crt_poly_length),
                )
                .for_each(|(x, y, z)| {
                    izip!(
                        x.chunks_exact(poly_length),
                        y.chunks_exact(poly_length),
                        z.chunks_exact_mut(poly_length),
                        moduli
                    )
                    .for_each(|(a, b, c, &modulus)| {
                        ArrayBase(a).sub_element_wise_inplace(
                            &ArrayBase(b),
                            &mut ArrayBase(c),
                            modulus,
                        );
                    });
                });
            }
        }
    };
}

macro_rules! impl_ntt {
    ($cipher:ident < $s:ident, $t:ident >,$ntt_cipher:ident) => {
        impl<$s, $t> $cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + DataMut,
            $t: UnsignedInteger,
        {
            /// Transforms `self` to ntt form.
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
            /// Transforms `self` to ntt form and stores in `result`.
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
            $s: RawData<Elem = $t> + DataMut,
            $t: UnsignedInteger,
        {
            /// Transforms `self` to coefficient form.
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
            /// Transforms `self` to coefficient form and stores in `result`.
            #[inline]
            pub fn to_coeff_form_inplace<Table, A>(
                &self,
                result: &mut $cipher<A>,
                ntt_table: &Table,
            ) where
                A: RawData<Elem = $t> + DataMut,
                Table: NttTable<ValueT = $t> + Ntt,
            {
                let poly_length = ntt_table.poly_length();
                result.data.copy_from_slice(self.data.as_ref());
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
            $s: RawData<Elem = $t> + DataMut,
            $t: UnsignedInteger,
        {
            /// Transforms `self` to ntt form.
            #[inline]
            pub fn into_ntt_form<Table>(self, table: &Table) -> $ntt_cipher<$s>
            where
                Table: DcrtTable<ValueT = $t> + Dcrt,
            {
                let crt_poly_length = table.crt_poly_length();
                let Self { mut data } = self;
                data.chunks_exact_mut(crt_poly_length).for_each(|crt_poly| {
                    table.transform_slice(crt_poly);
                });
                $ntt_cipher::new(data)
            }
        }

        impl<$s, $t> $cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + Data,
            $t: UnsignedInteger,
        {
            /// Transforms `self` to ntt form and stores in `result`.
            #[inline]
            pub fn to_ntt_form_inplace<Table, A>(&self, result: &mut $ntt_cipher<A>, table: &Table)
            where
                Table: DcrtTable<ValueT = $t> + Dcrt,
                A: RawData<Elem = $t> + DataMut,
            {
                let crt_poly_length = table.crt_poly_length();
                result.data.copy_from_slice(self.data.as_ref());
                result
                    .data
                    .chunks_exact_mut(crt_poly_length)
                    .for_each(|crt_poly| {
                        table.transform_slice(crt_poly);
                    });
            }
        }
    };
}

macro_rules! impl_crt_intt {
    ($ntt_cipher:ident < $s:ident, $t:ident >,$cipher:ident) => {
        impl<$s, $t> $ntt_cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + DataMut,
            $t: UnsignedInteger,
        {
            /// Transforms `self` to coefficient form.
            #[inline]
            pub fn into_coeff_form<Table>(self, table: &Table) -> $cipher<$s>
            where
                Table: DcrtTable<ValueT = $t> + Dcrt,
            {
                let crt_poly_length = table.crt_poly_length();
                let Self { mut data } = self;
                data.chunks_exact_mut(crt_poly_length).for_each(|crt_poly| {
                    table.inverse_transform_slice(crt_poly);
                });
                $cipher::new(data)
            }
        }

        impl<$s, $t> $ntt_cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + Data,
            $t: UnsignedInteger,
        {
            /// Transforms `self` to coefficient form and stores in `result`.
            #[inline]
            pub fn to_coeff_form_inplace<Table, A>(&self, result: &mut $cipher<A>, table: &Table)
            where
                Table: DcrtTable<ValueT = $t> + Dcrt,
                A: RawData<Elem = $t> + DataMut,
            {
                let crt_poly_length = table.crt_poly_length();
                result.data.copy_from_slice(self.data.as_ref());
                result
                    .data
                    .chunks_exact_mut(crt_poly_length)
                    .for_each(|crt_poly| {
                        table.inverse_transform_slice(crt_poly);
                    });
            }
        }
    };
}
