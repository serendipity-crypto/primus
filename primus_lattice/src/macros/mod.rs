macro_rules! impl_common {
    ($cipher:ident < $s:ident, $t:ident >) => {
        impl<$s, $t> $cipher<$s, $t>
        where
            $s: RawData<Elem = $t>,
            $t: UnsignedInteger,
        {
            #[doc = concat!(r" Creates a new [`",stringify!($cipher),"<",stringify!($s),", ",stringify!($t),">`].")]
            #[inline(always)]
            pub fn new(data: $s) -> Self {
                Self(data)
            }
        }

        impl<$s, $t> AsRef<[$t]> for $cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + Data,
            $t: UnsignedInteger,
        {
            #[inline(always)]
            fn as_ref(&self)->&[$t]{
                self.0.as_ref()
            }

        }

        impl<$s, $t> AsMut<[$t]> for $cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + DataMut,
            $t: UnsignedInteger,
        {
            #[inline(always)]
            fn as_mut(&mut self)->&mut [$t]{
                self.0.as_mut()
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

                Self(<$s>::from_slice(converted_data))
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

                self.0.copy_from_slice(converted_data);
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
                let converted_data: &[u8] = bytemuck::cast_slice(self.as_ref());

                converted_data.to_vec()
            }

            /// Converts `self` into bytes, stored in `data`.
            #[inline]
            pub fn to_bytes_inplace(&self, data: &mut [u8]) {
                let converted_data: &[u8] = bytemuck::cast_slice(self.as_ref());

                data.copy_from_slice(converted_data);
            }

            /// Returns the bytes count.
            #[inline]
            pub fn byte_count(&self) -> usize {
                self.0.byte_count()
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
            paste::paste! {
                #[doc = concat!(r" Creates a new [`",stringify!($cipher),"<",stringify!($s),", ",stringify!($t),">`] with all values or coefficients equal to zero.")]
                #[inline]
                pub fn zero([<$cipher:snake _len>]: usize) -> Self {
                    Self(<$s>::zero([<$cipher:snake _len>]))
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
                self.0.set_zero();
            }
        }
    };
}
macro_rules! impl_iters {
    ($cipher:ident) => {
        paste::paste! {
            pub struct [<$cipher Iter>]<'a, T>
            where
                T: UnsignedInteger,
            {
                pub(crate) iter: core::slice::ChunksExact<'a, T>
            }

            impl<'a, T: UnsignedInteger> [<$cipher Iter>]<'a, T> {
                #[inline]
                pub fn new(data:&'a [T], [<$cipher:snake _len>]:usize) -> Self{
                    Self {
                        iter: data.chunks_exact([<$cipher:snake _len>])
                    }
                }
            }

            impl<'a, T: UnsignedInteger> Iterator for [<$cipher Iter>]<'a, T> {
                type Item = $cipher<&'a [T], T>;

                #[inline]
                fn next(&mut self) -> Option<Self::Item> {
                    self.iter.next().map(|slice| $cipher(slice))
                }
            }

            impl<'a, T: UnsignedInteger> core::iter::FusedIterator for [<$cipher Iter>]<'a, T> {}
        }

        paste::paste! {
            pub struct [<$cipher IterMut>]<'a, T>
            where
                T: UnsignedInteger,
            {
                pub(crate) iter: core::slice::ChunksExactMut<'a, T>
            }

            impl<'a, T: UnsignedInteger> [<$cipher IterMut>]<'a, T> {
                #[inline]
                pub fn new(data:&'a mut [T], [<$cipher:snake _len>]:usize) -> Self{
                    Self {
                        iter: data.chunks_exact_mut([<$cipher:snake _len>])
                    }
                }
            }

            impl<'a, T: UnsignedInteger> Iterator for [<$cipher IterMut>]<'a, T> {
                type Item = $cipher<&'a mut [T], T>;

                #[inline]
                fn next(&mut self) -> Option<Self::Item> {
                    self.iter.next().map(|slice| $cipher(slice))
                }
            }

            impl<'a, T: UnsignedInteger> core::iter::FusedIterator for [<$cipher IterMut>]<'a, T> {}
        }
    };
}

macro_rules! impl_iter_sub_structure {
    ($cipher:ident < $s:ident, $t:ident >, $sub:ident) => {
        impl<$s, $t> $cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + Data,
            $t: UnsignedInteger,
        {
            paste::paste! {
                #[inline]
                pub fn [<iter_ $sub:snake>]<'a>(&'a self, [<$sub:snake _len>]: usize) -> [<$sub Iter>]<'a, $t> {
                    [<$sub Iter>] {
                        iter: self.0.chunks_exact([<$sub:snake _len>])
                    }

                }
            }
        }

        impl<$s, $t> $cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + DataMut,
            $t: UnsignedInteger,
        {
            paste::paste! {
                #[inline]
                pub fn [<iter_ $sub:snake _mut>]<'a>(
                    &'a mut self,
                    [<$sub:snake _len>]: usize,
                ) -> [<$sub IterMut>]<'a, $t> {
                    [<$sub IterMut>] {
                        iter: self.0.chunks_exact_mut([<$sub:snake _len>])
                    }
                }
            }
        }
    };
    ($cipher:ident < $s:ident, $t:ident >, $sub:ident, $sub_short:ident) => {
        impl<$s, $t> $cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + Data,
            $t: UnsignedInteger,
        {
            paste::paste! {
                #[inline]
                pub fn [<iter_ $sub_short>]<'a>(&'a self, [<$sub_short _len>]: usize) -> [<$sub Iter>]<'a, $t> {
                    [<$sub Iter>] {
                        iter: self.0.chunks_exact([<$sub_short _len>])
                    }

                }
            }
        }

        impl<$s, $t> $cipher<$s, $t>
        where
            $s: RawData<Elem = $t> + DataMut,
            $t: UnsignedInteger,
        {
            paste::paste! {
                #[inline]
                pub fn [<iter_ $sub_short _mut>]<'a>(
                    &'a mut self,
                    [<$sub_short _len>]: usize,
                ) -> [<$sub IterMut>]<'a, $t> {
                    [<$sub IterMut>] {
                        iter: self.0.chunks_exact_mut([<$sub_short _len>])
                    }
                }
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
                ArrayBase(self.as_mut()).add_element_wise_assign(&ArrayBase(rhs.as_ref()), modulus);
                self
            }

            /// Perform element-wise modular subtraction `self - rhs`.
            #[inline]
            pub fn sub_element_wise<M, A>(mut self, rhs: &$cipher<A>, modulus: M) -> Self
            where
                M: FieldContext<$t>,
                A: RawData<Elem = $t> + Data,
            {
                ArrayBase(self.as_mut()).sub_element_wise_assign(&ArrayBase(rhs.as_ref()), modulus);
                self
            }

            /// Performs an element-wise modular addition assignment `self += rhs`.
            #[inline]
            pub fn add_element_wise_assign<M, A>(&mut self, rhs: &$cipher<A>, modulus: M)
            where
                M: FieldContext<$t>,
                A: RawData<Elem = $t> + Data,
            {
                ArrayBase(self.as_mut()).add_element_wise_assign(&ArrayBase(rhs.as_ref()), modulus);
            }

            /// Performs an element-wise modular subtraction assignment `self -= rhs`
            #[inline]
            pub fn sub_element_wise_assign<M, A>(&mut self, rhs: &$cipher<A>, modulus: M)
            where
                M: FieldContext<$t>,
                A: RawData<Elem = $t> + Data,
            {
                ArrayBase(self.as_mut()).sub_element_wise_assign(&ArrayBase(rhs.as_ref()), modulus);
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
                ArrayBase(self.as_ref()).add_element_wise_inplace(
                    &ArrayBase(rhs.as_ref()),
                    &mut ArrayBase(result.as_mut()),
                    modulus,
                )
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
                ArrayBase(self.as_ref()).sub_element_wise_inplace(
                    &ArrayBase(rhs.as_ref()),
                    &mut ArrayBase(result.as_mut()),
                    modulus,
                )
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
                    self.0.chunks_exact_mut(crt_poly_length),
                    rhs.0.chunks_exact(crt_poly_length),
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
                    self.0.chunks_exact_mut(crt_poly_length),
                    rhs.0.chunks_exact(crt_poly_length),
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
                    self.0.chunks_exact(crt_poly_length),
                    rhs.0.chunks_exact(crt_poly_length),
                    result.0.chunks_exact_mut(crt_poly_length),
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
                    self.0.chunks_exact(crt_poly_length),
                    rhs.0.chunks_exact(crt_poly_length),
                    result.0.chunks_exact_mut(crt_poly_length),
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
                self.0.chunks_exact_mut(poly_length).for_each(|poly| {
                    ntt_table.transform_slice(poly);
                });
                $ntt_cipher::new(self.0)
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
                result.0.copy_from_slice(self.as_ref());
                result.0.chunks_exact_mut(poly_length).for_each(|poly| {
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
                self.0.chunks_exact_mut(poly_length).for_each(|poly| {
                    ntt_table.inverse_transform_slice(poly);
                });
                $cipher::new(self.0)
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
                result.0.copy_from_slice(self.as_ref());
                result.0.chunks_exact_mut(poly_length).for_each(|values| {
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
                let Self(mut data) = self;
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
                result.0.copy_from_slice(self.as_ref());
                result
                    .0
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
                let Self(mut data) = self;
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
                result.0.copy_from_slice(self.as_ref());
                result
                    .0
                    .chunks_exact_mut(crt_poly_length)
                    .for_each(|crt_poly| {
                        table.inverse_transform_slice(crt_poly);
                    });
            }
        }
    };
}
