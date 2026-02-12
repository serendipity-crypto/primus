use primus_factor::FactorMul;
use primus_integer::AsInto;
use primus_integer::{UnsignedInteger, izip};
use primus_modulo::Modulo;
use primus_modulo::MulModulo;
use primus_modulus::UintModulus;
use primus_reduce::FieldContext;

use crate::RNSBase;

#[derive(Clone)]
pub struct BaseConverter<T: UnsignedInteger, M: FieldContext<T>> {
    /// The base convert from.
    ibase: RNSBase<T, M>,
    /// The base convert into.
    obase: RNSBase<T, M>,
    /// Precomputed helper matrix, speed up conversion.
    base_change_matrix: Vec<T>,
}

impl<T: UnsignedInteger, M: FieldContext<T>> BaseConverter<T, M> {
    pub fn new(ibase: &RNSBase<T, M>, obase: &RNSBase<T, M>) -> Self {
        let ibase_moduli_count = ibase.moduli_count();
        let obase_moduli_count = obase.moduli_count();

        assert!(
            ibase_moduli_count.checked_mul(obase_moduli_count).is_some(),
            "the len can not be too large!"
        );

        let mut base_change_matrix = vec![T::ZERO; ibase_moduli_count * obase_moduli_count];

        for (row, &modulus) in base_change_matrix
            .chunks_exact_mut(ibase_moduli_count)
            .zip(obase.moduli())
        {
            for (ele, m_i) in row.iter_mut().zip(ibase.iter_punctured_product()) {
                *ele = m_i.modulo(modulus);
            }
        }

        Self {
            ibase: ibase.clone(),
            obase: obase.clone(),
            base_change_matrix,
        }
    }

    pub fn ibase(&self) -> &RNSBase<T, M> {
        &self.ibase
    }

    pub fn obase(&self) -> &RNSBase<T, M> {
        &self.obase
    }

    pub fn ibase_moduli_count(&self) -> usize {
        self.ibase.moduli_count()
    }

    pub fn obase_moduli_count(&self) -> usize {
        self.obase.moduli_count()
    }

    fn iter_base_change_matrix(&self) -> std::slice::ChunksExact<'_, T> {
        self.base_change_matrix
            .chunks_exact(self.ibase_moduli_count())
    }

    /// Convert residue numbers between different basis, output the result into `residues_out`.
    pub fn fast_convert(&self, residues_in: &[T], residues_out: &mut [T]) {
        debug_assert_eq!(residues_in.len(), self.ibase_moduli_count());
        debug_assert_eq!(residues_out.len(), self.obase_moduli_count());

        let temp: Vec<T> = izip!(
            residues_in,
            self.ibase.inv_punctured_product_mod_modulus(),
            self.ibase.moduli()
        )
        .map(|(&value, &inv, modulus)| inv.factor_mul_modulo(value, modulus.value_unchecked()))
        .collect();

        izip!(
            residues_out,
            self.iter_base_change_matrix(),
            self.obase.moduli()
        )
        .for_each(|(ele, base_chang_row, modulus)| {
            *ele = modulus.reduce_dot_product(&temp, base_chang_row);
        });
    }

    pub fn fast_convert_array(
        &self,
        crt_poly_in: &[T],
        crt_poly_out: &mut [T],
        poly_length: usize,
    ) {
        let ibase_moduli_count = self.ibase_moduli_count();

        let mut temp: Vec<T> = vec![T::ZERO; ibase_moduli_count * poly_length];

        izip!(
            crt_poly_in.chunks_exact(poly_length),
            self.ibase.inv_punctured_product_mod_modulus(),
            self.ibase.moduli()
        )
        .enumerate()
        .for_each(
            |(i, (poly, &inv_punctured_product_mod_modulus, &modulus))| {
                if inv_punctured_product_mod_modulus.value().is_one() {
                    izip!(poly, temp.iter_mut().skip(i).step_by(ibase_moduli_count)).for_each(
                        |(&x, ele)| {
                            *ele = x.modulo(modulus);
                        },
                    );
                } else {
                    izip!(poly, temp.iter_mut().skip(i).step_by(ibase_moduli_count)).for_each(
                        |(&x, ele)| {
                            *ele = x.mul_modulo(
                                inv_punctured_product_mod_modulus,
                                UintModulus(modulus.value_unchecked()),
                            );
                        },
                    );
                }
            },
        );

        izip!(
            crt_poly_out.chunks_exact_mut(poly_length),
            self.iter_base_change_matrix(),
            self.obase.moduli()
        )
        .for_each(|(poly, inv_punctured_product_mod_modulus, modulus)| {
            izip!(poly, temp.chunks_exact(ibase_moduli_count)).for_each(|(ele, product)| {
                *ele = modulus.reduce_dot_product(product, inv_punctured_product_mod_modulus);
            });
        });
    }

    pub fn exact_convert_array(
        &self,
        crt_poly_in: &[T],
        crt_poly_out: &mut [T],
        poly_length: usize,
    ) {
        let ibase_moduli_count = self.ibase_moduli_count();

        assert_eq!(
            self.obase_moduli_count(),
            1,
            "out base in exact_convert_array must be one."
        );

        let mut temp: Vec<T> = vec![T::ZERO; ibase_moduli_count * poly_length];
        let mut v: Vec<f64> = vec![0.0f64; ibase_moduli_count * poly_length];
        let mut aggregated_rounded_v: Vec<T> = vec![T::ZERO; poly_length];

        // Calculate [x_{i} * \hat{q_{i}}]_{q_{i}}
        izip!(
            crt_poly_in.chunks_exact(poly_length),
            self.ibase.inv_punctured_product_mod_modulus(),
            self.ibase.moduli()
        )
        .enumerate()
        .for_each(
            |(i, (poly, &inv_punctured_product_mod_modulus, &modulus))| {
                let divisor: f64 = modulus.value_unchecked().as_into();
                if inv_punctured_product_mod_modulus.value().is_one() {
                    // No multiplication needed
                    izip!(
                        poly,
                        temp.iter_mut().skip(i).step_by(ibase_moduli_count),
                        v.iter_mut().skip(i).step_by(ibase_moduli_count)
                    )
                    .for_each(|(&x, ele, fele)| {
                        // Reduce modulo ibase element
                        *ele = x.modulo(modulus);
                        let dividend: f64 = (*ele).as_into();
                        *fele = dividend / divisor;
                    });
                } else {
                    // Multiplication needed
                    izip!(
                        poly,
                        temp.iter_mut().skip(i).step_by(ibase_moduli_count),
                        v.iter_mut().skip(i).step_by(ibase_moduli_count)
                    )
                    .for_each(|(&x, ele, fele)| {
                        // Multiply coefficient of in with ibase_.inv_punctured_prod_mod_base_array_ element
                        *ele = x.mul_modulo(
                            inv_punctured_product_mod_modulus,
                            UintModulus(modulus.value_unchecked()),
                        );
                        let dividend: f64 = (*ele).as_into();
                        *fele = dividend / divisor;
                    });
                }
            },
        );

        // Aggrate v and rounding
        izip!(
            v.chunks_exact(ibase_moduli_count),
            aggregated_rounded_v.iter_mut()
        )
        .for_each(|(vi, ri)| {
            // Otherwise a memory space of the last execution will be used.
            let aggregated_v: f64 = vi.iter().sum();
            *ri = (aggregated_v + 0.5).as_into();
        });

        let p = self.obase.moduli()[0];
        let q_mod_p = self.ibase.moduli_product().0.modulo(p);
        let base_change_matrix_first = self.iter_base_change_matrix().next().unwrap();

        // Final multiplication
        izip!(
            crt_poly_out,
            temp.chunks_exact(ibase_moduli_count),
            aggregated_rounded_v,
        )
        .for_each(|(coeff, b, v)| {
            // Compute the base conversion sum modulo obase element
            let sum_mod_obase = p.reduce_dot_product(b, base_change_matrix_first);
            // Minus v*[q]_{p} mod p
            let v_q_mod_p = p.reduce_mul(v, q_mod_p);
            *coeff = p.reduce_sub(sum_mod_obase, v_q_mod_p);
        });
    }
}

#[cfg(test)]
mod tests {
    use primus_integer::BigUint;
    use primus_modulus::BarrettModulus;
    use rand::RngExt;

    use super::*;

    type ValueT = u32;

    #[test]
    #[ignore = "just for print"]
    fn test_base_convert() {
        let mut rng = rand::rng();

        let i = [31, 37, 41, 43];
        let o = [47, 53, 59, 61];
        let d = [i, o].concat();

        let basis: Vec<_> = i.map(<BarrettModulus<ValueT>>::new).into_iter().collect();
        let ibasis = RNSBase::new(&basis).unwrap();
        let in_len = ibasis.moduli_count();
        println!("{:?}", ibasis.moduli_product());

        let basis: Vec<_> = o.map(<BarrettModulus<ValueT>>::new).into_iter().collect();
        let obasis = RNSBase::new(&basis).unwrap();
        let out_len = obasis.moduli_count();
        println!("{:?}", obasis.moduli_product());

        let converter = BaseConverter::new(&ibasis, &obasis);

        let basis: Vec<_> = d
            .iter()
            .copied()
            .map(<BarrettModulus<ValueT>>::new)
            .collect();
        let dbasis = RNSBase::new(&basis).unwrap();
        println!("{:?}\n", dbasis.moduli_product());
        let mut ibasis_product = ibasis.moduli_product().0.to_vec();
        ibasis_product.resize(dbasis.moduli_product().len(), 0);

        let ibasis_product = BigUint(ibasis_product);

        for _ in 0..10 {
            let mut residues_in = Vec::with_capacity(in_len);
            for i in 0..in_len {
                residues_in.push(rng.random_range(0..converter.ibase.moduli()[i].value()));
            }

            let value = ibasis.compose(&residues_in);
            println!("{:?}", value);

            let mut residues_out = vec![0; out_len];
            converter.fast_convert(&residues_in, &mut residues_out);

            residues_in.extend_from_slice(&residues_out);
            let mut value = dbasis.compose(&residues_in);
            while value.cmp(&ibasis_product).is_ge() {
                let _ = value.sub_assign(&ibasis_product);
            }
            println!("{:?}\n", value);
        }
    }
}
