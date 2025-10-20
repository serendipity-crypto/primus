use primus_factor::FactorMul;
use primus_integer::{UnsignedInteger, izip};
use primus_reduce::FieldContext;

use crate::RNSBase;

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
        let input_len = ibase.moduli_count();
        let output_len = obase.moduli_count();

        assert!(
            input_len.checked_mul(output_len).is_some(),
            "the len can not be too large!"
        );

        let mut base_change_matrix = vec![T::ZERO; input_len * output_len];

        for (row, modulus) in base_change_matrix
            .chunks_exact_mut(input_len)
            .zip(obase.moduli())
        {
            for (ele, mi) in row.iter_mut().zip(ibase.punctured_product_iter()) {
                *ele = modulus.reduce(mi);
            }
        }

        Self {
            ibase: ibase.clone(),
            obase: obase.clone(),
            base_change_matrix,
        }
    }

    fn iter_base_change_matrix(&self) -> std::slice::ChunksExact<'_, T> {
        self.base_change_matrix
            .chunks_exact(self.ibase.moduli_count())
    }

    /// Convert residue numbers between different basis, output the result into `values_out`.
    pub fn fast_convert(&self, residues_in: &[T], residues_out: &mut [T]) {
        debug_assert_eq!(residues_in.len(), self.ibase.moduli_count());
        debug_assert_eq!(residues_out.len(), self.obase.moduli_count());

        let mut temp: Vec<T> = Vec::with_capacity(residues_in.len());
        for (&value, &inv, modulus) in izip!(
            residues_in,
            self.ibase.inv_punctured_product_mod_modulus(),
            self.ibase.moduli()
        ) {
            temp.push(inv.factor_mul_modulo(value, modulus.value_unchecked()));
        }

        for (ele, rhs, modulus) in izip!(
            residues_out.iter_mut(),
            self.iter_base_change_matrix(),
            self.obase.moduli()
        ) {
            *ele = modulus.reduce_dot_product(&temp, rhs);
        }
    }
}

#[cfg(test)]
mod tests {
    use primus_integer::BigIntegerOps;
    use primus_modulus::BarrettModulus;
    use rand::Rng;

    use super::*;

    type ValueT = u32;

    #[test]
    #[ignore = "just for pirint"]
    fn test_basis_convert() {
        let mut r = rand::rng();

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
            .into_iter()
            .collect();
        let dbasis = RNSBase::new(&basis).unwrap();
        println!("{:?}\n", dbasis.moduli_product());
        let mut ibasis_product = ibasis.moduli_product().to_vec();
        ibasis_product.resize(dbasis.moduli_product().len(), 0);

        for _ in 0..10 {
            let mut residues_in = Vec::with_capacity(in_len);
            for i in 0..in_len {
                residues_in.push(r.random_range(0..converter.ibase.moduli()[i].value()));
            }

            let value = ibasis.compose(&residues_in);
            println!("{:?}", value);

            let mut residues_out = vec![0; out_len];
            converter.fast_convert(&residues_in, &mut residues_out);

            residues_in.extend_from_slice(&residues_out);
            let mut value = dbasis.compose(&residues_in);
            while value.slice_cmp(&ibasis_product).is_ge() {
                let _ = value.slice_sub_assign(&ibasis_product);
            }
            println!("{:?}\n", value);
        }
    }
}
