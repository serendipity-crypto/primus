#[cfg(test)]
mod tests {
    use barrett::BarrettModulus;
    use integer::{BigIntegerOps, izip, multiply_many_values};
    use itertools::Itertools;
    use primus_decompose::big_integer::BigUintApproxSignedBasis;
    use primus_rns::RNSBase;
    use rand::Rng;
    use reduce::ops::ReduceMulAdd;

    type ValueT = u32;
    type SignedT = i64;

    #[test]
    fn test_single_decompose() {
        let mut rng = rand::rng();

        let a = BarrettModulus::new(134215681);
        let b = BarrettModulus::new(134176769);
        let moduli = &[a, b];

        let rns_base = RNSBase::new(moduli).unwrap();
        let composed_modulus: Vec<ValueT> = multiply_many_values(&[a.value(), b.value()]);
        let unused_bits = composed_modulus.last().unwrap().leading_zeros();
        let basis = BigUintApproxSignedBasis::new(&composed_modulus, 2, Some(25));

        println!("decompose_length:{}", basis.decompose_length());

        // make test simple
        assert!(basis.drop_bits() < ValueT::BITS);

        let basis_value = basis.basis_value();
        let log_basis = basis.log_basis() as usize;
        let mut decv = Vec::with_capacity(basis.decompose_length());

        let value = multiply_many_values(&[
            rng.random_range(0..a.value()),
            rng.random_range(0..b.value()),
        ]);

        let show = |value: &[ValueT]| {
            let mut value_str = String::new();
            let mut last = true;
            for &chunk in value.iter().rev() {
                if last {
                    value_str += &format!("{:01$b}", chunk, (ValueT::BITS - unused_bits) as usize);
                    last = false;
                } else {
                    value_str += &format!("{:01$b}", chunk, ValueT::BITS as usize);
                }
            }

            let (pre, end) = value_str.split_at(log_basis * basis.decompose_length());

            pre.chars().chunks(log_basis).into_iter().for_each(|v| {
                let str: String = v.collect();
                print!("{}|", str);
            });
            println!("{}", end);
        };

        println!("modulus: {:?}", composed_modulus);
        show(&composed_modulus);

        println!("value");
        show(&value);

        let (value_d, mut carry) = basis.init_value_carry(&value);

        println!("value_d");
        show(&value_d);

        for b in basis.decompose_iter() {
            let (di, ci) = b.decompose(&value_d, carry);
            decv.push(di);
            carry = ci;
        }

        let result =
            basis
                .scalar_iter()
                .zip(decv.iter())
                .fold(vec![0, 0], |mut acc, (scalar, dec)| {
                    let scalr_residue = rns_base.decompose(&scalar);
                    let dec_residue = rns_base.decompose(dec);
                    for (ac, s, d, m) in izip!(acc.iter_mut(), scalr_residue, dec_residue, moduli) {
                        *ac = m.reduce_mul_add(s, d, *ac);
                    }
                    acc
                });
        let result = rns_base.compose(&result);

        let mut cmp_value = vec![0; composed_modulus.len()];
        cmp_value[0] = basis_value / 2;
        for d in decv.iter().rev() {
            if basis_value > 2 {
                if d.slice_cmp(&cmp_value).is_ge() {
                    let mut signed = composed_modulus.clone();
                    let _ = signed.slice_sub_assign(&d);
                    print!("{:1$}|", -(signed[0] as SignedT), log_basis);
                } else {
                    print!("{:1$}|", d[0], log_basis);
                }
            } else {
                print!("{:1$}|", d[0], log_basis);
            }
        }
        println!();

        println!("value ={:?}", value);
        println!("result={:?}", result);

        let mut sub1 = result.clone();
        let _ = sub1.slice_sub_modulo_assign(&value, &composed_modulus);
        let mut sub2 = value.clone();
        let _ = sub2.slice_sub_modulo_assign(&result, &composed_modulus);
        if sub1.slice_cmp(&sub2).is_le() {
            println!("differ={:?}", sub1);
        } else {
            println!("differ={:?}", sub2);
        }
    }
}
