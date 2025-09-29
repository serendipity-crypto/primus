#[cfg(test)]
mod tests {
    use barrett::BarrettModulus;
    use itertools::Itertools;
    use primus_decompose::primitive::ApproxSignedBasis;
    use primus_integer::izip;
    use rand::{Rng, distr::Uniform};
    use reduce::ops::{ReduceMulAdd, ReduceSub};

    type ValueT = u32;
    type SignedT = i64;

    #[test]
    fn test_single_decompose() {
        let modulus_value: ValueT = 0b111_000_110;
        let modulus = <BarrettModulus<ValueT>>::new(modulus_value);
        let basis = ApproxSignedBasis::new(Some(modulus_value), 2, None);
        let modulus_bits = basis.value_bits() as usize;

        let difference_bound = basis.init_carry_mask().unwrap_or(0);

        let basis_value = basis.basis_value();
        let log_basis = basis.log_basis() as usize;
        let mut decv = Vec::with_capacity(basis.decompose_length());

        let mut value = 0b11_111_101;

        let show = |value: ValueT| {
            let value = if value > modulus_value.next_power_of_two() {
                value - modulus_value.next_power_of_two()
            } else {
                value
            };
            let value_str = format!("{:01$b}", value, modulus_bits);

            let (pre, end) = value_str.split_at(log_basis * basis.decompose_length());

            pre.chars().chunks(log_basis).into_iter().for_each(|v| {
                let str: String = v.collect();
                print!("{}|", str);
            });
            println!("{}", end);
        };

        println!("value");
        show(value);

        let (value_d, mut carry) = basis.init_value_carry(value);

        println!("value_d");
        show(value_d);

        for b in basis.decompose_iter() {
            let (di, ci) = b.decompose(value_d, carry);
            decv.push(di);
            carry = ci;
        }

        let result = basis
            .scalar_iter()
            .zip(decv.iter())
            .fold(0, |acc, (scalar, &dec)| {
                modulus.reduce_mul_add(scalar, dec, acc)
            });

        for &d in decv.iter().rev() {
            if basis_value > 2 {
                if d >= basis_value / 2 {
                    print!("{:1$}|", d as SignedT - modulus_value as SignedT, log_basis);
                } else {
                    print!("{:1$}|", d, log_basis);
                }
            } else {
                print!("{:1$}|", d, log_basis);
            }
        }
        println!();

        if value >= modulus_value {
            value -= modulus_value;
        }

        println!("value ={}", value);
        println!("result={}", result);

        let difference = modulus
            .reduce_sub(result, value)
            .min(modulus.reduce_sub(value, result));

        println!("difference={}", difference);
        println!("difference_bound={}", difference_bound);

        assert!(difference <= difference_bound);
    }

    #[test]
    fn test_approx_signed_decompose() {
        let mut rng = rand::rng();
        let modulus_value: ValueT = rng.random_range(512..(1 << 30));
        let modulus = <BarrettModulus<ValueT>>::new(modulus_value);
        let basis = ApproxSignedBasis::new(Some(modulus_value), 4, None);
        let modulus_bits = basis.value_bits() as usize;

        let difference_bound = basis.init_carry_mask().unwrap_or(0);

        let basis_value = basis.basis_value();
        let log_basis = basis.log_basis() as usize;
        let distr = Uniform::new(0, modulus_value).unwrap();

        let mut decv = Vec::with_capacity(basis.decompose_length());
        for value in rng.sample_iter(distr).take(1000) {
            decv.clear();

            let (value_d, mut carry) = basis.init_value_carry(value);
            for b in basis.decompose_iter() {
                let (di, ci) = b.decompose(value_d, carry);
                decv.push(di);
                carry = ci;
            }

            let result = basis
                .scalar_iter()
                .zip(decv.iter())
                .fold(0, |acc, (scalar, &dec)| {
                    modulus.reduce_mul_add(scalar, dec, acc)
                });

            let difference = modulus
                .reduce_sub(result, value)
                .min(modulus.reduce_sub(value, result));

            if difference > difference_bound {
                let show = |value: ValueT| {
                    let value_str = format!("{:01$b}", value, modulus_bits);

                    let (pre, end) = value_str.split_at(log_basis * basis.decompose_length());

                    pre.chars().chunks(log_basis).into_iter().for_each(|v| {
                        let str: String = v.collect();
                        print!("{}|", str);
                    });
                    println!("{}", end);
                };

                println!("value");
                show(value);

                for &d in decv.iter().rev() {
                    if basis_value > 2 {
                        if d >= basis_value / 2 {
                            print!("{:1$}|", d as SignedT - modulus_value as SignedT, log_basis);
                        } else {
                            print!("{:1$}|", d, log_basis);
                        }
                    } else {
                        print!("{:1$}|", d, log_basis);
                    }
                }
                println!();

                println!("value ={}", value);
                println!("result={}", result);
                println!("differ={}", difference);
                panic!("basis={}", basis_value)
            }
        }
    }

    #[test]
    fn test_decompose_slice() {
        const N: usize = 32;
        let mut rng = rand::rng();
        let modulus_value: ValueT = rng.random_range(128..(1 << 30));
        let modulus = <BarrettModulus<ValueT>>::new(modulus_value);
        let distr = Uniform::new(0, modulus_value).unwrap();

        let basis = ApproxSignedBasis::new(Some(modulus_value), 3, None);
        let differ_max = basis.init_carry_mask().unwrap_or(0);

        let input: Vec<ValueT> = rand::rng().sample_iter(distr).take(N).collect();

        let mut carries = vec![false; N];
        let mut adjust_input = input.clone();
        // basis.init_value_carry_slice_inplace(&mut adjust_input, &mut carries);
        basis.init_value_carry_slice(&input, &mut adjust_input, &mut carries);

        let mut output = vec![vec![0; N]; basis.decompose_length()];
        basis
            .decompose_iter()
            .zip(&mut output)
            .for_each(|(d, out)| d.decompose_slice_inplace(&adjust_input, out, &mut carries));

        let result = output.iter().zip(basis.scalar_iter()).fold(
            vec![0; N],
            |mut acc: Vec<ValueT>, (dec, scalar)| {
                acc.iter_mut().zip(dec.iter()).for_each(|(r, &d)| {
                    *r = modulus.reduce_mul_add(d, scalar, *r);
                });
                acc
            },
        );

        izip!(input, result).for_each(|(i, o)| {
            let difference = modulus.reduce_sub(i, o).min(modulus.reduce_sub(o, i));
            if difference > differ_max {
                println!("i ={}", i);
                println!("o ={}", o);
                println!("differ={}", difference);
            }
        });
    }
}
