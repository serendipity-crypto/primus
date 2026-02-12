#[cfg(test)]
mod tests {
    use itertools::Itertools;
    use primus_decompose::zama::ApproxSignedBasis;
    use primus_modulus::BarrettModulus;
    use primus_reduce::ops::{ReduceMulAdd, ReduceSub};
    use rand::{RngExt, distr::Uniform};

    type ValueT = u32;
    type SignedT = i64;

    #[test]
    fn test_single_decompose() {
        let modulus_value: ValueT = 0b111_000_110;
        let modulus = <BarrettModulus<ValueT>>::new(modulus_value);

        let basis = ApproxSignedBasis::new(modulus_value, 3, Some(2));
        let log_basis = basis.log_basis() as usize;

        let decomposed_bits = basis.decomposed_bits() as usize;
        let decompose_length = basis.decompose_length();

        let difference_bound = 1 << (basis.drop_length() * log_basis);

        let basis_value = basis.basis_value();
        let mut decv = Vec::with_capacity(basis.drop_length() + decompose_length);

        let value = rand::rng().random_range(0..modulus_value);

        let show = |value: ValueT| {
            let value_str = format!("{:01$b}", value, decomposed_bits);

            let (pre, end) = value_str.split_at(log_basis * decompose_length);

            pre.chars().chunks(log_basis).into_iter().for_each(|v| {
                let str: String = v.collect();
                print!("{}_", str);
            });

            end.chars().chunks(log_basis).into_iter().for_each(|v| {
                let str: String = v.collect();
                print!("{}_", str);
            });
            println!();
        };

        println!("value");
        show(value);

        let (value_d, sign) = basis.init_value_sign(value);

        println!("value_d");
        show(value_d);

        let mut carry = false;
        for b in basis.drop_iter() {
            let (di, ci) = b.decompose(value_d, sign, carry);
            decv.push(di);
            carry = ci;
        }

        for b in basis.decompose_iter() {
            let (di, ci) = b.decompose(value_d, sign, carry);
            decv.push(di);
            carry = ci;
        }

        let result = basis
            .scalar_iter()
            .zip(decv.iter().skip(basis.drop_length()))
            .fold(0, |acc, (scalar, &dec)| {
                modulus.reduce_mul_add(scalar, dec, acc)
            });

        for &d in decv.iter().rev() {
            if basis_value > 2 {
                if d > basis_value {
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
        let basis = ApproxSignedBasis::new(modulus_value, 3, None);
        let log_basis = basis.log_basis() as usize;

        let decomposed_bits = basis.decomposed_bits() as usize;
        let drop_length = basis.drop_length();
        let decompose_length = basis.decompose_length();

        let difference_bound = 1 << (drop_length * log_basis);

        let basis_value = basis.basis_value();
        let distr = Uniform::new(0, modulus_value).unwrap();

        let mut decv = Vec::with_capacity(drop_length + decompose_length);
        for value in rng.sample_iter(distr).take(1000) {
            decv.clear();

            let (value_d, sign) = basis.init_value_sign(value);

            let mut carry = false;

            for b in basis.drop_iter() {
                let (di, ci) = b.decompose(value_d, sign, carry);
                decv.push(di);
                carry = ci;
            }

            for b in basis.decompose_iter() {
                let (di, ci) = b.decompose(value_d, sign, carry);
                decv.push(di);
                carry = ci;
            }

            let result = basis
                .scalar_iter()
                .zip(decv.iter().skip(drop_length))
                .fold(0, |acc, (scalar, &dec)| {
                    modulus.reduce_mul_add(scalar, dec, acc)
                });

            let difference = modulus
                .reduce_sub(result, value)
                .min(modulus.reduce_sub(value, result));

            if difference > difference_bound {
                let show = |value: ValueT| {
                    let value_str = format!("{:01$b}", value, decomposed_bits);

                    let (pre, end) = value_str.split_at(log_basis * decompose_length);

                    pre.chars().chunks(log_basis).into_iter().for_each(|v| {
                        let str: String = v.collect();
                        print!("{}|", str);
                    });
                    println!("{}", end);
                };

                println!("value");
                show(value);

                println!("value_d");
                show(value_d);

                for &d in decv.iter().rev() {
                    if basis_value > 2 {
                        if d > basis_value {
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

    // #[test]
    // fn test_decompose_slice() {
    //     const N: usize = 32;
    //     let mut rng = rand::rng();
    //     let modulus_value: ValueT = rng.random_range(128..(1 << 30));
    //     let modulus = <BarrettModulus<ValueT>>::new(modulus_value);
    //     let distr = Uniform::new(0, modulus_value).unwrap();

    //     let basis = ApproxSignedBasis::new(Some(modulus_value), 3, None);
    //     let differ_max = basis.init_carry_mask().unwrap_or(0);

    //     let input: Vec<ValueT> = rand::rng().sample_iter(distr).take(N).collect();

    //     let mut carries = vec![false; N];
    //     let mut adjust_input = input.clone();
    //     // basis.init_value_carry_slice_inplace(&mut adjust_input, &mut carries);
    //     basis.init_value_carry_slice(&input, &mut adjust_input, &mut carries);

    //     let mut output = vec![vec![0; N]; basis.decompose_length()];
    //     basis
    //         .decompose_iter()
    //         .zip(&mut output)
    //         .for_each(|(d, out)| d.decompose_slice_inplace(&adjust_input, out, &mut carries));

    //     let result = output.iter().zip(basis.scalar_iter()).fold(
    //         vec![0; N],
    //         |mut acc: Vec<ValueT>, (dec, scalar)| {
    //             acc.iter_mut().zip(dec.iter()).for_each(|(r, &d)| {
    //                 *r = modulus.reduce_mul_add(d, scalar, *r);
    //             });
    //             acc
    //         },
    //     );

    //     izip!(input, result).for_each(|(i, o)| {
    //         let difference = modulus.reduce_sub(i, o).min(modulus.reduce_sub(o, i));
    //         if difference > differ_max {
    //             println!("i ={}", i);
    //             println!("o ={}", o);
    //             println!("differ={}", difference);
    //         }
    //     });
    // }
}
