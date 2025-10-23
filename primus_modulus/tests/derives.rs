macro_rules! test_modulus {
    () => {
        fn field_trait<M: FieldContext<ValueT>>(_modulus: M) {}

        #[test]
        fn test_modulus() {
            field_trait(Modulus);
            let modulus_value = Modulus::value();
            let distribution = Uniform::new(0, modulus_value).unwrap();

            let mut rng = rand::rng();
            let a = distribution.sample(&mut rng);
            let b = distribution.sample(&mut rng);
            let c = distribution.sample(&mut rng);

            let d = Modulus.reduce_add(a, b);
            assert_eq!(d, (a + b) % modulus_value);

            let d = Modulus.reduce_sub(a, b);
            assert_eq!(d, (a + modulus_value - b) % modulus_value);

            let d = Modulus.reduce_neg(a);
            assert_eq!(0, Modulus.reduce_add(a, d));

            let d = Modulus.reduce_mul(a, b);
            assert_eq!(
                d,
                (a as WideT * b as WideT % modulus_value as WideT) as ValueT
            );

            let d = Modulus.reduce_square(a);
            assert_eq!(
                d,
                (a as WideT * a as WideT % modulus_value as WideT) as ValueT
            );

            let d = Modulus.reduce_mul_add(a, b, c);
            assert_eq!(
                d,
                ((a as WideT * b as WideT + c as WideT) % modulus_value as WideT) as ValueT
            );

            if a != 0 {
                let d = Modulus.reduce_inv(a);
                assert_eq!(1, Modulus.reduce_mul(a, d));
            }

            if b != 0 {
                let d = Modulus.reduce_div(a, b);
                assert_eq!(a, Modulus.reduce_mul(b, d));
            }
        }
    };
}

#[cfg(all(test, feature = "derive"))]
mod u8tests {
    use primus_modulus::Barrett;
    use primus_reduce::FieldContext;
    use primus_reduce::ops::*;
    use rand::{distr::Uniform, prelude::*};

    #[derive(Barrett)]
    #[modulus(u8, value = 61)]
    struct Modulus;

    type ValueT = u8;
    type WideT = u16;

    test_modulus!();
}

#[cfg(all(test, feature = "derive"))]
mod u16tests {
    use primus_modulus::Barrett;
    use primus_reduce::FieldContext;
    use primus_reduce::ops::*;
    use rand::{distr::Uniform, prelude::*};

    #[derive(Barrett)]
    #[modulus(u16, value = 12289)]
    struct Modulus;

    type ValueT = u16;
    type WideT = u32;

    test_modulus!();
}

#[cfg(all(test, feature = "derive"))]
mod u32tests {
    use primus_modulus::Barrett;
    use primus_reduce::FieldContext;
    use primus_reduce::ops::*;
    use rand::{distr::Uniform, prelude::*};

    #[derive(Barrett)]
    #[modulus(u32, value = 536813569)]
    struct Modulus;

    type ValueT = u32;
    type WideT = u64;

    test_modulus!();
}

#[cfg(all(test, feature = "derive"))]
mod u64tests {
    use primus_modulus::Barrett;
    use primus_reduce::FieldContext;
    use primus_reduce::ops::*;
    use rand::{distr::Uniform, prelude::*};

    #[derive(Barrett)]
    #[modulus(u64, value = 4611686018427322369)]
    struct Modulus;

    type ValueT = u64;
    type WideT = u128;

    test_modulus!();
}
