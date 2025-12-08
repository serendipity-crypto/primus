use std::{marker::PhantomData, num::NonZero, sync::LazyLock};

use bigdecimal::{BigDecimal, Context, RoundingMode};
use num_traits::{FromPrimitive, One, ToPrimitive, Zero};
use primus_integer::{AsInto, Integer};
use rand::distr::{Distribution, Uniform};

/// cdt sampler
#[derive(Debug, Clone)]
pub struct CDTSampler<T: Integer> {
    std_dev: f64,
    cdt: Vec<u128>,
    phantom: PhantomData<T>,
}

impl<T: Integer> CDTSampler<T> {
    /// Generate a [`CDTSampler<T>`]
    pub fn new(std_dev: f64, tail_cut: f64) -> Self {
        let max_std_dev = std_dev * tail_cut;
        let mut length = max_std_dev.floor() as usize + 1;

        assert!(length <= 1024);
        if length <= 1 {
            length = 2;
        }

        let precision = if std_dev < 2.0 {
            128
        } else if std_dev < 5.0 {
            192
        } else {
            256
        };

        let context = Context::new(NonZero::new(precision).unwrap(), RoundingMode::HalfUp);

        let std_dev_hp = BigDecimal::from_f64(std_dev).unwrap();
        let var_hp = std_dev_hp.square();

        let minus_twice_variance_recip = -var_hp.double().inverse_with_context(&context);

        let mut pdf = vec![BigDecimal::default(); length];
        pdf[0] = BigDecimal::one().half();

        let mut pre = minus_twice_variance_recip.exp();
        pdf[1] = pre.clone();
        for i in 2..length {
            let factor = BigDecimal::from_usize(2 * i - 1).unwrap() * &minus_twice_variance_recip;
            pre = pre * factor.exp();
            pdf[i] = pre.clone();
        }

        let s = pdf.iter().fold(BigDecimal::zero(), |acc, v| acc + v);

        let mut cdt = Vec::with_capacity(length + 1);
        let mut acc = BigDecimal::zero();

        cdt.push(BigDecimal::zero());
        for p in pdf.iter() {
            acc += p;

            if acc < s {
                cdt.push(&acc / &s);
            } else {
                cdt.push(BigDecimal::one());
                break;
            }
        }
        assert_eq!(cdt.len(), length + 1);

        let cdt: Vec<u128> = cdt
            .into_iter()
            .map(|f| {
                (f * u128::MAX)
                    .with_scale_round(0, RoundingMode::HalfUp)
                    .to_u128()
                    .unwrap()
            })
            .collect();

        Self {
            std_dev,
            cdt,
            phantom: PhantomData::<T>,
        }
    }

    /// Returns the standard deviation of this [`CDTSampler<T>`].
    #[inline]
    pub fn std_dev(&self) -> f64 {
        self.std_dev
    }
}

static D: LazyLock<Uniform<u128>> = LazyLock::new(|| Uniform::new_inclusive(0, u128::MAX).unwrap());

impl<T: Integer> Distribution<T> for CDTSampler<T> {
    #[inline]
    fn sample<R: rand::Rng + ?Sized>(&self, rng: &mut R) -> T {
        let r: u128 = D.sample(rng);

        let positive = (r & 1) == 1;

        let idx = self.cdt.partition_point(|&x| x <= r) - 1;

        let v: T = idx.as_into();

        if v.is_zero() {
            return T::ZERO;
        }

        if positive { v } else { T::ZERO - v }
    }
}
