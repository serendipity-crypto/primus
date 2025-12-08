use primus_integer::{AsInto, UnsignedInteger};
use rand::distr::{Distribution, StandardUniform};
use rug::{Float, az::Cast};

const PRECISION: u32 = 512;

/// UnixCDTSampler
#[derive(Debug, Clone)]
pub struct UnixCDTSampler<T: UnsignedInteger> {
    std_dev: f64,
    modulus_minus_one: T,
    cdt: Vec<rug::Integer>,
}

impl<T: UnsignedInteger> UnixCDTSampler<T> {
    /// Generate UnixCDTSampler
    pub fn new(std_dev: f64, tail_cut: f64, modulus_minus_one: T) -> Self {
        let mut length = (std_dev * tail_cut).floor() as usize + 1;

        assert!(length <= 1024);
        if length <= 1 {
            length = 2;
        }

        let std_dev_hp = Float::with_val(PRECISION, std_dev);
        let var_hp = std_dev_hp.square();

        let minus_twice_variance_recip = -(var_hp * 2u32).recip();

        let mut pdf = vec![Float::new(PRECISION); length];
        pdf[0] = Float::with_val(PRECISION, 1) / 2;

        let mut pre = minus_twice_variance_recip.clone().exp();
        pdf[1] = pre.clone();
        for i in 2..length {
            let factor = Float::with_val(PRECISION, 2 * i - 1) * &minus_twice_variance_recip;
            pre = pre * factor.exp();
            pdf[i] = pre.clone();
        }

        let s = pdf.iter().fold(Float::new(PRECISION), |acc, v| acc + v);

        let mut cdt = Vec::with_capacity(length + 1);
        let mut acc = Float::new(PRECISION);

        cdt.push(Float::new(PRECISION));
        for p in pdf.iter() {
            acc += p;
            if acc < s {
                cdt.push(Float::with_val(PRECISION, &acc / &s));
            } else {
                cdt.push(Float::with_val(PRECISION, 1));
                break;
            }
        }
        assert_eq!(cdt.len(), length + 1);

        let scalar = rug::Integer::from(1) << 256;
        let cdt: Vec<rug::Integer> = cdt
            .into_iter()
            .map(|f| {
                let t: Float = f * &scalar;
                t.cast()
            })
            .collect();

        Self {
            std_dev,
            modulus_minus_one,
            cdt,
        }
    }

    /// Returns the standard deviation of this [`UnixCDTSampler<T>`].
    #[inline]
    pub fn std_dev(&self) -> f64 {
        self.std_dev
    }
}

impl<T: UnsignedInteger> Distribution<T> for UnixCDTSampler<T> {
    #[inline]
    fn sample<R: rand::Rng + ?Sized>(&self, rng: &mut R) -> T {
        let r: [u32; 8] = StandardUniform.sample(rng);
        let positive = (r[0] & 1) == 1;
        let r = rug::Integer::from_digits(&r, rug::integer::Order::Lsf);

        let idx = self.cdt.partition_point(|x| *x <= r) - 1;
        let v: T = idx.as_into();

        if v.is_zero() {
            return T::ZERO;
        }

        if positive {
            v
        } else {
            self.modulus_minus_one - v + T::ONE
        }
    }
}
