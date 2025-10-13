use std::f64::consts::{FRAC_1_SQRT_2, FRAC_2_SQRT_PI};

use primus_integer::Integer;
use rand::distr::{Distribution, Uniform};

/// Discrete Ziggurat
#[derive(Clone)]
pub struct DiscreteZiggurat<T: Integer> {
    std_dev: f64,
    x: Vec<f64>,
    y: Vec<f64>,
    sample_m: Uniform<usize>,
    sample_x: Vec<Uniform<T>>,
}

impl<T: Integer> DiscreteZiggurat<T> {
    /// Generate a [`DiscreteZiggurat<T>`]
    pub fn new(std_dev: f64, tail_cut: f64) -> Self {
        let x_m = (tail_cut * std_dev).floor();
        let sigma_square_mul_minus_two = std_dev * std_dev * (-2.0f64);

        let mut m = 3;
        'outer: loop {
            let mut x = Vec::with_capacity(m + 1);
            let mut y = Vec::with_capacity(m + 1);

            x.resize(m, 0.0);
            y.resize(m, 0.0);

            let initial_s = std_dev * FRAC_1_SQRT_2 * FRAC_2_SQRT_PI / (m as f64);

            let mut s = initial_s;
            loop {
                let mut pre_y = 0f64;
                let mut pre_x = x_m;
                for (y, x) in y.iter_mut().rev().zip(x.iter_mut().rev()) {
                    *y = s / (1.0f64 + pre_x) + pre_y;
                    *x = ((*y).ln() * sigma_square_mul_minus_two).sqrt().floor();
                    pre_y = *y;
                    pre_x = *x;
                }
                x[0] = 0f64;
                if y[0] > 1.0 {
                    break;
                }
                s += initial_s;
                if s > x_m + 1.0 {
                    m += m.next_power_of_two();
                    if m >= 512 {
                        panic!("error");
                    }
                    continue 'outer;
                }
            }
            x.push(x_m);
            y.push(0.0);
            let sample_x: Vec<Uniform<T>> = x
                .iter()
                .map(|&v| Uniform::new_inclusive(T::ZERO, T::as_from(v.floor())).unwrap())
                .collect();
            break Self {
                std_dev,
                x,
                y,
                sample_m: Uniform::new_inclusive(1, m).unwrap(),
                sample_x,
            };
        }
    }

    /// Returns the std dev of this [`DiscreteZiggurat<T>`].
    #[inline]
    pub fn std_dev(&self) -> f64 {
        self.std_dev
    }
}

impl<T: Integer> Distribution<T> for DiscreteZiggurat<T> {
    fn sample<R: rand::Rng + ?Sized>(&self, rng: &mut R) -> T {
        let pdf = |x: f64| ((x * x) / (-2.0 * self.std_dev * self.std_dev)).exp();
        let combine = |sign: bool, x: T| {
            if sign { x } else { T::ZERO - x }
        };

        loop {
            let i = self.sample_m.sample(rng);
            let sign: bool = rng.random();
            let x = self.sample_x[i].sample(rng);

            let x_f: f64 = x.as_into();

            if x_f <= self.x[i - 1] && x > T::ZERO {
                return combine(sign, x);
            } else if x == T::ZERO {
                if rng.random() {
                    return T::ZERO;
                } else {
                    continue;
                }
            } else {
                let mask = 2.0f64.powi(32);
                let y_prime = rng.next_u32();
                let y = (self.y[i - 1] - self.y[i]) * y_prime as f64;

                if self.x[i] + 1.0 <= self.std_dev {
                    if y <= mask
                        * s_line(i, self.x[i - 1], self.x[i], self.y[i - 1], self.y[i], x_f)
                        || y <= mask * (pdf(x_f) - self.y[i])
                    {
                        return combine(sign, x);
                    } else {
                        continue;
                    }
                } else if self.std_dev <= self.x[i - 1] {
                    if y >= mask
                        * s_line(
                            i,
                            self.x[i - 1],
                            self.x[i],
                            self.y[i - 1],
                            self.y[i],
                            x_f - 1.0,
                        )
                        || y > mask * (pdf(x_f) - self.y[i])
                    {
                        continue;
                    } else {
                        return combine(sign, x);
                    }
                } else if y <= mask * (pdf(x_f) - self.y[i]) {
                    return combine(sign, x);
                } else {
                    continue;
                }
            }
        }
    }
}

fn s_line(i: usize, x_i_minus_one: f64, x_i: f64, y_i_minus_one: f64, y_i: f64, x: f64) -> f64 {
    if x_i == x_i_minus_one {
        return -1.0;
    }
    if i == 1 {
        (y_i - 1.0) * (x - x_i) / (x_i - x_i_minus_one)
    } else {
        (y_i - y_i_minus_one) * (x - x_i) / (x_i - x_i_minus_one)
    }
}
