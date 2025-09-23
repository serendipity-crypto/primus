use integer::UnsignedInteger;
use primus_factor::{FactorMul, ShoupFactor};
use reduce::{FieldAdapter, ops::*};
use uint_modulus::UintModulus;

use crate::{NttError, reverse::ReverseLsbs, root::PrimitiveRoot};

use super::NttTable;

/// This struct store the pre-computed data for number theory transform and
/// inverse number theory transform.
///
/// ## The structure members meet the following conditions:
///
/// 1. `n = 1 << log_n`
/// 1. `root^{n} ≡ -1 (mod modulus)`
/// 1. `root * inv_root ≡ 1 (mod modulus)`
/// 1. `n * inv_n ≡ 1 (mod modulus)`
/// 1. `root_powers` holds 1~(n-1)-th powers of root in bit-reversed order, the 0-th power is left unset.
/// 1. `inv_root_powers` holds 1~(n-1)-th powers of inverse root in scrambled order, the 0-th power is left unset.
///
/// ## Compare three orders:
///
/// ```plain
/// normal order:        0  1  2  3  4  5  6  7
///
/// bit-reversed order:  0  4  2  6  1  5  3  7
///                         -  ----  ----------
/// scrambled order:     0  1  5  3  7  2  6  4
///                         ----------  ----  -
/// ```
pub struct UintNttTable<T: UnsignedInteger> {
    root: T,
    inv_root: T,
    modulus: T,
    log_n: u32,
    n: usize,
    inv_n: ShoupFactor<T>,
    root_powers: Vec<ShoupFactor<T>>,
    inv_root_powers: Vec<ShoupFactor<T>>,
    ordinal_root_powers: Vec<ShoupFactor<T>>,
    reverse_lsbs: Vec<usize>,
    // pool: Pool<Vec<T>>,
}

impl<T: UnsignedInteger> UintNttTable<T> {
    /// Returns the root of this [`TableWithShoupRoot<T>`].
    #[inline]
    pub fn root(&self) -> T {
        self.root
    }

    /// Returns the inverse element of the root of this [`TableWithShoupRoot<T>`].
    #[inline]
    pub fn inv_root(&self) -> T {
        self.inv_root
    }

    /// Returns the modulus of this [`TableWithShoupRoot<T>`].
    #[inline]
    pub fn modulus(&self) -> T {
        self.modulus
    }

    /// Returns the log n of this [`TableWithShoupRoot<T>`].
    #[inline]
    pub fn log_n(&self) -> u32 {
        self.log_n
    }

    /// Returns the n of this [`TableWithShoupRoot<T>`].
    #[inline]
    pub fn n(&self) -> usize {
        self.n
    }

    /// Returns the inverse element of the n of this [`TableWithShoupRoot<T>`].
    #[inline]
    pub fn inv_n(&self) -> ShoupFactor<T> {
        self.inv_n
    }

    /// Returns a reference to the root powers of this [`TableWithShoupRoot<T>`].
    #[inline]
    pub fn root_powers(&self) -> &[ShoupFactor<T>] {
        &self.root_powers
    }

    /// Returns a reference to the inverse elements of the root powers of this [`TableWithShoupRoot<T>`].
    #[inline]
    pub fn inv_root_powers(&self) -> &[ShoupFactor<T>] {
        &self.inv_root_powers
    }

    /// Returns a reference to the ordinal root powers of this [`TableWithShoupRoot<T>`].
    #[inline]
    pub fn ordinal_root_powers(&self) -> &[ShoupFactor<T>] {
        &self.ordinal_root_powers
    }

    /// Returns a reference to the reverse lsbs of this [`TableWithShoupRoot<T>`].
    #[inline]
    pub fn reverse_lsbs(&self) -> &[usize] {
        &self.reverse_lsbs
    }
}

impl<T: UnsignedInteger> NttTable for UintNttTable<T> {
    type ValueT = T;

    #[inline]
    fn new<M>(log_n: u32, modulus: M) -> Result<Self, NttError<Self::ValueT>>
    where
        M: FieldAdapter<Self::ValueT>,
    {
        let root = <T as PrimitiveRoot>::try_minimal_primitive_root(log_n + 1, modulus)?;

        let modulus = modulus.value_unchecked();

        let n = 1usize << log_n;
        let to_root_type = |x| -> ShoupFactor<T> { <ShoupFactor<T>>::new(x, modulus) };

        let root_one = to_root_type(T::ONE);
        let root_factor = to_root_type(root);

        let mut power = root;

        let mut ordinal_root_powers = vec![<ShoupFactor<T>>::default(); n * 2];
        let mut iter = ordinal_root_powers.iter_mut();
        *iter.next().unwrap() = root_one;
        *iter.next().unwrap() = root_factor;
        for root_power in iter {
            power = root_factor.factor_mul_modulo(power, modulus);
            *root_power = to_root_type(power);
        }

        let inv_root = ordinal_root_powers.last().unwrap().value();

        debug_assert_eq!(root_factor.factor_mul_modulo(inv_root, modulus), T::ONE);

        let reverse_lsbs: Vec<usize> = (0..n).map(|i| i.reverse_lsbs(log_n)).collect();

        let mut root_powers = vec![<ShoupFactor<T>>::default(); n];
        root_powers[0] = root_one;
        for (&root_power, &i) in ordinal_root_powers[0..n].iter().zip(reverse_lsbs.iter()) {
            root_powers[i] = root_power;
        }

        let mut inv_root_powers = vec![<ShoupFactor<T>>::default(); n];
        inv_root_powers[0] = root_one;
        for (&inv_root_power, &i) in ordinal_root_powers[n + 1..]
            .iter()
            .rev()
            .zip(reverse_lsbs.iter())
        {
            inv_root_powers[i + 1] = inv_root_power;
        }

        let n_cast =
            T::try_from(n).map_err(|_| NttError::DegreeConversionErr { degree: n, modulus })?;

        if n_cast >= modulus {
            return Err(NttError::DegreeTooLarge { degree: n, modulus });
        }

        let inv_n = to_root_type(UintModulus(modulus).reduce_inv(n_cast));

        // let pool = Pool::new_with(2, || vec![ConstZero::ZERO; n]);

        Ok(Self {
            root,
            inv_root,
            modulus,
            log_n,
            n,
            inv_n,
            root_powers,
            inv_root_powers,
            ordinal_root_powers,
            reverse_lsbs,
            // pool,
        })
    }

    #[inline(always)]
    fn poly_length(&self) -> usize {
        self.n
    }
}
