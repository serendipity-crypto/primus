use primus_integer::UnsignedInteger;

pub struct DcrtGlevContext<T: UnsignedInteger> {
    adjust_big_uint_values: Vec<T>,
    decomposed_unsigned_values: Vec<T>,
    carries: Vec<bool>,
    multi_residues: Vec<T>,
}

impl<T: UnsignedInteger> DcrtGlevContext<T> {
    pub fn new(poly_length: usize, crt_poly_len: usize, big_uint_poly_len: usize) -> Self {
        let adjust_big_uint_values = vec![T::ZERO; big_uint_poly_len];
        let decomposed_unsigned_values = vec![T::ZERO; poly_length];
        let carries = vec![false; poly_length];
        let multi_residues = vec![T::ZERO; crt_poly_len];

        Self {
            adjust_big_uint_values,
            decomposed_unsigned_values,
            carries,
            multi_residues,
        }
    }

    #[inline]
    pub fn as_mut(&mut self) -> (&mut [T], &mut [T], &mut [bool], &mut [T]) {
        (
            self.adjust_big_uint_values.as_mut(),
            self.decomposed_unsigned_values.as_mut(),
            self.carries.as_mut(),
            self.multi_residues.as_mut(),
        )
    }
}
