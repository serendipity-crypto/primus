use primus_integer::UnsignedInteger;

pub struct DcrtGlweKeySwitchingKey<T: UnsignedInteger> {
    key: Vec<T>,
}

impl<T: UnsignedInteger> DcrtGlweKeySwitchingKey<T> {
    // pub fn new_auto_key<R, M, Table>(
    //     sk: &CrtGlweSecretKey<T>,
    //     dcrt_sk: &DcrtGlweSecretKey<T>,
    //     degree: usize,
    //     gaussian: &SignedDiscreteGaussian<<T as UnsignedInteger>::SignedInteger>,
    //     moduli: &[M],
    //     table: &Table,
    //     rng: &mut R,
    // ) -> Self
    // where
    //     R: rand::Rng + rand::CryptoRng,
    //     M: FieldContext<T>,
    //     Table: DcrtTable<ValueT = T> + Dcrt,
    // {
    //     // let poly_length = sk.poly_length();
    //     // let dimension = sk.dimension();
    //     // let moduli_count = sk.moduli_count();

    //     todo!()
    // }
}
