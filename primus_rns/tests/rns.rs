use primus_barrett_modulus::BarrettModulus;
use primus_rns::RNSBase;

type ValueT = u16;

#[test]
fn test_rns() {
    let moduli = [3, 5, 7].map(BarrettModulus::<ValueT>::new);
    let base = RNSBase::new(&moduli).unwrap();

    let residues = &[2, 3, 2];
    let value = base.compose(residues);
    let dec = base.decompose(&value);
    assert_eq!(dec, residues);

    println!("Result: {:?}", value);
}
