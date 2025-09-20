use barrett::BarrettModulus;
use primus_rns::RNSBase;

type ValueT = u16;

#[test]
fn test_rns() {
    let moduli = [3, 5, 7].map(BarrettModulus::<ValueT>::new);
    let base = RNSBase::new(&moduli).unwrap();

    let result = base.compose(&[2, 3, 2]);

    println!("Result: {:?}", result);
}
