// Quick test for various edge case sigma values
use primus_distr::DiscreteZiggurat;

fn main() {
    // Test various edge cases
    let test_cases = vec![
        (0.5, "Very small sigma"),
        (1.0, "Unit sigma"),
        (1.9, "Small sigma"),
        (3.19, "Previously failing sigma"),
        (5.0, "Medium sigma"),
        (10.0, "Large sigma"),
        (100.0, "Very large sigma"),
    ];

    println!("Testing DiscreteZiggurat initialization with various sigma values:\n");

    for (sigma, desc) in test_cases {
        print!("{:.<30} (σ={:>6}): ", desc, sigma);
        match std::panic::catch_unwind(|| DiscreteZiggurat::<u64>::new(sigma, 12.0, u64::MAX)) {
            Ok(_) => println!("✓ Success"),
            Err(_) => println!("✗ Failed"),
        }
    }
}
