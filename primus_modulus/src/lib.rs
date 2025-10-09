pub use primus_integer as integer;
pub use primus_reduce as reduce;

mod barrett;
mod native;
mod power_of_two;
mod unsigned_integer;

#[cfg(feature = "derive")]
pub use primus_barrett_derive::Barrett;

pub use barrett::BarrettModulus;
pub use native::NativeModulus;
pub use power_of_two::PowOf2Modulus;
pub use unsigned_integer::UintModulus;
