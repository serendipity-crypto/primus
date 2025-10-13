/// Lwe Ciphertext
pub type LweCiphertext<T> = primus_lattice::Lwe<T>;

/// CmLwe Ciphertext
pub type MultiMsgLweCiphertext<T> = primus_lattice::MultiMsgLwe<T>;

/// Rlwe Ciphertext
pub type RlweCiphertext<T> = primus_lattice::Rlwe<T>;

/// Ntt version Rlwe Ciphertext
pub type NttRlweCiphertext<T> = primus_lattice::NttRlwe<T>;

/// Glwe Ciphertext
pub type GlweCiphertext<T> = primus_lattice::Glwe<T>;

/// Ntt version Glwe Ciphertext
pub type NttGlweCiphertext<T> = primus_lattice::NttGlwe<T>;

/// Glwe Ciphertext
pub type CrtGlweCiphertext<T> = primus_lattice::CrtGlwe<T>;

/// Ntt version Glwe Ciphertext
pub type DcrtGlweCiphertext<T> = primus_lattice::DcrtGlwe<T>;
