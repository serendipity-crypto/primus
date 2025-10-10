use primus_integer::UnsignedInteger;
use serde::{Deserialize, Serialize};

use crate::CrtRlwe;

/// A cryptographic structure for Ring Learning with Errors (RLWE).
/// This structure is used in advanced cryptographic systems and protocols, particularly
/// those that require efficient homomorphic encryption properties.
#[derive(Clone, Serialize, Deserialize)]
#[serde(bound(deserialize = "T: UnsignedInteger"))]
pub struct DcrtRlwe<T: UnsignedInteger> {
    pub(crate) data: Vec<T>,
}
