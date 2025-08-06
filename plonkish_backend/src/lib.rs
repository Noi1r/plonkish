#![allow(clippy::op_ref)]

pub mod anemoi_hash;
pub mod backend;
pub mod frontend;
pub mod pcs;
pub mod piop;
pub mod poly;
pub mod transform;
pub mod util;

pub use halo2_curves;

#[derive(Clone, Debug, PartialEq)]
pub enum Error {
    InvalidSumcheck(String),
    InvalidPcsParam(String),
    InvalidPcsOpen(String),
    InvalidSnark(String),
    Serialization(String),
    Transcript(std::io::ErrorKind, String),
    NotImplemented(String),
    InternalError(String),
    InvalidRotation(String),
    InvalidQuotient(String),
    InvalidInput(String),
}
