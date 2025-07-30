//! FN-DSA (Falcon) implementation
//!
//! This module provides a complete implementation of the FN-DSA signature scheme
//! based on the Falcon algorithm, designed for post-quantum security.

pub mod params;
pub mod flr;
pub mod fft;
pub mod poly;
pub mod ntru;
pub mod falcon_tree;
pub mod gaussian;
pub mod compression;
pub mod signature;
pub mod verification;

// Re-export key types and functions
pub use params::{FnDsaParams, FNDSA_512, FNDSA_1024};
pub use ntru::{NtruPrivateKey, NtruPublicKey, generate_ntru_keypair};
pub use falcon_tree::FalconTree;
pub use gaussian::{DiscreteGaussianSampler, GpvSampler};
pub use compression::CompressedSignature;
pub use signature::{FnDsaSignature, sign_message, sign_message_with_rng};
pub use verification::{verify_signature, batch_verify_signatures, verify_signature_detailed};