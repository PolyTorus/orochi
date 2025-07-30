#![cfg_attr(not(feature = "std"), no_std)]
#![doc = include_str!("../README.md")]
#![warn(missing_docs)]
#![warn(rust_2018_idioms)]
#![forbid(unsafe_code)]

//! Orochi: A no_std cryptographic library for post-quantum signatures
//!
//! This library provides implementations of various post-quantum signature schemes
//! with a focus on embedded and constrained environments through no_std support.

#[cfg(not(feature = "std"))]
extern crate alloc;

#[cfg(not(feature = "std"))]
use alloc::vec::Vec;

pub mod fndsa;

// Re-export main types and functions
pub use fndsa::{
    FnDsaParams, FNDSA_512, FNDSA_1024,
    NtruPrivateKey, NtruPublicKey, generate_ntru_keypair,
    FnDsaSignature, sign_message, sign_message_with_rng,
    verify_signature, batch_verify_signatures,
    FalconTree, DiscreteGaussianSampler, CompressedSignature,
};

/// Common error types for the library
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Error {
    /// Invalid parameter provided
    InvalidParameter,
    /// Key generation failed
    KeyGenerationFailed,
    /// Signature generation failed
    SigningFailed,
    /// Signature verification failed
    VerificationFailed,
    /// Invalid signature format
    InvalidSignature,
    /// Invalid public key format
    InvalidPublicKey,
    /// Invalid secret key format
    InvalidSecretKey,
    /// Random number generation failed
    RngError,
}

impl core::fmt::Display for Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Error::InvalidParameter => write!(f, "Invalid parameter"),
            Error::KeyGenerationFailed => write!(f, "Key generation failed"),
            Error::SigningFailed => write!(f, "Signing failed"),
            Error::VerificationFailed => write!(f, "Verification failed"),
            Error::InvalidSignature => write!(f, "Invalid signature"),
            Error::InvalidPublicKey => write!(f, "Invalid public key"),
            Error::InvalidSecretKey => write!(f, "Invalid secret key"),
            Error::RngError => write!(f, "Random number generation failed"),
        }
    }
}

#[cfg(feature = "std")]
impl std::error::Error for Error {}

/// Result type alias for operations that may fail
pub type Result<T> = core::result::Result<T, Error>;