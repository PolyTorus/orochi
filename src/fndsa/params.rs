//! FN-DSA parameter sets and constants
//! 
//! This module contains the exact parameter sets for FN-DSA as specified
//! in the NIST post-quantum cryptography standardization process.

/// Logarithmic degree for FN-DSA-512
pub const LOGN_512: usize = 9;

/// Logarithmic degree for FN-DSA-1024  
pub const LOGN_1024: usize = 10;

/// Prime modulus q used in both parameter sets
pub const Q: u32 = 12289;

/// FN-DSA parameter set specification
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FnDsaParams {
    /// Logarithmic degree (logn)
    pub logn: usize,
    /// Polynomial degree (n = 2^logn)
    pub n: usize,
    /// Prime modulus
    pub q: u32,
    /// Signature bound squared (β²)
    pub beta_squared: u64,
    /// NIST security level in bits
    pub security_level: usize,
    /// Public key size in bytes
    pub public_key_len: usize,
    /// Private key size in bytes (implementation dependent)
    pub private_key_len: usize,
    /// Maximum signature size in bytes
    pub signature_len: usize,
}

/// FN-DSA-512 parameter set (NIST security level 1)
pub const FNDSA_512: FnDsaParams = FnDsaParams {
    logn: LOGN_512,
    n: 1 << LOGN_512,        // 512
    q: Q,
    beta_squared: 34034726,   // From NIST specification
    security_level: 128,
    public_key_len: 897,
    private_key_len: 1281,
    signature_len: 666,
};

/// FN-DSA-1024 parameter set (NIST security level 5)
pub const FNDSA_1024: FnDsaParams = FnDsaParams {
    logn: LOGN_1024,
    n: 1 << LOGN_1024,       // 1024
    q: Q,
    beta_squared: 70265242,   // From NIST specification
    security_level: 256,
    public_key_len: 1793,
    private_key_len: 2305,
    signature_len: 1280,
};

/// Maximum supported degree
pub const MAX_LOGN: usize = 10;

/// Maximum polynomial degree
pub const MAX_N: usize = 1 << MAX_LOGN;

/// Salt size for hash-to-point in bytes
pub const SALT_LEN: usize = 40;

/// Hash output length factor (2 bytes per coefficient)
pub const HASH_BYTES_PER_COEFF: usize = 2;

/// Coefficient bound for signature compression
pub const COEFF_BOUND: i16 = 2047;

/// NTT primitive root for q = 12289
/// This is a primitive 2^11-th root of unity modulo q
pub const NTT_ROOT: u32 = 49;

/// Inverse of n modulo q for each supported degree
pub const INV_N_512: u32 = 12277;  // Precomputed: n^(-1) mod q for n=512
pub const INV_N_1024: u32 = 12265; // Precomputed: n^(-1) mod q for n=1024

/// Maximum attempts for key generation before giving up
pub const MAX_KEYGEN_ATTEMPTS: usize = 1000;

/// Gaussian sampling bound multiplier (6σ bound)
pub const GAUSSIAN_BOUND_MULTIPLIER: f64 = 6.0;

/// Hash domain separation for challenge generation
pub const DOMAIN_SEPARATOR: &[u8] = b"FN-DSA";

/// Implementation-specific constants
impl FnDsaParams {
    /// Get the appropriate parameter set for the given logarithmic degree
    pub const fn from_logn(logn: usize) -> Option<Self> {
        match logn {
            LOGN_512 => Some(FNDSA_512),
            LOGN_1024 => Some(FNDSA_1024),
            _ => None,
        }
    }
    
    /// Get the inverse of n modulo q for NTT operations
    pub const fn inv_n(self) -> u32 {
        match self.logn {
            LOGN_512 => INV_N_512,
            LOGN_1024 => INV_N_1024,
            _ => 0, // Invalid
        }
    }
    
    /// Get the signature bound β (square root of beta_squared)
    pub fn beta(self) -> f64 {
        (self.beta_squared as f64).sqrt()
    }
    
    /// Check if a given norm squared is within the signature bound
    pub fn norm_is_valid(self, norm_squared: u64) -> bool {
        norm_squared <= self.beta_squared
    }
    
    /// Calculate the expected hash output length
    pub const fn hash_output_len(self) -> usize {
        self.n * HASH_BYTES_PER_COEFF
    }
}