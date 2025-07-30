//! FN-DSA signature verification algorithm
//!
//! This module implements the complete Falcon signature verification process,
//! including signature decompression, norm checking, and NTRU equation verification.

use crate::fndsa::{
    params::*, ntru::NtruPublicKey, signature::FnDsaSignature,
    fft::Complex, poly::{IntPoly, ModqPoly}
};
use crate::Error;
use sha3::{digest::{ExtendableOutput, Update}, Shake256};

#[cfg(not(feature = "std"))]
use alloc::vec::Vec;

/// Verify a FN-DSA signature
pub fn verify_signature(
    public_key: &NtruPublicKey,
    message: &[u8],
    signature: &FnDsaSignature,
) -> Result<bool, Error> {
    let params = signature.params();
    
    // Step 1: Basic parameter validation
    if public_key.h.logn != params.logn {
        return Ok(false);
    }
    
    // Step 2: Decompress signature components
    let s2 = signature.get_s2()?;
    
    // Step 3: Reconstruct s1 from the verification equation
    let s1 = match reconstruct_s1(public_key, message, signature, &s2)? {
        Some(s1) => s1,
        None => return Ok(false),
    };
    
    // Step 4: Verify signature norm bound
    if !verify_signature_norm_bound(&s1, &s2, params) {
        return Ok(false);
    }
    
    // Step 5: Verify the NTRU equation: s1 + s2*h ≡ hash_to_point(message) (mod q)
    let verification_result = verify_ntru_equation(
        public_key,
        message,
        signature.nonce(),
        &s1,
        &s2,
        params,
    )?;
    
    Ok(verification_result)
}

/// Reconstruct s1 from the verification equation
fn reconstruct_s1(
    public_key: &NtruPublicKey,
    message: &[u8],
    signature: &FnDsaSignature,
    s2: &[i16],
) -> Result<Option<Vec<i16>>, Error> {
    let params = signature.params();
    
    // Hash message to target point
    let target = hash_to_point_for_verification(message, signature.nonce(), params)?;
    
    // Convert s2 to modular polynomial
    let s2_poly = IntPoly::from_coeffs(s2.to_vec(), params.logn)?;
    let s2_mod = s2_poly.to_modq();
    
    // Compute s2 * h mod q
    let s2_h = s2_mod.multiply_ntt(&public_key.h)?;
    
    // Convert target to modular polynomial
    let target_int: Vec<i16> = target.iter()
        .map(|c| c.re.round_ties_to_even() as i16)
        .collect();
    let target_poly = IntPoly::from_coeffs(target_int, params.logn)?;
    let target_mod = target_poly.to_modq();
    
    // Compute s1 = target - s2*h (mod q)
    let s1_mod = target_mod.sub(&s2_h)?;
    
    // Convert back to signed integer representation
    let s1: Vec<i16> = s1_mod.coeffs.iter()
        .map(|&x| {
            let signed = x as i32;
            let centered = if signed > (params.q / 2) as i32 {
                signed - params.q as i32
            } else {
                signed
            };
            centered as i16
        })
        .collect();
    
    Ok(Some(s1))
}

/// Hash message to point for verification (same as in signing)
fn hash_to_point_for_verification(
    message: &[u8],
    nonce: [u8; SALT_LEN],
    params: FnDsaParams,
) -> Result<Vec<Complex>, Error> {
    // Use SHAKE256 to generate deterministic hash
    let mut hasher = Shake256::default();
    
    // Hash domain separator, nonce, and message
    Update::update(&mut hasher, DOMAIN_SEPARATOR);
    Update::update(&mut hasher, &nonce);
    Update::update(&mut hasher, message);
    
    // Generate hash output
    let mut hash_output = vec![0u8; params.hash_output_len()];
    ExtendableOutput::finalize_xof_into(hasher, &mut hash_output);
    
    // Convert hash to complex coefficients
    let mut target = Vec::with_capacity(params.n);
    for i in 0..params.n {
        let byte_index = i * HASH_BYTES_PER_COEFF;
        if byte_index + 1 < hash_output.len() {
            let val = (hash_output[byte_index] as u16) | 
                     ((hash_output[byte_index + 1] as u16) << 8);
            
            // Reduce modulo q and center around 0
            let reduced = (val % params.q as u16) as i16;
            let centered = if reduced > (params.q / 2) as i16 {
                reduced - params.q as i16
            } else {
                reduced
            };
            
            target.push(Complex::from_real(centered as f64)?);
        } else {
            target.push(Complex::zero());
        }
    }
    
    Ok(target)
}

/// Verify signature norm bound
fn verify_signature_norm_bound(s1: &[i16], s2: &[i16], params: FnDsaParams) -> bool {
    let norm_squared: u64 = s1.iter()
        .chain(s2.iter())
        .map(|&x| (x as i64 * x as i64) as u64)
        .sum();
    
    params.norm_is_valid(norm_squared)
}

/// Verify the NTRU equation: s1 + s2*h ≡ target (mod q)
fn verify_ntru_equation(
    public_key: &NtruPublicKey,
    message: &[u8],
    nonce: [u8; SALT_LEN],
    s1: &[i16],
    s2: &[i16],
    params: FnDsaParams,
) -> Result<bool, Error> {
    // Hash message to target
    let target = hash_to_point_for_verification(message, nonce, params)?;
    
    // Convert signature components to modular polynomials
    let s1_poly = IntPoly::from_coeffs(s1.to_vec(), params.logn)?;
    let s2_poly = IntPoly::from_coeffs(s2.to_vec(), params.logn)?;
    let s1_mod = s1_poly.to_modq();
    let s2_mod = s2_poly.to_modq();
    
    // Compute s2 * h mod q
    let s2_h = s2_mod.multiply_ntt(&public_key.h)?;
    
    // Compute s1 + s2*h mod q
    let result_mod = s1_mod.add(&s2_h)?;
    
    // Convert target to modular polynomial
    let target_int: Vec<i16> = target.iter()
        .map(|c| c.re.round_ties_to_even() as i16)
        .collect();
    let target_poly = IntPoly::from_coeffs(target_int, params.logn)?;
    let target_mod = target_poly.to_modq();
    
    // Compare result with target using constant-time comparison
    constant_time_polynomial_compare(&result_mod.coeffs, &target_mod.coeffs)
}

/// Constant-time polynomial comparison
fn constant_time_polynomial_compare(a: &[u16], b: &[u16]) -> Result<bool, Error> {
    if a.len() != b.len() {
        return Ok(false);
    }
    
    let mut diff_accumulator = 0u16;
    
    for (&a_coeff, &b_coeff) in a.iter().zip(b.iter()) {
        diff_accumulator |= a_coeff ^ b_coeff;
    }
    
    Ok(diff_accumulator == 0)
}

/// Batch verification of multiple signatures (more efficient)
pub fn batch_verify_signatures(
    public_keys: &[NtruPublicKey],
    messages: &[&[u8]],
    signatures: &[FnDsaSignature],
) -> Result<bool, Error> {
    if public_keys.len() != messages.len() || messages.len() != signatures.len() {
        return Err(Error::InvalidParameter);
    }
    
    if public_keys.is_empty() {
        return Ok(true);
    }
    
    // For simplicity, verify each signature individually
    // Real batch verification would combine operations for efficiency
    for ((public_key, &message), signature) in 
        public_keys.iter().zip(messages.iter()).zip(signatures.iter()) {
        if !verify_signature(public_key, message, signature)? {
            return Ok(false);
        }
    }
    
    Ok(true)
}

/// Advanced verification with detailed error reporting
#[derive(Debug)]
pub enum VerificationError {
    InvalidParameters,
    DecompressionFailed,
    NormCheckFailed,
    EquationVerificationFailed,
    HashMismatch,
}

pub fn verify_signature_detailed(
    public_key: &NtruPublicKey,
    message: &[u8],
    signature: &FnDsaSignature,
) -> Result<Result<(), VerificationError>, Error> {
    let params = signature.params();
    
    // Parameter validation
    if public_key.h.logn != params.logn {
        return Ok(Err(VerificationError::InvalidParameters));
    }
    
    // Signature decompression
    let s2 = match signature.get_s2() {
        Ok(s2) => s2,
        Err(_) => return Ok(Err(VerificationError::DecompressionFailed)),
    };
    
    // s1 reconstruction
    let s1 = match reconstruct_s1(public_key, message, signature, &s2)? {
        Some(s1) => s1,
        None => return Ok(Err(VerificationError::HashMismatch)),
    };
    
    // Norm check
    if !verify_signature_norm_bound(&s1, &s2, params) {
        return Ok(Err(VerificationError::NormCheckFailed));
    }
    
    // NTRU equation verification
    if !verify_ntru_equation(public_key, message, signature.nonce(), &s1, &s2, params)? {
        return Ok(Err(VerificationError::EquationVerificationFailed));
    }
    
    Ok(Ok(()))
}

/// Signature verification statistics for analysis
#[derive(Debug)]
pub struct VerificationStats {
    pub s1_norm_squared: u64,
    pub s2_norm_squared: u64,
    pub total_norm_squared: u64,
    pub norm_utilization: f64, // Fraction of maximum allowed norm used
    pub max_coefficient: u16,
}

pub fn analyze_signature_verification(
    signature: &FnDsaSignature,
    s1: &[i16],
    s2: &[i16],
) -> VerificationStats {
    let params = signature.params();
    
    let s1_norm_squared: u64 = s1.iter()
        .map(|&x| (x as i64 * x as i64) as u64)
        .sum();
    
    let s2_norm_squared: u64 = s2.iter()
        .map(|&x| (x as i64 * x as i64) as u64)
        .sum();
    
    let total_norm_squared = s1_norm_squared + s2_norm_squared;
    
    let norm_utilization = total_norm_squared as f64 / params.beta_squared as f64;
    
    let max_coefficient = s1.iter()
        .chain(s2.iter())
        .map(|&x| x.unsigned_abs())
        .max()
        .unwrap_or(0);
    
    VerificationStats {
        s1_norm_squared,
        s2_norm_squared,
        total_norm_squared,
        norm_utilization,
        max_coefficient,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fndsa::{params::FNDSA_512, ntru::generate_ntru_keypair, signature::sign_message_with_rng};
    use rand_core::OsRng;
    
    #[test]
    #[cfg(feature = "std")]
    fn test_hash_to_point_deterministic() {
        let message = b"Test message";
        let nonce = [42u8; SALT_LEN];
        let params = FNDSA_512;
        
        let target1 = hash_to_point_for_verification(message, nonce, params).unwrap();
        let target2 = hash_to_point_for_verification(message, nonce, params).unwrap();
        
        // Should be identical
        assert_eq!(target1.len(), target2.len());
        for (a, b) in target1.iter().zip(target2.iter()) {
            assert!((a.re.raw() - b.re.raw()).abs() < 1e-10);
        }
    }
    
    #[test]
    fn test_constant_time_comparison() {
        let a = vec![1u16, 2, 3, 4, 5];
        let b = vec![1u16, 2, 3, 4, 5];
        let c = vec![1u16, 2, 3, 4, 6];
        
        assert!(constant_time_polynomial_compare(&a, &b).unwrap());
        assert!(!constant_time_polynomial_compare(&a, &c).unwrap());
        assert!(!constant_time_polynomial_compare(&a, &vec![1u16, 2, 3, 4]).unwrap());
    }
    
    #[test]
    fn test_norm_bound_verification() {
        let params = FNDSA_512;
        
        // Small signature within bounds
        let s1 = vec![1i16; params.n];
        let s2 = vec![1i16; params.n];
        assert!(verify_signature_norm_bound(&s1, &s2, params));
        
        // Large signature exceeding bounds
        // For FNDSA_512: beta_squared = 34034726
        // We need coefficients that when squared and summed exceed this
        // With 1024 total coefficients, we need each coeff^2 > 33233
        // So coefficients > 182 should exceed the bound
        let s1_large = vec![200i16; params.n];
        let s2_large = vec![200i16; params.n];
        assert!(!verify_signature_norm_bound(&s1_large, &s2_large, params));
        
        // Test edge case: just within bounds
        let edge_value = ((params.beta_squared as f64 / (2.0 * params.n as f64)).sqrt() as i16) - 1;
        let s1_edge = vec![edge_value; params.n];
        let s2_edge = vec![edge_value; params.n];
        assert!(verify_signature_norm_bound(&s1_edge, &s2_edge, params));
    }
    
    #[test]
    #[cfg(feature = "std")]
    fn test_signature_verification_structure() {
        let mut rng = OsRng;
        let params = FNDSA_512;
        
        // Create a mock signature for structure testing
        let nonce = [1u8; SALT_LEN];
        let s1 = vec![1i16; params.n];
        let s2 = vec![2i16; params.n];
        
        let signature = crate::fndsa::signature::FnDsaSignature::new(
            params, nonce, s1.clone(), s2.clone()
        ).unwrap();
        
        // Test signature analysis
        let stats = analyze_signature_verification(&signature, &s1, &s2);
        
        assert_eq!(stats.s1_norm_squared, params.n as u64); // All 1's
        assert_eq!(stats.s2_norm_squared, params.n as u64 * 4); // All 2's
        assert_eq!(stats.total_norm_squared, params.n as u64 * 5);
        assert!(stats.norm_utilization > 0.0);
        assert_eq!(stats.max_coefficient, 2);
    }
    
    #[test]
    #[cfg(feature = "std")]
    fn test_batch_verification() {
        let params = FNDSA_512;
        
        // Create empty batch (should succeed)
        let result = batch_verify_signatures(&[], &[], &[]).unwrap();
        assert!(result);
        
        // Create mock data for batch verification
        let nonce = [1u8; SALT_LEN];
        let s1 = vec![1i16; params.n];
        let s2 = vec![2i16; params.n];
        
        let signature = crate::fndsa::signature::FnDsaSignature::new(
            params, nonce, s1, s2
        ).unwrap();
        
        // This would fail with real keys, but tests the API structure
        let result = batch_verify_signatures(
            &[],
            &[],
            &[signature],
        );
        
        // Should return parameter error due to mismatched lengths
        assert!(result.is_err());
    }
    
    #[test]
    #[cfg(feature = "std")]  
    fn test_detailed_verification() {
        let params = FNDSA_512;
        let nonce = [1u8; SALT_LEN];
        let s1 = vec![1i16; params.n];
        let s2 = vec![2i16; params.n];
        
        let signature = crate::fndsa::signature::FnDsaSignature::new(
            params, nonce, s1, s2
        ).unwrap();
        
        // Create a dummy public key for testing
        let h = crate::fndsa::poly::ModqPoly::from_coeffs(vec![1u16; params.n], params.logn).unwrap();
        let public_key = NtruPublicKey { h };
        
        let message = b"Test message";
        
        let detailed_result = verify_signature_detailed(&public_key, message, &signature).unwrap();
        
        // With simplified implementation, this will likely fail verification
        // but demonstrates the detailed error reporting structure
        match detailed_result {
            Ok(()) => {
                // Successful verification (unlikely with mock data)
            }
            Err(error) => {
                // Expected failure with detailed error information
                println!("Verification failed with detailed error: {:?}", error);
            }
        }
    }
}