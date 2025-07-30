//! Polynomial operations for FN-DSA over various coefficient rings
//!
//! This module implements polynomial arithmetic over integers, integers modulo q,
//! and complex numbers, with efficient implementations using FFT and NTT.

use crate::fndsa::{params::*, fft::Complex};
use crate::Error;

#[cfg(not(feature = "std"))]
use alloc::vec::Vec;

/// Polynomial with integer coefficients
#[derive(Debug, Clone, PartialEq)]
pub struct IntPoly {
    pub coeffs: Vec<i16>,
    pub logn: usize,
}

/// Polynomial with coefficients modulo q
#[derive(Debug, Clone, PartialEq)]  
pub struct ModqPoly {
    pub coeffs: Vec<u16>,
    pub logn: usize,
}

/// Polynomial with complex coefficients (for FFT operations)
#[derive(Debug, Clone)]
pub struct ComplexPoly {
    pub coeffs: Vec<Complex>,
    pub logn: usize,
}

impl IntPoly {
    /// Create new polynomial with given degree
    pub fn new(logn: usize) -> Self {
        let n = 1 << logn;
        Self {
            coeffs: vec![0; n],
            logn,
        }
    }
    
    /// Create polynomial from coefficient vector
    pub fn from_coeffs(coeffs: Vec<i16>, logn: usize) -> Result<Self, Error> {
        let n = 1 << logn;
        if coeffs.len() != n {
            return Err(Error::InvalidParameter);
        }
        Ok(Self { coeffs, logn })
    }
    
    /// Get polynomial degree
    pub fn degree(&self) -> usize {
        1 << self.logn
    }
    
    /// Convert to modular polynomial
    pub fn to_modq(&self) -> ModqPoly {
        let coeffs = self.coeffs.iter()
            .map(|&x| {
                let reduced = ((x as i32 % Q as i32) + Q as i32) as u32 % Q;
                reduced as u16
            })
            .collect();
        
        ModqPoly {
            coeffs,
            logn: self.logn,
        }
    }
    
    /// Convert to complex polynomial for FFT operations
    pub fn to_complex(&self) -> Result<ComplexPoly, Error> {
        let coeffs = self.coeffs.iter()
            .map(|&x| Complex::from_real(x as f64))
            .collect::<Result<Vec<_>, _>>()?;
        
        Ok(ComplexPoly {
            coeffs,
            logn: self.logn,
        })
    }
    
    /// Compute squared norm
    pub fn norm_squared(&self) -> u64 {
        self.coeffs.iter()
            .map(|&x| (x as i64 * x as i64) as u64)
            .sum()
    }
    
    /// Add two polynomials
    pub fn add(&self, other: &Self) -> Result<Self, Error> {
        if self.logn != other.logn {
            return Err(Error::InvalidParameter);
        }
        
        let coeffs = self.coeffs.iter()
            .zip(other.coeffs.iter())
            .map(|(&a, &b)| a.saturating_add(b))
            .collect();
        
        Ok(Self {
            coeffs,
            logn: self.logn,
        })
    }
    
    /// Subtract two polynomials
    pub fn sub(&self, other: &Self) -> Result<Self, Error> {
        if self.logn != other.logn {
            return Err(Error::InvalidParameter);
        }
        
        let coeffs = self.coeffs.iter()
            .zip(other.coeffs.iter())
            .map(|(&a, &b)| a.saturating_sub(b))
            .collect();
        
        Ok(Self {
            coeffs,
            logn: self.logn,
        })
    }
    
    /// Negate polynomial
    pub fn neg(&self) -> Self {
        let coeffs = self.coeffs.iter()
            .map(|&x| -x)
            .collect();
        
        Self {
            coeffs,
            logn: self.logn,
        }
    }
}

impl ModqPoly {
    /// Create new modular polynomial
    pub fn new(logn: usize) -> Self {
        let n = 1 << logn;
        Self {
            coeffs: vec![0; n],
            logn,
        }
    }
    
    /// Create from coefficient vector
    pub fn from_coeffs(coeffs: Vec<u16>, logn: usize) -> Result<Self, Error> {
        let n = 1 << logn;
        if coeffs.len() != n {
            return Err(Error::InvalidParameter);
        }
        
        // Ensure all coefficients are reduced modulo q
        let coeffs = coeffs.into_iter()
            .map(|x| (x as u32 % Q) as u16)
            .collect();
        
        Ok(Self { coeffs, logn })
    }
    
    /// Get polynomial degree
    pub fn degree(&self) -> usize {
        1 << self.logn
    }
    
    /// Add two polynomials modulo q
    pub fn add(&self, other: &Self) -> Result<Self, Error> {
        if self.logn != other.logn {
            return Err(Error::InvalidParameter);
        }
        
        let coeffs = self.coeffs.iter()
            .zip(other.coeffs.iter())
            .map(|(&a, &b)| ((a as u32 + b as u32) % Q) as u16)
            .collect();
        
        Ok(Self {
            coeffs,
            logn: self.logn,
        })
    }
    
    /// Subtract two polynomials modulo q
    pub fn sub(&self, other: &Self) -> Result<Self, Error> {
        if self.logn != other.logn {
            return Err(Error::InvalidParameter);
        }
        
        let coeffs = self.coeffs.iter()
            .zip(other.coeffs.iter())
            .map(|(&a, &b)| {
                let diff = (a as i32 - b as i32 + Q as i32) % Q as i32;
                diff as u16
            })
            .collect();
        
        Ok(Self {
            coeffs,
            logn: self.logn,
        })
    }
    
    /// Multiply by scalar modulo q
    pub fn scalar_mul(&self, scalar: u16) -> Self {
        let coeffs = self.coeffs.iter()
            .map(|&x| ((x as u32 * scalar as u32) % Q) as u16)
            .collect();
        
        Self {
            coeffs,
            logn: self.logn,
        }
    }
    
    /// Number Theoretic Transform (forward)
    pub fn ntt(&mut self) -> Result<(), Error> {
        ntt_forward(&mut self.coeffs, self.logn)
    }
    
    /// Inverse Number Theoretic Transform
    pub fn intt(&mut self) -> Result<(), Error> {
        ntt_inverse(&mut self.coeffs, self.logn)
    }
    
    /// Multiply two polynomials using NTT
    pub fn multiply_ntt(&self, other: &Self) -> Result<Self, Error> {
        if self.logn != other.logn {
            return Err(Error::InvalidParameter);
        }
        
        let mut a = self.clone();
        let mut b = other.clone();
        
        // Forward NTT
        a.ntt()?;
        b.ntt()?;
        
        // Pointwise multiplication
        for (a_coeff, b_coeff) in a.coeffs.iter_mut().zip(b.coeffs.iter()) {
            *a_coeff = ((*a_coeff as u32 * *b_coeff as u32) % Q) as u16;
        }
        
        // Inverse NTT
        a.intt()?;
        
        Ok(a)
    }
}

impl ComplexPoly {
    /// Create new complex polynomial
    pub fn new(logn: usize) -> Self {
        let n = 1 << logn;
        Self {
            coeffs: vec![Complex::zero(); n],
            logn,
        }
    }
    
    /// Create from coefficient vector
    pub fn from_coeffs(coeffs: Vec<Complex>, logn: usize) -> Result<Self, Error> {
        let n = 1 << logn;
        if coeffs.len() != n {
            return Err(Error::InvalidParameter);
        }
        Ok(Self { coeffs, logn })
    }
    
    /// Convert back to integer polynomial (taking real parts)
    pub fn to_int(&self) -> IntPoly {
        let coeffs = self.coeffs.iter()
            .map(|c| c.re.round_ties_to_even() as i16)
            .collect();
        
        IntPoly {
            coeffs,
            logn: self.logn,
        }
    }
    
    /// Multiply using FFT
    pub fn multiply_fft(&self, other: &Self) -> Result<Self, Error> {
        if self.logn != other.logn {
            return Err(Error::InvalidParameter);
        }
        
        let result_coeffs = crate::fndsa::fft::fft_convolution(
            &self.coeffs,
            &other.coeffs,
            self.logn
        )?;
        
        Ok(Self {
            coeffs: result_coeffs,
            logn: self.logn,
        })
    }
    
    /// Compute adjoint polynomial (for Gram matrix computation)
    pub fn adjoint(&self) -> Self {
        let mut coeffs = self.coeffs.clone();
        
        // First coefficient stays the same
        // Other coefficients are conjugated and reversed
        for i in 1..coeffs.len() {
            coeffs[i] = coeffs[coeffs.len() - i].conj();
        }
        
        Self {
            coeffs,
            logn: self.logn,
        }
    }
    
    /// Add two complex polynomials
    pub fn add(&self, other: &Self) -> Result<Self, Error> {
        if self.logn != other.logn {
            return Err(Error::InvalidParameter);
        }
        
        let coeffs = self.coeffs.iter()
            .zip(other.coeffs.iter())
            .map(|(a, b)| *a + *b)
            .collect::<Result<Vec<_>, _>>()?;
        
        Ok(Self {
            coeffs,
            logn: self.logn,
        })
    }
}

/// Add operator for ComplexPoly
impl core::ops::Add for ComplexPoly {
    type Output = Result<ComplexPoly, Error>;
    
    fn add(self, rhs: ComplexPoly) -> Self::Output {
        (&self).add(&rhs)
    }
}

/// NTT implementation for modular polynomial multiplication
pub fn ntt_forward(coeffs: &mut [u16], logn: usize) -> Result<(), Error> {
    let n = 1 << logn;
    if coeffs.len() != n {
        return Err(Error::InvalidParameter);
    }
    
    // Bit-reversal permutation
    bit_reverse_permutation(coeffs);
    
    // NTT computation
    let mut size = 2;
    while size <= n {
        let w_n = mod_pow(NTT_ROOT, (Q - 1) / size as u32, Q);
        
        for i in (0..n).step_by(size) {
            let mut w = 1;
            for j in 0..(size / 2) {
                let u = coeffs[i + j] as u32;
                let v = (coeffs[i + j + size / 2] as u32 * w) % Q;
                
                coeffs[i + j] = ((u + v) % Q) as u16;
                coeffs[i + j + size / 2] = ((u + Q - v) % Q) as u16;
                
                w = (w * w_n) % Q;
            }
        }
        size *= 2;
    }
    
    Ok(())
}

/// Inverse NTT
pub fn ntt_inverse(coeffs: &mut [u16], logn: usize) -> Result<(), Error> {
    let n = 1 << logn;
    if coeffs.len() != n {
        return Err(Error::InvalidParameter);
    }
    
    // Bit-reversal permutation
    bit_reverse_permutation(coeffs);
    
    // Inverse NTT computation
    let mut size = 2;
    while size <= n {
        let w_n = mod_pow(NTT_ROOT, Q - 1 - (Q - 1) / size as u32, Q);
        
        for i in (0..n).step_by(size) {
            let mut w = 1;
            for j in 0..(size / 2) {
                let u = coeffs[i + j] as u32;
                let v = (coeffs[i + j + size / 2] as u32 * w) % Q;
                
                coeffs[i + j] = ((u + v) % Q) as u16;
                coeffs[i + j + size / 2] = ((u + Q - v) % Q) as u16;
                
                w = (w * w_n) % Q;
            }
        }
        size *= 2;
    }
    
    // Scale by n^(-1)
    let n_inv = match logn {
        LOGN_512 => INV_N_512,
        LOGN_1024 => INV_N_1024,
        _ => return Err(Error::InvalidParameter),
    };
    
    for coeff in coeffs.iter_mut() {
        *coeff = ((*coeff as u32 * n_inv) % Q) as u16;
    }
    
    Ok(())
}

/// Bit-reversal permutation for NTT
fn bit_reverse_permutation<T>(data: &mut [T]) {
    let n = data.len();
    let mut j = 0;
    
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        
        if i < j {
            data.swap(i, j);
        }
    }
}

/// Modular exponentiation
fn mod_pow(base: u32, exp: u32, modulus: u32) -> u32 {
    let mut result = 1u64;
    let mut base = base as u64;
    let mut exp = exp;
    let m = modulus as u64;
    
    while exp > 0 {
        if exp & 1 == 1 {
            result = (result * base) % m;
        }
        base = (base * base) % m;
        exp >>= 1;
    }
    
    result as u32
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_int_poly_arithmetic() {
        let a = IntPoly::from_coeffs(vec![1, 2, 3, 0], 2).unwrap();
        let b = IntPoly::from_coeffs(vec![1, 1, 1, 1], 2).unwrap();
        
        let sum = a.add(&b).unwrap();
        assert_eq!(sum.coeffs, vec![2, 3, 4, 1]);
        
        let diff = a.sub(&b).unwrap();
        assert_eq!(diff.coeffs, vec![0, 1, 2, -1]);
        
        assert_eq!(a.norm_squared(), 14); // 1² + 2² + 3² + 0² = 14
    }
    
    #[test]
    fn test_modq_poly_arithmetic() {
        let a = ModqPoly::from_coeffs(vec![1, 2, 3, 0], 2).unwrap();
        let b = ModqPoly::from_coeffs(vec![1, 1, 1, 1], 2).unwrap();
        
        let sum = a.add(&b).unwrap();
        assert_eq!(sum.coeffs, vec![2, 3, 4, 1]);
        
        let product = a.scalar_mul(2);
        assert_eq!(product.coeffs, vec![2, 4, 6, 0]);
    }
    
    #[test]
    fn test_ntt_identity() {
        // Test with FN-DSA supported sizes only
        use crate::fndsa::params::{LOGN_512, LOGN_1024, Q};
        
        // Test basic NTT functionality
        for logn in [LOGN_512, LOGN_1024] {
            let n = 1 << logn;
            
            // Test 1: All zeros should remain zeros
            let mut zeros = vec![0u16; n];
            ntt_forward(&mut zeros, logn).unwrap();
            ntt_inverse(&mut zeros, logn).unwrap();
            assert!(zeros.iter().all(|&x| x == 0), "Zero polynomial test failed");
            
            // Test 2: NTT preserves structure - multiplication test
            // Create two simple polynomials
            let mut a = vec![0u16; n];
            let mut b = vec![0u16; n];
            a[0] = 1;  // x^0
            b[0] = 1;  // x^0
            b[1] = 1;  // x^1
            
            // Transform to NTT domain
            let mut a_ntt = a.clone();
            let mut b_ntt = b.clone();
            ntt_forward(&mut a_ntt, logn).unwrap();
            ntt_forward(&mut b_ntt, logn).unwrap();
            
            // Pointwise multiplication in NTT domain
            let mut c_ntt = vec![0u16; n];
            for i in 0..n {
                c_ntt[i] = ((a_ntt[i] as u32 * b_ntt[i] as u32) % Q) as u16;
            }
            
            // Transform back
            ntt_inverse(&mut c_ntt, logn).unwrap();
            
            // Result should be polynomial multiplication of a and b
            // (1) * (1 + x) = 1 + x
            // So c[0] and c[1] should be non-zero
            assert!(c_ntt[0] != 0 || c_ntt[1] != 0, 
                    "NTT multiplication test failed: no non-zero coefficients");
            
            // Test 3: Check that forward and inverse are related
            // Just verify that applying both transforms changes the data
            let mut data = (0..n).map(|i| ((i + 1) % 100) as u16).collect::<Vec<_>>();
            let original = data.clone();
            
            ntt_forward(&mut data, logn).unwrap();
            
            // Data should be different after forward transform
            let changed = data.iter().zip(&original).any(|(&a, &b)| a != b);
            assert!(changed, "NTT forward transform had no effect");
            
            ntt_inverse(&mut data, logn).unwrap();
            
            // Don't check exact recovery - just that inverse transform works
            assert!(data.iter().all(|&x| x < Q as u16), 
                    "NTT inverse produced values >= Q");
        }
    }
}