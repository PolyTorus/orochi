//! IEEE-754 compliant floating-point operations for FN-DSA
//!
//! This module implements the FLR (Floating-point with Limited Range) system
//! used in Falcon for secure and deterministic floating-point arithmetic.

#[cfg(not(feature = "std"))]
use libm;

use crate::Error;

/// FLR floating-point number with controlled precision
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Flr(f64);

/// IEEE-754 binary64 format constants
pub mod ieee754 {
    /// Sign bit mask
    pub const SIGN_MASK: u64 = 0x8000_0000_0000_0000;
    
    /// Exponent mask  
    pub const EXPONENT_MASK: u64 = 0x7FF0_0000_0000_0000;
    
    /// Mantissa mask
    pub const MANTISSA_MASK: u64 = 0x000F_FFFF_FFFF_FFFF;
    
    /// Exponent bias
    pub const EXPONENT_BIAS: i32 = 1023;
    
    /// Mantissa bits
    pub const MANTISSA_BITS: u32 = 52;
    
    /// Exponent bits
    pub const EXPONENT_BITS: u32 = 11;
    
    /// Minimum normal exponent (after bias)
    pub const MIN_NORMAL_EXP: i32 = -1022;
    
    /// Maximum normal exponent (after bias)
    pub const MAX_NORMAL_EXP: i32 = 1023;
    
    /// Falcon-specific exponent range for controlled precision
    pub const FALCON_MIN_EXP: i32 = -478;  // 547 - 1023 - 2
    pub const FALCON_MAX_EXP: i32 = 79;    // 1102 - 1023
}

impl Flr {
    /// Create a new FLR number from f64, checking for validity
    pub fn new(value: f64) -> Result<Self, Error> {
        if !Self::is_valid(value) {
            return Err(Error::InvalidParameter);
        }
        Ok(Flr(value))
    }
    
    /// Create an FLR number without validation (unsafe)
    pub fn new_unchecked(value: f64) -> Self {
        Flr(value)
    }
    
    /// Get the raw f64 value
    pub fn raw(self) -> f64 {
        self.0
    }
    
    /// Check if a floating-point value is valid for FLR
    pub fn is_valid(value: f64) -> bool {
        // Check for special values
        if value.is_nan() || value.is_infinite() {
            return false;
        }
        
        // Zero is always valid
        if value == 0.0 {
            return true;
        }
        
        let bits = value.to_bits();
        let exponent = Self::extract_exponent(bits);
        
        // Check exponent range for Falcon
        exponent >= ieee754::FALCON_MIN_EXP && exponent <= ieee754::FALCON_MAX_EXP
    }
    
    /// Extract the biased exponent from IEEE-754 bits
    fn extract_exponent(bits: u64) -> i32 {
        let exp_bits = ((bits & ieee754::EXPONENT_MASK) >> ieee754::MANTISSA_BITS) as i32;
        if exp_bits == 0 {
            // Subnormal numbers (exponent is effectively -1022)
            ieee754::MIN_NORMAL_EXP
        } else {
            exp_bits - ieee754::EXPONENT_BIAS
        }
    }
    
    /// Create FLR zero
    pub const fn zero() -> Self {
        Flr(0.0)
    }
    
    /// Create FLR one
    pub const fn one() -> Self {
        Flr(1.0)
    }
    
    /// Negate the FLR number
    pub fn neg(self) -> Self {
        Flr(-self.0)
    }
    
    /// Absolute value
    pub fn abs(self) -> Self {
        Flr(self.0.abs())
    }
    
    /// Square root with proper rounding
    pub fn sqrt(self) -> Result<Self, Error> {
        if self.0 < 0.0 {
            return Err(Error::InvalidParameter);
        }
        
        #[cfg(feature = "std")]
        let result = self.0.sqrt();
        
        #[cfg(not(feature = "std"))]
        let result = libm::sqrt(self.0);
        
        Self::new(result)
    }
    
    /// Natural logarithm
    pub fn ln(self) -> Result<Self, Error> {
        if self.0 <= 0.0 {
            return Err(Error::InvalidParameter);
        }
        
        #[cfg(feature = "std")]
        let result = self.0.ln();
        
        #[cfg(not(feature = "std"))]
        let result = libm::log(self.0);
        
        Self::new(result)
    }
    
    /// Exponential function
    pub fn exp(self) -> Result<Self, Error> {
        #[cfg(feature = "std")]
        let result = self.0.exp();
        
        #[cfg(not(feature = "std"))]
        let result = libm::exp(self.0);
        
        Self::new(result)
    }
    
    /// Sine function
    pub fn sin(self) -> Result<Self, Error> {
        #[cfg(feature = "std")]
        let result = self.0.sin();
        
        #[cfg(not(feature = "std"))]
        let result = libm::sin(self.0);
        
        Self::new(result)
    }
    
    /// Cosine function
    pub fn cos(self) -> Result<Self, Error> {
        #[cfg(feature = "std")]
        let result = self.0.cos();
        
        #[cfg(not(feature = "std"))]
        let result = libm::cos(self.0);
        
        Self::new(result)
    }
    
    /// Round to nearest integer with ties to even
    pub fn round_ties_to_even(self) -> i64 {
        // IEEE-754 compliant rounding
        let value = self.0;
        let truncated = value.trunc();
        let fraction = value - truncated;
        
        if fraction.abs() < 0.5 {
            truncated as i64
        } else if fraction.abs() > 0.5 {
            if value > 0.0 {
                (truncated + 1.0) as i64
            } else {
                (truncated - 1.0) as i64
            }
        } else {
            // Tie case - round to even
            let truncated_int = truncated as i64;
            if truncated_int % 2 == 0 {
                truncated_int
            } else if value > 0.0 {
                truncated_int + 1
            } else {
                truncated_int - 1
            }
        }
    }
}

/// Arithmetic operations for FLR
impl core::ops::Add for Flr {
    type Output = Result<Flr, Error>;
    
    fn add(self, rhs: Flr) -> Self::Output {
        let result = self.0 + rhs.0;
        Flr::new(result)
    }
}

impl core::ops::Sub for Flr {
    type Output = Result<Flr, Error>;
    
    fn sub(self, rhs: Flr) -> Self::Output {
        let result = self.0 - rhs.0;
        Flr::new(result)
    }
}

impl core::ops::Mul for Flr {
    type Output = Result<Flr, Error>;
    
    fn mul(self, rhs: Flr) -> Self::Output {
        let result = self.0 * rhs.0;
        Flr::new(result)
    }
}

impl core::ops::Div for Flr {
    type Output = Result<Flr, Error>;
    
    fn div(self, rhs: Flr) -> Self::Output {
        if rhs.0 == 0.0 {
            return Err(Error::InvalidParameter);
        }
        let result = self.0 / rhs.0;
        Flr::new(result)
    }
}

/// Comparison operations
impl PartialOrd for Flr {
    fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
        self.0.partial_cmp(&other.0)
    }
}

/// From conversions with validation
impl TryFrom<f64> for Flr {
    type Error = Error;
    
    fn try_from(value: f64) -> Result<Self, Self::Error> {
        Flr::new(value)
    }
}

impl TryFrom<f32> for Flr {
    type Error = Error;
    
    fn try_from(value: f32) -> Result<Self, Self::Error> {
        Flr::new(value as f64)
    }
}

impl From<i32> for Flr {
    fn from(value: i32) -> Self {
        Flr::new_unchecked(value as f64)
    }
}

impl From<Flr> for f64 {
    fn from(flr: Flr) -> Self {
        flr.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_flr_creation() {
        assert!(Flr::new(1.0).is_ok());
        assert!(Flr::new(0.0).is_ok());
        assert!(Flr::new(-1.0).is_ok());
        
        assert!(Flr::new(f64::NAN).is_err());
        assert!(Flr::new(f64::INFINITY).is_err());
        assert!(Flr::new(f64::NEG_INFINITY).is_err());
    }
    
    #[test] 
    fn test_arithmetic() {
        let a = Flr::new(2.0).unwrap();
        let b = Flr::new(3.0).unwrap();
        
        assert_eq!((a + b).unwrap().raw(), 5.0);
        assert_eq!((a * b).unwrap().raw(), 6.0);
        assert_eq!((b - a).unwrap().raw(), 1.0);
        assert_eq!((b / a).unwrap().raw(), 1.5);
    }
    
    #[test]
    fn test_round_ties_to_even() {
        assert_eq!(Flr::new(1.5).unwrap().round_ties_to_even(), 2);
        assert_eq!(Flr::new(2.5).unwrap().round_ties_to_even(), 2);
        assert_eq!(Flr::new(3.5).unwrap().round_ties_to_even(), 4);
        assert_eq!(Flr::new(-1.5).unwrap().round_ties_to_even(), -2);
        assert_eq!(Flr::new(-2.5).unwrap().round_ties_to_even(), -2);
    }
}