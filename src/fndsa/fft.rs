//! Complex FFT implementation for FN-DSA polynomial operations
//!
//! This module implements the Fast Fourier Transform over complex numbers
//! required for Falcon's Fast Fourier Orthogonalization algorithm.

use crate::fndsa::flr::Flr;
use crate::Error;

#[cfg(not(feature = "std"))]
use alloc::vec::Vec;

/// Complex number using FLR floating-point arithmetic
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Complex {
    pub re: Flr,
    pub im: Flr,
}

impl Complex {
    /// Create a new complex number
    pub fn new(re: f64, im: f64) -> Result<Self, Error> {
        Ok(Complex {
            re: Flr::new(re)?,
            im: Flr::new(im)?,
        })
    }
    
    /// Create complex zero
    pub const fn zero() -> Self {
        Complex {
            re: Flr::zero(),
            im: Flr::zero(),
        }
    }
    
    /// Create complex one
    pub const fn one() -> Self {
        Complex {
            re: Flr::one(),
            im: Flr::zero(),
        }
    }
    
    /// Create complex from real part only
    pub fn from_real(re: f64) -> Result<Self, Error> {
        Ok(Complex {
            re: Flr::new(re)?,
            im: Flr::zero(),
        })
    }
    
    /// Complex conjugate
    pub fn conj(self) -> Self {
        Complex {
            re: self.re,
            im: self.im.neg(),
        }
    }
    
    /// Magnitude squared |z|²
    pub fn norm_squared(self) -> Result<Flr, Error> {
        let re_sq = (self.re * self.re)?;
        let im_sq = (self.im * self.im)?;
        re_sq + im_sq
    }
    
    /// Magnitude |z|
    pub fn magnitude(self) -> Result<Flr, Error> {
        let norm_sq = self.norm_squared()?;
        norm_sq.sqrt()
    }
    
    /// Argument (phase) of complex number
    pub fn arg(self) -> Result<Flr, Error> {
        #[cfg(feature = "std")]
        let result = self.im.raw().atan2(self.re.raw());
        
        #[cfg(not(feature = "std"))]
        let result = libm::atan2(self.im.raw(), self.re.raw());
        
        Flr::new(result)
    }
}

/// Addition
impl core::ops::Add for Complex {
    type Output = Result<Complex, Error>;
    
    fn add(self, rhs: Complex) -> Self::Output {
        Ok(Complex {
            re: (self.re + rhs.re)?,
            im: (self.im + rhs.im)?,
        })
    }
}

/// Subtraction
impl core::ops::Sub for Complex {
    type Output = Result<Complex, Error>;
    
    fn sub(self, rhs: Complex) -> Self::Output {
        Ok(Complex {
            re: (self.re - rhs.re)?,
            im: (self.im - rhs.im)?,
        })
    }
}

/// Multiplication
impl core::ops::Mul for Complex {
    type Output = Result<Complex, Error>;
    
    fn mul(self, rhs: Complex) -> Self::Output {
        // (a + bi)(c + di) = (ac - bd) + (ad + bc)i
        let ac = (self.re * rhs.re)?;
        let bd = (self.im * rhs.im)?;
        let ad = (self.re * rhs.im)?;
        let bc = (self.im * rhs.re)?;
        
        Ok(Complex {
            re: (ac - bd)?,
            im: (ad + bc)?,
        })
    }
}

/// Division
impl core::ops::Div for Complex {
    type Output = Result<Complex, Error>;
    
    fn div(self, rhs: Complex) -> Self::Output {
        // (a + bi)/(c + di) = [(a + bi)(c - di)] / (c² + d²)
        let denominator = rhs.norm_squared()?;
        let numerator = (self * rhs.conj())?;
        
        Ok(Complex {
            re: (numerator.re / denominator)?,
            im: (numerator.im / denominator)?,
        })
    }
}

/// Complex exponential: e^(i*θ) = cos(θ) + i*sin(θ)
pub fn complex_exp(theta: Flr) -> Result<Complex, Error> {
    let cos_theta = theta.cos()?;
    let sin_theta = theta.sin()?;
    
    Ok(Complex {
        re: cos_theta,
        im: sin_theta,
    })
}

/// FFT implementation using Cooley-Tukey algorithm
pub struct ComplexFFT {
    n: usize,
    logn: usize,
    omega: Vec<Complex>,      // Roots of unity
    inv_omega: Vec<Complex>,  // Inverse roots of unity
}

impl ComplexFFT {
    /// Create new FFT instance for given size
    pub fn new(logn: usize) -> Result<Self, Error> {
        let n = 1 << logn;
        let mut omega = Vec::with_capacity(n);
        let mut inv_omega = Vec::with_capacity(n);
        
        // Compute roots of unity: ω_n^k = e^(2πik/n)
        let pi = Flr::new(core::f64::consts::PI)?;
        let two = Flr::new(2.0)?;
        let n_flr = Flr::new(n as f64)?;
        let base_angle = ((two * pi)? / n_flr)?;
        
        for k in 0..n {
            let angle = (base_angle * Flr::new(k as f64)?)?;
            let w = complex_exp(angle)?;
            let w_inv = complex_exp(angle.neg())?;
            
            omega.push(w);
            inv_omega.push(w_inv);
        }
        
        Ok(ComplexFFT {
            n,
            logn,
            omega,
            inv_omega,
        })
    }
    
    /// Forward FFT (time to frequency domain)
    pub fn forward(&self, data: &mut [Complex]) -> Result<(), Error> {
        if data.len() != self.n {
            return Err(Error::InvalidParameter);
        }
        
        self.fft_internal(data, &self.omega)
    }
    
    /// Inverse FFT (frequency to time domain)
    pub fn inverse(&self, data: &mut [Complex]) -> Result<(), Error> {
        if data.len() != self.n {
            return Err(Error::InvalidParameter);
        }
        
        // Inverse FFT using conjugated roots
        self.fft_internal(data, &self.inv_omega)?;
        
        // Scale by 1/n
        let n_inv = Flr::new(1.0 / (self.n as f64))?;
        for x in data.iter_mut() {
            x.re = (x.re * n_inv)?;
            x.im = (x.im * n_inv)?;
        }
        
        Ok(())
    }
    
    /// Internal FFT implementation using Cooley-Tukey algorithm
    fn fft_internal(&self, data: &mut [Complex], roots: &[Complex]) -> Result<(), Error> {
        // Bit-reversal permutation
        self.bit_reverse_permutation(data);
        
        // FFT butterfly operations
        let mut size = 2;
        let mut step = self.n / 2;
        
        for _ in 0..self.logn {
            for i in (0..self.n).step_by(size) {
                for j in 0..(size / 2) {
                    let u = data[i + j];
                    let v = (data[i + j + size / 2] * roots[step * j])?;
                    
                    data[i + j] = (u + v)?;
                    data[i + j + size / 2] = (u - v)?;
                }
            }
            size *= 2;
            step /= 2;
        }
        
        Ok(())
    }
    
    /// Bit-reversal permutation for FFT
    fn bit_reverse_permutation(&self, data: &mut [Complex]) {
        let mut j = 0;
        for i in 1..self.n {
            let mut bit = self.n >> 1;
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
}

/// Convolution using FFT
pub fn fft_convolution(
    a: &[Complex], 
    b: &[Complex], 
    logn: usize
) -> Result<Vec<Complex>, Error> {
    let n = 1 << logn;
    if a.len() != n || b.len() != n {
        return Err(Error::InvalidParameter);
    }
    
    let fft = ComplexFFT::new(logn)?;
    
    let mut a_fft = a.to_vec();
    let mut b_fft = b.to_vec();
    
    // Forward FFT
    fft.forward(&mut a_fft)?;
    fft.forward(&mut b_fft)?;
    
    // Pointwise multiplication
    for (a_val, b_val) in a_fft.iter_mut().zip(b_fft.iter()) {
        *a_val = (*a_val * *b_val)?;
    }
    
    // Inverse FFT
    fft.inverse(&mut a_fft)?;
    
    Ok(a_fft)
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_complex_arithmetic() {
        let a = Complex::new(1.0, 2.0).unwrap();
        let b = Complex::new(3.0, 4.0).unwrap();
        
        let sum = (a + b).unwrap();
        assert!((sum.re.raw() - 4.0).abs() < 1e-10);
        assert!((sum.im.raw() - 6.0).abs() < 1e-10);
        
        let product = (a * b).unwrap();
        assert!((product.re.raw() - (-5.0)).abs() < 1e-10); // 1*3 - 2*4 = -5
        assert!((product.im.raw() - 10.0).abs() < 1e-10);   // 1*4 + 2*3 = 10
    }
    
    #[test]
    fn test_fft_identity() {
        let logn = 3;
        let n = 1 << logn;
        let fft = ComplexFFT::new(logn).unwrap();
        
        let mut data: Vec<Complex> = (0..n)
            .map(|i| Complex::from_real(i as f64).unwrap())
            .collect();
        let original = data.clone();
        
        // Forward then inverse should give identity
        fft.forward(&mut data).unwrap();
        fft.inverse(&mut data).unwrap();
        
        for (a, b) in data.iter().zip(original.iter()) {
            assert!((a.re.raw() - b.re.raw()).abs() < 1e-10);
            assert!((a.im.raw() - b.im.raw()).abs() < 1e-10);
        }
    }
}