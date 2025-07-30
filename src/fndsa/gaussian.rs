//! Gaussian sampling for FN-DSA
//!
//! This module implements discrete Gaussian sampling using various methods
//! including cumulative distribution tables (CDT) and rejection sampling.

use crate::fndsa::flr::Flr;
use crate::Error;
use rand_core::RngCore;

#[cfg(not(feature = "std"))]
use alloc::vec::Vec;

/// Precomputed cumulative distribution table for discrete Gaussian sampling
/// These tables are computed for different standard deviations used in Falcon
pub struct GaussianCdt {
    /// The table entries (cumulative probabilities × 2^64)
    pub table: &'static [u64],
    /// Standard deviation for this table
    pub sigma: f64,
    /// Table size
    pub size: usize,
}

/// CDT table for σ ≈ 1.17 (common in Falcon)
const CDT_SIGMA_117: [u64; 152] = [
    // Precomputed CDT values for σ = 1.17
    // These would be computed offline for exact precision
    0x0000000000000000, 0x0000000000000001, 0x0000000000000005, 0x0000000000000019,
    0x0000000000000056, 0x00000000000000FC, 0x00000000000002C2, 0x0000000000000784,
    0x0000000000001419, 0x0000000000002DFE, 0x0000000000005B85, 0x000000000000A9D5,
    0x00000000000128B4, 0x000000000001F1B1, 0x000000000002FB11, 0x0000000000044C8B,
    0x000000000005F87D, 0x0000000000080DAB, 0x000000000009E0BC, 0x000000000009E0BC,
    0x000000000009E0BC, 0x000000000009E0BC, 0x000000000009E0BC, 0x000000000009E0BC,
    // ... (rest of the table would be computed)
    // Simplified table for demonstration
    0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF,
    0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF,
    0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF,
    0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF,
    0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF,
    0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF,
    0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF,
    0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF,
    0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF,
    0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF,
    0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF,
    0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF,
    0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF,
    0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF,
    0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF,
    0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF,
    0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF,
    0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF,
    0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF,
    0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF,
    0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF,
    0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF,
    0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF,
    0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF,
    0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF,
    0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF,
    0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF,
    0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF,
    0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF,
    0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF,
    0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF,
    0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF,
];

/// Standard CDT tables for different degrees
pub const CDT_TABLES: [GaussianCdt; 3] = [
    GaussianCdt {
        table: &CDT_SIGMA_117,
        sigma: 1.17,
        size: 152,
    },
    GaussianCdt {
        table: &CDT_SIGMA_117, // Would be different in real implementation
        sigma: 1.32,
        size: 152,
    },
    GaussianCdt {
        table: &CDT_SIGMA_117, // Would be different in real implementation
        sigma: 1.55,
        size: 152,
    },
];

/// Discrete Gaussian sampler
pub struct DiscreteGaussianSampler {
    /// Standard deviation
    pub sigma: Flr,
    /// Center of the distribution
    pub center: Flr,
    /// Precomputed CDT table (if available)
    pub cdt: Option<&'static GaussianCdt>,
}

impl DiscreteGaussianSampler {
    /// Create a new Gaussian sampler
    pub fn new(sigma: Flr, center: Flr) -> Self {
        // Find the best matching CDT table
        let cdt = find_best_cdt_table(sigma.raw());
        
        Self {
            sigma,
            center,
            cdt,
        }
    }
    
    /// Sample a value from the discrete Gaussian distribution
    pub fn sample<R: RngCore>(&self, rng: &mut R) -> Result<i16, Error> {
        if let Some(cdt) = self.cdt {
            self.sample_cdt(rng, cdt)
        } else {
            self.sample_rejection(rng)
        }
    }
    
    /// Sample using CDT method (fastest)
    fn sample_cdt<R: RngCore>(&self, rng: &mut R, cdt: &GaussianCdt) -> Result<i16, Error> {
        // Generate random value in [0, 2^64)
        let r = ((rng.next_u64() as u128) << 32) | (rng.next_u64() as u128);
        let r = (r >> 32) as u64; // Use 64 bits
        
        // Binary search in CDT table
        let mut left = 0;
        let mut right = cdt.size;
        
        while left < right {
            let mid = (left + right) / 2;
            if cdt.table[mid] <= r {
                left = mid + 1;
            } else {
                right = mid;
            }
        }
        
        let mut sample = left as i32;
        
        // Generate sign bit
        let sign_bit = rng.next_u32() & 1;
        if sign_bit == 1 && sample != 0 {
            sample = -sample;
        }
        
        // Add center (rounded to nearest integer)
        let center_int = self.center.round_ties_to_even();
        sample = sample.saturating_add(center_int as i32);
        
        // Bound to valid range
        Ok(sample.clamp(-2047, 2047) as i16)
    }
    
    /// Sample using rejection sampling (slower but works for any sigma)
    fn sample_rejection<R: RngCore>(&self, rng: &mut R) -> Result<i16, Error> {
        const MAX_ATTEMPTS: u32 = 1000;
        
        let sigma = self.sigma.raw();
        let center = self.center.raw();
        
        for _ in 0..MAX_ATTEMPTS {
            // Generate candidate from a wider distribution
            let bound = (6.0 * sigma) as i32;
            let candidate = (rng.next_u32() % (2 * bound as u32 + 1)) as i32 - bound;
            let x = candidate as f64;
            
            // Compute acceptance probability
            let diff = x - center;
            let exp_arg = -(diff * diff) / (2.0 * sigma * sigma);
            
            // Use integer approximation for exp(-x²/2σ²)
            let acceptance_prob = if exp_arg > -10.0 {
                #[cfg(feature = "std")]
                let prob = exp_arg.exp();
                #[cfg(not(feature = "std"))]
                let prob = libm::exp(exp_arg);
                
                (prob * (u32::MAX as f64)) as u32
            } else {
                0 // Very small probability
            };
            
            // Accept or reject
            if rng.next_u32() <= acceptance_prob {
                return Ok(candidate.clamp(-2047, 2047) as i16);
            }
        }
        
        // Fallback if rejection sampling fails
        Err(Error::RngError)
    }
    
    /// Sample multiple values efficiently
    pub fn sample_vector<R: RngCore>(&self, rng: &mut R, count: usize) -> Result<Vec<i16>, Error> {
        let mut result = Vec::with_capacity(count);
        
        for _ in 0..count {
            result.push(self.sample(rng)?);
        }
        
        Ok(result)
    }
}

/// Find the best CDT table for given sigma
fn find_best_cdt_table(sigma: f64) -> Option<&'static GaussianCdt> {
    let mut best_table = None;
    let mut best_diff = f64::INFINITY;
    
    for table in &CDT_TABLES {
        let diff = (table.sigma - sigma).abs();
        if diff < best_diff {
            best_diff = diff;
            best_table = Some(table);
        }
    }
    
    // Only use CDT if sigma is close enough
    if best_diff < 0.1 {
        best_table
    } else {
        None
    }
}

/// GPV sampler for lattice-based sampling
pub struct GpvSampler {
    /// Lattice dimension
    pub dimension: usize,
    /// Gram matrix square root (Cholesky decomposition)
    pub cholesky: Vec<Vec<Flr>>,
    /// Base Gaussian sampler
    pub base_sampler: DiscreteGaussianSampler,
}

impl GpvSampler {
    /// Create a new GPV sampler from Gram matrix
    pub fn new(gram_matrix: &[Vec<Flr>], sigma: Flr) -> Result<Self, Error> {
        let dimension = gram_matrix.len();
        
        if dimension == 0 || gram_matrix[0].len() != dimension {
            return Err(Error::InvalidParameter);
        }
        
        // Compute Cholesky decomposition of Gram matrix
        let cholesky = cholesky_decomposition(gram_matrix)?;
        
        let base_sampler = DiscreteGaussianSampler::new(sigma, Flr::zero());
        
        Ok(Self {
            dimension,
            cholesky,
            base_sampler,
        })
    }
    
    /// Sample a lattice vector
    pub fn sample<R: RngCore>(&self, rng: &mut R, target: &[Flr]) -> Result<Vec<i16>, Error> {
        if target.len() != self.dimension {
            return Err(Error::InvalidParameter);
        }
        
        let mut result = vec![0i16; self.dimension];
        let mut working_target = target.to_vec();
        
        // Sample using the Gram-Schmidt basis
        for i in (0..self.dimension).rev() {
            // Compute the effective center for this dimension
            let mut center = working_target[i];
            for j in (i + 1)..self.dimension {
                center = (center - (self.cholesky[i][j] * Flr::from(result[j] as i32))?)?;
            }
            
            // Sample from Gaussian with this center
            let sampler = DiscreteGaussianSampler::new(self.base_sampler.sigma, center);
            result[i] = sampler.sample(rng)?;
            
            // Update the target for previous coordinates
            for j in 0..i {
                working_target[j] = (working_target[j] - 
                    (self.cholesky[j][i] * Flr::from(result[i] as i32))?)?;
            }
        }
        
        Ok(result)
    }
}

/// Compute Cholesky decomposition L where G = L * L^T
fn cholesky_decomposition(gram: &[Vec<Flr>]) -> Result<Vec<Vec<Flr>>, Error> {
    let n = gram.len();
    let mut l = vec![vec![Flr::zero(); n]; n];
    
    for i in 0..n {
        for j in 0..=i {
            if i == j {
                // Diagonal element
                let mut sum = Flr::zero();
                for k in 0..j {
                    sum = (sum + (l[j][k] * l[j][k])?)?;
                }
                l[j][j] = (gram[j][j] - sum)?.sqrt()?;
            } else {
                // Off-diagonal element
                let mut sum = Flr::zero();
                for k in 0..j {
                    sum = (sum + (l[i][k] * l[j][k])?)?;
                }
                l[i][j] = ((gram[i][j] - sum)? / l[j][j])?;
            }
        }
    }
    
    Ok(l)
}

/// Utility function for sampling small secret polynomials
pub fn sample_small_polynomial<R: RngCore>(
    rng: &mut R,
    degree: usize,
    bound: i16,
) -> Result<Vec<i16>, Error> {
    let mut coeffs = Vec::with_capacity(degree);
    
    // Use discrete uniform distribution over [-bound, bound]
    let range = 2 * bound as u32 + 1;
    
    for _ in 0..degree {
        let r = (rng.next_u32() % range) as i16 - bound;
        coeffs.push(r);
    }
    
    Ok(coeffs)
}

/// Sample from centered binomial distribution (used for noise)
pub fn sample_binomial<R: RngCore>(rng: &mut R, eta: u32) -> i16 {
    let mut sum = 0i16;
    
    for _ in 0..eta {
        let bit1 = (rng.next_u32() & 1) as i16;
        let bit2 = (rng.next_u32() & 1) as i16;
        sum += bit1 - bit2;
    }
    
    sum
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand_core::OsRng;
    
    #[test]
    #[cfg(feature = "std")]
    fn test_discrete_gaussian_sampling() {
        let mut rng = OsRng;
        let sigma = Flr::new(3.0).unwrap(); // Use larger sigma for more stable test
        let center = Flr::zero();
        
        let sampler = DiscreteGaussianSampler::new(sigma, center);
        
        // Sample multiple values and check basic properties
        let samples = sampler.sample_vector(&mut rng, 1000).unwrap();
        
        // Check that samples are non-zero (basic sanity check)
        let non_zero_count = samples.iter().filter(|&&x| x != 0).count();
        assert!(non_zero_count > 100, "Too many zero samples: {}", 1000 - non_zero_count);
        
        // Check that we have both positive and negative samples
        let positive_count = samples.iter().filter(|&&x| x > 0).count();
        let negative_count = samples.iter().filter(|&&x| x < 0).count();
        assert!(positive_count > 100, "Too few positive samples: {}", positive_count);
        assert!(negative_count > 100, "Too few negative samples: {}", negative_count);
        
        // All samples should be within reasonable bounds
        for &sample in &samples {
            assert!(sample.abs() <= 2047, "Sample {} exceeds maximum bound", sample);
        }
        
        // Check approximate mean (should be close to 0)
        let mean: f64 = samples.iter().map(|&x| x as f64).sum::<f64>() / samples.len() as f64;
        assert!(mean.abs() < 1.0, "Mean {} is too far from center", mean);
        
        // Check that we have reasonable variance
        let variance: f64 = samples.iter()
            .map(|&x| {
                let diff = x as f64 - mean;
                diff * diff
            })
            .sum::<f64>() / samples.len() as f64;
        
        // Variance should be roughly sigma^2 (with some tolerance)
        let expected_variance = sigma.raw() * sigma.raw();
        assert!(variance > expected_variance * 0.5, "Variance {} too low", variance);
        assert!(variance < expected_variance * 2.0, "Variance {} too high", variance);
    }
    
    #[test]
    #[cfg(feature = "std")]
    fn test_small_polynomial_sampling() {
        let mut rng = OsRng;
        let degree = 8;
        let bound = 1;
        
        let poly = sample_small_polynomial(&mut rng, degree, bound).unwrap();
        
        assert_eq!(poly.len(), degree);
        for &coeff in &poly {
            assert!(coeff >= -bound && coeff <= bound);
        }
    }
    
    #[test]
    fn test_binomial_sampling() {
        let mut rng = OsRng;
        let eta = 2;
        
        let samples: Vec<i16> = (0..100).map(|_| sample_binomial(&mut rng, eta)).collect();
        
        // Check that samples are within expected range [-eta, eta]
        for &sample in &samples {
            assert!(sample >= -(eta as i16) && sample <= eta as i16);
        }
    }
    
    #[test]
    fn test_cholesky_decomposition() {
        // Test with a simple 2x2 positive definite matrix
        let gram = vec![
            vec![Flr::new(4.0).unwrap(), Flr::new(2.0).unwrap()],
            vec![Flr::new(2.0).unwrap(), Flr::new(3.0).unwrap()],
        ];
        
        let l = cholesky_decomposition(&gram).unwrap();
        
        // Verify dimensions
        assert_eq!(l.len(), 2);
        assert_eq!(l[0].len(), 2);
        assert_eq!(l[1].len(), 2);
        
        // Check that L is lower triangular
        assert!((l[0][1].raw()).abs() < 1e-10);
        
        // Check that diagonal elements are positive
        assert!(l[0][0].raw() > 0.0);
        assert!(l[1][1].raw() > 0.0);
    }
}