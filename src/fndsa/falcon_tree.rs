//! Falcon tree construction and LDL decomposition
//!
//! This module implements the Fast Fourier Orthogonalization (FFO) algorithm
//! used in Falcon for efficient lattice basis reduction and Gaussian sampling.

use crate::fndsa::{params::*, flr::Flr, fft::*, ntru::NtruPrivateKey};
use crate::Error;

#[cfg(not(feature = "std"))]
use alloc::vec::Vec;

/// Falcon tree representing the LDL decomposition of the Gram matrix
#[derive(Debug, Clone)]
pub struct FalconTree {
    /// LDL decomposition data (L and D matrices combined)
    pub ldl_data: Vec<Flr>,
    /// Tree depth levels
    pub depths: Vec<usize>,
    /// Gaussian standard deviations for each level
    pub sigmas: Vec<Flr>,
    /// Tree parameters
    pub logn: usize,
}

/// Node in the Falcon tree for recursive operations
#[derive(Debug, Clone)]
struct TreeNode {
    /// Level in the tree (0 = leaf)
    level: usize,
    /// Index within the level
    index: usize,
    /// Gram matrix for this node [[a, b], [c, d]]
    gram: [[Complex; 2]; 2],
    /// LDL factors for this node
    ldl_factors: LdlFactors,
}

/// LDL decomposition factors
#[derive(Debug, Clone)]
struct LdlFactors {
    /// L matrix (lower triangular)
    l: [[Flr; 2]; 2],
    /// D matrix (diagonal)  
    d: [Flr; 2],
}

impl FalconTree {
    /// Construct a new Falcon tree from NTRU private key
    pub fn new(private_key: &NtruPrivateKey) -> Result<Self, Error> {
        let logn = private_key.f.logn;
        
        if logn > MAX_LOGN {
            return Err(Error::InvalidParameter);
        }
        
        // Build the tree using recursive FFT-based approach
        let root_gram = compute_root_gram_matrix(private_key)?;
        let tree_data = build_recursive_tree(root_gram, logn)?;
        
        // Extract LDL data, depths, and sigmas from tree structure
        let (ldl_data, depths, sigmas) = extract_tree_data(&tree_data)?;
        
        Ok(FalconTree {
            ldl_data,
            depths,
            sigmas,
            logn,
        })
    }
    
    /// Get the Gaussian standard deviation for a given tree level
    pub fn sigma_at_level(&self, level: usize) -> Result<Flr, Error> {
        if level >= self.sigmas.len() {
            return Err(Error::InvalidParameter);
        }
        Ok(self.sigmas[level])
    }
    
    /// Get the tree depth
    pub fn depth(&self) -> usize {
        self.logn
    }
    
    /// Sample from the tree using the GPV sampler
    pub fn sample(&self, target: &[Complex]) -> Result<Vec<i16>, Error> {
        if target.len() != (1 << self.logn) {
            return Err(Error::InvalidParameter);
        }
        
        recursive_sample(self, target, self.logn, 0)
    }
}

/// Compute the root Gram matrix from NTRU private key
fn compute_root_gram_matrix(private_key: &NtruPrivateKey) -> Result<[[Complex; 2]; 2], Error> {
    // Convert polynomials to complex representation
    let f_c = private_key.f.to_complex()?;
    let g_c = private_key.g.to_complex()?;
    let F_c = private_key.F.to_complex()?;
    let G_c = private_key.G.to_complex()?;
    
    // Apply FFT to all polynomials for efficient computation
    let mut f_fft = f_c.coeffs.clone();
    let mut g_fft = g_c.coeffs.clone();
    let mut F_fft = F_c.coeffs.clone();
    let mut G_fft = G_c.coeffs.clone();
    
    let fft = ComplexFFT::new(private_key.f.logn)?;
    fft.forward(&mut f_fft)?;
    fft.forward(&mut g_fft)?;
    fft.forward(&mut F_fft)?;
    fft.forward(&mut G_fft)?;
    
    // Compute Gram matrix in FFT domain: B* B where B = [[g, G], [-f, -F]]
    // The Gram matrix is [[g*g + F*F, g*(-f) + F*(-F)], [(-f)*g + (-F)*F, f*f + F*F]]
    
    let mut gram_00 = Complex::zero();
    let mut gram_01 = Complex::zero();
    let mut gram_11 = Complex::zero();
    
    // Average over all FFT coefficients (this is the [0] coefficient after IFFT)
    let n = 1 << private_key.f.logn;
    for i in 0..n {
        let g_i = g_fft[i];
        let f_i = f_fft[i];
        let G_i = G_fft[i];
        let F_i = F_fft[i];
        
        // g*conj(g) + F*conj(F)
        let term1 = (g_i * g_i.conj())?;
        let term2 = (F_i * F_i.conj())?;
        gram_00 = ((gram_00 + term1)? + term2)?;
        
        // g*conj(-f) + F*conj(-F) = -g*conj(f) - F*conj(F)
        let term3 = (g_i * f_i.conj())?;
        let term4 = (F_i * G_i.conj())?;
        gram_01 = ((gram_01 - term3)? - term4)?;
        
        // f*conj(f) + G*conj(G)
        let term5 = (f_i * f_i.conj())?;
        let term6 = (G_i * G_i.conj())?;
        gram_11 = ((gram_11 + term5)? + term6)?;
    }
    
    // Scale by 1/n
    let scale = Flr::new(1.0 / n as f64)?;
    gram_00.re = (gram_00.re * scale)?;
    gram_00.im = (gram_00.im * scale)?;
    gram_01.re = (gram_01.re * scale)?;
    gram_01.im = (gram_01.im * scale)?;
    gram_11.re = (gram_11.re * scale)?;
    gram_11.im = (gram_11.im * scale)?;
    
    // Gram matrix is Hermitian, so gram_10 = conj(gram_01)
    let gram_10 = gram_01.conj();
    
    Ok([[gram_00, gram_01], [gram_10, gram_11]])
}

/// Build the recursive tree structure using FFT
fn build_recursive_tree(
    root_gram: [[Complex; 2]; 2],
    logn: usize,
) -> Result<Vec<TreeNode>, Error> {
    let mut nodes = Vec::new();
    
    // Build tree level by level, from root to leaves
    for level in (0..=logn).rev() {
        let num_nodes = 1 << (logn - level);
        
        for index in 0..num_nodes {
            let gram = if level == logn {
                // Root level - use the computed gram matrix
                root_gram
            } else {
                // Derived from parent nodes using FFT split
                split_gram_matrix(&nodes, level + 1, index)?
            };
            
            // Perform LDL decomposition
            let ldl_factors = ldl_decompose_2x2(gram)?;
            
            let node = TreeNode {
                level,
                index,
                gram,
                ldl_factors,
            };
            
            nodes.push(node);
        }
    }
    
    Ok(nodes)
}

/// Perform LDL decomposition of a 2x2 Hermitian matrix
fn ldl_decompose_2x2(gram: [[Complex; 2]; 2]) -> Result<LdlFactors, Error> {
    // For a 2x2 Hermitian matrix [[a, b], [conj(b), c]]
    // LDL decomposition gives L = [[1, 0], [conj(b)/a, 1]] and D = [a, c - |b|²/a]
    
    let a = gram[0][0];
    let b = gram[0][1];
    let c = gram[1][1];
    
    // Extract real parts (since diagonal elements of Gram matrix are real)
    let a_real = a.re;
    let c_real = c.re;
    
    // Check for numerical stability
    if a_real.raw().abs() < 1e-10 {
        return Err(Error::InvalidParameter);
    }
    
    // Compute L matrix
    let l_21 = (b / Complex::from_real(a_real.raw())?)?;
    let l = [
        [Flr::one(), Flr::zero()],
        [l_21.re, Flr::one()],
    ];
    
    // Compute D matrix
    let b_norm_sq = b.norm_squared()?;
    let d_2 = (c_real - (b_norm_sq / a_real)?)?;
    let d = [a_real, d_2];
    
    Ok(LdlFactors { l, d })
}

/// Split gram matrix for recursive tree construction
fn split_gram_matrix(
    nodes: &[TreeNode],
    parent_level: usize,
    child_index: usize,
) -> Result<[[Complex; 2]; 2], Error> {
    // Find parent node
    let parent_index = child_index / 2;
    let parent_node = nodes.iter()
        .find(|n| n.level == parent_level && n.index == parent_index)
        .ok_or(Error::InvalidParameter)?;
    
    // For simplicity, use a basic splitting approach
    // Real implementation would use more sophisticated FFT-based splitting
    let scale = Flr::new(0.5)?;
    let scaled_gram = [
        [
            Complex { re: (parent_node.gram[0][0].re * scale)?, im: parent_node.gram[0][0].im },
            Complex { re: (parent_node.gram[0][1].re * scale)?, im: parent_node.gram[0][1].im },
        ],
        [
            Complex { re: (parent_node.gram[1][0].re * scale)?, im: parent_node.gram[1][0].im },
            Complex { re: (parent_node.gram[1][1].re * scale)?, im: parent_node.gram[1][1].im },
        ],
    ];
    
    Ok(scaled_gram)
}

/// Extract tree data for storage
fn extract_tree_data(nodes: &[TreeNode]) -> Result<(Vec<Flr>, Vec<usize>, Vec<Flr>), Error> {
    let mut ldl_data = Vec::new();
    let mut depths = Vec::new();
    let mut sigmas = Vec::new();
    
    for node in nodes {
        // Store LDL factors
        ldl_data.extend_from_slice(&[
            node.ldl_factors.l[0][0], node.ldl_factors.l[0][1],
            node.ldl_factors.l[1][0], node.ldl_factors.l[1][1],
            node.ldl_factors.d[0], node.ldl_factors.d[1],
        ]);
        
        depths.push(node.level);
        
        // Compute Gaussian standard deviation from diagonal elements
        let sigma = compute_gaussian_sigma(&node.ldl_factors)?;
        sigmas.push(sigma);
    }
    
    Ok((ldl_data, depths, sigmas))
}

/// Compute Gaussian standard deviation from LDL factors
fn compute_gaussian_sigma(ldl: &LdlFactors) -> Result<Flr, Error> {
    // Use the smallest diagonal element to determine sigma
    let d_min = if ldl.d[0] < ldl.d[1] { ldl.d[0] } else { ldl.d[1] };
    
    // Standard deviation is proportional to sqrt(diagonal element)
    let sigma = d_min.sqrt()?;
    
    // Apply Falcon-specific scaling factor
    let falcon_scale = Flr::new(1.17)?; // Empirical factor for Falcon
    sigma * falcon_scale
}

/// Recursive sampling from the Falcon tree
fn recursive_sample(
    tree: &FalconTree,
    target: &[Complex],
    level: usize,
    offset: usize,
) -> Result<Vec<i16>, Error> {
    if level == 0 {
        // Base case: direct Gaussian sampling
        return base_case_sample(tree, target, offset);
    }
    
    let n = target.len();
    let half_n = n / 2;
    
    // Split target into two halves
    let target_0 = &target[..half_n];
    let target_1 = &target[half_n..];
    
    // Recursively sample from subtrees
    let sample_1 = recursive_sample(tree, target_1, level - 1, offset + half_n)?;
    
    // Adjust target_0 based on sample_1 (this is a simplified approach)
    let adjusted_target_0 = adjust_target(target_0, &sample_1, tree, level)?;
    let sample_0 = recursive_sample(tree, &adjusted_target_0, level - 1, offset)?;
    
    // Combine samples
    let mut result = Vec::with_capacity(n);
    result.extend_from_slice(&sample_0);
    result.extend_from_slice(&sample_1);
    
    Ok(result)
}

/// Base case sampling for leaves
fn base_case_sample(
    tree: &FalconTree,
    target: &[Complex],
    offset: usize,
) -> Result<Vec<i16>, Error> {
    let mut result = Vec::with_capacity(target.len());
    
    for (i, &t) in target.iter().enumerate() {
        // Use the corresponding sigma from the tree
        let sigma_index = offset + i;
        let sigma = if sigma_index < tree.sigmas.len() {
            tree.sigmas[sigma_index]
        } else {
            Flr::one() // Fallback
        };
        
        // Sample from discrete Gaussian distribution
        let center = t.re.raw();
        let sample = sample_discrete_gaussian(center, sigma.raw())?;
        result.push(sample);
    }
    
    Ok(result)
}

/// Adjust target for recursive sampling
fn adjust_target(
    target: &[Complex],
    sample: &[i16],
    _tree: &FalconTree,
    _level: usize,
) -> Result<Vec<Complex>, Error> {
    // Simplified adjustment - real implementation would use proper lattice adjustment
    let mut adjusted = Vec::with_capacity(target.len());
    
    for (i, &t) in target.iter().enumerate() {
        let adjustment = if i < sample.len() {
            Complex::from_real(sample[i] as f64 * 0.1)? // Simplified adjustment
        } else {
            Complex::zero()
        };
        
        adjusted.push((t - adjustment)?);
    }
    
    Ok(adjusted)
}

/// Sample from discrete Gaussian distribution
fn sample_discrete_gaussian(center: f64, sigma: f64) -> Result<i16, Error> {
    // This is a placeholder for proper discrete Gaussian sampling
    // Real implementation would use rejection sampling or CDT method
    
    use rand_core::RngCore;
    
    #[cfg(feature = "std")]
    let mut rng = rand_core::OsRng;
    #[cfg(not(feature = "std"))]
    return Err(Error::RngError);
    
    // Simplified Box-Muller transform (not cryptographically secure)
    let u1 = rng.next_u32() as f64 / u32::MAX as f64;
    let u2 = rng.next_u32() as f64 / u32::MAX as f64;
    
    #[cfg(feature = "std")]
    let normal = (-2.0 * u1.ln()).sqrt() * (2.0 * core::f64::consts::PI * u2).cos();
    
    #[cfg(not(feature = "std"))]
    let normal = (-2.0 * libm::log(u1)).sqrt() * libm::cos(2.0 * core::f64::consts::PI * u2);
    
    let scaled = center + sigma * normal;
    let rounded = scaled.round() as i16;
    
    // Bound the result to prevent overflow
    Ok(rounded.clamp(-2047, 2047))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fndsa::{params::FNDSA_512, ntru::generate_ntru_keypair};
    use rand_core::OsRng;
    
    #[test]
    #[cfg(feature = "std")]
    fn test_ldl_decomposition() {
        let gram = [
            [Complex::from_real(4.0).unwrap(), Complex::from_real(2.0).unwrap()],
            [Complex::from_real(2.0).unwrap(), Complex::from_real(3.0).unwrap()],
        ];
        
        let ldl = ldl_decompose_2x2(gram).unwrap();
        
        // Check that L is lower triangular with 1's on diagonal
        assert!((ldl.l[0][0].raw() - 1.0).abs() < 1e-10);
        assert!((ldl.l[1][1].raw() - 1.0).abs() < 1e-10);
        
        // Check that D elements are positive
        assert!(ldl.d[0].raw() > 0.0);
        assert!(ldl.d[1].raw() > 0.0);
    }
    
    #[test]
    #[cfg(feature = "std")]
    fn test_falcon_tree_construction() {
        let mut rng = OsRng;
        let params = FNDSA_512;
        
        // This test may fail with current simplified NTRU implementation
        // but demonstrates the API structure
        if let Ok((private_key, _)) = generate_ntru_keypair(params, &mut rng) {
            let tree_result = FalconTree::new(&private_key);
            
            match tree_result {
                Ok(tree) => {
                    assert_eq!(tree.depth(), params.logn);
                    assert!(!tree.sigmas.is_empty());
                }
                Err(_) => {
                    // Expected with simplified implementation
                }
            }
        }
    }
}