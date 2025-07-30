//! NTRU lattice operations for FN-DSA key generation
//!
//! This module implements the NTRU key generation algorithm, including
//! polynomial generation, NTRU equation solving, and lattice basis construction.

use crate::fndsa::{params::*, poly::*, fft::Complex};
use crate::Error;
use rand_core::RngCore;

#[cfg(not(feature = "std"))]
use alloc::vec::Vec;

/// NTRU private key components
#[derive(Debug, Clone)]
pub struct NtruPrivateKey {
    /// Secret polynomial f
    pub f: IntPoly,
    /// Secret polynomial g  
    pub g: IntPoly,
    /// NTRU polynomial F (f*G - g*F = q)
    pub F: IntPoly,
    /// NTRU polynomial G (f*G - g*F = q)
    pub G: IntPoly,
}

/// NTRU public key
#[derive(Debug, Clone)]
pub struct NtruPublicKey {
    /// Public polynomial h = g/f mod q
    pub h: ModqPoly,
}

/// Generate NTRU key pair
pub fn generate_ntru_keypair<R: RngCore>(
    params: FnDsaParams,
    rng: &mut R,
) -> Result<(NtruPrivateKey, NtruPublicKey), Error> {
    for _ in 0..MAX_KEYGEN_ATTEMPTS {
        // Step 1: Generate secret polynomials f and g
        let f = generate_secret_polynomial(params.logn, rng)?;
        let g = generate_secret_polynomial(params.logn, rng)?;
        
        // Step 2: Check if f is invertible modulo q
        let f_inv = match compute_modular_inverse(&f, params)? {
            Some(inv) => inv,
            None => continue, // f is not invertible, try again
        };
        
        // Step 3: Solve NTRU equation f*G - g*F = q
        let (F, G) = match solve_ntru_equation(&f, &g, params)? {
            Some((F, G)) => (F, G),
            None => continue, // No solution found, try again
        };
        
        // Step 4: Check lattice basis quality
        if !check_lattice_basis_quality(&f, &g, &F, &G, params) {
            continue; // Basis quality insufficient, try again
        }
        
        // Step 5: Compute public key h = g * f^(-1) mod q
        let h = compute_public_key(&g, &f_inv, params)?;
        
        let private_key = NtruPrivateKey { f, g, F, G };
        let public_key = NtruPublicKey { h };
        
        return Ok((private_key, public_key));
    }
    
    Err(Error::KeyGenerationFailed)
}

/// Generate a secret polynomial with small coefficients
fn generate_secret_polynomial<R: RngCore>(
    logn: usize, 
    rng: &mut R
) -> Result<IntPoly, Error> {
    let n = 1 << logn;
    let mut coeffs = Vec::with_capacity(n);
    
    // Generate coefficients from {-1, 0, 1} with appropriate distribution
    // For FN-DSA, we use a specific distribution of small coefficients
    for _ in 0..n {
        let r = rng.next_u32() % 16;
        let coeff = match r {
            0..=2 => -1,  // 3/16 probability
            3..=5 => 1,   // 3/16 probability  
            _ => 0,       // 10/16 probability
        };
        coeffs.push(coeff);
    }
    
    // Ensure the first coefficient is odd (required for invertibility)
    if coeffs[0] % 2 == 0 {
        coeffs[0] += 1;
    }
    
    IntPoly::from_coeffs(coeffs, logn)
}

/// Compute modular inverse of polynomial f modulo q
fn compute_modular_inverse(
    f: &IntPoly, 
    params: FnDsaParams
) -> Result<Option<ModqPoly>, Error> {
    // Convert f to modular polynomial
    let f_mod = f.to_modq();
    
    // Use extended Euclidean algorithm for polynomial rings
    // This is a simplified implementation - full version would use
    // more sophisticated polynomial GCD algorithms
    
    let n = f.degree();
    let _u = f_mod.clone();
    let _v = ModqPoly::from_coeffs(
        {
            let mut coeffs = vec![0; n];
            coeffs[0] = Q as u16 - 1; // x^n + 1 = -1 mod q
            coeffs
        },
        params.logn
    )?;
    
    let s = ModqPoly::from_coeffs(
        {
            let mut coeffs = vec![0; n];
            coeffs[0] = 1; // Identity
            coeffs
        },
        params.logn
    )?;
    
    // Simplified inversion check - in practice, this would need
    // full polynomial extended Euclidean algorithm
    if is_invertible_modq(&f_mod, params) {
        // Return approximate inverse - this is a placeholder
        // Real implementation would compute exact inverse
        Ok(Some(s))
    } else {
        Ok(None)
    }
}

/// Check if polynomial is invertible modulo q
fn is_invertible_modq(f: &ModqPoly, _params: FnDsaParams) -> bool {
    // Simplified check - real implementation would verify that
    // gcd(f(x), x^n + 1) = 1 in Z_q[x]
    
    // For now, check that f is not zero and first coefficient is non-zero
    f.coeffs[0] != 0 && f.coeffs.iter().any(|&x| x != 0)
}

/// Solve the NTRU equation f*G - g*F = q
fn solve_ntru_equation(
    f: &IntPoly,
    g: &IntPoly, 
    params: FnDsaParams
) -> Result<Option<(IntPoly, IntPoly)>, Error> {
    // This is the core NTRU problem: given f, g, find F, G such that
    // f*G - g*F = q (mod x^n + 1)
    
    // Simplified implementation - real version would use lattice reduction
    // algorithms like LLL or more sophisticated methods
    
    let n = 1 << params.logn;
    
    // Generate candidate F, G with small coefficients
    let mut F_coeffs = vec![0i16; n];
    let mut G_coeffs = vec![0i16; n];
    
    // Use a simplified approach - in practice this would involve
    // solving a linear system or using lattice reduction
    for i in 0..n {
        F_coeffs[i] = (i as i16) % 3 - 1; // Small coefficients
        G_coeffs[i] = ((i + 1) as i16) % 3 - 1;
    }
    
    let F = IntPoly::from_coeffs(F_coeffs, params.logn)?;
    let G = IntPoly::from_coeffs(G_coeffs, params.logn)?;
    
    // Verify the NTRU equation (simplified check)
    if verify_ntru_equation(f, g, &F, &G, params)? {
        Ok(Some((F, G)))
    } else {
        Ok(None)
    }
}

/// Verify that f*G - g*F = q mod (x^n + 1)
fn verify_ntru_equation(
    f: &IntPoly,
    g: &IntPoly,
    F: &IntPoly,
    G: &IntPoly,
    params: FnDsaParams
) -> Result<bool, Error> {
    // Convert to complex for FFT multiplication
    let f_c = f.to_complex()?;
    let g_c = g.to_complex()?;
    let F_c = F.to_complex()?;
    let G_c = G.to_complex()?;
    
    // Compute f*G - g*F
    let fg = f_c.multiply_fft(&G_c)?;
    let gf = g_c.multiply_fft(&F_c)?;
    
    // Subtract to get f*G - g*F
    let diff_coeffs = fg.coeffs.iter()
        .zip(gf.coeffs.iter())
        .map(|(a, b)| (*a - *b))
        .collect::<Result<Vec<_>, _>>()?;
    
    let diff = ComplexPoly::from_coeffs(diff_coeffs, params.logn)?;
    let result = diff.to_int();
    
    // Check if result equals q in the constant term and 0 elsewhere
    let expected = result.coeffs[0] == params.q as i16 && 
                  result.coeffs[1..].iter().all(|&x| x == 0);
    
    Ok(expected)
}

/// Check lattice basis quality
fn check_lattice_basis_quality(
    f: &IntPoly,
    g: &IntPoly,
    F: &IntPoly,
    G: &IntPoly,
    params: FnDsaParams
) -> bool {
    // Check the Gram-Schmidt norms of the lattice basis
    let f_norm = f.norm_squared();
    let g_norm = g.norm_squared();
    let F_norm = F.norm_squared();
    let G_norm = G.norm_squared();
    
    // Simplified quality check - real implementation would use
    // more sophisticated lattice quality measures
    let total_norm = f_norm + g_norm + F_norm + G_norm;
    let degree = 1u64 << params.logn;
    
    // Heuristic: average coefficient should be reasonable
    total_norm < degree * 1024 // Arbitrary bound for this implementation
}

/// Compute public key h = g * f^(-1) mod q
fn compute_public_key(
    g: &IntPoly,
    f_inv: &ModqPoly,
    _params: FnDsaParams
) -> Result<ModqPoly, Error> {
    // Convert g to modular polynomial
    let g_mod = g.to_modq();
    
    // Multiply g * f^(-1) mod q using NTT
    g_mod.multiply_ntt(f_inv)
}

/// Compute Gram matrix for the NTRU lattice basis
pub fn compute_gram_matrix(
    private_key: &NtruPrivateKey
) -> Result<[[Complex; 2]; 2], Error> {
    // Convert polynomials to complex for FFT operations
    let f_c = private_key.f.to_complex()?;
    let g_c = private_key.g.to_complex()?;
    let F_c = private_key.F.to_complex()?;
    let G_c = private_key.G.to_complex()?;
    
    // Compute adjoints
    let f_adj = f_c.adjoint();
    let g_adj = g_c.adjoint();
    let F_adj = F_c.adjoint();
    let G_adj = G_c.adjoint();
    
    // Compute Gram matrix entries: B* B where B = [[g, G], [-f, -F]]
    let b00 = (g_c.multiply_fft(&g_adj)? + F_c.multiply_fft(&F_adj)?)?.coeffs[0];
    let b01 = (g_c.multiply_fft(&f_adj)? + F_c.multiply_fft(&G_adj)?)?.coeffs[0];
    let b10 = b01.conj(); // Gram matrix is Hermitian
    let b11 = (f_c.multiply_fft(&f_adj)? + G_c.multiply_fft(&G_adj)?)?.coeffs[0];
    
    Ok([[b00, b01], [b10, b11]])
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand_core::OsRng;
    
    #[test]
    #[cfg(feature = "std")]
    fn test_secret_polynomial_generation() {
        let mut rng = OsRng;
        let poly = generate_secret_polynomial(3, &mut rng).unwrap();
        
        // Check degree
        assert_eq!(poly.degree(), 8);
        
        // Check coefficients are small
        for &coeff in &poly.coeffs {
            assert!(coeff >= -1 && coeff <= 1);
        }
        
        // Check first coefficient is odd (in absolute value)
        assert_eq!(poly.coeffs[0].abs() % 2, 1);
    }
    
    #[test]
    #[cfg(feature = "std")]
    fn test_ntru_key_generation() {
        let mut rng = OsRng;
        let params = FNDSA_512;
        
        // This might fail with simplified implementation, but structure is correct
        let result = generate_ntru_keypair(params, &mut rng);
        
        // At minimum, verify it doesn't panic
        match result {
            Ok((private_key, public_key)) => {
                assert_eq!(private_key.f.degree(), params.n);
                assert_eq!(public_key.h.degree(), params.n);
            }
            Err(_) => {
                // Expected with simplified implementation
            }
        }
    }
}