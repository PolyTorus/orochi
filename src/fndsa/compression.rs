//! Signature compression and decompression for FN-DSA
//!
//! This module implements the Golomb-Rice encoding used to compress
//! FN-DSA signatures to their minimal size.

use crate::fndsa::params::*;
use crate::Error;

#[cfg(not(feature = "std"))]
use alloc::vec::Vec;

/// Signature format header
const SIGNATURE_HEADER_512: u8 = 0x30 + LOGN_512 as u8;
const SIGNATURE_HEADER_1024: u8 = 0x30 + LOGN_1024 as u8;

/// FN-DSA signature binary format
#[derive(Debug, Clone)]
pub struct CompressedSignature {
    /// Parameter set identifier
    pub params: FnDsaParams,
    /// Random nonce (salt) for hash-to-point
    pub nonce: [u8; SALT_LEN],
    /// Compressed signature polynomial s2
    pub compressed_s2: Vec<u8>,
}

impl CompressedSignature {
    /// Create a new compressed signature
    pub fn new(params: FnDsaParams, nonce: [u8; SALT_LEN], s2: &[i16]) -> Result<Self, Error> {
        let compressed_s2 = compress_signature_polynomial(s2)?;
        
        Ok(Self {
            params,
            nonce,
            compressed_s2,
        })
    }
    
    /// Serialize to binary format
    pub fn to_bytes(&self) -> Vec<u8> {
        let header = match self.params.logn {
            LOGN_512 => SIGNATURE_HEADER_512,
            LOGN_1024 => SIGNATURE_HEADER_1024,
            _ => 0x30, // Fallback
        };
        
        let mut bytes = Vec::with_capacity(self.params.signature_len);
        
        // Header (1 byte)
        bytes.push(header);
        
        // Nonce (40 bytes)
        bytes.extend_from_slice(&self.nonce);
        
        // Compressed signature
        bytes.extend_from_slice(&self.compressed_s2);
        
        // Pad to exact signature length
        bytes.resize(self.params.signature_len, 0);
        
        bytes
    }
    
    /// Deserialize from binary format
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, Error> {
        if bytes.len() < 1 + SALT_LEN {
            return Err(Error::InvalidSignature);
        }
        
        // Parse header
        let header = bytes[0];
        let logn = (header - 0x30) as usize;
        let params = FnDsaParams::from_logn(logn).ok_or(Error::InvalidSignature)?;
        
        if bytes.len() != params.signature_len {
            return Err(Error::InvalidSignature);
        }
        
        // Extract nonce
        let mut nonce = [0u8; SALT_LEN];
        nonce.copy_from_slice(&bytes[1..1 + SALT_LEN]);
        
        // Extract compressed signature
        let compressed_s2 = bytes[1 + SALT_LEN..].to_vec();
        
        Ok(Self {
            params,
            nonce,
            compressed_s2,
        })
    }
    
    /// Decompress the signature polynomial
    pub fn decompress_s2(&self) -> Result<Vec<i16>, Error> {
        decompress_signature_polynomial(&self.compressed_s2, self.params.n)
    }
}

/// Compress signature polynomial using Golomb-Rice encoding
pub fn compress_signature_polynomial(s2: &[i16]) -> Result<Vec<u8>, Error> {
    let mut compressed = Vec::new();
    let mut bit_buffer = 0u32;
    let mut bit_count = 0u8;
    
    for &coeff in s2 {
        // Check coefficient bounds
        if coeff < -COEFF_BOUND || coeff > COEFF_BOUND {
            return Err(Error::InvalidParameter);
        }
        
        // Encode coefficient using Golomb-Rice method
        encode_coefficient(coeff, &mut bit_buffer, &mut bit_count, &mut compressed)?;
    }
    
    // Flush remaining bits
    if bit_count > 0 {
        compressed.push((bit_buffer << (8 - bit_count)) as u8);
    }
    
    Ok(compressed)
}

/// Decompress signature polynomial from Golomb-Rice encoding
pub fn decompress_signature_polynomial(data: &[u8], n: usize) -> Result<Vec<i16>, Error> {
    let mut result = Vec::with_capacity(n);
    let mut byte_index = 0;
    let mut bit_buffer = 0u32;
    let mut bit_count = 0u8;
    
    for _ in 0..n {
        // Decode next coefficient
        let coeff = decode_coefficient(data, &mut byte_index, &mut bit_buffer, &mut bit_count)?;
        result.push(coeff);
    }
    
    // Verify that remaining bits are zero (padding)
    verify_padding(data, byte_index, bit_buffer, bit_count)?;
    
    Ok(result)
}

/// Encode a single coefficient using Golomb-Rice encoding
fn encode_coefficient(
    coeff: i16,
    bit_buffer: &mut u32,
    bit_count: &mut u8,
    output: &mut Vec<u8>,
) -> Result<(), Error> {
    // Extract sign and absolute value
    let sign = if coeff < 0 { 1u32 } else { 0u32 };
    let abs_val = coeff.unsigned_abs() as u32;
    
    // Encode sign bit + low 7 bits (8 bits total)
    let low_bits = (sign << 7) | (abs_val & 0x7F);
    write_bits(low_bits, 8, bit_buffer, bit_count, output)?;
    
    // Encode high bits using unary encoding
    let high_bits = abs_val >> 7;
    
    // Write 'high_bits' zero bits followed by one '1' bit
    for _ in 0..high_bits {
        write_bits(0, 1, bit_buffer, bit_count, output)?;
    }
    write_bits(1, 1, bit_buffer, bit_count, output)?; // Stop bit
    
    Ok(())
}

/// Decode a single coefficient from Golomb-Rice encoding
fn decode_coefficient(
    data: &[u8],
    byte_index: &mut usize,
    bit_buffer: &mut u32,
    bit_count: &mut u8,
) -> Result<i16, Error> {
    // Read sign bit and low 7 bits
    let low_data = read_bits(data, 8, byte_index, bit_buffer, bit_count)?;
    let sign = (low_data >> 7) & 1;
    let mut abs_val = low_data & 0x7F;
    
    // Decode unary-encoded high bits
    loop {
        let bit = read_bits(data, 1, byte_index, bit_buffer, bit_count)?;
        if bit == 1 {
            break; // Stop bit found
        }
        abs_val += 0x80; // Add 128 for each zero bit
        
        if abs_val > COEFF_BOUND as u32 {
            return Err(Error::InvalidSignature);
        }
    }
    
    // Reject invalid "-0" encoding
    if sign == 1 && abs_val == 0 {
        return Err(Error::InvalidSignature);
    }
    
    // Apply sign
    let result = if sign == 1 {
        -(abs_val as i16)
    } else {
        abs_val as i16
    };
    
    Ok(result)
}

/// Write bits to output buffer
fn write_bits(
    value: u32,
    num_bits: u8,
    bit_buffer: &mut u32,
    bit_count: &mut u8,
    output: &mut Vec<u8>,
) -> Result<(), Error> {
    if num_bits > 32 {
        return Err(Error::InvalidParameter);
    }
    
    let mask = if num_bits == 32 { 0xFFFFFFFF } else { (1u32 << num_bits) - 1 };
    let value = value & mask;
    
    *bit_buffer = (*bit_buffer << num_bits) | value;
    *bit_count += num_bits;
    
    // Output complete bytes
    while *bit_count >= 8 {
        *bit_count -= 8;
        let byte = (*bit_buffer >> *bit_count) as u8;
        output.push(byte);
        *bit_buffer &= (1u32 << *bit_count) - 1;
    }
    
    Ok(())
}

/// Read bits from input buffer
fn read_bits(
    data: &[u8],
    num_bits: u8,
    byte_index: &mut usize,
    bit_buffer: &mut u32,
    bit_count: &mut u8,
) -> Result<u32, Error> {
    if num_bits > 32 {
        return Err(Error::InvalidParameter);
    }
    
    // Ensure we have enough bits in the buffer
    while *bit_count < num_bits {
        if *byte_index >= data.len() {
            return Err(Error::InvalidSignature);
        }
        
        *bit_buffer = (*bit_buffer << 8) | (data[*byte_index] as u32);
        *byte_index += 1;
        *bit_count += 8;
    }
    
    // Extract the requested bits
    *bit_count -= num_bits;
    let mask = if num_bits == 32 { 0xFFFFFFFF } else { (1u32 << num_bits) - 1 };
    let result = (*bit_buffer >> *bit_count) & mask;
    *bit_buffer &= (1u32 << *bit_count) - 1;
    
    Ok(result)
}

/// Verify that padding bits are zero
fn verify_padding(
    data: &[u8],
    byte_index: usize,
    bit_buffer: u32,
    bit_count: u8,
) -> Result<(), Error> {
    // Check remaining bits in buffer are zero
    if bit_count > 0 && bit_buffer != 0 {
        return Err(Error::InvalidSignature);
    }
    
    // Check remaining bytes are zero
    for i in byte_index..data.len() {
        if data[i] != 0 {
            return Err(Error::InvalidSignature);
        }
    }
    
    Ok(())
}

/// Calculate the expected compressed size for a signature
pub fn estimate_compressed_size(s2: &[i16]) -> usize {
    let mut total_bits = 0;
    
    for &coeff in s2 {
        let abs_val = coeff.unsigned_abs() as u32;
        
        // 8 bits for sign + low 7 bits
        total_bits += 8;
        
        // Unary encoding for high bits
        let high_bits = abs_val >> 7;
        total_bits += high_bits + 1; // high_bits zeros + 1 stop bit
    }
    
    // Convert to bytes (round up)
    ((total_bits + 7) / 8) as usize
}

/// Signature statistics for analysis
#[derive(Debug)]
pub struct SignatureStats {
    pub mean_abs_coeff: f64,
    pub max_abs_coeff: u16,
    pub compressed_size: usize,
    pub compression_ratio: f64,
}

/// Analyze signature compression efficiency
pub fn analyze_signature(s2: &[i16]) -> SignatureStats {
    let sum: u32 = s2.iter().map(|&x| x.unsigned_abs() as u32).sum();
    let mean_abs_coeff = sum as f64 / s2.len() as f64;
    
    let max_abs_coeff = s2.iter()
        .map(|&x| x.unsigned_abs())
        .max()
        .unwrap_or(0);
    
    let compressed_size = estimate_compressed_size(s2);
    let uncompressed_size = s2.len() * 2; // 2 bytes per coefficient
    let compression_ratio = uncompressed_size as f64 / compressed_size as f64;
    
    SignatureStats {
        mean_abs_coeff,
        max_abs_coeff,
        compressed_size,
        compression_ratio,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_compress_decompress_simple() {
        let s2 = vec![1, -2, 0, 127, -128, 3, -1, 0];
        
        let compressed = compress_signature_polynomial(&s2).unwrap();
        let decompressed = decompress_signature_polynomial(&compressed, s2.len()).unwrap();
        
        assert_eq!(s2, decompressed);
    }
    
    #[test]
    fn test_compress_decompress_boundary_values() {
        let s2 = vec![COEFF_BOUND, -COEFF_BOUND, 0, 1, -1];
        
        let compressed = compress_signature_polynomial(&s2).unwrap();
        let decompressed = decompress_signature_polynomial(&compressed, s2.len()).unwrap();
        
        assert_eq!(s2, decompressed);
    }
    
    #[test]
    fn test_invalid_coefficient_bounds() {
        let s2 = vec![COEFF_BOUND + 1]; // Too large
        let result = compress_signature_polynomial(&s2);
        assert!(result.is_err());
        
        let s2 = vec![-COEFF_BOUND - 1]; // Too small
        let result = compress_signature_polynomial(&s2);
        assert!(result.is_err());
    }
    
    #[test]
    fn test_signature_format_roundtrip() {
        let params = FNDSA_512;
        let nonce = [42u8; SALT_LEN];
        let s2: Vec<i16> = (0..params.n).map(|i| (i % 17) as i16 - 8).collect();
        
        let compressed_sig = CompressedSignature::new(params, nonce, &s2).unwrap();
        let bytes = compressed_sig.to_bytes();
        
        assert_eq!(bytes.len(), params.signature_len);
        assert_eq!(bytes[0], SIGNATURE_HEADER_512);
        
        let recovered = CompressedSignature::from_bytes(&bytes).unwrap();
        assert_eq!(recovered.params.logn, params.logn);
        assert_eq!(recovered.nonce, nonce);
        
        let recovered_s2 = recovered.decompress_s2().unwrap();
        assert_eq!(recovered_s2, s2);
    }
    
    #[test]
    fn test_compression_statistics() {
        let s2: Vec<i16> = (0..512).map(|i| (i % 7) as i16 - 3).collect();
        let stats = analyze_signature(&s2);
        
        assert!(stats.mean_abs_coeff >= 0.0);
        assert!(stats.max_abs_coeff <= 3);
        assert!(stats.compression_ratio > 1.0); // Should achieve some compression
        assert_eq!(stats.compressed_size, estimate_compressed_size(&s2));
    }
    
    #[test]
    fn test_bit_operations() {
        let mut bit_buffer = 0u32;
        let mut bit_count = 0u8;
        let mut output = Vec::new();
        
        // Write some bits
        write_bits(0b101, 3, &mut bit_buffer, &mut bit_count, &mut output).unwrap();
        write_bits(0b11000, 5, &mut bit_buffer, &mut bit_count, &mut output).unwrap();
        
        // Flush
        if bit_count > 0 {
            output.push((bit_buffer << (8 - bit_count)) as u8);
        }
        
        // Read back
        let mut byte_index = 0;
        let mut bit_buffer = 0u32;
        let mut bit_count = 0u8;
        
        let val1 = read_bits(&output, 3, &mut byte_index, &mut bit_buffer, &mut bit_count).unwrap();
        let val2 = read_bits(&output, 5, &mut byte_index, &mut bit_buffer, &mut bit_count).unwrap();
        
        assert_eq!(val1, 0b101);
        assert_eq!(val2, 0b11000);
    }
}