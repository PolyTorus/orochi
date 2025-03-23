pub fn trim_i8_encode(f: &[i8], nbits: u32, d: &mut [u8]) -> usize {
    let mut k = 0;
    let mut acc = 0;
    let mut acc_len = 0;
    let mask = (8u32 << nbits) - 1;

    for i in 0..f.len() {
        acc = (acc << nbits) | (((f[i] as u8) as u32) & mask);
        acc_len += nbits;

        while acc_len >= 8 {
            acc_len -= 8;
            d[k] = (acc >> acc_len) as u8;
            k += 1;
        }
    }

    if acc_len > 0 {
        d[k] = (acc << (8 - acc_len)) as u8;
        k += 1;
    }

    k
}

pub fn modq_encode(h: &[u16], d: &mut [u8]) -> usize {
    assert!((h.len() & 3) == 0);
    let mut j = 0;
    for i in 0..(h.len() >> 2) {
        let x0 = h[4 * i + 0] as u64;
        let x1 = h[4 * i + 1] as u64;
        let x2 = h[4 * i + 2] as u64;
        let x3 = h[4 * i + 3] as u64;
        let x = (x0 << 42) | (x1 << 28) | (x2 << 14) | x3;
        d[j..(j + 7)].copy_from_slice(&x.to_be_bytes()[1..8]);
        j += 7;
    }
    j
}

pub fn modq_decode(d: &[u8], h: &mut [u16]) -> Option<usize> {
    let n = h.len();
    if n == 0 {
        return Some(0);
    }
    assert!((n & 3) == 0);
    let needed = 7 * (n >> 2);
    if d.len() != needed {
        return None;
    }
    let mut ov = 0xFFFF;
    let x = ((d[0] as u64) << 48)
        | ((d[1] as u64) << 40)
        | ((d[2] as u64) << 32)
        | ((d[3] as u64) << 24)
        | ((d[4] as u64) << 16)
        | ((d[5] as u64) << 8)
        | (d[6] as u64);
    let h0 = ((x >> 42) as u32) & 0x3FFF;
    let h1 = ((x >> 28) as u32) & 0x3FFF;
    let h2 = ((x >> 14) as u32) & 0x3FFF;
    let h3 = (x as u32) & 0x3FFF;
    ov &= h0.wrapping_sub(12289);
    ov &= h1.wrapping_sub(12289);
    ov &= h2.wrapping_sub(12289);
    ov &= h3.wrapping_sub(12289);
    h[0] = h0 as u16;
    h[1] = h1 as u16;
    h[2] = h2 as u16;
    h[3] = h3 as u16;
    for i in 1..(n >> 2) {
        let x = u64::from_be_bytes(
            *<&[u8; 8]>::try_from(&d[(7 * i - 1)..(7 * i + 7)]).unwrap());
        let h0 = ((x >> 42) as u32) & 0x3FFF;
        let h1 = ((x >> 28) as u32) & 0x3FFF;
        let h2 = ((x >> 14) as u32) & 0x3FFF;
        let h3 = (x as u32) & 0x3FFF;
        ov &= h0.wrapping_sub(12289);
        ov &= h1.wrapping_sub(12289);
        ov &= h2.wrapping_sub(12289);
        ov &= h3.wrapping_sub(12289);
        h[4 * i + 0] = h0 as u16;
        h[4 * i + 1] = h1 as u16;
        h[4 * i + 2] = h2 as u16;
        h[4 * i + 3] = h3 as u16;
    }
    if (ov & 0x8000) == 0 {
        return None;
    }
    Some(needed)
}


pub fn comp_encode(s: &[i16], d: &mut [u8]) -> bool {
    let mut acc = 0;
    let mut acc_len = 0;
    let mut j = 0;
    for i in 0..s.len() {
        // Invariant: acc_len <= 7 at the beginning of each iteration.

        let x = s[i] as i32;
        if x < -2047 || x > 2047 {
            return false;
        }

        // Get sign and absolute value.
        let sw = (x >> 16) as u32;
        let w = ((x as u32) ^ sw).wrapping_sub(sw);

        // Encode sign bit then low 7 bits of the absolute value.
        acc <<= 8;
        acc |= sw & 0x80;
        acc |= w & 0x7F;
        acc_len += 8;

        // Encode the high bits. Since |x| <= 2047, the value in the high
        // bits is at most 15.
        let wh = w >> 7;
        acc <<= wh + 1;
        acc |= 1;
        acc_len += wh + 1;

        // We appended at most 8 + 15 + 1 = 24 bits, so the total number of
        // bits still fits in the 32-bit accumulator. We output complete
        // bytes.
        while acc_len >= 8 {
            acc_len -= 8;
            if j >= d.len() {
                return false;
            }
            d[j] = (acc >> acc_len) as u8;
            j += 1;
        }
    }

    // Flush remaining bits (if any).
    if acc_len > 0 {
        if j >= d.len() {
            return false;
        }
        d[j] = (acc << (8 - acc_len)) as u8;
        j += 1;
    }

    // Pad with zeros.
    for k in j..d.len() {
        d[k] = 0;
    }
    true
}

pub fn comp_decode(d: &[u8], v: &mut [i16]) -> bool {
    let mut i = 0;
    let mut acc = 0;
    let mut acc_len = 0;
    for j in 0..v.len() {
        // Invariant: acc_len <= 7 at the beginning of each iteration.

        // Get next 8 bits and split them into sign bit (s) and low bits
        // of the absolute value (m).
        if i >= d.len() {
            return false;
        }
        acc = (acc << 8) | (d[i] as u32);
        i += 1;
        let s = (acc >> (acc_len + 7)) & 1;
        let mut m = (acc >> acc_len) & 0x7F;

        // Get next bits until a 1 is reached.
        loop {
            if acc_len == 0 {
                if i >= d.len() {
                    return false;
                }
                acc = (acc << 8) | (d[i] as u32);
                i += 1;
                acc_len = 8;
            }
            acc_len -= 1;
            if ((acc >> acc_len) & 1) != 0 {
                break;
            }
            m += 0x80;
            if m > 2047 {
                return false;
            }
        }

        // Reject "-0" (invalid encoding).
        if (s & (m.wrapping_sub(1) >> 31)) != 0 {
            return false;
        }

        // Apply the sign to get the value.
        let sw = s.wrapping_neg();
        let w = (m ^ sw).wrapping_sub(sw);
        v[j] = w as i16;
    }

    // Check that unused bits are all zero.
    if acc_len > 0 {
        if (acc & ((1 << acc_len) - 1)) != 0 {
            return false;
        }
    }
    for k in i..d.len() {
        if d[k] != 0 {
            return false;
        }
    }
    true
}
