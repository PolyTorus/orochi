#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

pub fn mqpoly_small_is_invertible(logn: u32,
    f: &[i8], tmp: &mut [u16]) -> bool
{
    let n = 1usize << logn;
    mqpoly_small_to_int(logn, f, tmp);
    mqpoly_int_to_NTT(logn, tmp);
    let mut r = 0xFFFFFFFF;
    for i in 0..n {
        r &= (tmp[i] as u32).wrapping_sub(Q);
    }
    (r >> 16) != 0
}


pub fn mqpoly_div_small(logn: u32,
    f: &[i8], g: &[i8], h: &mut [u16], tmp: &mut [u16])
{
    let n = 1usize << logn;
    mqpoly_small_to_int(logn, f, tmp);
    mqpoly_small_to_int(logn, g, h);
    mqpoly_int_to_NTT(logn, tmp);
    mqpoly_int_to_NTT(logn, h);
    for i in 0..n {
        h[i] = mq_div(h[i] as u32, tmp[i] as u32) as u16;
    }
    mqpoly_NTT_to_int(logn, h);
    mqpoly_int_to_ext(logn, h);
}

pub const SQBETA: [u32; 11] = [
    0,        // unused
    101498,
    208714,
    428865,
    892039,
    1852696,
    3842630,
    7959734,
    16468416,
    34034726,
    70265242,
];

const Q: u32 = 12289;

// -1/q mod 2^32
const Q1I: u32 = 4143984639;

// 2^64 mod q
const R2: u32 = 5664;


pub fn mqpoly_signed_to_ext(logn: u32, v: &[i16], d: &mut [u16]) {
    for i in 0..(1usize << logn) {
        let x = (v[i] as i32) as u32;
        d[i] = x.wrapping_add((x >> 16) & Q) as u16;
    }
}

#[inline(always)]
fn mq_add(x: u32, y: u32) -> u32 {
    // a = q - (x + y)
    // -q <= a <= q - 2  (represented as u32)
    let a = Q.wrapping_sub(x + y);

    // If a < 0, add q.
    // b = -(x + y) mod q
    // 0 <= b <= q - 1
    let b = a.wrapping_add(Q & (a >> 16));

    // q - b = x + y mod q
    // 1 <= q - b <= q
    Q - b
}

// Subtraction modulo q (internal representation).
#[inline(always)]
fn mq_sub(x: u32, y: u32) -> u32 {
    // -(q - 1) <= a <= q - 1
    let a = y.wrapping_sub(x);
    // 0 <= b <= q - 1
    let b = a.wrapping_add(Q & (a >> 16));
    Q - b
}

// Halving modulo q (internal representation).
#[inline(always)]
fn mq_half(x: u32) -> u32 {
    (x + ((x & 1).wrapping_neg() & Q)) >> 1
}
