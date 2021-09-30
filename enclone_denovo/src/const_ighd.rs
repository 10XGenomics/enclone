// Copyright (c) 2020 10x Genomics, Inc. All rights reserved.

// Code to help find IGHD constant regions.

use std::arch::x86_64::{
    __m128, __m128i, __m256, _mm256_add_epi32, _mm256_add_ps, _mm256_castps256_ps128,
    _mm256_cvtepu8_epi32, _mm256_extractf128_ps, _mm256_i32gather_ps, _mm256_mullo_epi32,
    _mm256_set1_epi32, _mm256_set1_ps, _mm256_setr_epi32, _mm_add_ps, _mm_add_ss, _mm_cvtss_f32,
    _mm_loadu_si128, _mm_movehdup_ps, _mm_movehl_ps,
};

// Score relative to IGHD PWM.

/// this version detects whether it can use the fast version at runtime.  

pub fn ighd_score(x: &[u8], ighd_pwm2: &Vec<f32>) -> f32 {
    let sub_vec_size = 5;
    let logp;
    if is_x86_feature_detected!("avx2") {
        logp = unsafe { s32_avx2(x, ighd_pwm2, sub_vec_size) };
    } else {
        logp = ighd_score_orig(x, ighd_pwm2);
    }
    logp
}

fn ighd_score_orig(x: &[u8], u: &Vec<f32>) -> f32 {
    let sub_vec_size = 5;
    let n = u.len() / sub_vec_size;
    let mut a = 0.0;
    for m in 0..n {
        a += u[m * sub_vec_size + x[m] as usize];
    }
    a
}

#[target_feature(enable = "avx2")]
unsafe fn s32_avx2(x: &[u8], u: &Vec<f32>, sub_vec_size: usize) -> f32 {
    let n = u.len() / sub_vec_size;
    let n0 = n - n % 8;
    assert!(u.len() < 1 << 16);
    // we will accumulate the sum in the 8 slots of this vector
    let mut sum_lanes: __m256 = _mm256_set1_ps(0.0);
    let mut comb = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
    comb = _mm256_mullo_epi32(comb, _mm256_set1_epi32(sub_vec_size as i32));
    let step = _mm256_set1_epi32((sub_vec_size * 8) as i32);
    // compute 8 slots at a time
    for m in (0..n0).step_by(8) {
        // Add on the x-offsets into each subvector
        let x_values = _mm_loadu_si128(x[m..(m + 8)].as_ptr() as *const __m128i);
        let x_values_epi32 = _mm256_cvtepu8_epi32(x_values);
        let offsets = _mm256_add_epi32(comb, x_values_epi32);
        // extract 8 values from u, at the final 'offsets'
        let extracted: __m256 = _mm256_i32gather_ps(u.as_ptr(), offsets, 4);
        // increment the parallel sums
        sum_lanes = _mm256_add_ps(sum_lanes, extracted);
        // move the comb over to the next group of 8 subvecs
        comb = _mm256_add_epi32(comb, step);
    }
    // add up the parallel sum values into a single number
    let mut a = hsum256_ps_avx(sum_lanes);
    // do the remaining work
    for m in n0..n {
        a += u[m * sub_vec_size + x[m] as usize];
    }
    a
}

/// Sum the slots in 8xf32 vector

#[target_feature(enable = "avx2")]
unsafe fn hsum256_ps_avx(vin: __m256) -> f32 {
    let vlow: __m128 = _mm256_castps256_ps128(vin);
    let vhigh: __m128 = _mm256_extractf128_ps(vin, 1); // high 128
    let v = _mm_add_ps(vlow, vhigh); // add the low 128
    let mut shuf = _mm_movehdup_ps(v); // broadcast elements 3,1 to 2,0
    let mut sums = _mm_add_ps(v, shuf);
    shuf = _mm_movehl_ps(shuf, sums); // high half -> low half
    sums = _mm_add_ss(sums, shuf);

    _mm_cvtss_f32(sums)
}

pub fn ighd_score2(
    x: &[u8],
    ighd_pwm: &Vec<f32>,
    max_ins_len: usize,
    ins_start: usize,
    ins_pos: &mut usize,
    ins_len: &mut usize,
) -> f32 {
    let sub_vec_size = 5;
    let n = ighd_pwm.len() / sub_vec_size;
    let mut best_score = ighd_score(x, ighd_pwm);
    let mut y = vec![0_u8; n];
    const GAP_PENALTY: f32 = 0.0_f32; // was 0.25
    for ins in (3..=3 * max_ins_len).step_by(3) {
        for i in 0..ins_start + 3 {
            y[i] = x[i];
        }
        for i in 0..ins {
            y[i + ins_start + 3] = sub_vec_size as u8 - 1;
        }
        for i in ins_start + 3..n - ins {
            y[i + ins] = x[i];
        }
        for m in 0..n - ins - 3 - ins_start {
            if m % 3 == 0 {
                let score2 = ighd_score(&y, ighd_pwm) + (ins as f32 * GAP_PENALTY);
                if score2 < best_score {
                    *ins_pos = ins_start + m + 3;
                    *ins_len = ins;
                    best_score = score2;
                }
            }
            y.swap(ins_start + m + 3, ins_start + m + ins + 3);
        }
    }
    best_score
}
