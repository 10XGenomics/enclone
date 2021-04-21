//! Wavelet Matrix data structure for DNA alphabet.
//! The implementation is based on the paper
//! [Claude Francisco and Gonzalo Navarro. The wavelet matrix. SPIRE (2012)](https://doi.org/10.1007/978-3-642-34109-0_18)
//!
//! # Example
//!
//! ```
//! use bio::data_structures::wavelet_matrix::WaveletMatrix;
//! let text = b"AANGGT$ACCNTT$";
//! let wm = WaveletMatrix::new(text);
//! assert_eq!(wm.rank(b'A', 0), 1);
//! assert_eq!(wm.rank(b'G', 9), 2);
//! assert_eq!(wm.rank(b'T', 13), 3);
//! ```

use crate::data_structures::rank_select::RankSelect;
use bv::BitVec;
use bv::BitsMut;

const DNA2INT: [u8; 128] = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //  0
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 10
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 20
    0, 0, 0, 0, 0, 0, 5, 0, 0, 0, // 30
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, // 40
    2, 3, 4, 5, 6, 7, 0, 0, 0, 0, // 50
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, // 60
    0, 2, 0, 0, 0, 0, 0, 0, 4, 0, // 70
    0, 0, 0, 0, 3, 0, 0, 0, 0, 0, // 80
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, // 90
    0, 0, 0, 2, 0, 0, 0, 0, 0, 0, // 100
    4, 0, 0, 0, 0, 0, 3, 0, 0, 0, // 110
    0, 0, 0, 0, 0, 0, 0, 0,
]; // 120

#[derive(Serialize, Deserialize)]
pub struct WaveletMatrix {
    width: usize,  // levels[0].len()
    height: usize, // zeros.len() == levels.len()
    zeros: Vec<u64>,
    levels: Vec<RankSelect>,
}

fn build_partlevel(
    vals: &Vec<u8>,
    shift: u8,
    next_zeros: &mut Vec<u8>,
    next_ones: &mut Vec<u8>,
    bits: &mut BitVec<u8>,
    prev_bits: u64,
) {
    let mut p = prev_bits;
    for val in vals {
        let bit = ((DNA2INT[usize::from(*val)] >> shift) & 1) == 1; // get shifted lsb
        bits.set_bit(p, bit);
        p += 1;
        if bit {
            next_ones.push(*val);
        } else {
            next_zeros.push(*val);
        }
    }
}

impl WaveletMatrix {
    /// Construct a new instance of the wavelet matrix of given text of length n (DNA alphabet plus sentinel symbol).
    /// Complexity: O(n).
    pub fn new(text: &[u8]) -> Self {
        let width = text.len();
        let height: usize = 3; // hardcoded for alphabet size <= 8 (ACGTN$)

        let mut curr_zeros: Vec<u8> = text.to_vec().clone();
        let mut curr_ones: Vec<u8> = Vec::new();

        let mut zeros: Vec<u64> = Vec::new();
        let mut levels: Vec<RankSelect> = Vec::new();

        for level in 0..height {
            let mut next_zeros: Vec<u8> = Vec::with_capacity(width);
            let mut next_ones: Vec<u8> = Vec::with_capacity(width);
            let mut curr_bits: BitVec<u8> = BitVec::new_fill(false, width as u64);
            let shift = (height - level - 1) as u8;
            build_partlevel(
                &curr_zeros,
                shift,
                &mut next_zeros,
                &mut next_ones,
                &mut curr_bits,
                0,
            );
            build_partlevel(
                &curr_ones,
                shift,
                &mut next_zeros,
                &mut next_ones,
                &mut curr_bits,
                curr_zeros.len() as u64,
            );

            curr_zeros = next_zeros;
            curr_ones = next_ones;

            let level = RankSelect::new(curr_bits, 1);
            levels.push(level);
            zeros.push(curr_zeros.len() as u64);
        }

        WaveletMatrix {
            width,
            height,
            zeros,
            levels,
        }
    }

    fn check_overflow(&self, p: u64) -> bool {
        p >= self.width as u64
    }

    fn prank(&self, level: usize, p: u64, val: u8) -> u64 {
        if p == 0 {
            0
        } else {
            if val == 0 {
                self.levels[level].rank_0(p - 1).unwrap()
            } else {
                self.levels[level].rank_1(p - 1).unwrap()
            }
        }
    }

    /// Compute the number of occurrences of symbol val in the original text up to position p (inclusive).
    /// Complexity O(1).
    pub fn rank(&self, val: u8, p: u64) -> u64 {
        if self.check_overflow(p) {
            panic!("Invalid p (it must be in range 0..wm_size-1");
        }
        let height = self.height as usize;
        let mut spos = 0;
        let mut epos = p + 1;
        for level in 0..height {
            let shift = height - level - 1;
            let bit = ((DNA2INT[val as usize] >> shift) & 1) == 1; // get shifted lsb
            if bit {
                spos = self.prank(level, spos, 1) + self.zeros[level];
                epos = self.prank(level, epos, 1) + self.zeros[level];
            } else {
                spos = self.prank(level, spos, 0);
                epos = self.prank(level, epos, 0);
            }
        }
        epos - spos
    }
}
