// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Unoptimized functions for packing and unpacking some data structures.

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn u32_bytes(x: usize) -> Vec<u8> {
    (x as u32).to_le_bytes().to_vec()
}

pub fn u32_from_bytes(x: &[u8]) -> u32 {
    u32::from_le_bytes([x[0], x[1], x[2], x[3]])
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn save_vec_string(x: &Vec<String>) -> Vec<u8> {
    let mut bytes = Vec::<u8>::new();
    bytes.append(&mut u32_bytes(x.len()));
    for i in 0..x.len() {
        bytes.append(&mut u32_bytes(x[i].len()));
        bytes.append(&mut x[i].as_bytes().to_vec());
    }
    bytes
}

pub fn restore_vec_string(x: &Vec<u8>, pos: &mut usize) -> Result<Vec<String>, ()> {
    if *pos + 4 > x.len() {
        return Err(());
    }
    let n = u32_from_bytes(&x[*pos..*pos + 4]) as usize;
    *pos += 4;
    let mut y = vec![String::new(); n];
    for j in 0..n {
        if *pos + 4 > x.len() {
            return Err(());
        }
        let k = u32_from_bytes(&x[*pos..*pos + 4]) as usize;
        *pos += 4;
        if *pos + k > x.len() {
            return Err(());
        }
        let s = String::from_utf8(x[*pos..*pos + k].to_vec());
        if s.is_err() {
            return Err(());
        }
        *pos += k;
        y[j] = s.unwrap();
    }
    Ok(y)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn save_vec_vec_u8(x: &Vec<Vec<u8>>) -> Vec<u8> {
    let mut bytes = Vec::<u8>::new();
    bytes.append(&mut u32_bytes(x.len()));
    for i in 0..x.len() {
        bytes.append(&mut u32_bytes(x[i].len()));
        bytes.append(&mut x[i].clone());
    }
    bytes
}

pub fn restore_vec_vec_u8(x: &Vec<u8>, pos: &mut usize) -> Result<Vec<Vec<u8>>, ()> {
    if *pos + 4 > x.len() {
        return Err(());
    }
    let n = u32_from_bytes(&x[*pos..*pos + 4]) as usize;
    *pos += 4;
    let mut y = vec![Vec::<u8>::new(); n];
    for j in 0..n {
        if *pos + 4 > x.len() {
            return Err(());
        }
        let k = u32_from_bytes(&x[*pos..*pos + 4]) as usize;
        *pos += 4;
        if *pos + k > x.len() {
            return Err(());
        }
        y[j] = x[*pos..*pos + k].to_vec();
        *pos += k;
    }
    Ok(y)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn save_vec_vec_u32(x: &Vec<Vec<u32>>) -> Vec<u8> {
    let mut bytes = Vec::<u8>::new();
    bytes.append(&mut u32_bytes(x.len()));
    for i in 0..x.len() {
        bytes.append(&mut u32_bytes(x[i].len()));
        for j in 0..x[i].len() {
            bytes.append(&mut x[i][j].to_le_bytes().to_vec());
        }
    }
    bytes
}

pub fn restore_vec_vec_u32(x: &Vec<u8>, pos: &mut usize) -> Result<Vec<Vec<u32>>, ()> {
    if *pos + 4 > x.len() {
        return Err(());
    }
    let n = u32_from_bytes(&x[*pos..*pos + 4]) as usize;
    *pos += 4;
    let mut y = vec![Vec::<u32>::new(); n];
    for j in 0..n {
        if *pos + 4 > x.len() {
            return Err(());
        }
        let k = u32_from_bytes(&x[*pos..*pos + 4]) as usize;
        *pos += 4;
        if *pos + 4 * k > x.len() {
            return Err(());
        }
        for _ in 0..k {
            y[j].push(u32_from_bytes(&x[*pos..*pos + 4]));
            *pos += 4;
        }
    }
    Ok(y)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn save_vec_bool(x: &Vec<bool>) -> Vec<u8> {
    let mut bytes = Vec::<u8>::new();
    bytes.append(&mut u32_bytes(x.len()));
    for i in 0..x.len() {
        bytes.push(if x[i] { 1 } else { 0 });
    }
    bytes
}

pub fn restore_vec_bool(x: &Vec<u8>, pos: &mut usize) -> Result<Vec<bool>, ()> {
    if *pos + 4 > x.len() {
        return Err(());
    }
    let n = u32_from_bytes(&x[*pos..*pos + 4]) as usize;
    *pos += 4;
    if *pos + n > x.len() {
        return Err(());
    }
    let mut y = vec![false; n];
    for j in 0..n {
        y[j] = if x[*pos] == 1 { true } else { false };
        *pos += 1;
    }
    Ok(y)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn save_vec_u32(x: &Vec<u32>) -> Vec<u8> {
    let mut bytes = Vec::<u8>::new();
    bytes.append(&mut u32_bytes(x.len()));
    for i in 0..x.len() {
        bytes.append(&mut x[i].to_le_bytes().to_vec());
    }
    bytes
}

pub fn restore_vec_u32(x: &Vec<u8>, pos: &mut usize) -> Result<Vec<u32>, ()> {
    if *pos + 4 > x.len() {
        return Err(());
    }
    let n = u32_from_bytes(&x[*pos..*pos + 4]) as usize;
    *pos += 4;
    if *pos + 4 * n > x.len() {
        return Err(());
    }
    let mut y = vec![0; n];
    for j in 0..n {
        y[j] = u32_from_bytes(&x[*pos..*pos + 4]);
        *pos += 4;
    }
    Ok(y)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn save_u32(x: u32) -> Vec<u8> {
    x.to_le_bytes().to_vec()
}

pub fn restore_u32(x: &Vec<u8>, pos: &mut usize) -> Result<u32, ()> {
    if *pos + 4 > x.len() {
        return Err(());
    }
    let y = u32_from_bytes(&x[*pos..*pos + 4]);
    *pos += 4;
    Ok(y)
}
