// Copyright 2014-2015 Johannes Köster, Peer Aramillo Irizar.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::borrow::Borrow;
use std::mem;

use bit_set::BitSet;
use vec_map::VecMap;

pub type SymbolRanks = VecMap<u8>;

/// Representation of an alphabet.
#[derive(Debug, PartialEq)]
pub struct Alphabet {
    pub symbols: BitSet,
}

impl Alphabet {
    pub fn new<C, T>(symbols: T) -> Self
    where
        C: Borrow<u8>,
        T: IntoIterator<Item = C>,
    {
        let mut s = BitSet::new();
        s.extend(symbols.into_iter().map(|c| *c.borrow() as usize));

        Alphabet { symbols: s }
    }

    pub fn insert(&mut self, a: u8) {
        self.symbols.insert(a as usize);
    }

    pub fn is_word<C, T>(&self, text: T) -> bool
    where
        C: Borrow<u8>,
        T: IntoIterator<Item = C>,
    {
        text.into_iter()
            .all(|c| self.symbols.contains(*c.borrow() as usize))
    }

    pub fn max_symbol(&self) -> Option<u8> {
        self.symbols.iter().max().map(|a| a as u8)
    }

    pub fn len(&self) -> usize {
        self.symbols.len()
    }

    pub fn is_empty(&self) -> bool {
        self.symbols.is_empty()
    }

    pub fn intersection(&self, other: &Alphabet) -> Self {
        return Alphabet {
            symbols: self.symbols.intersection(&other.symbols).collect(),
        };
    }

    pub fn difference(&self, other: &Alphabet) -> Self {
        return Alphabet {
            symbols: self.symbols.difference(&other.symbols).collect(),
        };
    }

    pub fn union(&self, other: &Alphabet) -> Self {
        return Alphabet {
            symbols: self.symbols.union(&other.symbols).collect(),
        };
    }
}

#[derive(Serialize, Deserialize)]
pub struct RankTransform {
    pub ranks: SymbolRanks,
}

impl RankTransform {
    pub fn new(alphabet: &Alphabet) -> Self {
        let mut ranks = VecMap::new();
        for (r, c) in alphabet.symbols.iter().enumerate() {
            ranks.insert(c, r as u8);
        }

        RankTransform { ranks }
    }

    pub fn get(&self, a: u8) -> u8 {
        *self.ranks.get(a as usize).expect("Unexpected character.")
    }

    pub fn transform<C, T>(&self, text: T) -> Vec<u8>
    where
        C: Borrow<u8>,
        T: IntoIterator<Item = C>,
    {
        text.into_iter()
            .map(|c| {
                *self
                    .ranks
                    .get(*c.borrow() as usize)
                    .expect("Unexpected character in text.")
            })
            .collect()
    }

    pub fn qgrams<C, T>(&self, q: u32, text: T) -> QGrams<'_, C, T::IntoIter>
    where
        C: Borrow<u8>,
        T: IntoIterator<Item = C>,
    {
        let bits = (self.ranks.len() as f32).log2().ceil() as u32;
        assert!(
            (bits * q) as usize <= mem::size_of::<usize>() * 8,
            "Expecting q to be smaller than usize / log2(|A|)"
        );

        let mut qgrams = QGrams {
            text: text.into_iter(),
            ranks: self,
            bits,
            mask: (1 << (q * bits)) - 1,
            qgram: 0,
        };

        for _ in 0..q - 1 {
            qgrams.next();
        }

        qgrams
    }

    pub fn alphabet(&self) -> Alphabet {
        let mut symbols = BitSet::with_capacity(self.ranks.len());
        symbols.extend(self.ranks.keys());
        Alphabet { symbols }
    }

    pub fn get_width(&self) -> usize {
        (self.ranks.len() as f32).log2().ceil() as usize
    }
}

/// Iterator over q-grams.
pub struct QGrams<'a, C, T>
where
    C: Borrow<u8>,
    T: Iterator<Item = C>,
{
    text: T,
    ranks: &'a RankTransform,
    bits: u32,
    mask: usize,
    qgram: usize,
}

impl<'a, C, T> QGrams<'a, C, T>
where
    C: Borrow<u8>,
    T: Iterator<Item = C>,
{
    /// Push a new character into the current qgram.
    fn qgram_push(&mut self, a: u8) {
        self.qgram <<= self.bits;
        self.qgram |= a as usize;
        self.qgram &= self.mask;
    }
}

impl<'a, C, T> Iterator for QGrams<'a, C, T>
where
    C: Borrow<u8>,
    T: Iterator<Item = C>,
{
    type Item = usize;

    fn next(&mut self) -> Option<usize> {
        match self.text.next() {
            Some(a) => {
                let b = self.ranks.get(*a.borrow());
                self.qgram_push(b);
                Some(self.qgram)
            }
            None => None,
        }
    }
}

/// Returns the english ascii lower case alphabet.
pub fn english_ascii_lower_alphabet() -> Alphabet {
    Alphabet::new(&b"abcdefghijklmnopqrstuvwxyz"[..])
}

/// Returns the english ascii upper case alphabet.
pub fn english_ascii_upper_alphabet() -> Alphabet {
    Alphabet::new(&b"ABCDEFGHIJKLMNOPQRSTUVWXYZ"[..])
}
