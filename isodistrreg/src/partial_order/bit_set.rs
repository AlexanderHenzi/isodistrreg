use std::array;
use std::fmt::{Debug, Display, Formatter};
use std::ops::{BitAnd, BitAndAssign, BitOr, BitOrAssign, BitXor, BitXorAssign, Not};

type Word = u64;

#[derive(Clone, Eq, PartialEq, Hash)]
pub struct BitSet<const B: usize> {
    words: [Word; B],
}

impl<const B: usize> BitSet<B> {
    #[must_use]
    #[inline]
    pub fn new() -> Self {
        Self { words: [0; B] }
    }
    #[must_use]
    #[inline]
    pub fn singleton(index: usize) -> Self {
        debug_assert!(index < Self::capacity());

        let mut words = [0; B];
        let (word_index, bit) = Self::address(index);
        words[word_index] = 1 << bit;
        Self { words }
    }
    #[must_use]
    pub fn fill(n: usize) -> Self {
        debug_assert!(n <= Self::capacity());

        let (word_index, bit) = Self::address(n);
        let mut words = [0; B];

        words[..word_index].fill(Word::MAX);
        if bit > 0 {
            words[word_index] = (1_u64 << bit) - 1;
        }

        Self { words }
    }
    #[must_use]
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.words.iter().all(|&w| w == 0)
    }
    #[inline(always)]
    pub fn add(&mut self, index: usize) {
        debug_assert!(index < Self::capacity());
        let (word_index, bit) = Self::address(index);
        self.words[word_index] |= 1 << bit;
    }
    #[inline(always)]
    pub fn remove(&mut self, index: usize) {
        debug_assert!(index < Self::capacity());
        let (word_index, bit) = Self::address(index);
        self.words[word_index] &= !(1 << bit);
    }
    #[must_use]
    #[inline]
    pub fn contains(&self, index: &usize) -> bool {
        if *index < Self::capacity() {
            let (word_index, bit) = Self::address(*index);
            self.words[word_index] & (1 << bit) != 0
        } else {
            false
        }
    }
    #[must_use]
    #[inline]
    pub fn min(&self) -> Option<usize> {
        for (word_index, word) in self.words.iter().enumerate() {
            if *word != 0 {
                return Some(word_index * Self::bits_per_word() + word.trailing_zeros() as usize);
            }
        }
        None
    }
    #[must_use]
    #[inline]
    pub fn max(&self) -> Option<usize> {
        for (word_index, word) in self.words.iter().enumerate().rev() {
            if *word != 0 {
                return Some(
                    (word_index + 1) * Self::bits_per_word() - word.leading_zeros() as usize,
                );
            }
        }
        None
    }
    #[must_use]
    #[inline]
    pub fn len(&self) -> usize {
        self.words
            .iter()
            .map(|word| word.count_ones() as usize)
            .sum()
    }
    #[must_use]
    pub fn iter(&self) -> BitSetIter<'_, B> {
        BitSetIter {
            word_index: 0,
            base: 0,
            word_state: self.words[0],
            ptr: self,
        }
    }
    #[must_use]
    #[inline(always)]
    pub const fn address(index: usize) -> (usize, usize) {
        (index / Self::bits_per_word(), index % Self::bits_per_word())
    }
    #[must_use]
    #[inline(always)]
    pub const fn capacity() -> usize {
        B * Self::bits_per_word()
    }
    #[must_use]
    #[inline(always)]
    pub const fn bits_per_word() -> usize {
        Word::BITS as usize
    }
}
impl<const B: usize> Default for BitSet<B> {
    fn default() -> Self {
        Self::new()
    }
}
impl<'a, const B: usize> IntoIterator for &'a BitSet<B> {
    type Item = usize;
    type IntoIter = BitSetIter<'a, B>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

impl<const B: usize> BitSet<B> {
    #[must_use]
    #[inline(always)]
    pub fn complement(&self) -> Self {
        Self {
            words: array::from_fn(|i| !self.words[i]),
        }
    }
    #[must_use]
    #[inline(always)]
    pub fn difference(&self, other: &Self) -> Self {
        self.bin_op(other, |l, r| l & !r)
    }
    #[must_use]
    #[inline(always)]
    pub fn intersection(&self, other: &Self) -> Self {
        self.bin_op(other, |l, r| l & r)
    }
    #[must_use]
    #[inline(always)]
    pub fn union(&self, other: &Self) -> Self {
        self.bin_op(other, |l, r| l | r)
    }
    #[must_use]
    #[inline(always)]
    pub fn with_singleton(&self, index: usize) -> Self {
        let mut with = self.clone();
        with.add(index);
        with
    }
    #[must_use]
    #[inline(always)]
    pub fn symmetric_difference(&self, other: &Self) -> Self {
        self.bin_op(other, |l, r| l ^ r)
    }
    #[inline(always)]
    fn bin_op(&self, other: &Self, bin_op: impl Fn(Word, Word) -> Word) -> Self {
        Self {
            words: array::from_fn(|i| bin_op(self.words[i], other.words[i])),
        }
    }
    #[must_use]
    pub fn is_subset_of(&self, other: &Self) -> bool {
        self.words
            .iter()
            .zip(other.words.iter())
            .all(|(&sub, &sup)| sub & sup == sub)
    }
    #[must_use]
    pub fn is_superset_of(&self, other: &Self) -> bool {
        // or sub & sup == sub
        self.words
            .iter()
            .zip(other.words.iter())
            .all(|(&sup, &sub)| sup | sub == sup)
    }
}

impl<const B: usize> Not for BitSet<B> {
    type Output = Self;

    #[inline(always)]
    fn not(self) -> Self::Output {
        self.complement()
    }
}

impl<const B: usize> BitOr for BitSet<B> {
    type Output = Self;

    #[inline(always)]
    fn bitor(mut self, rhs: Self) -> Self::Output {
        BitOrAssign::bitor_assign(&mut self, rhs);
        self
    }
}

impl<const B: usize> BitOrAssign for BitSet<B> {
    #[inline(always)]
    fn bitor_assign(&mut self, rhs: Self) {
        *self = self.union(&rhs);
    }
}

impl<const B: usize> BitAnd for BitSet<B> {
    type Output = Self;

    #[inline(always)]
    fn bitand(mut self, rhs: Self) -> Self::Output {
        BitAndAssign::bitand_assign(&mut self, rhs);
        self
    }
}

impl<const B: usize> BitAndAssign for BitSet<B> {
    #[inline(always)]
    fn bitand_assign(&mut self, rhs: Self) {
        *self = self.intersection(&rhs);
    }
}

impl<const B: usize> BitXor for BitSet<B> {
    type Output = Self;

    #[inline(always)]
    fn bitxor(mut self, rhs: Self) -> Self::Output {
        BitXorAssign::bitxor_assign(&mut self, rhs);
        self
    }
}

impl<const B: usize> BitXorAssign for BitSet<B> {
    #[inline(always)]
    fn bitxor_assign(&mut self, rhs: Self) {
        *self = self.symmetric_difference(&rhs);
    }
}

impl<const B: usize> FromIterator<usize> for BitSet<B> {
    fn from_iter<T: IntoIterator<Item = usize>>(iter: T) -> Self {
        let mut container = BitSet::new();
        for item in iter {
            container.add(item);
        }
        container
    }
}

impl<const B: usize> Display for BitSet<B> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.write_str("{")?;
        match self.len() {
            0 => {}
            1 => {
                let only = self.iter().next().unwrap();
                Debug::fmt(&only, f)?;
            }
            _ => {
                let mut items = self.iter();
                let first = items.next().unwrap();
                Debug::fmt(&first, f)?;
                for item in items {
                    f.write_str(", ")?;
                    Debug::fmt(&item, f)?;
                }
            }
        }
        f.write_str("}")?;

        Ok(())
    }
}

impl<const B: usize> Debug for BitSet<B> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.write_str("BitSet<")?;
        Debug::fmt(&B, f)?;
        f.write_str(" x ")?;
        Debug::fmt(&Word::BITS, f)?;
        f.write_str(" = ")?;
        Debug::fmt(&Self::capacity(), f)?;
        f.write_str("> ")?;
        Display::fmt(self, f)?;

        Ok(())
    }
}

#[derive(Clone, Debug)]
pub struct BitSetIter<'a, const B: usize> {
    // Index of the current word in `state`.
    word_index: usize,
    // Base bit offset for the current word: word_index * BitSet::<B>::bits_per_word().
    base: usize,
    // Remaining bits in the current word.
    word_state: Word,
    ptr: &'a BitSet<B>,
}

impl<const B: usize> Iterator for BitSetIter<'_, B> {
    // TODO: Return u32 instead?
    type Item = usize;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if self.word_state != 0 {
                let bit = self.word_state.trailing_zeros() as usize;
                self.word_state &= self.word_state - 1;
                return Some(self.base + bit);
            }

            self.word_index += 1;
            if self.word_index >= B {
                return None;
            }
            self.base += BitSet::<B>::bits_per_word();
            self.word_state = self.ptr.words[self.word_index];
        }
    }
}

#[cfg(test)]
mod test {
    use crate::partial_order::BitSet;

    #[test]
    fn test_fill() {
        for n in [0, 1, BitSet::<1>::capacity() - 1, BitSet::<1>::capacity()] {
            let filled = BitSet::<1>::fill(n);
            assert_eq!(filled.len(), n);
        }
    }
}
