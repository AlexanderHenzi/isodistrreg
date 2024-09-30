use crate::structures::Direction;
use crate::total_order::routines::pool_partitions_from_right;
use crate::total_order::structures::WeightedPartition;
use std::iter::{Flatten, RepeatN, Scan, repeat_n};
use std::mem;
use std::vec::IntoIter;

pub fn algorithm<D: Direction>(
    data: impl IntoIterator<Item = (f64, f64, f64)>,
) -> PartitionIterator {
    let mut allocated: Vec<_> = data.into_iter().collect();
    allocated.sort_unstable_by(|l, r| {
        // Sort by covariate and response
        l.0.total_cmp(&r.0).then(l.1.total_cmp(&r.1))
    });

    let deduplicated = allocated.chunk_by(|l, r| l.0 == r.0).map(|chunk| {
        let (weighted_sum, weight) = chunk.iter().fold(
            (0.0, 0.0),
            |(acc_sum_prod, acc_weight), &(_, response, weight)| {
                (acc_sum_prod + weight * response, acc_weight + weight)
            },
        );
        (weighted_sum / weight, weight)
    });

    algorithm_pre_sorted_deduplicated::<D>(deduplicated)
}

pub fn algorithm_pre_sorted_deduplicated<D: Direction>(
    data: impl Iterator<Item = (f64, f64)>,
) -> PartitionIterator {
    let mut partitions = Vec::new();

    for (index, (response, weight)) in data.enumerate() {
        if weight > 0.0 {
            partitions.push(WeightedPartition {
                index: index + 1,
                weight,
                value: response,
            });
            pool_partitions_from_right::<D>(&mut partitions);
        } else if weight == 0.0 {
            // the value might not be meaningful
            if let Some(last) = partitions.last_mut() {
                last.index += 1;
            }
        } else {
            panic!("can't perform isotonic regression with negative or NAN weights");
        }
    }

    let n = partitions.last().map(|p| p.index).unwrap_or(0);
    PartitionIterator {
        inner: partitions
            .into_iter()
            .scan(0, expand_partition as ScanFn)
            .flatten(),
        // input data corresponds to unique covariate indices
        n,
    }
}

pub struct PartitionIterator {
    inner: Flatten<Scan<IntoIter<WeightedPartition>, usize, ScanFn>>,
    /// Number of elements this iterator will produce
    n: usize,
}

impl Iterator for PartitionIterator {
    type Item = f64;

    fn next(&mut self) -> Option<Self::Item> {
        let value = self.inner.next();
        if value.is_some() {
            self.n -= 1;
        }
        value
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.n, Some(self.n))
    }
}

impl ExactSizeIterator for PartitionIterator {
    fn len(&self) -> usize {
        self.n
    }
}

fn expand_partition(start_index: &mut usize, partition: WeightedPartition) -> Option<RepeatN<f64>> {
    let previous_index = mem::replace(start_index, partition.index);
    Some(repeat_n(partition.value, partition.index - previous_index))
}
type ScanFn = fn(&mut usize, WeightedPartition) -> Option<RepeatN<f64>>;

#[cfg(test)]
mod test {
    use crate::structures::{Decreasing, Increasing};
    use crate::test::is_relative_eq_vec;
    use crate::total_order::tonic_regression::algorithm_pre_sorted_deduplicated;
    use std::iter::repeat;

    #[test]
    fn test_isotone() {
        assert!(is_relative_eq_vec(
            &algorithm_pre_sorted_deduplicated::<Increasing>(
                [1.0, 2.0, 3.0].into_iter().zip(repeat(1.0))
            )
            .collect::<Vec<_>>(),
            &[1.0, 2.0, 3.0],
        ));
    }

    #[test]
    fn test_antitone() {
        assert!(is_relative_eq_vec(
            &algorithm_pre_sorted_deduplicated::<Decreasing>(
                [1.0, 2.0, 3.0].into_iter().zip(repeat(1.0))
            )
            .collect::<Vec<_>>(),
            &[2.0, 2.0, 2.0],
        ));
    }

    #[test]
    fn test_antitone_weighted() {
        assert!(is_relative_eq_vec(
            &algorithm_pre_sorted_deduplicated::<Decreasing>(
                [1.0, 2.0, 3.0].into_iter().zip([1.0, 2.0, 4.0].into_iter())
            )
            .collect::<Vec<_>>(),
            &[17.0 / 7.0; 3],
        ));
    }

    #[test]
    fn test_zero_weighted() {
        assert!(is_relative_eq_vec(
            &algorithm_pre_sorted_deduplicated::<Decreasing>(
                [1.0, 2.0, 99.0, 3.0]
                    .into_iter()
                    .zip([1.0, 2.0, 0.0, 4.0].into_iter())
            )
            .collect::<Vec<_>>(),
            &[17.0 / 7.0; 4],
        ));
    }
}
