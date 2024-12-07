import sys
import os
import random
import itertools
from statistics import mean, stdev
from EnhancedKMerTracker import EnhancedKMerTracker
from DynamicKMerTracker import DynamicKMerTracker
from KMerTracker import KMerTracker


def generate_sequence_with_frequencies(length, k, max_frequency):
    """Generate a sequence with controlled k-mer frequencies."""
    num_kmers = length // k
    kmers = [''.join(random.choices('ACGT', k=k)) for _ in range(num_kmers)]
    kmer_frequencies = {kmer: random.randint(1, max_frequency) for kmer in kmers}
    sequence = ''.join(itertools.chain.from_iterable([kmer * freq for kmer, freq in kmer_frequencies.items()]))
    random.shuffle(list(sequence))
    return sequence, kmer_frequencies


def profile_frequency_tracking(tracker, kmer_frequencies, margin=1):
    """Profile the frequency tracking capabilities of a k-mer tracker."""
    true_positives = 0
    false_positives = 0
    false_negatives = 0
    total_kmers = len(kmer_frequencies)
    
    for kmer, freq in kmer_frequencies.items():
        for _ in range(freq):
            tracker.add_kmer(kmer)
    
    for kmer, true_freq in kmer_frequencies.items():
        predicted_freq = tracker.get_kmer_frequency(kmer)
        if abs(predicted_freq - true_freq) <= margin:
            true_positives += 1
        elif predicted_freq > true_freq:
            false_positives += 1
        elif predicted_freq < true_freq:
            false_negatives += 1
    
    accuracy = true_positives / total_kmers
    fpr = false_positives / (false_positives + total_kmers - true_positives) if total_kmers > true_positives else 0
    fnr = false_negatives / (false_negatives + true_positives) if true_positives > 0 else 1
    
    return {
        "accuracy": accuracy,
        "fpr": fpr,
        "fnr": fnr
    }


def compare_frequency_tracking_tools(sequence_lengths, k_values, max_frequencies):
    """Compare different tools for frequency tracking over varied datasets."""
    results = {tool: [] for tool in ["KMerTracker", "EnhancedKMerTracker", "DynamicKMerTracker"]}

    for seq_len, k, max_freq in itertools.product(sequence_lengths, k_values, max_frequencies):
        print(f"\nSequence Length: {seq_len}, k: {k}, Max Frequency: {max_freq}")
        sequence, kmer_frequencies = generate_sequence_with_frequencies(seq_len, k, max_freq)

        trackers = {
            "KMerTracker": KMerTracker(cuckoo_bucket_size=4, cuckoo_num_buckets=50000, cms_width=1000, cms_depth=5),
            "EnhancedKMerTracker": EnhancedKMerTracker(cuckoo_bucket_size=4, cuckoo_num_buckets=50000, cms_width=1000, cms_depth=5, freq_threshold=10),
            "DynamicKMerTracker": DynamicKMerTracker(cuckoo_bucket_size=4, cuckoo_num_buckets=50000, cms_width=1000, cms_depth=5, initial_freq_threshold=10, use_threshold=0.75, smoothing_window=10),
        }

        for name, tracker in trackers.items():
            metrics = profile_frequency_tracking(tracker, kmer_frequencies)
            results[name].append(metrics)
            print(f"{name} Metrics:")
            for metric, value in metrics.items():
                print(f"  {metric.capitalize()}: {value:.4f}")

    return results


def main():
    """Main execution function."""
    random.seed(42)
    
    # Define varied datasets
    sequence_lengths = [5000, 10000, 20000]  # Reduced to avoid long runtime
    k_values = [21]
    max_frequencies = [10, 20, 50]

    # Compare tools
    print("\nComparing Frequency Tracking Tools...")
    compare_frequency_tracking_tools(sequence_lengths, k_values, max_frequencies)


if __name__ == "__main__":
    main()