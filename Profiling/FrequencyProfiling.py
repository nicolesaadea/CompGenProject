import sys
import os
import random
import itertools
import matplotlib.pyplot as plt
from Bio import SeqIO
from Implementations.CuckooFilter import CuckooFilter
from Implementations.EnhancedKMerTracker import EnhancedKMerTracker
from Implementations.DynamicKMerTracker import DynamicKMerTracker

def get_kmers(sequence, k):
    """Get all k-mers from a sequence."""
    for i in range(len(sequence) - k + 1):
        yield sequence[i:i + k]

def generate_sequence_with_frequencies(length, k, max_frequency):
    """Generate a sequence with controlled k-mer frequencies."""
    num_kmers = length // k
    kmers = [''.join(random.choices('ACGT', k=k)) for _ in range(num_kmers)]
    kmer_frequencies = {kmer: random.randint(1, max_frequency) for kmer in kmers}
    sequence = ''.join(itertools.chain.from_iterable([kmer * freq for kmer, freq in kmer_frequencies.items()]))
    random.shuffle(list(sequence))
    return sequence, kmer_frequencies

def profile_frequency_tracking(tracker, kmer_frequencies, margin=1):
    """
    Profile the frequency tracking capabilities of a k-mer tracker.
    
    Args:
        tracker (object): The k-mer tracker to test (e.g., CuckooFilter, EnhancedKMerTracker).
        kmer_frequencies (dict): Ground truth k-mer frequencies {kmer: frequency}.
        margin (int): Acceptable margin of error for true positives.
    
    Returns:
        dict: Accuracy metrics including FPR, FNR, and overall accuracy.
    """
    true_positives = 0
    false_positives = 0
    false_negatives = 0
    total_kmers = len(kmer_frequencies)
    
    # Insert k-mers into the tracker
    for kmer, freq in kmer_frequencies.items():
        for _ in range(freq):
            tracker.add_kmer(kmer)
    
    # Query and compare frequencies
    for kmer, true_freq in kmer_frequencies.items():
        predicted_freq = tracker.get_kmer_frequency(kmer)
        if abs(predicted_freq - true_freq) <= margin:
            true_positives += 1
        elif predicted_freq > true_freq:
            false_positives += 1
        elif predicted_freq < true_freq:
            false_negatives += 1
    
    # Calculate metrics
    accuracy = true_positives / total_kmers
    fpr = false_positives / (false_positives + total_kmers - true_positives) if total_kmers > true_positives else 0
    fnr = false_negatives / (false_negatives + true_positives) if true_positives > 0 else 1
    
    print(f"Accuracy: {accuracy:.4f}")
    print(f"False Positive Rate: {fpr:.4f}")
    print(f"False Negative Rate: {fnr:.4f}")
    
    return {
        "accuracy": accuracy,
        "fpr": fpr,
        "fnr": fnr,
        "true_positives": true_positives,
        "false_positives": false_positives,
        "false_negatives": false_negatives
    }

def compare_frequency_tracking_tools(sequence_length, k, max_frequency):
    """
    Compare different tools for their frequency tracking capabilities.
    
    Args:
        sequence_length (int): Length of the sequence.
        k (int): k-mer size.
        max_frequency (int): Maximum frequency for k-mers.
    """
    sequence, kmer_frequencies = generate_sequence_with_frequencies(sequence_length, k, max_frequency)
    
    trackers = {
        "CuckooFilter": CuckooFilter(bucket_size=4, num_buckets=50000),
        "EnhancedKMerTracker": EnhancedKMerTracker(
            cuckoo_bucket_size=4,
            cuckoo_num_buckets=50000,
            cms_width=1000,
            cms_depth=5,
            freq_threshold=10
        ),
        "DynamicKMerTracker": DynamicKMerTracker(
            cuckoo_bucket_size=4,
            cuckoo_num_buckets=50000,
            cms_width=1000,
            cms_depth=5,
            initial_freq_threshold=10,
            use_threshold=0.75,
            smoothing_window=10
        )
    }
    
    results = {}
    for name, tracker in trackers.items():
        print(f"\nTesting {name} for frequency tracking...")
        results[name] = profile_frequency_tracking(tracker, kmer_frequencies)
    
    # Display results
    print("\nFrequency Tracking Results:")
    for name, metrics in results.items():
        print(f"{name}:")
        print(f"  Accuracy: {metrics['accuracy']:.4f}")
        print(f"  False Positive Rate: {metrics['fpr']:.4f}")
        print(f"  False Negative Rate: {metrics['fnr']:.4f}")

    # Visualization of metrics
    names = list(results.keys())
    accuracies = [metrics['accuracy'] for metrics in results.values()]
    fprs = [metrics['fpr'] for metrics in results.values()]
    fnrs = [metrics['fnr'] for metrics in results.values()]

    x = range(len(names))
    plt.figure(figsize=(12, 6))
    plt.bar(x, accuracies, color='green', label='Accuracy')
    plt.bar(x, fprs, bottom=accuracies, color='red', label='False Positive Rate')
    plt.bar(x, fnrs, bottom=[a + f for a, f in zip(accuracies, fprs)], color='blue', label='False Negative Rate')
    plt.xticks(x, names)
    plt.ylabel("Proportion")
    plt.title("Frequency Tracking Metrics")
    plt.legend()
    plt.grid(True)
    plt.show()

def main():
    """Main execution function."""
    random.seed(42)
    
    # Frequency tracking testing.
    print("\nComparing Frequency Tracking Tools:")
    compare_frequency_tracking_tools(sequence_length=10000, k=21, max_frequency=10)

if __name__ == "__main__":
    main()
