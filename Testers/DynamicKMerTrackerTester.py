from Implementations.DynamicKMerTracker import DynamicKMerTracker

def get_kmers(sequence, k):
    """Generate all k-mers from a given sequence."""
    for i in range(len(sequence) - k + 1):
        yield sequence[i:i + k]

# Initialize the DynamicKMerTracker
cuckoo_bucket_size = 4
cuckoo_num_buckets = 50000
cms_width = 1000
cms_depth = 5
initial_freq_threshold = 3
use_threshold = 0.75
tracker = DynamicKMerTracker(cuckoo_bucket_size, cuckoo_num_buckets, cms_width, cms_depth, initial_freq_threshold, use_threshold)

print("Inserting k-mers into the tracker:")
low_frequency_kmers = ["ACTG", "TGAC"]
high_frequency_kmers = ["CTGA", "GACT"]

# Insert low-frequency k-mers
for kmer in low_frequency_kmers:
    tracker.add_kmer(kmer)
    frequency = tracker.get_kmer_frequency(kmer)
    print(f"Inserted {kmer} once, frequency: {frequency}")
    tracker.add_kmer(kmer)
    frequency = tracker.get_kmer_frequency(kmer)
    print(f"Inserted {kmer} twice, frequency: {frequency}")

# Insert high-frequency k-mers
for kmer in high_frequency_kmers:
    for _ in range(5):
        tracker.add_kmer(kmer)
    frequency = tracker.get_kmer_frequency(kmer)
    print(f"Inserted {kmer} five times, frequency: {frequency}")

# Lookup all k-mers to verify frequencies
print("\nLooking up all k-mers:")
all_kmers = low_frequency_kmers + high_frequency_kmers
for kmer in all_kmers:
    frequency = tracker.get_kmer_frequency(kmer)
    print(f"{kmer}: frequency {frequency}")

# Test dynamic adjustment during insertion
test_sequence = "ACTGCTGACTGACTGACTGCTGACTG"
print("\nTesting dynamic adjustment:")
for kmer in get_kmers(test_sequence, k=4):
    tracker.add_kmer(kmer)
    print(f"Inserted {kmer}, current threshold: {tracker.freq_threshold}, current use rate: {tracker.use_rate():.2f}")