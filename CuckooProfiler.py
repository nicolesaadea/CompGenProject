import time
from memory_profiler import memory_usage
import random
from CuckooFilterImpl import CuckooFilter

def get_kmers(sequence, k):
    """Generate all k-mers from a given sequence."""
    for i in range(len(sequence) - k + 1):
        yield sequence[i:i + k]


# Profiler for a specifically defined cuckoo filter
def profile_cuckoo_filter(cf, sequence, k):
    inserted_kmers = get_kmers(sequence, k)

    # Store starting time and memory
    start_time = time.time()
    initial_memory = memory_usage()[0]

    # Insert all kmers
    for kmer in inserted_kmers:
        cf.insert(kmer)

    # Measure time and memory taken to store kmers
    end_time = time.time()
    final_memory = memory_usage()[0]
    elapsed_time = end_time - start_time
    memory_used = final_memory - initial_memory
    print(f"Insertion Time: {elapsed_time:.4f} seconds")
    print(f"Memory Used for Insertion: {memory_used:.4f} MB")

    # Lookup a random subset of 5 kmers and test time
    lookup_kmers = random.sample(inserted_kmers, min(5, len(inserted_kmers)))
    start_time = time.time()
    for kmer in lookup_kmers:
        print(f"Lookup {kmer}: {cf.lookup(kmer)}")
    end_time = time.time()
    print(f"Lookup Time: {end_time - start_time:.4f} seconds")


# Generate a sequence with a set duplication rate of kmers
def generate_sequence_with_duplication(length, k, duplication_rate):
    unique_kmers = set()
    sequence = ""
    while len(sequence) < length:
        kmer = ''.join(random.choices('ACGT', k=k))

        # Add kmer if it is a duplicate and withing the duplicate limit
        if kmer in unique_kmers and random.random() < duplication_rate:
            sequence += kmer
        else:
            unique_kmers.add(kmer)
            sequence += kmer
    return sequence

# Test filters with a specified configuration over different duplication rates
def test_duplication_rates(sequence_length, k, duplication_rates, bucket_size, num_buckets):
    for rate in duplication_rates:
        print(f"\nTesting duplication rate: {rate}")
        sequence = generate_sequence_with_duplication(sequence_length, k, rate)
        cf = CuckooFilter(bucket_size=bucket_size, num_buckets=num_buckets)
        profile_cuckoo_filter(cf, sequence, k)

# Test filters with a specified configuration over different duplication rates
def test_kmer_sizes(sequence, kmer_sizes, bucket_size, num_buckets):
    for k in kmer_sizes:
        print(f"\nTesting k-mer size: {k}")
        cf = CuckooFilter(bucket_size=bucket_size, num_buckets=num_buckets)
        profile_cuckoo_filter(cf, sequence, k)

# Test filters with a variety of configurations
def test_filter_configurations(sequence, k, configurations):
    for config in configurations:
        bucket_size, num_buckets, max_kicks = config
        print(f"Testing configuration: bucket_size={bucket_size}, num_buckets={num_buckets}, max_kicks={max_kicks}")
        cf = CuckooFilter(bucket_size=bucket_size, num_buckets=num_buckets, max_kicks=max_kicks)
        profile_cuckoo_filter(cf, sequence, k)