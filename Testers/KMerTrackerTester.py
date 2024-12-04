from Implementations.KMerTracker import KMerTracker

def get_kmers(sequence, k):
    """Generate all k-mers from a given sequence."""
    for i in range(len(sequence) - k + 1):
        yield sequence[i:i + k]

# Initialize KMerTracker
cuckoo_bucket_size = 4
cuckoo_num_buckets = 50000
cms_width = 1000
cms_depth = 5
tracker = KMerTracker(cuckoo_bucket_size, cuckoo_num_buckets, cms_width, cms_depth)

# Example sequence and k-mer size
sequence = "ACTGACTGACTGACTGACTG"
k = 4

# Insert k-mers into the tracker
print("Inserting k-mers into the tracker:")
for kmer in get_kmers(sequence, k):
    tracker.add_kmer(kmer)
    count = tracker.get_count(kmer)
    print(f"Inserted {kmer}, estimated frequency: {count}")

# Test lookups for existing k-mers
print("\nLooking up k-mers that exist:")
existing_kmers = ["ACTG", "TGAC", "CTGA"]
for kmer in existing_kmers:
    count = tracker.get_count(kmer)
    print(f"Lookup for {kmer}: estimated frequency {count}")

# Test lookups for non-existing k-mers
print("\nLooking up k-mers that do not exist:")
non_existing_kmers = ["GGGG", "CCCC", "AAAA"]
for kmer in non_existing_kmers:
    count = tracker.get_count(kmer)
    print(f"Lookup for {kmer}: estimated frequency {count}")

# Increment counts for existing k-mers
print("\nIncrementing counts for existing k-mers:")
increment_kmers = ["ACTG", "TGAC", "ACTG", "CTGA"]
for kmer in increment_kmers:
    tracker.add_kmer(kmer)
    count = tracker.get_count(kmer)
    print(f"Incremented {kmer}, new estimated frequency: {count}")
