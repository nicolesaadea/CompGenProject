from Implementations.CuckooFilter import CuckooFilter

def get_kmers(sequence, k):
    """Generate all k-mers from a given sequence."""
    for i in range(len(sequence) - k + 1):
        yield sequence[i:i + k]

cf = CuckooFilter(bucket_size=4, num_buckets=50000)

# Example sequence
sequence = "ACTGACTGACTGACTGACTG"
k = 4 

# Insert k-mers into the filter
print("Inserting k-mers into the filter:")
for kmer in get_kmers(sequence, k):
    result = cf.insert(kmer)
    print(f"Inserted {kmer}: {result}")

# Test lookups for k-mers that should exist
print("\nLooking up k-mers that exist:")
existing_kmers = ["ACTG", "TGAC", "CTGA"]
for kmer in existing_kmers:
    result = cf.lookup(kmer)
    print(f"Lookup for {kmer}: {'Found' if result else 'Not Found'}")

# Test lookups for k-mers that do not exist
print("\nLooking up k-mers that do not exist:")
non_existing_kmers = ["GGGG", "CCCC", "AAAA"]
for kmer in non_existing_kmers:
    result = cf.lookup(kmer)
    print(f"Lookup for {kmer}: {'Found' if result else 'Not Found'}")

# Delete some k-mers
print("\nDeleting k-mers:")
kmers_to_delete = ["ACTG", "TGAC"]
for kmer in kmers_to_delete:
    result = cf.delete(kmer)
    print(f"Deleted {kmer}: {'Success' if result else 'Failed'}")

# Verify deletion
print("\nLooking up k-mers after deletion:")
for kmer in existing_kmers:
    result = cf.lookup(kmer)
    print(f"Lookup for {kmer} after deletion: {'Found' if result else 'Not Found'}")
