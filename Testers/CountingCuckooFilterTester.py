from Implementations.CountingCuckooFilter import CountingCuckooFilter

# Now you can use CountingCuckooFilter
filter = CountingCuckooFilter(bucket_size=4, num_buckets=20)
print(filter)


def get_kmers(sequence, k):
    """Generate all k-mers from a given sequence."""
    for i in range(len(sequence) - k + 1):
        yield sequence[i:i + k]

cf = CountingCuckooFilter(bucket_size=4, num_buckets=50000)

# Example sequence
sequence = "ACTGACTGACTGACTGACTG"
k = 4

# Insert k-mers into the filter
print("Inserting k-mers into the filter:")
for kmer in get_kmers(sequence, k):
    cf.insert(kmer)
    count = cf.get_count(kmer)
    print(f"Inserted {kmer}, count: {count}")

# Test lookups for k-mers that should exist
print("\nLooking up k-mers that exist:")
existing_kmers = ["ACTG", "TGAC", "CTGA"]
for kmer in existing_kmers:
    result = cf.lookup(kmer)
    count = cf.get_count(kmer)
    print(f"Lookup for {kmer}: {'Found' if result else 'Not Found'}, count: {count}")

# Test lookups for k-mers that do not exist
print("\nLooking up k-mers that do not exist:")
non_existing_kmers = ["GGGG", "CCCC", "AAAA"]
for kmer in non_existing_kmers:
    result = cf.lookup(kmer)
    count = cf.get_count(kmer)
    print(f"Lookup for {kmer}: {'Found' if result else 'Not Found'}, count: {count}")

# Increment counts for some k-mers
print("\nIncrementing counts for existing k-mers:")
increment_kmers = ["ACTG", "TGAC", "ACTG", "CTGA"]
for kmer in increment_kmers:
    cf.insert(kmer)
    count = cf.get_count(kmer)
    print(f"Incremented {kmer}, new count: {count}")

# Delete some k-mers
print("\nDeleting k-mers:")
kmers_to_delete = ["ACTG", "TGAC"]
for kmer in kmers_to_delete:
    cf.delete(kmer)
    count = cf.get_count(kmer)
    print(f"Deleted {kmer}, remaining count: {count}")

# Verify deletion
print("\nLooking up k-mers after deletion:")
for kmer in existing_kmers:
    result = cf.lookup(kmer)
    count = cf.get_count(kmer)
    print(f"Lookup for {kmer} after deletion: {'Found' if result else 'Not Found'}, count: {count}")
