from Implementations.CuckooFilterWithCMS import CuckooFilterWithCMS

bucket_size = 4
num_buckets = 50000
cms_width = 1000
cms_depth = 5
max_kicks = 500
cf = CuckooFilterWithCMS(bucket_size, num_buckets, cms_width, cms_depth, max_kicks)

# Example sequence
sequence = "ACTGACTGACTGACTGACTG"
k = 4

def get_kmers(sequence, k):
    """Generate all k-mers from a given sequence."""
    for i in range(len(sequence) - k + 1):
        yield sequence[i:i + k]

# Insert k-mers into the filter
print("Inserting k-mers into the filter:")
for kmer in get_kmers(sequence, k):
    result = cf.insert(kmer)
    frequency = cf.get_estimated_frequency(kmer)
    print(f"Inserted {kmer}: {'Success' if result else 'Failed'}, estimated frequency: {frequency}")

# Lookup existing k-mers
print("\nLooking up k-mers that exist:")
existing_kmers = ["ACTG", "TGAC", "CTGA"]
for kmer in existing_kmers:
    found, frequency = cf.lookup(kmer)
    print(f"Lookup for {kmer}: {'Found' if found else 'Not Found'}, estimated frequency: {frequency}")

# Lookup non-existing k-mers
print("\nLooking up k-mers that do not exist:")
non_existing_kmers = ["GGGG", "CCCC", "AAAA"]
for kmer in non_existing_kmers:
    found, frequency = cf.lookup(kmer)
    print(f"Lookup for {kmer}: {'Found' if found else 'Not Found'}, estimated frequency: {frequency}")

# Increment counts for some k-mers
print("\nIncrementing counts for existing k-mers:")
increment_kmers = ["ACTG", "TGAC", "ACTG", "CTGA"]
for kmer in increment_kmers:
    cf.insert(kmer)
    frequency = cf.get_estimated_frequency(kmer)
    print(f"Incremented {kmer}, new estimated frequency: {frequency}")

# Delete some k-mers
print("\nDeleting k-mers:")
kmers_to_delete = ["ACTG", "TGAC"]
for kmer in kmers_to_delete:
    result = cf.delete(kmer)
    frequency = cf.get_estimated_frequency(kmer)
    print(f"Deleted {kmer}: {'Success' if result else 'Failed'}, remaining estimated frequency: {frequency}")

# Verify deletion
print("\nLooking up k-mers after deletion:")
for kmer in existing_kmers:
    found, frequency = cf.lookup(kmer)
    print(f"Lookup for {kmer} after deletion: {'Found' if found else 'Not Found'}, estimated frequency: {frequency}")
