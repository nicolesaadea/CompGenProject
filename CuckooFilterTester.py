from CuckooFilterImpl import CuckooFilter

""""""

def get_kmers(sequence, k):
    """Generate all k-mers from a given sequence."""
    for i in range(len(sequence) - k + 1):
        yield sequence[i:i + k]

# Initialize Cuckoo Filter for genomic data
cf = CuckooFilter(bucket_size=4, num_buckets=50000)  # Adjust size as needed

# Example sequence and k-mer size
sequence = "ACTGACTGACTGACTGACTG"
k = 4  # You might choose larger k for real applications

# Insert k-mers into the filter
for kmer in get_kmers(sequence, k):
    cf.insert(kmer)

# Lookup specific k-mers
print(cf.lookup("ACTG"))  # Should return True if "ACTG" is in the sequence
print(cf.lookup("TGAC"))  # Check for presence of another k-mer