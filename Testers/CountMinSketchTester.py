
from Implementations.CountMinSketch import CountMinSketch

def get_kmers(sequence, k):
    """Generate all k-mers from a given sequence."""
    for i in range(len(sequence) - k + 1):
        yield sequence[i:i + k]

# Initialize the Count-Min Sketch
cms_width = 100
cms_depth = 5
cms = CountMinSketch(cms_width, cms_depth)

# Example sequence and k-mer size
sequence = "ACTGACTGACTGACTGACTG"
k = 4

# Insert k-mers into the Count-Min Sketch
print("Inserting k-mers into Count-Min Sketch:")
for kmer in get_kmers(sequence, k):
    cms.add(kmer)
    estimated_count = cms.estimate(kmer)
    print(f"Inserted {kmer}, estimated frequency: {estimated_count}")

# Test frequency estimation for specific k-mers
print("\nEstimating frequencies for k-mers:")
test_kmers = ["ACTG", "TGAC", "CTGA", "GGGG", "CCCC"]
for kmer in test_kmers:
    estimated_count = cms.estimate(kmer)
    print(f"Estimated frequency of {kmer}: {estimated_count}")

# Increment counts for some k-mers
print("\nIncrementing frequencies for some k-mers:")
increment_kmers = ["ACTG", "TGAC", "ACTG", "CTGA", "ACTG"]
for kmer in increment_kmers:
    cms.add(kmer)
    estimated_count = cms.estimate(kmer)
    print(f"Incremented {kmer}, new estimated frequency: {estimated_count}")
