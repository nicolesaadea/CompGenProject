
from CountMinSketchImp import KMerTracker

kmer_tracker = KMerTracker(cuckoo_bucket_size=4, cuckoo_num_buckets=1000, cms_width=1000, cms_depth=5)

# Simulating adding k-mers
kmer_tracker.add_kmer("ACGT")
kmer_tracker.add_kmer("ACGT")
kmer_tracker.add_kmer("ACGA")
kmer_tracker.add_kmer("ACGT")  # Duplicate, should be ignored by Cuckoo Filter

print("Frequency of k-mer 'ACGT':", kmer_tracker.get_count("ACGT"))  # Expected output: 3
print("Frequency of k-mer 'ACGA':", kmer_tracker.get_count("ACGA"))  # Expected output: 1