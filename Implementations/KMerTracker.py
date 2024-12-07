from CuckooFilter import CuckooFilter
from CountMinSketch import CountMinSketch

class KMerTracker:
    """
    k-mer Frequency Tracker combining Cuckoo Filter and Count-Min Sketch.
    
    Attributes:
        cuckoo_filter (CuckooFilter): Tracks unique k-mers
        cms (CountMinSketch): Estimates the frequency of k-mers
    """

    def __init__(self, cuckoo_bucket_size, cuckoo_num_buckets, cms_width, cms_depth):
        """
        Initialize the KMerTracker with specified parameters.

        Args:
            cuckoo_bucket_size (int): Maximum number of items per bucket in the Cuckoo Filter
            cuckoo_num_buckets (int): Total number of buckets in the Cuckoo Filter
            cms_width (int): Number of buckets per hash function row in Count-Min Sketch
            cms_depth (int): Number of hash functions in Count-Min Sketch
        """
        self.cuckoo_filter = CuckooFilter(cuckoo_bucket_size, cuckoo_num_buckets)
        self.cms = CountMinSketch(cms_width, cms_depth)

    def add_kmer(self, kmer):
        """
        Add a k-mer to the tracker, ensuring its presence in the Cuckoo Filter
        and updating its frequency in the Count-Min Sketch.

        Args:
            kmer (str): k-mer to add or update
        """
        # Check if the k-mer exists in the Cuckoo Filter
        if not self.cuckoo_filter.lookup(kmer):
            # Insert the k-mer into the Cuckoo Filter
            self.cuckoo_filter.insert(kmer)
        
        # Update frequency in the Count-Min Sketch
        self.cms.add(kmer)

    def get_kmer_frequency(self, kmer):
        """
        Retrieve the estimated frequency of a k-mer using the Count-Min Sketch.

        Args:
            kmer (str): k-mer to query
        
        Returns:
            int: Estimated frequency of the k-mer
        """
        return self.cms.estimate(kmer)
