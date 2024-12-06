from CountingCuckooFilter import CountingCuckooFilter
from CountMinSketch import CountMinSketch

class EnhancedKMerTracker:
    """
    k-mer Frequency Tracker combining Counting Cuckoo Filter and Count-Min Sketch. 
    Prioritizes accuracy for low-frequency k-mers and scalability for high-frequency k-mers.
    
    Fields:
        cuckoo_filter (CountingCuckooFilter): Tracks exact frequencies for low-frequency k-mers
        cms (CountMinSketch): Tracks approximate frequencies for high-frequency k-mers
        freq_threshold (int): Frequency threshold for moving k-mers from the Cuckoo Filter to the Count-Min Sketch
    """

    def __init__(self, cuckoo_bucket_size, cuckoo_num_buckets, cms_width, cms_depth, freq_threshold):
        """
        Initialize the EnhancedKMerTracker with specified parameters.
        
        Args:
            cuckoo_bucket_size (int): Maximum number of items each bucket in the Cuckoo Filter
            cuckoo_num_buckets (int): Total number of buckets in Cuckoo Filter
            cms_width (int): Number of buckets per hash function row in Count-Min Sketch
            cms_depth (int): Number of hash functions in Count-Min Sketch
            freq_threshold (int): Frequency threshold for moving k-mers to Count-Min Sketch
        """
        self.cuckoo_filter = CountingCuckooFilter(cuckoo_bucket_size, cuckoo_num_buckets)
        self.cms = CountMinSketch(cms_width, cms_depth)
        self.freq_threshold = freq_threshold

    def add_kmer(self, kmer):
        """
        Add a k-mer to the appropriate data structure or increment its frequency.

        Args:
            kmer (str): k-mer to add or increment.
        """
        # Check if the k-mer is already in the Count-Min Sketch
        if self.cms.estimate(kmer) > 0:
            # Increment frequency in the Count-Min Sketch
            self.cms.add(kmer)
            return
        
        # If not in the Count-Min Sketch check the Cuckoo Filter
        if not self.cuckoo_filter.lookup(kmer):
            # Insert the k-mer into the Cuckoo Filter
            self.cuckoo_filter.insert(kmer)
        
        elif self.cuckoo_filter.get_count(kmer) >= self.freq_threshold:
            # Migrate to Count-Min Sketch if above threshold
            count = self.cuckoo_filter.get_count(kmer)
            self.cms.add(kmer)
            self.cms.add(kmer, count - 1)
            # Remove k-mer from Cuckoo Filter
            self.cuckoo_filter.delete(kmer)
        
        else:
            # Increment the count in the Cuckoo Filter
            self.cuckoo_filter.insert(kmer)

    def get_kmer_frequency(self, kmer):
        """
        Get the frequency of a k-mer from the data structure

        Args:
            kmer (str): k-mer to retrieve the frequency for
        
        Returns:
            int: frequency of the k-mer or 0 if not found.
        """
        # Check if the k-mer is in Count-Min Sketch
        estimate = self.cms.estimate(kmer)
        if estimate > 0:
            return estimate
        
        # Otherwise check the Cuckoo Filter
        if self.cuckoo_filter.lookup(kmer):
            return self.cuckoo_filter.get_count(kmer)
        
        return 0
