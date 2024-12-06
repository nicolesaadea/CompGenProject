from CountingCuckooFilter import CountingCuckooFilter
from CountMinSketch import CountMinSketch
from collections import deque


class DynamicKMerTracker:
    """
    K-mer Frequency Tracker with dynamic thresholding. Combines Counting Cuckoo Filter and 
    Count-Min Sketch to prioritize accuracy for low-frequency k-mers and scalability for high-frequency k-mers.
    Dynamic thresholding adjusts the frequency threshold based on Cuckoo Filter use.
    
    Fields:
        cuckoo_filter (CountingCuckooFilter): Tracks exact frequencies for low-frequency k-mers
        cms (CountMinSketch): Tracks approximate frequencies for high-frequency k-mers
        freq_threshold (int): Frequency threshold for moving k-mers from Cuckoo Filter to Count-Min Sketch
        use_threshold (float): Maximum allowed use of the Cuckoo Filter before adjusting `freq_threshold`
    """

    def __init__(self, cuckoo_bucket_size, cuckoo_num_buckets, cms_width, cms_depth, initial_freq_threshold=50, use_threshold=0.75, smoothing_window=10):
        """
        Initialize the DynamicKMerTracker with specified parameters.
        
        Args:
            cuckoo_bucket_size (int): Maximum number of items each bucket in the Cuckoo Filter
            cuckoo_num_buckets (int): Total number of buckets in Cuckoo Filter
            cms_width (int): Number of buckets per hash function row in Count-Min Sketch
            cms_depth (int): Number of hash functions in Count-Min Sketch
            initial_freq_threshold (int): Initial frequency threshold for transitioning to Count-Min Sketch
            use_threshold (float): Max use rate of the Cuckoo Filter before dynamic adjustment
            smoothing_window (int): Number of recent use_rate values to average
        """
        self.cuckoo_filter = CountingCuckooFilter(cuckoo_bucket_size, cuckoo_num_buckets)
        self.cms = CountMinSketch(cms_width, cms_depth)
        self.freq_threshold = initial_freq_threshold
        self.use_threshold = use_threshold
        self.smoothing_window = smoothing_window
        self.use_rate_history = deque(maxlen=smoothing_window)

    def add_kmer(self, kmer):
        """
        Add a k-mer to the appropriate data structure or increment its frequency.

        Args:
            kmer (str): k-mer to add or increment.
        """
        # Adjust the frequency threshold dynamically based on use
        self.adjust_threshold()

        # Check if the k-mer is already in the Count-Min Sketch
        if self.cms.estimate(kmer) > 0:
            # Increment frequency in the Count-Min Sketch
            self.cms.add(kmer)
            return
        
        # If not in the Count-Min Sketch, check the Cuckoo Filter
        if not self.cuckoo_filter.lookup(kmer):
            # Insert the k-mer into the Cuckoo Filter
            self.cuckoo_filter.insert(kmer)
        
        elif self.cuckoo_filter.get_count(kmer) >= self.freq_threshold:
            # Migrate to Count-Min Sketch if frequency exceeds the threshold
            count = self.cuckoo_filter.get_count(kmer)
            self.cms.add(kmer, count)  # Add cumulative frequency to Count-Min Sketch
            self.cuckoo_filter.delete(kmer)  # Remove the k-mer from the Cuckoo Filter
        
        else:
            # Increment the count in the Cuckoo Filter
            self.cuckoo_filter.insert(kmer)

    def adjust_threshold(self):
        """
        Adjust the frequency threshold based on the smoothed use rate.
        """
        # Calculate current use rate
        current_use_rate = self.use_rate()
        self.use_rate_history.append(current_use_rate)  # Add to history
        
        # Compute running average of use rate
        smoothed_use_rate = sum(self.use_rate_history) / len(self.use_rate_history)

        # Adjust the frequency threshold based on average use rate
        if smoothed_use_rate > self.use_threshold:
            self.freq_threshold = max(1, self.freq_threshold - 1)
        elif smoothed_use_rate < self.use_threshold * 0.5:
            self.freq_threshold += 1

    def use_rate(self):
        """
        Calculate the use rate of the Cuckoo Filter.

        Returns:
            float: use rate as a fraction of the filter's total capacity.
        """
        total_capacity = self.cuckoo_filter.num_buckets * self.cuckoo_filter.bucket_size
        filled_slots = sum(len(bucket) for bucket in self.cuckoo_filter.buckets)
        return filled_slots / total_capacity

    def get_kmer_frequency(self, kmer):
        """
        Get the frequency of a k-mer from the data structure.

        Args:
            kmer (str): k-mer to retrieve the frequency for.
        
        Returns:
            int: Frequency of the k-mer or 0 if not found.
        """
        # Check if the k-mer is in Count-Min Sketch
        estimate = self.cms.estimate(kmer)
        if estimate > 0:
            return estimate
        
        # Otherwise check the Cuckoo Filter
        if self.cuckoo_filter.lookup(kmer):
            return self.cuckoo_filter.get_count(kmer)
        
        return 0
