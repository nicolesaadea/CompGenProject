import numpy as np
import hashlib

from CuckooFilterImpl import CuckooFilter

class CountMinSketch:
    def __init__(self, width, depth):
        # Initialize the Count-Min Sketch with specified width and depth
        self.width = width  
        self.depth = depth  
        self.table = np.zeros((depth, width), dtype=int)  # Create a 2D count table
        self.hash_functions = [self.hash_function(i) for i in range(depth)]  # Generate hash functions

    def hash_function(self, seed):
        # Return a hash function that produces an index for the given item based on the seed
        def hash_fn(item):
            return (int(hashlib.md5((str(seed) + item).encode()).hexdigest(), 16) % self.width)
        return hash_fn

    def add(self, item):
        # Increment the count for the item in the Count-Min Sketch for each hash function
        for i in range(self.depth):
            self.table[i][self.hash_functions[i](item)] += 1

    def estimate(self, item):
        # Return the minimum count estimate for the item across all hash functions
        return min(self.table[i][self.hash_functions[i](item)] for i in range(self.depth))
    
class KMerTracker:
    def __init__(self, cuckoo_bucket_size, cuckoo_num_buckets, cms_width, cms_depth):
        # Initialize KMerTracker with Cuckoo Filter and Count-Min Sketch instances
        self.cuckoo_filter = CuckooFilter(cuckoo_bucket_size, cuckoo_num_buckets)
        self.cms = CountMinSketch(cms_width, cms_depth)

    def add_kmer(self, kmer):
        # Insert the k-mer into the Cuckoo Filter if it's not a duplicate and add it to Count-Min Sketch
        if not self.cuckoo_filter.lookup(kmer):  
            self.cuckoo_filter.insert(kmer)       
        self.cms.add(kmer)                       

    def get_count(self, kmer):
        # Return the estimated frequency of the k-mer from the Count-Min Sketch
        return self.cms.estimate(kmer)