import numpy as np
import hashlib

class CountMinSketch:
    """
    Count-Min Sketch for approximate frequency tracking.
    
    Fields:
        width (int): Number of buckets in each hash function row
        depth (int): Number of hash function rows
        table (2D array): 2D array for storing counts
        hash_functions (list): List of hash functions used for indexing items
    """

    def __init__(self, width, depth):
        """
        Initialize the Count-Min Sketch with specified width and depth.
        
        Args:
            width (int): Number of buckets per hash function row
            depth (int): Number of hash functions
        """
        self.width = width
        self.depth = depth
        self.table = np.zeros((depth, width), dtype=int)
        self.hash_functions = [self._hash_function(i) for i in range(depth)]

    def _hash_function(self, seed):
        """
        Generate a hash function based on a certain seed.
        
        Args:
            seed (int): Seed used to generate a hash function
        
        Returns:
            function: Hash function that computes an index for a given item.
        """
        def hash_fn(item):
            """
            Hash an item to produce an index in the Count-Min Sketch.
            
            Args:
                item (str): Item to be hashed
            
            Returns:
                int: Index to insert in
            """
            hash_value = hashlib.md5((str(seed) + item).encode()).hexdigest()
            return int(hash_value, 16) % self.width

        return hash_fn

    def add(self, item, count=1):
        """
        Increment the count for an item across all hash functions.
        
        Args:
            item (str): Item to be added to the Count-Min Sketch
            count (int): Frequency increment for the item (default 1)
        """
        for i in range(self.depth):
            self.table[i][self.hash_functions[i](item)] += count

    def estimate(self, item):
        """
        Estimate the frequency of an item in the Count-Min Sketch.
        
        Args:
            item (str): Item to estimate frequency for
        
        Returns:
            int: Estimated frequency of item (upper bound)
        """
        return min(self.table[i][self.hash_functions[i](item)] for i in range(self.depth))