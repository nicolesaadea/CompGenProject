import hashlib
import random
import numpy as np

class CuckooFilterWithCMS:
    """
    Cuckoo Filter integrated with Count-Min Sketch for combined membership testing 
    and frequency estimation. Provides efficient item insertion, lookup, deletion, 
    and frequency tracking.

    Attributes:
        bucket_size (int): Maximum number of items per bucket in the Cuckoo Filter
        num_buckets (int): Total number of buckets in the Cuckoo Filter
        max_kicks (int): Maximum number of rehash attempts for Cuckoo Filter insertions
        buckets (list): Buckets for the Cuckoo Filter
        cms_width (int): Width of the Count-Min Sketch (number of hash buckets per row)
        cms_depth (int): Depth of the Count-Min Sketch (number of hash functions)
        cms_table (np.ndarray): Table for frequency estimation in the Count-Min Sketch
        hash_functions (list): List of hash functions for Count-Min Sketch
    """

    def __init__(self, bucket_size, num_buckets, cms_width, cms_depth, max_kicks=500):
        """
        Initialize the Cuckoo Filter and Count-Min Sketch with specified parameters.

        Args:
            bucket_size (int): Maximum number of items per bucket in the Cuckoo Filter
            num_buckets (int): Total number of buckets in the Cuckoo Filter
            cms_width (int): Width of the Count-Min Sketch
            cms_depth (int): Depth of the Count-Min Sketch
            max_kicks (int): Maximum number of rehash attempts for Cuckoo Filter
        """
        # Initialize Cuckoo Filter
        self.bucket_size = bucket_size
        self.num_buckets = num_buckets
        self.max_kicks = max_kicks
        self.buckets = [[] for _ in range(num_buckets)]

        # Initialize Count-Min Sketch
        self.cms_width = cms_width
        self.cms_depth = cms_depth
        self.cms_table = np.zeros((cms_depth, cms_width), dtype=int)
        self.hash_functions = [self._hash_function(i) for i in range(cms_depth)]

    def _hash(self, item):
        """
        Generate a 16-bit hash for an item.

        Args:
            item (str): Item to hash
        
        Returns:
            int: 16-bit fingerprint of the item
        """
        return int(hashlib.sha256(item.encode()).hexdigest(), 16) % (2 ** 16)

    def _get_index(self, fingerprint, i):
        """
        Compute the bucket index for a given fingerprint and offset.

        Args:
            fingerprint (int): Fingerprint of the item
            i (int): Offset for bucket computation
        
        Returns:
            int: Index of the bucket
        """
        return (fingerprint + i) % self.num_buckets

    def _hash_function(self, seed):
        """
        Generate a hash function for Count-Min Sketch based on a seed.

        Args:
            seed (int): Seed for the hash function
        
        Returns:
            function: Hash function for Count-Min Sketch
        """
        def hash_fn(item):
            hash_value = hashlib.md5((str(seed) + item).encode()).hexdigest()
            return int(hash_value, 16) % self.cms_width
        return hash_fn

    def insert(self, item):
        """
        Insert an item into the Cuckoo Filter and update its frequency in the Count-Min Sketch.

        Args:
            item (str): Item to insert
        
        Returns:
            bool: True if insertion succeeds, False otherwise
        """

        # Get indexes to try insertion
        fingerprint = self._hash(item)
        index1 = self._get_index(fingerprint, 0)
        index2 = self._get_index(fingerprint, index1)

        # Attempt insertion into the buckets
        if len(self.buckets[index1]) < self.bucket_size:
            self.buckets[index1].append(fingerprint)
        elif len(self.buckets[index2]) < self.bucket_size:
            self.buckets[index2].append(fingerprint)
        else:
            # If both buckets are full randomly evict an item and insert
            index = random.choice([index1, index2])
            for _ in range(self.max_kicks):
                evicted_fingerprint = self.buckets[index][random.randint(0, self.bucket_size - 1)]
                self.buckets[index].remove(evicted_fingerprint)
                self.buckets[index].append(fingerprint)

                # Recompute the bucket index for the evicted item and attempt insertion
                index = self._get_index(evicted_fingerprint, index)
                if len(self.buckets[index]) < self.bucket_size:
                    self.buckets[index].append(evicted_fingerprint)
                    break

                fingerprint = evicted_fingerprint
            else:
                return False

        # Update the frequency of the item in the Count-Min Sketch
        for i in range(self.cms_depth):
            index = self.hash_functions[i](item)
            self.cms_table[i][index] += 1

        return True


    def lookup(self, item):
        """
        Check if an item exists in the Cuckoo Filter and retrieve its frequency.

        Args:
            item (str): Item to lookup
        
        Returns:
            tuple: (bool, int or None) Presence in the Cuckoo Filter and frequency from Count-Min Sketch
        """
        fingerprint = self._hash(item)
        index1 = self._get_index(fingerprint, 0)
        index2 = self._get_index(fingerprint, index1)

        if fingerprint in self.buckets[index1] or fingerprint in self.buckets[index2]:
            frequency = min(self.cms_table[i][self.hash_functions[i](item)] for i in range(self.cms_depth))
            return True, frequency
        return False, None

    def delete(self, item):
        """
        Delete an item from the Cuckoo Filter and update the Count-Min Sketch.

        Args:
            item (str): Item to delete
        
        Returns:
            bool: True if the item was deleted, False otherwise
        """
        fingerprint = self._hash(item)
        index1 = self._get_index(fingerprint, 0)
        index2 = self._get_index(fingerprint, index1)

        if fingerprint in self.buckets[index1]:
            self.buckets[index1].remove(fingerprint)
        elif fingerprint in self.buckets[index2]:
            self.buckets[index2].remove(fingerprint)
        else:
            return False

        # Update Count-Min Sketch
        for i in range(self.cms_depth):
            index = self.hash_functions[i](item)
            self.cms_table[i][index] = max(0, self.cms_table[i][index] - 1)

        return True

    def get_estimated_frequency(self, item):
        """
        Retrieve the estimated frequency of an item from the Count-Min Sketch.

        Args:
            item (str): Item to query
        
        Returns:
            int: Estimated frequency of the item
        """
        return min(self.cms_table[i][self.hash_functions[i](item)] for i in range(self.cms_depth))