import hashlib
import random

class CountingCuckooFilter:
    """
    A Counting Cuckoo Filter for approximate membership testing and frequency tracking.
    
    Attributes:
        bucket_size (int): Maximum number of items per bucket
        num_buckets (int): Total number of buckets
        max_kicks (int): Maximum number of rehash attempts before insertion fails
        buckets (list): List of buckets, each storing fingerprints and counts
    """

    def __init__(self, bucket_size, num_buckets, max_kicks=500):
        """
        Initialize the Counting Cuckoo Filter with specified parameters.
        
        Args:
            bucket_size (int): Maximum number of items per bucket
            num_buckets (int): Total number of buckets
            max_kicks (int): Maximum number of rehash attempts during insertion
        """
        self.bucket_size = bucket_size
        self.num_buckets = num_buckets
        self.max_kicks = max_kicks
        self.buckets = [{} for _ in range(num_buckets)]

    def _hash(self, item):
        """
        Generate a hash-based fingerprint for an item.
        
        Args:
            item (str): Item to hash
        
        Returns:
            int: 16-bit fingerprint of the item
        """
        return int(hashlib.sha256(item.encode()).hexdigest(), 16) % (2 ** 16)

    def _get_index(self, fingerprint, i):
        """
        Calculate the index of a bucket for a given fingerprint.
        
        Args:
            fingerprint (int): Item's fingerprint
            i (int): Offset (to get secondary index)
        
        Returns:
            int: Index of the bucket in the filter
        """
        return (fingerprint + i) % self.num_buckets

    def insert(self, item):
        """
        Insert an item into the Counting Cuckoo Filter, incrementing its count if it already exists.
        
        Args:
            item (str): Item to insert
        """
        fingerprint = self._hash(item)
        index1 = self._get_index(fingerprint, 0)
        index2 = self._get_index(fingerprint, index1)

        # Increment count in the first available bucket or insert if new item
        for index in [index1, index2]:
            if len(self.buckets[index]) < self.bucket_size or fingerprint in self.buckets[index]:
                self.buckets[index][fingerprint] = self.buckets[index].get(fingerprint, 0) + 1
                return True

        # If both buckets are full, evict an element and reinsert it (max_kicks attempts)
        index = random.choice([index1, index2])
        for _ in range(self.max_kicks):
            # Evict a random fingerprint
            evicted_fingerprint, count = random.choice(list(self.buckets[index].items()))
            del self.buckets[index][evicted_fingerprint]

             # Insert the new fingerprint
            self.buckets[index][fingerprint] = 1

            index = self._get_index(evicted_fingerprint, index)
            if len(self.buckets[index]) < self.bucket_size:
                # Reinsert the evicted fingerprint
                self.buckets[index][evicted_fingerprint] = count
                return

            fingerprint = evicted_fingerprint

    def lookup(self, item):
        """
        Check if an item exists in the Counting Cuckoo Filter.
        
        Args:
            item (str): Item to look up
        
        Returns:
            bool: True if the item is found
        """
        fingerprint = self._hash(item)
        index1 = self._get_index(fingerprint, 0)
        index2 = self._get_index(fingerprint, index1)
        return fingerprint in self.buckets[index1] or fingerprint in self.buckets[index2]

    def get_count(self, item):
        """
        Get the count of an item in the Counting Cuckoo Filter.
        
        Args:
            item (str): Item to retrieve the count for
        
        Returns:
            int: Count of the item or 0 if not found
        """
        fingerprint = self._hash(item)
        index1 = self._get_index(fingerprint, 0)
        index2 = self._get_index(fingerprint, index1)

        for index in [index1, index2]:
            if fingerprint in self.buckets[index]:
                return self.buckets[index][fingerprint]
        return 0

    def delete(self, item):
        """
        Remove an item from the Counting Cuckoo Filter
        
        Args:
            item (str): Item to remove
        """
        fingerprint = self._hash(item)
        index1 = self._get_index(fingerprint, 0)
        index2 = self._get_index(fingerprint, index1)

        for index in [index1, index2]:
            if fingerprint in self.buckets[index]:
                del self.buckets[index][fingerprint]
                return
