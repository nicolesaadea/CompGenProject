import hashlib
import random
import numpy as np

class CuckooFilterWithCMS:
    def __init__(self, bucket_size, num_buckets, cms_width, cms_depth, max_kicks=500):
        # Cuckoo filter settings
        self.bucket_size = bucket_size
        self.num_buckets = num_buckets
        self.max_kicks = max_kicks
        self.buckets = [[] for _ in range(num_buckets)]
        
        # Count-Min Sketch settings
        self.cms_width = cms_width
        self.cms_depth = cms_depth
        self.cms_table = np.zeros((cms_depth, cms_width), dtype=int)
        self.hash_functions = [self._hash_function(i) for i in range(cms_depth)]

    def _hash(self, item):
        #Return a hash of the item.
        return int(hashlib.sha256(item.encode()).hexdigest(), 16) % (2 ** 16)

    def _get_index(self, fingerprint, i):
        #Get the index in the bucket based on the fingerprint and an index i.
        return (fingerprint + i) % self.num_buckets

    def _hash_function(self, seed):
        #Generate a hash function for Count-Min Sketch based on a seed.
        def hash_fn(item):
            hash_value = hashlib.md5((str(seed) + item).encode()).hexdigest()
            return int(hash_value, 16) % self.cms_width
        return hash_fn

    def insert(self, item):
        #Insert an item into the Cuckoo Filter and update its frequency in the CMS.
        # Cuckoo filter insert logic
        fingerprint = self._hash(item)
        index1 = self._get_index(fingerprint, 0)
        index2 = self._get_index(fingerprint, index1)

        # Try to insert into the first available index
        if len(self.buckets[index1]) < self.bucket_size:
            self.buckets[index1].append(fingerprint)
        elif len(self.buckets[index2]) < self.bucket_size:
            self.buckets[index2].append(fingerprint)
        else:
            # Cuckoo filtering with kicking out an existing item if both are full
            index = random.choice([index1, index2])
            for _ in range(self.max_kicks):
                evicted_fingerprint = self.buckets[index][random.randint(0, self.bucket_size - 1)]
                self.buckets[index].remove(evicted_fingerprint)
                self.buckets[index].append(fingerprint)

                index = self._get_index(evicted_fingerprint, index)
                if len(self.buckets[index]) < self.bucket_size:
                    self.buckets[index].append(evicted_fingerprint)
                    break

                fingerprint = evicted_fingerprint
            else:
                return False  # Failed to insert after max kicks

        # Update the frequency of the item in the CMS
        for i in range(self.cms_depth):
            index = self.hash_functions[i](item)
            self.cms_table[i][index] += 1

        return True

    def lookup(self, item):
        #Check if an item exists in the Cuckoo Filter
        fingerprint = self._hash(item)
        index1 = self._get_index(fingerprint, 0)
        index2 = self._get_index(fingerprint, index1)

        if fingerprint in self.buckets[index1] or fingerprint in self.buckets[index2]:
            # Get frequency estimate from the CMS
            min_frequency = min(self.cms_table[i][self.hash_functions[i](item)] for i in range(self.cms_depth))
            return True, min_frequency
        return False, None

    def delete(self, item):
        #Delete an item from the Cuckoo Filter
        fingerprint = self._hash(item)
        index1 = self._get_index(fingerprint, 0)
        index2 = self._get_index(fingerprint, index1)

        if fingerprint in self.buckets[index1]:
            self.buckets[index1].remove(fingerprint)
            # Reduce the count in CMS
            for i in range(self.cms_depth):
                index = self.hash_functions[i](item)
                self.cms_table[i][index] = max(0, self.cms_table[i][index] - 1)
            return True
        elif fingerprint in self.buckets[index2]:
            self.buckets[index2].remove(fingerprint)
            # Reduce the count in CMS
            for i in range(self.cms_depth):
                index = self.hash_functions[i](item)
                self.cms_table[i][index] = max(0, self.cms_table[i][index] - 1)
            return True
        return False

    def get_estimated_frequency(self, item):
        #Return the frequency estimate for an item from the CMS
        min_frequency = min(self.cms_table[i][self.hash_functions[i](item)] for i in range(self.cms_depth))
        return min_frequency
