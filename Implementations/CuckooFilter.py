import hashlib
import random

class CuckooFilter:
    def __init__(self, bucket_size, num_buckets, max_kicks=500):
        self.bucket_size = bucket_size
        self.num_buckets = num_buckets
        self.max_kicks = max_kicks
        self.buckets = [[] for _ in range(num_buckets)]

    def _hash(self, item):
        return int(hashlib.sha256(item.encode()).hexdigest(), 16) % (2 ** 16)

    def _get_index(self, fingerprint, i):
        return (fingerprint + i) % self.num_buckets

    def insert(self, item):
        fingerprint = self._hash(item)
        index1 = self._get_index(fingerprint, 0)
        index2 = self._get_index(fingerprint, index1)

        if len(self.buckets[index1]) < self.bucket_size:
            self.buckets[index1].append(fingerprint)
            return True
        elif len(self.buckets[index2]) < self.bucket_size:
            self.buckets[index2].append(fingerprint)
            return True

        index = random.choice([index1, index2])
        for _ in range(self.max_kicks):
            evicted_fingerprint = self.buckets[index][random.randint(0, self.bucket_size - 1)]
            self.buckets[index].remove(evicted_fingerprint)
            self.buckets[index].append(fingerprint)

            index = self._get_index(evicted_fingerprint, index)
            if len(self.buckets[index]) < self.bucket_size:
                self.buckets[index].append(evicted_fingerprint)
                return True

            fingerprint = evicted_fingerprint

        return False

    def lookup(self, item):
        fingerprint = self._hash(item)
        index1 = self._get_index(fingerprint, 0)
        index2 = self._get_index(fingerprint, index1)

        return fingerprint in self.buckets[index1] or fingerprint in self.buckets[index2]

    def delete(self, item):
        fingerprint = self._hash(item)
        index1 = self._get_index(fingerprint, 0)
        index2 = self._get_index(fingerprint, index1)

        if fingerprint in self.buckets[index1]:
            self.buckets[index1].remove(fingerprint)
            return True
        elif fingerprint in self.buckets[index2]:
            self.buckets[index2].remove(fingerprint)
            return True
        return False
