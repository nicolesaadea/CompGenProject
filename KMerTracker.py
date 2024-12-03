from CuckooFilterImpl import CuckooFilter
from CountMinSketchImp import CountMinSketch

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