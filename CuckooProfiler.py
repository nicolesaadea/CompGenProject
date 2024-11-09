import time
from memory_profiler import memory_usage
import random
from Bio import SeqIO
from CuckooFilterImpl import CuckooFilter
import subprocess

def profile_memory(func, *args, **kwargs):
    """Profile memory usage of a function"""
    mem_usage = memory_usage((func, args, kwargs), interval=0.1)
    return max(mem_usage) - min(mem_usage)

def get_kmers(sequence, k):
    """Get all kmers for a sequence"""
    for i in range(len(sequence) - k + 1):
        yield sequence[i:i + k]

def profile_cuckoo_filter(cf, sequence, k):
    """Profile a cuckoo filter's performance"""
    inserted_kmers = list(get_kmers(sequence, k))

    def insert_all_kmers():
        for kmer in inserted_kmers:
            cf.insert(kmer)

    start_time = time.time()
    memory_used = profile_memory(insert_all_kmers)
    elapsed_time = time.time() - start_time

    print(f"Insertion Time: {elapsed_time:.4f} seconds")
    print(f"Memory Used for Insertion: {memory_used:.4f} MB")

    # Lookup testing
    lookup_kmers = random.sample(inserted_kmers, min(5, len(inserted_kmers)))
    start_time = time.time()
    for kmer in lookup_kmers:
        print(f"Lookup {kmer}: {cf.lookup(kmer)}")
    end_time = time.time()
    print(f"Lookup Time: {end_time - start_time:.4f} seconds")

def generate_sequence_with_duplication(length, k, duplication_rate):
    """Generate a sequence with controlled duplication rate"""
    num_kmers = length // k
    num_unique = int(num_kmers * (1 - duplication_rate))
    num_duplicates = num_kmers - num_unique

    unique_kmers = [''.join(random.choices('ACGT', k=k)) for _ in range(num_unique)]
    duplicate_kmers = random.choices(unique_kmers, k=num_duplicates)
    all_kmers = unique_kmers + duplicate_kmers
    random.shuffle(all_kmers)
    return ''.join(all_kmers)

def run_ngsreads_treatment(input_fastq, output_fastq):
    """Run NGSReadsTreatment and profile its performance"""
    start_time = time.time()
    initial_memory = memory_usage()[0]
    
    command = f"NGSReadsTreatment -i {input_fastq} -o {output_fastq}"
    subprocess.run(command, shell=True)
    
    end_time = time.time()
    final_memory = memory_usage()[0]
    
    return end_time - start_time, final_memory - initial_memory

def run_cuckoo_filter_deduplication(input_fastq, output_fastq, bucket_size, num_buckets, k):
    """Run Cuckoo filter deduplication and profile its performance"""
    cf = CuckooFilter(bucket_size=bucket_size, num_buckets=num_buckets)
    
    def deduplicate_reads():
        seen_kmers = set()  # For verification
        with open(output_fastq, "w") as out_handle:
            for record in SeqIO.parse(input_fastq, "fastq"):
                kmer = str(record.seq[:k])
                if not cf.lookup(kmer):
                    cf.insert(kmer)
                    seen_kmers.add(kmer)  # For verification
                    SeqIO.write(record, out_handle, "fastq")
        return len(seen_kmers)

    start_time = time.time()
    # Properly scope the memory profiling
    memory_used = profile_memory(deduplicate_reads)
    elapsed_time = time.time() - start_time

    return elapsed_time, memory_used

def count_reads(fastq_file, unique=False):
    """Count total or unique reads in a FASTQ file"""
    if unique:
        return len(set(str(record.seq) for record in SeqIO.parse(fastq_file, "fastq")))
    return sum(1 for _ in SeqIO.parse(fastq_file, "fastq"))

def compare_deduplication_tools(input_fastq, output_fastq, bucket_size, num_buckets, k):
    """Compare different deduplication tools"""
    results = []
    total_reads = count_reads(input_fastq)

    # Test NGSReadsTreatment
    time_taken, memory_used = run_ngsreads_treatment(input_fastq, f"{output_fastq}_ngsreads.fastq")
    unique_reads_ngs = count_reads(f"{output_fastq}_ngsreads.fastq", unique=True)
    dedup_percentage_ngs = (1 - unique_reads_ngs / total_reads) * 100
    results.append(("NGSReadsTreatment", time_taken, memory_used, dedup_percentage_ngs))


    time_taken, memory_used = run_cuckoo_filter_deduplication(input_fastq, f"{output_fastq}_cuckoo.fastq", 
                                                            bucket_size, num_buckets, k)
    unique_reads_cuckoo = count_reads(f"{output_fastq}_cuckoo.fastq", unique=True)
    dedup_percentage_cuckoo = (1 - unique_reads_cuckoo / total_reads) * 100
    results.append(("Cuckoo Filter", time_taken, memory_used, dedup_percentage_cuckoo))


    print("\nComparison Results:")
    for tool, time_taken, memory_used, dedup_percentage in results:
        print(f"{tool}:")
        print(f"  Time: {time_taken:.4f} seconds")
        print(f"  Memory: {memory_used:.4f} MB")
        print(f"  Deduplication: {dedup_percentage:.2f}%")

def test_configurations(sequence_length=10000, k=21, bucket_size=4, num_buckets=50000):
    """Run various configuration tests"""
    # Test duplication rates
    print("\nTesting Different Duplication Rates")
    for rate in [0.1, 0.25, 0.5, 0.75]:
        print(f"\nDuplication rate: {rate}")
        sequence = generate_sequence_with_duplication(sequence_length, k, rate)
        cf = CuckooFilter(bucket_size=bucket_size, num_buckets=num_buckets)
        profile_cuckoo_filter(cf, sequence, k)

    # Test k-mer sizes
    print("\nTesting Different k-mer Sizes")
    for ksize in [15, 21, 25]:
        print(f"\nk-mer size: {ksize}")
        sequence = generate_sequence_with_duplication(sequence_length, ksize, 0.5)
        cf = CuckooFilter(bucket_size=bucket_size, num_buckets=num_buckets)
        profile_cuckoo_filter(cf, sequence, ksize)

    # Test filter configurations
    print("\nTesting Different Filter Configurations")
    configs = [(4, 50000, 500), (4, 100000, 1000), (8, 50000, 500)]
    for b_size, b_num, max_kicks in configs:
        print(f"\nConfiguration: bucket_size={b_size}, num_buckets={b_num}, max_kicks={max_kicks}")
        sequence = generate_sequence_with_duplication(sequence_length, k, 0.5)
        cf = CuckooFilter(bucket_size=b_size, num_buckets=b_num, max_kicks=max_kicks)
        profile_cuckoo_filter(cf, sequence, k)

def main():
    """Main execution function"""
    test_configurations()

    # Run deduplication comparison
    input_fastq = "input_reads.fastq"
    output_fastq = "output"
    print("\nRunning Deduplication Comparison")
    compare_deduplication_tools(input_fastq, output_fastq, 
                              bucket_size=4, num_buckets=50000, k=21)

if __name__ == "__main__":
    main()