import time
from memory_profiler import memory_usage
import matplotlib.pyplot as plt
import tracemalloc
import os
import psutil
import time
import subprocess
import random
import gc
from Bio import SeqIO
from CuckooFilterImpl import CuckooFilter
import subprocess

# def profile_memory(func, *args, **kwargs):
#     """Profile net change in memory usage of a function."""
#     gc.collect()
#     # Get memory usage before function execution
#     start_mem = memory_usage()[0]
#     func(*args, **kwargs)
#     gc.collect()
#     # Get memory usage after function execution
#     end_mem = memory_usage()[0]
#     # Calculate net change in memory usage
#     mem_diff = end_mem - start_mem
#     return mem_diff

def get_kmers(sequence, k):
    """Get all kmers for a sequence"""
    for i in range(len(sequence) - k + 1):
        yield sequence[i:i + k]

import psutil
import time

def profile_cuckoo_filter(cf, sequence, k):
    """Profile a cuckoo filter's performance using psutil."""
    inserted_kmers = list(get_kmers(sequence, k))
    
    process = psutil.Process(os.getpid())
    max_memory = process.memory_info().rss / (1024 * 1024)  # Initial memory in MB
    
    # Time and memory for insertion
    start_time = time.perf_counter()
    for kmer in inserted_kmers:
        cf.insert(kmer)
        mem_usage_mb = process.memory_info().rss / (1024 * 1024)
        if mem_usage_mb > max_memory:
            max_memory = mem_usage_mb
    insertion_time = time.perf_counter() - start_time
    insertion_memory = max_memory
    
    # Lookup operations (as before)
    lookup_kmers_present = random.sample(inserted_kmers, min(1000, len(inserted_kmers)))
    lookup_kmers_absent = [''.join(random.choices('ACGT', k=k)) for _ in range(1000)]
    total_lookups = lookup_kmers_present + lookup_kmers_absent
    random.shuffle(total_lookups)
    
    start_time = time.perf_counter()
    lookup_results = [cf.lookup(kmer) for kmer in total_lookups]
    lookup_time = time.perf_counter() - start_time
    
    true_positives = sum(1 for i in range(len(lookup_kmers_present)) if lookup_results[i])
    false_positives = sum(1 for i in range(len(lookup_kmers_present), len(total_lookups)) if lookup_results[i])
    
    tpr = true_positives / len(lookup_kmers_present)
    fpr = false_positives / len(lookup_kmers_absent)
    
    print(f"Insertion Time: {insertion_time:.4f} seconds")
    print(f"Peak Memory Used for Insertion: {insertion_memory:.4f} MB")
    print(f"Lookup Time for {len(total_lookups)} k-mers: {lookup_time:.4f} seconds")
    print(f"True Positive Rate: {tpr:.4f}")
    print(f"False Positive Rate: {fpr:.4f}")
    
    return {
        "insertion_time": insertion_time,
        "insertion_memory": insertion_memory,
        "lookup_time": lookup_time,
        "true_positive_rate": tpr,
        "false_positive_rate": fpr
    }

    
    
    
def generate_synthetic_fastq(filename, num_reads, read_length, duplication_rate):
    """Generate a synthetic FASTQ file with specified parameters"""
    sequences = []
    num_unique = int(num_reads * (1 - duplication_rate))
    num_duplicates = num_reads - num_unique

    # Generate unique reads
    for _ in range(num_unique):
        seq = ''.join(random.choices('ACGT', k=read_length))
        sequences.append(seq)

    # Generate duplicate reads
    duplicates = random.choices(sequences, k=num_duplicates)
    sequences.extend(duplicates)
    random.shuffle(sequences)

    # Write to FASTQ file
    with open(filename, "w") as fh:
        for i, seq in enumerate(sequences):
            fh.write(f"@read_{i}\n{seq}\n+\n{'~'*read_length}\n")

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


def run_ngsreads_treatment(input_fastq):
    """Run NGSReadsTreatment and profile its performance using psutil."""
    command = f"java -jar NGSReadsTreatment_v1.3.jar {input_fastq} 2"
    
    # Start the subprocess
    process = subprocess.Popen(command, shell=True)
    ps_process = psutil.Process(process.pid)
    
    start_time = time.time()
    max_memory = 0
    
    # Monitor the process until it finishes
    while process.poll() is None:
        try:
            # Get the memory usage of the subprocess
            mem_info = ps_process.memory_info()
            memory_usage_mb = mem_info.rss / (1024 * 1024)  # Convert to MB
            if memory_usage_mb > max_memory:
                max_memory = memory_usage_mb
        except psutil.NoSuchProcess:
            # The process might have already finished
            break
        time.sleep(0.1)  # Poll every 100 ms
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    
    return elapsed_time, max_memory



def run_cuckoo_filter_deduplication(input_fastq, output_fastq, bucket_size, num_buckets, k):
    """Run Cuckoo filter deduplication and profile its performance using psutil."""
    gc.collect()
    process = psutil.Process(os.getpid())
    baseline_memory = process.memory_info().rss / (1024 * 1024)  # Initial memory in MB
    max_memory = baseline_memory
    last_sample_time = time.time()
    sample_interval = 0.1  # Sample every 0.1 seconds

    cf = CuckooFilter(bucket_size=bucket_size, num_buckets=num_buckets)

    def deduplicate_reads():
        nonlocal max_memory, last_sample_time  # Declare variables as nonlocal
        with open(output_fastq, "w") as out_handle:
            for idx, record in enumerate(SeqIO.parse(input_fastq, "fastq")):
                kmer = str(record.seq[:k])
                if not cf.lookup(kmer):
                    cf.insert(kmer)
                    SeqIO.write(record, out_handle, "fastq")
                # Monitor memory usage every 0.1 seconds
                current_time = time.time()
                if current_time - last_sample_time >= sample_interval:
                    mem_usage_mb = process.memory_info().rss / (1024 * 1024)
                    max_memory = max(max_memory, mem_usage_mb)
                    last_sample_time = current_time

    start_time = time.time()
    deduplicate_reads()
    elapsed_time = time.time() - start_time

    net_memory_usage = max_memory - baseline_memory
    return elapsed_time, net_memory_usage


def count_reads(fastq_file, unique=False):
    """Count total or unique reads in a FASTQ file"""
    if unique:
        unique_sequences = set()
        for record in SeqIO.parse(fastq_file, "fastq"):
            unique_sequences.add(str(record.seq))
        return len(unique_sequences)
    return sum(1 for _ in SeqIO.parse(fastq_file, "fastq"))

def compare_deduplication_tools(input_fastq, output_fastq, bucket_size, num_buckets, k):
    """Compare different deduplication tools using consistent memory profiling."""
    results = []
    total_reads = count_reads(input_fastq)

    # Test NGSReadsTreatment
    time_taken, memory_used = run_ngsreads_treatment(input_fastq)
    new_input_fastq = input_fastq.split(".")[0]
    unique_reads_ngs = count_reads(f"{new_input_fastq}_1_trated.fastq", unique=True)
    dedup_percentage_ngs = (1 - unique_reads_ngs / total_reads) * 100
    results.append(("NGSReadsTreatment", time_taken, memory_used, dedup_percentage_ngs))

    # Test Cuckoo Filter Deduplication
    time_taken, memory_used = run_cuckoo_filter_deduplication(
        input_fastq, f"{output_fastq}_cuckoo.fastq",
        bucket_size, num_buckets, k
    )
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
    duplication_rates = [0.1, 0.25, 0.5, 0.75]
    tpr_results = []
    fpr_results = []
    insertion_times = []
    lookup_times = []
    memory_usages = []
    
    print("\nTesting Different Duplication Rates")
    for rate in duplication_rates:
        print(f"\nDuplication rate: {rate}")
        sequence = generate_sequence_with_duplication(sequence_length, k, rate)
        cf = CuckooFilter(bucket_size=bucket_size, num_buckets=num_buckets)
        metrics = profile_cuckoo_filter(cf, sequence, k)
        
        tpr_results.append(metrics["true_positive_rate"])
        fpr_results.append(metrics["false_positive_rate"])
        insertion_times.append(metrics["insertion_time"])
        lookup_times.append(metrics["lookup_time"])
        memory_usages.append(metrics["insertion_memory"])
        
    # Visualize results
    plt.figure(figsize=(12, 8))
    plt.plot(duplication_rates, tpr_results, label='True Positive Rate', marker='o')
    plt.plot(duplication_rates, fpr_results, label='False Positive Rate', marker='o')
    plt.xlabel('Duplication Rate')
    plt.ylabel('Rate')
    plt.title('TPR and FPR vs. Duplication Rate')
    plt.legend()
    plt.grid(True)
    plt.show()
    
    # Visualization of Insertion and Lookup Times vs Duplication Rate
    plt.figure(figsize=(10, 6))
    plt.plot(duplication_rates, insertion_times, label='Insertion Time', marker='o')
    plt.plot(duplication_rates, lookup_times, label='Lookup Time', marker='o')
    plt.xlabel('Duplication Rate')
    plt.ylabel('Time (seconds)')
    plt.title('Insertion and Lookup Times vs. Duplication Rate')
    plt.legend()
    plt.grid(True)
    plt.show()

    # Visualization of Memory Usage vs Duplication Rate
    plt.figure(figsize=(10, 6))
    plt.plot(duplication_rates, memory_usages, label='Memory Usage', marker='o')
    plt.xlabel('Duplication Rate')
    plt.ylabel('Memory Usage (MB)')
    plt.title('Memory Usage vs. Duplication Rate')
    plt.legend()
    plt.grid(True)
    plt.show()

    kmer_sizes = [15, 21, 25]
    tpr_results = []
    fpr_results = []
    insertion_times = []
    lookup_times = []
    memory_usages = []
    # Test k-mer sizes
    print("\nTesting Different k-mer Sizes")
    for ksize in kmer_sizes:
        print(f"\nk-mer size: {ksize}")
        sequence = generate_sequence_with_duplication(sequence_length, ksize, 0.5)
        cf = CuckooFilter(bucket_size=bucket_size, num_buckets=num_buckets)
        metrics = profile_cuckoo_filter(cf, sequence, ksize)
        
        tpr_results.append(metrics["true_positive_rate"])
        fpr_results.append(metrics["false_positive_rate"])
        insertion_times.append(metrics["insertion_time"])
        lookup_times.append(metrics["lookup_time"])
        memory_usages.append(metrics["insertion_memory"])
        
    # Visualize results for k-mer sizes
    plt.figure(figsize=(10, 6))
    plt.plot(kmer_sizes, tpr_results, label='True Positive Rate', marker='o')
    plt.plot(kmer_sizes, tpr_results, label='False Positive Rate', marker='o')
    plt.xlabel('k-mer Size')
    plt.ylabel('Rate')
    plt.title('TPR and FPR vs. k-mer Size')
    plt.legend()
    plt.grid(True)
    plt.show()
        
    # Visualization of Insertion and Lookup Times vs k-mer Size
    plt.figure(figsize=(10, 6))
    plt.plot(kmer_sizes, insertion_times, label='Insertion Time', marker='o')
    plt.plot(kmer_sizes, lookup_times, label='Lookup Time', marker='o')
    plt.xlabel('k-mer Size')
    plt.ylabel('Time (seconds)')
    plt.title('Insertion and Lookup Times vs. k-mer Size')
    plt.legend()
    plt.grid(True)
    plt.show()
    
    # Visualization of Memory Usage vs k-mer Size
    plt.figure(figsize=(10, 6))
    plt.plot(kmer_sizes, memory_usages, label='Memory Usage', marker='o')
    plt.xlabel('k-mer Size')
    plt.ylabel('Memory Usage (MB)')
    plt.title('Memory Usage vs. k-mer Size')
    plt.legend()
    plt.grid(True)
    plt.show()
    
        
    configs = [
        {'bucket_size': 4, 'num_buckets': 50000, 'max_kicks': 500},
        {'bucket_size': 4, 'num_buckets': 100000, 'max_kicks': 1000},
        {'bucket_size': 8, 'num_buckets': 50000, 'max_kicks': 500}
    ]
    config_labels = []
    # Test filter configurations
    print("\nTesting Different Filter Configurations")
    for config in configs:
        b_size = config['bucket_size']
        b_num = config['num_buckets']
        max_kicks = config['max_kicks']
        label = f"b_size={b_size}, b_num={b_num}, max_kicks={max_kicks}"
        config_labels.append(label)
        print(f"\nConfiguration: {label}")
        sequence = generate_sequence_with_duplication(sequence_length, k, 0.5)
        cf = CuckooFilter(bucket_size=b_size, num_buckets=b_num, max_kicks=max_kicks)
        profile_cuckoo_filter(cf, sequence, k)

def main():
    """Main execution function"""
    random.seed(42)
    test_configurations()

    # Run deduplication comparison
    input_fastq = "input_reads.fastq"
    if not os.path.exists(input_fastq):
        generate_synthetic_fastq(input_fastq, num_reads=10000, read_length=100, duplication_rate=0.5)
    output_fastq = "output"
    print("\nRunning Deduplication Comparison")
    compare_deduplication_tools(input_fastq, output_fastq, 
                              bucket_size=4, num_buckets=50000, k=21)

if __name__ == "__main__":
    main()