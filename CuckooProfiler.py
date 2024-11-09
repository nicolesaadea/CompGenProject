import time
from memory_profiler import memory_usage
import random
from CuckooFilterImpl import CuckooFilter
import subprocess

# SeqIO is used to parse sequencing reads as this is an ifficient implementation and is not the focus of the project
from Bio import SeqIO

# Get all kmers for a sequence
def get_kmers(sequence, k):
    for i in range(len(sequence) - k + 1):
        yield sequence[i:i + k]


# Profiler for a specifically defined cuckoo filter
def profile_cuckoo_filter(cf, sequence, k):
    inserted_kmers = get_kmers(sequence, k)

    # Store starting time and memory
    start_time = time.time()
    initial_memory = memory_usage()[0]

    # Insert all kmers
    for kmer in inserted_kmers:
        cf.insert(kmer)

    # Measure time and memory taken to store kmers
    end_time = time.time()
    final_memory = memory_usage()[0]
    elapsed_time = end_time - start_time
    memory_used = final_memory - initial_memory
    print(f"Insertion Time: {elapsed_time:.4f} seconds")
    print(f"Memory Used for Insertion: {memory_used:.4f} MB")

    # Lookup a random subset of 5 kmers and test time
    lookup_kmers = random.sample(inserted_kmers, min(5, len(inserted_kmers)))
    start_time = time.time()
    for kmer in lookup_kmers:
        print(f"Lookup {kmer}: {cf.lookup(kmer)}")
    end_time = time.time()
    print(f"Lookup Time: {end_time - start_time:.4f} seconds")



# Generate a sequence with a set duplication rate of kmers
def generate_sequence_with_duplication(length, k, duplication_rate):
    unique_kmers = set()
    sequence = ""
    while len(sequence) < length:
        kmer = ''.join(random.choices('ACGT', k=k))

        # Add kmer if it is a duplicate and withing the duplicate limit
        if kmer in unique_kmers and random.random() < duplication_rate:
            sequence += kmer
        else:
            unique_kmers.add(kmer)
            sequence += kmer
    return sequence



# Test filters with a specified configuration over different duplication rates
def test_duplication_rates(sequence_length, k, duplication_rates, bucket_size, num_buckets):
    for rate in duplication_rates:
        print(f"\nTesting duplication rate: {rate}")
        sequence = generate_sequence_with_duplication(sequence_length, k, rate)
        cf = CuckooFilter(bucket_size=bucket_size, num_buckets=num_buckets)
        profile_cuckoo_filter(cf, sequence, k)



# Test filters with a specified configuration over different kmer sizes
def test_kmer_sizes(sequence_length, ks , duplication_rate, bucket_size, num_buckets, ):
    for k in ks:
        print(f"\nTesting k-mer size: {k}")
        sequence = generate_sequence_with_duplication(sequence_length, k, duplication_rate)
        cf = CuckooFilter(bucket_size=bucket_size, num_buckets=num_buckets)
        profile_cuckoo_filter(cf, sequence, k)



# Test filters with a variety of configurations
def test_filter_configurations(sequence_length, k, duplication_rate, configurations):
    for config in configurations:
        bucket_size, num_buckets, max_kicks = config
        print(f"Testing configuration: bucket_size={bucket_size}, num_buckets={num_buckets}, max_kicks={max_kicks}")
        sequence = generate_sequence_with_duplication(sequence_length, k, duplication_rate)
        cf = CuckooFilter(bucket_size=bucket_size, num_buckets=num_buckets, max_kicks=max_kicks)
        profile_cuckoo_filter(cf, sequence, k)



# Profile NGSReadsTreatment for comparison (MUST ADD DOWNLOAD NGSReadsTreatment AND ADD TO PATH)
def run_ngsreads_treatment(input_fastq, output_fastq):
    start_time = time.time()
    initial_memory = memory_usage()[0]
    
    command = f"NGSReadsTreatment -i {input_fastq} -o {output_fastq}"
    subprocess.run(command, shell=True)
    
    end_time = time.time()
    final_memory = memory_usage()[0]
    
    elapsed_time = end_time - start_time
    memory_used = final_memory - initial_memory
    
    return elapsed_time, memory_used

# Profile our Cuckoo filter for comparison
def run_cuckoo_filter_deduplication(input_fastq, output_fastq, bucket_size, num_buckets, k):
    cf = CuckooFilter(bucket_size=bucket_size, num_buckets=num_buckets)
    start_time = time.time()
    initial_memory = memory_usage()[0]
    
    with open(output_fastq, "w") as out_handle:
        for record in SeqIO.parse(input_fastq, "fastq"):
            kmer = str(record.seq[:k])
            if not cf.lookup(kmer):
                cf.insert(kmer)
                SeqIO.write(record, out_handle, "fastq")

    end_time = time.time()
    final_memory = memory_usage()[0]
    
    elapsed_time = end_time - start_time
    memory_used = final_memory - initial_memory
    
    return elapsed_time, memory_used

def count_total_reads(fastq_file):
    return sum(1 for _ in SeqIO.parse(fastq_file, "fastq"))

# Helper function to count unique reads in a FASTQ file
def count_unique_reads(fastq_file):
    unique_sequences = set()
    for record in SeqIO.parse(fastq_file, "fastq"):
        unique_sequences.add(str(record.seq))
    return len(unique_sequences)

# Compare all results for deduplication efficiency
def compare_deduplication_tools(input_fastq, output_fastq, bucket_size, num_buckets, k):
    results = []
    total_reads = count_total_reads(input_fastq)

    #FastUniq could also be tested

    # Run NGSReadsTreatment
    time_taken, memory_used = run_ngsreads_treatment(input_fastq, f"{output_fastq}_ngsreads.fastq")
    unique_reads_ngs = count_unique_reads(f"{output_fastq}_ngsreads.fastq")
    deduplication_percentage_ngs = (1 - unique_reads_ngs / total_reads) * 100
    results.append(("NGSReadsTreatment", time_taken, memory_used, deduplication_percentage_ngs))

    # Run Cuckoo Filter
    time_taken, memory_used = run_cuckoo_filter_deduplication(input_fastq, f"{output_fastq}_cuckoo.fastq", bucket_size, num_buckets, k)
    unique_reads_cuckoo = count_unique_reads(f"{output_fastq}_cuckoo.fastq")
    deduplication_percentage_cuckoo = (1 - unique_reads_cuckoo / total_reads) * 100
    results.append(("Cuckoo Filter", time_taken, memory_used, deduplication_percentage_cuckoo))

    # Print Comparisons
    print("\nComparison Results:")
    for tool, time_taken, memory_used, deduplication_percentage in results:
        print(f"{tool}: Time = {time_taken:.4f} seconds, Memory Used = {memory_used:.4f} MB, Deduplication Percentage = {deduplication_percentage:.2f}%")


def main():

    # Default Configurations
    sequence_length = 10000
    k = 21
    bucket_size = 4
    num_buckets = 50000
    
    duplication_rates = [0.1, 0.25, 0.5, 0.75]
    print("\nTesting Different Duplication Rates")
    test_duplication_rates(sequence_length, k, duplication_rates, bucket_size, num_buckets)

    kmer_sizes = [15, 21, 25]
    print("\nTesting Different k-mer Sizes")
    test_kmer_sizes(sequence_length, kmer_sizes, duplication_rate=0.5, bucket_size=bucket_size, num_buckets=num_buckets)

    configurations = [(4, 50000, 500), (4, 100000, 1000), (8, 50000, 500)]
    print("\nTesting Different Filter Configurations")
    test_filter_configurations(sequence_length, k, duplication_rate=0.5, configurations=configurations)


    input_fastq = "input_reads.fastq"
    output_fastq = "output"
    print("\nRunning Deduplication Comparison")
    compare_deduplication_tools(input_fastq, output_fastq, bucket_size, num_buckets, k)


if __name__ == "__main__":
    main()