# profile_cuckoo.py

import sys
import os
import time
import psutil
import gc
import random
import matplotlib.pyplot as plt
from Bio import SeqIO
from pathlib import Path

# Add the parent directory to sys.path (if necessary for imports)
current_dir = Path(__file__).resolve().parent
parent_dir = current_dir.parent
sys.path.insert(0, str(parent_dir))

from Implementations.CuckooFilter import CuckooFilter

def generate_synthetic_fastq(filename, num_reads, read_length, duplication_rate):
    """Generate a synthetic FASTQ file with specified parameters."""
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
    print(f"Generated synthetic FASTQ file: {filename}")

def run_cuckoo_filter_deduplication(input_fastq, output_fastq, bucket_size, num_buckets, k):
    """
    Run Cuckoo filter deduplication and profile its performance using psutil.

    Args:
        input_fastq (Path): Path to the input FASTQ file.
        output_fastq (Path): Path to the output deduplicated FASTQ file.
        bucket_size (int): Bucket size parameter for the Cuckoo Filter.
        num_buckets (int): Number of buckets for the Cuckoo Filter.
        k (int): k-mer size for deduplication.

    Returns:
        dict: Dictionary containing profiling metrics.
    """
    gc.collect()
    process = psutil.Process(os.getpid())
    baseline_memory = process.memory_info().rss / (1024 * 1024)  # Initial memory in MB
    max_memory = baseline_memory

    cf = CuckooFilter(bucket_size=bucket_size, num_buckets=num_buckets)

    total_reads = 0
    unique_reads = 0

    start_time = time.time()
    max_memory_during_insertion = baseline_memory

    # Open input and output FASTQ files
    with open(input_fastq, "r") as infile, open(output_fastq, "w") as outfile:
        for record in SeqIO.parse(infile, "fastq"):
            total_reads += 1
            kmer = str(record.seq[:k])

            if not cf.lookup(kmer):
                cf.insert(kmer)
                SeqIO.write(record, outfile, "fastq")
                unique_reads += 1

            # Update memory usage
            current_memory = process.memory_info().rss / (1024 * 1024)  # in MB
            if current_memory > max_memory_during_insertion:
                max_memory_during_insertion = current_memory

    end_time = time.time()
    elapsed_time = end_time - start_time
    peak_memory = max_memory_during_insertion

    dedup_percentage = (1 - unique_reads / total_reads) * 100 if total_reads else 0.0

    return {
        "time_taken": elapsed_time,
        "peak_memory": peak_memory,
        "dedup_percentage": dedup_percentage
    }

def count_reads(fastq_file, unique=False):
    """Count total or unique reads in a FASTQ file."""
    if unique:
        unique_sequences = set()
        for record in SeqIO.parse(fastq_file, "fastq"):
            unique_sequences.add(str(record.seq))
        return len(unique_sequences)
    return sum(1 for _ in SeqIO.parse(fastq_file, "fastq"))

def compare_cuckoo(input_fastq, output_fastq, bucket_size, num_buckets, k):
    """Compare Cuckoo Filter deduplication."""
    print(f"Counting total reads in {input_fastq}...")
    total_reads = count_reads(input_fastq)

    print(f"Running Cuckoo Filter deduplication on {input_fastq}...")
    metrics = run_cuckoo_filter_deduplication(
        input_fastq, output_fastq,
        bucket_size, num_buckets, k
    )

    print("\nCuckoo Filter Profiling Results:")
    print(f"Time Taken: {metrics['time_taken']:.4f} seconds")
    print(f"Peak Memory Used: {metrics['peak_memory']:.4f} MB")
    print(f"Deduplication Percentage: {metrics['dedup_percentage']:.2f}%")

def main():
    """Main execution function for Cuckoo Filter profiling."""
    random.seed(42)

    # Define duplication rates and corresponding filenames
    duplication_rates = [0.25, 0.5, 0.75, 0.9]
    filenames = {
        0.25: "input_reads_25_percent.fastq",
        0.5: "input_reads_50_percent.fastq",
        0.75: "input_reads_75_percent.fastq",
        0.9: "input_reads_90_percent.fastq"
    }

    # Generate synthetic FASTQ files if they don't exist
    for rate in duplication_rates:
        input_fastq = Path(filenames[rate])
        if not input_fastq.exists():
            print(f"Generating synthetic FASTQ file: {input_fastq} with duplication rate {rate*100}%")
            generate_synthetic_fastq(
                input_fastq, num_reads=100000,
                read_length=100, duplication_rate=rate
            )
        else:
            print(f"FASTQ file already exists: {input_fastq}")

    # Run deduplication profiling for each synthetic FASTQ file
    bucket_size = 4
    num_buckets = 50000
    k = 21
    for rate in duplication_rates:
        input_fastq = Path(filenames[rate])
        output_fastq = Path(f"output_cuckoo_{int(rate*100)}_percent.fastq")
        print(f"\nRunning Cuckoo Filter on {input_fastq}...")
        compare_cuckoo(input_fastq, output_fastq, bucket_size, num_buckets, k)

if __name__ == "__main__":
    main()
