# profile_ngsreads.py

import sys
import os
import time
import psutil
import subprocess
import random
import gc
import matplotlib.pyplot as plt
from Bio import SeqIO

# Optional: Use pathlib for better path handling (recommended)
from pathlib import Path

# Add the parent directory to sys.path (if necessary for imports)
current_dir = Path(__file__).resolve().parent
parent_dir = current_dir.parent
sys.path.insert(0, str(parent_dir))

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

def run_ngsreads_treatment(input_fastq):
    """Run NGSReadsTreatment and profile its performance using psutil."""
    # Correctly define the jar file path using a raw string or forward slashes
    jar_file = r"Profiling\NgsReadsTreatment_v1.3.jar"  # Option 1: Raw string
    # Alternatively, use forward slashes
    # jar_file = "Profiling/NgsReadsTreatment_v1.3.jar"

    # Verify that the JAR file exists
    if not os.path.exists(jar_file):
        raise FileNotFoundError(f"The JAR file was not found: {jar_file}")

    # Construct the command
    command = f"java -jar {jar_file} {input_fastq} 2"

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

def format_fastq(input_fastq, output_fastq):
    """
    Format the FASTQ file by removing lines before the first '@' and ensuring proper FASTQ format.

    Args:
        input_fastq (str): Path to the input FASTQ file.
        output_fastq (str): Path to the output formatted FASTQ file.
    """
    found_first_at = False
    with open(input_fastq, 'r') as infile, open(output_fastq, 'w') as outfile:
        for line in infile:
            if not found_first_at:
                if line.startswith('@'):
                    found_first_at = True
                    outfile.write(line)
            else:
                outfile.write(line)

    # Optional: Validate the FASTQ format
    with open(output_fastq, 'r') as fh:
        lines = fh.readlines()
        if len(lines) % 4 != 0:
            raise ValueError(f"Formatted FASTQ file {output_fastq} has incomplete records.")
        for i in range(0, len(lines), 4):
            if not lines[i].startswith('@'):
                raise ValueError(f"Record starting at line {i+1} in {output_fastq} does not start with '@'.")
            if not lines[i+2].startswith('+'):
                raise ValueError(f"Record starting at line {i+1} in {output_fastq} does not have '+' separator.")

    print(f"Formatted FASTQ file: {output_fastq}")

def count_reads(fastq_file, unique=False):
    """Count total or unique reads in a FASTQ file."""
    if unique:
        unique_sequences = set()
        for record in SeqIO.parse(fastq_file, "fastq"):
            unique_sequences.add(str(record.seq))
        return len(unique_sequences)
    return sum(1 for _ in SeqIO.parse(fastq_file, "fastq"))

def compare_ngsreads(input_fastq, output_fastq):
    """Run and profile NGSReadsTreatment deduplication."""
    total_reads = count_reads(input_fastq)

    time_taken, memory_used = run_ngsreads_treatment(input_fastq)


    treated_fastq = input_fastq.replace(".fastq", "_1_trated.fastq")

    if not os.path.exists(treated_fastq):
        print(f"Warning: Expected output file {treated_fastq} not found.")
        unique_reads_ngs = 0
        dedup_percentage_ngs = 0.0
    else:
        # Define a temporary formatted FASTQ file
        formatted_fastq = input_fastq.replace(".fastq", "_1_trated_formatted.fastq")
        
        # Format the treated FASTQ file
        try:
            format_fastq(treated_fastq, formatted_fastq)
        except ValueError as ve:
            print(f"Error formatting FASTQ file {treated_fastq}: {ve}")
            unique_reads_ngs = 0
            dedup_percentage_ngs = 0.0
        else:
            # Count unique reads from the formatted FASTQ file
            unique_reads_ngs = count_reads(formatted_fastq, unique=True)
            dedup_percentage_ngs = (1 - unique_reads_ngs / total_reads) * 100

    print("\nNGSReadsTreatment Profiling Results:")
    print(f"Time Taken: {time_taken:.4f} seconds")
    print(f"Peak Memory Used: {memory_used:.4f} MB")
    print(f"Deduplication Percentage: {dedup_percentage_ngs:.2f}%")

def main():
    """Main execution function for NGSReadsTreatment profiling."""
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
        input_fastq = filenames[rate]
        if not os.path.exists(input_fastq):
            print(f"Generating synthetic FASTQ file: {input_fastq} with duplication rate {rate*100}%")
            generate_synthetic_fastq(input_fastq, num_reads=100000, read_length=100, duplication_rate=rate)
        else:
            print(f"FASTQ file already exists: {input_fastq}")

    # Run deduplication profiling for each synthetic FASTQ file
    output_fastq = "output"
    print("\nRunning NGSReadsTreatment Deduplication Profiling")
    for rate in duplication_rates:
        input_fastq = filenames[rate]
        print(f"\nRunning NGSReadsTreatment on {input_fastq}")
        compare_ngsreads(input_fastq, output_fastq)

if __name__ == "__main__":
    main()
