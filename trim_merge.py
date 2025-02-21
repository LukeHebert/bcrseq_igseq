#!/usr/bin/env python3
"""
Takes one forward and one reverse read file from Illumina platform paired-end
sequencing of B cell receptors and:
1) generates fastQC .html files to manually evaluate quality metrics of sequencing data
2) trims library prep adapter sequences & trims poor quality 3' ends of reads
3) merges the paired forward & reverse reads

A log file will be created in your current working directory.

Example use:
python trim_merge.py r1.fastq.gz r2.fastq.gz
"""

__author__ = "Luke S Hebert"
__version__ = "0.1.0"
__license__ = "MIT"

import datetime
import argparse
import subprocess
import os

def main(args):
    """ Main entry point of the module """
    
    # Create log file with current date and time in its name
    now = datetime.datetime.now()
    directory = os.path.dirname(args.R1)
    log_name = f'log_QC-trimming-merging_{now.strftime("%Y-%m-%d_%H-%M-%S")}.txt'
    log_path = os.path.join(directory, log_name)

    # Run FastQC on initial input files
    qc_args = run_fastqc([args.R1, args.R2], directory, log_path)
    
    if not args.notrim:
        # Trim adapters & poor quality read ends
        trimmed_r1, trimmed_r2, trim_args = call_cutadapt(args.R1, args.R2, directory, log_path)
    else:
        # Or use original files if user specifies to skip trimming
        trimmed_r1, trimmed_r2 = args.R1, args.R2
        trim_args = "Trimming skipped due to --notrim flag"
    
    # Call "PEAR" Paired End reAd mergeR to consolidate forward R1 & reverse R2 reads
    merge_args = call_pear(trimmed_r1, trimmed_r2, log_path)

    # Run FastQC on the output assembled file
    out_base = trimmed_r1.replace('_trimmed.fastq.gz', '').replace('_R1','')
    assembled_file = f"{out_base}.assembled.fastq"
    run_fastqc([assembled_file], directory, log_path)
    
    # Append the tool-calling command and the time it took to log file
    elapsed_str = time_passed(now)
    with open(log_path, 'a') as log:
        log.write(f'\n\nQUALITY CHECK COMMAND:\n{qc_args}'
                  f'\n\nADAPTER & PHRED TRIMMING COMMAND:\n{trim_args}'
                  f'\n\nMERGE COMMAND:\n{merge_args}'
                  f'\n\nTOTAL SCRIPT RUN TIME\n(HR:MIN:SEC):\n{elapsed_str}')

def run_fastqc(files, output_dir, log_pathway):
    """ Runs FastQC on a list of files and outputs the reports to the specified directory """
    fastqc_path = '/stor/work/Georgiou/Sharing_Folder/FastQC/fastqc'
    command = f"{fastqc_path} -o {output_dir} " + " ".join(files) + f" 2>&1 | tee -a {log_pathway}"
    subprocess.run(command, shell=True, check=True)
    return command

def call_cutadapt(r1, r2, output_dir, log_pathway):
    """ Trims adapters and low-quality ends from reads using cutadapt """
    #These are the TruSeq standard adapter sequences
    adapter1 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'
    adapter2 = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
    quality_cutoff = '20'  # Quality score cutoff for 3' ends
    trimmed_r1 = os.path.join(output_dir, os.path.basename(r1).replace('.fastq.gz', '_trimmed.fastq.gz'))
    trimmed_r2 = os.path.join(output_dir, os.path.basename(r2).replace('.fastq.gz', '_trimmed.fastq.gz'))
    cutadapt_path = '/usr/local/bin/cutadapt'
    # Note: -m 50 (keep minimum length cutoff) is needed to prevent downstream PEAR merging errors with sequence records of length 0
    # Note: the --pair-filter=both is needed to prevent uneven output file lengths which also causes a downstream PEAR error
    cmd = f"{cutadapt_path} -q {quality_cutoff} -a {adapter1} -A {adapter2} -o {trimmed_r1} -p {trimmed_r2} {r1} {r2} -m 50 -j 0 2>&1 | tee -a {log_pathway}"
    subprocess.run(cmd, shell=True, check=True)
    subprocess.run(cmd, shell=True, check=True)
    return trimmed_r1, trimmed_r2, cmd

def call_pear(r1, r2, log_pathway):
    """ Calls PEAR read merger to find forward-reverse read couples based on
    overlapping sequence ends """
    
    out = r1.replace('_trimmed.fastq.gz', '').replace('_R1','')
    pear = ('/stor/work/Georgiou/Sharing_Folder/PEAR_0.9.11/'
            'pear-0.9.11-linux-x86_64/bin/pear')
    command = (f"{pear} -f {r1} -r {r2} -o {out} -v 10 -m 700 -n 50 -u 1 -j 20 "
               f"2>&1 | tee -a {log_pathway}")
    subprocess.run(command, shell=True, check=True)
    return command

def time_passed(start_time):
    """ Makes a human-readable string of elapsed time from start_time to 
    when this function is called"""
    elapsed_time = datetime.datetime.now() - start_time
    elapsed_hr = int(elapsed_time.total_seconds() // 3600)
    elapsed_min = int((elapsed_time.total_seconds() % 3600) // 60)
    elapsed_sec = int(elapsed_time.total_seconds() % 60)
    return f"{elapsed_hr:02}:{elapsed_min:02}:{elapsed_sec:02}"

if __name__ == "__main__":
    """ This is executed when run from the command line """
    
    # Get the arguments from terminal
    parser = argparse.ArgumentParser(
        description='Generates .html files for quality assessment, trims adapters and low quality 3 prime read ends, and merges paired forward & reverse reads.')
    parser.add_argument('R1', metavar='forward_reads', type=str, 
        help='R1 fastq.gz file')
    parser.add_argument('R2', metavar='reverse_reads', type=str, 
        help='R2 fastq.gz file')
    parser.add_argument('--notrim', action='store_true', 
                        help='If set, skips the adapter & quality score trimming step before merging paired reads.')
    args = parser.parse_args()
    
    # Pass the arguments to the main pairing function
    main(args)
