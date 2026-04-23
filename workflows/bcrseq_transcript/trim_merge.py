#!/usr/bin/env python3
"""
Takes one forward and one reverse read file from Illumina platform paired-end
sequencing of B cell receptors and:
1) generates fastQC .html files to manually evaluate quality metrics of sequencing data
2) trims library prep adapter sequences & trims poor quality 3' ends of reads
3) merges the paired forward & reverse reads

The primary output file "basename.assembled.fastq" will be written to the same 
directory as the input R1 & R2 files, and other outputs & a log file will be 
saved to a new 'trim_merge' (quality control) subdirectory.

Use `python trim_merge.py --help` for instructions
"""

import datetime
import argparse
import subprocess
import os

def parse_arguments():
    """ Parse command-line arguments """
    parser = argparse.ArgumentParser(
        description='Generates .html files for quality assessment, trims adapters and low quality 3 prime read ends, and merges paired forward & reverse reads.')
    parser.add_argument('R1', metavar='forward_reads', type=str, 
        help='R1 fastq.gz file')
    parser.add_argument('R2', metavar='reverse_reads', type=str, 
        help='R2 fastq.gz file')
    parser.add_argument('--notrim', action='store_true', 
                        help='If set, skips the adapter & quality score trimming step before merging paired reads.')
    return parser.parse_args()

def run_fastqc(files, output_dir, log_pathway):
    """ Runs FastQC on a list of files and outputs the reports to the specified directory """
    fastqc_path = '/stor/work/Georgiou/Sharing_Folder/FastQC/fastqc'
    command = f"{fastqc_path} -o {output_dir} " + " ".join(files) + f" 2>&1 | tee -a {log_pathway}"
    subprocess.run(command, shell=True, check=True)
    return command

def call_cutadapt(r1, r2, output_dir, log_pathway):
    """ Trims adapters and low-quality ends from reads using cutadapt """
    # These are the TruSeq standard adapter sequences
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
    return trimmed_r1, trimmed_r2, cmd

def call_pear(r1, r2, output_dir, log_pathway):
    """ Calls PEAR read merger to find forward-reverse read couples based on
    overlapping sequence ends """
    
    base = os.path.basename(r1).replace('_trimmed.fastq.gz', '').replace('_R1','')
    out = os.path.join(output_dir, base)
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

def main(args):
    """ Main entry point of the module """
    
    # Create log file with current date and time in its name
    now = datetime.datetime.now()
    input_dir = os.path.dirname(args.R1)
    out_dir = os.path.join(input_dir, "trim_merge")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    log_name = f'log_QC-trimming-merging_{now.strftime("%Y-%m-%d_%H-%M-%S")}.txt'
    log_path = os.path.join(out_dir, log_name)

    # Run FastQC on initial input files
    qc_args = run_fastqc([args.R1, args.R2], out_dir, log_path)
    
    if not args.notrim:
        # Trim adapters & poor quality read ends
        trimmed_r1, trimmed_r2, trim_args = call_cutadapt(args.R1, args.R2, out_dir, log_path)
    else:
        # Or use original files if user specifies to skip trimming
        # In this case, copy the original files to the qc directory to keep generated files together
        trimmed_r1 = os.path.join(out_dir, os.path.basename(args.R1))
        trimmed_r2 = os.path.join(out_dir, os.path.basename(args.R2))
        subprocess.run(f"cp {args.R1} {trimmed_r1}", shell=True, check=True)
        subprocess.run(f"cp {args.R2} {trimmed_r2}", shell=True, check=True)
        trim_args = "Trimming skipped due to --notrim flag"
    
    # Call "PEAR" Paired End reAd mergeR to consolidate forward R1 & reverse R2 reads
    merge_args = call_pear(trimmed_r1, trimmed_r2, out_dir, log_path)

    # Determine base name for the assembled file
    base = os.path.basename(trimmed_r1).replace('_trimmed.fastq.gz', '').replace('_R1','')
    # The assembled file is originally created in the 'trim_merge' directory
    old_assembled_file = os.path.join(out_dir, f"{base}.assembled.fastq")
    # Move the assembled file to the input directory
    assembled_file = os.path.join(input_dir, f"{base}.assembled.fastq")
    subprocess.run(f"mv {old_assembled_file} {assembled_file}", shell=True, check=True)
    
    # Run FastQC on the output assembled file
    run_fastqc([assembled_file], out_dir, log_path)
    
    # Append the tool-calling command and the time it took to log file
    elapsed_str = time_passed(now)
    with open(log_path, 'a') as log:
        log.write(f'\n\nQUALITY CHECK COMMAND:\n{qc_args}'
                  f'\n\nADAPTER & PHRED TRIMMING COMMAND:\n{trim_args}'
                  f'\n\nMERGE COMMAND:\n{merge_args}'
                  f'\n\nTOTAL SCRIPT RUN TIME\n(HR:MIN:SEC):\n{elapsed_str}')

if __name__ == "__main__":
    """ This is executed when run from the command line """
    args = parse_arguments()
    main(args)
