#!/usr/bin/env python3 
"""
This script takes a ".assembled.fastq" file from (merged, if needed) BCRseq 
sequences as input and calls a tool of the user's choice to map the V, D, J, 
and C genes within reads to their germline sequences. The tool also finds 
frameworks and  complementarity determining regions, simulates translation, 
and calculates somatic hypermutation rate.

Use `identify_genes.py --help` for instructions
"""

import datetime
import argparse
import subprocess
import os
from Bio import SeqIO

def call_MiXCR(infile, species, left, right, nano, log_file_name, output_dir):
    """ Calls MiLab's MiXCR tool for mapping and annotating gene elements from
    BCR mRNA sequences."""
    
    # To parse the user's species argument
    spp = {'human':'hsa', 'hsa':'hsa', 'mouse':'mmu', 'mmu':'mmu', 'dog':'dog'}
    
    # Determine output file names based on input file type and save them to output_dir
    base = os.path.basename(infile)
    if infile.lower().endswith(('.fastq', '.fq')):
        alignment_filename = base.replace('assembled.fastq', f'{spp[species.lower()].upper()}alignment.vdjca')
        not_aligned_filename = base.replace('.fastq','_na.fastq')
    elif infile.lower().endswith(('.fasta', '.fa')):
        alignment_filename = base.replace('assembled.fasta', f'{spp[species.lower()].upper()}alignment.vdjca')
        not_aligned_filename = base.replace('.fasta','_na.fasta')
    else:
        raise ValueError("Input file must be in FASTA or FASTQ format.")
    
    alignment = os.path.join(output_dir, alignment_filename)
    not_aligned = os.path.join(output_dir, not_aligned_filename)

    # For the final TSV file, keep it in the same directory as the input file
    input_dir = os.path.dirname(os.path.abspath(infile))
    final_tsv = os.path.join(input_dir, os.path.basename(alignment).replace('alignment.vdjca','MiXCR.tsv'))
    
    # Call MiXCR to align
    mixcr = ('/stor/work/Georgiou/Sharing_Folder/MiXCR_4.3.2/mixcr-4.3.2/mixcr')
    # for atypical species, an IMGT-based reference library is needed
    if spp[species.lower()] not in ['hsa','mmu']:
        if nano:
            command1 = (f"{mixcr} -Xmx30g align "
                f"--preset ont-rna-seq-vdj-full-length "  # nanopore data preset
                f"--keep-non-CDR3-alignments "  # very important for amplicons that aren't meant to include all of the CDR3
                f"--rna -s {spp[species.lower()]} "
                f"--library imgt.202312-3.sv8.json "
                f"--not-aligned-R1 {not_aligned} "
                f"-f {infile} {alignment} "
                f"2>&1 | tee -a {log_file_name}")      
        else:
            command1 = (f"{mixcr} -Xmx30g align --preset generic-bcr-amplicon "
                f"--floating-left-alignment-boundary  {left} "
                f"--floating-right-alignment-boundary {right} "
                f"--rna -s {spp[species.lower()]} "
                f"--library imgt.202312-3.sv8.json "
                f"--not-aligned-R1 {not_aligned} "
                f"-f {infile} {alignment} "
                f"2>&1 | tee -a {log_file_name}")      
    # for human and mouse, MiXCR does not need an external library
    else:
        if nano:
            command1 = (f"{mixcr} -Xmx30g align "
                f"--preset ont-rna-seq-vdj-full-length "  # nanopore data preset
                f"--keep-non-CDR3-alignments "  # very important for amplicons that aren't meant to include all of the CDR3
                f"--rna -s {spp[species.lower()]} -f {infile} {alignment} "
                f"2>&1 | tee -a {log_file_name}")
        else:
            command1 = (f"{mixcr} -Xmx30g align --preset generic-bcr-amplicon "
                f"--floating-left-alignment-boundary  {left} "
                f"--floating-right-alignment-boundary {right} "
                f"--rna -s {spp[species.lower()]} -f {infile} {alignment} "
                f"2>&1 | tee -a {log_file_name}")
    print('\tCalling MiXCR for alignment...')
    subprocess.run(command1, shell=True, check=True)
    
    # Call MiXCR to export alignments as a .tsv in the input file's directory
    command2 = (f"/stor/work/Georgiou/Sharing_Folder/MiXCR_4.3.2/mixcr-4.3.2/"
        f"mixcr -Xmx30g exportAlignments "
        f"-f -targetSequences -readIds "
        f"-vGene -dGene -jGene -cGene "
        f"-nFeature FR1 -nFeature CDR1 -nFeature FR2 -nFeature CDR2 "
        f"-nFeature FR3 -nFeature CDR3 -nFeature FR4 "
        f"-aaFeature FR1 -aaFeature CDR1 -aaFeature FR2 -aaFeature CDR2 "
        f"-aaFeature FR3 -aaFeature CDR3 -aaFeature FR4 "
        f"-vBestIdentityPercent -dBestIdentityPercent -jBestIdentityPercent "
        f"-cBestIdentityPercent -positionOf CBegin -positionOf CExon1End -positionOf CEnd "
        f"{alignment} {final_tsv} 2>&1 | tee -a {log_file_name}")
    print('\tCalling MiXCR to export alignments...')
    subprocess.run(command2, shell=True, check=True)
    
    return '\n\n'.join([command1, command2])

def call_igblast(infile, species, log_file_name, output_dir):
    """ Calls BLAST-based NCBI tool IgBLAST for mapping and annotating gene 
    elements from BCR mRNA sequences."""
    
    # If input is FASTQ, convert to FASTA and save the converted file to output_dir.
    # If already FASTA, use as is.
    if infile.lower().endswith(('.fastq', '.fq')):
        print('\tConverting .fastq to .fasta...')
        base = os.path.basename(infile)
        infasta_filename = base.replace('.fastq','.fasta')
        infasta = os.path.join(output_dir, infasta_filename)
        sequences = SeqIO.parse(infile, "fastq")
        SeqIO.write(sequences, infasta, "fasta")
    elif infile.lower().endswith(('.fasta', '.fa')):
        print('\tInput file is in FASTA format, no conversion needed.')
        infasta = infile
    else:
        raise ValueError("Input file must be in FASTA or FASTQ format.")

    # Call IgBLAST & output the annotation .tsv in AIRR format ("19")
    basename = os.path.splitext(os.path.basename(infasta))[0]
    outnametsv = f"{basename}_IgBLAST.tsv"

    igblast_dir = ('/stor/work/Georgiou/Sharing_Folder/IgBLAST_1.21.0/'
                   'ncbi-igblast-1.21.0')
    igblastn = ('/stor/work/Georgiou/Sharing_Folder/IgBLAST_1.21.0/'
                'ncbi-igblast-1.21.0/bin/igblastn')
    # Use the directory of the input file for the final TSV output
    input_dir = os.path.dirname(os.path.abspath(infile))
    command = (f"cd {igblast_dir}; "
               f"{igblastn} "
               f"-germline_db_V internal_data/{species}/{species}_V "
               f"-germline_db_J internal_data/{species}/{species}_J "
               f"-germline_db_D internal_data/{species}/{species}_D "
               f"-c_region_db internal_data/{species}/{species}_C "
               f"-auxiliary_data optional_file/{species}_gl.aux "
               f"-organism {species} "
               f"-query {infasta} "
               f"-outfmt 19 "
               f"-num_threads 40 "
               f"-out {input_dir}/{outnametsv} "
               f"2>&1 | tee -a {log_file_name} "
               f"2>/dev/null")
    print('\tCalling IgBLAST...')
    subprocess.run(command, shell=True, check=True)
    return command

def time_passed(start_time):
    """ Makes a human-readable string of ellapsed time from start_time to 
    when this function is called"""
    elapsed_time = datetime.datetime.now() - start_time
    elapsed_hr = int(elapsed_time.total_seconds() // 3600)
    elapsed_min = int((elapsed_time.total_seconds() % 3600) // 60)
    elapsed_sec = int(elapsed_time.total_seconds() % 60)
    return f"{elapsed_hr:02}:{elapsed_min:02}:{elapsed_sec:02}"

def parse_arguments():
    """Handles argument parsing for the script."""
    parser = argparse.ArgumentParser(description=('Calls a mapping+alignment'
        '+annotation tool with preset parameters to annotate each read\'s '
        'genetic elements.'))
    parser.add_argument('reads', metavar='input_reads', type=str, 
        help='Input (assembled aka merged) FASTQ or FASTA file name.') 
    parser.add_argument('sp', metavar='species', type=str, 
        help=('Species of sample donor e.g. "human", "mouse", "ferret"'
            ' (if using -t igblast); don\'t include quotations.'))    
    parser.add_argument("-t", "--tool", default="igblast", 
        choices=["mixcr", "igblast"], help=("Software that performs the gene "
        "mapping and annotations."))
    
    # MiXCR-specific optional arguments
    parser.add_argument('--leftbound', type=str, default='FR1Begin', 
        help=('If using MiXCR: leftmost boundary that will determine alignment'
            ' area. Default is FR1Begin. See MiXCR anchorpoints webpage for more options.'))
    parser.add_argument('--rightbound', type=str, default='FR4End', 
        help=('If using MiXCR: rightmost boundary that will determine alignment'
            ' area. Default is FR4End. See MiXCR anchorpoints webpage for more options.'))
    parser.add_argument('--nanopore', action='store_true', 
        help=('If specifically calling MiXCR to annotate sequences from Oxford Nanopore Technology platform.'))
    
    return parser.parse_args()

def main(args):
    """ Main entry point of the app """
    # Determine the directory of the input file
    input_dir = os.path.dirname(os.path.abspath(args.reads))
    
    # Create a new output directory (inside the input file's directory)
    output_dir = os.path.join(input_dir, "identify_genes")
    os.makedirs(output_dir, exist_ok=True)
    
    # Create log file with current date and time in its name in the new directory
    now = datetime.datetime.now()
    log_name = f'log_mapping_{now.strftime("%Y-%m-%d_%H-%M-%S")}.txt'
    log_path = os.path.join(output_dir, log_name)

    # Call the specified (or default: IgBLAST) gene mapping & annotation tool
    if args.tool == 'mixcr':
        args_str = call_MiXCR(args.reads, args.sp, args.leftbound, 
            args.rightbound, args.nanopore, log_path, output_dir)
    elif args.tool == 'igblast':
        args_str = call_igblast(args.reads, args.sp, log_path, output_dir)
    else:
        print('Invalid mapping tool name was used as a -t/--tool argument.')
        exit()

    # Append the tool-calling command and the time it took to log file
    elapsed_str = time_passed(now)
    with open(log_path, 'a') as log:
        log.write(f'\n\nCOMMAND(s):\n{args_str}\n\n'
            f'MAPPING-ANNOTATING TIME\n(HR:MIN:SEC):\n{elapsed_str}')

if __name__ == "__main__":
    """ This is executed when run from the command line """
    args = parse_arguments()
    main(args)
