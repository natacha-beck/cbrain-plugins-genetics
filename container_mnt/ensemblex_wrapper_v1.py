#!/usr/bin/python3

# This script is a wrapper for the ensemblex command line interface.
# It is used to run run the ensemblex pipeline from the command line.

import argparse
import subprocess
import sys
import os
import glob
import time
import types
import re
import shutil

# Add current directory to sys.path to ensure we can import wrapper_utils
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from wrapper_utils import print_header, print_log, is_dir, is_file, run_command

def handle_boolean(value, style=0):
    '''Converts a boolean value to the appropriate string representation for the config file.
       style=0: 'yes'/'no'
       style=1: 'TRUE'/'FALSE'
       style=2: 'T'/'F'
    '''
    if style == 0:
        return "yes"  if value else "no"
    elif style == 1:
        return "TRUE" if value else "FALSE"
    elif style == 2:
        return "T"    if value else "F"
    else:
        raise ValueError("Invalid style for boolean conversion.")

def handle_boolean_mixiture(args):
    args.PAR_vireo_forcelearnGT                                 = handle_boolean(args.PAR_vireo_forcelearnGT, 2)
    args.PAR_vartrix_umi                                        = handle_boolean(args.PAR_vartrix_umi, 1)
    args.PAR_ensemblex_merge_constituents                       = handle_boolean(args.PAR_ensemblex_merge_constituents, 0)
    args.PAR_ensemblex_probabilistic_weighted_ensemble          = handle_boolean(args.PAR_ensemblex_probabilistic_weighted_ensemble, 0)
    args.PAR_ensemblex_preliminary_parameter_sweep              = handle_boolean(args.PAR_ensemblex_preliminary_parameter_sweep, 0)
    args.PAR_ensemblex_graph_based_doublet_detection            = handle_boolean(args.PAR_ensemblex_graph_based_doublet_detection, 0)
    args.PAR_ensemblex_preliminary_ensemble_independent_doublet = handle_boolean(args.PAR_ensemblex_preliminary_ensemble_independent_doublet, 0)
    args.PAR_ensemblex_ensemble_independent_doublet                      = handle_boolean(args.PAR_ensemblex_ensemble_independent_doublet, 0)
    args.PAR_ensemblex_doublet_Demuxalot_threshold              = handle_boolean(args.PAR_ensemblex_doublet_Demuxalot_threshold, 0)
    args.PAR_ensemblex_doublet_Demuxalot_no_threshold           = handle_boolean(args.PAR_ensemblex_doublet_Demuxalot_no_threshold, 0)
    args.PAR_ensemblex_doublet_Demuxlet_threshold               = handle_boolean(args.PAR_ensemblex_doublet_Demuxlet_threshold, 0)
    args.PAR_ensemblex_doublet_Demuxlet_no_threshold            = handle_boolean(args.PAR_ensemblex_doublet_Demuxlet_no_threshold, 0)
    args.PAR_ensemblex_doublet_Souporcell_threshold             = handle_boolean(args.PAR_ensemblex_doublet_Souporcell_threshold, 0)
    args.PAR_ensemblex_doublet_Souporcell_no_threshold          = handle_boolean(args.PAR_ensemblex_doublet_Souporcell_no_threshold, 0)
    args.PAR_ensemblex_doublet_Vireo_threshold                  = handle_boolean(args.PAR_ensemblex_doublet_Vireo_threshold, 0)
    args.PAR_ensemblex_doublet_Vireo_no_threshold               = handle_boolean(args.PAR_ensemblex_doublet_Vireo_no_threshold, 0)
    args.PAR_ensemblex_compute_singlet_confidence               = handle_boolean(args.PAR_ensemblex_compute_singlet_confidence, 0)
    args.PAR_ensemblex_doublet_Freemuxlet_threshold             = handle_boolean(args.PAR_ensemblex_doublet_Freemuxlet_threshold, 0)
    args.PAR_ensemblex_doublet_Freemuxlet_no_threshold          = handle_boolean(args.PAR_ensemblex_doublet_Freemuxlet_no_threshold, 0)

    return args

def parse_args_list():
    parser = argparse.ArgumentParser(description="ensemblex v1 wrapper script")

    ##############
    # Setup args #
    ##############
    parser.add_argument("--method",      required=True, help="Pipeline method to execute")
    parser.add_argument("--working_dir", default="working_directory", required=False, type=str, help="Working directory for ensemblex (default: working_directory)")
    parser.add_argument("--output_dir",  default="ensemblex_output",  required=True,  type=str, help="Output directory for ensemblex results (default: ensemblex_output)")

    ##############
    # Files args #
    ##############
    parser.add_argument("--step",                 type=str,            help="Step of the ensemblex pipeline to execute (e.g., 'setup', 'demuxalot', 'freemuxlet', 'vireo', 'soupercell', 'ensemblexing').")
    parser.add_argument("--gene_expression",      type=is_file,        help="Gene expression bam file of the pooled samples (e.g., 10X Genomics possorted_genome_bam.bam)")
    parser.add_argument("--gene_expression_bai",  type=is_file,        help="Gene expression bam index file of the pooled samples (e.g., 10X Genomics possorted_genome_bam.bam.bai)")
    parser.add_argument("--barcodes",             type=is_file,        help="Barcodes tsv file of the pooled cells (e.g., 10X Genomics barcodes.tsv)")
    parser.add_argument("--pooled_samples",       type=is_file,        help="VCF file describing the genotypes of the pooled samples (optional, required for demultiplexing tools with prior genotype information)")
    parser.add_argument("--genome_reference",     type=str,            help="Genome reference fasta file (e.g., 10X Genomics)")
    parser.add_argument("--genome_reference_fai", type=str,            help="Genome reference fasta index file (e.g., 10X Genomics)")
    parser.add_argument("--genotype_reference",   type=is_file,        help="Population reference vcf file (e.g., 1000 Genomes Project)")

    ##################
    # Demuxalot args #
    ##################
    parser.add_argument("--PAR_demuxalot_genotype_names", default=None, type=str, help="List of Sample ID's in the sample VCF file (e.g., 'Sample_1,Sample_2,Sample_3').")
    parser.add_argument("--PAR_demuxalot_prior_strength", default=100, type=int, help="Read prior strength (default: 100).")
    parser.add_argument("--PAR_demuxalot_minimum_coverage", default=200, type=int, help="Minimum read coverage (default: 200).")
    parser.add_argument("--PAR_demuxalot_minimum_alternative_coverage", default=10, type=int, help="Minimum alternative read coverage (default: 10).")
    parser.add_argument("--PAR_demuxalot_n_best_snps_per_donor", default=100, type=int, help="Number of best SNPs for each donor to use for demultiplexing (default: 100).")
    parser.add_argument("--PAR_demuxalot_genotypes_prior_strength", default=1, type=int, help="Genotype prior strength (default: 1.0).")
    parser.add_argument("--PAR_demuxalot_doublet_prior", default=0.25, type=float, help="Doublet prior strength (default: 0.25).")

    #############################################################
    #                    Demuxlet args                          #
    # (Note: only for demuxlet with prior genotype information) #
    #############################################################
    parser.add_argument("--PAR_demuxlet_field", default="GT", type=str, help="Field to extract the genotypes (GT), genotype likelihood (PL), or posterior probability (GP) from the sample .vcf file (default: GT).")

    ##################################################################
    # Freemuxlet args                                                #
    # (Note: only for freemuxlet without prior genotype information) #
    ##################################################################
    parser.add_argument("--PAR_freemuxlet_nsample", type=int, help="Number of pooled samples.")

    ##############
    # Vireo args #
    ##############
    parser.add_argument("--PAR_vireo_N", type=int, default=None, help="Number of pooled samples.")
    parser.add_argument("--PAR_vireo_minMAF", type=float, default=0.1, help="Minimum minor allele frequency (default: 0.1).")
    parser.add_argument("--PAR_vireo_minCOUNT", type=int, default=20, help="Minimum aggregated count (default: 20).")
    parser.add_argument("--PAR_vireo_forcelearnGT", action='store_true', help="Whether or not to treat donor GT as prior only (default: False).")
    # Vireo args for prior genotype information
    parser.add_argument("--PAR_vireo_type", default="GT", type=str, help="Field to extract the genotypes (GT), genotype likelihood (PL), or posterior probability (GP) from the sample .vcf file (default: GT).")

    ###################
    # Souporcell args #
    ###################
    parser.add_argument("--PAR_minimap2", type=str, default=None, help="For information regarding the minimap2 parameters.")
    parser.add_argument("--PAR_freebayes", type=str, default=None, help="For information regarding the freebayes parameters.")
    parser.add_argument("--PAR_vartrix_umi", action='store_true', help="Whether or no to consider UMI information when populating coverage matrices.")
    parser.add_argument("--PAR_vartrix_mapq", type=int, default=30, help="Minimum read mapping quality.")
    parser.add_argument("--PAR_vartrix_threads", type=int, default=8, help="Number of threads for computing.")
    parser.add_argument("--PAR_souporcell_k", type=int, default=None, help="Number of pooled samples.")
    parser.add_argument("--PAR_souporcell_t", type=int, default=8, help="Number of threads for computing.")

    ##################
    # Ensemblex args #
    ##################
    # Pool parameters
    parser.add_argument("--PAR_ensemblex_sample_size", type=int, default=None, help="Number of samples multiplexed in the pool.")
    parser.add_argument("--PAR_ensemblex_expected_doublet_rate", type=float, default=None, help="Expected doublet rate for the pool, If using 10X Genomics, the expected doublet rate can be estimated based on the number of recovered cells.")

    # Set up parameters
    parser.add_argument("--PAR_ensemblex_merge_constituents", action='store_true', help="Whether or not to merge the output files of the constituent demultiplexing tools.")

    # Step 1 parameters: Probabilistic-weighted ensemble
    # PAR_ensemblex_probabilistic_weighted_ensemble	Yes	Whether or not to perform Step 1: Probabilistic-weighted ensemble. If running Ensemblex on a pool for the first time, this parameter should be set to "Yes". Subsequent runs of Ensemblex (e.g., parameter optimization) can have this parameter set to "No" as the pipeline will automatically detect the previously generated Step 1 output file.
    parser.add_argument("--PAR_ensemblex_probabilistic_weighted_ensemble", action='store_true', help="Whether or not to perform Step 1: Probabilistic-weighted ensemble.")

    # Step 2 parameters: Graph-based doublet detection
    parser.add_argument("--PAR_ensemblex_preliminary_parameter_sweep", action='store_true', help="Whether or not to perform a preliminary parameter sweep for Step 2: Graph-based doublet detection.")
    parser.add_argument("--PAR_ensemblex_nCD", type=int, default=None, help="Manually defined number of confident doublets in the pool (nCD).")
    parser.add_argument("--PAR_ensemblex_pT", type=float, default=None, help="Manually defined percentile threshold of the nearest neighour frequency (pT).")
    parser.add_argument("--PAR_ensemblex_graph_based_doublet_detection", action='store_true', help="Whether or not to perform Step 2: Graph-based doublet detection. If PAR_ensemblex_nCD and PAR_ensemblex_pT are not defined by the user (NULL), ensemblex will automatically determine the optimal parameter values using an unsupervised parameter sweep.\nIf PAR_ensemblex_nCD and PAR_ensemblex_pT are defined by the user, graph-based doublet detection will be performed with the user-defined values.")

    # Step 3 parameters: Ensemble-independent doublet detection
    parser.add_argument("--PAR_ensemblex_preliminary_ensemble_independent_doublet", action='store_true', help="Whether or not to perform a preliminary parameter sweep for Step 3: Ensemble-independent doublet detection.")
    parser.add_argument("--PAR_ensemblex_ensemble_independent_doublet", action='store_true', help="Whether or not to perform Step 3: Ensemble-independent doublet detection.")
    parser.add_argument("--PAR_ensemblex_doublet_Demuxalot_threshold", action='store_true', help="Whether or not to label doublets identified by Demuxalot as doublets.\nOnly doublets with assignment probabilities exceeding Demuxalot's recommended probability threshold will be labeled as doublets by ensemblex.")
    parser.add_argument("--PAR_ensemblex_doublet_Demuxalot_no_threshold", action='store_true', help="Whether or not to label doublets identified by Demuxalot as doublets, regardless of the corresponding assignment probability.")
    parser.add_argument("--PAR_ensemblex_doublet_Souporcell_threshold", action='store_true', help="Whether or not to label doublets identified by Souporcell as doublets.\nOnly doublets with assignment probabilities exceeding Souporcell's recommended probability threshold will be labeled as ensemblex.")
    parser.add_argument("--PAR_ensemblex_doublet_Souporcell_no_threshold", action='store_true', help="Whether or not to label doublets identified by Souporcell as doublets, regardless of the corresponding assignment probability.")
    parser.add_argument("--PAR_ensemblex_doublet_Vireo_threshold", action='store_true', help="Whether or not to label doublets identified by Vireo as doublets.\nOnly doublets with assignment probabilities exceeding Vireo's recommended probability threshold will be labeled as doublets by ensemblex.")
    parser.add_argument("--PAR_ensemblex_doublet_Vireo_no_threshold", action='store_true', help="Whether or not to label doublets identified by Vireo as doublets, regardless of the corresponding assignment probability.")

    # Step 3 specific parameters for demuxlet
    parser.add_argument("--PAR_ensemblex_doublet_Demuxlet_threshold", action='store_true', help="Whether or not to label doublets identified by Demuxlet as doublets.\nOnly doublets with assignment probabilities exceeding Demuxlet's recommended probability threshold will be labeled as doublets by ensemblex.")
    parser.add_argument("--PAR_ensemblex_doublet_Demuxlet_no_threshold", action='store_true', help="Whether or not to label doublets identified by Demuxlet as doublets, regardless of the corresponding assignment probability.")

    # Step 3 specific parameters for freemuxlet
    parser.add_argument("--PAR_ensemblex_doublet_Freemuxlet_threshold", action='store_true', help="Whether or not to label doublets identified by Freemuxlet as doublets.\nOnly doublets with assignment probabilities exceeding Freemuxlet's recommended probability threshold will be labeled as doublets by ensemblex.")
    parser.add_argument("--PAR_ensemblex_doublet_Freemuxlet_no_threshold", action='store_true', help="Whether or not to label doublets identified by Freemuxlet as doublets, regardless of the corresponding assignment probability.")


    # Confidence score parameters
    parser.add_argument("--PAR_ensemblex_compute_singlet_confidence", action='store_true', help="Whether or not to compute Ensemblex's singlet confidence score.\nThis will define low confidence assignments which should be removed from downstream analyses.")

    return handle_boolean_mixiture(parser.parse_args())

def build_file_io_ctx(args):
    file_io_ctx                            = types.SimpleNamespace()
    base_workdir                           = os.path.join(os.getcwd(), args.working_dir)
    file_io_ctx.working_dir                = base_workdir
    base_workdir_inputfiles                = os.path.join(base_workdir, "input_files")
    file_io_ctx.config_ensemblex_ini       = os.path.join(base_workdir, "job_info", "configs", "ensemblex_config.ini")
    file_io_ctx.pooled_bam_input           = os.path.join(base_workdir_inputfiles, "pooled_bam.bam")
    file_io_ctx.genotype_reference_input   = os.path.join(base_workdir_inputfiles, "reference.vcf")
    file_io_ctx.pooled_samples_input       = os.path.join(base_workdir_inputfiles, "pooled_samples.vcf")
    file_io_ctx.gene_expression_bai_input  = os.path.join(base_workdir_inputfiles, "pooled_bam.bam.bai")
    file_io_ctx.barcodes_input             = os.path.join(base_workdir_inputfiles, "pooled_barcodes.tsv")
    file_io_ctx.genome_reference_input     = os.path.join(base_workdir_inputfiles, "reference.fa")
    file_io_ctx.genome_reference_fai_input = os.path.join(base_workdir_inputfiles, "reference.fa.fai")

    return file_io_ctx

def create_option_dict(args, config_file_path):
    # read the config file and create a dictionary of the options
    method = args.method
    option_dict = {}
    option_dict[method] = {}
    with open(config_file_path, "r") as f:
        for line in f:
            if re.match(r"^##\s+.*\s+##", line):
                current_section = re.findall(r"^##\s+(.*)\s+##", line)[0].split()[0].lower()
                option_dict[method][current_section] = {}
            if re.match(r"^([A-Za-z0-9_]+)=(.*)", line):
                key, value = re.findall(r"^([A-Za-z0-9_]+)=(.*)", line)[0]
                if key == "METHOD":
                    value = value.strip().lower()
                    if value != args.method.lower():
                        print_log(f"Warning: METHOD in config file ({value}) does not match the method specified in the command line arguments ({args.method}).")
                        exit(1)
                else:
                    option_dict[method][current_section][key] = value
            if re.match(r"^#(PAR_[A-Za-z0-9_]+)=(.*)", line):
                key, value = re.findall(r"^#(PAR_[A-Za-z0-9_]+)=(.*)", line)[0]
                option_dict[method][current_section][key] = value

    return option_dict

def overwrite_ensemblex_config(args, file_io_ctx, option_dict):
    new_config_lines = []
    method           = args.method
    with open(file_io_ctx.config_ensemblex_ini, "r") as f:
        current_section = None
        for line in f:
            if re.match(r"^##\s+.*\s+##", line):
                current_section = re.findall(r"^##\s+(.*)\s+##", line)[0].split()[0].lower()
                new_config_lines.append(line)
            elif re.match(r"^([A-Za-z0-9_]+)=(.*)", line) or re.match(r"^#(PAR_[A-Za-z0-9_]+)=(.*)", line):
                if re.match(r"^([A-Za-z0-9_]+)=(.*)", line):
                    key, value = re.findall(r"^([A-Za-z0-9_]+)=(.*)", line)[0]
                else:
                    key, value = re.findall(r"^#(PAR_[A-Za-z0-9_]+)=(.*)", line)[0]

                if key == "METHOD":
                    value = value.strip().lower()
                    if value != args.method.lower():
                        print_log(f"Warning: METHOD in config file ({value}) does not match the method specified in the command line arguments ({args.method}).")
                        exit(1)
                    else:
                        new_config_lines.append(line)
                # Special case to handle CONTAINER_CMD
                elif key == "CONTAINER_CMD":
                    line = f"CONTAINER_CMD='{os.environ.get('WRAPPER_CONTAINER_CMD', 'singularity')}'\n"
                    new_config_lines.append(line)
                # Overwrite the parameters in the config file with the values from the command line arguments
                else:
                    if vars(args).get(key) is not None:
                        value = vars(args).get(key)
                        if isinstance(value, str) and re.search(r"\s", value):
                            value = f"'{value}'"
                        line = f"{key}={value}\n"
                        new_config_lines.append(line)
                    else:
                        new_config_lines.append(line)
            else:
                new_config_lines.append(line)

    # Move the old config file to a backup file
    backup_config_path = file_io_ctx.config_ensemblex_ini + ".bk." + time.strftime("%Y%m%d-%H%M%S")
    os.rename(file_io_ctx.config_ensemblex_ini, backup_config_path)
    # Write the new config file
    with open(file_io_ctx.config_ensemblex_ini, "w") as f:
        f.writelines(new_config_lines)

def files_preparation(args, file_io_ctx):
    '''Prepares input files for ensemblex.'''
    workdir_input_files_dir = os.path.join(os.getcwd(), file_io_ctx.working_dir, "input_files")

    # --gene_expression
    gene_expression_file      = args.gene_expression
    is_file(gene_expression_file)
    if os.path.exists(file_io_ctx.pooled_bam_input):
        print_log(f"file exists: {os.path.exists(file_io_ctx.pooled_bam_input)}")
    else:
        run_command(["cp", gene_expression_file, file_io_ctx.pooled_bam_input], shell=False)

    # --genotype_reference
    if args.genotype_reference:
        genotype_reference_file   = args.genotype_reference
        is_file(genotype_reference_file)
        if not os.path.exists(file_io_ctx.genotype_reference_input):
            run_command(["cp", genotype_reference_file, file_io_ctx.genotype_reference_input], shell=False)

    # --pooled_samples (optional)
    if args.pooled_samples:
        pooled_samples_file       = args.pooled_samples
        is_file(pooled_samples_file)
        if os.path.exists(file_io_ctx.pooled_samples_input):
            print_log(f"file exists: {os.path.exists(file_io_ctx.pooled_samples_input)}")
        else:
            run_command(["cp", pooled_samples_file, file_io_ctx.pooled_samples_input], shell=False)

    # --gene_expression_bai
    gene_expression_bai_file  = args.gene_expression_bai
    is_file(gene_expression_bai_file)
    if not os.path.exists(file_io_ctx.gene_expression_bai_input):
        run_command(["cp", gene_expression_bai_file, file_io_ctx.gene_expression_bai_input], shell=False)

    # --barcodes
    barcodes_file             = args.barcodes
    is_file(barcodes_file)
    if not os.path.exists(file_io_ctx.barcodes_input):
        run_command(["cp", barcodes_file, file_io_ctx.barcodes_input], shell=False)


    # --genome_reference
    ref_data_dir = os.environ.get("REF_DATA_DIR", "")
    genome_reference_file     = ref_data_dir + "/" + args.genome_reference
    is_file(genome_reference_file)
    if not os.path.exists(file_io_ctx.genome_reference_input):
        run_command(["cp", genome_reference_file, file_io_ctx.genome_reference_input], shell=False)

    # --genome_reference_fai
    genome_reference_fai_file = ref_data_dir + "/" + args.genome_reference_fai
    is_file(genome_reference_fai_file)
    if not os.path.exists(file_io_ctx.genome_reference_fai_input):
        run_command(["cp", genome_reference_fai_file, file_io_ctx.genome_reference_fai_input], shell=False)

def sort_vcf_files(args, file_io_ctx):
    '''Prepares input files for the demultiplexing step of ensemblex.
        - Sort pooled samples .vcf file
        - Sort reference .vcf file
    '''

    if args.pooled_samples:
        print_log("Performing files preparation for demultiplexing (sorting vcf files)...")
        launch_sort_pooled_cmd    = ["bash", os.path.expandvars("$ensemblex_HOME") + "/launch_ensemblex.sh", "-d", f"{file_io_ctx.working_dir}", "--step", "sort", "--vcf", file_io_ctx.pooled_samples_input, "--bam", file_io_ctx.pooled_bam_input, "--sortout", "pooled_samples_sorted.vcf"]
        run_command(launch_sort_pooled_cmd, False)
        print_log(f"Moving sorted pooled samples vcf to {file_io_ctx.pooled_samples_input}...")
        os.rename("pooled_samples_sorted.vcf", file_io_ctx.pooled_samples_input)

    launch_sort_reference_cmd = ["bash", os.path.expandvars("$ensemblex_HOME") + "/launch_ensemblex.sh", "-d", f"{file_io_ctx.working_dir}", "--step", "sort", "--vcf", file_io_ctx.genotype_reference_input, "--bam", file_io_ctx.pooled_bam_input, "--sortout", "reference_sorted.vcf"]
    run_command(launch_sort_reference_cmd, False)
    os.rename("reference_sorted.vcf", file_io_ctx.genotype_reference_input)

def verify_step_completion(args, file_io_ctx):
    '''Verifies the completion of the specified ensemblex step by checking for the existence of the expected output file.'''
    step_output_files = {
        "GT": {
            "setup":      "input_files/reference.vcf",
            "demuxalot":  "demuxalot/Demuxalot_result.csv",
            "demuxlet":   "demuxlet/outs.best",
            "souporcell": "souporcell/clusters.tsv",
            "vireo-GT":   "vireo_gt/donor_ids.tsv",
            "ensemblexing": "ensemblex_gt/confidence/ensemblex_final_cell_assignment.csv"
        },
        "noGT": {
            "setup":      "input_files/reference.vcf",
            "demuxalot":  "demuxalot/Demuxalot_result.csv",
            "freemuxlet": "freemuxlet/outs.clust1.samples*",
            "vireo":      "vireo/donor_ids.tsv",
            "souporcell": "souporcell/clusters.tsv",
            "ensemblexing": "ensemblex/confidence/ensemblex_final_cell_assignment.csv"
        }
    }

    expected_relative = step_output_files.get(args.method, {}).get(args.step)
    if not expected_relative:
        print_log(f"Error: No expected output file defined for method '{args.method}', step '{args.step}'.")
        exit(1)

    matches_expected_output_file = [file for file in glob.glob(os.path.join(file_io_ctx.working_dir, expected_relative)) if os.path.getsize(file) > 0]
    if matches_expected_output_file:
        expected_output_file = matches_expected_output_file[0]
        print_log(f"Step '{args.step}' completed successfully. Output file: {expected_output_file}")
    else:
        print_log(f"Error: Expected output files for step '{args.step}' not found.")
        exit(1)


def main():
    tool_name       = "ensemblex"
    tool_version    = "v1"
    wrapper_version = "1.0.0"
    wrapper_author  = "Natacha Beck: nbeck@mcin.ca"
    print_header(tool_name, tool_version, wrapper_version, wrapper_author)

    args          = parse_args_list()
    file_io_ctx   = build_file_io_ctx(args)
    option_dict   = None

    if args.step == "setup":
        if args.method == "GT":
            launch_ensemblex_init_withGT_cmd = ["bash", os.path.expandvars("$ensemblex_HOME") + "/launch_ensemblex.sh", "-d", file_io_ctx.working_dir, "--step", "init-GT"]
            run_command(launch_ensemblex_init_withGT_cmd, False)
            option_dict = create_option_dict(args, file_io_ctx.config_ensemblex_ini)
        elif args.method == "noGT":
            launch_ensemblex_init_withoutGT_cmd = ["bash", os.path.expandvars("$ensemblex_HOME") + "/launch_ensemblex.sh", "-d", file_io_ctx.working_dir, "--step", "init-noGT"]
            run_command(launch_ensemblex_init_withoutGT_cmd, False)
            option_dict = create_option_dict(args, file_io_ctx.config_ensemblex_ini)

        overwrite_ensemblex_config(args, file_io_ctx, option_dict)

        files_preparation(args, file_io_ctx)
        sort_vcf_files(args, file_io_ctx)
        verify_step_completion(args, file_io_ctx)
    else:
        if file_io_ctx.working_dir is None:
            print_log("Error: --working_dir must be specified for steps other than 'setup'.")
            exit(1)
        else:
            option_dict = create_option_dict(args, file_io_ctx.config_ensemblex_ini)

            overwrite_ensemblex_config(args, file_io_ctx, option_dict)
            launch_ensemblex_cmd = ["bash", os.path.expandvars("$ensemblex_HOME") + "/launch_ensemblex.sh", "-d", file_io_ctx.working_dir, "--step", args.step]
            run_command(launch_ensemblex_cmd, False)

            verify_step_completion(args, file_io_ctx)



    # After running ensemblex, move the output files to the output directory specified by the user
    try:
        os.rename(file_io_ctx.working_dir, args.output_dir)
        print_log(f"Output directory {args.output_dir} created successfully.")
    except FileExistsError:
        print_log(f"Error: Output directory {args.output_dir} already exists.")

if __name__ == "__main__":
    main()
