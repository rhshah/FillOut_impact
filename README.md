# FillOut_dmp

VCF FILES DESCRIPTION:
input_vcf_full_adepth_merged.vcf is an example ouput of the script when a full vcf file (vcf v4.2) is given.
input_vcf_8col_adepth_merged.vcf is an example ouput of the script when a 8 column vcf file (vcf v4.2) is given.


NAME

allele_depth_to_vcf - Reformat and merge variants from multiple alleledepth files into a single vcf.

SYNOPSIS

 perl allele_depth_to_vcf.pl --help
 perl allele_depth_to_vcf.pl --i_bam list_of_bams_files.txt --i_vcf file.vcf --qsub /common/sge/bin/lx-amd64/qsub --queue test.q

 Sample call for SGE:
 qsub -q test.q -V -N jobname -wd working_dir_path -e jobname.stderr -o jobname.stdout -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y perl allele_depth_to_vcf.pl --i_bam my_test_bams.txt --i_vcf input.vcf --o_dir output_dir_path --qsub qsub_executable_path --queue example.q --o_vcf name_for_output_vcf
 
 Sample call for LSF:
 bsub -q examplq_queue -J jobname -cwd working_dir_path -e jobname.stderr -o jobname.stdout -We 02:00 -R "rusage[mem=2]" -M 4 -n 1 "perl allele_depth_to_vcf.pl --i_vcf input.vcf --i_bam files_of_bam_paths.txt --bsub bsub_executable_path --queue example_queue --o_dir output_dir_path --o_vcf name_for_output_vcf"

OPTIONS

 --i_vcf          Path to input vcf file in vcf4.2 format. vcf4.2 file with eight columns is also accepted. (required)
 --o_dir     	  Path to output directory for all the output files of this program (optional)
 --o_vcf          Path to output multi-sample VCF [<o_dir>/allele_depth_merged.vcf] (optional)
 --g_allele       Path to dmp_genotype_allele.pl script DMP IMPACT pipeline [/dmp/resources/prod/software/dmp-impact-res/production/bin/dmp_genotype_allele.pl] (optional)
 --samtools       Path to samtools [default] (optional)
 --bedtools       Path to bedtools [default] (optional)
 --perl           Path to the perl binary that you want to use to run this program [/usr/bin/env perl] (optional)
 --qsub           Path to qsub executable for SGE (required if --bsub is not specified)
 --bsub           Path to qsub executable for LSF (required if --qsub is not specified)
 --queue          Name of the SGE / LSF queue where this software should run (required)
 --i_bam          Comma-separated absolute paths to bam files or a text file, with a .txt suffix, of all the absolute paths to bam files with line breaks (required)
 --RefFile        Path to reference genome file [default] (optional)
 --remove         Remove intermediate mpileup.alledepth and other files (optional)
 --help           Print help message and quit (optional)
 --man            Print detailed message and usage instructions and quit (optional)

DESCRIPTION

 allele_depth_to_vcf v1.0. This script genotypes the bam files provided using the MSK-IMPACT dmp_genotype_allele.pl script and merges all the resulting mpileup.alleledepth into a single VCF, in which each line contains a unique mutation variant with allele depth and variant frequency information for all the bam files.
 
AUTHORS

 Gowtham Jayakumaran ( jayakumg <at> mskcc <dot> org )

=cut
