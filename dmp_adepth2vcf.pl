#!/usr/bin/env perl

# dmp_adepth2vcf.pl - Genotype one or more bams; reformat and merge allelic info of all variants from the the resulting multiple alleledepth files into a single vcf.
# version: v1.2
# Author: Gowtham Jayakumaran ( jayakumg@mskcc.org )
# Adapted from IMPACT pipeline.
# Date: 25/Apr/2016
# Last Modified: 05/May/2016
# Modification timeline ( for major changes )
# 28/Apr/2016 --- v1.1 - Added support for batch processing and the option to handle either qsub or bsub.
# 05/May/2016 --- v1.2 - Updated qsub/bsub wait option. Made compatible with the impact pipeline.

use strict;
use warnings;
use IO::File;
use Getopt::Long qw( GetOptions );
use Pod::Usage qw( pod2usage );
use Cwd qw( abs_path );
use MSKCC_DMP_Logger;

# Initiate run time
my $now = time;

# parse options and print error, if any, or print usage, if help or man requested.
my( $help, $man, $delete, $i_vcf, $g_allele, $samtools, $bedtools, $RefFile, $i_bam, $mmq, $mbq, $o_dir, $identifier, $o_vcf, $perl, $qsub, $bsub, $queue );
my $logger = MSKCC_DMP_Logger->get_logger('ADEPTH_TO_VCF');
$logger->start_local();

if( @ARGV < 1 or
	$ARGV[0] !~ m/^-/ or
	!GetOptions(
		'help!' => \$help,
		'man!' => \$man,
		'delete!' => \$delete,
		'i_vcf=s' => \$i_vcf,
		'g_allele=s' => \$g_allele,
		'samtools=s' => \$samtools,
		'bedtools=s' => \$bedtools,
		'RefFile=s' => \$RefFile,
		'i_bam=s' => \$i_bam,
		'mmq=s' => \$mmq,
		'mbq=s' => \$mbq,
		'o_dir=s' => \$o_dir,
		'identifier=s' => \$identifier,
		'o_vcf=s' => \$o_vcf,
		'perl=s' => \$perl,
		'qsub=s' => \$qsub,
		'bsub=s' => \$bsub,
		'queue=s' => \$queue
		) )	{
	pod2usage( -verbose => 1, -input => \*DATA, -exitval => 2 );
}
pod2usage( -verbose => 1, -input => \*DATA, -exitval => 0 ) if( $help );
pod2usage( -verbose => 2, -input => \*DATA, -exitval => 0 ) if( $man );

$logger->info( "Starting $0: Setting default paths and constants ..." );
# obtain absolute paths for the input and output files
for my $file ( $i_vcf, $g_allele, $samtools, $bedtools, $RefFile, $i_bam, $o_dir, $o_vcf, $perl, $qsub, $bsub )	{
	$file = abs_path( $file ) if defined $file;
}

# Setting house keeping parameters
if( ( !defined $delete ) or ( $delete != 1 ) or ( $delete != 2) ) 	{
	$delete = 2;
	$logger->warn( "No or invalid parameter given for the argument --delete. Default is set and house keeping will be performed if $0 completes execution without errors." );
	}else	{
		$logger->info( "House keeping will be performed if $0 completes without any errors." ) if( $delete == 2 );
		$logger->info( "House keeping will not be performed." ) if( $delete == 1 );
	}

# set default paths and constants
$g_allele = "/dmp/resources/prod/software/dmp-impact-res/VERSIONS/dmp-impact-res-v118/bin/dmp_genotype_allele.pl" unless defined $g_allele;
$samtools = "/dmp/resources/prod/tools/bio/samtools/production/samtools" unless defined $samtools;
$perl = "$^X" unless defined $perl and ( -e $perl ); # /usr/bin/env perl
$bedtools = "/home/shahr2/Software/BEDTools/current/bin/" unless defined $bedtools;
$RefFile = "/home/shahr2/Resources/IMPACT/References/Homo_sapiens_assembly19_formatted.fasta" unless defined $RefFile;

# default parameters and constants for genotype_allele script. Arguments for mpileUpOutFile and outFile (alleledepth) will be set to default by the genotype_allele script.
if( !defined $mbq )
{
    $logger->warn( "Minimum base quality (-mbq) is not provided. It will be set to 5 for genotyping all the bam files." );
    $mbq = 5;
}
if( !defined $mmq )
{
    $logger->warn( "Minimum mapping quality (-mmq) is not provided. It will be set to 5 for genotyping all the bam files." );
    $mmq = 5;
}

die $logger->fatal( "Please provide either --qsub or --bsub, and not both." ) if ( ( $qsub and $bsub ) or ( !$qsub and !$bsub ) );
die $logger->fatal( "Please provide the -queue on which you want this program to run." ) unless defined $queue;
die $logger->fatal( "Please provide the path for your genotype script." ) unless ( -e $g_allele );
die $logger->fatal( "Please provide a vcf file with the .vcf suffix." ) unless defined $i_vcf and ( -e $i_vcf ) and ( $i_vcf =~ m/(\.)?(vcf)?$/ );
die $logger->fatal( "Please provide either a name for the output vcf file or a unique identifier (eg. Patient MRN), and not both." ) if ( ( $o_vcf and $identifier ) or ( !$o_vcf and !$identifier ) );
die $logger->fatal( "Please provide a proper argument for --i_bam." ) unless defined $i_bam and ( -e $i_bam );
die $logger->fatal( "Please provide the path to samtools." ) unless defined $samtools and ( -e $samtools );
die $logger->fatal( "Please provide the path to bedtools." ) unless defined $bedtools and ( -e $bedtools );
die $logger->fatal( "Please provide the path to a reference genome hg19.fasta." ) unless defined $RefFile and ( -e $RefFile );
die $logger->fatal( "Please provide the path to perl." ) unless defined $perl and ( -e $perl );

# set output dir
unless( defined $o_dir and ( -e $o_dir ) )	{
	$o_dir = `pwd`;
	chomp $o_dir;
	$logger->warn( "Path to output directory not provided. Current working directory will be used." );
}

# set output .vcf file. To be standardized later
if( defined $o_vcf )	{
	map{ s/(.+\/)?/$o_dir\//; s/(\.)?(maf|tsv|csv|txt|vcf)?$/.vcf/ } $o_vcf;
	$o_vcf .= ".vcf" unless( $o_vcf =~ m/(\.)?(vcf)?$/ );
}else	{
	$identifier =~ s/^\s+|\s+$|\r|\n//g;
	$o_vcf = $o_dir."/".$identifier."_allele_depth_merged.vcf";
}
$logger->info( "Output vcf will be printed to $o_vcf." );

# parse the arguments to --i_bam depending on whether its a comma-separated string of bam file paths or a text file with bam file paths
my @bam = ();
if ($i_bam =~ /.txt$/ )	{
	my $bam_file = IO::File->new( $i_bam ) or die $logger->fatal( "Cannot open file $i_bam." );
	$logger->info( "Bam files are provided as a file of files, will one file path per line." );
	while ( <$bam_file> )	{
		next if ( $_ =~ /^#/ );
		$_ =~ s/^\s+|\s+$|\r|\n//g;
		push @bam, $_;
	}
	$bam_file->close;
}else {
	@bam = map{ s/^\s+|\s+$|\r|\n//g; $_ } split( /,/, $i_bam );
	$logger->warn( "Bam files are provided as comma-separated string of bam file paths." );
}

# check whether the vcf file provided is in the proper format. A more stringent check can be implemented, if necessary.
my $vcf_header = `grep '^#' $i_vcf`;
$vcf_header =~ s/\n//g;
map{ ( m/^\#\#fileformat=VCFv4.2.+\#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO/ ) or 
	die $logger->fatal( "ERROR: $i_vcf is not in the vcf4.2 format. Please check https://samtools.github.io/hts-specs/VCFv4.2.pdf for more info." ) } $vcf_header;

# Extracting sample names from bam file, to be included in the VCF
my @bam_sample_names = ();
foreach my $bam_file ( @bam )	{
	die $logger->fatal( "ERROR: $bam_file does not exist. ERROR: $!" ) unless ( -e $bam_file );
	my $cmd = "$samtools view -H $bam_file";
	$logger->debug( "COMMAND: $cmd" );
	my $header = `$samtools view -H $bam_file`;
	map{ s/\n//g; s/.+\@RG\sID:(\S+)\s.+/$1/; } $header;
	die $logger->fatal( "Could not get sample name by parsing the header of $bam_file using $samtools view -H." ) unless $header;
	push @bam_sample_names, $header;
}

# set the absolute path and bam file names of all mpileup.alleledepth output files in $o_dir
my @alleledepth_files = @bam;
my @alleledepth_file_names = @bam;
map{ s/(.+\/)?/$o_dir\//; s/(\.bam)?$/_mpileup.alleledepth/; } @alleledepth_files;
map{ s/(.+\/)?//; s/(\.bam)?$/_mpileup.alleledepth/; } @alleledepth_file_names;

# Make a dummy wait file to run after process
MakeCSH( $o_dir );

# run dmp_genotype_allele.pl on each of the input bam files and direct the ouput mpileup.alleledepth files to $o_dir
 my @notifyNames = ();
 my @deletefiles = ();

my $bam_file_count = 0;
while( my $bam_file = shift @bam )	{
		die $logger->fatal( "$bam_file does not exist. ERROR: $!" ) unless ( -e $bam_file );
		my $job_name = ++$bam_file_count;
		$logger->info( "Genotyping $bam_file using variant information from $i_vcf." );
		eval {
		if( $qsub )	{
			my $cmd = "$qsub -q $queue -V -N genotype.$job_name.$$ -wd $o_dir -e genotype.$job_name.$$.stderr -o genotype.$job_name.$$.stdout -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y $perl $g_allele -o $o_dir -fmv $i_vcf -bam $bam_file -rf $RefFile -s $samtools -o $o_dir -b $bedtools -mmq $mmq -mbq $mbq -qsub $qsub -q $queue"; # For arguments --outFile|-of, --mpileUpOutFile|-mof, and --bamId|-bi in $g_allele script, no parameters is supplied here and the default will be set.
			$logger->debug( "COMMAND: $cmd" );
			`$cmd`;
			my $cmd2 = "$qsub -q $queue -V -wd $o_dir -hold_jid genotype.$job_name.$$ -N Notify_genotype.$job_name.$$ -e Notify_genotype.$job_name.$$.stderr -o Notify_genotype.$job_name.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y \"$o_dir/Notify.csh\"";
			$logger->debug( "COMMAND: $cmd2" );
			`$cmd2`;
		}else	{
			my $cmd = "$bsub -q $queue -J genotype.$job_name.$$ -cwd $o_dir -e genotype.$job_name.$$.stderr -o genotype.$job_name.$$.stdout -We 24:00 -R \"rusage[mem=2]\" -M 4 -n 1 \"$perl $g_allele -fmv $i_vcf -bam $bam_file -rf $RefFile -s $samtools -o $o_dir -q $queue -b $bedtools -mmq $mmq -mbq $mbq -bsub $bsub\"";# For arguments --outFile|-of, --mpileUpOutFile|-mof, and --bamId|-bi in $g_allele script, no parameters is supplied here and the default will be set.
			$logger->debug( "COMMAND: $cmd" );
			`$cmd`;
			my $cmd2 = "$bsub -q $queue -cwd $o_dir -w \"post_done(genotype.$job_name.$$)\" -J Notify_genotype.$job_name.$$ -e Notify_genotype.$job_name.$$.stderr -o Notify_genotype.$job_name.$$.stat -We 24:00 -R \"rusage[mem=2]\" -M 4 -n 1 \"$o_dir/Notify.csh\"";
			$logger->debug( "COMMAND: $cmd2" );
			`$cmd2`;
		}
	};
die $logger->fatal( "Genotype bam files: Job submission failed. ERROR: $@" ) if ($@);
push( @notifyNames, "Notify_genotype.".$job_name.".".$$.".stat" );
push( @deletefiles, "genotype.".$job_name.".".$$.".stderr" );
push( @deletefiles, "genotype.".$job_name.".".$$.".stdout" );
push( @deletefiles, "Notify_genotype.".$job_name.".".$$.".stderr" );
push( @deletefiles, "Notify_genotype.".$job_name.".".$$.".stat" );
}
push( @deletefiles, "Notify.csh" );

# Waiting for all the bam files to be genotyped and the resulting mpileup.alleledepth files generated.
WaitToFinish( $o_dir, @notifyNames );

$logger->info( "Genotyping is done for all the given bam files." );
undef $!;

# parse all the mpileup.alleledepth files in $o_dir and merge the output into a single VCF
my $vcf_out = IO::File->new( $o_vcf, ">" ) or die logger->fatal( "Failed to create file $o_vcf." );

# vcf headers and comments
$vcf_out->print( "##fileformat=VCFv4.2\n" );
$vcf_out->print( "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">\n" );
$vcf_out->print( "##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Depth matching reference (REF) allele\">\n" );
$vcf_out->print( "##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Depth matching alternate (ALT) allele\">\n" );
$vcf_out->print( "##FORMAT=<ID=VF,Number=1,Type=Float,Description=\"Variant frequence (AD/DP)\">\n" );
$vcf_out->print( "##FORMAT=<ID=DPP,Number=1,Type=Integer,Description=\"Depth on postitive strand\">\n" );
$vcf_out->print( "##FORMAT=<ID=DPN,Number=1,Type=Integer,Description=\"Depth on negative strand\">\n" );
$vcf_out->print( "##FORMAT=<ID=RDP,Number=1,Type=Integer,Description=\"Reference depth on postitive strand\">\n" );
$vcf_out->print( "##FORMAT=<ID=RDN,Number=1,Type=Integer,Description=\"Reference depth on negative strand\">\n" );
$vcf_out->print( "##FORMAT=<ID=ADP,Number=1,Type=Integer,Description=\"Alternate depth on postitive strand\">\n" );
$vcf_out->print( "##FORMAT=<ID=ADN,Number=1,Type=Integer,Description=\"Alternate depth on negative strand\">\n" );
$vcf_out->print( "##INFO=<analysis_type=mpileup.alleledepth_files_merged_into_vcf,comma_separated_input_file_names=[".join( ",", @alleledepth_file_names )."]>\n" );
$vcf_out->print( "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t".join( "\t", @bam_sample_names )."\n");

my %variant_metrics = (); # hash of all unique variants and their coordinates as keys and their correspoding allelic info as values
while( my $ad_file_path = shift @alleledepth_files )	{
	chomp $ad_file_path;
	die $logger->fatal( "ERROR: $ad_file_path does not exist." ) unless ( -e $ad_file_path );
	my $ad_file = IO::File->new( $ad_file_path ) or logger->fatal( "Could not open file $ad_file_path" );
	chomp( my $header = $ad_file->getline );
	die $logger->fatal( "$ad_file is not in the proper mpileup.alleledepth file format." ) unless ( $header =~ m/^Ref_BAM\tSample\tChrom\tPOS\tRef\tAlt\tTotal_Depth\tRef_Depth\tAlt_Counts\tAlt_Freq\tRef_Forward\tRef_Reverse\tAlt_Forward\tAlt_Reverse$/ );
	$logger->info( "Gathering allele depth and vf info from $ad_file_path ..." );
	while( my $line = $ad_file->getline )	{
		chomp $line;
		my( $variant_info, $allelic_info ) = ad_to_vcf( $line );
		( $variant_metrics{$variant_info} ) ? ( $variant_metrics{$variant_info} .= "\t$allelic_info" ) : ( $variant_metrics{$variant_info} = "$allelic_info" );
	}
	$ad_file->close;
}

my $allelic_format = qw( DP:RD:AD:VF:DPP:DPN:RDP:RDN:ADP:ADN );
my @variant_infos = sort keys %variant_metrics;

# print the data in %variant_metrics to the .vcf file
while ( my $curr_variant = shift @variant_infos ) {
	$vcf_out->print( join( "\t", $curr_variant, ".\t.\t.", $allelic_format, $variant_metrics{$curr_variant} )."\n" );
}

$vcf_out->close;

$logger->info( "Completed merging all the relevant alleledepth files into a vcf." );
undef %variant_metrics;

# Execute housekeeping, print runtime, and exit message.
Housekeeping( $o_dir, @deletefiles );
$now = time - $now;
$now = sprintf( "%02d:%02d:%02d", int( $now / 3600 ), int( ( $now % 3600 ) / 60) , int( $now % 60 ) );
$logger->info( "Total running time: $now\n" );
$logger->info( "Done! $0 has completed execution.\n" );

exit;
####################################################

# ad_to_vcf reformats the allelic info from the mpileup.alleledepth into the .vcf format 
sub ad_to_vcf {
	my( $variant_line ) = @_;
	my @var = split( /\t/, $variant_line);
	map{ s/^\s+|\s+$|\r|\n//g; $_ } @var;
	$var[9] = 0 if( ( @var == 14 ) and ( $var[9] =~ /^[0.]+$/ ) ); # Simplify vf values of 0.0000 to 0.
	my $variant_info = join( "\t", @var[2..3], ".", @var[4..5] );
	my $allelic_info = qw ( 0:0:0:0:0:0:0:0:0:0 ); # for variants in the input vcf that are not genotyped in a given bam file
	$allelic_info = join( ":", @var[6..9], $var[10]+$var[12], $var[11]+$var[13], @var[10..13] ) if( @var == 14 );
	return( $variant_info, $allelic_info );
}
#####################

# Make Notification file --- adapted from RunIlluminaProcess.pl
sub MakeCSH	{
	my( $outdir ) = @_;
	my $filename = $outdir."/Notify.csh";
	my $ntmp = IO::File->new( $filename, ">" ) or die $logger->warn( "Failed to create $filename." );
	if( $ntmp )	{
		$ntmp->print( "#!/bin/csh\n" );
		$ntmp->print( "#Notification File\n" );
		$ntmp->print( "echo", " This is Done", "\n" );
		$ntmp->close;
		eval { `chmod +x $filename`; };
		$logger->fatal( "Cannot change permissions for $filename. Error:$@" ) if ($@);
	}
	return;
}
#####################
# Waiting for the process to finish --- adapted from RunIlluminaProcess.pl
sub WaitToFinish{
	my ( $outdir, @waitfilenames ) = @_;
	$logger->info( "Waiting for the Process to finish..." );
	foreach my $wfile( @waitfilenames )	{
		next if ( $wfile eq "NULL" );
		sleep 10 while ( !( -e "$outdir/$wfile" ) );
		while( -e "$outdir/$wfile" )	{
			open( FH, "<", "$outdir/$wfile" );
			while( <FH> )	{
				if( $_ =~ /This is Done/ig )	{
					#print "\nFinished: $wfile\n";
					last;
				}else	{
					wait;
				}
			}
			last;
		}
		close(FH);
	}

	foreach my $wfile ( @waitfilenames )	{
		next if ( $wfile eq "NULL" );
		eval { `rm $outdir/$wfile`; };
		if ($@) { $logger->warn( "Cannot remove $outdir/$wfile. Error:$@" ); }
	}
	return;
}

######################
# Delete Unwanted Files --- adapted from RunIlluminaProcess.pl
sub Housekeeping	{
    my( $outdir, @list ) = @_;
    $logger->info( "Removing Unwanted files...\n" );
    foreach my $file ( @list )    {
    	$file = $outdir."/".$file;
    	if( -e $file )	{
    		eval{ `rm $file`; };
        	if( $@ )	{
        		$logger->fatal( "Cannot remove $outdir/$file. Error:$@" );
        	}
        }
    }
    return;
}

######################

__DATA__

=head1 NAME

dmp_adepth2vcf.pl v1.2 - Genotype one or more bams; reformat and merge allelic info of all variants from the the resulting multiple alleledepth files into a single vcf.

=head1 SYNOPSIS

 perl dmp_adepth2vcf.pl --help

 SGE: perl dmp_adepth2vcf.pl path_to_dmp_adepth2vcf.pl --i_vcf path_to_input_vcf_file.vcf --i_bam path_to_bams.txt --RefFile path_to_hg19.fasta --samtools path_to_samtools_binary --bedtools path_to_dir_with_bedtools --bsub path_to_bsub --queue example_queue --o_dir path_to_out_dir --identifier name_for_output_vcf
 
 LSF: perl dmp_adepth2vcf.pl path_to_dmp_adepth2vcf.pl --g_allele path_to_dmp_genotype_allele.pl --i_vcf path_to_input_vcf_file.vcf --i_bam path_to_bams.txt --RefFile path_to_hg19.fasta --samtools path_to_samtools_binary --bedtools path_to_dir_with_bedtools --bsub path_to_bsub --queue example_queue --o_dir path_to_out_dir --identifier name_for_output_vcf


=head1 OPTIONS

 --i_vcf          Path to input vcf file in vcf4.2 format. vcf4.2 file with eight columns is also accepted. (required)
 --o_dir     	  Path to output directory for all the output files of this program (optional)
 --o_vcf          Path to output multi-sample VCF [<o_dir>/allele_depth_merged.vcf] (optional)
 --g_allele       Path to dmp_genotype_allele.pl script DMP IMPACT pipeline [/dmp/resources/prod/software/dmp-impact-res/production/bin/dmp_genotype_allele.pl] (optional)
 --samtools       Path to samtools [/dmp/resources/prod/tools/bio/samtools/production/samtools] (optional)
 --bedtools       Path to bedtools [/home/shahr2/Software/BEDTools/current/bin/] (optional)
 --perl           Path to the perl binary that you want to use to run this program [/usr/bin/env perl] (optional)
 --qsub           Path to qsub executable for SGE (required if --bsub is not specified)
 --bsub           Path to qsub executable for LSF (required if --qsub is not specified)
 --queue          Name of the SGE / LSF queue where this software should run (required)
 --i_bam          Comma-separated absolute paths to bam files or a text file, with a .txt suffix, of all the absolute paths to bam files with line breaks (required)
 --identifier     Unique identifier (eg. Patient MRN) to be used in the name of the output merged.vcf (required unless --o_vcf is specified)
 --RefFile        Path to reference genome file [/home/shahr2/Resources/IMPACT/References/Homo_sapiens_assembly19_formatted.fasta] (optional)
 --mbq            Min. Base Quality Threshold to be set for genotyping the bam files [default:5] (optional)
 --mmq            Min. Mapping Quality Threshold to be set for genotyping the bam files [default:5] (optional)
 --delete         2=>To delete files 1=> To keep files [default:2 --- stdout, stderr, stat, and csh files generated by this script are deleted if the script completes without any errors] (optional)
 --help           Print help message and quit (optional)
 --man            Print detailed message and usage instructions and quit (optional)

=head1 DESCRIPTION

 dmp_adepth2vcf.pl v1.2. This script genotypes the bam files provided using the MSK-IMPACT dmp_genotype_allele.pl script and merges all the resulting mpileup.alleledepth into a single VCF, in which each line contains a unique variant with allele depth and variant frequency information for all the samples.

 Example usage on SGE or LSG:
 Sample call for SGE:
 qsub_path -q queue_name -V -N jobname -wd working_dir_path -e stderr_name -o stdout_name -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y perl_path path_to_dmp_adepth2vcf.pl --i_vcf path_to_input_vcf_file.vcf --i_bam path_to_bams.txt --RefFile path_to_hg19.fasta --samtools path_to_samtools_binary --bedtools path_to_dir_with_bedtools --bsub path_to_bsub --queue example_queue --o_dir path_to_out_dir --identifier name_for_output_vcf
 
 Sample call for LSF:
 bsub_path -q queue_name -J jobname -cwd working_dir_path -e stderr_name -o stdout_name -We 24:00 -R "rusage[mem=2]" -M 4 -n 1 "perl_path path_to_dmp_adepth2vcf.pl --g_allele path_to_dmp_genotype_allele.pl --i_vcf path_to_input_vcf_file.vcf --i_bam path_to_bams.txt --RefFile path_to_hg19.fasta --samtools path_to_samtools_binary --bedtools path_to_dir_with_bedtools --bsub path_to_bsub --queue example_queue --o_dir path_to_out_dir --identifier name_for_output_vcf"

 
=head1 AUTHORS

 Gowtham Jayakumaran ( jayakumg@mskcc.org )

=cut
