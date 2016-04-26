#!/usr/bin/env perl

# allele_depth_to_vcf - Reformat and merge variants from multiple alleledepth files into a single vcf.
# version: v1.0
# Author: Gowtham Jayakumaran ( jayakumg <at> mskcc <dot> org )
# Date: 25/Apr/2016
# Last Modified: 25/Apr/2016

use strict;
use warnings;
use IO::File;
use Getopt::Long qw( GetOptions );
use Pod::Usage qw( pod2usage );
use feature qw( say );
use Cwd qw( abs_path );


# parse options and print error, if any, or print usage, if help requested.
my( $help, $man, $keep, $i_vcf, $g_allele, $samtools, $i_bam, $o_dir, $o_vcf );
if( @ARGV < 1 or
	$ARGV[0] !~ m/^-/ or
	!GetOptions(
		'help!' => \$help,
		'man!' => \$man,
		'keep!' => \$keep,
		'i_vcf=s' => \$i_vcf,
		'g_allele=s' => \$g_allele,
		'samtools=s' => \$samtools,
		'i_bam=s' => \$i_bam,
		'o_dir=s' => \$o_dir,
		'o_vcf=s' => \$o_vcf
		) )	{
	pod2usage( -verbose => 1, -input => \*DATA, -exitval => 2 );
}
pod2usage( -verbose => 1, -input => \*DATA, -exitval => 0 ) if( $help );
pod2usage( -verbose => 2, -input => \*DATA, -exitval => 0 ) if( $man );

# set default paths and constants
$g_allele = "/dmp/resources/prod/software/dmp-impact-res/VERSIONS/dmp-impact-res-v118/bin/dmp_genotype_allele.pl" unless defined $g_allele;
$samtools = "/dmp/resources/prod/tools/bio/samtools/production/samtools" unless defined $samtools;
$samtools = `which samtools` unless defined $samtools and ( -e $samtools );
chomp ( $samtools );
( $samtools and -e $samtools ) or die "ERROR: Please install samtools and include it in your PATH.\n";
die "ERROR: Please provide the path for your genotype script. Try perl $0 --help for more info.\n" unless ( -e $g_allele ) ;
die "ERROR: Please provide a vcf file. Try perl $0 --help for more info.\n" unless defined $i_vcf and ( -e $i_vcf ) ;
die "ERROR: Please provide a proper argument for --i_bam. Try perl $0 --help for more info.\n" unless defined $i_bam and ( -e $i_bam ) ;

# set output dir and files
$o_dir = `pwd` unless defined $o_dir and ( -e $o_dir );
chomp $o_dir;

# set output .vcf file
if( defined $o_vcf )	{
	$o_vcf =~ s/(\.)?(maf|tsv|csv|txt)?$/.vcf/;
	$o_vcf .= ".vcf" unless( $o_vcf =~ m/(\.)?(vcf)?$/ );
}else	{
	$o_vcf = $o_dir."/allele_depth_merged.vcf" unless defined $o_vcf;
}

# obtain absolute paths for all the input and output files
for my $file ( $g_allele, $i_vcf, $o_vcf, $i_bam )	{
	$file = abs_path( $file );
}

# parse the arguments to --i_bam depending on whether its a comma-separated string of bam file paths or a text file with bam file paths;
my @bam = ();
if ($i_bam =~ /.txt$/ )	{
	my $bam_file = IO::File->new( $i_bam ) or die "ERROR: Cannot open file $i_bam\n";
	while ( <$bam_file> )	{
		next if ( $_ =~ /^#/ );
		$_ =~ s/^\s+|\s+$|\r|\n//g;
		push @bam, $_;
	}
	$bam_file->close;
}else {
	@bam = map{ s/^\s+|\s+$|\r|\n//g; $_ } split( /,/, $i_bam );
}

# check whether the vcf file provided is in the proper format. A more stringent check can be implemented, if necessary.
my $vcf_header = `grep '^#' $i_vcf`;
$vcf_header =~ s/\n//g;
map{ ( m/^\#\#fileformat=VCFv4.2.+\#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT/ ) or 
	die "ERROR: $i_vcf is not in the vcf4.2 format. Please check https://samtools.github.io/hts-specs/VCFv4.2.pdf for more info.\n" } $vcf_header;

# Extracting sample names from bam file, to be included in the VCF
my @bam_sample_name = ();
foreach my $bam_file ( @bam )	{
	die "ERROR: $bam_file does not exist.\n" unless ( -e $bam_file );
	my $header = `$samtools view -H $bam_file`;
	map{ s/\n//g; s/.+\@RG\sID:(\S+)\s.+/$1/; } $header;
	push @bam_sample_name, $header;
}

# get the absolute path and bam file names of all mpileup.alleledepth output files in $o_dir
my @alleledepth_files = @bam;
my @alleledepth_file_names = @bam;
map{ s/(.+\/)?/$o_dir\//; s/(\.bam)?$/_mpileup.alleledepth/; } @alleledepth_files;
map{ s/(.+\/)?//; s/(\.bam)?$/_mpileup.alleledepth/; } @alleledepth_file_names;


# run dmp_genotype_allele.pl on each of the input bam files and direct the ouput mpileup.alleledepth files to $o_dir
while ( my $bam_file = shift @bam )	{
	die "ERROR: $bam_file does not exist.\n" unless ( -e $bam_file );
	say "Genotyping $bam_file ...";
	my $exitstatus = system ( "/usr/bin/env perl $g_allele --outdir $o_dir -fmv $i_vcf -bam $bam_file --qsub /common/sge/bin/lx-amd64/qsub -q test.q 1>/dev/null" );
	die "ERROR: $g_allele did not execute successfully. Please check your arguments and path of your input files.\n" unless( $exitstatus == 0 );
}


# parse all the mpileup.alleledepth files in $o_dir and merge the output into a single VCF
my $vcf_out = IO::File->new( $o_vcf, ">" ) or die "ERROR: Failed to create file $o_vcf\n";

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
$vcf_out->print( "##FILE_TYPE=<analysis_type=merged_mpileup.alleledepth_files input_file=[".join( ", ", @alleledepth_file_names )."]>\n" );
$vcf_out->print( "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t".join( "\t", @bam_sample_name )."\n");

my %variant_metrics = (); # hash of all unique variants and their coordinates as keys and their correspoding allelic info as values
while( my $ad_file_path = shift @alleledepth_files )	{
	chomp $ad_file_path;
	die "ERROR: $ad_file_path does not exist.\n" unless ( -e $ad_file_path );
	my $ad_file = IO::File->new( $ad_file_path ) or die "ERROR: Cannot open file: $ad_file_path\n";
	chomp( my $header = $ad_file->getline );
	die "ERROR: $ad_file is not in the proper mpileup.alleledepth file format" unless ( $header =~ m/^Ref_BAM\tSample\tChrom\tPOS\tRef\tAlt\tTotal_Depth\tRef_Depth\tAlt_Counts\tAlt_Freq\tRef_Forward\tRef_Reverse\tAlt_Forward\tAlt_Reverse$/ );
	say "Parsing $ad_file_path ...";
	while( my $line = $ad_file->getline )	{
		chomp $line;
		my( $variant_info, $allelic_info ) = ad_to_vcf( $line );
		( $variant_metrics{$variant_info} ) ? ( $variant_metrics{$variant_info} .= "\t$allelic_info" ) : ( $variant_metrics{$variant_info} = "$allelic_info" );
	}
	$ad_file->close;
	system( "rm $ad_file_path $o_dir/*stderr $o_dir/*stat $o_dir/Notify.csh $o_dir/Impact_local.log 1>/dev/null 2>/dev/null" ) unless defined $keep;
}

my $allelic_format = qw( DP:RD:AD:VF:DPP:DPN:RDP:RDN:ADP:ADN );
my @variant_infos = sort keys %variant_metrics;

# print the data in %variant_metrics to the .vcf file
while ( my $curr_variant = shift @variant_infos ) {
	$vcf_out->print( join( "\t", $curr_variant, ".\t.\t.", $allelic_format, $variant_metrics{$curr_variant} )."\n" );
}

$vcf_out->close;

undef %variant_metrics;

# function ad_to_vcf to reformats the allelic info from the mpileip.alleledepth into the .vcf format 
sub ad_to_vcf {
	my ( $variant_line ) = @_;
	my @var = split( /\t/, $variant_line);
	map{ s/^\s+|\s+$|\r|\n//g; $_ } @var;
	$var[9] = 0 if( ( @var == 14 ) and ( $var[9] =~ /^[0.]+$/ ) );
	my $variant_info = join( "\t", @var[2..3], ".", @var[4..5] );
	my $allelic_info = qw ( 0:0:0:0:0:0:0:0:0:0 ); # for variants in the input vcf that are not genotyped in a given bam file
	$allelic_info = join( ":", @var[6..9], $var[10]+$var[12], $var[11]+$var[13], @var[10..13] ) if( @var == 14 );
	return( $variant_info, $allelic_info );
}

__DATA__

=head1 NAME

allele_depth_to_vcf - Reformat and merge variants from multiple alleledepth files into a single vcf.

=head1 SYNOPSIS

 perl allele_depth_to_vcf.pl --help
 perl allele_depth_to_vcf.pl --i_bam list_of_bams_files.txt --i_vcf file.vcf --keep

=head1 OPTIONS

 --i_vcf          Path to input vcf file in vcf4.2 format. vcf4.2 file eight columns is also accepted. (required)
 --o_dir     	  Path to output directory for all the output files of this program (optional)
 --o_vcf          Path to output multi-sample VCF [<o_dir>/allele_depth_merged.vcf] (optional)
 --g_allele       Path to dmp_genotype_allele.pl script DMP IMPACT pipeline [/dmp/resources/prod/software/dmp-impact-res/production/bin/dmp_genotype_allele.pl] (optional)
 --samtools       Path to samtools [/dmp/resources/prod/tools/bio/samtools/production/samtools] (optional)
 --i_bam          Comma-separated absolute paths to bam files or a text file, with a .txt suffix, of all the absolute paths to bam files with line breaks (required)
 --keep           Keep intermediate mpileup.alledepth and other files (optional)
 --help           Print help message and quit (optional)
 --man            Print detailed message and usage instructions and quit (optional)

=head1 DESCRIPTION

 allele_depth_to_vcf v1.0. This script genotypes the bam files provided using the MSK-IMPACT dmp_genotype_allele.pl script and merges all the resulting mpileup.alleledepth into a single VCF, in which each line contains a unique mutation variant with allele depth and variant frequency information for all the bam files.
 
=head1 AUTHORS

 Gowtham Jayakumaran ( jayakumg <at> mskcc <dot> org )

=head1 LICENSE

 Apache-2.0 | Apache License, Version 2.0 | https://www.apache.org/licenses/LICENSE-2.0

=cut
