#!/usr/bin/perl
use strict 'vars';
use vars qw(@files %o_geno $f $i $j $dirname $good $lastname  @sort_snp $last_chrom $no_snps %data @fields @subfields %rsname %pos %chrom $snp_count %fam %indiv %mom %dad %affect %gender 
%good_sample $snp_comp $sample_comp $total_snps @fields @files @filess $first_sample $last_sample %mz_name_to_ha_name );

$dirname = ".";
use 5.10.0;

if(@ARGV != 2) 
{
  print "\n Usage: ${0} pvar_file outname \n\n "; 
  exit(1);
}

open(oFILE,">$ARGV[1]") || die "\nCan not open $ARGV[1] for writing\n"; 
open(FILE,"$ARGV[0]") || die "\nCan not open $ARGV[0] for reading which should contain the pvar file\n"; 
my $still_header = 1;
while($still_header)
{
	$_ = <FILE>;
	if(/^#CHROM/)
	{
		$still_header = 0;
	}
	print oFILE;
}
# 1-55051215-G-GA
while(<FILE>)
{
	chomp;
	@fields = split('\t');
	$fields[2] = "$fields[0]\-$fields[1]\-$fields[3]\-$fields[4]";
 	local $" = "\t";
	print oFILE "@fields\n";
}
