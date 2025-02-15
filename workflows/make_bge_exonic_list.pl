#!/usr/bin/perl
use strict 'vars';
use vars qw(%names @fields);

use 5.10.0;

if(@ARGV != 3) 
{
  print "\n Usage: ${0} pvar_file annotation_file_for_additional_exonic outname \n\n "; 
  exit(1);
}
open(oFILE,">$ARGV[2]") || die "\nCan not open $ARGV[2] for writing\n"; 
open(FILE,"$ARGV[0]") || die "\nCan not open $ARGV[0] for reading which should contain the pvar file\n"; 
my $still_header = 1;
while($still_header)
{
	$_ = <FILE>;
	if(/^#CHROM/)
	{
		$still_header = 0;
	}
	#print oFILE;
}
while(<FILE>)
{
	chomp;
	@fields = split('\t');
	$names{$fields[2]} = 1;
}
close(FILE);
open(FILE,"$ARGV[1]") || die "\nCan not open $ARGV[1] for reading which should contain the pvar file\n"; 
$_ = <FILE>;
chomp;
@fields = split('\t');
my $s_id = 0;
my $s_alleles = 0;
my $ref_id = 0;
my $vcf_pos = 0;
my $alt_id = 0;
for(my $i=0;$i<@fields;$i++)
{
  #print "\n trying $_";
  if($fields[$i] eq 'dbSNP.id')
  {
    $s_id = $i;
  }
  if($fields[$i] eq 'ref')
  {
    $ref_id = $i;
  }
  if($ref_id <= 0)
  {
  	if($fields[$i] eq 'inputRef')
  	{
	    $ref_id = $i;
  	}
  }
  if($fields[$i] eq 'alt')
  {
    $alt_id = $i;
  }
  if($fields[$i] eq 'vcfPos')
  {
    $vcf_pos = $i;
  }

}
print "\n the s_id is $s_id with ref = $ref_id and alt = $alt_id and vcf_pos = $vcf_pos\n";
while(<FILE>)
{
	chomp;
	@fields = split('\t');
	if($fields[6] eq 'exonic')
	{
		$names{$fields[$s_id]} = 1;
	  	$_ = $fields[$ref_id]; s/\|//g; $fields[$ref_id] = $_;
	  	$_ = $fields[$alt_id]; 
 	 	if(/^\+/)
  		{
			my $j = substr $fields[$ref_id], 0, 1;
			$fields[$ref_id] = $j;
			s/\|//g;  s/\+//g; $_ = "$j$_" ;$fields[$alt_id] = $_;
  		}
		my $this = "$fields[0]\_$fields[$vcf_pos]";
		$_ = $fields[$s_id];
		if(/^rs/)
		{
			$_ = $fields[$s_id];
			my @sfields = split('[\|;]');
			$names{$sfields[0]} = 1;
  		}
		else
		{
			$names{$this} = 1;
		}
	}
}
foreach my $i (keys %names)
{
	print oFILE "$i\n";
}
