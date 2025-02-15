#!/usr/bin/perl
use strict 'vars';
use vars qw(@files %o_geno $f $i $j $dirname $good $lastname  @sort_snp $last_chrom $no_snps %data @fields @subfields %rsname %pos %chrom $snp_count %fam %indiv %mom %dad %affect %gender 
%good_sample $snp_comp $sample_comp $total_snps @fields @files @filess $first_sample $last_sample %mz_name_to_ha_name );

$dirname = ".";
use 5.10.0;

if(@ARGV != 3) 
{
  print "\n Usage: ${0} annotation_file bimfile outname \n\n "; 
  exit(1);
}
#pedigree info file
open(FILE,"$ARGV[0]") || die "\nCan't open file $ARGV[0] which should contain annotations\n";
$_ = <FILE>;
chomp;
@fields = split('\t');
my $s_id = 0;
my $s_alleles = 0;
my $ref_id = 0;
my $vcf_pos = 0;
my $alt_id = 0;
for($i=0;$i<@fields;$i++)
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
  $_ = $fields[$ref_id]; s/\|//g; $fields[$ref_id] = $_;
  $_ = $fields[$alt_id]; 
  if(/^\+/)
  {
	my $j = substr $fields[$ref_id], 0, 1;
	$fields[$ref_id] = $j;
	s/\|//g;  s/\+//g; $_ = "$j$_" ;$fields[$alt_id] = $_;
  }
  my $this = "$fields[0]\_$fields[$vcf_pos]\_$fields[$ref_id]\_$fields[$alt_id]";
  $this = "$fields[0]\_$fields[$vcf_pos]";
  $_ = $this;
  if(!/^chr/)
  {
	$this = "chr$this";
  }
  $_ = $fields[$s_id];
  if(/^rs/)
  {
	$_ = $fields[$s_id];
	my @sfields = split('[\|;]');
	$rsname{$this} = $sfields[0];
  }
  else
  {
	my @sfields = split('_',$this);
	$rsname{$this} = "$sfields[0]\_$sfields[1]";
  }
  #print "\n $this $rsname{$this}";
}
close(FILE);

open(oFILE,">$ARGV[2]") || die "\nCan not open $ARGV[2] for writing\n"; 
open(FILE,"$ARGV[1]") || die "\nCan not open $ARGV[1] for reading which should contain the pvar file\n"; 
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
while(<FILE>)
{
	chomp;
	@fields = split('\t');
	my $j = length $fields[5];
  if($fields[0] > 22)
  {
	if($fields[0] == 23)
	{
		$fields[0] = "X";
	}
	elsif($fields[0] == 24)
	{
		$fields[0] = "Y";
	}
	elsif($fields[0] == 25)
	{
		$fields[0] = "M";
	}
  }
	# my $this = "$fields[0]\_$fields[3]\_$fields[4]\_$fields[5]";
	my $this = "$fields[0]\_$fields[1]";
	$_ = $this;
	if(!/^chr/)
	{
		$this = "chr$this";
	}
	#print "\n In BIM $this  .$rsname{$this}.";
	my $name = $this;
	if(exists($rsname{$this}))
	{
		$name = $rsname{$this};
	}
	$fields[2] = $name;
 	local $" = "\t";
	print oFILE "@fields\n";
}
