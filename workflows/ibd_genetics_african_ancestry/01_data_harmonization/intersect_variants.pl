#!/usr/bin/perl
use strict 'vars';
use vars qw(%vars %pos @fields %exon);

use 5.10.0;

if(@ARGV != 4) 
{
  print "\n Usage: ${0} pvar_file1_without_extension annotation_file1 pvar_file_no_exten anntation_file2 \n\n "; 
  exit(1);
}

my @var_names;
for(my $i = 0; $i < 3; $i+=2)
{
	open(FILE,"$ARGV[$i].pvar") || die "\nCan not open $ARGV[$i].pvar for reading which should contain the pvar file\n"; 
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
	my $k = $i;
	if($k > 0)
	{
		$k--;
	}
	my $pvar_count = 0;
	my $f = $ARGV[$i];
	print "\n Working on file $f \n";
	while(<FILE>)
	{
		chomp;
		@fields = split('\t');
		my $j = $fields[0] + 0;
		$pos{$fields[2]} = "$fields[0]\-$fields[1]";
		$var_names[$k][$pvar_count] = $fields[2];
		if( ($j >= 1 ) && ($j <= 22) && ($fields[3] ne "\*") && ($fields[4] ne "\*") )
		{
			$vars{$f}{$fields[2]} = 1;
		}
		else
		{
			$vars{$f}{$fields[2]} = 0;
		}
		$pvar_count++;
		if($pvar_count % 1000000 == 0)
		{
			print "\n\t\tRead $pvar_count lines";
			# close(FILE);
		}
		# print oFILE "$fields[2]\n";
	}
	close(FILE);
	open(FILE,"$ARGV[$i+1]") || die "\nCan not open $ARGV[$i+1] for reading which should contain annotation\n"; 
	$_ = <FILE>;
	print "\n Starting to read annotation file \n";
	$pvar_count = 0;
	while(<FILE>)
	{
		chomp;
		@fields = split('\t');
		$_ = $fields[0];
		s/^chr//;
		$fields[0] = $_;
		my $this = "$fields[0]\-$fields[4]";
		$_ = $fields[6];
		if(/exonic/)
		{
			$exon{$this} = 1;
		}
		$pvar_count++;
		if($pvar_count % 1000000 == 0)
		{
			print "\n\t\tRead $pvar_count annotation lines";
			# close(FILE);
		}
	}
	close(FILE);
}
print "\n Time to write output \n\n";
print "File\tTotal\tKept\tRemoved\tExon_Tot\tEx_Kept\tEx_Removed\n";
for(my $i = 0; $i < 3; $i+=2)
{
	my $f = $ARGV[$i];
	my $of = $ARGV[2];
	if($i == 2)
	{
		$of = $ARGV[0];
	}

	open(oFILE,">$f.keep.txt") || die "\nCan not open $f.keep.txt writing\n"; 
	open(oFILE2,">$f.remove.txt") || die "\nCan not open $f.remove.txt writing\n"; 
	my @tot_count;
	my @ex_count;
	$tot_count[0] = $tot_count[1] = $tot_count[2] = 0;
	$ex_count[0] = $ex_count[1] = $ex_count[2] = 0;
	my $k = $i;
	if($k > 0)
	{
		$k--;
	}
	foreach my $j (@{$var_names[$k]})
	{
		$tot_count[0]++;
		if(exists($vars{$of}{$j}))
		{
			$vars{$f}{$j} += $vars{$of}{$j};
		}
		if(exists($exon{$pos{$j}}))
		{
			$vars{$f}{$j} += $exon{$pos{$j}};
			$ex_count[0]++;
			if($vars{$f}{$j} >= 2)
			{
				$ex_count[2]++;
			}
			else
			{
				$ex_count[1]++;
			}
		}
		if($vars{$f}{$j} >= 2)
		{
			$tot_count[2]++;
			print oFILE "$j\n";
		}
		else
		{
			$tot_count[1]++;
			print oFILE2 "$j\n";
		}
	}
	print "$f\t$tot_count[0]\t$tot_count[2]\t$tot_count[1]";
	print "\t$ex_count[0]\t$ex_count[2]\t$ex_count[1]\n";
	close(oFILE);
	close(oFILE2);
}
