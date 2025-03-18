#!/usr/bin/perl
use IO::Compress::Gzip qw($GzipError);
use IO::Uncompress::Gunzip qw($GunzipError);
use strict 'vars';
use vars qw(@fields @f1_temp @f2_temp);

if(@ARGV != 3)
{
	die "\n Usage: ${0} merge_header first_vcffile second_vcffile\n";
}
open(FILE,"$ARGV[0]") || die "\n Can not open $ARGV[0], which should be the merged_header for reading\n";
while(<FILE>)
{
	print ;
}
close(FILE);

my $fh1 = IO::Uncompress::Gunzip->new($ARGV[1]) or die "For file $ARGV[1] gzip failed: $GunzipError\n";
$_ = $fh1->getline();
while(!/^#CHROM	POS/)
{
	$_ = $fh1->getline();
}
chomp;
@fields = split('\t');
@f1_temp = @fields;
print "$fields[0]";
for(my $i =1;$i < 9;$i++)
{
	print "\t$fields[$i]";
}
for(my $i =9;$i < @fields;$i++)
{
	print "\t$fields[$i]";
}
my $fh2 = IO::Uncompress::Gunzip->new($ARGV[2]) or die "For file $ARGV[2] gzip failed: $GunzipError\n";
$_ = $fh2->getline();
while(!/^#CHROM	POS/)
{
	$_ = $fh2->getline();
}
chomp;
@fields = split('\t');
@f2_temp = @fields;
for(my $i =9;$i < @fields;$i++)
{
	print "\t$fields[$i]";
}
print "\n";
my $l1 =  $fh1->getline();
chomp $l1;
my $l2 =  $fh2->getline();
chomp $l2;
my @f1 = split('\t',$l1);
my @f2 = split('\t',$l2);
my $old_c1 = $f1[0];
my $old_c2 = $f1[0];
while( (@f1 > 9) || (@f2 > 9) )
{
	my $pos1 = $f1[1];
	my $pos2 = $f2[1];
	my $c1 = $f1[0];
	my $c2 = $f2[0];
	$f1[7] = ".";
	$f2[7] = ".";
	if($c1 eq $c2) 
	{
		if ($pos1 == $pos2) 
		{
			my @al;
			my @f2_map;
			$f2_map['.'] = '.';
		        $al[0] = $f1[3];
			my $ac = 1;
			my @ssf1 = split("\,",$f1[4]);
			for(my $i = 0; $i<@ssf1;$i++)
			{
				$al[$ac++] = $ssf1[$i];
			}
			my $found = 0;
			for(my $i = 0; $i < $ac; $i++)
			{
				if($al[$i] eq $f2[3])
				{
					$f2_map[0] = $i;
					$found = 1;
					$i = $ac+1;
				}
			}	
			if(!$found)
			{
				$al[$ac] = $f2[3];
				$f2_map[0] = $ac++;
			}
			my @ssf2 = split("\,",$f2[4]);
			for(my $j=0;$j<@ssf2;$j++)
			{
				my $found = 0;
				for(my $i = 0; $i < $ac; $i++)
				{
					if($al[$i] eq $ssf2[$j])
					{
						$f2_map[$j+1] = $i;
						$found = 1;
						$i = $ac+1;
					}
				}	
				if(!$found)
				{
					$al[$ac] = $ssf2[$j];
					$f2_map[$j+1] = $ac++;
				}
			}
			$f1[4] = $al[1];
			for(my $i = 2;$i < @al;$i++)
			{
				$f1[4] .= ",$al[$i]";
			}
			print $f1[0];
			for(my $i=1;$i<@f1;$i++)
			{
				print "\t$f1[$i]";
			}
			for(my $i=9;$i<@f2;$i++)
			{
				my @ss = split("[\|\\\/]",$f2[$i]);
				print "\t$f2_map[$ss[0]]/$f2_map[$ss[1]]";
			}
			print "\n";


			$l1 =  $fh1->getline();
			chomp $l1;
			@f1 = split('\t',$l1);
			$l2 =  $fh2->getline();
			chomp $l2;
			@f2 = split('\t',$l2);
		}
		elsif($pos1 < $pos2)
		{
			print $f1[0];
			for(my $i=1;$i<@f1;$i++)
			{
				print "\t$f1[$i]";
			}
			for(my $i=9;$i<@f2_temp;$i++)
			{
				print "\t0/0";
			}
			print "\n";
				
			$l1 =  $fh1->getline();
			chomp $l1;
			@f1 = split('\t',$l1);
		}
		else
		{
			print $f2[0];
			for(my $i=1;$i<9;$i++)
			{
				print "\t$f2[$i]";
			}
			for(my $i=9;$i<@f1_temp;$i++)
			{
				print "\t0/0";
			}
			for(my $i=9;$i<@f2;$i++)
			{
				print "\t$f2[$i]";
			}
			print "\n";
			$l2 =  $fh2->getline();
			chomp $l2;
			@f2 = split('\t',$l2);
		}
		$old_c1 = $c1;
		$old_c2 = $c2;
	}
	elsif($c1 eq $old_c1)
	{
		print $f1[0];
		for(my $i=1;$i<@f1;$i++)
		{
			print "\t$f1[$i]";
		}
		for(my $i=9;$i<@f2_temp;$i++)
		{
			print "\t0/0";
		}
		print "\n";
		$l1 =  $fh1->getline();
		chomp $l1;
		@f1 = split('\t',$l1);
	}
	else
	{
		print $f2[0];
		for(my $i=1;$i<9;$i++)
		{
			print "\t$f2[$i]";
		}
		for(my $i=9;$i<@f1_temp;$i++)
		{
			print "\t0/0";
		}
		for(my $i=9;$i<@f2;$i++)
		{
				print "\t$f2[$i]";
		}
		print "\n";
		$l2 =  $fh2->getline();
		chomp $l2;
		@f2 = split('\t',$l2);
	}
}

