#!/usr/bin/perl
use IO::Compress::Gzip qw($GzipError);
use IO::Uncompress::Gunzip qw($GunzipError);
use strict 'vars';
use vars qw(@fields $f1_len $f2_len);

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
$f1_len = @fields;
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
$f2_len = @fields;
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
my $old_c2 = $f2[0];
while( ((length $l1) > 20) || ((length $l2) > 20) )
{
	my $pos1 = $f1[1];
	my $pos2 = $f2[1];
	my $c1 = $f1[0];
	my $c2 = $f2[0];
	$f1[7] = ".";
	$f2[7] = ".";
	$l1 = $l2 = "";
	if($c1 eq $c2) 
	{
		if($pos1 == $pos2) 
		{
			my @al;
			my %f2_map;
			my %f1_map;
			$f1_map{'.'} = '.';
			$f2_map{'.'} = '.';
			$f1_map{0} = 0;
			$f2_map{0} = 0;
			my $ac = 1;
			my $suff1 = "";
			my $suff2 = "";
			if( (length $f1[3]) >= (length $f2[3]) )
			{
		        	$al[0] = $f1[3];
				if($f2[3] ne $f1[3])
				{
					my $jj = length $f2[3];
					$suff2 = substr $f1[3],$jj;
					
				}
			}
			else
			{
				$al[0] = $f2[3];
				my $jj = length $f1[3];
				$suff1 = substr $f2[3],$jj;
			}
			my @ssf1 = split("\,",$f1[4]);
			for(my $j = 0; $j<@ssf1;$j++)
			{
				my $found = 0;
				$ssf1[$j] .= $suff1;
				for(my $i = 0; $i < $ac; $i++)
				{
					if($al[$i] eq $ssf1[$j])
					{
						$f1_map{$j+1} = $i;
						$found = 1;
						$i = $ac+1;
					}
				}	
				if(!$found)
				{
					$al[$ac] = $ssf1[$j];
					$f1_map{$j+1} = $ac++;
				}
			}
			my @ssf2 = split("\,",$f2[4]);
			for(my $j=0;$j<@ssf2;$j++)
			{
				my $found = 0;
				$ssf2[$j] .= $suff2;
				for(my $i = 0; $i < $ac; $i++)
				{
					if($al[$i] eq $ssf2[$j])
					{
						$f2_map{$j+1} = $i;
						$found = 1;
						$i = $ac+1;
					}
				}	
				if(!$found)
				{
					$al[$ac] = $ssf2[$j];
					$f2_map{$j+1} = $ac++;
				}
			}
			$f1[3] = $f2[3] = $al[0];
			$f1[4] = $al[1];
			for(my $i = 2;$i < @al;$i++)
			{
				$f1[4] .= ",$al[$i]";
			}
			print $f1[0];
			for(my $i=1;$i<9;$i++)
			{
				print "\t$f1[$i]";
			}
			if($ac < 10)
			{
				for(my $i=9;$i<$f1_len;$i++)
				{
					my $all_1 = substr $f1[$i],0,1;
					my $all_2 = substr $f1[$i],2,1;
					print "\t$f1_map{$all_1}/$f1_map{$all_2}";
				}
				for(my $i=9;$i<$f2_len;$i++)
				{
					my $all_1 = substr $f2[$i],0,1;
					my $all_2 = substr $f2[$i],2,1;
					print "\t$f2_map{$all_1}/$f2_map{$all_2}";
				}
			}
			else
			{
				for(my $i=9;$i<$f1_len;$i++)
				{
					my @ssf = split('[\\\/|]',$f1[$i]);
					my $all_1 = $ssf[0];
					my $all_2 = $ssf[1];
					print "\t$f1_map{$all_1}/$f1_map{$all_2}";
				}
				for(my $i=9;$i<$f2_len;$i++)
				{
					my @ssf = split('[\\\/|]',$f2[$i]);
					my $all_1 = $ssf[0];
					my $all_2 = $ssf[1];
					print "\t$f2_map{$all_1}/$f2_map{$all_2}";
				}
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
			for(my $i=1;$i<$f1_len;$i++)
			{
				print "\t$f1[$i]";
			}
			for(my $i=9;$i<$f2_len;$i++)
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
			for(my $i=9;$i<$f1_len;$i++)
			{
				print "\t0/0";
			}
			for(my $i=9;$i<$f2_len;$i++)
			{
				print "\t$f2[$i]";
			}
			print "\n";
			$l2 =  $fh2->getline();
			chomp $l2;
			@f2 = split('\t',$l2);
		}
	}
	elsif( ($c1 le $c2) && (@f1 > 10) )
	{
		print $f1[0];
		for(my $i=1;$i<$f1_len;$i++)
		{
			print "\t$f1[$i]";
		}
		for(my $i=9;$i<$f2_len;$i++)
		{
			print "\t0/0";
		}
		print "\n";
		$l1 =  $fh1->getline();
		chomp $l1;
		@f1 = split('\t',$l1);
	}
	elsif(@f2 > 10)
	{
		print $f2[0];
		for(my $i=1;$i<9;$i++)
		{
			print "\t$f2[$i]";
		}
		for(my $i=9;$i<$f1_len;$i++)
		{
			print "\t0/0";
		}
		for(my $i=9;$i<$f2_len;$i++)
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

