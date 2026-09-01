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
while( ((length $l1) > 20) || ((length $l2) > 20) )
{
	my $pos1 = $f1[1];
	my $pos2 = $f2[1];
	my $c1 = $f1[0];
	my $c2 = $f2[0];
	$f1[7] = ".";
	$f2[7] = ".";
	$l1 = $l2 = "";
	my $v1 = $f1[2];
	my $v2 = $f2[2];
	if($v1 eq $v2)
	{
		print $f1[0];
		for(my $i=1;$i<$f1_len;$i++)
		{
			print "\t$f1[$i]";
		}
		for(my $i=9;$i<$f2_len;$i++)
		{
			print "\t$f2[$i]";
		}
		print "\n";
		$l1 =  $fh1->getline();
		chomp $l1;
		@f1 = split('\t',$l1);
		$l2 =  $fh2->getline();
		chomp $l2;
		@f2 = split('\t',$l2);
	}
	elsif( (@f1 > 8) && (($c1 < $c2) || ( ($c1 == $c2) && ($pos1 < $pos2)) || ( ($c1 == $c2) && ($pos1 == $pos2)  &&($v1 le $v2))))
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
	elsif(@f2 > 8)
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

