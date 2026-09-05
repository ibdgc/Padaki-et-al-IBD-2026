#!/usr/bin/perl
use IO::Zlib;
use strict 'vars';
use vars qw(@fields %keep_snp @keep_which);

if(@ARGV != 6)
{
	die "\n Usage: ${0} ANNOTATION_FILE PVAR_FILE MAX_FREQ MIN_CADD REPLACEMENT_ONLY[y,n] LOF_ONLY[y,n]\n";
}
my $top_freq = 1.0 - $ARGV[2];
my $bottom_freq = $ARGV[2];

if( ($top_freq < 0.5) || ($bottom_freq < 1e-8) )
{
	die "\n Impossible Maximum frequency of $ARGV[2] \n";
}

my $REPLACE_only = 0;
my $LOF_only = 0;
$_ = $ARGV[4];
if(/[Yy]/)
{
	$REPLACE_only = 1;
}
elsif(!/[Nn]/)
{
	die "\n $ARGV[4] is not a yes or no answer to REPLACEMENT_ONLY \n";
}
$_ = $ARGV[5];
if(/[Yy]/)
{
	$LOF_only = 1;
}
elsif(!/[Nn]/)
{
	die "\n $ARGV[4] is not a yes or no answer to LOF_ONLY \n";
}



tie *FILE, 'IO::Zlib', "$ARGV[0]", "rb";
my $cadd_thres = $ARGV[3];
if($cadd_thres > 200)
{
	die "\n Silly CADD threshold of $ARGV[3] \n";
}
$_ = <FILE>;
chomp;
@fields = split('\t');

my $vcf_pos = -1;
my $cadd_no = -1;
my $gene_name_no = -1;
my $gene_symbol_no = -1;
my $rs_name_no = -1;
my $func1 = -1;
my $func2 = -1;
my $a_freq = -1;

chomp;
@fields = split('\t');
for(my $i = 2; $i < @fields;$i++)
{	
	my $j = $fields[$i];
	if($j eq "refSeq.siteType")
	{
		$func1 = $i;
	} 
	elsif($j eq "refSeq.exonicAlleleFunction")
	{
		$func2 = $i;
	}
	elsif($j eq "vcfPos")
	{
		$vcf_pos = $i;
	}
	elsif($j eq "nearest.refSeq.name")
	{
		$gene_name_no = $i;
	}
	elsif($j eq "nearest.refSeq.name2")
	{
		$gene_symbol_no = $i;
	}
	elsif($j eq "cadd")
	{
		$cadd_no = $i;
	}
	elsif($j eq "id")
	{
		$rs_name_no = $i;
	}
	elsif($j eq "gnomad.joint.AF_joint_afr")
	{
		$a_freq = $i;
	}
}
	
if($vcf_pos < 0) 
{
	die "\n Can't find vcfPos field \n";
}
if($cadd_no < 1) 
{
	die "\n Can't find cadd field \n";
}
if($gene_name_no < 1) 
{
	die "\n Can't find nearest gene name field \n";
}
if($gene_symbol_no < 1) 
{
	die "\n Can't find nearest gene symbol field \n";
}
if($rs_name_no < 1) 
{
	die "\n Can't find rs_name \n";
}
if($func1 < 1) 
{
	die "\n Can't find refseq.sitetype \n";
}
if($func2 < 1) 
{
	die "\n Can't find exonic function \n";
}
if($a_freq < 1) 
{
	die "\n Can't find allele frequencies\n";
}
while(<FILE>)
{
	chomp;
	@fields = split('\t');
	for(my $i=0;$i<@keep_which;$i++)
	{
		$keep_which[$i] = 0;
	}
	my $keep_it = 0;
	$_ = $fields[$cadd_no];
	my @sfields = split('[\|\;]');	
	if(@sfields+0 >= 1)
	{
		foreach my $i (@sfields)
		{
			if($i >= $cadd_thres)
			{
				#print "\n $fields[0] $fields[$vcf_pos] $fields[2] $i $fields[$cadd_no] with threshold $cadd_thres";
				$keep_it = 1;
			}
		}
	}
	else
	{
		if( ($fields[2] eq "DEL") || ($fields[2] eq "INS") )
		{
			$keep_it = 1;
		}
	}

	if($keep_it && ($LOF_only || $REPLACE_only) )
	{
		my @sfields = split('[\|\;]',"$fields[$func1]");
		# print "\n just split $fields[$func2] \n";
		$keep_it = 0;
		my $which = 0;
		foreach my $k (@sfields)
		{
			$_ = $k;
			#print "\n testing $k";
			if(/non[Ss]yn/)
			{
				if(!$LOF_only)
				{
					$keep_it = 1;
					$keep_which[$which] = 1;
				}
			}
			elsif(/[sS]top/)
			{
				$keep_it = 1;
				$keep_which[$which] = 1;
			}
			elsif(/[[sS]plice/)
			{
				$keep_it = 1;
				$keep_which[$which] = 1;
			}
			elsif(/indel-frameshift/)
			{
				$_ = $fields[$func1];
				my @sss = split('[\|\;]');
				foreach my $ss (@sss)
				{
					if($ss eq "exonic")
					{
						$keep_it = 1;
						$keep_which[$which] = 1;
					}		
				}
			}
			$which++;
		}
		my @sfields = split('[\|\;]',"$fields[$func2]");
		# print "\n just split $fields[$func2] \n";
		$which = 0;
		foreach my $k (@sfields)
		{
			$_ = $k;
			#print "\n testing $k";
			if(/non[Ss]yn/)
			{
				if(!$LOF_only)
				{
					$keep_it = 1;
					$keep_which[$which] = 1;
				}
			}
			elsif(/[sS]top/)
			{
				$keep_it = 1;
				$keep_which[$which] = 1;
			}
			elsif(/[[sS]plice/)
			{
				$keep_it = 1;
				$keep_which[$which] = 1;
			}
			elsif(/indel-frameshift/)
			{
				$_ = $fields[$func1];
				my @sss = split('[\|\;]');
				foreach my $ss (@sss)
				{
					if($ss eq "exonic")
					{
						$keep_it = 1;
						$keep_which[$which] = 1;
					}		
				}
			}
			$which++;
		}
		
	}
	elsif($keep_it)
	{
		
		my @sfields = split('[\;\|]',$fields[$gene_symbol_no]);
		for(my $i =0; $i < @sfields; $i++)
		{
			$keep_which[$i] = 1;
		}
	}

	if($keep_it)
	{
		my @sfields = split('[\|\;]',$fields[$a_freq]);
		foreach my $i (@sfields)
		{
			# print "\n $i $bottom_freq $top_freq\n";
			if( ($i > $bottom_freq) && ($i < $top_freq) )
			{
				$keep_it = 0;
			}  
		}
	}

	if($keep_it)
	{
		my $gene = "!!!!";
		$_ = $fields[$gene_symbol_no];
		my @sfields = split('[\;\|]');
		my %gc = undef;
		my $best = $gene;
		$gc{$best} = -1;
		my $which = 0;
		
		foreach my $j (@sfields)
		{
			if( ($j ne "\\") && ($keep_which[$which] > 0) )
			{
				$gc{$j}++;
			}
			$which++;
		}
		foreach my $j (keys %gc)
		{
			if($gc{$j} > $gc{$best})
			{
				$best = $j;
				$gene = $best;
			}
		}
		if($gc{$best} < 1)
		{
			$gene = "!!!!";
			$_ = $fields[$gene_name_no];
			@sfields = split('[\;\|]');
			%gc = undef;
			$best = $gene;
			$gc{$best} = -1;
			$which = 0;
			foreach my $j (@sfields)
			{
				if( ($j ne "\\") && ($keep_which[$which] > 0))
				{
					$gc{$j}++;
				}
				$which++;
			}
			foreach my $j (keys %gc)
			{
				if($gc{$j} > $gc{$best})
				{
					$best = $j;
					$gene = $best;
				}
			}
			if($gc{$best} < 1)
			{
				$gene = "NA";
			}
		}
		
		
		$_ = $fields[0];
		# if(/^chr/)
		# {
		#	s/^chr//;
		#	$fields[0] = $_;
		# }
		if( ($gene ne "NA") || ($gc{$best} > 0))
		{
			# print "\nThe snp appears to be $fields[0] $fields[$vcf_pos] $fields[$rs_name_no] $best\n";
			$keep_snp{"$fields[0]\_$fields[$vcf_pos]"} = $gene;
			$keep_snp{$fields[$rs_name_no]} = $gene;
			#print "\t\t\t$fields[0]\_$fields[$vcf_pos] $fields[$rs_name_no] $gene\n";
		}
	}
}
close(FILE);
open(FILE2,"$ARGV[1]") || die "\n Could not open $ARGV[1] which should be the pvar_file \n";
while(<FILE2>)
{
	@fields = split('\s');
	if(exists($keep_snp{$fields[2]}))
	{
		print "$keep_snp{$fields[2]}\t$fields[2]\n";
	}
}
