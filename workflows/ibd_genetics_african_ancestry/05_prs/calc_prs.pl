#!/usr/bin/perl
use strict 'vars';
use Math::Gauss ':all';
use vars qw(@fields %sex %affect $prev $thres %sam_anc $pthres $dist_thres  $set_count);


if(@ARGV < 5) 
{
	print "\n Usage: ${0} ancestry_file phenotype_famfile effects_file annotation_file prevalence [pvalue_threshold] [dist_thres]\n";
	exit(1);
}

$pthres = 1.0;
if( ($ARGV[5] > 0) && ($ARGV[5] < 1.0))
{
	$pthres = $ARGV[5];
}
$dist_thres = 0;
if( ($ARGV[6] > 0) )
{
	$dist_thres = $ARGV[6];
}
if( ($ARGV[4] <= 0) || ($ARGV[4] >= 1.0))
{
	die "\n Prevalence is outside a reasonable range $ARGV[4] \n";
}
$prev = $ARGV[4];
$thres = inv_cdf(1.0 - $prev);

open(FILE,"$ARGV[0]") || die "\n Can't read ancestry file\n\n";;
$_ = <FILE>;
s/\"//g;
chomp;
@fields = split(',');
my $name_col = -1;
my $afr_col = -1;
my $amr_col = -1;
my $eas_col = -1;
my $eur_col = -1;
my $sas_col = -1;
for(my $i = 0; $i < @fields; $i++)
{
	# print "\n.$fields[$i].\n";
	$_ = $fields[$i];
	if(/Sample/)
	{
		# print "\n Found name column at $i \n";
		$name_col = $i;
	}
	elsif($fields[$i] eq "AFR")
	{
		$afr_col = $i;
	}
	elsif($fields[$i] eq "AMR")
	{
		$amr_col = $i;
	}
	elsif($fields[$i] eq "EAS")
	{
		$eas_col = $i;
	}
	elsif($fields[$i] eq "EUR")
	{
		$eur_col = $i;
	}
	elsif($fields[$i] eq "SAS")
	{
		$sas_col = $i;
	}
}
if($name_col < 0) { die "\n Can't find Sample ID \n";}
if($afr_col < 0) { die "\n Can't find afr col \n";}
if($amr_col < 0) { die "\n Can't find amr col \n";}
if($eas_col < 0) { die "\n Can't find eas col \n";}
if($eur_col < 0) { die "\n Can't find eur col \n";}
if($sas_col < 0) { die "\n Can't find sas col \n";}
while(<FILE>)
{
	s/\"//g;
	chomp;
	@fields = split(',');
	my $this = $fields[$name_col];
	$sam_anc{$this}{"afr"} = $fields[$afr_col];
	$sam_anc{$this}{"amr"} = $fields[$amr_col];
	$sam_anc{$this}{"eas"} = $fields[$eas_col];
	$sam_anc{$this}{"nfe"} = $fields[$eur_col];
	$sam_anc{$this}{"sas"} = $fields[$sas_col];
}
close(FILE);

my @super_pop = qw(afr amr eas nfe sas);

open(FILE,"$ARGV[1]") || die "\n Can't read phenotype file $ARGV[1]\n\n";
while(<FILE>)
{
	chomp;
	@fields = split('\s');
	$sex{$fields[1]} = $fields[4];
	$affect{$fields[1]} = $fields[5];
}
close(FILE);
open(FILE,"$ARGV[2]") || die "\n Can't read snp effects file $ARGV[2]\n\n";
$_ = <FILE>;
chomp;
@fields = split('\t');

my $snp_col = -1;
my $or_col = -1;
my $freq_col = -1;
my $pvalue_col = -1;
my $pos_col = -1;

my %or;
my %beta;
my %ref_freq;

for(my $i = 0; $i <@fields;$i++)
{
	# print "\n.$fields[$i].\n";
	if($fields[$i] eq "hm_rsid")
	{
		$snp_col = $i;
	}
	elsif($fields[$i] eq "hm_pos")
	{
		$pos_col = $i;
	}
	elsif($fields[$i] eq "hm_beta")
	{
		$or_col = $i;
	}
	elsif($fields[$i] eq "hm_effect_allele_frequency")
	{
		$freq_col = $i;
	}
	elsif($fields[$i] eq "p_value")
	{
		$pvalue_col = $i;
	}
}

if($snp_col < 0) { die "\n Can't find snp name col \n";}
if($or_col < 0) { die "\n Can't find OR name col \n";}
if($freq_col < 0) { die "\n Can't find Beta name col \n";}
if($pvalue_col < 0) { die "\n Can't find pvalue name col \n";}
if($pos_col < 0) { die "\n Can't find pos name col \n";}


my $last_pos = -$dist_thres;
$set_count = -1;
my @snp_sets;
my $this_set_count = 0;
my %pvalues;
# print "\n About to start reading effects \n\n";
while(<FILE>)
{
	chomp;
	@fields = split('\t');
	my $this_p = $fields[$pvalue_col]+0;
	my $this_pos = $fields[$pos_col]+0;
	if($this_p <= $pthres)
	{
		my $name = $fields[$snp_col];
		$_ = $name;
		if(/^rs/)
		{
			if(abs($this_pos - $last_pos) > $dist_thres)
			{
				$set_count++;
				$this_set_count = 0;
			}	
			$pvalues{$name} = $this_p;
			# print "\n Storing $name in set = $set_count in slot $this_set_count";
			$snp_sets[$set_count][$this_set_count] = $name;
			$this_set_count++;
			$last_pos = $this_pos;
			my $kk = abs($fields[$or_col]);
			if( ($kk > 0.01) && ($kk < 3) )
			{
				$or{$name} = $fields[$or_col];
				$ref_freq{$name} = $fields[$freq_col];
				# print "\n found $name with or = .$or{$name}. and beta = .$beta{$name}. ";
				my $pene_effect = $prev/($ref_freq{$name} + (1-$ref_freq{$name})/exp($or{$name}));
				my $pene_ref = ($prev -$pene_effect*$ref_freq{$name})/(1-$ref_freq{$name});
				my $alpha_effect = $thres - inv_cdf(1.0 - $pene_effect); 
				my $alpha_ref = $thres - inv_cdf(1.0 - $pene_ref); 
				$beta{$name} = $alpha_effect - $alpha_ref;
				# print "\n $name $this_pos $ref_freq{$name} $pene_effect $pene_ref $beta{$name}"; 
			}
		}
	}
}
close(FILE);
# print "\n About to start thinning effects \n\n";
foreach my $i (@snp_sets)
{
	my $best = @$i[0];
	foreach my $j (@$i) 
	{
		# print "\n \t comparing $j and $best \n";
		if($pvalues{$j} < $pvalues{$best})
		{
				$best = $j;
		}
	}
	# print "\n best pvalue is $best $pvalues{$best} \n";
	foreach my $j (@$i) 
	{
		if($j ne $best)
		{
			delete $beta{$j};
			delete $or{$j};
		}
	}
}
# print "\n Finished thinning \n";
open(FILE,"$ARGV[3]") || die "\n Can't read annotations file $ARGV[3]\n\n";
$_ = <FILE>;
chomp;
@fields = split('\t');
my $het_col = -1;
my $hom_col = -1;
my $missing_col = -1;
$snp_col = -1;
my %anc_col;
foreach my $pop (@super_pop)
{
	$anc_col{$pop} = -1;
}

for(my $i = 0; $i < @fields ; $i++)
{
	if($fields[$i] eq "heterozygotes")
	{
		$het_col = $i;
	}
	elsif($fields[$i] eq "homozygotes")
	{
		$hom_col = $i;
	}
	elsif($fields[$i] eq "missingGenos")
	{
		$missing_col = $i;
	}
	elsif($fields[$i] eq "gnomad.genomes.id")
	{
		$snp_col = $i;
	}
	else
	{
		foreach my $pop (@super_pop)
		{
			my $j = "gnomad.genomes.AF\_joint\_$pop";
			if($fields[$i] eq $j)
			{
				$anc_col{$pop} = $i;
			}
		}
	}
}

if($snp_col < 0) { die "\n Could not find a field for variant name \n"; };
if($het_col < 0) { die "\n Could not find a field for the heterozygotes \n"; };
if($hom_col < 0) { die "\n Could not find a field for the homozygotes \n"; };
if($missing_col < 0) { die "\n Could not find a field for the missing \n"; };
foreach my $pop (@super_pop)
{
	if($anc_col{$pop} < 0)
	{
		die "\n Could not find the allele frequency column for $pop \n";
	}
}

my $ave_Fst =0;
my %prs_data;
my %prs_naive;
my %or_data;
my @samples = (sort (keys %sam_anc));
my $N = @samples+0;
my $snp_count = 0;
my $Va = 0;
my $ref_or_mean = 0;
my $ref_or_var = 0;
my $ref_prs_mean = 0;
my $ref_prs_var = 0;
while(<FILE>)
{
	chomp;
	@fields = split('\t');
	my @this_mis;
	my @this_het; 
	my @this_hom;
	my $this;
	$this =  $fields[$snp_col];
	if(exists($beta{$this}))
	{
		my $this_or = $or{$this};
		my $this_beta = $beta{$this};
		# print "\n$this\t$this_or\t$this_beta";
		@this_mis = split('\;',$fields[$missing_col]);
		@this_het = split('\;',$fields[$het_col]);
		@this_hom = split('\;',$fields[$hom_col]);
		my $mean_q;
		my $var_q;
		my %sample_q;
		my %maf;
		my $this_mean_var;
		foreach my $s (@samples)
		{
			$sample_q{$s} = 0;
			$maf{$s} = 0;
			my $this_var_q =0 ;
			foreach my $pop (@super_pop)
			{
				$sample_q{$s} += $fields[$anc_col{$pop}]*$sam_anc{$s}{$pop};
				$this_var_q += $fields[$anc_col{$pop}]*$fields[$anc_col{$pop}]*$sam_anc{$s}{$pop};
			}
			my $this_q = $sample_q{$s};
			$this_var_q -= $this_q*$this_q;
			$this_mean_var += $this_var_q;
			$mean_q += $this_q;
			$var_q += $this_q*$this_q;
		}
		$mean_q /= $N;
		$var_q /= $N;
		$this_mean_var /= $N;
		$var_q -= $mean_q*$mean_q;
		if($mean_q > 0)
		{
			$snp_count++;
			$var_q +=  $this_mean_var;
			$ave_Fst += $var_q / ( $mean_q*(1-$mean_q));
		}
		foreach my $i (@this_mis)
		{
			$maf{$i} = -1;
		}
		foreach my $i (@this_het)
		{
			$maf{$i} = 1;
		}
		foreach my $i (@this_hom)
		{
			$maf{$i} = 2;
		}
		$Va += $this_beta*$this_beta*2*$mean_q*(1-$mean_q);
		$ref_or_mean += $this_or * 2.0 * $ref_freq{$this};
		$ref_prs_mean += $this_beta * 2.0 * $ref_freq{$this};
		$ref_prs_var += $this_beta*$this_beta*2*$ref_freq{$this}*(1-$ref_freq{$this});
		$ref_or_var += $this_or*$this_or*2*$ref_freq{$this}*(1-$ref_freq{$this});
		
		foreach my $s (@samples)
		{
			if($maf{$s} >= 0)
			{
				$prs_data{$s} += $this_beta*($maf{$s} - 2.0*$sample_q{$s});
				if($maf{$s} >0)
				{
					$prs_naive{$s} += $this_beta*($maf{$s});
					$or_data{$s} += $this_or*($maf{$s});
				}
				# print "$s $this_beta $this_or $maf{$s} $mean_q\n";
			}
		}
	}
}
close(FILE);
$ave_Fst /= $snp_count;
my %mean_or;
my %mean_prs_or;
my %mean_afr;
my %var_afr;
my %mean_prs_naive;
my %var_prs_naive;
my %mean_prs;
my %var_or;
my %var_prs_or;
my %var_prs;
my %cov_or;
my %cov_prs_naive;
my %cov_prs;
my %scount;
my %mean_log_or;
my %var_log_or;
my $or_prev = (1.0 - $prev) / $prev;
my $Ve = sqrt(1.0 - $Va);
foreach my $s (@samples)
{
	if(exists($sex{$s}))
	{
		$mean_afr{$affect{$s}} += $sam_anc{$s}{"afr"};
		$var_afr{$affect{$s}} += $sam_anc{$s}{"afr"}*$sam_anc{$s}{"afr"};

		$mean_log_or{$affect{$s}} += $or_data{$s};
		$var_log_or{$affect{$s}} += $or_data{$s}*$or_data{$s};

		my $kk = exp($or_data{$s});

		$mean_or{$affect{$s}} += $kk;
		$var_or{$affect{$s}} += $kk*$kk;
		$cov_or{$affect{$s}} += $sam_anc{$s}{"afr"}*$kk;

		my $this_prev = 1.0 - cdf($thres,$prs_data{$s},$Ve);
		my $real_or = $or_prev * $this_prev/ (1.0-$this_prev);

		$mean_prs_or{$affect{$s}} += $real_or;
		$var_prs_or{$affect{$s}} += $real_or*$real_or;

		$mean_prs{$affect{$s}} += $prs_data{$s};
		$mean_prs_naive{$affect{$s}} += $prs_naive{$s};
		$var_prs{$affect{$s}} += $prs_data{$s}*$prs_data{$s};
		$var_prs_naive{$affect{$s}} += $prs_naive{$s}*$prs_naive{$s};
		$cov_prs{$affect{$s}} += $sam_anc{$s}{"afr"}*$prs_data{$s};
		$cov_prs_naive{$affect{$s}} += $sam_anc{$s}{"afr"}*$prs_naive{$s};

		$scount{$affect{$s}}++;

		$mean_afr{"All"} += $sam_anc{$s}{"afr"};
		$var_afr{"All"} += $sam_anc{$s}{"afr"}*$sam_anc{$s}{"afr"};

		$mean_log_or{"All"} += $or_data{$s};
		$var_log_or{"All"} += $or_data{$s}*$or_data{$s};

		$mean_or{"All"} += $kk;
		$var_or{"All"} += $kk*$kk;
		$cov_or{"All"} += $sam_anc{$s}{"afr"}*$kk;
		$mean_prs_or{"All"} += $real_or;
		$var_prs_or{"All"} += $real_or*$real_or;


		$mean_prs{"All"} += $prs_data{$s};
		$mean_prs_naive{"All"} += $prs_naive{$s};
		$var_prs{"All"} += $prs_data{$s}*$prs_data{$s};
		$var_prs_naive{"All"} += $prs_naive{$s}*$prs_naive{$s};
		$cov_prs{"All"} += $sam_anc{$s}{"afr"}*$prs_data{$s};
		$cov_prs_naive{"All"} += $sam_anc{$s}{"afr"}*$prs_naive{$s};

		$scount{"All"}++;
	}
}

my $k = exp($ref_or_mean + $ref_or_var/2);
print "Total Snps\tFst\tVa\tTraining_log(OR)_Mean\tTraining_log(OR)_Var\tTraining_OR_mean\tTraining_PRS_Mean\tTraining_PRS_Var\n";
print "$snp_count\t$ave_Fst\t$Va\t$ref_or_mean\t$ref_or_var\t$k\t$ref_prs_mean\t$ref_prs_var\n\n";
print "Affectation\tCount\tMean_Naive_log(OR)\tVar_Naive_log(OR)\tMean_Naive(OR)\tVar_Naive(OR)\tMean_Or\tVar_OR\tMean_PRS\tVar_PRS\tMean_PRS_Naive\tVar_PRS_Naive\tMean_afr\tCorr[Naive_Or,Afr]\tCorr[PRS,Afr]\tCorr[PRS_naive,Afr]\n";
foreach my $i (sort (keys %scount))
{
	if($scount{$i} > 0)
	{
		$mean_or{$i} /= $scount{$i};
		$mean_log_or{$i} /= $scount{$i};
		$mean_prs_or{$i} /= $scount{$i};
		$mean_prs{$i} /= $scount{$i};
		$mean_prs_naive{$i} /= $scount{$i};
		$mean_afr{$i} /= $scount{$i};
		$var_or{$i} /= $scount{$i};
		$var_log_or{$i} /= $scount{$i};
		$var_prs_or{$i} /= $scount{$i};
		$var_prs{$i} /= $scount{$i};
		$var_prs_naive{$i} /= $scount{$i};
		$var_afr{$i} /= $scount{$i};
		$cov_or{$i} /= $scount{$i};
		$cov_prs{$i} /= $scount{$i};
		$cov_prs_naive{$i} /= $scount{$i};

		$var_or{$i} -= $mean_or{$i}*$mean_or{$i};
		$var_log_or{$i} -= $mean_log_or{$i}*$mean_log_or{$i};
		$var_prs_or{$i} -= $mean_prs_or{$i}*$mean_prs_or{$i};
		$var_prs{$i} -= $mean_prs{$i}*$mean_prs{$i};
		$var_prs_naive{$i} -= $mean_prs_naive{$i}*$mean_prs_naive{$i};
		$var_afr{$i} -= $mean_afr{$i}*$mean_afr{$i};
		$cov_or{$i} -= $mean_or{$i}*$mean_afr{$i};
		$cov_prs{$i} -= $mean_prs{$i}*$mean_afr{$i};
		$cov_prs_naive{$i} -= $mean_prs_naive{$i}*$mean_afr{$i};
		$cov_or{$i} /= sqrt($var_or{$i}*$var_afr{$i});
		$cov_prs{$i} /= sqrt($var_prs{$i}*$var_afr{$i});
		$cov_prs_naive{$i} /= sqrt($var_prs_naive{$i}*$var_afr{$i});
	
		print "$i\t$scount{$i}\t$mean_log_or{$i}\t$var_log_or{$i}\t$mean_or{$i}\t$var_or{$i}\t$mean_prs_or{$i}\t$var_prs_or{$i}\t$mean_prs{$i}\t$var_prs{$i}\t$mean_prs_naive{$i}\t$var_prs_naive{$i}\t$mean_afr{$i}\t$cov_or{$i}\t$cov_prs{$i}\t$cov_prs_naive{$i}\n";
	}
}
print "\nSample\tSex\tAffectation\tPercent_Afr\tNaive_Or\tPRS\tCorrected_Or\n";
foreach my $s (@samples)
{
	my $this_or = exp($or_data{$s});
	my $this_prev = 1.0 - cdf($thres,$prs_data{$s},$Ve);
	my $real_or = $or_prev * $this_prev/ (1.0-$this_prev);
	my $j = "afr";
	if(exists($sex{$s}))
	{
		print "$s\t$sex{$s}\t$affect{$s}\t$sam_anc{$s}{$j}\t$this_or\t$prs_data{$s}\t$real_or\n";
	}
}
