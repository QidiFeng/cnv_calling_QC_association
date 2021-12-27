##combine LRR_SD, BAF_SD, GCWF, CNVnum, CNVprop, Aneuploidy outlier. as well as unknown labels in AFFT in phe.

#!usr/bin/perl -w
use strict;

my ($sample_outlier, $aneuploidy, $phe, $all_outlier) = @ARGV;
my %outlier;
open OUTLIER1, "$sample_outlier";
while(<OUTLIER1>){
	next if $. == 1;
	my ($cohort, $sample) = (split /\s+/, $_)[0, 1];
	$outlier{$cohort}{$sample}=1;
}
close OUTLIER1;

open OUTLIER2, "$aneuploidy";
while(<OUTLIER2>){
	next if $. == 1;
	my ($cohort, $sample) = (split /\s+/, $_)[0, 1];
	$outlier{$cohort}{$sample}=1;
}
##remove samples with AFF=-9 
open OUTLIER3, "$phe";
while(<OUTLIER3>){
	next if $. == 1;
	my ($fid, $aff) = (split)[0, 2];
	if(($aff ne "1")and($aff ne "2")){
		my ($cohort, $sample) = split /\*/, $fid;
		$outlier{$cohort}{$sample}=1;
	}
}
close OUTLIER3;

open OUT, ">$all_outlier";
for my $cohort (sort keys %outlier){
	for my $sample (sort keys %{$outlier{$cohort}}){
			print OUT "$cohort\t$sample\n";
	}
}
close OUT;
