#remove LRR_SD, GCWF, BAF_SD, cnv_num, cnv_prop and aneuploidy outlier.
#remove what phe do not have, which exclude the GenotypeQC filtered samples.

#!usr/bin/perl -w
use strict;

my ($all_sample_outlier, $phe_file, $annealed_cnv, $cnv_out, $fam_out) = @ARGV;

my %outlier;
open OUTLIER, "$all_sample_outlier";
while(<OUTLIER>){
	chomp; my ($cohort, $sample) = split;
	$outlier{$cohort}{$sample}=1;
}
close OUTLIER;

open CNV, "$annealed_cnv";
my %kept_samples;
while(<CNV>){
	if($. == 1){
		next;
	}
	my @entry = split /\s+/, $_;
	my ($cohort, $sample)=split /\*/, $entry[0];
	next if exists $outlier{$cohort}{$sample};
	$kept_samples{$cohort}{$sample}=1;
}
close CNV;

my %phe_samples;
open CNV_PHE, "$phe_file";
open FAM, ">$fam_out";
while(<CNV_PHE>){
	next if $. == 1;
	my @entry = split /\s+/, $_;
	my ($fid, $iid, $aff, $sex) = @entry[0..3];
	my ($cohort, $sample)=split /\*/, $entry[0];
	next if exists $outlier{$cohort}{$sample};
	next if not exists $kept_samples{$cohort}{$sample};
	$phe_samples{$cohort}{$sample}=1;
	print FAM "$fid $iid 0 0 $sex $aff\n";
}
close CNV_PHE;

open CNV, "$annealed_cnv";
open CNV_OUT, ">$cnv_out";
while(<CNV>){
	if($. == 1){
		print CNV_OUT "$_";
		next;
	}
	my @entry = split /\s+/, $_;
	my ($cohort, $sample) = split /\*/, $entry[0];
	print CNV_OUT "$_" if exists $phe_samples{$cohort}{$sample};
}
close CNV;
close CNV_OUT;
