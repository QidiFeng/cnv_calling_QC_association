#!usr/bin/perl -w
use strict;
my ($cohortlist, $ipn_dir, $ipn_out) = @ARGV;
$ipn_dir =~ s/\/$//g;
my @cohorts;
open COHORTLIST, "$cohortlist";
while(<COHORTLIST>){
	chomp; push @cohorts, $_;
}
close COHORTLIST;

open OUT, ">$ipn_out";
print OUT "Cohort\tSample\tCNV_state\tChr\tStart\tEnd\n";
for my $cohort (@cohorts){
	print "$cohort\n";
	my @rawcnvs = `ls $ipn_dir/ipn/$cohort/out/*all_calls.txt`;
	for my $rawcnv(@rawcnvs){
		print "$cohort, $rawcnv\n";
	open RAWCNV, "$rawcnv";
	while(<RAWCNV>){
		next if /^#/;
		my @entry = split /\s+/, $_;
		my ($state, $chr, $start, $end, $sample) = @entry[0, 1, 2, 3, 11];
		next if (($state ne "Gain")and($state ne "Loss"));
		$sample = (split /X/, $sample)[-1];
		$sample =~ s/\./\-/g;
		print OUT "$cohort\t$sample\t$state\t$chr\t$start\t$end\n";
	}
	}
}
close OUT;
