#!usr/bin/perl -w
use strict;

my ($cohortlist, $penn_dir, $ipn_penn_intersect, $out) = @ARGV;
$penn_dir =~ s/\/$//g;
my @cohorts;
open COHORTLIST, "$cohortlist";
while(<COHORTLIST>){
	chomp; push @cohorts, $_;
}
close COHORTLIST;

open OUT, ">$out";
print OUT "Cohort\tSample\tLRR_SD\tBAF_SD\tGCWF\tCNV_num\tCNV_total_length\n";
my %info1;
for my $cohort (@cohorts){
	print "$cohort\n";
	my $logfile = "$penn_dir/penn/$cohort/out/sample.log";
	print "log file is $logfile\n";
	open LOG, "$logfile";
	while(<LOG>){
		if(/quality summary/){
			$_ =~ /.*\/(.*)\: LRR_mean=(.*) LRR_median=(.*) LRR_SD=(.*) BAF_mean=(.*) BAF_median=(.*) BAF_SD=(.*) BAF_DRIFT=(.*) WF=(.*) GCWF=(.*)\n/;
			my ($sample, $lrr_mean, $lrr_median, $lrr_sd, $baf_mean, $baf_median, $baf_sd, $baf_drift, $wf, $gcwf) = ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10);
			$info1{$cohort}{$sample} = join "\t", ($lrr_sd, $baf_sd, $gcwf);
		}
	}
	close LOG;
}

my %cnv_num;
my %total_length;
open RAWCNV, "$ipn_penn_intersect";
while(<RAWCNV>){
	next if $. == 1;
	chomp;
	my ($fid, $iid, $chr, $start, $end, $type) = (split)[0, 1, 2, 3, 4, 5];
	my ($cohort, $sample) = split /\*/, $fid;
	my $length=$end-$start;
	my @entry = split /\s+/, $_;
	$cnv_num{$cohort}{$sample} += 1;
	$total_length{$cohort}{$sample} += $length;
}
close RAWCNV;

for my $cohort (keys %cnv_num){
	for my $sample(keys %{$cnv_num{$cohort}}){
			print OUT "$cohort\t$sample\t$info1{$cohort}{$sample}\t$cnv_num{$cohort}{$sample}\t$total_length{$cohort}{$sample}\n";
	}
}
close OUT;
