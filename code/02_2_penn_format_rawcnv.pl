#!usr/bin/perl -w
use strict;

my ($cohortlist, $penn_dir, $penn_out) = @ARGV;

$penn_dir =~ s/\/$//g;
my @cohorts;
open COHORTLIST, "$cohortlist";
while(<COHORTLIST>){
	chomp; push @cohorts, $_;
}
close COHORTLIST;

open PLINK_CNV, ">$penn_out";
print PLINK_CNV "FID\tIID\tCHR\tBP1\tBP2\tTYPE\tSCORE\tSITES\n";
for my $cohort (@cohorts){
	print "processing $cohort\n";
	my $cnvfile = "$penn_dir/penn/$cohort/out/sample.rawcnv";
	open CNV, "$cnvfile";
	while(<CNV>){
		chomp;
		my @entry = split /\s+/, $_;
		$entry[0] =~ /chr(.*):(.*)-(.*)/;
		my ($chr, $bp1, $bp2) = ($1, $2, $3);
		my $sites = (split /\=/, $entry[1])[-1];
		my $type = (split /\=/, $entry[3])[-1];
		my $score = (split /\=/, $entry[-1])[-1];
		my $ind = (split /\//, $entry[4])[-1];
		my $ind2 = $cohort."*".$ind;
		print PLINK_CNV "$ind2\t$ind2\t$chr\t$bp1\t$bp2\t$type\t$score\t$sites\n";
	}
	close CNV;
}
close PLINK_CNV;

