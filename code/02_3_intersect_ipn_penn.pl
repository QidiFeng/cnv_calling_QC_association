#!usr/bin/perl -w
use strict;
use lib './';
use Annot qw( intersect_func);

my ($ipn_rawcnv, $penn_rawcnv, $out) = @ARGV;
##for usage of intersect_func
open IPN, "$ipn_rawcnv";
my %ipn;
while(<IPN>){
	next if $. == 1;
	chomp;
	$_ =~ s/Gain/3/g;
	$_ =~ s/Loss/1/g;
	my ($cohort, $sample, $type, $chr, $start, $end) = split;
	$ipn{$cohort}{$sample}{$type}{$chr}{$start} = $end;
}
close IPN;

##main script starts here
open PEN, "$penn_rawcnv";
open OUT, ">$out";
print OUT "FID\tIID\tCHR\tBP1\tBP2\tTYPE\n";
while(<PEN>){
	next if $. == 1;
	chomp;
	my ($fid, $iid, $chr, $start, $end, $type) = (split)[0, 1, 2, 3, 4, 5];
	next unless (($type eq '1')|($type eq '3'));
	my ($cohort, $sample) = split /\*/, $fid;
	my ($bn1_bn2) = intersect_func(\%ipn, $cohort, $sample, $chr, $start, $end, $type);
	my %bn1_bn2 = %{$bn1_bn2};
	my @tmp_bn1_bn2 = keys %bn1_bn2;
	next if $#tmp_bn1_bn2 == 0;
	delete $bn1_bn2{"NA"};
	for my $bd1 (keys %bn1_bn2){
		my $bd2 = $bn1_bn2{$bd1};
		print OUT "$fid\t$iid\t$chr\t$bd1\t$bd2\t$type\n";
	}
}
close OUT;
close PEN;

