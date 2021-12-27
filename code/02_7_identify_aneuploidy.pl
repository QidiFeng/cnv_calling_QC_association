#!usr/bin/perl
use strict;

my ($annealed_in, $hg19_chr_length, $aneuploidy_info_out, $aneuploidy_outlier) = @ARGV;

my %total_length;
open ANNEALED, "$annealed_in";
while(<ANNEALED>){
	next if $. == 1;
	my ($fid, $iid, $chr, $bp1, $bp2) = (split /\s+/, $_)[0, 1, 2, 3, 4];
	$total_length{$fid}{$chr} += ($bp2-$bp1);
}
close ANNEALED;

###save chromosome length
open INFO, "$hg19_chr_length";
my %chr_length;
while(<INFO>){
	my ($chr, $length) = (split /\s+/,$_)[0,1];
	$chr_length{$chr}=$length;
}
close INFO;

open OUT, ">$aneuploidy_info_out";
print OUT "Cohort\tSample\tChr\tcnv_prop\n";
open OUT1, ">$aneuploidy_outlier";
print OUT1 "Cohort\tSample\tChr\tcnv_prop\tOUTLIER\n";
for my $fid (keys %total_length){
	my ($cohort, $sample) = split /\*/, $fid;
		for my $chr(sort {$a<=>$b} keys %{$total_length{$fid}}){
			my $prop = $total_length{$fid}{$chr}/$chr_length{$chr};
			print OUT "$cohort\t$sample\t$chr\t$prop\n";
			print OUT1 "$cohort\t$sample\t$chr\t$prop\tAneuploidy\n" if $prop>0.1;
		}
}
close OUT;
close OUT1;
