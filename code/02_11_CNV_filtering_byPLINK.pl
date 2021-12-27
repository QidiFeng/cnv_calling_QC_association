#!usr/bin/perl -w
use strict;

my ($cohortlist, $centromere_telomere, $cnv_outlierremoved, $cnv_plinked) = @ARGV;
my @cohorts;
open COHORTLIST, "$cohortlist";
while(<COHORTLIST>){
	chomp; push @cohorts, $_;
}
close COHORTLIST;

##creat map file for outlierremoved.cnv
`plink --noweb --cnv-list $cnv_outlierremoved.cnv --cnv-make-map --out $cnv_outlierremoved`;

##cohort specific maf 
for my $cohort(@cohorts){
	`grep $cohort $cnv_outlierremoved.fam > tmp.$cohort.fam`;
}
for my $cohort(@cohorts){
	my $sample =  `wc -l tmp.$cohort.fam`;
	my $number = (split /\s+/, $sample)[0];
	print "sample is $number\n";
	my $freq = int($number/100);
	`plink --noweb --cfile $cnv_outlierremoved --keep tmp.$cohort.fam --cnv-freq-exclude-above $freq --cnv-write --out tmp2.$cohort`;
}
###combine all cohorts for .cnv file and .fam file
open OUTCNV, ">tmp3.cnv";
open OUTFAM, ">tmp3.fam";
for my $i(0..$#cohorts){
	my $cohort=$cohorts[$i];
	open FILE, "tmp2.$cohort.cnv";
	while(<FILE>){
		if($. == 1){
			if($i==0){
				print OUTCNV "$_";next;
			}else{next;}
		}
		print OUTCNV "$_";
	}
	close FILE;
	
	open FILE, "tmp2.$cohort.fam";
	while(<FILE>){
		print OUTFAM "$_";
	}
	close FILE;
}
close OUTCNV;
close OUTFAM;

##creat map file for combined cohorts.
`plink --noweb --cnv-list tmp3.cnv --cnv-make-map --out tmp3`;

## global maf, exclude specific regions, cnv sites >10, cnvkb > 20
	my $totalsample =  `wc -l tmp3.fam`;
	my $totalnumber = (split /\s+/, $totalsample)[0];
	print "sample is $totalnumber\n";
	my $totalfreq = int($totalnumber/100);
`plink --noweb --cfile tmp3 --cnv-overlap 0.50 --cnv-exclude $centromere_telomere --cnv-kb 20 --cnv-sites 10 --cnv-freq-exclude-above $totalfreq --cnv-write --out $cnv_plinked`;
`plink --noweb --cnv-list $cnv_plinked.cnv --cnv-make-map --out $cnv_plinked`;

##remove tmp files
`rm tmp*`;
