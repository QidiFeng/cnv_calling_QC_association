#!usr/bin/perl -w
use strict;
use lib './';
use Annot qw(numgene_exon_func numloci_func anneal_func offer_snp offer_exon);

my ($ipn_penn_intersect, $penn_dir, $cohortlist, $exon_gene_list, $ipn_penn_intersect_sorted, $annealed_out) = @ARGV;

##read rawcnv and save as cohort-sample-type-chr-start-end
open RAW, "$ipn_penn_intersect";
my %raw;
while(<RAW>){
	next if $. == 1;
	chomp;
	my ($fid, $iid, $chr, $bp1, $bp2, $type) = split /\t/, $_;
	my ($cohort, $sample)=split /\*/, $fid;
	$raw{$cohort}{$sample}{$type}{$chr}{$bp1}=$bp2;
}
close RAW;
print "saving raw hash\n";

open RAW_TMP, ">$ipn_penn_intersect_sorted";
	for my $cohort(sort keys %raw){
		for my $sample(sort keys %{$raw{$cohort}}){
			for my $type(sort keys %{$raw{$cohort}{$sample}}){
				for my $chr(sort keys %{$raw{$cohort}{$sample}{$type}}){
					for my $pos (sort {$a<=>$b} keys %{$raw{$cohort}{$sample}{$type}{$chr}}){
						print RAW_TMP "$cohort\t$sample\t$type\t$chr\t$pos\t$raw{$cohort}{$sample}{$type}{$chr}{$pos}\n";
				}
			}
		}
	}
}
close RAW_TMP;

my $annealed = anneal_func(\%raw);
my %annealed=%{$annealed};

##get snp annot and gene annot
my $snp = offer_snp($penn_dir, $cohortlist);
my %snp = %{$snp};
my ($exon_annot_end, $exon_annot_gene, $start_all, $end_all) = offer_exon($exon_gene_list);
my %exon_annot_end = %{$exon_annot_end};
my %exon_annot_gene = %{$exon_annot_gene};
my %start_all = %{$start_all};
my %end_all = %{$end_all};

##output annealed and do gene and snp annotation.
open OUT, ">$annealed_out";
print OUT "FID\tIID\tCHR\tBP1\tBP2\tTYPE\tSCORE\tSITES\n";
for my $cohort (sort keys %annealed){
	print "$cohort\n";
	for my $sample (sort keys %{$annealed{$cohort}}){
		my $fid = $cohort."*".$sample;
		for my $type (sort keys %{$annealed{$cohort}{$sample}}){
			for my $chr (sort {$a<=>$b} keys %{$annealed{$cohort}{$sample}{$type}}){
				for my $start (sort {$a<=>$b} keys %{$annealed{$cohort}{$sample}{$type}{$chr}}){
					my $end = $annealed{$cohort}{$sample}{$type}{$chr}{$start};
					my $genenum = numgene_exon_func(\%start_all, \%exon_annot_end, \%exon_annot_gene, $chr, $start, $end);
					my $locinum = numloci_func(\%snp, $cohort, $chr, $start, $end);
					my $length=$end - $start;
					my $exp_site = $length/20000;
					if($exp_site>$locinum){
						next;
					}
					else{
						print OUT "$fid\t$fid\t$chr\t$start\t$end\t$type\t$genenum\t$locinum\n";
					}
				}
			}
		}
	}
}
close OUT;
