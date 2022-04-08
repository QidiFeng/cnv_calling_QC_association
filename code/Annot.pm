package Annot;
use strict;
use warnings;
use Exporter;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(numgene_func numgene_exon_func numloci_func offer_snp offer_exon offer_gene merge_func anneal_func intersect_func);

##anneal segement of gam < 20%
sub anneal_func {
	my ($raw) =@_;
	my %raw = %{$raw};
	my %annealed;
	for my $cohort(keys %raw){
		for my $sample(keys %{$raw{$cohort}}){
			for my $type(keys %{$raw{$cohort}{$sample}}){
				for my $chr(keys %{$raw{$cohort}{$sample}{$type}}){
					my @pos = sort {$a<=>$b} keys %{$raw{$cohort}{$sample}{$type}{$chr}};
					my ($start, $end) = ($pos[0],$raw{$cohort}{$sample}{$type}{$chr}{$pos[0]});
					for my $i (@pos[1..$#pos]){
						my $j = $raw{$cohort}{$sample}{$type}{$chr}{$i};
						my $total_length = $j-$start;
						my $gap = $i-$end;
						my $gap_prop = $gap/$total_length;
						if($gap_prop<0.3){
							$start=$start;
							$end=$j;
						}
						else{
							$annealed{$cohort}{$sample}{$type}{$chr}{$start}=$end;
							$start=$i;
							$end=$j;
						}
					}
					$annealed{$cohort}{$sample}{$type}{$chr}{$start}=$end;
				}
			}
		}
	}
	return(\%annealed);
}



##number of genes##
sub numgene_func {
	my ($start_all, $gene_annot, $chr, $pos1, $pos2) = @_;
	my %start_all = %{$start_all};
	my %gene_annot = %{$gene_annot};
	my $count=0;
	for my $start(@{$start_all{$chr}}){
		if($start>$pos2){
			last;
		}
		my $end = $gene_annot{$chr}{$start};
		if((($start>=$pos1)and($start<=$pos2)) | (($end>=$pos1)and($end<=$pos2)) | (($start<=$pos1)and($end>=$pos2))){
			$count++;
			}
	}
	return($count);
}

##number of genes by only accounting for overlapping exons##
sub numgene_exon_func {
	my ($start_all, $exon_annot_end, $exon_annot_gene, $chr, $pos1, $pos2) = @_;
	my %start_all = %{$start_all};
	my %exon_annot_end = %{$exon_annot_end};
	my %exon_annot_gene = %{$exon_annot_gene};
	my %count;
	for my $start(@{$start_all{$chr}}){
		if($start>$pos2){
			last;
		}
		my $end = $exon_annot_end{$chr}{$start};
		if((($start>=$pos1)and($start<=$pos2)) | (($end>=$pos1)and($end<=$pos2)) | (($start<=$pos1)and($end>=$pos2))){
			my $gene = $exon_annot_gene{$chr}{$start}{$end};
			$count{$gene} = 1;
			}
	}
	my @counts = keys %count;
	my $count_num = $#counts+1;
	return($count_num);
}


##number of loci##
sub numloci_func {
	my ($snp, $cohort, $chr, $start, $end)=@_;
	my %snp = %{$snp};
	my $snpnum=0;
	for my $pos (sort {$a<=>$b} keys %{$snp{$cohort}{$chr}}){
		if(($start<=$pos) and ($pos<=$end)){
			$snpnum++;		
		}
	}
	return($snpnum);
}

##for useage of numloci_func
sub offer_snp {
	my ($penn_dir, $cohortlist) = @_;
	$penn_dir =~ s/\/$//;
	my @cohorts;
	open COHORTLIST, "$cohortlist";
	while(<COHORTLIST>){
		chomp; push @cohorts, $_;
	}
	close COHORTLIST;
	my %snp;
	for my $cohort(@cohorts){
	open SNPLIST, "$penn_dir/penn/$cohort/data_aux/snplist.txt";
	while(<SNPLIST>){
		next if $. == 1;
		chomp;
		my ($rsid, $chr, $pos) = split /\s+/, $_;
		$snp{$cohort}{$chr}{$pos}=$rsid;
	}
	close SNPLIST;
	}
	return(\%snp);
}

##for useage of numgene_func
sub offer_gene {
	my %gene_annot;
	open GENE, "/mnt/disks/sdb/1_sc_asia/240_annotation_file/hg19_refGene_plink.txt";
	while(<GENE>){
		chomp;
		my ($chr, $start, $end) = (split /\s+/, $_)[0, 1, 2];
		$gene_annot{$chr}{$start}=$end;
	}
	close GENE;
	my (%start_all, %end_all);
	for my $chr (sort keys %gene_annot){
		for my $start (sort {$a<=>$b} keys %{$gene_annot{$chr}}){
			push @{$start_all{$chr}}, $start;
			push @{$end_all{$chr}}, $gene_annot{$chr}{$start};
		}
	}
	return(\%gene_annot,\%start_all,\%end_all);
}
##for useage of numgene_exon_func
sub offer_exon {
	my ($exon_gene_list) = @_;
	my %exon_annot_end;
	my %exon_annot_gene;
	open GENE, "$exon_gene_list";
	while(<GENE>){
		next if $. == 1;
		next unless /cmpl/;
		chomp;
		my ($chr, $start, $end, $genename) = (split /\s+/, $_)[2, 6, 7, 12];
		$chr =~ /chr(.*)/;
		$chr = $1;
		$exon_annot_end{$chr}{$start}=$end;
		$exon_annot_gene{$chr}{$start}{$end}=$genename;
	}
	close GENE;

	my (%start_all, %end_all);
	for my $chr (sort keys %exon_annot_end){
		for my $start (sort {$a<=>$b} keys %{$exon_annot_end{$chr}}){
			push @{$start_all{$chr}}, $start;
			push @{$end_all{$chr}}, $exon_annot_end{$chr}{$start};
		}
	}
	return(\%exon_annot_end,\%exon_annot_gene,\%start_all,\%end_all);
}
##merg ipn and penncnv##
sub merge_func {
	my ($ipn, $cohort, $sample, $chr, $start, $end, $type) = @_;
	my %ipn = %{$ipn};
	my ($bn1, $bn2) = ("NA", "NA");
	my $intersect=0;
	my @points = ($start, $end);
	for my $i (sort {$a<=>$b} keys %{$ipn{$cohort}{$sample}{$type}{$chr}}){
		my $j = $ipn{$cohort}{$sample}{$type}{$chr}{$i};
		if((($i>=$start)and($i<=$end)) | (($j>=$start)and($j<=$end)) | (($i<=$start)and($j>=$end))){
			$intersect=1;
			push @points, ($i, $j);
		}
	}
	if($intersect==1){
	@points = sort {$a<=>$b} @points;
	($bn1, $bn2) = @points[0,-1];
	}
	return($bn1, $bn2);
}
##intersect ipn and penncnv##
sub intersect_func {
	my ($ipn, $cohort, $sample, $chr, $start, $end, $type) = @_;
	my %ipn = %{$ipn};
	my $intersect=0;
	my %bn1_bn2;
	$bn1_bn2{"NA"}="NA";
	for my $i (sort {$a<=>$b} keys %{$ipn{$cohort}{$sample}{$type}{$chr}}){
		my $j = $ipn{$cohort}{$sample}{$type}{$chr}{$i};
		if((($i>=$start)and($i<=$end)) | (($j>=$start)and($j<=$end)) | (($i<=$start)and($j>=$end))){
			my @points;
			push @points, ($i, $j, $start, $end);
			@points = sort {$a<=>$b} @points;
			my ($bn1, $bn2) = @points[1,2];
			$bn1_bn2{$bn1}=$bn2;
		}
	}
	return(\%bn1_bn2);
}
