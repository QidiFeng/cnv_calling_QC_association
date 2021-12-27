#Notice: please use the same $ipn_dir and $penn_dir as before in the first part.
#$chromosome_length file, exon_gene_list file and centromere_telomere file can be downloaded from github together with codes.

ipn_dir=$1
penn_dir=$2
out_dir=$3
cohort_list=$4
chromosome_length=$5
phe_file=$6
exon_gene_list=$7
centromere_telomere=$8

mkdir $out_dir
echo "formatting ipattern raw results"
perl 02_1_ipn_format_rawcnv.pl $cohort_list $ipn_dir $out_dir/02_1_ipn_rawcnv.txt

echo "formatting penncnv raw results"
perl 02_2_penn_format_rawcnv.pl $cohort_list $penn_dir $out_dir/02_2_penn_rawcnv.txt

echo "Getting intersection between ipattern and penncnv"
perl 02_3_intersect_ipn_penn.pl $out_dir/02_1_ipn_rawcnv.txt $out_dir/02_2_penn_rawcnv.txt $out_dir/02_3_ipn_penn_intersection.txt

echo "Get information of LRRSD, BAFSD, GCFSD, GCWF, CNV proportion, CNV number for each sample"
perl 02_4_extract_lrrsd_bafsd_gcfsd_gcwf_cnvprop_cnvnum.pl $cohort_list $penn_dir $out_dir/02_3_ipn_penn_intersection.txt $out_dir/02_4_rawcnv_quality_summary.txt

echo "Get sample outlier due to LRRSD, BAFSD, GCFSD, GCWF, CNV proportion and CNV number"
Rscript 02_5_get_sample_outlier.r $out_dir/02_4_rawcnv_quality_summary.txt $cohort_list $out_dir/02_5_sample_outlier.txt $out_dir/02_5_sample_outlier_criteria.table.txt

echo "Annealing CNV"
perl 02_6_anneal_rawcnv.pl $out_dir/02_3_ipn_penn_intersection.txt $penn_dir $cohort_list $exon_gene_list $out_dir/02_6_ipn_penn_intersection.sorted.txt $out_dir/02_6_annealed.txt

echo "Identify Aneuploidy outlier"
perl 02_7_identify_aneuploidy.pl $out_dir/02_6_annealed.txt $chromosome_length $out_dir/02_7_aneuploidy_info.txt $out_dir/02_7_aneuploidy_outlier.txt

echo "Summarizing all sample outliers"
perl 02_9_summarize_allsampleoutliers.pl $out_dir/02_5_sample_outlier.txt  $out_dir/02_7_aneuploidy_outlier.txt $phe_file $out_dir/02_8_all_sample_outliers.txt

echo "Removing all sample outliers"
perl 02_10_remove_sample_outlier.pl $out_dir/02_8_all_sample_outliers.txt $phe_file $out_dir/02_6_annealed.txt $out_dir/02_9_sampleoutlier_removed.cnv $out_dir/02_9_sampleoutlier_removed.fam

echo "Filtering CNVs"
perl 02_11_CNV_filtering_byPLINK.pl $cohort_list $centromere_telomere $out_dir/02_9_sampleoutlier_removed $out_dir/02_11_plinked
