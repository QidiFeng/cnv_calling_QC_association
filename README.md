# Prerequisite
1. iPattern and PennCNV installed.
2. Finalreport file. Finalreport format has a header which contains "SNP Name", "Sample ID", "Allele1 - Forward", "Allele2 - Forward", "Chr", "Position", "X", "Y", "B Allele Freq" and "Log R Ratio"; And they were seperated by tab. 
3. Phenofile. $phenofile should contain FID, IID, AFF, SEX, ancestry, CNV\_platform, C1, C2, C3, C4, C5 columns.Please use "cohort\*sampleid" as FID and IID. FID could be the same with IID.DO NOT use "\*" or "-" within cohort names or sampleid (please replace "\*" or "-" with "\_").
4. Cohortlist. One cohort name each line. Cohort name should be consistent with cohort name in $phenofile.
5. Please keep all codes in ./code directory and conduct commands under ./code directory. 


# Part1: PennCNV and iPattern input preparation######

`perl 01_01_signalDensity_samplelist.pl $cohort_name $finalreport $phenofile $penndir $ipndir`

notice: This command process one cohort each time. If you have several cohorts, please repeatedly use this command to produce inputs for all cohorts. Please use full directories for penncnv\_input\_dir and ipn\_input\_dir. This script will produce signalfiles for penncnv under $penndir/$cohort/data, and snplist.txt, samplelist under $penndir/$cohort/data\_aux; Also it will produce all files needed for ipn under $ipndir/$cohort/data and $ipndir/$cohort/data\_aux

example: `perl 01_01_signalDensity_samplelist.pl Ma_xajd finalreport.txt MA.phe /mnt/disks/sdb/1_sc_asia/pipeline/03_for_github/ /mnt/disks/sdb/1_sc_asia/pipeline/03_for_github/`

2. `cd $penndir/$cohort/data_aux/`

3. `perl $penn_package_dir/compile_pfb.pl -listfile samplelist -snpposfile snplist.txt -output pfbfile`

notice: This script will produce pfbfile for Penncnv. compile\_pfb.pl script is in PennCNV package.

example: `perl /home/fengqidi/software/PennCNV-1.0.5/compile_pfb.pl -listfile samplelist -snpposfile snplist.txt -output pfbfile`

4. `cd $penndir/$cohort/data_aux/`

5. `perl $penn_package_dir/cal_gc_snp.pl gc5Base.txt snplist.txt -output gcmodel`

notice: This script will produce gcmodel for Penncnv. cal\_gc\_snp.pl script is in PennCNV package.

example: `perl /home/fengqidi/software/PennCNV-1.0.5/cal_gc_snp.pl gc5Base.txt snplist.txt -output gcmodel`


# Part2: Running PennCNV and iPattern######

## running iPattern
1. `cd $ipn_package_dir/ipn_0.582/preprocess/ilmn/`

2. `bash ilmn.sh $ipn_dir/$cohort/data_aux/parameter.file`

example: `bash ilmn.sh /mnt/disks/sdb/1_sc_asia/pipeline/03_for_github/ipn/Ma_xajd/data_aux/conf.sublistaa.txt`

## running penncnv

1. `perl $penn_package_dir/detect_cnv.pl -test -hmm $penn_package_dir/lib/hhall.hmm -pfb $penndir/$cohort/data_aux/pfbfile -gcmodel $penndir/$cohort/data_aux/gcmodel -listfile $penndir/$cohort/data_aux/samplelist -confidence -log $penndir/$cohort/out/sample.log -out $penndir/$cohort/out/sample.rawcnv`

example: `perl /home/fengqidi/software/PennCNV-1.0.5/detect_cnv.pl -test -hmm /home/fengqidi/software/PennCNV-1.0.5/lib/hhall.hmm -pfb pfbfile -gcmodel gcmodel -listfile samplelist -confidence -log ../out/sample.log -out ../out/sample.rawcnv`


# Part3: Intersection of Raw CNV (PennCNV and iPattern) and QC######

1. `bash 02_0_CNVcombine_QC.sh $ipn_dir $penn_dir $out_dir $cohort_list $chromosome_length $phenofile $exon_gene_list $centromere_telomere`

notice: 02\_0\_CNVcombine\_QC.sh is a pipeline bash script, which carry out all intersection and QC operations.
Please use the same $ipn\_dir and $penn\_dir as before in the first part. $chromosome\_length file, $exon\_gene\_list file and $centromere\_telomere file can be downloaded in ./example directory.

example: `bash 02_0_CNVcombine_QC.sh /mnt/disks/sdb/1_sc_asia/pipeline/03_for_github/ /mnt/disks/sdb/1_sc_asia/pipeline/03_for_github/ ./QC_ed`


# Part4: CNVburden calculation######
`Rscript 03_1_cnvburden.r $working_dir $out_dir/02_11_plinked $gene_file $output $phenofile`

