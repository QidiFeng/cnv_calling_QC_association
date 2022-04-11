#!usr/bin/perl -w
use strict;
use lib './';
use Func qw( which_col );

my ($cohort, $finalreport, $phenofile, $penndir, $ipndir) = @ARGV;

$penndir =~ s/\/$//g;
$ipndir =~ s/\/$//g;

##making directories for inputs.
my ($penn_data, $penn_dataaux, $penn_out, $ipn_data, $ipn_dataaux, $ipn_out) = ($penndir.'/penn/'.$cohort.'/data', $penndir.'/penn/'.$cohort.'/data_aux', $penndir.'/penn/'.$cohort.'/out', $ipndir.'/ipn/'.$cohort.'/data', $ipndir.'/ipn/'.$cohort.'/data_aux', $ipndir.'/ipn/'.$cohort.'/out');
my @dirs= ($penn_data,$penn_dataaux,$ipn_data,$ipn_dataaux, $penn_out, $ipn_out);
for my $dir (@dirs){
	if(-d $dir){
		print "$dir already exists, good, processing...\n";
	}else{
		`mkdir -p $dir`;
		print "making $dir\n";
	}
}


##producing inputs.
open FINALREPORT, "$finalreport";
my $turn=0;
my $sample0="INIT"; my $printsnp=0; my %exist_samples;
my ($snpidcol, $samplecol, $allel1forcol, $allele2forcol, $chrcol, $poscol, $xcol, $ycol, $ballelecol, $logrcol);
			print "outputing signal files, snp files...\n";
while(<FINALREPORT>){
	$_ =~ s/^(\s+)//g; $_ =~ s/(\s+)$//g;
	##skip the Header
	if(/SNP Name/ and /Sample ID/){
		($snpidcol, $samplecol, $allel1forcol, $allele2forcol, $chrcol, $poscol, $xcol, $ycol, $ballelecol, $logrcol) = which_col($_);
		$turn=1; next;
	}
	##read the important part.
	if($turn==1){
		my @entry = split /\t/, $_;
		my ($snpid, $sample, $allele1for, $allele2for, $chr, $pos, $x, $y, $ballele, $logr) = @entry[$snpidcol, $samplecol, $allel1forcol, $allele2forcol, $chrcol, $poscol, $xcol, $ycol, $ballelecol, $logrcol];
		if($sample ne $sample0){
			if($sample0 eq "INIT"){
				my $pennsnp = "snplist.txt";
				my $ipnsnp = "probe.txt";
				$printsnp = 1;
				open PENSNP, ">$penn_dataaux/$pennsnp";
				print PENSNP "Name\tChr\tPosition\n";
				open IPNSNP, ">$ipn_dataaux/$ipnsnp";
				print IPNSNP "SNP Name\tChr\tPosition\n";
			}else{
				$printsnp = 0;
				close PENSNP; close IPNSNP;
				close PENNSIG; close IPNSIG;
			}
			$exist_samples{$sample}=1;
			my $pennsig = $sample;
			my $ipnsig = $sample.".txt";
			open PENNSIG, ">$penn_data/$pennsig";
			print PENNSIG "Name\t$sample.Log R Ratio\t$sample.B Allele Freq\n";
			open IPNSIG, ">$ipn_data/$ipnsig";
			print IPNSIG "[Header]\nSample ID\tSNP Name\tChr\tPosition\tAllele1 - Forward\tAllele2 - Forward\tX\tY\tB Allele Freq\tLog R Ratio\n";
		}
			print PENNSIG "$snpid\t$logr\t$ballele\n";
			print IPNSIG "$sample\t$snpid\t$chr\t$pos\t$allele1for\t$allele2for\t$x\t$y\t$ballele\t$logr\n";
			if($printsnp==1){
				print PENSNP "$snpid\t$chr\t$pos\n";
				print IPNSNP "$snpid\t$chr\t$pos\n";
			}
			$sample0 = $sample;
	}
}
close FINALREPORT;

##Producing sample info files.
open PHENO, "$phenofile";
open OUT_PEN_SAMPLE, ">$penn_dataaux/samplelist";
open OUT_IPN_SAMPLE, ">$ipn_dataaux/samplelist.txt";
open OUT_IPN_GENDER, ">$ipn_dataaux/gender.txt";
open OUT_IPN_BADSAMPLE, ">$ipn_dataaux/badsample.txt";
close OUT_IPN_BADSAMPLE;
print "outputing sample list...\n";
while(<PHENO>){
	chomp; next if $.==1;
	my ($fid, $gender) = (split /\t/, $_)[0, 3];
	$gender =~ tr/1/M/; $gender =~ tr/2/F/;
	my ($tmpcohort, $sample) = split /\*/, $fid;
	next unless $tmpcohort eq $cohort;
	if(exists $exist_samples{$sample}){
		print OUT_PEN_SAMPLE "$penn_data/$sample\n";
		print OUT_IPN_SAMPLE "$ipn_data/$sample.txt\n";
		print OUT_IPN_GENDER "$sample\t$gender\n";
	}
}
close PHENO;
close OUT_PEN_SAMPLE;
close OUT_IPN_SAMPLE;
close OUT_IPN_GENDER;

##split ipn samplelist if sample number is >300;
`split -l 200 $ipn_dataaux/samplelist.txt $ipn_dataaux/sublist`;

##produce ipn parameter file;
	my @batchs = `ls $ipn_dataaux/sublist*`;
	for my $i (@batchs){
		$i =~ s/^(\s+)//g; $i =~ s/(\s+)$//g;
		my $batch = (split /\//, $i)[-1];
	my $conf = "$ipn_dataaux/conf.$batch.txt";
	open CONF, ">$conf";
	print CONF "GENDERFILE=$ipn_dataaux/gender.txt\n";
	print CONF "BADNAME=$ipn_dataaux/badsample.txt\n";
	print CONF "DATAFILE=$ipn_dataaux/$batch\n";
	print CONF "EXPERIMENT=$cohort.$batch\n";
	print CONF "OUTDIR=$ipn_out\n";
	print CONF "PROBEFILE=$ipn_dataaux/probe.txt\n";
	print CONF "TEMPDIR=$ipn_out/\n";
	print CONF "DATADIR=$ipn_out/\n";
	print CONF "CALLDIR=$ipn_out/\n";
	print CONF "DOLOG=TRUE\nDOCLEANUP=TRUE\n";
	print CONF "NOQSUB=YES\n";
	}
