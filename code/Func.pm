package Func;
use strict;
use warnings;
use Exporter;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(which_col para_file);

##identify the column number for each factor
sub which_col {
	my ($header) = @_;
	my ($snpidcol, $samplecol, $allele1forcol, $allele2forcol, $chrcol, $poscol, $xcol, $ycol, $ballelecol, $logrcol);
	my @headers = split /\t/, $header;
	for my $i (0..$#headers){
		if($headers[$i] eq 'SNP Name'){
			$snpidcol=$i;
		}elsif($headers[$i] eq 'Sample ID'){
			$samplecol=$i;
		}elsif($headers[$i] eq 'Allele1 - Forward'){
			$allele1forcol = $i;
		}elsif($headers[$i] eq 'Allele2 - Forward'){
			$allele2forcol = $i;
		}elsif($headers[$i] eq 'Chr'){
			$chrcol=$i;
		}elsif($headers[$i] eq 'Position'){
			$poscol=$i;
		}elsif($headers[$i] eq 'X'){
			$xcol=$i;
		}elsif($headers[$i] eq 'Y'){
			$ycol=$i;
		}elsif($headers[$i] eq 'B Allele Freq'){
			$ballelecol=$i;
		}elsif($headers[$i] eq 'Log R Ratio'){
			$logrcol=$i;
		}
	}
	return($snpidcol, $samplecol, $allele1forcol, $allele2forcol, $chrcol, $poscol, $xcol, $ycol, $ballelecol, $logrcol);
}

