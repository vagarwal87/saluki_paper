#!/usr/bin/perl

use allfxns;

$exprMat = shift;
$spec = ($exprMat =~ /mouse/)? 1 : 0; #species is mouse or human?

sub readFasta{
	local $fasta = shift;
	local %fasta = ();
	if ($spec){
		open DNA, "zcat $fasta |" || die "Could not open fasta file for $fasta\n";
	}
	else{
		open DNA, "<$fasta" || die "Could not open fasta file for $fasta\n";
	}
	while ($line = <DNA>){ chomp $line;
		if ($line =~ /^>\s?(\w+\d)\.?\d*/){ $header = $1; }
		else { $fasta{$header} .= $line; }
	}
	close DNA;
	return \%fasta;
}

if ($spec){
	open IN, "zcat Mus_musculus.GRCm3.90.chosenTranscript.gtf.gz | ";
}
else{
	open IN, "<Homo_sapiens.GRCh38.83.chosenTranscript.gtf";
}
while(<IN>){
	($region, $start, $stop, $last) = (split /\t/)[2,3,4,-1];
	($id) = ($last =~ /(ENS\w*GR?[\d|\.]+)/);
	$lengths{$id}{$region} += ($stop-$start);
	$cdsexoncount{$id}++ if $region eq 'CDS';
}
close IN;

if ($spec){
	%utr5p = %{ readFasta("mm10_ensembl90_5utrs.fa.gz") };
	%orfs = %{ readFasta("mm10_ensembl90_orfs.fa.gz") };
	%utr3p = %{ readFasta("mm10_ensembl90_3utrs.fa.gz") };
}
else{
	%utr5p = %{ readFasta("Homo_sapiens.GRCh38.83.chosenTranscript.5pUTRs.fa") };
	%orfs = %{ readFasta("Homo_sapiens.GRCh38.83.chosenTranscript.ORFs.fa") };
	%utr3p = %{ readFasta("Homo_sapiens.GRCh38.83.chosenTranscript.3pUTRs.fa") };
}

print join("\t", "ENSID", "HALFLIFE", "UTR5LEN", "CDSLEN", "INTRONLEN", "UTR3LEN", "UTR5GC", "CDSGC", "UTR3GC", "ORFEXONDENSITY", "5UTR", "ORF", "3UTR"), "\n";
open IN, "<$exprMat";
while(<IN>){ chomp;
	@a=split /\t/;
	$id = $a[0];
	if ($lengths{$id}{"five_prime_utr"} >= 10 && $lengths{$id}{"three_prime_utr"} >= 10 && $lengths{$id}{"CDS"} >= 10){
		print join("\t", $id, $a[2], int($lengths{$id}{"five_prime_utr"}), int($lengths{$id}{"CDS"}),
			int($lengths{$id}{"transcript"})-(int($lengths{$id}{"three_prime_utr"})+int($lengths{$id}{"CDS"})+int($lengths{$id}{"five_prime_utr"})),
			int($lengths{$id}{"three_prime_utr"}), gcContent($utr5p{$id}), gcContent($orfs{$id}), gcContent($utr3p{$id}),
			sprintf("%.2f", $cdsexoncount{$id}*1000/$lengths{$id}{"CDS"}), $utr5p{$id}, $orfs{$id}, $utr3p{$id}), "\n";
	}
	else{ $count++; }
}

print STDERR "$count IDs with region missing or too short\n";
