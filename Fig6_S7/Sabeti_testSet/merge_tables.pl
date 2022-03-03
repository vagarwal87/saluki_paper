#!/usr/bin/perl

use allfxns;

open IN, "<ENCFF770UJN.fasta";
while(<IN>){ chomp;
  $id = substr($_,1);
  $seq = <IN>;
  chomp $seq;
  @a = split /\//, $id;
  foreach $id (@a){
    if($id =~ /ref/){
      $var = "ref";
      $id =~ s/_ref//;
    } else {
      $var = "alt";
      $id =~ s/_alt//;
    }
    $gc{$id}{$var}=gcContent($seq);
  }
}
close IN;

open IN, "<ENCFF090JTW.fasta";
while(<IN>){ chomp;
  $id = substr($_,1);
  $seq = <IN>;
  chomp $seq;
  @a = split /\//, $id;
  foreach $id (@a){
    if($id =~ /ref/){
      $var = "ref";
      $id =~ s/_ref//;
    } else {
      $var = "alt";
      $id =~ s/_alt//;
    }
    $gc{$id}{$var}=gcContent($seq);
  }
}
close IN;

open IN, "<fastUTR_MPRA_predictions_ENCFF770UJN.txt";
<IN>;
while(<IN>){ chomp;
  ($id, $val) = split /\t/;
  @a = split /\//, $id;
  foreach $id (@a){
    if($id =~ /ref/){
      $var = "ref";
      $id =~ s/_ref//;
    } else {
      $var = "alt";
      $id =~ s/_alt//;
    }
    $vals{$id}{$var}=$val;
  }
}
close IN;

open IN, "<fastUTR_MPRA_predictions_ENCFF090JTW.txt";
<IN>;
while(<IN>){ chomp;
  ($id, $val) = split /\t/;
  @a = split /\//, $id;
  foreach $id (@a){
    if($id =~ /ref/){
      $var = "ref";
      $id =~ s/_ref//;
    } else {
      $var = "alt";
      $id =~ s/_alt//;
    }
    $vals{$id}{$var}=$val;
  }
}
close IN;

open IN, "<fastUTR_mpra.txt";
$header = <IN>;
print "Ref_pred\tAlt_pred\tRef_GC\tAlt_GC\t$header";
while(<IN>){
    $id = (split /\t/)[0];
    if (!defined($vals{$id}{"ref"}) || !defined($vals{$id}{"alt"})){
      print STDERR "$id skipped\n";
    } else{
      print join("\t", $vals{$id}{"ref"}, $vals{$id}{"alt"}, $gc{$id}{"ref"}, $gc{$id}{"alt"}, $_);
    }
}
close IN;
