#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Cwd qw(abs_path);

my ($idlist, $repseq, $table, $outdir) = @ARGV;

GetOptions(
		'id:s'      =>  \$idlist ,
		'rep:s'      =>  \$repseq,
		'tab:s'    =>  \$table,
		'outdir:s' =>  \$outdir,
);

$outdir ||= './';
-d $outdir || mkdir $outdir;
$outdir = abs_path($outdir);

$idlist && -s $idlist  || die "perl $0 --id id.list --rep rep-seq.fna --tab feature-table.tsv --outdir ./\n";
$repseq && -s $repseq  || die "perl $0 --id id.list --rep rep-seq.fna --tab feature-table.tsv --outdir ./\n";
$table && -s $table  || die "perl $0 --id id.list --rep rep-seq.fna --tab feature-table.tsv --outdir ./\n";

my $outtab =  $outdir . '/rename-feature-table.tsv';
my $outrepseq =  $outdir . '/rename-rep-seq.fna';

my %old_new = id_hash($idlist);

open TAB, $table;
open OUT1, ">$outtab";
while(<TAB>){
  /^#/ && print OUT1 "$_";
  /^#/ && next;
  chomp;
  my @temp = split /\t/;
  my $old_id = shift @temp;
  $old_new{$old_id} || next;
  print OUT1 "$old_new{$old_id}\t", join("\t", @temp), "\n";
  #print OUT1 "$old_new{$old_id}\t", join("\t", @temp), "\n";
}
close TAB;
close OUT1;

open OUT2, ">$outrepseq";
open REP, $repseq;
$/ = ">";
<REP>;
while(<REP>){
  chomp;
  my ($old_id, $seq) = split /\s+/;
  #print $old_id, "\t", $seq, "\n";
  $old_new{$old_id} || next;
  print OUT2 ">$old_new{$old_id}\n$seq\n";
}
$/ = "\n";
close REP;
close OUT2;

sub id_hash{
  my ($f) = @_;
  open IN, $f;
  my %hash;
  while(<IN>){
    chomp;
    my @temp = split /\t/;
    $hash{$temp[0]} = $temp[1];
  }
  close IN;
  return %hash;
}
