#!/usr/bin/perl -w
use strict;

@ARGV == 1 || die "Usage: perl $0 genome.coverage.bed > bed.coverage.stat.xls\n";

my $bed_stat = shift;

#scaffold1	0	1	1155	0.000865801
#seqID	depth	number_of_bases_that_have_the_depth	seq_length	fraction_of_bases
my %seq_cov;
open BED, $bed_stat || die $!;
while(<BED>){
	chomp;
	my ($id, $depth, $num, $length, $perc_cov) = split /\t/;
	$depth == 0 && next;
	$seq_cov{$id}->{'mean_cov'} += $depth * $perc_cov;
	$seq_cov{$id}->{"len"} = $length;
	$seq_cov{$id}->{"cov_len"} += $num;
}
close BED;
print "ID\tLength\tCov_Length\tMeanCoverage\n";
for my $id (sort { $seq_cov{$a}->{'len'} <=> $seq_cov{$b}->{'len'} } keys %seq_cov){
	print "$id\t$seq_cov{$id}->{len}\t$seq_cov{$id}->{cov_len}\t$seq_cov{$id}->{mean_cov}\n";
}

