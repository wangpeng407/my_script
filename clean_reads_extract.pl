#!/usr/bin/perl -w
#use strict;
use warnings;
use autodie qw(open close);
use FileCache;
use Cwd qw(abs_path);
use Getopt::Long;
use PerlIO::gzip;
#use FileCache maxopen => 10000;

my $samtools = "~/samtools/samtools-1.3/samtools";

#geneID-readsID     ReadsMapping.bowtie.bam
#mlgID-geneID       mgs.cluster
#readsID-sequences  fq1.gz & fq2.gz
my ($reads1, $reads2, $bamfile, $mgscluster, $outdir);
GetOptions(
		'1:s'      =>  \$reads1,
		'2:s'      =>  \$reads2,
		'bam:s'    =>  \$bamfile,
		'mgs:s'    =>  \$mgscluster,
		'outdir:s' =>  \$outdir
);

($reads1 && -s $reads1 && $reads2 && -s $reads2 && $mgscluster && -s $mgscluster) ||  &usage();

$outdir ||= './';
$outdir = abs_path($outdir);

$| = 1;

print STDERR datestr(), ": Program starting .\n\n";

my %ReadGene = readsID2geneID($bamfile);
print STDERR datestr(), ": Reading bamfile finished!\n\n";

my ($GeneMGS_ref, $MGSs) = geneID2mgsID($mgscluster);
print STDERR datestr(), ": Building Gene_ID and MGS_ID .\n\n";

my %GeneMGS = %{$GeneMGS_ref};
my @MGS = @{$MGSs};

$reads1 =~ /gz$/ ? (open FQ1, "<:gzip", $reads1 || die "Can not open $reads1.\n") : (open FQ1,  $reads1 || die "Can not open $reads1..\n") ;
$reads2 =~ /gz$/ ? (open FQ2, "<:gzip", $reads2 || die "Can not open $reads2.\n") : (open FQ2,  $reads2 || die  "Can not open $reads2..\n");
my %hash;
while(my $head1 = <FQ1>){
	my $seq1 = <FQ1>;  my $stand1 = <FQ1>; my $que1 = <FQ1>;
	my $head2 = <FQ2>; my $seq2 = <FQ2>; my $stand2 = <FQ2>; my $que2 = <FQ2>;
	my $readsid=(split /\s+/,$head1)[0];
	$readsid =~ s/^\@//;
	$ReadGene{$readsid} || next;
	my $cagid = $GeneMGS{$ReadGene{$readsid}};
#	print "readsid:$readsid\tgeneid:$ReadGene{$readsid}\tcagid:$cagid\n";
	my $handle1 = $cagid . 'OUT1';
	my $handle2 = $cagid . 'OUT2';
#	print $handle1 $head1,$seq1,$stand1,$que1;
#	print $handle2 $head2,$seq2,$stand2,$que2;
	$hash{$cagid}->{$readsid}->{$handle1} = $head1 . $seq1 . $stand1 . $que1;
	$hash{$cagid}->{$readsid}->{$handle2} = $head2 . $seq2 . $stand2 . $que2;
}

close FQ1;
close FQ2;
print STDERR datestr(), ": Reading readfiles finished .\n\n";

for my $tempid (@MGS){
    my $path = "$outdir/$tempid";
    -d $path || mkdir $path;
    my $handle1 = $tempid . 'OUT1';
    my $handle2 = $tempid . 'OUT2';
    my $f1 = $path . '/temp.' . $tempid . '.fq1.gz';
    my $f2 = $path . '/temp.' . $tempid . '.fq2.gz';
    -s $f1 && `rm $f1`;
    -s $f2 && `rm $f2`;
    open $handle1, ">:gzip", $f1 || die $!;
	open $handle2, ">:gzip", $f2 || die $!;
#	$handle2  = cacheout( ">:gzip", $f2) || die $!;
	for my $readID (keys %{$hash{$tempid}}){
		print $handle1 $hash{$tempid}->{$readID}->{$handle1};
		print $handle2 $hash{$tempid}->{$readID}->{$handle2};
	}
	print STDERR datestr(), ": Extracting $tempid reads finished .\n";
	close $handle1;
	close $handle2;
}

print STDERR datestr(), ": Extracting all MLGs reads finished .\n\n";
sub readsID2geneID{
	my ($bamfile) = @_;
	(-B "$bamfile") ? (open BAM, "$samtools view $bamfile | " || die $!) :
	(open BAM,  $bamfile || die $!);
	my (%reads_gene);
	while(<BAM>){
		/^@/ && next;
		chomp;
		my ($readsid, $geneid) = (split /\s+/)[0, 2];
		$geneid || next;
		$geneid eq '*' && next;
		$reads_gene{$readsid} = $geneid;
	}
	close BAM;
	return 	%reads_gene;
}

sub geneID2mgsID{
	my ($mgsfile) = @_;
	open MGS, $mgsfile || die $!;
	my %genes_mlg; my @mgs;
	<MGS>;
	while(<MGS>){
		chomp;
		my ($mgsid, $genes) = (split /\s+/)[0, -1];
		my @Genes = split /,/, $genes;
		push @mgs, $mgsid;
		for my $id (@Genes){
			$genes_mlg{$id} = $mgsid;
		}
	}
	close MGS;
	return (\%genes_mlg, \@mgs);
}

sub usage{
	print "\nUsage: script for extracting paired sequences for MGS/MLG according to bowtie mapping results\n";
	print "\n|---|      |---|      |---| 2019-11-28 |---|\n";
	print "\nExample: perl $0 -1 sample.fq1.gz -2 sample.fq2.gz -mgs mgs.cluster -bam readsmapping.sam -outdir ./\n\n";
	exit;
}

sub datestr{
    use POSIX qw(strftime);
    return strftime "%Y-%m-%d %H:%M:%S", localtime;
}
