#!/usr/bin/perl -w
use strict;
use warnings;
use List::Util qw(sum);
use Getopt::Long;

my %opt = ("num_each_group" => 1);

GetOptions(\%opt, 'num_each_group|n:s');

@ARGV == 2 || die &usage();

my ($group, $pcfile) = @ARGV;

my %sam_group = sample_group($group);

my %group_sample_data;

-s $pcfile || die "No $pcfile, exit...\n";

open PC, $pcfile;
<PC>;
while(<PC>){
	chomp;
	my @ll = split /\t/;
	my $samID = shift @ll;
	$sam_group{$samID} || warn "No group assigned for $samID, escape it.\n";
	$sam_group{$samID} || next;
	$group_sample_data{$sam_group{$samID}}->{$samID} = \@ll;
}
close PC;

for my $g (sort {$a cmp $b} keys %group_sample_data){
	my %temp_g_sample_data = %{$group_sample_data{$g}};
	my @temp_data;
	for my $s (keys %temp_g_sample_data){
		push @temp_data, $temp_g_sample_data{$s};
	}
	my $mid_point_value = middle_point(@temp_data);
	my %temp_sample_middle_dist;
	for my $s (keys %temp_g_sample_data){
		$temp_sample_middle_dist{$s} = euclidean_dist($mid_point_value, $temp_g_sample_data{$s});
		print "$s-$g\t$temp_sample_middle_dist{$s}\n";
	}
	my @min_value = (sort {$a <=> $b } values %temp_sample_middle_dist )[0..$opt{num_each_group}-1];
	my @min_value_sam = ();
	for my $m (@min_value){
		my @temp_min_value_sam = grep { $temp_sample_middle_dist{$_} == $m } keys %temp_sample_middle_dist;
		@min_value_sam = (@min_value_sam, @temp_min_value_sam);
	}
	
	print "# ", join(",", @min_value_sam), "\n\n";
}




sub middle_point{
	my (@points) = @_;
	my $dim_num = scalar(@{$points[0]});
	my @middle_value;
	for my $i (0 .. $dim_num-1){
		my $d = mean( map {$_->[$i]} @points );
		push @middle_value, $d;
	}
	return \@middle_value;
}

sub euclidean_dist{
	my ($p1, $p2) = @_;
	my @p11 = @{$p1}; my @p22 = @{$p2};
	my $squar_sum;
	for my $i (0..$#p11){
		$squar_sum += abs($p11[$i] - $p22[$i]) ** 2;
	}
	return sqrt($squar_sum);
}

sub mean{
	my (@num) = @_;
	return sum(@num) / scalar(@num);
}

sub sample_group{
	my ($gl) = @_;
	open GL, $gl || die $!;
	my %sample_group;
	while(<GL>){
		chomp;
		my @ll = split;
		$sample_group{$ll[0]} = $ll[1];
	}
	return %sample_group;
}


sub usage{
	print "Picking closest point between sample and middle, according to euclidean distance\n\n";
	print "perl $0 group.list matrix.txt -n 10 > selected.mat.txt\n\n";
	print "--num_each_group <integer>   sample number selected in a group\n\n";
	print "group.list: sample_name\\tgroup\n\n";
	exit;
}
