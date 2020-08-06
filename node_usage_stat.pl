#!/usr/bin/perl -w
use strict;
use warnings;

my @nodes = &get_node_name();

printf("%-30s %-30s %-40s %-25s\n", ('#node', 'left | {total} [CPU]', 'virtual [theory] (used) {total} [mem]', 'queue'));
for my $n (@nodes){
	printf("%-30s %-30s %-40s %-25s\n", node_usage_stat($n));
}

sub get_node_name{
	my $cont = `qselect -U $ENV{USER}`;
	my @nds;
	my @queueNode = split /\s+/, $cont;
	for my $s (@queueNode){
		my $n = (split /\@/, $s)[-1];
		push @nds, $n;
	}
	return @nds;
}

sub node_usage_stat{
	my ($node) = @_;
	open CONT, "qhost -q -F vf,p,m_core -h $node |";
	<CONT>; <CONT>; <CONT>;
	my $temp = <CONT>;
	my ($nodename, $total_mem, $used_mem) = (split /\s+/, $temp)[0, -4, -3];
	$used_mem = $used_mem eq '-' ? '0G' : $used_mem;
	$total_mem = transform_num($total_mem);
	$used_mem = transform_num($used_mem);
	my $left_mem = sprintf("%.1f", $total_mem - $used_mem);
	my ($left_cpu_num, $virtual_free, $queue, $total_cpu, $warn);
	while(<CONT>){
		s/^\s+//g;
		if(/num_proc=(\S+)/){
			$left_cpu_num = sprintf("%.0f", $1);
		}elsif(/virtual_free=(\S+)/){
			$virtual_free = transform_num($1);
		#tjsmp07_1024.q       BIP   0/10/64
		}elsif(/BIP/){
			/admin.q/ && next;
			my @temp = split /\s+/;
			$total_cpu = (split /\//, $temp[2])[-1];
			$total_cpu <= 1 && next;
			$queue = $temp[0];
			$warn = scalar(@temp) > 3 ? $temp[-1] : "";
			$nodename = $warn ? $nodename . '[' .  $warn . ']' : $nodename;
		}
	}
	#NodeName\tQueue\tTotMem\tUsedMem\tVirMem\tNCPU\tLCPU\n
	$virtual_free ||= 0;
	my $mem = "$virtual_free [$left_mem] ($used_mem) {$total_mem}";
	my $nod = "$left_cpu_num | {$total_cpu}";
	return($nodename, $nod, $mem, $queue);

}


sub transform_num{
	my ($strs) = @_;
	if($strs =~ /(\S+)(G)/i){
		return sprintf("%.1f", $1);
	}elsif($strs =~ /(\S+)M/i){
		return sprintf('%.1f', $1/1024);
	}
}

