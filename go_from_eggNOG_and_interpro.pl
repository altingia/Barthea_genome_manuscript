#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 eggNOG_hmmer.emapper.annotations interpro.tsv > go.annot

USAGE
if (@ARGV==0) { die $usage }

open IN, '<', $ARGV[0] or die "Can not open file $ARGV[0], $!\n";
my (%go, %go_eggNOG, %go_interpro);
while (<IN>) {
    my $gene_id = $1 if m/^(\S+)/;
    if (my @go = m/(GO:\d+)/g) {
        foreach (@go) {
            $go{$gene_id}{$_} = 1;
            $go_eggNOG{$gene_id}{$_} = 1;
        }
    }
}
close IN;

open IN, '<', $ARGV[1] or die "Can not open file $ARGV[1], $!\n";
while (<IN>) {
    my $gene_id = $1 if m/^(\S+)/;
    if (my @go = m/(GO:\d+)/g) {
        foreach (@go) {
            $go{$gene_id}{$_} = 1;
            $go_interpro{$gene_id}{$_} = 1;
        }
    }
}
close IN;

my ($num_eggNOG, $num_interpro, $num_common, $num_total);
my ($number_eggNOG, $number_interpro, $number_common, $number_total);
foreach my $gene_id (sort keys %go) {
    $number_total ++;
    $number_eggNOG ++ if exists $go_eggNOG{$gene_id};
    $number_interpro ++ if exists $go_interpro{$gene_id};
    $number_common ++ if (exists $go_eggNOG{$gene_id} && exists $go_interpro{$gene_id});
    foreach (sort keys %{$go{$gene_id}}) {
        $num_total ++;
        $num_eggNOG ++ if $go_eggNOG{$gene_id}{$_};
        $num_interpro ++ if $go_interpro{$gene_id}{$_};
        $num_common ++ if ($go_eggNOG{$gene_id}{$_} && $go_interpro{$gene_id}{$_});
        print "$gene_id\t$_\n";
    }
}
print STDERR "\tTotal\teggNOG\tInterpro\tCommon\n";
print STDERR "GO:\t$num_total\t$num_eggNOG\t$num_interpro\t$num_common\n";
print STDERR "Gene:\t$number_total\t$number_eggNOG\t$number_interpro\t$number_common\n";
