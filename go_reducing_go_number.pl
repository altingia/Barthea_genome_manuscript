#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    $0 gene_ontology_edit.obo go.annot > go_reduced.annot

    程序通过gene_ontology_edit.obo文件解析GO term之间的is_a关系；然后对单个基因的GO注释编号之间的关系进行解析，去除上层的GO编号，保留注释更精细的子级编号。

    程序读取输入的go.annot文件的前两列数据。第一列，基因ID；第二列，GO编号。若基因有多个GO编号注释，则用多列表示。程序同样输出两列，并额外输出第三列对GO的描述性信息。

USAGE
if (@ARGV==0){die $usage}

# 分析 gene_ontology_edit.obo 文件
$/ = "\n\n";
open IN, $ARGV[0] or die "Can not open file $ARGV[0], $!\n";
my (%alias, %isa, %name);
while (<IN>) {
    next if m/is_obsolete: true/;
    if (m/id: (GO:\d+).*\nname: (.*)/) {
        my ($id, $name) = ($1, $2);
        $name{$id} = $name;

        my @alias =  m/alt_id: (GO:\d+)/g;
        foreach (@alias) {
            $alias{$_} = $id;
            $name{$_} = $name;
        }

        my @is_a = m/is_a: (GO:\d+)/g;
        foreach (@is_a) {
            $isa{$_}{$id} = 1;
        }
    }
}
close IN;
$/ = "\n";

open IN, $ARGV[1] or die "Can not open file $ARGV[1], $!\n";
my ($gene, $annot);
while (<IN>) {
    @_ = split /\t/;
    if ($_[0] eq $gene) {
        $annot .= $_;
    }
    else {
        my @out = &reduce_go($annot);
        foreach (@out) {
            print "$gene\t$_\t$name{$_}\n";
        }
        $gene = $_[0];
        $annot = $_;
    }
}
close IN;
my @out = &reduce_go($annot);
foreach (@out) {
    print "$gene\t$_\t$name{$_}\n";
}


sub reduce_go {
    my @line = split /\n/, $_[0];
    my %go;
    foreach (@line) {
        @_ = split /\t/;
        if (exists $alias{$_[1]}) {
            $go{$alias{$_[1]}} = 1;
        }
        else {
            $go{$_[1]} = 1 if exists $name{$_[1]};
        }
    }

    my @go = keys %go;
    foreach my $go (@go) {
        my %progeny = &get_all_progeny($go);
        foreach (keys %progeny) {
            if (exists $go{$_}) {
                delete $go{$go};
                last;
            }
        }
    }

    my @out = sort keys %go;
    return @out;
}

sub get_all_progeny {
    my $query_id = shift @_;
    my %target;
    my @list = keys %{$isa{$query_id}};
    while (@list) {
        my $one = shift @list;
        $target{$one} = 1;
        push @list, keys %{$isa{$one}};
    }
    return %target;
}
