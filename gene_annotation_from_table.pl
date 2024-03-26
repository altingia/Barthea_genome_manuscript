#!/usr/bin/perl
use strict;
use Getopt::Long;

my $usage = <<USAGE;
Usage:
    $0 [options] annot.tab > annot.txt

    程序用提取表格格式的注释结果，将注释结果转换为2列式的结果：第一列，geneID；第二列，注释结果，多个注释结果之间用“分号空格”分割。
    程序输入信息有三列：第一列，geneID；第二列，注释数据库的编号；第三列，对编号的描述。

    程序运行示例：
    cut -f 1,3,4 cog.annot | $0 - > COG.txt
    cut -f 2,5,9 gene_association.gaf2 | $0 - > GO.txt
    grep IPR interpro.tsv | cut -f 1,12,13 | $0 - > Interpro.txt

USAGE
if (@ARGV==0){die $usage}

my ($gene, $annot);
while (<>) {
    next if m/^#/;
    @_ = split /\t/;
    next if @_ ne 3;

    if ($_[0] eq $gene) {
        $annot .= $_;
    }
    else {
        my $out = &get_annot($annot);
        print "$gene\t$out\n" if $gene;
        $gene = $_[0];
        $annot = $_;
    }
}
my $out = &get_annot($annot);
print "$gene\t$out\n";


sub get_annot {
    my %line;
       foreach (split /\n/, $_[0]) {
        $line{$_} = 1;
    }

    my @out;
    foreach (sort keys %line) {
        @_ = split /\t/;
        $_[2] =~ s/;\s+/;/g;
        push @out, "$_[1]: $_[2]";
    }

    my $out = join "; ", @out;
    return $out;
}
