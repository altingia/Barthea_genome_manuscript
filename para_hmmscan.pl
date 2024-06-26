#!/usr/bin/perl
use strict;
use Getopt::Long;

my $usage = <<USAGE;
Usage:
    perl $0 [options] protein.fasta > out.domtbl

    程序调用hmmscan将protein.fasta文件中的蛋白序列比对到hmm数据库中，并将hmmscan的domtbl结果进行过滤并输出到标准输出。
    调用hmmscan进行比对时，使用--cut_ga参数，即根据hmm文件中自带的GA类型的bit scores阈值进行过滤。在hmmscan文档中提示为：GA thresholds are generally considered to be the reliable curated thresholds defining family membership。使用--cut_ga方法能对不同的hmm模型使用不同的阈值进行过滤，但此时hmmscan不能使用evalue方法进行过滤。
    在得到hmmscan结果后，进一步使用evalue和coverage进行过滤。若一个蛋白序列和一个hmm模型有多个匹配区域，则计算多个匹配区域去冗余后的长度对hmm模型的覆盖率。

    --outformat
    当添加该参数时，输出简单的7列表格结果：（1）GeneID；（2）HMM accession，例如Pfam编号；（3）HMM Name；（4）E-value；（5）Score；（6）Coverage；（7）Description。默认情况下，输出的表格格式和hmmscan的domtbl一致且含有大量注释行信息。

    --evalue1 <float>    default: 1e-5
    --evalue2 <float>    default: 1e-3
    --hmm_length <int>    default: 80
    设置hmmsacn命令的evalue阈值。当hmm模型长度低于--hmm_length参数值时，使用--evalue2阈值；当hmm模型长度不低于--hmm_length参数值时，使用--evalue1阈值。此外，此处的E-value阈值对象是hmmscan结果中full sequence，而不是domain。

    --coverage <float>    default: 0.25
    设置蛋白序列对hmm模型的覆盖度阈值。

    --no_cut_ga    default: none
    默认设置下，运行hmmscan命令时使用了--cut_ga参数。但当输入的HMM数据库并不是PFam-A时，其HMM数据库中并不包含GA信息，会导致hmmscan运行失败。此时，需要添加该参数，则不会在hmmscan命令中使用--cut_ga参数，从而避免程序运行失败。

    --chunk <in>    default: 100
    设置每个数据块的protein序列条数。程序将protein.fasta序列从头到尾的分割成多份，每100条相邻的序列分配到一个fasta文件中；每个fasta文件写出一条hmmscan命令。

    --cpu <int>    default: 4
    程序调用ParaFly对hmmscan命令进行并行化计算，此参数传递给ParaFly，表示并行运行数目。

    --hmmscan_cpu <int>    default: 2
    传递给hmmscan命令的--cpu参数，用于设置单个hmmscan命令消耗的CPU线程。

    --hmm_db <string>    default: "/opt/biosoft/hmmer-3.3.2/Pfam-AB.hmm"
    设置hmmscan比对的数据库路径。

    --tmp_prefix <string>    default: "hmmscan"
    设置临时文件或文件夹前缀。默认设置下，程序生成hmmscan.command.list, hmmscan.tmp/等临时文件或目录。

USAGE
if (@ARGV==0){die $usage}

my ($outfmat, $evalue1, $evalue2, $hmm_length, $coverage, $no_cut_ga, $chunk, $cpu, $hmmscan_cpu, $hmm_db, $tmp_prefix);
GetOptions(
    "outformat!" => \$outfmat,
    "evalue1:f" => \$evalue1,
    "evalue2:f" => \$evalue2,
    "hmm_length:i" => \$hmm_length,
    "coverage:f" => \$coverage,
    "no_cut_ga!" => \$no_cut_ga,
    "chunk:i" => \$chunk,
    "cpu:i" => \$cpu,
    "hmmscan_cpu:i" => \$hmmscan_cpu,
    "hmm_db:s" => \$hmm_db,
    "tmp_prefix:s" => \$tmp_prefix,
);
$evalue1 ||= 1e-5;
$evalue2 ||= 1e-3;
$hmm_length ||= 80;
$coverage ||= 0.25;
$chunk ||= 100;
$cpu ||= 4;
$hmmscan_cpu ||= 2;
$hmm_db ||= "/opt/biosoft/hmmer-3.3.2/Pfam-AB.hmm";
$tmp_prefix ||= "hmmscan";

mkdir "$tmp_prefix.tmp" unless -e "$tmp_prefix.tmp";
open IN, $ARGV[0] or die "Can not open the input file: $ARGV[0]\n$!\n";
open CMD, ">", "$tmp_prefix.command.list"  or die "Cannot create file $tmp_prefix.command.list\n$!\n";
my ($fasta, $number, $chunk_number, @chunk);
while (<IN>) {
    if (m/^>/) {
        $number ++;
        if ($number > $chunk) {
            $chunk_number ++;
            push @chunk, "$tmp_prefix.tmp/chunk.$chunk_number";
            open OUT, ">", "$tmp_prefix.tmp/chunk.$chunk_number.fasta" or die "Can not create file $tmp_prefix.tmp/chunk.$chunk_number.fasta\n$!\n";
            print OUT $fasta;
            my $out = "hmmscan --cpu $hmmscan_cpu --cut_ga -o $tmp_prefix.tmp/chunk.$chunk_number.txt --tblout $tmp_prefix.tmp/chunk.$chunk_number.tbl --domtblout $tmp_prefix.tmp/chunk.$chunk_number.domtbl --pfamtblout $tmp_prefix.tmp/chunk.$chunk_number.pfamtbl $hmm_db $tmp_prefix.tmp/chunk.$chunk_number.fasta\n";
            $out =~ s/--cut_ga// if $no_cut_ga;
            print CMD $out;
            close OUT;
            $number  = 1;
            $fasta = "";
        }
    }
    $fasta .= $_;
}
close IN;
if ($fasta) {
    $chunk_number ++;
    push @chunk, "$tmp_prefix.tmp/chunk.$chunk_number";
    open OUT, ">", "$tmp_prefix.tmp/chunk.$chunk_number.fasta" or die "Can not create file $tmp_prefix.tmp/chunk.$chunk_number.fasta\n$!\n";
    my $out = "hmmscan --cpu $hmmscan_cpu --cut_ga -o $tmp_prefix.tmp/chunk.$chunk_number.txt --tblout $tmp_prefix.tmp/chunk.$chunk_number.tbl --domtblout $tmp_prefix.tmp/chunk.$chunk_number.domtbl --pfamtblout $tmp_prefix.tmp/chunk.$chunk_number.pfamtbl $hmm_db $tmp_prefix.tmp/chunk.$chunk_number.fasta\n";
    $out =~ s/--cut_ga// if $no_cut_ga;
    print CMD $out;
    print OUT $fasta;
    close OUT;
}
close CMD;

my $cmdString = "ParaFly -c $tmp_prefix.command.list -CPU $cpu &> /dev/null";
print STDERR "CMD: $cmdString\n";
(system $cmdString) == 0 or die "Failed to execute: $cmdString\n";

print "GeneID\tHMM_Accession\tHMM_Name\tE-value\tScore\tCoverage\tDescription\n" if $outfmat;
foreach (@chunk) {
    open IN, "$_.domtbl" or die "Can not open file $_.domtbl\n$!\n";
    my ($out_head, $out_tail, %info);
    $out_head = <IN>; $out_head .= <IN>; $out_head .= <IN>;
    while (<IN>) {
        if (m/^#/) {
            $out_tail .= $_;
        }
        else {
            @_ = split /\s+/;
            $info{$_[3]}{$_[0]} .= $_;
        }
    }
    close IN;

    my %annot;
    foreach my $geneID (keys %info) {
        foreach my $name (keys %{$info{$geneID}}) {
            my $line = $info{$geneID}{$name};
            chomp($line);
            my (@match_region_of_hmm, $hmm_len, $evalue, $accession, $score, $description);
            foreach (split /\n/, $line) {
                chomp;
                @_ = split /\s+/;
                push @match_region_of_hmm, "$_[15]\t$_[16]";
                $hmm_len = $_[2];
                $evalue = $_[6];
                $accession = $_[1];
                $score = $_[7];
                foreach (1..22) { shift @_; }
                $description = join " ", @_;
            }
            my $match_length = &match_length(@match_region_of_hmm);
            my $hmm_coverage = $match_length / $hmm_len;
            next if $hmm_coverage < $coverage;
            next if $evalue > $evalue2;
            next if ($hmm_len >= $hmm_length && $evalue > $evalue1);
            $annot{$geneID}{$name}{"evalue"} = $evalue;
            $hmm_coverage = int($hmm_coverage * 10000 + 0.5) / 100;
            $annot{$geneID}{$name}{"out"} = "$geneID\t$accession\t$name\t$evalue\t$score\t$hmm_coverage\%\t$description";
        }
    }

    print $out_head unless $outfmat;
    foreach my $geneID (sort keys %annot) {
        foreach my $name (sort {$annot{$geneID}{$a}{"evalue"} <=> $annot{$geneID}{$b}{"evalue"}} keys %{$annot{$geneID}}) {
            if ($outfmat) {
                print $annot{$geneID}{$name}{"out"} . "\n";
            }
            else {
                print $info{$geneID}{$name};
            }
        }
    }
    print $out_tail unless $outfmat;
}

sub match_length {
    my @inter_sorted_site;
    foreach (@_) {
        my @aaa = $_ =~ m/(\d+)/g;
        @aaa = sort { $a <=> $b } @aaa;
        push @inter_sorted_site, "$aaa[0]\t$aaa[1]";
    }
    @inter_sorted_site = sort { $a <=> $b } @inter_sorted_site;

    my $out_site_number;
    my $former_region = shift @inter_sorted_site;
    my @aaa = $former_region =~ m/(\d+)/g;
    $out_site_number += ($aaa[1] - $aaa[0] + 1);
    foreach (@inter_sorted_site) {
        my @former_region = $former_region =~ m/(\d+)/g;
        my @present_region = $_ =~ m/(\d+)/g;

        if ($present_region[0] > $former_region[1]) {
            $out_site_number += ($present_region[1] - $present_region[0] + 1);
            $former_region = $_;
        }
        elsif ($present_region[1] > $former_region[1]) {
            $out_site_number += ($present_region[1] - $former_region[1]);
            $former_region = $_;
        }
        else {
            next
        }
    }
    return $out_site_number;
}
