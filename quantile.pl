#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

my $version = 1.0;
my ($nohead, $file, $output, $help, $sep, $get_version);

GetOptions( 'nohead|n'   => \$nohead,
            'file|f=s'   => \$file,
            'output|o=s' => \$output,
            'help|h|?'   => \$help,
            'sep|s'      => \$sep,
            'version|v'  => \$get_version,
);

pod2usage(1) if $help;

if ($get_version) {
    print "quantile version $version\n";
    exit;
}

if (@ARGV) {
    $file = shift @ARGV unless $file;
    if (@ARGV) {
        $output = shift @ARGV unless $output;
    }
    die "Unknow argument(s) @ARGV given!\n" if @ARGV;
}

pod2usage(1) unless ($file);

$sep = '\t' unless ($sep);

unless ($output) {
    $output = $file;
    if ($output =~ /\./) { $output =~ s/\.[^.]+$/.quantile/ }
    else { $output .= '.quantile' }
}

my (@samples, @genes, %data, %genes);
my $ncol; # number of columns for data

open my $in, '<', $file || die "Can't open $file for reading. $!\n";
if ($nohead) {
    chomp(my $line = <$in>);
    my @t = split /$sep/, $line;
    $ncol = $#t;
    @samples = (1 .. $ncol); # indexes as default sample name
    push @genes, $t[0];
    $genes{$t[0]} = undef;
    for (1 .. $ncol) {
        push @{$data{$_}}, $t[$_]; # use sample names as keys
    }
}
else {
    chomp(my $head = <$in>);
    @samples = split /$sep/, $head;
    chomp(my $line = <$in>);
    my @t = split /$sep/, $line;
    $ncol = $#t;
    shift @samples if @samples > $ncol;
    my $gene = shift @t;
    push @genes, $gene;
    $genes{$gene} = undef;
    for (@samples) {
        my $d = shift @t;
        push @{$data{$_}}, $d;
    }
}

while (<$in>) {
    chomp;
    my @t = split /$sep/;
    die "Each row should have the same number of columns.\n"
      . "$ncol columns (exclude gene name) detected for most rows,\n"
      . "but $#t cloumns detected for line $.:\n$_\n"
      if $#t != $ncol;
    my $gene = shift @t;
    next if exists $genes{$gene};
    $genes{$gene} = undef;
    push @genes, $gene;
    for (@samples) {
        my $d = shift @t;
        push @{$data{$_}}, $d;
    }
}
close $in;

my $nrow = $#genes;

my %index;
my @sum; # sum for each row
for my $k (@samples) {
    my @t = @{$data{$k}};
    my %sort_index;
    my $i = 0;
    my @sort = sort {$a <=> $b} @t;
    for (0 .. $nrow) {
        if (exists $sort_index{$sort[$_]}) { $i++ }
        else { $sort_index{$sort[$_]} = $i++; }
        $sum[$_] += $sort[$_]; # summary each row
    }
    my @new_index = @sort_index{@t};
    for (0 .. $nrow) {
        if (exists $sort_index{$sort[$_]}) { $i++ }
        else { $sort_index{$sort[$_]} = $i++; }
        $sum[$_] += $sort[$_]; # summary each row
    }
    my @new_index = @sort_index{@t};
    for (@genes) { # set the new index to each gene as a single row
        my $i = shift @new_index;
        push @{$index{$_}}, $i;
    }
}

@sum = map { sprintf("%.2f", $_ / $ncol) } @sum; # mean value

open STDOUT, '>', $output || die "Can't write to $output. $!\n";

print join("\t", '', @samples) . "\n" unless $nohead;
for (@genes) {
    print join("\t", $_, @sum[@{$index{$_}}]), "\n";
}


=pod

=head1 NAME

quantile    Do quantile normalization on given file data.

=head1 SYNOPSIS

=item B<quantile> [options] input [output]

B<quantile> can do quantile normalization for a tabular file, which typically a
macroarray expression profile file, which file row are array IDs (sample IDs),
followed by rows begin with the gene IDs. Each row should contain the same
number of columns, which should equal the number of sample IDs (exclude gene ID
itself).

=head2 OPTIONS

=item B<-file, -f>

Specify the tabular file to do quantile normalization. Needed. You can omit
B<-f>, in this case, the first unpaired argument will be treated as file.

=item B<-nohead, -n>

If set, no header line (sample IDs) contained in the file. Default is not set.

=item B<-output, -o>

Specify the output file name. Default output file name is derived from B<-f>,
which will replace the suffix with '.quantile'. You can omit B<-o>, in this
case, the first unpaired argument after file is set will be treated as output.

=item B<-sep, -s>

Specify the separator for each column. Default is '\t'.

=item B<-help, -h, -?>

Display this help information.

=head2 B<EXAMPLES>

quantile input.txt output.quantile # set input and output name
quantile input.txt -sep , # set input and use default output, comma-separated

=head2 B<AUTHOR>

                    Zhijun HAN (<hangeneral@126.com>)
                            Jan 5th, 2013
                    Group of Epigenome Biology, PICB
                           Shanghai, China
