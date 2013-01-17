quantile
========

quantile normalization script in Perl

Usage:

quantile [options] input [output]

  -file, -f
  
    Specify the input tabular file, which each row must contain the same number of columns.
    The first column is treated as unique ID for each row (except the first row, since it maybe
    the header line).
    
  -nohead, -n
  
    If set, the first line will not be treated as sample names line.
