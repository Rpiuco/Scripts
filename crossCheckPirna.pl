#!usr/bin/perl

use strict;
use warnings;
use Getopt::Std;

my %options=();
getopts("f:", \%options);

my $pirnaDBFile = 'hsa_pirna.fa';
my $tcgaSequenceFile = $options{'f'}.'.bam.q20.fa';

my %pirnaDBData;
my $pirnaExistTcga;

$/ = '>';
open(FHPDBF, $pirnaDBFile);
<FHPDBF>;

while(my $rowData = <FHPDBF>) 
{

    chomp($rowData);
    my @dataArrayOne = split("\n", $rowData);

    $pirnaDBData{$dataArrayOne[1]} = $dataArrayOne[0];

}
close(FHPDBF);

$/ = '>';
open(FHTSF, $tcgaSequenceFile);
<FHTSF>;

while(my $rowDataTcga = <FHTSF>) 
{

    chomp($rowDataTcga);
    my @dataArrayOne = split("\n", $rowDataTcga);

    if (defined($pirnaDBData{$dataArrayOne[1]})) {
        $pirnaExistTcga .= $pirnaDBData{$dataArrayOne[1]}."\n";
    }

}
close(FHTSF);

open(FHWET, ">pirnaExistTcga.".$options{'f'}.".txt");
print FHWET $pirnaExistTcga;
close(FHWET);