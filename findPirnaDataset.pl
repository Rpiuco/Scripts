########################################################
# FILE: findPirnaDataset
# USAGE: 
# DESCRIPTION: Find piRNA sequences of fasta file on .sam files
# OPTIONS: - s: the SAM alignment file
#          - f: the FASTA sequence file
# REQUIREMENTS: piRNA sequences in fasta format
# BUGS: --- 
# NOTES: ---  
# AUTHOR: Ricardo Piuco 
# COMPANY: Bioinformatics Lab - Hospital Sirio Libanes 
# VERSION: 1.0 
# CREATED: Mon Set 05 10:19:20 BRT 2016 
# REVISION: --- 
# LICENSE: GPL 
########################################################

#!usr/bin/perl

use warnings;
use strict;
use Getopt::Std;

# Get options
my %options=();
getopts("f:s:", \%options);

# Check if options are available
if (defined($options{'s'}) && defined($options{'f'})) {

    my $samArchive = $options{'s'};
    my $fastaFile = $options{'f'};

    # Fill a hash with fasta file data
    $/ = '>';
    open(FFD, $fastaFile) or die "Problem Fast File";
    <FFD>;

    my %fastaData; 
    my %samData;
    my $listPirnaExist; 

    while (my $rowFastaData = <FFD>) {

        chomp($rowFastaData);

        my ($pirnaCode, $pirnaSequence) = split("\n", $rowFastaData);

        $fastaData{$pirnaSequence} = $pirnaCode;

        undef($pirnaCode);
        undef($pirnaSequence);

    }

    close(FFD);

    # Fill a hash with sam data file
    $/ = "\n";
    open(FSD, $samArchive) or die "Problem Sam File";

    while (my $rowSamData = <FSD>) {

        if (substr($rowSamData, 0, 1) eq '@') {
            next;
        }

        chomp($rowSamData);

        my @samDataRaw = split("\t", $rowSamData);

        if (!exists($samData{$samDataRaw[9]})) {
            $samData{$samDataRaw[9]} = $samDataRaw[0];
        }
       
        undef(@samDataRaw);
    }

    close(FSD);

    # Check if the piRNA sequence is present on sam file
    foreach my $pirnaExists (keys %fastaData) {

        if (exists($samData{$pirnaExists})) {
            $listPirnaExist .= $fastaData{$pirnaExists}."\n".$pirnaExists."\n";
        }

    }

    # save the piRNA list that appeared on the SAM file
    open(my $newFile, '>pirnaExistsDataset.fa');
    print $newFile $listPirnaExist;
    close($newFile)

} else {
    print "You need to specify the fasta file (-f) and the sam file (-s).";
}


