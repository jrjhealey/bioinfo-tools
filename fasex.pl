#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;

my $string = pop;
my $seqout = Bio::SeqIO->new( -format => 'Fasta', -fh => \*STDOUT );

for my $inFile (@ARGV) {
    my $seqin = Bio::SeqIO->new( -format => 'Fasta', -file => $inFile );

    while ( my $seq = $seqin->next_seq() ) {

        if (   $seq->id =~ m/\Q$string\E/i
            or $seq->description =~ m/\Q$string\E/i )
        {
            $seqout->write_seq($seq);
        }
    }
}
