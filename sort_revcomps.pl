#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Bio::SearchIO;
use Bio::SeqIO;
use Cwd;

my $help = 0;
my $targetsDir;
my $targetFasta;
my $outDir;

GetOptions  ("targets=s" => \$targetsDir,
             "targetFasta=s" => \$targetFasta,
             "out=s"      => \$outDir,
             "help|man" => \$help) || die "Couldn't get options with GetOpt::Long: $!\n";

if (!$targetsDir or !$targetFasta or !$outDir or $help) {
    die "Must supply --targets and --targetFasta and --out.\n";
}



my %targetFastaHash;
my $targetIn = Bio::SeqIO->new(-file=>$targetFasta,
                               -format=>'fasta');
while (my $seq = $targetIn->next_seq()) {
    $targetFastaHash{$seq->display_id()} = $seq;
}


my $startingDir = getcwd();
opendir(my $targetsDirFH, $targetsDir);
my @targetFastas = readdir($targetsDirFH);
closedir($targetsDirFH);

chdir ($targetsDir);
unless (-d $outDir) {
    mkdir $outDir;
}

foreach my $targetFasta (@targetFastas) {
    next if ($targetFasta !~ /\.fasta/);
    my $targetName;
    if ($targetFasta =~ /(.*)\.fasta/) {
        $targetName = $1;
    }
    my $indSeqOut = Bio::SeqIO->new(-file => ">targetCurrentlyProcessing.fasta",
                                    -format => 'fasta');
    $indSeqOut->write_seq($targetFastaHash{$targetName});
    
    print "Processing $targetFasta:\n";
    # First make a hash of all the sequences in that fasta file, and find the longest one
    my $seqIn = Bio::SeqIO->new(-file => $targetFasta,
                                -format => 'fasta');
    my %seqs;
    while (my $seq = $seqIn->next_seq()) {
        $seqs{$seq->display_id()} = $seq;
    }

    system("makeblastdb -in targetCurrentlyProcessing.fasta -dbtype nucl -out targetCurrentlyProcessing");
    system("blastn -db targetCurrentlyProcessing -task blastn -query $targetFasta -out all_bl2_targetCurrentlyProcessing.blast");
    
    # Manually delete database just to be safe
    system("rm targetCurrentlyProcessing*");

    
    
    # Now process the blast results file
    my $blastIn = Bio::SearchIO->new(-file => 'all_bl2_targetCurrentlyProcessing.blast',
                                     -format => 'blast');
    print "Strand results:\n";
    while (my $result = $blastIn->next_result()) {
        if ($result->num_hits() == 0) {
            print "No blast hit found for " . $result->query_name() . " in $targetFasta. MUST CHECK THIS ONE MANUALLY\n";
            next;
        }
        
        my $hit = $result->next_hit();
        my @strands = $hit->strand(); # Should return a two-member array, like (1,1) or (1,-1);
        # foreach my $strand (@strands) {
        #     print $strand . " ";
        # }
        # print "\n";
        if ($strands[0] == -1) {
            print "STRAND is minus for query!!\n";
            #code
        }
        if ($strands[1] == -1) {
            print "Revcomping " . $result->query_name() . " for " . $hit->name() . "\n";
            $seqs{$result->query_name} = $seqs{$result->query_name}->revcom();
        }

        
        
    }
    my $outSeqName = $outDir . "/$targetName" . "_revcomped.fasta";
    my $seqOut = Bio::SeqIO->new(-file   => ">$outSeqName",
                                 -format => 'fasta');
    foreach my $seqName (sort keys %seqs) {
        $seqOut->write_seq($seqs{$seqName});
    }
    
    unlink("all_bl2_targetCurrentlyProcessing.blast");
    print "\n\n";
}






















