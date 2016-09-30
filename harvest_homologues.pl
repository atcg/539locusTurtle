#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;
use Cwd;
use Bio::Tools::Run::Alignment::Muscle;
use Data::Dumper;

my $help = 0;
my $outDir;
my $RBBHdir;
my $alignmentDir;

GetOptions  ("RBBHs=s"      => \$RBBHdir,
             "alignments=s" => \$alignmentDir,
             "out=s"        => \$outDir,
             "help|man"     => \$help) || die "Couldn't get options with GetOpt::Long: $!\n";

if (!$RBBHdir or !$alignmentDir or !$outDir or $help) {
    die "Must supply --RBBHs and --alignments and --out.\n";
}

my $startingDir = getcwd();
opendir (my $RBBHdirFH, $RBBHdir);
my @RBBHfiles = readdir $RBBHdirFH or die "Couldn't read files in $RBBHdir: $!\n";
closedir $RBBHdirFH;

chdir $RBBHdir or die "Couldn't change directories into $RBBHdir:$!\n";

my %allSeqsHash;
my @targetHomologues;

foreach my $RBBHfile (@RBBHfiles) {
    if ($RBBHfile =~ /abyss_51_(.*)_RBBHs.fasta/) {
        # $1 is the sample name
        my $sampleName = $1;
        my $seqIn = Bio::SeqIO->new(-file => $RBBHfile,
                                    -format => 'fasta');
        while (my $seq = $seqIn->next_seq()) {
            if ($seq->display_id =~ /\d+_(\S+)/) {
                # $1 should now be the name of the homologous target
                my $newName = $sampleName . "_" . $1;
                $seq->display_id($newName);
                push (@targetHomologues, $1);
                $allSeqsHash{$1}{$sampleName} = $seq;
            }
        }
    }
}
chdir $startingDir or die "Couldn't change back to $startingDir: $!\n";
unless (-d $outDir) {
    mkdir $outDir or die "Couldn't make $outDir: $!\n";
}


chdir $outDir or die "Couldn't change into $outDir: $!\n";
foreach my $targetHomologue (@targetHomologues) {
    my $fileName = "$targetHomologue.fasta";
    my $seqOut = Bio::SeqIO->new(-file => ">$fileName",
                                 -format => 'fasta');
    foreach my $sampleTarget (sort keys %{$allSeqsHash{$targetHomologue}}) {
        # print $allSeqsHash{$targetHomologue}{$sampleTarget}->display_id() . "\n";
        # $sampleTarget here should be all the hash keys for the target, meaning the names of all the biological samples
        # $allSeqsHash{$targetHomologue}{$sampleTarget} should be the relevant Bio::Seq object
        $seqOut->write_seq($allSeqsHash{$targetHomologue}{$sampleTarget});        
    }
}
# print Dumper(\%allSeqsHash);


























