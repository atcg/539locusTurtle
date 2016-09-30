#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Cwd;

my $help = 0;
my $RBBHsDir;
my $readsDir;
my $outDir;

GetOptions  ("RBBHs=s"      => \$RBBHsDir,
             "reads=s"      => \$readsDir,
             "out=s"        => \$outDir,
             "help|man" => \$help) || die "Couldn't get options with GetOpt::Long: $!\n";

if (!$RBBHsDir or !$readsDir or !$outDir or $help) {
    die "Must supply --RBBHs and --reads and --out.\n";
}
# perl ../../../turtlebaitspipeline/map_reads_to_RBBHs.pl --RBBHs ./ --reads constituentReads/ --out readMapping > readmappingSummary 2>&1
my $startingDir = getcwd();

opendir (my $RBBHdirFH, $RBBHsDir);
my @RBBHs = readdir($RBBHdirFH);
closedir ($RBBHdirFH);

unless (-d $outDir) {
    mkdir $outDir;
}

chdir $outDir;

foreach my $RBBHfile (@RBBHs) {
    next unless ($RBBHfile =~ /RBBHs\.fasta/);
    my $reference = "../$RBBHfile";
    my $sample;
    my $pooling;
    if ($RBBHfile =~ /abyss.*\_([1248]+x)\_(.*)\_RBBHs.fasta/) {
        #$1 should be something like 1x or 4x or 1248x
        #$2 should be something like HBS_108583_GCAGGCGT or USNM_520644
        $pooling = $1;
        $sample = $2;
    }
    my $reads1 = "../$readsDir" . $pooling . "_" . $sample . ".un1.fastq";
    my $reads2 = "../$readsDir" . $pooling . "_" . $sample . ".un2.fastq";
    my $readsSingles = "../$readsDir" . $pooling . "_" . $sample . "_joined_and_both_singles.fastq";
    my $singlesSamFile = $pooling . "_" . $sample . ".singles.sam";
    my $pairedSamFile = $pooling . "_" . $sample . ".paired.sam";
    my $singlesBamFile = $pooling . "_" . $sample . ".singles.bam";
    my $pairedBamFile = $pooling . "_" . $sample . ".paired.bam";
    my $mergedBamFile = $pooling . "_" . $sample . ".merged.bam";
    system("bwa index $reference");
    system("bwa mem -t 12 $reference $readsSingles > $singlesSamFile");
    system("bwa mem -t 12 $reference $reads1 $reads2 > $pairedSamFile");
    system("samtools view -b -S $singlesSamFile > $singlesBamFile");
    system("samtools view -b -S $pairedSamFile > $pairedBamFile");
    system("samtools merge $mergedBamFile $singlesBamFile $pairedBamFile");
    
    # Mark duplicates and use mpileup
    my $cleanedBam = $pooling . "_" . $sample . ".merged.cleaned.bam";
    my $sortedBam = $pooling . "_" . $sample . ".merged.cleaned.sorted.bam";
    my $markDupsBam = $pooling . "_" . $sample . ".merged.cleaned.sorted.markDups.bam";
    my $markDupsMetrics = $pooling . "_" . $sample . ".merged.sorted.cleaned.markDups.metrics";
    my $pileupFile = $pooling . "_" . $sample . ".mpileup";
    system("java -jar ~/bin/picard/picard-tools/CleanSam.jar I=$mergedBamFile O=$cleanedBam");    
    system("java -jar ~/bin/picard/picard-tools/AddOrReplaceReadGroups.jar I=$cleanedBam O=$sortedBam SORT_ORDER=coordinate RGPL=illumina RGPU=Test RGLB=Lib1 RGID=$sample RGSM=$sample VALIDATION_STRINGENCY=LENIENT");
    system("java -jar ~/bin/picard/picard-tools/MarkDuplicates.jar I=$sortedBam O=$markDupsBam METRICS_FILE=$markDupsMetrics MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=250 ASSUME_SORTED=true REMOVE_DUPLICATES=false");
    system("samtools mpileup $markDupsBam > $pileupFile");

    print "Stats for $RBBHfile:\n";
    system("samtools flagstat $markDupsBam");
    print "\n\n\n";
    unlink ($singlesSamFile, $pairedSamFile);
}










