#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Text::CSV;
use Data::Dumper;
use Bio::SeqIO;
use Bio::Seq;

my $help = 0;
my $pileupsDir;
my $contigsDir;
my $outDir;
my $windowSize;
my $minimumDepth;
my $minimumContigLength; # Default to 200


# As it stands now both the pileups files and contigs files must be in the current working
# directory when the script is called. I can fix this by appending $pileupsDir and $contigsDir
# when actually opening the files (readdir does not append the file path before the file name)
GetOptions  ("pileups=s"         => \$pileupsDir,
             "contigs=s"         => \$contigsDir,
             "out=s"             => \$outDir,
             "window=i"          => \$windowSize,
             "mindepth=i"        => \$minimumDepth,
             "mincontiglength=i" => \$minimumContigLength,
             "help|man"          => \$help) || die "Couldn't get options with GetOpt::Long: $!\n";

if (!$pileupsDir or !$contigsDir or !$outDir or !$windowSize or !$minimumDepth or !$minimumContigLength or $help) {
    die "Must supply --pileups and --contigs and --out and --window and --mindepth and --mincontiglength.\n";
}
opendir (my $pileupsDH, $pileupsDir);
my @pileupFiles = readdir($pileupsDH);
closedir($pileupsDH);


my %trimLocations;
# %trimLocations will eventually look like:
#      Contig20743 = 1x_HBS_121002 => (35,984),
#                    2x_HBS_19863 => (71, 874)
#     Contig7721 = 1x_HBS_121002 => (121,764),
#                  2x_HBS_19863 => (29, 1201)
#
# That would mean trim everything before bp 35 in 1x_HBS_121002 for contig20743 and after bp 984, etc...

my %pileupHash;
# %pileupHash will look like:
#     Contig20743 = 1x_HBS_121002 => (3,3,3,4,4,4,5,6,5,6,6,6,6,7,8,9,8,8,8,8,7,4,3,2),
#                   2x_HBS_19863  => (4,6,7,8,9,9,10,10,10,11,12,12,11,7,5,4,3,3,3,3,3,1,2,3,2,1,2)

my %targetRange;
# %targetRange will look like:
#     Contig20743 = 1x_HBS_121002 => (1,973),
#                   2x_HBS_19863  => (40, 1294);

my $csv = Text::CSV->new ( { binary => 1, sep_char => "\t", quote_char => "~"} )  # should set binary attribute.
                or die "Cannot use CSV: ".Text::CSV->error_diag ();

foreach my $inFile (@pileupFiles) {
    
    open my $inFileFH, "<", "$inFile" or die "$inFile: $!";
    my $sampleName;
    if ($inFile =~ /(.*)\.mpileup/) {
        $sampleName = $1;
    } else {
        next;
    }
    
    while (my $row = $csv->getline($inFileFH)) {
        #print $row->[3] . " ";
        my $targetName;

        if ($row->[0] =~ /\d+\_(.*)/) {
            $targetName = $1;
        } else {
            die "Target name doesn't start with [digit(s)]_. Exiting because we assume they all do.\n";
        }
        unless (exists $targetRange{$targetName}{$sampleName}) {
            @{$targetRange{$targetName}{$sampleName}} = ($row->[1],0); # When we first come across a target, we initialize it into the hash, and take the first row bp position as the lower bound
        }
        
        if ($row->[1] =~ /(\d+)/) {
            if ($1 > ${$targetRange{$targetName}}{$sampleName}[1]) {
                ${$targetRange{$targetName}{$sampleName}}[1] = $1; # Update the upper range of the target every time we come across a higher bp
            }            
        }
        push @{$pileupHash{$targetName}{$sampleName}}, $row->[3];        
    }
    
    print "Finished processing $sampleName\n";
    $csv->eof or $csv->error_diag();
    close $inFileFH;    
}
#print Dumper(\%pileupHash);
#print Dumper(\%targetRange);

my %realCutSites;
# Now we want to scan through the arrays for each sample within each target in %pileupHash
foreach my $target (sort keys %pileupHash) {
    #print "$target\n";
    foreach my $sample (sort keys %{$pileupHash{$target}}) {
        #print "$sample\n";
        my $leftCut;
        my $rightCut;
        my @meetsMinimum;
        # @{$pileupHash{$target}{$sample}} is the array of read depths for that target for that individual
        # It will contain as many values as there are positions in the mpileup file. So if the first position
        # for a locus in the mpileup file is 42 and the last position is 50, it will contain 9 elements. If we are
        # running a 2bp sliding window, we'd want the last starting position to be 49, so the bottom would read:
        # foreach my $windowStart (0 .. (9-2))
        foreach my $windowStart (0 .. scalar(@{$pileupHash{$target}{$sample}} - $windowSize)) {
            # Example, if we're starting at base pair one, and the window is 20bp long, it should run from 1:20
            my @windowDepths = @{$pileupHash{$target}{$sample}}[$windowStart .. ($windowStart + ($windowSize-1))];
            if ($sample eq '2x_HBS_119093_TGACGTCG' and $target eq 'ACM4') {
                foreach my $element (@windowDepths) {
                    print "$element ";
                }
                print "\n";
            }
            
            
            my $windowAverage = average(\@windowDepths);
            if ($sample eq '2x_HBS_119093_TGACGTCG' and $target eq 'ACM4') {
                #print $windowAverage . "\n";
            }
            
            if ($windowAverage >= $minimumDepth) {
                push(@meetsMinimum, 1);
            } else {
                push(@meetsMinimum, 0);
            }            
        }
        my $longestRun = 0;
        my $longestRunEndBeginningBase = 0;
        my $run = 0;
        my $counter = 1;
        foreach my $pass (@meetsMinimum) {
            if ($pass == 0) {
                $run = 0;
                $counter++;
            } else {
                $run++;
                if ($run > $longestRun) {
                    $longestRun = $run;
                    $longestRunEndBeginningBase = $counter;
                }
                $counter++;
            }
        }
        if ($target eq 'ACM4' and $sample =~ '2x_HBS_119093') {
            print "Longest run for $target and $sample: $longestRun. Elements in array: " . scalar(@meetsMinimum) . "\n";
        }
        
        # We only want contigs that are at least 200bp in length, so don't intialize targets that don't have a run at least that long
        next if ($longestRun <= $minimumContigLength);
        
        # $longestRunEnd now is indexed in terms of the order of bases present in the pileup file. We need to
        # change it so that it starts with bp 1
        $longestRunEndBeginningBase = $longestRunEndBeginningBase + $targetRange{$target}{$sample}[0] - 1;
        my $starting = $longestRunEndBeginningBase - $longestRun;
        $leftCut = $longestRunEndBeginningBase - $longestRun;
        
        # The sliding window stops sliding when it is n base pairs away from the end of the sequence, where n is the
        # length of the window. So the last "good" base is the position of the last window starting position plus the
        # size of the window itself. Let's say the binary pass/fail array looks like 0 0 1 1 1 1 1 0 for an 8bp sequence.
        # 
        
        $rightCut = $longestRunEndBeginningBase + $windowSize - 1;
        
        #print "Interval for $target and $sample: " . $starting . " to " . $longestRunEnd . "\n";
        #   #Since array counting starts at 0 and not 1, we want to add 1 to the $leftCut and $rightCut
        #   #This will translate to removing the base where the transition happened. For instance, if the
        #   #depths are (1,1,2,2,3,3,4) and we have a 2bp window with minimum depth of 2, the remaining bases
        #   #will have depths of (2,3,3,4)
        #   $leftCut = $leftCut + 1;
        #   $rightCut = $rightCut + 1;
        
        # WE USE THE ABOVE SCAN TO TELL US WHICH POSITION IN THE ARRAY DROPS BELOW THE SLIDING WINDOW AVERAGE THRESHOLD.
        # WE THEN MUST TRANSLATE THIS TO THE BEGINNING AND END OF THE ACTUAL CONTIG (SOMETIMES THE PILEUP FILE DOES NOT
        # BEGIN WITH BP1 IN ALL CONTIGS)
        #
        # We made %targetRange above for precisely this purpose.
        #print "Starting of reads for $target and $sample: $targetRange{$target}{$sample}[0]\n";
        # $leftCut is the bp where we want to trim everything to the left of it, so instead of just starting at bp1 automatically,
        # let's start at $leftCut + $targetRange{$target}{$sample}[0];
        $leftCut = $leftCut + $targetRange{$target}{$sample}[0];
        
        # Similarly, we want to cut off everything to the right of $rightCut, but $rightCut is currently in the context of position in
        # the array of depths in the pileup file. Sometimes these depths don't start at bp1
        #$rightCut = $rightCut + $targetRange{$target}{$sample}[0];
        if ($target eq 'ACM4' and $sample =~ '2x_HBS_119093') {
            print "LongestRunEndBeginningBase : $longestRunEndBeginningBase\n";
            print "Total possible target range for $target and $sample: $targetRange{$target}{$sample}[0] to $targetRange{$target}{$sample}[1]\n";
            print "Interval for $target and $sample: " . $leftCut . " to " . $rightCut . "\n";
        }
        $realCutSites{$target}{$sample}[0] = $leftCut;
        $realCutSites{$target}{$sample}[1] = $rightCut;
        #print "Left cut for $target and $sample : $leftCut\nRight cut for $target and $sample : $rightCut\n";
    }
}
# print Dumper(\%realCutSites);
# Contig5853 for sample 2x_HBS_119093_TGACGTCG is a good test case. There is a bad beginning, followed by a passing segment then a failing segment. Then a long run passing segment. It should take the longer (second) passing segment.

# OK, now we've got all of the coordinates that we want. Now we want to apply these
# trimming actions to the contigs. So we're going to make a hash of the actual sequence
# files
opendir (my $contigsDH, $contigsDir);
my @contigFiles = readdir($contigsDH);
closedir($contigsDH);

print "Putting all targets into a hash so we can slice them up\n";
my %contigsHash;
foreach my $contigFile (@contigFiles) {
    if ($contigFile =~ /abyss\_\d+\_(.*)\_RBBHs\.fasta/) {
        my $sample = $1;
        my $seqIn = Bio::SeqIO->new(-file => $contigFile,
                                    -format => 'fasta');
        while (my $seq = $seqIn->next_seq()) {
            if ($seq->display_id() =~ /\d+\_(.*)/) {
                my $target = $1;
                my $sample = 
                $contigsHash{$target}{$sample} = $seq;
            }
        }
    }
}
print "Finished putting all targets into a hash so we can slice them up\n";

#print Dumper(\%contigsHash);


# Now we're going to iterate through all of the cut sites in %realCutSites to physically trim the Bio::Seq objects
# in %contigsHash based on the sliding window trimmer.

my %seqOutHash;

foreach my $target (sort keys %realCutSites) {
    foreach my $sample (sort keys %{$realCutSites{$target}}) {
        print "Target: $target, Sample: $sample\n";
        unless (exists $seqOutHash{$sample}) {
            my $seqOutName = $sample . "_window" . $windowSize . "bp_mindepth" . $minimumDepth . "_trimmed.fasta";
            $seqOutHash{$sample} = Bio::SeqIO->new(-file => ">$seqOutName",
                                                   -format => 'fasta');
        }
        #print "Processing $target for $sample\n";
        my $newSeqName = $contigsHash{$target}{$sample}->display_id() . " " . $contigsHash{$target}{$sample}->desc();
        #print $newSeqName . "\n";
        #$seqOutHash{$sample}->write_seq('seq');
        #my $trimmedSeq = Bio::Seq->new()
        my $trimmedSeq = Bio::Seq->new(-display_id => $contigsHash{$target}{$sample}->display_id(),
                                       -seq => $contigsHash{$target}{$sample}->subseq($realCutSites{$target}{$sample}[0],$realCutSites{$target}{$sample}[1]),
                                       -alphabet => 'dna',
                                       -desc => $contigsHash{$target}{$sample}->desc());
        $seqOutHash{$sample}->write_seq($trimmedSeq);
    }
}




sub average{
        my($data) = @_;
        if (not @$data) {
                die("Empty array\n");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}

















