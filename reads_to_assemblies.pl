#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Cwd;
use Parallel::ForkManager;

my $help = 0;
my $outDir;
my $logFile;
my $readsDir;
my $adaptersDir;
my $kmer = 51;
my $threadsMax = 1;

GetOptions  ("out=s"    => \$outDir,
             "reads=s"  => \$readsDir,
             "adapters=s" => \$adaptersDir,
             "log=s"    => \$logFile,
             "threads=i"=> \$threadsMax,
             "kmer=i"   => \$kmer,
             "help|man" => \$help) || pod2usage(2);

if (!$outDir or !$logFile or !$readsDir or $help) {
    pod2usage(-exitval => 0, -verbose => 2, -noperldoc => 1);
}

my $startingDir = getcwd();
unless (-d $outDir) {
    mkdir $outDir or die "Couldn't make output directory $outDir: $!\n";
}


open(my $logFH, ">", $logFile) or die "Couldn't open log file $logFile for writing: $!\n";

# First gather up all the file names of the files in the specified reads directory
my %readFilesHash;
my %sampleNamesHash;
chdir($readsDir);
opendir my $readsDirFH, "./";
my @readsFiles = readdir $readsDirFH or die "Couldn't readdir $readsDirFH: $!\n";
closedir $readsDirFH;
chdir $startingDir;
my $trimmomaticDir = $outDir . "/trimmomatic";
unless (-d $trimmomaticDir) {
    mkdir $trimmomaticDir or die "Couldn't make directory $trimmomaticDir: $!\n";
}

# Now do sequence QC using Trimmomatic
{
my @trimmomaticCommands;
print $logFH "Generating trimmomatic commands on all read files.\n";
foreach my $file (@readsFiles) {
    if ($file =~ /HSGM001_(.*)_(.*)_L002_R1_(00[12]).fastq.gz/) { # Only looking for the R1s here, because we only want one trimmomatic run per library (or in the case of too many reads for a single fastq.gz file, do a trimmomatic run for each pair of read files)
        # $1 = sample name
        # $2 = index sequence
        # $3 = reads group 1 or 2 (the only ones that will get 2 are the ones that had too many reads to fit in a single fastq file)
        print $logFH "--------------------------------------------------\n";
        print $logFH "Sequence file found: HSGM001_$1\_$2\_L002_R1\_$3\.fastq\.gz\n";
        my $R1File = $readsDir . "HSGM001_$1\_$2\_L002_R1\_$3\.fastq\.gz";
        print $logFH "R1 file: $readsDir" . "HSGM001_$1\_$2\_L002_R1\_$3\.fastq\.gz\n";
        my $R2File = $readsDir . "HSGM001_$1\_$2\_L002_R2\_$3\.fastq\.gz";
        print $logFH "R2 file: $readsDir" . "HSGM001_$1\_$2\_L002_R2\_$3\.fastq\.gz\n";
        my $readGroupName = $1 . "_" . $2 . "_" . "$3";
        
        my $adaptersFile = $adaptersDir . "individuals_$1";
        print $logFH "Adapters file to be used for the sequence group: " . $adaptersFile . "\n";
        
        my $R1OutFilePaired = "$trimmomaticDir/$1\_$2\_L002\_R1\_$3\_paired_trimmed.fastq.gz";
        print $logFH "R1 paired trimmomatic output file: $trimmomaticDir/$1\_$2\_L002\_R1\_$3\_paired_trimmed.fastq.gz\n";
        my $R1OutFileSingles = "$trimmomaticDir/$1\_$2\_L002\_R1\_$3\_singles_trimmed.fastq.gz";
        print $logFH "R1 singles trimmomatic output file: $trimmomaticDir/$1\_$2\_L002\_R1\_$3\_singles_trimmed.fastq.gz\n";
        my $R2OutFilePaired = "$trimmomaticDir/$1\_$2\_L002\_R2\_$3\_paired_trimmed.fastq.gz";
        print $logFH "R2 paired trimmomatic output file: $trimmomaticDir/$1\_$2\_L002\_R2\_$3\_paired_trimmed.fastq.gz\n";
        my $R2OutFileSingles = "$trimmomaticDir/$1\_$2\_L002\_R2\_$3\_singles_trimmed.fastq.gz";
        print $logFH "R2 singles trimmomatic output file: $trimmomaticDir/$1\_$2\_L002\_R2\_$3\_singles_trimmed.fastq.gz\n";
        $sampleNamesHash{$readGroupName}{'R1_paired_trimmed'} = $R1OutFilePaired;
        $sampleNamesHash{$readGroupName}{'R1_singles_trimmed'} = $R1OutFileSingles;
        $sampleNamesHash{$readGroupName}{'R2_paired_trimmed'} = $R2OutFilePaired;
        $sampleNamesHash{$readGroupName}{'R2_singles_trimmed'} = $R2OutFileSingles;        
        push (@trimmomaticCommands, "java -Xmx8G -jar ~/bin/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 2 -phred33 $R1File $R2File $R1OutFilePaired $R1OutFileSingles $R2OutFilePaired $R2OutFileSingles ILLUMINACLIP:$adaptersFile:2:30:10 LEADING:5 TRAILING:15 SLIDINGWINDOW:4:20 MINLEN:40");
    }
}
print $logFH "--------------------------------------------------\n";
print $logFH "Finished generating trimmomatic commands on all read files.\n\n\n";

print $logFH "Running all trimmomatic commands\n";
my $counter = 0;
my $forkManager = new Parallel::ForkManager($threadsMax);
foreach my $trimCommand (@trimmomaticCommands) {
    $counter++;
    print $logFH "--------------------------------------------------\n";
    print $logFH "Trimmomatic command $counter: \n\t"; # Indent the next line to make it easier to find the commands in the text
    print $logFH $trimCommand . "\n";
    sleep 10;
    print "\n";
    $forkManager->start and next;
    print "\n";
    system("$trimCommand");
    print "Finished running the following:\n\t$trimCommand\n\n";
    $forkManager->finish;
}
$forkManager->wait_all_children;
print $logFH "--------------------------------------------------\n";
print $logFH "Finished running all trimmomatic commands\n";
print $logFH "--------------------------------------------------\n\n";
}

# Now we can iterate through the %sampleNamesHash that we populated earlier to do the following:
#  1. Find read groups that have both 001 and 002, and cat the 002s onto the end of the 001s
#  1b. Delete the 002 files
#  2. Run fastq-join on the output R1 and R2 paired files
#  3. Merge the R1 and R2 singles files together into one singles file, and add the fastq-join joined reads with the "true" single reads
my @filesToZip;
my @assemblyFiles;

my $forkManager2 = new Parallel::ForkManager(2); # We want to be able to run in parallel two biological samples for the gunzip, fastq-join, gzip, and assembly stages

foreach my $readGroup (sort keys %sampleNamesHash) {
    sleep 10;
    print $logFH "\n\nStarting to process $readGroup through the fastq-join and assembly stages\n\n";
    $forkManager2->start and next;
    # $readGroup should be something like "1x_HBS_119014_GGTACGCT_001"
    if ($readGroup =~ /(.*)_001$/) {
        my $readGroup2 = $1 . "_002";
        if (exists $sampleNamesHash{$readGroup2}) {
            print $logFH "\n--------------------------------------------------\n";
            print $logFH "\nAppending all reads from 002 readgroups to the end of the corresponding 001s for $1 (and deleting the 002 files)\n";
            print $logFH "\nConcatenation commands:\n\t";
            print $logFH "\ncat $sampleNamesHash{$readGroup2}{'R1_paired_trimmed'} >> $sampleNamesHash{$readGroup}{'R1_paired_trimmed'}\n\t";
            print $logFH "\ncat $sampleNamesHash{$readGroup2}{'R1_singles_trimmed'} >> $sampleNamesHash{$readGroup}{'R1_singles_trimmed'}\n\t";
            print $logFH "\ncat $sampleNamesHash{$readGroup2}{'R2_paired_trimmed'} >> $sampleNamesHash{$readGroup}{'R2_paired_trimmed'}\n\t";
            print $logFH "\ncat $sampleNamesHash{$readGroup2}{'R2_singles_trimmed'} >> $sampleNamesHash{$readGroup}{'R2_singles_trimmed'}\n\n";
            system("cat $sampleNamesHash{$readGroup2}{'R1_paired_trimmed'} >> $sampleNamesHash{$readGroup}{'R1_paired_trimmed'}");
            system("cat $sampleNamesHash{$readGroup2}{'R1_singles_trimmed'} >> $sampleNamesHash{$readGroup}{'R1_singles_trimmed'}");
            system("cat $sampleNamesHash{$readGroup2}{'R2_paired_trimmed'} >> $sampleNamesHash{$readGroup}{'R2_paired_trimmed'}");
            system("cat $sampleNamesHash{$readGroup2}{'R2_singles_trimmed'} >> $sampleNamesHash{$readGroup}{'R2_singles_trimmed'}");
            # Now that we've appended the 002 readgroups to the end of the 001s, we can delete the 002 files:
            unlink($sampleNamesHash{$readGroup2}{'R1_paired_trimmed'},
                   $sampleNamesHash{$readGroup2}{'R1_singles_trimmed'},
                   $sampleNamesHash{$readGroup2}{'R2_paired_trimmed'},
                   $sampleNamesHash{$readGroup2}{'R2_singles_trimmed'});
            print $logFH "\nFinished appending all reads from 002 readgroups to the end of the corresponding 001s for $1 (and deleting the 002 files)\n\n";
            print $logFH "\n--------------------------------------------------\n\n";

        }
        print $logFH "\n--------------------------------------------------\n";
        print $logFH "\nUnzipping files for $1 in the trimmomatic folder for input into fastq-join\n";
        system("gunzip $sampleNamesHash{$readGroup}{'R1_paired_trimmed'}");
        system("gunzip $sampleNamesHash{$readGroup}{'R1_singles_trimmed'}");
        system("gunzip $sampleNamesHash{$readGroup}{'R2_paired_trimmed'}");
        system("gunzip $sampleNamesHash{$readGroup}{'R2_singles_trimmed'}");
        print $logFH "\nFinished unzipping files for $1 in the trimmomatic folder for input into fastq-join\n\n";
        print $logFH "\n--------------------------------------------------\n\n";


        my $fastqJoinOutDir = $outDir . "/fastq-join/";
        unless (-d $fastqJoinOutDir) {
            mkdir $fastqJoinOutDir;
        }
        chdir $fastqJoinOutDir;
        print $logFH "\n--------------------------------------------------\n";
        print $logFH "\nRunning fastq-join on $1\n";
        print $logFH "\nfastq-join command used:\n\t";
        my $fastqJoinOutputName = $1 . ".%.fastq";
        print $logFH "\nfastq-join -v ' ' -m 8 ../trimmomatic/$1\_L002\_R1\_001_paired_trimmed.fastq ../trimmomatic/$1\_L002\_R2\_001_paired_trimmed.fastq -o $fastqJoinOutputName";
        system("fastq-join -v ' ' -m 8 ../trimmomatic/$1\_L002\_R1\_001_paired_trimmed.fastq ../trimmomatic/$1\_L002\_R2\_001_paired_trimmed.fastq -o $fastqJoinOutputName");
        print $logFH "\nFinished running fastq-join on $1\n";
        print $logFH "\n--------------------------------------------------\n\n";
        
        my $joinedFile = $1 . ".join.fastq";
        my $newR1file = $1 . ".un1.fastq";
        my $newR2file = $1 . ".un2.fastq";

        my $singlesFile1 = "../trimmomatic/$1\_L002_R1_001_singles_trimmed.fastq";
        my $singlesFile2 = "../trimmomatic/$1\_L002_R2_001_singles_trimmed.fastq";
        my $joinedAndSinglesFile = $1 . "_joined_and_both_singles.fastq";
        
        print $logFH "\n--------------------------------------------------\n";
        print $logFH "\nCombining fastq-joined reads with the singletons Trimmomatic made from R1 and R2.\n";
        print $logFH "\nConcatenation command:\n\t";
        print $logFH "\ncat $joinedFile $singlesFile1 $singlesFile2 > $joinedAndSinglesFile\n";
        system("cat $joinedFile $singlesFile1 $singlesFile2 > $joinedAndSinglesFile");
        print $logFH "\nRemoving the joined and singles files\n";
        print $logFH "\nFile removal command:\n\t";
        print $logFH "\nunlink($joinedFile)\n";
        unlink($joinedFile);
        print $logFH "\nFinished combining joined reads with singletons and removing the consitutuent files.\n";
        print $logFH "\n--------------------------------------------------\n\n";
        chdir $startingDir;
        
        my $assemblyDir = $outDir . "/assembly/";
        unless (-d $assemblyDir) {
            mkdir $assemblyDir;
        }
        my $sampleAssemblyDir = $assemblyDir . "k$kmer" . "_$1/"; # $1 should be something like "1x_HBS_119014_GGTACGCT" from the regex at the top of the loop
        unless (-d $sampleAssemblyDir) {
            mkdir $sampleAssemblyDir;
        }
        my $abyssOutputPrefix = "abyss_$kmer\_$1";
        my $abyssLogFile = $abyssOutputPrefix . ".log";
        
        # Run the assembly
        print $logFH "--------------------------------------------------\n";
        print $logFH "\nRunning Abyss using ../../fastq-join/$newR1file, ../../fastq-join/$newR2file, and ../../fastq-join/$joinedAndSinglesFile\n";
        print $logFH "\nAbyss command (after changing the working directory to $sampleAssemblyDir):\n\t";
        print $logFH "\nabyss-pe k=$kmer name=$abyssOutputPrefix in='../../fastq-join/$newR1file ../../fastq-join/$newR2file' se='../../fastq-join/$joinedAndSinglesFile' > $abyssLogFile 2>&1\n";
        chdir $sampleAssemblyDir;
        system("abyss-pe k=$kmer name=$abyssOutputPrefix in='../../fastq-join/$newR1file ../../fastq-join/$newR2file' se='../../fastq-join/$joinedAndSinglesFile' > $abyssLogFile 2>&1");
        chdir $startingDir;
        print $logFH "\nFinished running Abyss using $newR1file, $newR2file, and $joinedAndSinglesFile\n\n";
        print $logFH "\n--------------------------------------------------\n\n";

        
        my $originalR1pairedTrimmed = $1 if ($sampleNamesHash{$readGroup}{'R1_paired_trimmed'} =~ /(.*)\.gz/);
        my $originalR1singlesTrimmed = $1 if ($sampleNamesHash{$readGroup}{'R1_singles_trimmed'} =~ /(.*)\.gz/);
        my $originalR2pairedTrimmed = $1 if ($sampleNamesHash{$readGroup}{'R2_paired_trimmed'} =~ /(.*)\.gz/);
        my $originalR2singlesTrimmed = $1 if ($sampleNamesHash{$readGroup}{'R2_singles_trimmed'} =~ /(.*)\.gz/);        
        
        push (@filesToZip, $fastqJoinOutDir.$joinedAndSinglesFile, $fastqJoinOutDir.$newR1file, $fastqJoinOutDir.$newR2file, $trimmomaticDir.$originalR1pairedTrimmed, $trimmomaticDir.$originalR1singlesTrimmed, $trimmomaticDir.$originalR2pairedTrimmed, $trimmomaticDir.$originalR2singlesTrimmed);
        my $assembly = $sampleAssemblyDir . $abyssOutputPrefix . "-contigs.fa";
        push (@assemblyFiles, $assembly);
    }
    $forkManager2->finish;

    print $logFH "--------------------------------------------------\n";
    print $logFH "Finished running all fastq-join and assembly commands\n";
    print $logFH "--------------------------------------------------\n\n";
}
$forkManager2->wait_all_children;

# GZip files to save some file space
foreach my $file (@filesToZip) {
    print $logFH "gzipping $file\n";
    system("gzip $file");
}
print $logFH "The following is a list of all the assembly files:\n";
foreach my $file (@assemblyFiles) {
    print $logFH "$file\n";
}

        # # GZip operations now to save some file space
        # print $logFH "--------------------------------------------------\n";
        # print $logFH "\nZipping the files in the fastq-join folder back up to save space\n";
        # print $logFH "\nZip command:\n\t";
        # print $logFH "\ngzip $joinedAndSinglesFile\n\t";
        # print $logFH "\ngzip $newR1file\n\t";
        # print $logFH "\ngzip $newR2file\n";
        # system("gzip $joinedAndSinglesFile");
        # system("gzip $newR1file");
        # system("gzip $newR2file");
        # print $logFH "\nFinished zipping the files back up to save space\n\n";
        # print $logFH "--------------------------------------------------\n";
        # 
        # print $logFH "--------------------------------------------------\n";
        # print $logFH "Zipping the files in the trimmomatic folder back up to save space\n";
        # 
        # print $logFH "Zip command:\n\t";
        # print $logFH "gzip $originalR1pairedTrimmed\n\t";
        # print $logFH "gzip $originalR1singlesTrimmed\n\t";
        # print $logFH "gzip $originalR2pairedTrimmed\n\t";
        # print $logFH "gzip $originalR2singlesTrimmed\n";
        # system("gzip $originalR1pairedTrimmed");
        # system("gzip $originalR1singlesTrimmed");
        # system("gzip $originalR2pairedTrimmed");
        # system("gzip $originalR2singlesTrimmed");
        # print $logFH "Finished zipping the files in the trimmomatic folder back up to save space\n\n";
        # print $logFH "--------------------------------------------------\n";








#Documentation
__END__

=head1 NAME

reads_to_assemblies.pl ##CHANGE

=head1 SYNOPSIS 

perl reads_to_assemblies.pl --out <directory_name> --log <file> --threads <integer> --reads <directory> --kmer <integer>

 Options:
   -in=s            Input filename
   -out=s           Name of directory where all output will be saved
   -reads=s         Name of directory where all raw reads are
   -log=s           Log filename
   -threads=i       Maximum number of threads (default 1)
   -kmer=i          kmer value to use for Abyss assemblies
   -help|man        Prints out documentation


=head1 DESCRIPTION

This program takes raw read files in fastq.gz format and does sequence quality
trimming and adapter contamination removal. It then merges overlapping paired
end reads, and combines the single orphaned reads from the qc process with the
joined reads to create the "single-end" library for each sample.

It then does a de novo assembly using Abyss for each biological sample.

=cut









