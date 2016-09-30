#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;
use Bio::SearchIO;
use Cwd;

my $help = 0;
my $allAssemblyDirs;
my $targets;

GetOptions  ("targets=s"    => \$targets,
             "assemblies=s" => \$allAssemblyDirs,
             "help|man"     => \$help) || die "Couldn't get options with GetOpt::Long: $!\n";

if (!$targets or !$allAssemblyDirs or $help) {
    die "Must supply --in and --out and --targets and --assemblies.\n";
}

# This script goes through all the assembly folders in a directory and does the following:
#  1. make a blast database from the final assembly and keep it in that folder with the contigs.fa file
#  2. blast the targets against this blast database
#  3. blast the assembly against the target blast database
#  4. finds reciprocal best blast hits for all targets and then:
#    4a. rename the assembly contig to represent which target it is homologous to
#    4b. write the matching renamed assembled contigs to an output file which contains all the putative target regions for that library

# Get the names of all the read files in the rawReads directory
opendir my $assembliesDirFH, $allAssemblyDirs or die "cannot open dir $allAssemblyDirs: $!";
my @readDirAssemblyDirs= readdir $assembliesDirFH; # This should contain a bunch of assembly directories with names like "8x_HBS_33307_CCGTAAGA_51"
closedir $assembliesDirFH;

# Make sure we have a blast database of the targets ready:
system("makeblastdb -in $targets -dbtype nucl -out targets");
my %targetHash;
my $targetIn = Bio::SeqIO->new(-file => $targets,
                            -format => 'fasta');
while (my $seq = $targetIn->next_seq()) {
    $targetHash{$seq->display_id()} = $seq;
}

# Now we actually need to go into each of these directories to do our analyses
my $startingDir = getcwd;
foreach my $assemblyDir (@readDirAssemblyDirs) {
    chdir "$allAssemblyDirs/$assemblyDir";
    opendir my $individualAssemblyDirFH, "./";
    my @individualAssemblyFiles = readdir $individualAssemblyDirFH or die "Couldn't readdir $individualAssemblyDirFH: $!\n";
    closedir $individualAssemblyDirFH;
    foreach my $file (@individualAssemblyFiles) {
        if ($file =~ /(.*)\-contigs\.fa/) {
            # Make the blast database:
            system("makeblastdb -in $file -dbtype nucl -out assembly");
            
            # Blast the targets against the assembly database we just created (only show the best hit)
            system("blastn -db assembly -query $startingDir/$targets -outfmt 6 -max_target_seqs 1 -out targets_bl2_assembly.blast");
            
            # Blast the assembly against the targets database
            system("blastn -db $startingDir/targets -query $file -outfmt 6 -max_target_seqs 1 -out assembly_bl2_targets.blast");
            
            # Now find the reciprocal best blast hits
            my %assemblyHash;
            my $assemblyIn = Bio::SeqIO->new(-file => $file,
                                            -format => 'fasta');
            while (my $seq = $assemblyIn->next_seq()) {
                $assemblyHash{$seq->display_id()} = $seq;
            }

            my $assembly2targetsBlast = Bio::SearchIO->new(-file => 'assembly_bl2_targets.blast',
                                                        -format => 'blasttable');
            
            my %assem2targBlastResults;
            while (my $result = $assembly2targetsBlast->next_result()) {
                $assem2targBlastResults{$result->query_name()} = ($result->next_hit())->name(); # This finds the name of the first hit and makes that the value for the key that is the name of the query
            }
            
            my $targ2assembBlast = Bio::SearchIO->new(-file => 'targets_bl2_assembly.blast',
                                                    -format => 'blasttable');
            my %targ2assemBlastResults;
            while (my $result = $targ2assembBlast->next_result()) {
                # print $result->query_name() . "\n";
                $targ2assemBlastResults{$result->query_name()} = ($result->next_hit())->name();  # This finds the name of the first hit and makes that the value for the key that is the name of the query
            }
            my $RBBfile = $1 . "_RBBHs.fasta";
            my $seqOut = Bio::SeqIO->new(-file=>">$RBBfile",
                                        -format => 'fasta');
            # That should be good the the names of all the hits for both blast reports
            my $counter = 0;
            foreach my $seqName (sort keys %targ2assemBlastResults) {
                # $seqName should be something like "Contig63"
                # $targ2assemBlastResults{$seqName} will be something like "6" (just a single number)
                # We want to see if $assem2targBlastResults{[the name of the contig in the assembly that was best match for target]} is the same as the name of the contig in the targets that matched the assembly
                if ($assem2targBlastResults{$targ2assemBlastResults{$seqName}} eq $seqName) {
                    #print "Match!\n";
                    $counter++;
                    #print "seqname: $seqName\n";
                    #print "Name: $targ2assemBlastResults{$seqName}\n";
                    my $seqDisplayName = $assemblyHash{$targ2assemBlastResults{$seqName}}->display_id();
                    my $newSeqName = $seqDisplayName . "_$seqName";
                    $assemblyHash{$targ2assemBlastResults{$seqName}}->display_id($newSeqName);
                    $seqOut->write_seq($assemblyHash{$targ2assemBlastResults{$seqName}});
                }
            }
            print "$counter total RBB hits found for assembly $file\n";
        }
    }    
    chdir $startingDir;    
}











