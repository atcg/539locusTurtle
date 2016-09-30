#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Cwd;

my $help = 0;
my $inDir;
my $outDir;

GetOptions  ("indir=s"      => \$inDir,
             "outdir=s"      => \$outDir,
             "help|man" => \$help) || die "Couldn't get options with GetOpt::Long: $!\n";

if (!$inDir or !$outDir or $help) {
    die "Must supply --indir and --outdir.\n";
}

my $startingDir = getcwd();
opendir(my $inDirFH, $inDir);
my @dirFiles = readdir($inDirFH);
closedir($inDirFH);
chdir($inDir);

unless(-d $outDir) {
    mkdir $outDir;
}

my %translateTable = ("1x_HBS_108574_ACGCAACT" => "sternon_baurii_1x_HBS_108574_ACGCAACT",
"1x_HBS_112648_ATTGACCG" => "Chrysemys_picta_1x_HBS_112648_ATTGACCG",
"1x_HBS_112754_GTTCGGTA" => "Pelusios_subniger_1x_HBS_112754_GTTCGGTA",
"1x_HBS_119014_GGTACGCT" => "Phrynops_hilarii_1x_HBS_119014_GGTACGCT",
"1x_HBS_121002_CCGGTGGA" => "Acanthochelys_spixii_1x_HBS_121002_CCGGTGGA",
"1x_USNM_520644_AGCTCAAC" => "Lissemys_scutata_1x_USNM_520644_AGCTCAAC",
"2x_CAS_196664_TAAGAGCG" => "Pelomedusa_subrufa_2x_CAS_196664_TAAGAGCG",
"2x_FMNH_256544_TGTTACCT" => "Cuora_galbinifrons_2x_FMNH_256544_TGTTACCT",
"2x_HBS_103217_GAGAGTAC" => "Erymnochelys_madagascariensis_2x_HBS_103217_GAGAGTAC",
"2x_HBS_108583_GCAGGCGT" => "Macrochelys_temminckii_2x_HBS_108583_GCAGGCGT",
"2x_HBS_112449_GATGGAGA" => "Manouria_emys_emys_2x_HBS_112449_GATGGAGA",
"2x_HBS_112648_TTCCAAGG" => "Chrysemys_picta_2x_HBS_112648_TTCCAAGG",
"2x_HBS_112687_ACCGCCTA" => "Eretmochelys_imbricata_2x_HBS_112687_ACCGCCTA",
"2x_HBS_112754_TGTAGCCA" => "Pelusios_subniger_2x_HBS_112754_TGTAGCCA",
"2x_HBS_116234_TTCTCCTT" => "Emys_blandingii_2x_HBS_116234_TTCTCCTT",
"2x_HBS_116867_AACCGAGT" => "Dermatemys_mawii_2x_HBS_116867_AACCGAGT",
"2x_HBS_117914_AGTATAGC" => "Staurotypus_triporcatus_2x_HBS_117914_AGTATAGC",
"2x_HBS_117925_CTAACGAT" => "Podocnemis_sextuberculata_2x_HBS_117925_CTAACGAT",
"2x_HBS_119093_TGACGTCG" => "Carretochelys_insculpta_2x_HBS_119093_TGACGTCG",
"2x_HBS_119201_AATACTTC" => "Emys_marmorata_2x_HBS_119201_AATACTTC",
"2x_HBS_121002_GCCATAAC" => "Acanthochelys_spixii_2x_HBS_121002_GCCATAAC",
"2x_HBS_123633_AAGCTTAT" => "Rhinoclemmys_punctularia_2x_HBS_123633_AAGCTTAT",
"2x_HBS_124036_TCCAGGAA" => "Homopus_areolatus_2x_HBS_124036_TCCAGGAA",
"2x_HBS_16255_TTGTTGGC" => "Platysternon_megacephalum_2x_HBS_16255_TTGTTGGC",
"2x_HBS_23551_GGCGAAGG" => "Chelydra_serpentina_2x_HBS_23551_GGCGAAGG",
"2x_MVZ_149847_CCAGCACC" => "Dermochelys_coriacea_2x_MVZ_149847_CCAGCACC",
"2x_USNM_520644_CCAAGACT" => "Lissemys_scutata_2x_USNM_520644_CCAAGACT",
"2x_USNM_572079_CGCTTGGA" => "Dogania_subplana_2x_USNM_572079_CGCTTGGA",
"4x_HBS_112648_GGTAGTGT" => "Chrysemys_picta_4x_HBS_112648_GGTAGTGT",
"4x_HBS_112754_TGCGAACT" => "Pelusios_subniger_4x_HBS_112754_TGCGAACT",
"4x_HBS_121002_GAGTCTCT" => "Acanthochelys_spixii_4x_HBS_121002_GAGTCTCT",
"4x_USNM_520644_ATTGCGTG" => "Lissemys_scutata_4x_USNM_520644_ATTGCGTG",
"8x_HBS_100002_ACTGAGGT" => "Chrysemys_picta_8x_HBS_100002_ACTGAGGT",
"8x_HBS_100088_CACTGACA" => "Chrysemys_picta_8x_HBS_100088_CACTGACA",
"8x_HBS_100115_CAGTCCAA" => "Chrysemys_picta_8x_HBS_100115_CAGTCCAA",
"8x_HBS_100121_TCGACATC" => "Chrysemys_picta_8x_HBS_100121_TCGACATC",
"8x_HBS_112648_AGAATGCC" => "Chrysemys_picta_8x_HBS_112648_AGAATGCC",
"8x_HBS_112754_ATGCACGA" => "Pelusios_subniger_8x_HBS_112754_ATGCACGA",
"8x_HBS_121002_TACGGTTG" => "Acanthochelys_spixii_8x_HBS_121002_TACGGTTG",
"8x_HBS_23676_AAGGCGTT" => "Chrysemys_picta_8x_HBS_23676_AAGGCGTT",
"8x_HBS_26917_GTCCACAT" => "Chrysemys_picta_8x_HBS_26917_GTCCACAT",
"8x_HBS_27136_TACACGCT" => "Chrysemys_picta_8x_HBS_27136_TACACGCT",
"8x_HBS_27180_ATCTGTCC" => "Chrysemys_picta_8x_HBS_27180_ATCTGTCC",
"8x_HBS_27254_GAGGACTT" => "Chrysemys_picta_8x_HBS_27254_GAGGACTT",
"8x_HBS_27287_GTTGTTCG" => "Chrysemys_picta_8x_HBS_27287_GTTGTTCG",
"8x_HBS_27379_AGACCGTA" => "Chrysemys_picta_8x_HBS_27379_AGACCGTA",
"8x_HBS_27682_ATGGCGAA" => "Chrysemys_picta_8x_HBS_27682_ATGGCGAA",
"8x_HBS_27748_GAATCCGA" => "Chrysemys_picta_8x_HBS_27748_GAATCCGA",
"8x_HBS_27885_ACCTTCTC" => "Chrysemys_picta_8x_HBS_27885_ACCTTCTC",
"8x_HBS_28003_GAAGGAAG" => "Chrysemys_picta_8x_HBS_28003_GAAGGAAG",
"8x_HBS_28659_TCGAGTGA" => "Chrysemys_picta_8x_HBS_28659_TCGAGTGA",
"8x_HBS_31496_CACCTGTT" => "Chrysemys_picta_8x_HBS_31496_CACCTGTT",
"8x_HBS_31527_GTCACTGT" => "Chrysemys_picta_8x_HBS_31527_GTCACTGT",
"8x_HBS_33307_CCGTAAGA" => "Chrysemys_picta_8x_HBS_33307_CCGTAAGA",
"8x_HBS_35486_GTCTTGCA" => "Chrysemys_picta_8x_HBS_35486_GTCTTGCA",
"8x_USNM_520644_CAAGTGCA" => "Lissemys_scutata_8x_USNM_520644_CAAGTGCA");


foreach my $file(@dirFiles) {
    # Translate all instances of every sample number with sampleName_sample_number
    if ($file =~ /(.*)\.aligned/) {
        open(my $FH, "<", $file) or die "Couldn't open $file for reading: $!\n";
        my $outFileName = $outDir . "/" . $1 . "_rn.aligned";
        open(my $outFH, ">", $outFileName) or die "Couldn't open $outFileName for writing: $!\n";
        while (my $line = <$FH>) {
            if ($line =~ /^\>/) {
                $line =~ s/(1x_HBS_108574_ACGCAACT|1x_HBS_112648_ATTGACCG|1x_HBS_112754_GTTCGGTA|1x_HBS_119014_GGTACGCT|1x_HBS_121002_CCGGTGGA|1x_USNM_520644_AGCTCAAC|2x_CAS_196664_TAAGAGCG|2x_FMNH_256544_TGTTACCT|2x_HBS_103217_GAGAGTAC|2x_HBS_108583_GCAGGCGT|2x_HBS_112449_GATGGAGA|2x_HBS_112648_TTCCAAGG|2x_HBS_112687_ACCGCCTA|2x_HBS_112754_TGTAGCCA|2x_HBS_116234_TTCTCCTT|2x_HBS_116867_AACCGAGT|2x_HBS_117914_AGTATAGC|2x_HBS_117925_CTAACGAT|2x_HBS_119093_TGACGTCG|2x_HBS_119201_AATACTTC|2x_HBS_121002_GCCATAAC|2x_HBS_123633_AAGCTTAT|2x_HBS_124036_TCCAGGAA|2x_HBS_16255_TTGTTGGC|2x_HBS_23551_GGCGAAGG|2x_MVZ_149847_CCAGCACC|2x_USNM_520644_CCAAGACT|2x_USNM_572079_CGCTTGGA|4x_HBS_112648_GGTAGTGT|4x_HBS_112754_TGCGAACT|4x_HBS_121002_GAGTCTCT|4x_USNM_520644_ATTGCGTG|8x_HBS_100002_ACTGAGGT|8x_HBS_100088_CACTGACA|8x_HBS_100115_CAGTCCAA|8x_HBS_100121_TCGACATC|8x_HBS_112648_AGAATGCC|8x_HBS_112754_ATGCACGA|8x_HBS_121002_TACGGTTG|8x_HBS_23676_AAGGCGTT|8x_HBS_26917_GTCCACAT|8x_HBS_27136_TACACGCT|8x_HBS_27180_ATCTGTCC|8x_HBS_27254_GAGGACTT|8x_HBS_27287_GTTGTTCG|8x_HBS_27379_AGACCGTA|8x_HBS_27682_ATGGCGAA|8x_HBS_27748_GAATCCGA|8x_HBS_27885_ACCTTCTC|8x_HBS_28003_GAAGGAAG|8x_HBS_28659_TCGAGTGA|8x_HBS_31496_CACCTGTT|8x_HBS_31527_GTCACTGT|8x_HBS_33307_CCGTAAGA|8x_HBS_35486_GTCTTGCA|8x_USNM_520644_CAAGTGCA)/$translateTable{$1} || $1/e;
                $line =~ s/(.*)\_\S+\s.*$/$1/;
            }
        print $outFH $line;            
        }
    }
}




































