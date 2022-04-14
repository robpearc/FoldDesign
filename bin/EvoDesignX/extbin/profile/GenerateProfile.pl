#!/usr/bin/perl

###############################################################################
# GenerateProfile.pl
#                              
# Prepares profile file of a given target protein for subsequent design process
################################################################################

use strict;
use warnings;

if(-e "prf.txt"){
  printf "PSSM file prf.txt exists, do not run GenerateProfile.pl\n";
  exit;
}

######################################
use Getopt::Long qw(GetOptions);
    Getopt::Long::Configure qw(gnu_getopt);
use Data::Dumper;

# PROGRAM PATH
use File::Basename;
use Cwd 'abs_path';
my $bin_path = dirname(abs_path(__FILE__));

#################################################
use Benchmark qw(:hireswallclock);
my $starttime = Benchmark->new;
my $finishtime;
my $timespent;

######################################
# INIT.VARS
my $hostName = `hostname`;
my $timeStart = localtime();

######################################
# INPUT-OPTIONS-DEFAULTS
my $cutoff=0.7;

######################################
# PDB-DATABASE
my $pdbset="/nfs/amino-library/PDB";
if($hostName =~ /comet/){
  $pdbset="/oasis/projects/nsf/mia181/zhanglab/library/PDB";
}

######################################
# PARSE INPUT OPTIONS with GETOPT
GetOptions(
    'pdblib|l=s' => \$pdbset,
    'tmcutoff|t=f' => \$cutoff,
) or die "Illegal arguments are used for GenerateProfile.pl";

######################################
# PARSE INPUT ARGUMENTS
my $querypdb = $ARGV[0];
my $TargetID = $querypdb;
$TargetID =~ s/\.pdb//g;

print "#starting GenerateProfile.pl at $timeStart\n";
print "host: $hostName";

################################################################################
####  GENERATE SEQUENCES PROFILE 
################################################################################
print "structure alignment will be done against PDB structures available in: $pdbset\n";
my @decoys = glob("$pdbset/*.pdb");

# Total REDUNDANT PDBs in database 
my $ntmp = scalar @decoys;
printf "number of redundant pdb templates:  %d\n", $ntmp;
$ntmp = 0; #Restart Protein Counter

print "Starting the structure alignment process, with TM-Score cut-off: $cutoff \n";
print "TM-Scores with each template are given in first column...\n";

open TMFILE,  ">msa_from_tmalign.txt" or die "Error in creating msa_from_tmalign.txt\n";
open TMFASTA, ">msa_from_tmalign.fas" or die "Error in creating msa_from_tmalign.fas\n";

for(my $i = 0; $i<@decoys; ++$i){
  chomp(my $dec = $decoys[$i]);
  my @pdbpath=split(/\//, $dec);
  my $PDBID = $pdbpath[$#pdbpath];
  $PDBID =~s/\.pdb//g;
  my $idLength = length $PDBID;
  
  #Run TMAlign, for full chains only; Fragmented domains are skipped;
  my @tmo;
  if($idLength<6){
    $ntmp++; #Count real number of PDBs searched
    @tmo = `$bin_path/TMalign $dec $querypdb `;
  }

  my $tmScore;
  my $homology = 0;
  my $seqid;
  for(my $j=0; $j<@tmo; $j++) {
    if($tmo[$j] =~ /TM-score=/) {
      # Extract TM-Scores from TM-Align output
      my @dat1=split(/,/, $tmo[$j]);
      my @dat2=split(/=/, $dat1[2]);
      # Printing TM-Scores for each decoy
      $tmScore = $dat2[1];
      print "$dat2[1]\t$dec\n";
      `echo "$dat2[1]\t$PDBID" >> tmscore.txt`;

      # print the sequence identity
      my @dat3=split(/=/, $dat1[3]);
      chomp($seqid = $dat3[1]);
      `echo "$seqid\t$PDBID" >> seqid.txt`;

      if($dat2[1]>$cutoff){
        $homology=1;
      }
      else{
        $homology=0;
      }
      last;
    }
  }
 
  #Skip-NonHomologues
  if($homology==0) {next;}
  `echo "$seqid\t$PDBID" >> msa.seqid`;

  my ($dali,$tali);
  #Process-Homologues
  for(my $j=0; $j<@tmo; $j++) {
    if($tmo[$j] =~ /denotes residue/) {
      #Template/Decoy-Sequence
      chomp($dali = $tmo[$j+1]);
      #Target-Sequence
      chomp($tali = $tmo[$j+3]);
      last;
    }
  }

  #my ($lena,$ti);
  #my @ali;
  ##Save-Aligned-Decoy
  #$lena = length $tali;
  #for(my $j=0; $j<$lena; $j++) {
  #  if(substr($tali, $j, 1) ne '-') {
  #    $ali[$i][$ti++] = substr($dali, $j, 1);
  #  }
  #}

  my ($lena,$ti);
  my @ali;
  #Save-Aligned-Decoy
  $lena = length $tali;
  my $sameaa = 0;
  my $nongapchar = 0;
  for(my $j=0; $j<$lena; $j++) {
    my $char_t = substr($tali, $j, 1);
    if($char_t ne '-') {
      my $char_d = substr($dali, $j, 1);
      $ali[$i][$ti++] = $char_d;
      if($char_d ne '-'){
        $nongapchar++;
      }
      if($char_t eq $char_d){
        $sameaa++;
      }
    }
  }

  ##set filters to get high-quality template structural homologs
  ##print "Sequence identity (comparted to query): %f\n", $sameaa/$ti;
  ##1) reject sequence with too many gaps
  if($nongapchar/$ti < 0.7){
    next;
  }
  ##2) reject sequences nearly identical to the query sequence
  #elsif($nongapchar/$ti>0.99 && $sameaa/$nongapchar>0.99){
  #  next;
  #}

  #Print-Decoy-Sequences 
  print TMFASTA "> $PDBID;tmSCO=$tmScore \n";
  for(my $j=0; $j<$ti; $j++) {
    print TMFILE $ali[$i][$j];
    print TMFASTA $ali[$i][$j];
  }
  print TMFILE "\n";
  print TMFASTA "\n";
}
    
close TMFILE;
close TMFASTA;

# Total NON-REDUNDANT PDBs in database 
print "\n";
printf "Number of Non-Redundat PDB templates searched: %d\n", $ntmp;

# Sort from Highest TM-Align Score 
`sort -k1 -n -r tmscore.txt > sorted_tmscore.txt`;
`head -100 sorted_tmscore.txt >./top_tmscore.txt`;

#Sort from Highest sequence identity
`sort -k1 -n -r seqid.txt > sorted_seqid.txt`;
`head -100 sorted_seqid.txt > ./top_seqid.txt`;

# PyPlot: TM-Scores Histogram
`$bin_path/bin_plot/./prfHistogram.py top_tmscore.txt`;

# Print Suggestion for an "Ideal" TM-CutOff Score for Target Protein
open STMFILE, "<./sorted_tmscore.txt" or print "File sorted_tmscore.txt not found! \n";

my @linevar = <STMFILE>;
my @colvar;

print "\n";
print "Top-100 protein with highest similarity histogram has been generated...\n";
print "\n";
print "*******************************************************************\n";
print "Top-N proteins with global fold similarity with the target protein,\n";
print "has a minimum TM-Score value of:\n";
print "\n";

@colvar = split(/\s/, $linevar[9]);
print "Top-10: "; print "$colvar[0]\n";

@colvar = split(/\s/, $linevar[19]);
print "Top-20: "; print "$colvar[0]\n";

@colvar = split(/\s/, $linevar[29]);
print "Top-30: "; print "$colvar[0]\n";

@colvar = split(/\s/, $linevar[39]);
print "Top-40: "; print "$colvar[0]\n";

@colvar = split(/\s/, $linevar[49]);
print "Top-50: "; print "$colvar[0]\n";

close STMFILE;

print "\n";
print "*******************************************************************\n";

###############################################################################
# if the structural analogs obtained by TM-align is not enough, the profile 
# generated will be inaccurate for protein design; in this case, psi-blast is 
# used to search for sequence analogs to generate msa instead.
###############################################################################

my $analog_num = 0;
my $MIN_ANALOG_NUM = 10;
$analog_num = `cat msa_from_tmalign.txt | wc -l`;

if($analog_num < $MIN_ANALOG_NUM ){
  my $PDB2FAS = "$bin_path/PDB2FAS";
  my $database = "/nfs/amino-library/nr/nr";
  if($hostName =~ /comet/){
    $database = "/oasis/projects/nsf/mia181/zhanglab/library/nr/nr";
  }
  # convert pdb into fasta format;
  system "$PDB2FAS $querypdb > query.fasta";
  my $arguments = "-query query.fasta -out blast.xml -db $database -outfmt 5 -evalue 1e-4";
  # do psi-blast to search for sequence analogs and write in .xml format
  system "$bin_path/ncbi-blast-2.7.1+/bin/psiblast $arguments";
  system "$bin_path/xml2msa.pl  blast.xml > msa_from_psiblast.txt";

  # combine the msa obtained from tmalign & psi-blast
  my (@msa1, @msa2);
  open MSA1, "<msa_from_tmalign.txt";
  @msa1 = <MSA1>; chomp(@msa1);
  close MSA1;
  open MSA2, "<msa_from_psiblast.txt";
  @msa2 = <MSA2>; chomp(@msa2);
  close MSA2;
  system "cat msa_from_tmalign.txt > msa.txt";
  for(my $i = 0; $i < @msa2; $i++){
    my $seq2 = $msa2[$i];
    my $same = 0;
    for(my $j = 0; $j < @msa1; $j++){
      my $seq1 = $msa1[$j];
      if($seq1 eq $seq2){
        $same = 1;
        next;
      }
    }
    if($same == 0){
      system "echo $seq2 >> msa.txt";
    }
  }

  system "rm query.fasta";
}
else{
  `cp msa_from_tmalign.txt msa.txt`;
}

##############################################################################
# normally, we make profile using msa generated by TMalign, but if msa is not 
# enough, we make profile using msa generated by psi-blast
##############################################################################
# Make-Profile 
print "generating position specific scoring matrix... \n";
system "$bin_path/mkprf msa.txt > prf.txt";

# timestamping to end program
$finishtime = Benchmark->new;
$timespent = timediff($finishtime,$starttime);
print "generating Profile, computational time spent: ". timestr($timespent);


################################################################################
#### END 
################################################################################

my $timeEnd = localtime();

print "#normal ending of GenerateProfile.pl at $timeEnd\n";
print "\n";

