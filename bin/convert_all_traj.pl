#!/usr/bin/env perl

my $datadir="$ARGV[0]";
my $bindir="$ARGV[1]";

my @decoys = qw(
    remc30.pdb
    remc31.pdb
    remc32.pdb
    remc33.pdb
    remc34.pdb
    remc35.pdb
    remc36.pdb
    remc37.pdb
    remc38.pdb
    remc39.pdb
);


system("mkdir $datadir/all_decoys_intermediate");
system("mkdir $datadir/all_decoys");
my $cycle=31;
my $energy;
open(FILE_ALL,">$datadir/all_decoys/all_decoys.pdb");
foreach my $decoy (@decoys){
    my $flagFirst = "true";
    my $num=0;
    my @decoy_coord=();
    chomp(my @all_decoys = `cat $datadir/$decoy`);
    foreach my $line (@all_decoys){
        chomp(my @splitLine = split(/\s+/,$line));
        if($splitLine[0]){
            if($flagFirst ne "true"){
                open(FILE,">$datadir/all_decoys_intermediate/decoy_$cycle\_$num");
                foreach my $decoy_line (@decoy_coord){
                    print FILE "$decoy_line\n";
                }
                close(FILE);
                system("$bindir/convert.py $datadir/all_decoys_intermediate/decoy_$cycle\_$num $datadir/all_decoys_intermediate/decoy_$cycle\_$num.pdb");
                system("$bindir/pulchra $datadir/all_decoys_intermediate/decoy_$cycle\_$num.pdb");
                chomp(my @full_atom_coord = `tail -n +2 $datadir/all_decoys_intermediate/decoy_$cycle\_$num.rebuilt.pdb`);
                print FILE_ALL "REMARK Decoy from REMC Replica $cycle Cycle $num with Energy $energy\n";
                foreach my $f_decoy_line (@full_atom_coord){
                    print FILE_ALL "$f_decoy_line\n";
                }
                @decoy_coord=();
            }
            else{
                $flagFirst="false";
            }
            $energy="$splitLine[1]";
            $num++;
        }
        else{
            push(@decoy_coord,$line);
        }
    }
    $cycle++;
}
close(FILE_ALL);

system("rm -rf $datadir/all_decoys_intermediate");
system("cd $datadir/all_decoys && tar -cvf all_decoys.pdb.tar all_decoys.pdb");
system("bzip2 $datadir/all_decoys/all_decoys.pdb.tar");
system("cp $datadir/all_decoys/all_decoys.pdb.tar.bz2 $datadir/final_designs/");
system("rm -rf $datadir/all_decoys");
