#!/usr/bin/perl

use warnings;
use strict;

#grep -v ^"#" true_positives.txt | perl -lane 'if($F[1] eq "k-turn"){print "$F[0]\tk-turn-1\t$F[2]\n$F[0]\tk-turn-2\t$F[2]"}elsif($F[1] eq "sarcin-ricin"){print "$F[0]\tsarcin-ricin-1\t$F[2]\n$F[0]\tsarcin-ricin-2\t$F[2]"}elsif($F[1] eq "Terminator"){print "$F[0]\tterminator1\t$F[2]\n$F[0]\tterminator2\t$F[2]"}else{$F=join("\t", @F); print $F}' | sort -d > true_positives_processed.txt

my $totalTrues=0;
my $totalFalse=2208*10; #number of shuffled alignments
my %totalTruesB2;
my %rmfam;
my %seenMatch;
my $trueFile = "true_positives.txt"; 
open(TF, "< $trueFile"); 
my (%fam2mots, %fammotcat);
while(my $t=<TF>){
    
    my @t = split(/\s+/, $t); 
    if(scalar(@t)>1){
	push(@{ $fam2mots{$t[0]} },$t[1]);
	$fammotcat{"$t[0]:$t[1]"}=1;
	$seenMatch{"$t[0]:$t[1]"}=0;
	$totalTrues++;
	my $rmfam = $t[1];
	if (not defined($totalTruesB2{$rmfam})){
	    $totalTruesB2{$rmfam}=0; 
	}
	$totalTruesB2{$rmfam}++;
	$rmfam{$rmfam}=1;

    }
    
}
close(TF);

my $cmFile = "../scripts/RMfam.cm"; 
open(CF, "< $cmFile"); 
open(UT, "> thresholds.dat"); 
print UT "NAME\tGA\n";
my ($motName, $motGA, %motifs);
while(my $cf=<CF>){
    if ($cf=~/^NAME\s+(\S+)/){
	$motName=$1; 
    }
    elsif($cf=~/^GA\s+(\S+)/){
	$motGA=$1;
	print UT "$motName\t$motGA\n";
	$motifs{$motName}=$motGA;
    }
}
close(CF);
close(UT);


######################################################################
#ALIGNMENTS

my (@mappings);
#my $netFile  = "network.dat";
my $netFile  = "../annotations/Rfam.tab";
#array of arrays
open(NF, "sort -k5nr < $netFile |"); 
open(UTT, "> true-predictions-alignments.dat"); 
open(UTF, "> false-predictions-alignments.dat"); 
$totalTrues=0; 
while(my $n=<NF>){
    
    my @n=split(/\s+/, $n);
    if($n=~/RF\d{5}\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/){
	my ($rfam,$rmfam,$n1,$n2,$n3)=($1,$2,$3,$4,$5);
	if( defined($fammotcat{"$rfam:$rmfam"}) ){
	    $seenMatch{"$rfam:$rmfam"}++;
	    #TRUE!
	    push(@mappings, [ $rfam,$rmfam,$n1,$n2,$n3, "true" ]); 
	    $totalTrues++;
	    $rmfam{$rmfam}=1;
	    print UTT "$n";
	}
	else{
	    print UTF "$n";
	    #print "FALSE: $n";
	}
    }
}
close(NF); 
close(UTT); 
close(UTF); 

open(UT, "> false-negative-alignments.dat");
foreach my $k ( sort { $a cmp $b } keys %seenMatch){
    my @k = split(/:/, $k);
    print UT "$k[0]\t$k[1]\n" if ($seenMatch{$k} == 0);
}
close(UT);

my $netSFile = "../annotations/Rfam-shuffle.tab";
$totalFalse=0;
my %totalFalseB2;
open(NF, "< $netSFile"); 
while(my $n=<NF>){
    
    my @n=split(/\s+/, $n);
    push(@mappings, [ @n[1..5], "false" ]); 
    $totalFalse++;
    
    $totalFalseB2{$n[2]}=0 if (not defined($totalFalseB2{$n[2]})); 
    $totalFalseB2{$n[2]}++;
    
}
close(NF); 


#TP, TN, FP, FN

my (%tpB2,%tnB2,%fpB2,%fnB2); 
foreach my $rm (keys %rmfam){
    $totalFalseB2{$rm}=0 if (not defined($totalFalseB2{$rm}));
    $totalTruesB2{$rm} =0 if (not defined($totalTruesB2{$rm}));
    ($tpB2{$rm},$tnB2{$rm},$fpB2{$rm},$fnB2{$rm})=(0,2208*10,0,$totalTruesB2{$rm}); #$totalFalseB2{$rm}
    #print "motif:[$rm] tF[$totalFalseB2{$rm}] tT[$totalTruesB2{$rm}]\n"; 
}

#sorted on sum-bits
my @labels=("fracSeqs", "sumBits", "weightedSumBits"); 
open(UT2, "> ROC_data_sumBits_benchmark2\.dat");
print UT2 "bits\tTP\tTN\tFP\tFN\tMCC\tSENS\tSPEC\tPPV\tNPV\tACC\tFDR\tmotif\n";
my(%mccB2, %sensB2, %specB2, %ppvB2, %npvB2, %accB2, %fdrB2);
for(my $i=0; $i<scalar(@labels); $i++){
    my $arrayCoord=$i+2; 
    @mappings = sort { $b->[$arrayCoord] <=> $a->[$arrayCoord] } @mappings; 
    my ($tp,$tn,$fp,$fn)=(0,$totalFalse,0,$totalTrues); #$totalFalse
    
    open(UT, "> ROC_data_$labels[$i]\.dat");
    
    print UT "$labels[$i]\tTP\tTN\tFP\tFN\tMCC\tSENS\tSPEC\tPPV\tNPV\tACC\tFDR\n";
    foreach my $row (@mappings){
	
	my $rmfam=$row->[1];
	my $str=join("][", @{$row});
	if($row->[5] eq "true"){
	    $tp++;
	    $fn--;
	    if($labels[$i] eq "sumBits"){
		$tpB2{$rmfam}++;
		$fnB2{$rmfam}--;	    
	    }
	}
	else{
	    $fp++; 
	    $tn--;
	    if($labels[$i] eq "sumBits"){
		$fpB2{$rmfam}++;
		$tnB2{$rmfam}--;
	    }
	}
	
	my($mcc, $sens, $spec, $ppv, $npv, $acc, $fdr)=acc($tp,$tn,$fp,$fn); 
	printf UT "$row->[$arrayCoord]\t$tp\t$tn\t$fp\t$fn\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\n", $mcc, $sens, $spec, $ppv, $npv, $acc, $fdr; 
	
	if($labels[$i] eq "sumBits" && defined($tpB2{$rmfam}) && defined($tnB2{$rmfam}) && defined($fpB2{$rmfam}) && defined($fnB2{$rmfam}) ){
	    ($mccB2{$rmfam}, $sensB2{$rmfam}, $specB2{$rmfam}, $ppvB2{$rmfam}, $npvB2{$rmfam}, $accB2{$rmfam}, $fdrB2{$rmfam})=acc($tpB2{$rmfam},$tnB2{$rmfam},$fpB2{$rmfam},$fnB2{$rmfam}); 
	    printf UT2 "$row->[$arrayCoord]\t$tpB2{$rmfam}\t$tnB2{$rmfam}\t$fpB2{$rmfam}\t$fnB2{$rmfam}\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t$rmfam\n", $mccB2{$rmfam}, $sensB2{$rmfam}, $specB2{$rmfam}, $ppvB2{$rmfam}, $npvB2{$rmfam}, $accB2{$rmfam}, $fdrB2{$rmfam}; 
	}
    }
    
    close(UT); 
}
close(UT2); 

######################################################################
#SEQUENCES
@mappings=();
$totalTrues=0;
my $file="../annotations/Rfam-sequences.tab";
my %totalTruesSeq;
if(-e $file){
    open(IN, "< $file"); 
    while(my $in=<IN>){
	chomp($in); 
	my @in=split(/\t/, $in);
	if( defined($fammotcat{"$in[0]:$in[1]"}) ){
	    #TRUE!
	    push(@mappings, [ $in[2], "true", $in[0], $in[1] ]); 
	    $totalTrues++;
	    
	    $totalTruesSeq{$in[1]}=0 if (not defined($totalTruesSeq{$in[1]}));
	    $totalTruesSeq{$in[1]}++;
	}
    }
    close(IN); 
}

my %totalFalseSeq;
$totalFalse=0;
$file="../annotations/Rfam-shuffle-sequences.tab"; 
if(-e $file){
    open(IN, "< $file"); 
    while(my $in=<IN>){
	chomp($in); 
	my @in=split(/\t/, $in);
	
	if( scalar(@in) == 3 && $in[2] =~ /^\d+\.\d+$/){
	    #DAS IST FALSCH!
	    push(@mappings, [ $in[2], "false", $in[0], $in[1] ]); 
	    $totalFalse++;
	    $totalFalseSeq{$in[1]}=0 if (not defined($totalFalseSeq{$in[1]}));
	    $totalFalseSeq{$in[1]}++;
	}
	
    }
    close(IN); 
}


my (%tpSeq,%tnSeq,%fpSeq,%fnSeq); 
foreach my $rm (keys %motifs){
    $totalFalseSeq{$rm}=0 if (not defined($totalFalseSeq{$rm}));
    $totalTruesSeq{$rm}=0 if (not defined($totalTruesSeq{$rm}));
    ($tpSeq{$rm},$tnSeq{$rm},$fpSeq{$rm},$fnSeq{$rm})=(0,$totalFalseSeq{$rm},0,$totalTruesSeq{$rm}); 
    print "$rm:[$tpSeq{$rm},$tnSeq{$rm},$fpSeq{$rm},$fnSeq{$rm}]\n";
}

my ($tp,$tn,$fp,$fn)=(0,$totalFalse,0,$totalTrues); 
@mappings = sort { $b->[0] <=> $a->[0] } @mappings; 

open(UT,  "> ROC_data_single_seq_sumBits.dat");
open(UTS, "> ROC_data_single_seq_sumBits-motifPerMotif.dat");

print UT  "sngl_seq_sumBits\tTP\tTN\tFP\tFN\tMCC\tSENS\tSPEC\tPPV\tNPV\tACC\tFDR\n";
print UTS "sngl_seq_sumBits\tTP\tTN\tFP\tFN\tMCC\tSENS\tSPEC\tPPV\tNPV\tACC\tFDR\tmotif\n";
foreach my $row (@mappings){
    
    my $rmfam=$row->[3];
    if($row->[1] eq "true"){
	$tp++;
	$fn--;
	$tpSeq{$rmfam}++;
	$fnSeq{$rmfam}--;
    }
    else{
	$fp++; 
	$tn--;
	$fpSeq{$rmfam}++;
	$tnSeq{$rmfam}--;
    }
    
    my ($mcc, $sens, $spec, $ppv, $npv, $acc, $fdr)=acc($tp,$tn,$fp,$fn); 
    printf UT "$row->[0]\t$tp\t$tn\t$fp\t$fn\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\n", $mcc, $sens, $spec, $ppv, $npv, $acc, $fdr; 
    
    ($mcc, $sens, $spec, $ppv, $npv, $acc, $fdr)=acc($tpSeq{$rmfam},$tnSeq{$rmfam},$fpSeq{$rmfam},$fnSeq{$rmfam}); 
    printf UTS "$row->[0]\t$tpSeq{$rmfam}\t$tnSeq{$rmfam}\t$fpSeq{$rmfam}\t$fnSeq{$rmfam}\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t$rmfam\n", $mcc, $sens, $spec, $ppv, $npv, $acc, $fdr; 
    
}
close(UT); 
close(UTS); 

######################################################################
#Benchmark 2: benchmarking on the training data: Rmfam seqs = true, shuffled negative control sequences
#echo -e "db\tseq\tbits\teval\tmotif" > RMfam/benchmark/benchmark2.txt && cat */seq-scores.dat | egrep '^SEED|^PDB-shuffle' | sort -k3nr >> RMfam/benchmark/benchmark2.txt

open(B2, "< benchmark2.txt");
my @b2=<B2>;
close(B2);
#my ($tp,$tn,$fp,$fn)=(0,$totalFalse,0,$totalTrues); 
my (%totalFalse, %totalTrues);
($totalFalse,$totalTrues)=(0,0);
foreach my $b2 (@b2){
    next if $b2 =~ /^db/; 
    chomp($b2);
    
    my @b=split(/\t/, $b2);
    ($totalFalse{$b[4]}, $totalTrues{$b[4]})=(0,0) if (not defined($totalTrues{$b[4]})); 
    
    if($b2=~/^SEED-shuffled/ or $b2=~/^PDB-shuffle/){
	$totalFalse{$b[4]}++;	
	$totalFalse++;
    }
    elsif($b2=~/^SEED\t/){
	$totalTrues{$b[4]}++;
	$totalTrues++;
    }
}

($tp,$tn,$fp,$fn)=(0,$totalFalse,0,$totalTrues); 
my (%tp,%tn,%fp,%fn);
foreach my $k (keys %totalTrues){
    ($tp{$k},$tn{$k},$fp{$k},$fn{$k})=(0,$totalFalse{$k},0,$totalTrues{$k}); 
}


open(UT,  "> ROC_benchmark2.dat");
open(UTA, "> ROC_benchmark2-all.dat");
print UT  "bits\tTP\tTN\tFP\tFN\tMCC\tSENS\tSPEC\tPPV\tNPV\tACC\tFDR\tmotif\n";
print UTA "bits\tTP\tTN\tFP\tFN\tMCC\tSENS\tSPEC\tPPV\tNPV\tACC\tFDR\n";
foreach my $row (@b2){
    chomp($row);
    my @row = split(/\t/, $row); 
    next if $row =~ /^db/;
    next if (not defined($row[4]));
    if($row[0] eq "SEED"){
	$tp{$row[4]}++;
	$fn{$row[4]}--;
	$fn{$row[4]} = 0 if $fn{$row[4]}<0;
	$tp++;
	$fn--;
    }
    elsif($row[0] eq "SEED-shuffled" or $row[0] eq "PDB-shuffle"){
	$fp{$row[4]}++; 
	$tn{$row[4]}--;
	$tn{$row[4]} = 0 if $tn{$row[4]}<0;
	$fp++;
	$tn--;
    }
    
    my($mcc, $sens, $spec, $ppv, $npv, $acc, $fdr)=acc($tp{$row[4]},$tn{$row[4]},$fp{$row[4]},$fn{$row[4]}); 
    printf UT  "$row[2]\t$tp{$row[4]}\t$tn{$row[4]}\t$fp{$row[4]}\t$fn{$row[4]}\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t$row[4]\n", $mcc, $sens, $spec, $ppv, $npv, $acc, $fdr; 
    ($mcc, $sens, $spec, $ppv, $npv, $acc, $fdr)=acc($tp,$tn,$fp,$fn); 
    printf UTA "$row[2]\t$tp\t$tn\t$fp\t$fn\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\n", $mcc, $sens, $spec, $ppv, $npv, $acc, $fdr; 
    
}
close(UT);
close(UTA);

######################################################################
#PLOT RESULTS

system("R CMD BATCH --no-save plotROC.R");

exit(0); 

######################################################################
sub acc{
    
    my ($tp,$tn,$fp,$fn)=@_;
    
    my $mcc=0.0; 
    my $denom=($tp+$fp)*($tp+$fn)*($tn+$fp)*($tn+$fn);
    $mcc=($tp*$tn-$fp*$fn)/sqrt($denom) if ($denom>0);
    
    my $sens=0.0; 
    $sens = $tp/($tp+$fn) if(($tp+$fn)>0);
    
    my $spec=0.0; 
    $spec = $tn/($tn+$fp) if(($tn+$fp)>0);
    
    my $ppv=0.0; 
    $ppv = $tp/($tp+$fp) if(($tp+$fp)>0);
    
    my $npv=0.0; 
    $npv = $tn/($tn+$fn) if(($tn+$fn)>0);
    
    my $acc=0.0;
    $acc = ($tp+$tn)/($tp+$tn+$fp+$fn) if(($tp+$tn+$fp+$fn)>0); 
    
    my $fdr=0.0;
    $fdr = ($fp)/($tp+$fp) if(($tp+$fp)>0); 
        
    return($mcc, $sens, $spec, $ppv, $npv, $acc, $fdr); 
    
}
