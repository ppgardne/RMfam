#!/usr/bin/perl

use warnings;
use strict;

#grep -v ^"#" true_positives.txt | perl -lane 'if($F[1] eq "k-turn"){print "$F[0]\tk-turn-1\t$F[2]\n$F[0]\tk-turn-2\t$F[2]"}elsif($F[1] eq "sarcin-ricin"){print "$F[0]\tsarcin-ricin-1\t$F[2]\n$F[0]\tsarcin-ricin-2\t$F[2]"}elsif($F[1] eq "Terminator"){print "$F[0]\tterminator1\t$F[2]\n$F[0]\tterminator2\t$F[2]"}else{$F=join("\t", @F); print $F}' | sort -d > true_positives_processed.txt

my $trueFile = "true_positives.txt"; 
open(TF, "< $trueFile"); 
my (%fam2mots, %fammotcat);
while(my $t=<TF>){
    
    my @t = split(/\s+/, $t); 
    if(scalar(@t)>1){
	push(@{ $fam2mots{$t[0]} },$t[1]);
	$fammotcat{"$t[0]:$t[1]"}=1;
    }
    
}
close(TF);

######################################################################
#ALIGNMENTS

my (@mappings);
#my $netFile  = "network.dat";
my $netFile  = "../annotations/Rfam.tab";
my $totalTrues=0;
#array of arrays
open(NF, "< $netFile"); 
while(my $n=<NF>){
    
    my @n=split(/\s+/, $n);
    if($n=~/RF\d{5}\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/){
	my ($rfam,$rmfam,$n1,$n2,$n3)=($1,$2,$3,$4,$5);
	if( defined($fammotcat{"$rfam:$rmfam"}) ){
	    #TRUE!
	    push(@mappings, [ $rfam,$rmfam,$n1,$n2,$n3, "true" ]); 
	    $totalTrues++;
	    #print "TRUE: $n";
	}
	else{
	    print "FALSE: $n";
	}
    }
}
close(NF); 

#exit(0);

my $netSFile = "../annotations/Rfam-shuffle.tab";
my $totalFalses=0;
open(NF, "< $netSFile"); 
while(my $n=<NF>){
    
    my @n=split(/\s+/, $n);
    push(@mappings, [ @n[1..5], "false" ]); 
    $totalFalses++;
    
}
close(NF); 


#TP, TN, FP, FN

#sorted on sum-bits
my @labels=("fracSeqs", "sumBits", "weightedSumBits"); 
for(my $i=0; $i<scalar(@labels); $i++){
    my $arrayCoord=$i+2; 
    @mappings = sort { $b->[$arrayCoord] <=> $a->[$arrayCoord] } @mappings; 
    my ($tp,$tn,$fp,$fn)=(0,$totalFalses,0,$totalTrues); 
    
    open(UT, "> ROC_data_$labels[$i]\.dat");
    
    print UT "$labels[$i]\tTP\tTN\tFP\tFN\tMCC\tSENS\tSPEC\tPPV\tNPV\tACC\tFDR\n";
    foreach my $row (@mappings){
	
	my $str=join("][", @{$row});
	if($row->[5] eq "true"){
	    $tp++;
	    $fn--;
	}
	else{
	    $fp++; 
	    $tn--;
	}
	
	my($mcc, $sens, $spec, $ppv, $npv, $acc, $fdr)=acc($tp,$tn,$fp,$fn); 
	printf UT "$row->[$arrayCoord]\t$tp\t$tn\t$fp\t$fn\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\n", $mcc, $sens, $spec, $ppv, $npv, $acc, $fdr; 
	
    }
    
    close(UT); 
}

######################################################################
#SEQUENCES
@mappings=();
$totalTrues=0;
foreach my $f (keys %fam2mots){
    my $file="../data/sequence-annotations/rfam/" . $f  . "-seed.stk.selectn5.network";
    
    if(-e $file){
	open(IN, "< $file"); 
	while(my $in=<IN>){
	    chomp($in); 
	    my @in=split(/\t/, $in);
	    if( defined($fammotcat{"$f:$in[1]"}) ){
		#TRUE!
		push(@mappings, [ $in[2], "true" ]); 
		$totalTrues++;
	    }	    
	}
	close(IN); 
    }
}

$totalFalses=0;
my @files=glob("../data/sequence-annotations/rfam-shuffle/*-shuffle-selectn5.fa.network"); 
foreach my $file (@files){
    if(-e $file){
	open(IN, "< $file"); 
	while(my $in=<IN>){
	    chomp($in); 
	    my @in=split(/\t/, $in);
	    
	    if( defined($in[1]) && $file=~/rfam-shuffle\/(\S+)-selectn5.fa.network/){
		#FALSE!
		push(@mappings, [ $in[2], "false" ]); 
		$totalFalses++;
	    }
	    
	}
	close(IN); 
    }
}


my ($tp,$tn,$fp,$fn)=(0,$totalFalses,0,$totalTrues); 
@mappings = sort { $b->[0] <=> $a->[0] } @mappings; 

open(UT, "> ROC_data_single_seq_sumBits\.dat");

print UT "sngl_seq_sumBits\tTP\tTN\tFP\tFN\tMCC\tSENS\tSPEC\tPPV\tNPV\tACC\tFDR\n";
foreach my $row (@mappings){
    
    my $str=join("][", @{$row});
    if($row->[1] eq "true"){
	$tp++;
	$fn--;
    }
    else{
	$fp++; 
	$tn--;
    }
    
    my($mcc, $sens, $spec, $ppv, $npv, $acc, $fdr)=acc($tp,$tn,$fp,$fn); 
    printf UT "$row->[0]\t$tp\t$tn\t$fp\t$fn\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\n", $mcc, $sens, $spec, $ppv, $npv, $acc, $fdr; 
    
}
close(UT); 

######################################################################
#Benchmark 2:
#echo -e "db\tseq\tbits\teval\tmotif" > RMfam/benchmark/benchmark2.txt && cat */seq-scores.dat | egrep '^SEED|^PDB-shuffle' | sort -k3nr >> RMfam/benchmark/benchmark2.txt

open(B2, "< benchmark2.txt");
my @b2=<B2>;
close(B2);
#my ($tp,$tn,$fp,$fn)=(0,$totalFalses,0,$totalTrues); 
my (%totalFalses, %totalTrues);
foreach my $b2 (@b2){
    next if $b2 =~ /^db/; 
    chomp($b2);
    
    my @b=split(/\t/, $b2);
    ($totalFalses{$b[4]}, $totalTrues{$b[4]})=(0,0) if (not defined($totalTrues{$b[4]})); 
    
    if($b2=~/^SEED-shuffled/){
	$totalFalses{$b[4]}++;	
    }
    elsif($b2=~/^SEED\t/){
	$totalTrues{$b[4]}++;
    }
}

my (%tp,%tn,%fp,%fn);
foreach my $k (keys %totalTrues){
    ($tp{$k},$tn{$k},$fp{$k},$fn{$k})=(0,$totalFalses{$k},0,$totalTrues{$k}); 
}

open(UT, "> ROC_benchmark2\.dat");
print UT "bits\tTP\tTN\tFP\tFN\tMCC\tSENS\tSPEC\tPPV\tNPV\tACC\tFDR\tmotif\n";
foreach my $row (@b2){
    chomp($row);
    my @row = split(/\t/, $row); 
    next if $row =~ /^db/;
    next if (not defined($row[4]));
    if($row[0] eq "SEED"){
	$tp{$row[4]}++;
	$fn{$row[4]}--;
    }
    else{
	$fp{$row[4]}++; 
	$tn{$row[4]}--;
    }
    
    my($mcc, $sens, $spec, $ppv, $npv, $acc, $fdr)=acc($tp{$row[4]},$tn{$row[4]},$fp{$row[4]},$fn{$row[4]}); 
    printf UT "$row[2]\t$tp{$row[4]}\t$tn{$row[4]}\t$fp{$row[4]}\t$fn{$row[4]}\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t$row[4]\n", $mcc, $sens, $spec, $ppv, $npv, $acc, $fdr; 
     
}
close(UT);

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
