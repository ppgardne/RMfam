#!/usr/bin/perl

#search the motif CMs against a database of alignments or sequences. 

use warnings;
use strict;
use Getopt::Long;
use IO::File;
use Data::Dumper; #for debug

# Examples:
# k-turns:
# http://www.dundee.ac.uk/biocentre/nasg/kturn/kturns_known.php
# U4, SSU, LSU, SAM, Lysine, T-box, cyclic di-GMP-II, G2nA SAM, 
# possibly RNases P & MRP, metazoan SRP
# Csr-Rsm:
# RF00018-RsmB; RF00084-CsrC; RF00166-RsmZ; RF00195-RsmY; RF02144-RsmX
# sarcin-ricin:
# bacterial 23S ribosomal RNA, archeael 23S ribosomal RNA, archeal 5S ribosomal RNA, 
# Bacillus subtilis RNase P

my( $local, 
    $global,
    $thresh,
    $evalueThresh,
    $outfile,
    $fastaIn,
    $emblOut, 
    $gffOut, 
    $alnOut,
    $netOut,
    $pid,
    $clean,
    $verbose,
    $help
    );

my $VERSION = '0.01';
my $cmdline = $0.' '.join( ' ', @ARGV );
my $starttime = `date`;
chomp $starttime;
my $idf               = 0.9; #value for esl-weight: filter sequences by $idf fraction identity 
my $fractionMotifs    = 0.1; #fraction of sequences in a filtered alignment covered by a motif before inclusion
my $minNumberHits     = 2;   #minimum number of sequences covered by a motif before inclusion
&GetOptions( "local"         => \$local,
	     "global"        => \$global,
	     "t=s"           => \$thresh,
	     "e=s"           => \$evalueThresh,
	     "fm=s"          => \$fractionMotifs,
	     "min=s"         => \$minNumberHits,
	     "idf=s"         => \$idf,
	     "f|fasta"       => \$fastaIn,
	     "e|embl"        => \$emblOut,
	     "g|gff"         => \$gffOut,
	     "a|aln"         => \$alnOut,
	     "n|net"         => \$netOut,
	     "pid=s"         => \$pid,
	     "o=s"           => \$outfile,
	     "c|clean"       => \$clean,
	     "v"             => \$verbose,
	     "h"             => \$help
	     );

my $cmfile = shift;
my $infile = shift;

if( $help ) {
    &help();
    exit(0);
}
elsif (not -e $cmfile or not -e $infile){    
    print "MISSING an essential CM [$cmfile] or fasta/Stockholm [$infile] file!\n";
    &help();
    exit(1);        
}

$alnOut = 1 if ((not defined $gffOut) && (not defined $emblOut) && (not defined $alnOut));

my ($ffafile,$resfile, $cmsfile,$weights,$sumWeights);
if (not defined $fastaIn){
#IF STOCKHOLM FILE:
    ($ffafile,$resfile, $cmsfile) = ($pid . '.filtered.fasta', $pid . '.tblout', $pid . '.cmscan' ) if (defined $pid);
    print STDERR "Filter and reformat alignment [$infile]\n"          if( $verbose );
    ($ffafile,$weights,$sumWeights) = stockholm2filteredfasta($infile,$idf) if (not defined $pid);
}
else{
#IF FASTA FILE:
    $ffafile = $infile;
}

my $noSeqs = compute_number_of_seqs($ffafile);

print STDERR "run infernal search [$cmfile] against [$ffafile]\n"                         if( $verbose );
($resfile, $cmsfile) = run_infernal_search( $cmfile, $ffafile, $thresh, $evalueThresh )   if( not defined $pid);
print STDERR "parse infernal results [$resfile]\n"                                        if( $verbose );
my ($features,$motifLabels, $idCounts, $sumBits, $weightedSumBits) = parse_infernal_table( $resfile, $noSeqs, $fractionMotifs, $minNumberHits, $weights,$sumWeights, $fastaIn );

###############

if( defined $gffOut ) {
    print_gff($infile,$features,$starttime,$VERSION,$cmdline,$cmfile);
}

if( defined $emblOut ) {
    print "Working on EMBL still!\n";
}

if(defined $alnOut && not defined $fastaIn){    
    my $aOut = print_annotated_alignment($infile,$features,$motifLabels);
    system("esl-reformat pfam $aOut > $infile\.annotated.stk");
    print "created annotated alignment [$infile\.annotated.stk]\n" if( $verbose );
}

if( defined $netOut ) {
    #Print data for displaying as a network:
    #--alignment id, rmfam IDs, 
    #--scores: fraction seqs in alignment, sumBits, 
    my $nOut = print_network($infile,$features,$idCounts,$sumBits,$weightedSumBits, $noSeqs, $fractionMotifs, $minNumberHits, $sumWeights); 
    print "Network data written to [$nOut]\n" if( $verbose );
}

#CLEANUP FILES!
if(defined($clean)){
    unlink( $ffafile ) if (-s $ffafile);
    unlink( $resfile ) if (-s $resfile);
    unlink( $cmsfile ) if (-s $cmsfile);
    unlink( "$$.annotated.stk" ) if (-s "$$.annotated.stk");
}

exit(0);

sub stockholm2filteredfasta {
    my ($infile,$idf) = @_;
    
    system "esl-weight -f --idf $idf $infile > $$.filtered.stk && esl-reformat -r -u fasta $$.filtered.stk > $$.filtered.fasta"
	and die "FATAL: failed to run [esl-weight -f --idf $idf $infile > $$.filtered.stk && esl-reformat -r -u fasta $$.filtered.stk > $$.filtered.fasta]!\n[$!]";
    
    my ($weights,$sum) = stockholm2weights("$$.filtered.stk"); 
    
    return ("$$.filtered.fasta",$weights,$sum);
    
}

#stockholm2weights: Use esl-weights to obtain Gerstein/Sonnhammer/Chothia tree weights for each sequence
sub stockholm2weights {
    my $infile = shift;
    my $sum=0;
    my %weights;
    open(F, "esl-weight -g $infile | ") or die "FATAL: could not open pipe for [esl-weight -g $infile]\n[$!]";
    while(my $w = <F>){
	if($w=~/^#=GS\s+(\S+)\s+WT\s+(\S+)/){
	    $weights{$1}=$2; 
	    $sum+=$2;
	}
    }

return \%weights, $sum;
    
}


sub run_infernal_search {
    my ($cmfile, $fafile, $thresh, $evalueThresh) = @_;
    my $options = "";

    if( defined $thresh ) {
	$options = " -T $thresh ";
    }
    elsif ( defined $evalueThresh ){
	$options = " -E $evalueThresh ";
    }
    else {
	$options = " --cut_ga ";
    }
    
    if( $global ) {
	$options .= " -g";
    }
    elsif( $local ) {
	# default in infernal 1.0
    }
#    my $exe = "cmsearch --toponly --tabfile $$.tabfile $options $cmfile $fafile > $$.cmsearch";
#    my $exe = "cmsearch --toponly --tblout $$.tabfile $options $cmfile $fafile > $$.cmsearch";
    my $exe = "cmscan --max --toponly --tblout $$.tblout $options $cmfile $fafile > $$.cmscan";
 
    print "$exe\n" if (defined $verbose);
    system "$exe" and die "FATAL: failed to execute [$exe]\n[$!]";
    return ("$$.tblout", "$$.cmscan");
}

sub parse_infernal_table {
    my ($file, $noSeqs, $fractionMotifs, $minNumberHits, $weights,$sum,$fastaIn) = @_;
    my $fh;
    my %f2;
    my %idCounts;
    my (%sumBits,%weightedSumBits);
    my %seenSeqidRmfam;
    my %motifLabels;
    my %taken; 
    
    #Accounting loop:
    #(only really needed for annotating alignments)
    $fh = IO::File->new( $file );
    while(<$fh>) {
#	if( /^\#\s+CM:\s+(\S+)/ ) {
#	    $rmfamid = $1;
#	}
	
	next if( /^\#/ );
		
	my( $rmfamid, $model, $seqid, $start, $end, $modst, $moden, $bits, $evalue, $gc );
#	if( (( $model, $seqid, $start, $end, $modst, $moden, $bits, $evalue, $gc ) =
#	    /^\s*(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)$/) || 
#	            (( $seqid, $start, $end, $modst, $moden, $bits, $evalue, $gc ) =
#	            /^\s*(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)$/) ) {
	
	my @fields = split(/\s+/, $_);
	if(defined($fields[15]) && $fields[15]){
#1                   2         3                          4         5   6          7      8          9      10     11    12     13  14    15      16      17  18                   
#target name         accession query name                 accession mdl mdl from   mdl to seq from   seq to strand trunc pass   gc  bias  score   E-value inc description of target
#sarcin-ricin-2       -         CP000002.3/3366091-3365915 -          cm        1       61       18       82      +    no    1 0.55   0.0   19.5    0.0035 !   Sarcin-ricin motif 2
	    ($rmfamid,$seqid,$bits)=($fields[0],$fields[2],$fields[14]);
	    $sumBits{$rmfamid}+=$bits;
	    $weights->{$seqid} = 1.0 if (not defined $weights->{$seqid});
	    $weightedSumBits{$rmfamid}+=$bits*$weights->{$seqid};
	    
	    $idCounts{$rmfamid}++ if (not defined $seenSeqidRmfam{$seqid}{$rmfamid});
	    $seenSeqidRmfam{$seqid}{$rmfamid}++;
	    
	}
    }
    $fh->close;
    
    $fh = IO::File->new( $file );
    while(<$fh>) {
	# if( /^\#\s+CM:\s+(\S+)/ ) {
	#     $rmfamid = $1;
	# }
	
	next if( /^\#/ );
	#filters (only really needed for annotating alignments):
	
	
	my @fields = split(/\s+/, $_);

	my( $rmfamid, $model, $seqid, $start, $end, $strand, $modst, $moden, $bits, $evalue, $printEval, $gc );
#	if( (( $model, $seqid, $start, $end, $modst, $moden, $bits, $evalue, $gc ) =
#	    /^\s*(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)$/) || 
#	            (( $seqid, $start, $end, $modst, $moden, $bits, $evalue, $gc ) =
#	            /^\s*(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)$/) ) {
	if(defined($fields[15]) && $fields[15]){
#1                   2         3                          4         5   6          7      8          9      10     11    12     13  14    15      16      17  18                   
#target name         accession query name                 accession mdl mdl from   mdl to seq from   seq to strand trunc pass   gc  bias  score   E-value inc description of target
#sarcin-ricin-2       -         CP000002.3/3366091-3365915 -          cm        1       61       18       82      +    no    1 0.55   0.0   19.5    0.0035 !   Sarcin-ricin motif 2
	    
	    ($rmfamid, $seqid, $start, $end, $strand, $bits, $evalue)=($fields[0],$fields[2],$fields[7],$fields[8],$fields[9],$fields[14],$fields[15]);
	    
	    next if (($idCounts{$rmfamid} < $noSeqs*$fractionMotifs) && (not defined $fastaIn));
	    next if (($idCounts{$rmfamid} < $minNumberHits) && (not defined $fastaIn)); 
	    $motifLabels{$rmfamid} = assign_motif_label($rmfamid,\%taken) if not defined $motifLabels{$rmfamid};
	    $taken{$motifLabels{$rmfamid}}=1;
	    
	    my $strand = 1;
	    if( $end < $start ) {
		( $start, $end ) = ( $end, $start );
		$strand = -1;
	    }
	    
	    my $printEval = 'NA';
	    if( $evalue =~ /[0-9]/ ) {
		$printEval = $evalue;
	    }

	    my %f = ( seqid       => $seqid,
		      start       => $start,
		      end         => $end,
		      strand      => $strand,
		      score       => $bits,
		      evalue      => $printEval,
		      rmfamid     => $rmfamid,
		      label       => $motifLabels{$rmfamid}
		);

	    push( @{ $f2{$seqid} }, \%f ); #hash of arrays of hashes is a better structure for later (annotating alignments)
	    
	    #print "f:[[[[[[" . Dumper(%f) . "]]]]]]\n" if (defined $verbose);

	    
	}
    }
    $fh->close;
    
#    print "f2:[[[" . Dumper(%f2) . "]]]";
    
    return (\%f2, \%motifLabels, \%idCounts, \%sumBits, \%weightedSumBits);
}

sub assign_motif_label {
    
    my ($rmfamid,$taken)=@_;
    my @chars = qw(A B C D E F G H I J K L M N O P Q R S T U V W X Y Z a b c d e f g h i j k l m n o p q r s t u v w x y z 1 2 3 4 5 6 7 8 9 0);
    my @rmfamid = split(//,$rmfamid);
    foreach my $c (@rmfamid,@chars){
	next if ($c !~ /\w|\d/);
	return $c if (not defined $taken->{$c});
    }   
}

sub print_gff {
    my ($infile,$features,$starttime,$VERSION,$cmdline,$cmfile) = @_;
    my $outfile = $infile . ".gff"; 
    
    my $endtime = `date`;
    chomp $endtime;
    
    open(F, "> $outfile");
    print F "##gff-version 3
# rmfam_scan.pl (v$VERSION)
# command line:     $cmdline
# CM file:          $cmfile
# query file:       $infile
# start time:       $starttime
# end time:         $endtime\n";

    foreach my $seqid (keys %{$features}){
	foreach my $f ( @{$features->{$seqid}} ){
	    my ($start,$end, $char, $rmfamid, $score, $evalue) = ($f->{'start'}, $f->{'end'}, $f->{'label'}, $f->{'rmfamid'}, $f->{'score'}, $f->{'evalue'});
	    my $strand = '+';
	    
	    if ($start > $end){
		$strand = '-';
		($start,$end) = ($end,$start);
	    }	
	    my $eString='';
	    $eString="Evalue=$evalue;" if (defined($evalue) && $evalue ne '-' && $evalue ne 'NA');
	    print F "$seqid\tCMSCAN102\tmotif\t$start\t$end\t$score\t$strand\t.\tName=$rmfamid;$eString\n";	    
	}
    }
    close(F);
    return $outfile;
}

sub print_annotated_alignment{
    my ($infile, $features, $motifLabels) = @_;
    
    #1. CONVERT TO PFAM & READ SEQS INTO AN ARRAY (storing sequence positions in hashes)
    my $alnLength = compute_length_of_alignment($infile);
    open(F, "esl-reformat -ru --mingap pfam $infile | ") or die "FATAL: could not open pipe for reading $infile\n[$!]";
    my @stk = <F>;
    my %positions2seqid; 
    my %seqid2positions; 
    my %motifLines; #hash of arrays: seq1:MT.1,MT.2,... seq2:MT.1,MT.2,...
    my %motiffedSeqLineNumbers;
    my $cnt=0;
    my $firstSeqLine;

    foreach my $stk (@stk){
	#next if ($stk=~/^#/ or $stk=~/^\/\//);
	if($stk=~/^(\S+)\s+\S+$/){
	    $positions2seqid{$cnt}=$1;
	    $seqid2positions{$1}=$cnt;
	    $firstSeqLine=$cnt if not defined $firstSeqLine;
	}
	$cnt++;
    }
    print "Read [$cnt] sequences from [$infile]\n" if ($verbose);
    
    #2. CREATE MT LINES FOR EACH SEQ WITH AN ANNOTATION (TAKE CARE OF OVERLAPS)
    my @seqCoords2alnCoords;
    foreach my $seqid (keys %seqid2positions){
	if(defined $features->{$seqid}){
	    my $fpos=0;

	    foreach my $f ( @{$features->{$seqid}} ){
		
		#make a map from sequence to alignment coordinate space
		my $alnSeq=$stk[$seqid2positions{$seqid}];
		
		if($alnSeq=~/^(\S+)\s+(\S+)$/){
		    my ($lid,$lseq)=($1,$2);
		    #sanity check:
		    if($lid ne $seqid){
			print "WARNING: seqId:[$seqid] and alnSeq:[$alnSeq] don't match!\n";
			next;
		    }
		    
		    $alnSeq=$lseq; 
		    my @alnSeq = split(//, $alnSeq);
		    #sanity check:
		    if(scalar(@alnSeq) != $alnLength){
			printf "WARNING: the lengths [$alnLength]!=[%d] computed from seqId [$seqid] and alnSeq:\n[$alnSeq]\ndon't match! [CHECK: could be a gap-only column]\n", scalar(@alnSeq);
			#next;
		    }
		    
		    my ($aCnt,$sCnt)=(0,0); 
		    foreach my $as (@alnSeq){
			if(is_nucleotide($as)){
			    $seqCoords2alnCoords[$sCnt]=$aCnt;
			    $sCnt++;
			}
			$aCnt++;
		    }
		}
		else {
		    printf "WARNING: line number [%d] is meant to correspond to $seqid! Check the formatting.\n", $seqid2positions{$seqid};
		}
		
		#Choose a motif line for the annotation (avoiding overlaps)
		#--the motif line selection approach is to check for overlaps 
		#  with the seen motifs, incrementing the line by 1 for each overlap.
		my $mtCnt=0;
		my ($start,$end, $char, $rmfamid) = ($f->{'start'}, $f->{'end'}, $f->{'label'}, $f->{'rmfamid'});
		
		for (my $gpos=0; $gpos<$fpos; $gpos++){
		    my $g = ${ $features->{$seqid} }[$gpos];
		    $mtCnt++ if ( overlap($start,$end, $g->{'start'},$g->{'end'}) );
		}
		
		$motifLines{$seqid}[$mtCnt] = '.' x $alnLength if not defined $motifLines{$seqid}[$mtCnt];
		#my $lalnSeq=$stk[$seqid2positions{$seqid}];
		#if($lalnSeq =~ /^(\S+)\s+(\S+)$/){
		#    $lalnSeq=$2;
		#}
		#print "lalnSeq:[[[$lalnSeq]]]\n";
		#my @lalnSeq=split(//,$lalnSeq);
		#system("sfetch -d $pid.filtered.fasta -f $start -t $end -r \42$seqid s:$start e:$end as:$seqCoords2alnCoords[$start] ae:$seqCoords2alnCoords[$end] $mtCnt $char:$rmfamid\42 \42$seqid\42");
		for(my $mpos=$start-1; $mpos<$end; $mpos++){
		    my $aCoord = $seqCoords2alnCoords[$mpos];
		    #print $lalnSeq[$aCoord];
		    #if(is_nucleotide($lalnSeq[$aCoord])){
		    substr($motifLines{$seqid}[$mtCnt],$aCoord,1)=$char;
		    #}
		}
		$fpos++;
	    }
	    $motiffedSeqLineNumbers{$seqid2positions{$seqid}}=$seqid if defined($motifLines{$seqid}[0]); 
	}
    }
    
    
    #3. PRINT
    my $outFile="$$.annotated.stk";
    my $fh = IO::File->new(  );
    $fh->open("> $outFile");
    for(my $ii=0; $ii<scalar(@stk); $ii++){
	
	if ($ii == $firstSeqLine-1){
	    #print 1 letter codes to rmfamid mapping 
	    foreach my $l (sort {$motifLabels->{$a} cmp $motifLabels->{$b}} keys %{$motifLabels}){
		printf $fh "#=GF MT.%s   %s\n", $motifLabels->{$l}, $l; 
	    }
	}
	
	print $fh $stk[$ii];
	
	if (defined $motiffedSeqLineNumbers{$ii} && defined $positions2seqid{$ii}){
	    my $mCnt=0;
	    
	    foreach my $mt (@{$motifLines{$positions2seqid{$ii}}}){
		printf $fh "#=GR %s MT.$mCnt %s\n", $motiffedSeqLineNumbers{$ii}, $mt if (defined $mt);
		$mCnt++;
	    }
	}
	
    }
    
    return $outFile;
}

sub print_network {
    my ($infile, $features, $idCounts, $sumBits, $weightedSumBits, $noSeqs, $fractionMotifs, $minNumberHits, $sumWeights) = @_;
    my $outfile = $infile . ".network";
    my $aId = extract_id($infile);
    my $aAc = extract_acc($infile);
    open(F, "> $outfile");
    
    #becomes a mean if weights is not defined:
    $sumWeights = $noSeqs if (not defined $sumWeights);
    $sumWeights = 1.0     if ($sumWeights < 1.0);
    foreach my $fId (keys %{$weightedSumBits}){
	$weightedSumBits->{$fId} = $weightedSumBits->{$fId}/$sumWeights; 	
    }
    
    foreach my $mId (sort keys %{$idCounts}){
	next if ( $idCounts->{$mId}/$noSeqs < $fractionMotifs );
	next if ( $idCounts->{$mId}         < $minNumberHits  );
	printf F "%10s %20s %20s " . " "x10 . "%0.2f" . " "x10 . " %0.2f" . " "x10 . " %0.2f\n", 
	$aAc,$aId,$mId, $idCounts->{$mId}/$noSeqs, $sumBits->{$mId}, $weightedSumBits->{$mId};
    }
    close(F);
    
    return $outfile;
}




######################################################################
#Utilities:
######################################################################
sub compute_number_of_seqs {    
    my $file = shift;
    my $numberOfSeqs;
    open(ALI,"esl-seqstat --rna $file |") or die "FATAL: could not open [esl-seqstat $file] pipe:[$!]";
    #Grab the fields out of seqstat
    while(<ALI>) {
	if (/^Number of sequences:\s+(\d+)/){
	    $numberOfSeqs=$1;
	}
    }
    close(ALI);
    print "WARNING: numberOfSeqs is undefined for [esl-seqstat $file]!" if not defined $numberOfSeqs;
    return $numberOfSeqs;
}

######################################################################
sub compute_length_of_alignment {    
    my $file = shift;
    my $alnLength;
    open(ALI,"esl-alistat --rna $file |") or die "FATAL: could not open [esl-seqstat $file] pipe:[$!]";
    #Grab the fields out of seqstat
    while(<ALI>) {
	if (/^Alignment length:\s+(\d+)/){
	    $alnLength=$1;
	}
    }
    close(ALI);
    print "WARNING: alnLength is undefined for [esl-seqstat $file]!" if not defined $alnLength;
    return $alnLength;
}

######################################################################
#returns true if input character is a nucleotide (IUPAC codes):
sub is_nucleotide {
    my $a = shift;
    
    if (defined($a)){
	$a =~ tr/a-z/A-Z/;
    }
    
    if (defined($a) && length($a) && ($a =~ /[ACGUTRYWSMKBDHVN]/) ){
	return 1;
    }
    else {
	return 0;
    }
    
}

######################################################################
#Returns true if the coordinates for two regions ($x1, $y1) and ($x2, $y2) overlap:
# - assumes that $x1 < $y1 and $x2 < $y2.
sub overlap {
    my($x1, $y1, $x2, $y2) = @_;
    
    if ( ($x1<=$x2 && $x2<=$y1) || ($x1<=$y2 && $y2<=$y1) || ($x2<=$x1 && $x1<=$y2) || ($x2<=$y1 && $y1<=$y2)  ){
        return 1;
    }
    else {
        return 0;
    }
}

######################################################################
#Extract an identifier from a file, if not found then use filename
sub extract_id {
    my $aFile=shift;
    open(F, "< $aFile");
    
    while(my $l = <F>){
	if($l=~/\#=GF\s+ID\s+(\S+)/){
	    return $1;
	}
    }
    return $aFile;
}

######################################################################
#Extract an accession from a file, if not found then use filename
sub extract_acc {
    my $aFile=shift;
    open(F, "< $aFile");
    
    while(my $l = <F>){
	if($l=~/\#=GF\s+AC\s+(\S+)/){
	    return $1;
	}
    }
    return $aFile;
}

######################################################################
#Extract a GA threshold from a file, if not found then use filename
sub extract_ga {
    my $aFile=shift;
    open(F, "< $aFile");
    
    while(my $l = <F>){
	if($l=~/\#=GF\s+GA\s+(\S+)/){
	    return $1;
	}
    }
    return 0.0;
}


sub help {
    print STDERR <<EOF;

$0: search a fasta or stockholm file with RMfam covariance models

Usage: $0 <options> cm_file Stockholm_file/fasta_file
    Options
        -h              : show this help
	-v              : verbose - prints lots of largely unimportant information
	-c|--clean      : delete unnecessary results files
    Cmscan options:
	-t <bits>       : specify cutoff in bits      [DEFAULT is to use curated GA thresholds]
	--local         : perform local mode search   [DEFAULT]
	--global        : perform global mode search

    RMfam options:
	-f|--fasta      : search database is in fasta format
	For searching alignments:
	-fm  <num>      : fraction of sequences a motif must hit before inclusion [DEFAULT: $fractionMotifs]
	-min <num>      : number   of sequences a motif must hit before inclusion [DEFAULT: $minNumberHits]
	-idf <num>      : filter out seqs w/ fractional ident > num [DEFAULT: $idf]
	
    Output options:
	-o              : output file
	-e|--embl       : output in EMBL format
	-g|--gff        : output in GFF format                 [DEFAULT for fasta input]
	-a|--aln        : output in annotated Stockholm format [DEFAULT for Stockholm input]
	-n|--net        : output in tabular format for network visualisation
	
    Miscellaneous options:
	--pid           : restart a job using precomputed cmscan results
	
      TODO:
	--try a tree weighting on the sum of bit scores metric
	--allow a threshold on the sumBits score
	--add some more caveats to what gets summed for the "fm" option <a fuzzy alignment approach - window can be proportional to the specificity of the motif model>
	--record model specificity in the alignment?
	--add support for HMM and-or PWMs 
        --add secondary structure information to the annotated alignment?
            --try to include motif secondary structures to annotations 
              (using the cmscan outfile)
	--add an unstranded search option (removes --toponly from cmscan) &/or 
          search the reverse complement
	--for fasta files present a better summary. Eg. sumBits (not clans), list motifs, ...
	
	--test unfiltered runs
	--try multiple CMs for each family with eg. base-triples, pseudoknots, ...
	--only report the highest scoring terminator motif (eg. compete "clans")
	
	--add a test for concordance with the consensus structure

	--update code for Infernal 1.rc2
	--add support for Elena-s alignment annotation format

	--where are the missing pseudoknot annotations?
EOF
}





