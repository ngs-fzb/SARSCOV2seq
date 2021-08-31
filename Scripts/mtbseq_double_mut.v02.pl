#!/usr/bin/perl

################################################################
# Script for changing variant files so that multiple changes within a codon are called correctly
# and insertions/deletion are now a one line output, not a one line per allele output

#v.02 The same as v.06 of multi codon change script
################################################################

use strict;
use warnings;
no warnings "numeric";
no warnings "uninitialized";
use 5.010;
use Data::Dumper;
use IO::Handle;
use Getopt::Std;
use File::Basename;

my %codon = ('TCA'=>'Ser','TCC'=>'Ser','TCG'=>'Ser','TCT'=>'Ser',
			'TTC'=>'Phe','TTT'=>'Phe',
			'TTA'=>'Leu','TTG'=>'Leu',
			'TAC'=>'Tyr','TAT'=>'Tyr',
			'TAA'=>'_',
			'TAG'=>'_',
			'TGC'=>'Cys','TGT'=>'Cys',
			'TGA'=>'_',
			'TGG'=>'Trp',
			'CTA'=>'Leu','CTC'=>'Leu','CTG'=>'Leu','CTT'=>'Leu',
			'CCA'=>'Pro','CCC'=>'Pro','CCG'=>'Pro','CCT'=>'Pro',
			'CAC'=>'His','CAT'=>'His',
			'CAA'=>'Gln','CAG'=>'Gln',
			'CGA'=>'Arg','CGC'=>'Arg','CGG'=>'Arg','CGT'=>'Arg',
			'ATA'=>'Ile','ATC'=>'Ile','ATT'=>'Ile',
			'ATG'=>'Met',
			'ACA'=>'Thr','ACC'=>'Thr','ACG'=>'Thr','ACT'=>'Thr',
			'AAC'=>'Asn','AAT'=>'Asn',
			'AAA'=>'Lys','AAG'=>'Lys',
			'AGC'=>'Ser','AGT'=>'Ser',
			'AGA'=>'Arg','AGG'=>'Arg',
			'GTA'=>'Val','GTC'=>'Val','GTG'=>'Val','GTT'=>'Val',
			'GCA'=>'Ala','GCC'=>'Ala','GCG'=>'Ala','GCT'=>'Ala',
			'GAC'=>'Asp','GAT'=>'Asp',
			'GAA'=>'Glu','GAG'=>'Glu',
			'GGA'=>'Gly','GGC'=>'Gly','GGG'=>'Gly','GGT'=>'Gly'
			);

# switches autoflush on
$| = 1;

my @input_files = @ARGV;
my $header = "Pos\tRef\tType\tAllel\tFreq";
$header .= "\tSubst\tGene\tGeneName\tProduct\tResistanceSNP";
$header .= "\tPhyloSNP\tInterestingRegion";


foreach my $input_file (@input_files) {
  my $filename = basename($input_file, ".tsv");
  #my @names = split(/_/, $filename);
  #my $fullID = join "_",$names[0], $names[1];
  #my $output_file = $fullID . "_codon_change_table.tsv";
	my $output_file = $filename . "_true_codon.tsv";
  open(OUT,">$output_file") or die $!;
  print OUT "$header\n";
  parse_input($input_file);
  close OUT;
}
exit(0);

################################################
# subroutines
################################################


sub parse_input {
  my($input_file) = @_;
  my $header="";

  open FILE, $input_file;
  my @content = (<FILE>);
  chomp @content;
  close FILE;

  shift @content;    #remove the header of the file

  my @lines = ();    #array of arrays, each input array is a line from the input file
  my @tab = ();      #one array from @lines that will be modified and printed to output
  my @holder = ();   #an array that is used for formating the Subst column output
  my @ins_allel = ();#array of alleles within an insertion
	my @ins_order = ();#array which will be used to set the correct order of alleles within an insertion
  my @del_ref = ();  #array of reference alleles which are removed by a deletion
  my $complement = 0;#Check if the gene is on the complement strand
  my $count = 0; 	 #counter for each line in the input file
  my $cc1 = 0;		 #First nucleotide in a codon
  my $cc2 = 0;		 #Second nucleotide in a codon
  my $cc3 = 0;		 #Third nucleotide in a codon
  my $no_change = 0; #Ref alleles of the codon
  my $change = 0;	 #Right side of the Subst column
  my $del_start = 0; #First position of a deletion
  my $del = 0;   	 #Deletion length counter, switch to check whether the loop is already in a deletion
  my $ins = 0;   	 #Insertion length counter, switch to check whether the loop is already in a insertion
  my $multi = 0; 	 #T/F switch for seeing whether the prior line is already a codon change
  my $aa = 0;    	 #Final amino acid after the multicodon change
  my $aa_pos = 0;	 #Left side of the Subst column
  my $tab_line = 0;  #Line to print OUT
  my $pos = 0;	     #Position column value when multicodon/ins/del
  my $ref = 0;		 #Reference -||-
  my $allel = 0;	 #Allele -||-
  my $freq = 0;      #Allele frequency
	my $covf = 0;    #Coverage forward
	my $covr = 0;    #Coverage reverse


  foreach my $line (@content) {
	chomp($line);
    # this will take care of different line break characters
    $line =~ s/\015?\012?$//;

    next unless $line;
    $lines[$count] = [split /\t/, $line];
	$count += 1;
  }
  for my $i (0 .. $#lines) {                                                                             #loop over the variant file
	if ($lines[$i][3] eq "SNP") {
	$complement = substr $lines[$i][11], -1;
	if ($complement ne "c"){
		if (!$multi) {																					 #check to see whether we are looping within a multichange codon
			my $cc = (split ' ', $lines[$i][10])[1];
			if (!$cc) {
				@tab = @{$lines[$i]};
				$tab_line = join "\t", $tab[0],@tab[2..4],$tab[8],@tab[10..16];
				print OUT "$tab_line\n";
				next;
			}
			$cc = (split '/', substr( $cc, 1, (length($cc) - 2) ))[1];
			if ((split '',$cc)[2] =~ /[A,T,C,G]/) {                                                      # if aaA codon type print it out
				@tab = @{$lines[$i]};
				$tab_line = join "\t",$tab[0],@tab[2..4],$tab[8],@tab[10..16];
				print OUT "$tab_line\n";
			} elsif ((split '',$cc)[1] =~ /[A,T,C,G]/) {
				if ($lines[$i][0] == ($lines[$i+1][0]-1)) {                                            #if aAA codon type print it out and set the counter to skip the next iteration
					$cc2 = (split ' ', $lines[$i][10])[1];
					$cc3 = (split ' ', $lines[$i+1][10])[1];
					$no_change = (split '/', substr( $cc2, 1, (length($cc2) - 2) ))[0];
					$cc1 = substr((split '/', substr( $cc2, 1, (length($cc2) - 2) ))[1],0,-2);
					$cc2 = substr((split '/', substr( $cc2, 1, (length($cc2) - 2) ))[1],1,-1);
					$cc3 = substr((split '/', substr( $cc3, 1, (length($cc3) - 2) ))[1],2);
					$cc = uc(join "",$cc1,$cc2,$cc3);
					$change = join "/",$no_change,(join "",$cc1,$cc2,$cc3);
					@holder = split "",$change; push @holder, ")"; unshift @holder, "(";
					$change = join "", @holder;
					$aa_pos = substr((split ' ', $lines[$i][10])[0],0,-3);
					$aa = $codon{$cc};
					$aa_pos = join "", $aa_pos,$aa;
					#$pos = join "-", $lines[$i][0],$lines[$i+1][0]; #start and stop
					$pos = $lines[$i][0];
					$ref = join "", $lines[$i][2],$lines[$i+1][2];
					$allel = join "", $lines[$i][4],$lines[$i+1][4];
					if ($lines[$i][5] <= $lines[$i+1][5]) {
						$covf = $lines[$i][5];
					} else {
						$covf = $lines[$i+1][5];
					}
					if ($lines[$i][6] <= $lines[$i+1][6]) {
						$covr = $lines[$i][6];
					} else {
						$covr = $lines[$i+1][6];
					}
					if ($lines[$i][8] <= $lines[$i+1][8]) {
						$freq = $lines[$i][8];
					} else {
						$freq = $lines[$i+1][8];
					}
					@tab = @{$lines[$i]};
					$tab_line = join "\t", $pos,$ref,$tab[3],$allel,$freq,(join " ",$aa_pos,$change),@tab[11..16];
					print OUT "$tab_line\n";
					$multi += 1;
				} else{																					  #if aAa codon type print it out
					@tab = @{$lines[$i]};
					$tab_line = join "\t", $tab[0],@tab[2..4],$tab[8],@tab[10..16];
					print OUT "$tab_line\n";
				}
			} else {
				if ($lines[$i][0] == ($lines[$i+1][0]-1) && $lines[$i+1][2] eq "SNP") {	                                              #AA.
					if ($lines[$i][0] == ($lines[$i+2][0]-2) && $lines[$i+2][2] eq "SNP") {                                           #if AAA codon type print it out and set the counter to skip the next iteration
						$cc1 = (split ' ', $lines[$i][10])[1];
						$cc2 = (split ' ', $lines[$i+1][10])[1];
						$cc3 = (split ' ', $lines[$i+2][10])[1];
						$no_change = (split '/', substr( $cc2, 1, (length($cc2) - 2) ))[0];
						$cc1 = substr((split '/', substr( $cc1, 1, (length($cc1) - 2) ))[1],0,-2);
						$cc2 = substr((split '/', substr( $cc2, 1, (length($cc2) - 2) ))[1],1,-1);
						$cc3 = substr((split '/', substr( $cc3, 1, (length($cc3) - 2) ))[1],2);
						$cc = join "",$cc1,$cc2,$cc3;
						$change = join "/",$no_change,(join "",$cc1,$cc2,$cc3);
						@holder = split "",$change; push @holder, ")"; unshift @holder, "(";
						$change = join "", @holder;
						$aa_pos = substr((split ' ', $lines[$i][10])[0],0,-3);
						$aa = $codon{$cc};
						$aa_pos = join "", $aa_pos,$aa;
						#$pos = join "-", $lines[$i][0],$lines[$i+2][0]; #start and stop
						$pos = $lines[$i][0];
						$ref = join "", $lines[$i][2],$lines[$i+1][2],$lines[$i+2][2];
						$allel = join "", $lines[$i][4],$lines[$i+1][4],$lines[$i+2][4];
						if (($lines[$i][5] <= $lines[$i+1][5]) and ($lines[$i][5] <= $lines[$i+2][5])) {
							$covf = $lines[$i][5];
						} elsif ($lines[$i+1][5] <= $lines[$i+2][5]) {
							$covf = $lines[$i+1][5];
						} else {
							$covf = $lines[$i+2][5];
						}
						if (($lines[$i][6] <= $lines[$i+1][6]) and ($lines[$i][6] <= $lines[$i+2][6])) {
							$covr = $lines[$i][6];
						} elsif ($lines[$i+1][6] <= $lines[$i+2][6]) {
							$covr = $lines[$i+1][6];
						} else {
							$covr = $lines[$i+2][6];
						}
						if (($lines[$i][8] <= $lines[$i+1][8]) and ($lines[$i][8] <= $lines[$i+2][8])) {
							$freq = $lines[$i][8];
						} elsif ($lines[$i+1][8] <= $lines[$i+2][8]) {
							$freq = $lines[$i+1][8];
						} else {
							$freq = $lines[$i+2][8];
						}
						@tab = @{$lines[$i]};
						$tab_line = join "\t",$pos,$ref,$tab[3],$allel,$freq,(join " ",$aa_pos,$change),@tab[11..16];
						print OUT "$tab_line\n";
						$multi += 1;
					} else {																			  #if AAa codon type print it out and set the counter to skip the next iteration
						$cc1 = (split ' ', $lines[$i][10])[1];
						$cc2 = (split ' ', $lines[$i+1][10])[1];
						$cc3 = (split ' ', $lines[$i+1][10])[1];
						$no_change = (split '/', substr( $cc2, 1, (length($cc2) - 2) ))[0];
						$cc1 = substr((split '/', substr( $cc1, 1, (length($cc1) - 2) ))[1],0,-2);
						$cc2 = substr((split '/', substr( $cc2, 1, (length($cc2) - 2) ))[1],1,-1);
						$cc3 = substr((split '/', substr( $cc3, 1, (length($cc3) - 2) ))[1],2);
						$cc = uc(join "",$cc1,$cc2,$cc3);
						$change = join "/",$no_change,(join "",$cc1,$cc2,$cc3);
						@holder = split "",$change; push @holder, ")"; unshift @holder, "(";
						$change = join "", @holder;
						$aa_pos = substr((split ' ', $lines[$i][10])[0],0,-3);
						$aa = $codon{$cc};
						$aa_pos = join "", $aa_pos,$aa;
						#$pos = join "-", $lines[$i][0],$lines[$i+1][0]; #start and stop
						$pos = $lines[$i][0];
						$ref = join "", $lines[$i][2],$lines[$i+1][2];
						$allel = join "", $lines[$i][4],$lines[$i+1][4];
						if ($lines[$i][5] <= $lines[$i+1][5]) {
							$covf = $lines[$i][5];
						} else {
							$covf = $lines[$i+1][5];
						}
						if ($lines[$i][6] <= $lines[$i+1][6]) {
							$covr = $lines[$i][6];
						} else {
							$covr = $lines[$i+1][6];
						}
						if ($lines[$i][8] <= $lines[$i+1][8]) {
							$freq = $lines[$i][8];
						} else {
							$freq = $lines[$i+1][8];
						}
						@tab = @{$lines[$i]};
						$tab_line = join "\t",$pos,$ref,$tab[3],$allel,$freq,(join " ",$aa_pos,$change),@tab[11..16];
						print OUT "$tab_line\n";
						$multi += 1;
					}
				} elsif ($lines[$i][0] == ($lines[$i+1][0]-2) && $lines[$i+1][2] eq "SNP") {                                          #if AaA codon type print it out and set the counter to skip the next iteration
						$cc1 = (split ' ', $lines[$i][10])[1];
						$cc2 = (split ' ', $lines[$i][10])[1];
						$cc3 = (split ' ', $lines[$i+1][10])[1];
						$no_change = (split '/', substr( $cc1, 1, (length($cc1) - 2) ))[0];
						$cc1 = substr((split '/', substr( $cc1, 1, (length($cc1) - 2) ))[1],0,-2);
						$cc2 = substr((split '/', substr( $cc2, 1, (length($cc2) - 2) ))[1],1,-1);
						$cc3 = substr((split '/', substr( $cc3, 1, (length($cc3) - 2) ))[1],2);
						$cc = uc(join "",$cc1,$cc2,$cc3);
						$change = join "/",$no_change,(join "",$cc1,$cc2,$cc3);
						@holder = split "",$change; push @holder, ")"; unshift @holder, "(";
						$change = join "", @holder;
						$aa_pos = substr((split ' ', $lines[$i][10])[0],0,-3);
						$aa = $codon{$cc};
						$aa_pos = join "", $aa_pos,$aa;
						#$pos = join "-", $lines[$i][0],$lines[$i+1][0]; #start and stop
						$pos = $lines[$i][0];
						$ref = join "", $lines[$i][2],$lines[$i+1][2];
						$allel = join "", $lines[$i][4],$lines[$i+1][4];
						if ($lines[$i][5] <= $lines[$i+1][5]) {
							$covf = $lines[$i][5];
						} else {
							$covf = $lines[$i+1][5];
						}
						if ($lines[$i][6] <= $lines[$i+1][6]) {
							$covr = $lines[$i][6];
						} else {
							$covr = $lines[$i+1][6];
						}
						if ($lines[$i][8] <= $lines[$i+1][8]) {
							$freq = $lines[$i][8];
						} else {
							$freq = $lines[$i+1][8];
						}
						@tab = @{$lines[$i]};
						$tab_line = join "\t",$pos,$ref,$tab[3],$allel,$freq,(join " ",$aa_pos,$change),@tab[11..16];
						print OUT "$tab_line\n";
						$multi += 1;
				} else {																				  #if Aaa codon type print it out
					@tab = @{$lines[$i]};
					$tab_line = join "\t", $tab[0],@tab[2..4],$tab[8],@tab[10..16];
					print OUT "$tab_line\n";
				}

			}
		} else {
			if (($lines[$i][0] == ($lines[$i+1][0]-1))
			and (substr((split ' ', $lines[$i][10])[0],0,-3) eq substr((split ' ', $lines[$i+1][10])[0],0,-3))) { #check if in a middle of multi change codon and skip it
				next
			} else{																								  #set the multi codon change to 0 when moving onto next one
				$multi = 0;
			}
		}
	} else {                                                                                                      #complement SNPs basically reverse the previous loops
		if (!$multi) {																					 #check to see whether we are looping within a multichange codon
			my $cc = (split ' ', $lines[$i][10])[1];
			if (!$cc) {
				@tab = @{$lines[$i]};
				$tab_line = join "\t", $tab[0],@tab[2..4],$tab[8],@tab[10..16];
				print OUT "$tab_line\n";
				next;
			}
			$cc = (split '/', substr( $cc, 1, (length($cc) - 2) ))[1];
			if ((split '',$cc)[0] =~ /[A,T,C,G]/) {                                                      # if Aaa codon type print it out
				@tab = @{$lines[$i]};
				$tab_line = join "\t",$tab[0],@tab[2..4],$tab[8],@tab[10..16];
				print OUT "$tab_line\n";
			} elsif ((split '',$cc)[1] =~ /[A,T,C,G]/) {
				if ($lines[$i][0] == ($lines[$i+1][0]-1) && $lines[$i+1][2] eq "SNP") {                                            #if AAa codon type print it out and set the counter to skip the next iteration
					$cc2 = (split ' ', $lines[$i][10])[1];
					$cc1 = (split ' ', $lines[$i+1][10])[1];
					$no_change = (split '/', substr( $cc2, 1, (length($cc2) - 2) ))[0];
					$cc3 = substr((split '/', substr( $cc2, 1, (length($cc2) - 2) ))[1],2);
					$cc1 = substr((split '/', substr( $cc1, 1, (length($cc2) - 2) ))[1],0,-2);
					$cc2 = substr((split '/', substr( $cc2, 1, (length($cc2) - 2) ))[1],1,-1);
					$cc = uc(join "",$cc1,$cc2,$cc3);
					$change = join "/",$no_change,(join "",$cc1,$cc2,$cc3);
					@holder = split "",$change; push @holder, ")"; unshift @holder, "(";
					$change = join "", @holder;
					$aa_pos = substr((split ' ', $lines[$i][10])[0],0,-3);
					$aa = $codon{$cc};
					$aa_pos = join "", $aa_pos,$aa;
					#$pos = join "-", $lines[$i][0],$lines[$i+1][0]; #start and stop
					$pos = $lines[$i][0];
					$ref = join "", $lines[$i][2],$lines[$i+1][2];
					$allel = join "", $lines[$i][4],$lines[$i+1][4];
					if ($lines[$i][5] <= $lines[$i+1][5]) {
						$covf = $lines[$i][5];
					} else {
						$covf = $lines[$i+1][5];
					}
					if ($lines[$i][6] <= $lines[$i+1][6]) {
						$covr = $lines[$i][6];
					} else {
						$covr = $lines[$i+1][6];
					}
					if ($lines[$i][8] <= $lines[$i+1][8]) {
						$freq = $lines[$i][8];
					} else {
						$freq = $lines[$i+1][8];
					}
					@tab = @{$lines[$i]};
					$tab_line = join "\t", $pos,$ref,$tab[3],$allel,$freq,(join " ",$aa_pos,$change),@tab[11..16];
					print OUT "$tab_line\n";
					$multi += 1;
				} else{																					  #if aAa codon type print it out
					@tab = @{$lines[$i]};
					$tab_line = join "\t", $tab[0],@tab[2..4],$tab[8],@tab[10..16];
					print OUT "$tab_line\n";
				}
			} else {
				if ($lines[$i][0] == ($lines[$i+1][0]-1) && $lines[$i+1][2] eq "SNP") {	                                              #.AA
					if ($lines[$i][0] == ($lines[$i+2][0]-2) && $lines[$i+2][2] eq "SNP") {                                           #if AAA codon type print it out and set the counter to skip the next iteration
						$cc3 = (split ' ', $lines[$i][10])[1];
						$cc2 = (split ' ', $lines[$i+1][10])[1];
						$cc1 = (split ' ', $lines[$i+2][10])[1];
						$no_change = (split '/', substr( $cc2, 1, (length($cc2) - 2) ))[0];
						$cc1 = substr((split '/', substr( $cc1, 1, (length($cc1) - 2) ))[1],0,-2);
						$cc2 = substr((split '/', substr( $cc2, 1, (length($cc2) - 2) ))[1],1,-1);
						$cc3 = substr((split '/', substr( $cc3, 1, (length($cc3) - 2) ))[1],2);
						$cc = join "",$cc1,$cc2,$cc3;
						$change = join "/",$no_change,(join "",$cc1,$cc2,$cc3);
						@holder = split "",$change; push @holder, ")"; unshift @holder, "(";
						$change = join "", @holder;
						$aa_pos = substr((split ' ', $lines[$i][10])[0],0,-3);
						$aa = $codon{$cc};
						$aa_pos = join "", $aa_pos,$aa;
						#$pos = join "-", $lines[$i][0],$lines[$i+2][0]; #start and stop
						$pos = $lines[$i][0];
						$ref = join "", $lines[$i][2],$lines[$i+1][2],$lines[$i+2][2];
						$allel = join "", $lines[$i][4],$lines[$i+1][4],$lines[$i+2][4];
						if (($lines[$i][5] <= $lines[$i+1][5]) and ($lines[$i][5] <= $lines[$i+2][5])) {
							$covf = $lines[$i][5];
						} elsif ($lines[$i+1][5] <= $lines[$i+2][5]) {
							$covf = $lines[$i+1][5];
						} else {
							$covf = $lines[$i+2][5];
						}
						if (($lines[$i][6] <= $lines[$i+1][6]) and ($lines[$i][6] <= $lines[$i+2][6])) {
							$covr = $lines[$i][6];
						} elsif ($lines[$i+1][6] <= $lines[$i+2][6]) {
							$covr = $lines[$i+1][6];
						} else {
							$covr = $lines[$i+2][6];
						}
						if (($lines[$i][8] <= $lines[$i+1][8]) and ($lines[$i][8] <= $lines[$i+2][8])) {
							$freq = $lines[$i][8];
						} elsif ($lines[$i+1][8] <= $lines[$i+2][8]) {
							$freq = $lines[$i+1][8];
						} else {
							$freq = $lines[$i+2][8];
						}
						@tab = @{$lines[$i]};
						$tab_line = join "\t",$pos,$ref,$tab[3],$allel,$freq,(join " ",$aa_pos,$change),@tab[11..16];
						print OUT "$tab_line\n";
						$multi += 1;
					} else {																			  #if aAA codon type print it out and set the counter to skip the next iteration
						$cc3 = (split ' ', $lines[$i][10])[1];
						$cc2 = (split ' ', $lines[$i+1][10])[1];
						$cc1 = (split ' ', $lines[$i+1][10])[1];
						$no_change = (split '/', substr( $cc2, 1, (length($cc2) - 2) ))[0];
						$cc1 = substr((split '/', substr( $cc1, 1, (length($cc1) - 2) ))[1],0,-2);
						$cc2 = substr((split '/', substr( $cc2, 1, (length($cc2) - 2) ))[1],1,-1);
						$cc3 = substr((split '/', substr( $cc3, 1, (length($cc3) - 2) ))[1],2);
						$cc = uc(join "",$cc1,$cc2,$cc3);
						$change = join "/",$no_change,(join "",$cc1,$cc2,$cc3);
						@holder = split "",$change; push @holder, ")"; unshift @holder, "(";
						$change = join "", @holder;
						$aa_pos = substr((split ' ', $lines[$i][10])[0],0,-3);
						$aa = $codon{$cc};
						$aa_pos = join "", $aa_pos,$aa;
						#$pos = join "-", $lines[$i][0],$lines[$i+1][0]; #start and stop
						$pos = $lines[$i][0];
						$ref = join "", $lines[$i][2],$lines[$i+1][2];
						$allel = join "", $lines[$i][4],$lines[$i+1][4];
						if ($lines[$i][5] <= $lines[$i+1][5]) {
							$covf = $lines[$i][5];
						} else {
							$covf = $lines[$i+1][5];
						}
						if ($lines[$i][6] <= $lines[$i+1][6]) {
							$covr = $lines[$i][6];
						} else {
							$covr = $lines[$i+1][6];
						}
						if ($lines[$i][8] <= $lines[$i+1][8]) {
							$freq = $lines[$i][8];
						} else {
							$freq = $lines[$i+1][8];
						}
						@tab = @{$lines[$i]};
						$tab_line = join "\t",$pos,$ref,$tab[3],$allel,$freq,(join " ",$aa_pos,$change),@tab[11..16];
						print OUT "$tab_line\n";
						$multi += 1;
					}
					#} elsif ($lines[$i][0] == ($lines[$i+1][0]-2) && $lines[$i+1][2] eq "SNP") {
				} elsif ($lines[$i][0] == ($lines[$i+1][0]-2) && $lines[$i+1][2] eq "SNP") {                                          #if AaA codon type print it out and set the counter to skip the next iteration
						$cc3 = (split ' ', $lines[$i][10])[1];
						$cc2 = (split ' ', $lines[$i][10])[1];
						$cc1 = (split ' ', $lines[$i+1][10])[1];
						$no_change = (split '/', substr( $cc1, 1, (length($cc1) - 2) ))[0];
						$cc1 = substr((split '/', substr( $cc1, 1, (length($cc1) - 2) ))[1],0,-2);
						$cc2 = substr((split '/', substr( $cc2, 1, (length($cc2) - 2) ))[1],1,-1);
						$cc3 = substr((split '/', substr( $cc3, 1, (length($cc3) - 2) ))[1],2);
						$cc = uc(join "",$cc1,$cc2,$cc3);
						$change = join "/",$no_change,(join "",$cc1,$cc2,$cc3);
						@holder = split "",$change; push @holder, ")"; unshift @holder, "(";
						$change = join "", @holder;
						$aa_pos = substr((split ' ', $lines[$i][10])[0],0,-3);
						$aa = $codon{$cc};
						$aa_pos = join "", $aa_pos,$aa;
						#$pos = join "-", $lines[$i][0],$lines[$i+1][0]; #start and stop
						$pos = $lines[$i][0];
						$ref = join "", $lines[$i][2],$lines[$i+1][2];
						$allel = join "", $lines[$i][4],$lines[$i+1][4];
						if ($lines[$i][5] <= $lines[$i+1][5]) {
							$covf = $lines[$i][5];
						} else {
							$covf = $lines[$i+1][5];
						}
						if ($lines[$i][6] <= $lines[$i+1][6]) {
							$covr = $lines[$i][6];
						} else {
							$covr = $lines[$i+1][6];
						}
						if ($lines[$i][8] <= $lines[$i+1][8]) {
							$freq = $lines[$i][8];
						} else {
							$freq = $lines[$i+1][8];
						}
						@tab = @{$lines[$i]};
						$tab_line = join "\t",$pos,$ref,$tab[3],$allel,$freq,(join " ",$aa_pos,$change),@tab[11..16];
						print OUT "$tab_line\n";
						$multi += 1;
				} else {																				  #if Aaa codon type print it out
					@tab = @{$lines[$i]};
					$tab_line = join "\t", $tab[0],@tab[2..4],$tab[8],@tab[10..16];
					print OUT "$tab_line\n";
				}

			}
		} else {
			if (($lines[$i][0] == ($lines[$i+1][0]-1))
			and (substr((split ' ', $lines[$i][10])[0],0,-3) eq substr((split ' ', $lines[$i+1][10])[0],0,-3))) { #check if in a middle of multi change codon and skip it
				next
			} else{																								  #set the multi codon change to 0 when moving onto next one
				$multi = 0;
			}
		}

	}







	} elsif ($lines[$i][3] eq "Ins"){																	  #check if insertion
		if ($lines[$i][0] == $lines[$i+1][0]){                                                            #check if next line is still an insertion, and grow the existing one if yes
			$ins += 1;
			if ($ins == 1) {
				$covf = $lines[$i][5];
				$covr = $lines[$i][6];
				$freq = $lines[$i][8];
			}
			push @ins_allel, $lines[$i][4];
			push @ins_order, $lines[$i][1];
			if ($covf >= $lines[$i+1][5]) {
				$covf = $lines[$i+1][5];
			}
			if ($covr >= $lines[$i+1][6]) {
				$covr = $lines[$i+1][6];
			}
			if ($freq >= $lines[$i+1][8]) {
				$freq = $lines[$i+1][8];
			}
		} else {															#if this is the last ins line add it to ones before and print it all OUT
			if ($ins == 0) {
				$covf = $lines[$i][5];
				$covr = $lines[$i][6];
				$freq = $lines[$i][8];
			}																						  #if this is the last ins line add it to ones before and print it all OUT
			push @ins_allel, $lines[$i][4];
			push @ins_order, $lines[$i][1];
			no warnings "numeric";
			my @idx = sort { $ins_order[$a] <=> $ins_order[$b] } 0 .. $#ins_order;
			@ins_order = @ins_order[@idx];
			@ins_allel = @ins_allel[@idx];
			$pos = join "+",$lines[$i][0],(scalar @ins_allel);
			$ref = $lines[$i][2];
			$allel = join "", @ins_allel;
			@tab = @{$lines[$i]};
			$tab_line = join "\t",$pos,$ref,$tab[3],$allel,$freq,@tab[10..16];
			print OUT "$tab_line\n";
			$ins = 0;
			@ins_order = ();
			@ins_allel = ();
		}
	} else {																							  #the same as above but for deletions
		$del += 1;
		if ($del == 1) {
				$del_start = $lines[$i][0];
				$covf = $lines[$i][5];
				$covr = $lines[$i][6];
				$freq = $lines[$i][8];
			}
		if ($lines[$i+1][3] eq "Del"  && ($lines[$i][0] == ($lines[$i+1][0]-1))){
			push @del_ref, $lines[$i][2];
			if ($covf >= $lines[$i+1][5]) {
				$covf = $lines[$i+1][5];
			}
			if ($covr >= $lines[$i+1][6]) {
				$covr = $lines[$i+1][6];
			}
			if ($freq >= $lines[$i+1][8]) {
				$freq = $lines[$i+1][8];
			}
		} else {
			push @del_ref, $lines[$i][2];
			#$pos = join "-", $del_start, $lines[$i][0]; #start and stop
			$pos = $del_start;
			$ref = join "", @del_ref;
			$allel = $lines[$i][4];
			@tab = @{$lines[$i]};
			$tab_line = join "\t",$pos,$ref,$tab[3],$allel,$freq,@tab[10..16];
			print OUT "$tab_line\n";
			$del = 0;
			@del_ref = ();
		}
	}
  }
}
return(1);
