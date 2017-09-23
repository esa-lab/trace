#!/usr/bin/env perl   

#  -*- perl -*-

# RNA-trace
# written by Ebbe S. Andersen <esa@inano.au.dk>, 2013
#
# Updated by Cody Geary, April 2015
#  Added support for NGACURYKMSWVHBDT
#  Added support for -,|,+ in plain ascii, which then converts to special characters
#  Removed space-padding requirement at the end of most lines (only the first line must be padded enough)
#
#  Converts T to U
#
#  June20, 2015 edit: bug fixes in sequence read-in parser
#
#  Fixed support for "!" extra base pairs.
#  Added strand-path diagram as extra output at the end
#
#  June30 2017 - Update to generate target.txt outputs  use: perl trace_pattern.pl pattern.txt > target.txt



use strict;
use warnings FATAL => qw ( all );
use Data::Dumper;
use Encode;
use utf8;


# >>>>>>>>>>>>>>>>>> RUN PROGRAM <<<<<<<<<<<<<<<<<<<<

my ( $file1, $file2, $line, @cols, $seq );

( $file1, $file2 ) = @ARGV;

my @pri = ( );
if ( defined $file2 ) {
    open FILE, "< $file2";
    while ( $line = <FILE> ) {
        print "$line";
        print "\n";
        @pri = split(//, $line);
    }
} else {
    $file2 = "no";
}

if ( not open FILE, "< $file1" ) { die "file not found!"; }

# translate unicode to symbols:
# ----------------------------------------
# symbol unicode       unique description
# ----------------------------------------
# _      \302\240      240    space       
# -      \342\224\200  200    straight    
# i      \342\224\202  202    up-down
# p      \342\224\212  212    pair        
# x      \342\224\274  274    cross       
# L      \342\225\255  255    down-right
# J      \342\225\256  256    down-left
# 7      \342\225\257  257    up-left
# r      \342\225\260  260    up-right
# b                           base-pair vertical
# ----------------------------------------

my @m = ( );
my @n = ( );
my @t = ( );
@cols = ( );
my $i = 0;
my $j = 0;
my $cols = 0;
my $l = 0;
my $k = 0;
my $maxlength = 0;
my $endpadding = 0;
my $numtees = 0;
my $name = 'Untitled';
my $KL_pattern = '';
my $nt=0;

#vars for the pattern_checker subroutine
my @map = ();
my @repeat_map =(); 
my $wobble_zones=0;
my $wobble_window=10;
my $pattern_repeats=0;
my $poly_repeats=0;
my $restriction_sites =0;
my @trial_sol = ();


my @a = ( );
my @b = ( );
my @p = ( );
my @seq = ( );
my $print_pattern = 0;

	my $find = "";
	my $replace = "";


	while ( $line = <FILE> ) {
		$l = length $line;
		
		$find = "I";
		$replace = "!";
		$line =~ s/$find/$replace/g;   #replace all the "I" with "!" -- I dunno why this didn't work just searching for $cols eq "I" below...

		if ($j<1){$maxlength = $l;}
		$k = 0;
		for ($i=0;$i<$l;$i++) {
			$cols = substr("$line", $i, 1);
			if ( $cols eq "#" ) { $i=$i+$maxlength; }
			 if ( $cols eq ">" ) { $name = substr("$line", $i+1, 32); $i=$i+$maxlength; } #Extract the File Name Here with a max 32 chars, skip to the next line
			 if ( $cols eq "@" ) { $KL_pattern = substr("$line", $i+1, 32);  $i=$i+$maxlength; } #Extract the KL matching pattern with a max 32 chars, skip to the next line
			if ( $cols eq " " ) { $m[$j][$k] = " "; $k++; }     
			if ( $cols eq "*" ) { $m[$j][$k] = "*"; $k++; }        
			if ( $cols eq "\257" ) { $m[$j][$k] = "J"; $k++; }
			if ( $cols eq "\/" &&  substr("$line", $i+1, 1)!~ /[+xNXGACURYKMSWVHBDT35-]/ ) { $m[$j][$k] = "J"; $k++; }     
			if ( $cols eq "\200" ) { $m[$j][$k] = "-"; $k++; }
			if ( $cols eq "-" ) { $m[$j][$k] = "-"; $k++; }     
			if ( $cols eq "\260" ) { $m[$j][$k] = "L"; $k++; }
			if ( $cols eq "\\" &&  substr("$line", $i+1, 1)=~ /[+xNXGACURYKMSWVHBDT35-]/ ) { $m[$j][$k] = "L"; $k++; }      
			if ( $cols eq "\255" ) { $m[$j][$k] = "r"; $k++; }
			if ( $cols eq "\/" &&  substr("$line", $i+1, 1)=~ /[+xNXGACURYKMSWVHBDT35-]/ ) { $m[$j][$k] = "r"; $k++; }     
			if ( $cols eq "\256" ) { $m[$j][$k] = "7"; $k++; } 
			if ( $cols eq "\\" &&  substr("$line", $i+1, 1)!~ /[+xNXGACURYKMSWVHBDT35-]/ ) { $m[$j][$k] = "7"; $k++; } 
			if ( $cols eq "\202" ) { $m[$j][$k] = "i"; $k++; }
			if ( $cols eq "|" ) { $m[$j][$k] = "i"; $k++; }   
			if ( $cols eq "^" ) { $m[$j][$k] = "i"; $k++; }     
			if ( $cols eq "\212" ) { $m[$j][$k] = "p"; $k++; } 
			if ( $cols eq ":" ) { $m[$j][$k] = "p"; $k++; } 
			if ( $cols eq "\;" ) { $m[$j][$k] = "p"; $k++; }    
			if ( $cols eq "!" ) { $m[$j][$k] = "!"; $k++; }          
			if ( $cols eq "\274" ) { $m[$j][$k] = "x"; $k++; }
			if ( $cols eq "+" ) { $m[$j][$k] = "x"; $k++; }
			if ( $cols eq "="    ) { $m[$j][$k] = "b"; $k++; }     
			if ( $cols =~ /\w/ ) { $m[$j][$k] = "$cols"; $k++; }
		}
		if ($k<$maxlength){$endpadding=$maxlength-$k;}

		for ($i=0;$i<$endpadding;$i++) {
			$m[$j][$k] = " "; $k++;
		}
		$m[$j][$k] = "\n"; $k++;

		$j++;
	}


	#scrub filename and KL-pattern inputs
	$find = " ";
	$replace = "";
	$name =~ s/$find/$replace/g;
	$KL_pattern =~ s/$find/$replace/g;
	my @KL_pat = split(//,$KL_pattern);
	my $current_KL=0;


	# Find 5 prime end

	my $r = 0;
	my $c = 0;
	my $d = "left";
	for ($i=0;$i<1000;$i++) {
		for ($j=0;$j<1000;$j++) {
			if ( defined $m[$i][$j] ) { 
				if ($m[$i][$j] =~ /\d+/ ) {
					if ( scalar $m[$i][$j] && $m[$i][$j] == 5 ) { 
						$r = $i;
						$c = $j;
						if ( $m[$r][$c+1] =~ /[NXGACURYKMSWVHBDT-]/ ) { $d = "right"; }
						if ( $m[$r][$c-1] =~ /[NXGACURYKMSWVHBDT-]/ ) { $d = "left"; }
						if ( $m[$r+1][$c] =~ /[NXGACURYKMSWVHBDTi]/ ) { $d = "down"; }
						if ( $m[$r-1][$c] =~ /[NXGACURYKMSWVHBDTi]/ ) { $d = "up"; }
						#print "\nThe 5p end is found at row $r, column $c and is running $d.\n"; 
					}
				}
			}
		}
	}

	# Find number of nucleotides in blueprint

	$nt = 0;
	for ($i=0;$i<1000;$i++) {
		for ($j=0;$j<1000;$j++) {
			if ( defined $m[$i][$j] ) {
				if ($m[$i][$j] =~ /[NXGACURYKMSWVHBDT]/ ) {
					$nt++;
				}
			}
		} 
	}


	#Trace the structure

	my $r2 = scalar ($r);
	my $c2 = scalar ($c);
	my $d2 = $d;
	my $num = 0;
	@seq = ( );
	my $test = "test";
	for ($k=0;$k<$nt+10000;$k++) { 
		# TRACE HORIZONTAL STRAND
		if ( $file2 eq "no" ) {
			if ( $d eq "right" && $m[$r][$c+1] =~ /[xNXGACURYKMSWVHBDT-]/ ) { 
				if ( $m[$r][$c+1] =~ /[NXGACURYKMSWVHBDT]/ ) { $num++; $n[$r][$c+1] = $num; push @seq, $m[$r][$c+1]; } 
				$c++;
			}
			if ( $d eq "left"  && $m[$r][$c-1] =~ /[xNXGACURYKMSWVHBDT-]/ ) { 
				if ( $m[$r][$c-1] =~ /[NXGACURYKMSWVHBDT]/ ) { $num++; $n[$r][$c-1] = $num; push @seq, $m[$r][$c-1]; }
				$c--;
			}

			if ( $d eq "up" && $m[$r-1][$c] =~ /[xNXGACURYKMSWVHBDTi]/ ) { 
				if ( $m[$r-1][$c] =~ /[NXGACURYKMSWVHBDT]/ ) { $num++; $n[$r-1][$c] = $num; push @seq, $m[$r-1][$c]; } 
				$r--; 
			}
			if ( $d eq "down" && $m[$r+1][$c] =~ /[xNXGACURYKMSWVHBDTi]/ ) { 
				if ( $m[$r+1][$c] =~ /[NXGACURYKMSWVHBDT]/ ) { $num++; $n[$r+1][$c] = $num; push @seq, $m[$r+1][$c]; } 
				$r++; 
			}



		} else { # when a sequence is available as file2 then we add it on the m-grid here
			if ( $d eq "right" && $m[$r][$c+1] =~ /[xNXGACURYKMSWVHBDT-]/ ) { 
				if ( $m[$r][$c+1] =~ /[NXGACURYKMSWVHBDT]/ ) { $m[$r][$c+1] = $pri[$num]; $num++; $n[$r][$c+1] = $num; }
				$c++;
			}
			if ( $d eq "left"  && $m[$r][$c-1] =~ /[xNXGACURYKMSWVHBDT-]/ ) { 
				if ( $m[$r][$c-1] =~ /[NXGACURYKMSWVHBDT]/ ) { $m[$r][$c-1] = $pri[$num]; $num++; $n[$r][$c-1] = $num; }
				$c--; 
			}


			if ( $d eq "up" && $m[$r-1][$c] =~ /[xNXGACURYKMSWVHBDTi]/ ) { 
				if ( $m[$r-1][$c] =~ /[NXGACURYKMSWVHBDT]/ ) { $m[$r-1][$c] = $pri[$num]; $num++; $n[$r-1][$c] = $num; }
				$r--;
			}
			if ( $d eq "down" && $m[$r+1][$c] =~ /[xNXGACURYKMSWVHBDTi]/ ) { 
				if ( $m[$r+1][$c] =~ /[NXGACURYKMSWVHBDT]/ ) { $m[$r+1][$c] = $pri[$num]; $num++; $n[$r+1][$c] = $num; }
				$r++; 
			}


		}
		# CROSS-OVER
		if ( $d eq "right" ) {
			if ( $m[$r][$c+1] eq "7" ) { $d = "down"; $c++;}
			if ( $m[$r][$c+1] eq "J" ) { $d = "up"; $c++;}
		}    
		if ( $d eq "left" ) {
			if ( $m[$r][$c-1] eq "L" ) { $d = "up"; $c--;}
			if ( $m[$r][$c-1] eq "r" ) { $d = "down";  $c--;}
		}

		if ( $d eq "down" ) {
			if ( $m[$r+1][$c] eq "L" ) { $d = "right"; $r++;}
			if ( $m[$r+1][$c] eq "J" ) { $d = "left"; $r++;}
		}    
		if ( $d eq "up" ) {
			if ( $m[$r-1][$c] eq "7" ) { $d = "left"; $r--;}
			if ( $m[$r-1][$c] eq "r" ) { $d = "right"; $r--;}
		}


		if ( $m[$r][$c+1] =~ /\d/ ) { if ( $m[$r][$c+1] == 3 ) { $test = "success"; last; } }
		if ( $m[$r][$c-1] =~ /\d/ ) { if ( $m[$r][$c-1] == 3 ) { $test = "success"; last; } }

		if ( $m[$r+1][$c] =~ /\d/ ) { if ( $m[$r+1][$c] == 3 ) { $test = "success"; last; } }
		if ( $m[$r-1][$c] =~ /\d/ ) { if ( $m[$r-1][$c] == 3 ) { $test = "success"; last; } }


	}
	if ( $test eq "success" ) {
		#print "The structure has been successfully traced (3p end found).\n";
	} else {
		print "The trace through the structure failed (3p end not found). Ended at row $r, column $c.\n";
	}
	#print "The are " . scalar (@seq) . " nts from 5p to 3p.\n";

	# Find base pairs

	$r = scalar ($r2); # reset row
	$c = scalar ($c2); # reset column
	$d = $d2; # restore 5p direction
	@a = ( );
	@b = ( );
	@p = ( );
	for ($k=0;$k<$nt+10000;$k++) { 
		# TRACE HORIZONTAL STRAND
		if ( $d eq "right" && $m[$r][$c+1] =~ /[xNXGACURYKMSWVHBDT-]/ ) { 
			if ( $m[$r][$c+1] =~ /[NXGACURYKMSWVHBDT]/ ) { 
				if ( $m[$r+1][$c+1] =~ /[!p\*]/ ) {  #check above for pairs
					push @a, $n[$r][$c+1];
					push @b, $n[$r+2][$c+1]; 
					push @p, $m[$r+1][$c+1]; 
				} elsif ( $m[$r-1][$c+1] =~ /[!p\*]/ ) {  #check below for pairs
					push @a, $n[$r][$c+1];
					push @b, $n[$r-2][$c+1]; 
					push @p, $m[$r-1][$c+1]; 
				} else {
					if ( $m[$r][$c+1] =~ /[X]/) {
						push @a, $n[$r][$c+1];
						push @b, 0;
						push @p, $KL_pat[$current_KL];
						if ($m[$r][$c+2] !~ /[X]/) {$current_KL++;}
					} else {
						push @a, $n[$r][$c+1];
						push @b, 0; 
						push @p, "-";
					}               
				}
			}
			$c++;
		} 
	
	
	
		if ( $d eq "left"  && $m[$r][$c-1] =~ /[xNXGACURYKMSWVHBDT-]/ ) { 
			if ( $m[$r][$c-1] =~ /[NXGACURYKMSWVHBDT]/ ) { 
				if ( $m[$r+1][$c-1] =~ /[!p\*]/ ) {
					push @a, $n[$r][$c-1];
					push @b, $n[$r+2][$c-1]; 
					push @p, $m[$r+1][$c-1]; 
				} elsif ( $m[$r-1][$c-1] =~ /[!p\*]/ ) {
					push @a, $n[$r][$c-1];
					push @b, $n[$r-2][$c-1]; 
					push @p, $m[$r-1][$c-1]; 
				} else {
					if ( $m[$r][$c-1] =~ /[X]/) {
						push @a, $n[$r][$c-1];
						push @b, 0;
						push @p, $KL_pat[$current_KL];
						if ($m[$r][$c-2] !~ /[X]/) {$current_KL++;}
					} else {
						push @a, $n[$r][$c-1];
						push @b, 0; 
						push @p, "-";
					}               
				}
			}
			$c--; 
		}

		if ( $d eq "down" && $m[$r+1][$c] =~ /[xNXGACURYKMSWVHBDTi]/ ) { 
			if ( $m[$r+1][$c] =~ /[NXGACURYKMSWVHBDT]/ ) { 
				if ( $m[$r+1][$c+1] =~ /[b\*]/ ) {
					push @a, $n[$r+1][$c];
					push @b, $n[$r+1][$c+2]; 
					push @p, $m[$r+1][$c+1];  
				} elsif ( $m[$r+1][$c-1] =~ /[b\*]/ ) {
					push @a, $n[$r+1][$c];
					push @b, $n[$r+1][$c-2]; 
					push @p, $m[$r+1][$c-1]; 
				} else {
					push @a, $n[$r+1][$c];
					push @b, 0; 
					push @p, "i";
				}
			}
			$r++;
		}

		if ( $d eq "up" && $m[$r-1][$c] =~ /[xNXGACURYKMSWVHBDTi]/ ) { 
			if ( $m[$r-1][$c] =~ /[NXGACURYKMSWVHBDT]/ ) { 
				if ( $m[$r-1][$c+1] =~ /[b\*]/ ) {
					push @a, $n[$r-1][$c];
					push @b, $n[$r-1][$c+2]; 
					push @p, $m[$r-1][$c+1]; 
				} elsif ( $m[$r-1][$c-1] =~ /[b\*]/ ) {
					push @a, $n[$r-1][$c];
					push @b, $n[$r-1][$c-2]; 
					push @p, $m[$r-1][$c-1]; 
				} else {
					push @a, $n[$r-1][$c];
					push @b, 0; 
					push @p, "i";
				}
			}
			$r--;
		}


		if ( $d eq "right" && $m[$r][$c+1] =~ /[x-]/ ) { $c++; }
		if ( $d eq "left"  && $m[$r][$c-1] =~ /[x-]/ ) { $c--; }

		if ( $d eq "down" && $m[$r+1][$c] =~ /[xi]/ ) { $r++; }
		if ( $d eq "up" && $m[$r-1][$c] =~ /[xi]/ ) { $r--; }




		# CROSS-OVER

		if ( $d eq "right" ) {
			if ( $m[$r][$c+1] eq "7" ) { $d = "down"; $c++;}
			if ( $m[$r][$c+1] eq "J" ) { $d = "up"; $c++;}
		}    
		if ( $d eq "left" ) {
			if ( $m[$r][$c-1] eq "L" ) { $d = "up"; $c--;}
			if ( $m[$r][$c-1] eq "r" ) { $d = "down";  $c--;}
		}

		if ( $d eq "down" ) {
			if ( $m[$r+1][$c] eq "L" ) { $d = "right"; $r++;}
			if ( $m[$r+1][$c] eq "J" ) { $d = "left"; $r++;}
		}    
		if ( $d eq "up" ) {
			if ( $m[$r-1][$c] eq "7" ) { $d = "left"; $r--;}
			if ( $m[$r-1][$c] eq "r" ) { $d = "right"; $r--;}
		}
	}




# print filename
print "$name\n";

# print structure

my $a = 0;
my $b = 0;
$i = 0;

my $structure_map = "";

foreach $a ( @a ) {
    if ( defined $b[$i] && defined $a ) {
        if ( $a > $b[$i] and $b[$i] != 0 ) { 
            if ( $p[$i] eq "p" ) { print ")"; $structure_map = $structure_map.")"; } 
            if ( $p[$i] eq "b" ) { print ")"; $structure_map = $structure_map.")";} 
            if ( $p[$i] eq "!" ) { print "}"; $structure_map = $structure_map.")";} 
            if ( $p[$i] eq "*" ) { print "]"; $structure_map = $structure_map."]";}
        }
        if ( $a < $b[$i] and $b[$i] != 0 ) { 
            if ( $p[$i] eq "p" ) { print "("; $structure_map = $structure_map."(";} 
            if ( $p[$i] eq "b" ) { print "("; $structure_map = $structure_map."(";} 
            if ( $p[$i] eq "!" ) { print "{"; $structure_map = $structure_map."(";} 
            if ( $p[$i] eq "*" ) { print "["; $structure_map = $structure_map."[";} 
        }
        if ( $b[$i] == 0 ) {
        	if ( $p[$i] =~ /[ABCDEFGHIJ1234567890]/ ) { print $p[$i]; $structure_map = $structure_map.".";} else { print "."; $structure_map = $structure_map.".";}
        }
    } else {
        if ( $p[$i] =~ /[ABCDEFGHIJ1234567890]/ ) { print $p[$i]; $structure_map = $structure_map.".";} else { print "."; $structure_map = $structure_map.".";}
    } 
    $i++;
}
print "\n";

my $sequence = "";
# print sequence
foreach $seq ( @seq ) {
	if ( $seq =~ /[NGACURYKMSWVHBD]/ ) { print "$seq"; $sequence=$sequence.$seq;}
	if ( $seq eq "T" ) { print "U"; $sequence = $sequence."U";}
	if ( $seq eq "X" ) { print "N"; $sequence = $sequence."N";}
}
print "\n\n\n";

print "My Structure map:  \n";
print "$structure_map \n\n";

&map($structure_map);  #calculate the map function


# print "\n\n2D diagram with sequence\n";
for ($i=0;$i<1000;$i++) {
	for ($j=0;$j<1000;$j++) {
		if ( defined $m[$i][$j] ) {
			if ( $m[$i][$j] =~ /[NXGACURYKMSWVHBXD35\*]/ ) { print "$m[$i][$j]"; } 
			if ( $m[$i][$j] eq "T" ) { print "U"; } 
			if ( $m[$i][$j] eq "\n" ) { print "$m[$i][$j]"; }
			if ( $m[$i][$j] eq " " ) { print "$m[$i][$j]"; } 
			if ( $m[$i][$j] eq "7" ) { print "\342\225\256"; }
			if ( $m[$i][$j] eq "-" ) { print "\342\224\200"; }
			if ( $m[$i][$j] eq "r" ) { print "\342\225\255"; }
			if ( $m[$i][$j] eq "L" ) { print "\342\225\260"; }
			if ( $m[$i][$j] eq "J" ) { print "\342\225\257"; }
			if ( $m[$i][$j] eq "i" ) { print "\342\224\202"; }
			if ( $m[$i][$j] eq "p" ) { print "\342\224\212"; }
			if ( $m[$i][$j] eq "b" ) { print "="; }
			if ( $m[$i][$j] eq "!" ) { print "I"; }           
			if ( $m[$i][$j] eq "x" ) { print "\342\224\274"; }
		}
	}
}

print "\n\nStrand Path\n";
    for ($i=0;$i<1000;$i++) {
        for ($j=0;$j<1000;$j++) {
            if ( defined $m[$i][$j] ) { 
                if ( $m[$i][$j] =~ /[35]/ ) { print "$m[$i][$j]"; } 
                if ( $m[$i][$j] =~ /[NXGACURYKMSWVHBD]/ ) { print "\342\224\200"; } 
                if ( $m[$i][$j] eq "T" ) { print "\342\224\200"; } 
                if ( $m[$i][$j] eq "\n" ) { print "$m[$i][$j]"; }
                if ( $m[$i][$j] eq " " ) { print "$m[$i][$j]"; } 
                if ( $m[$i][$j] eq "7" ) { print "\342\225\256"; }
                if ( $m[$i][$j] eq "-" ) { print "\342\224\200"; }
                if ( $m[$i][$j] eq "r" ) { print "\342\225\255"; }
                if ( $m[$i][$j] eq "L" ) { print "\342\225\260"; }
                if ( $m[$i][$j] eq "J" ) { print "\342\225\257"; }
                if ( $m[$i][$j] eq "i" ) { print "\342\224\202"; }
                if ( $m[$i][$j] eq "p" ) { print "\342\224\202"; }
                if ( $m[$i][$j] eq "b" ) { print "\342\224\202"; }
                if ( $m[$i][$j] eq "!" ) { print "\342\224\202"; } 
                if ( $m[$i][$j] eq "\*" ) { print "\342\224\202"; }           
                if ( $m[$i][$j] eq "x" ) { print "\342\224\274"; }
            }
        }
    }


my @sequence_array = ();
@sequence_array = split(//,$sequence);
my $sequence_check = 1;
@trial_sol=@sequence_array;

for ($i=0; $i<$nt; $i++){   #$nt is the length of the strand, measured earlier
	if ( $sequence_array[$i] ne "A" && $sequence_array[$i] ne "U" && $sequence_array[$i] ne "C" && $sequence_array[$i] ne "G" ){$sequence_check=0;}
}



if ($sequence_check==0){
	print"Sequence not suitable for pattern search...\n";
} else {
	print "Sequence good for pattern search..\n";
	
	&countrepeats;	
	
	print "\n\nHighlighting Repeat Sequences \n";
	
	print " 8S|8W Repeats: $pattern_repeats	  5G|5C|5A|5U Repeats: $poly_repeats   Restriction: $restriction_sites   Complements: $wobble_zones \n";
	

	my $printer_counter = 0;
    for ($i=0;$i<1000;$i++) {
        for ($j=0;$j<1000;$j++) {
            if ( defined $m[$i][$j] ) { 
                if ( $m[$i][$j] =~ /[35]/ ) { print "$m[$i][$j]"; } 
                if ( $m[$i][$j] =~ /[NXACGURYKMSWVHBDTPWS]/ ) {
                	if ($repeat_map[$n[$i][$j]-1] eq "-") {
                		print "\342\227\246";
                	} else {
                		print "$repeat_map[$n[$i][$j]-1]";
                	}
                }                 
                
                if ( $m[$i][$j] eq "\n" ) { print "$m[$i][$j]"; }
                if ( $m[$i][$j] eq " " ) { print "$m[$i][$j]"; } 
                if ( $m[$i][$j] eq "7" ) { print "\342\225\256"; }
                if ( $m[$i][$j] eq "-" ) { print "\342\224\200"; }
                if ( $m[$i][$j] eq "r" ) { print "\342\225\255"; }
                if ( $m[$i][$j] eq "L" ) { print "\342\225\260"; }
                if ( $m[$i][$j] eq "J" ) { print "\342\225\257"; }
                if ( $m[$i][$j] eq "i" ) { print "\342\224\202"; }
                if ( $m[$i][$j] eq "p" ) { print "\342\224\202"; }
                if ( $m[$i][$j] eq "b" ) { print "\342\224\202"; }
                if ( $m[$i][$j] eq "!" ) { print "\342\224\202"; } 
                if ( $m[$i][$j] eq "\*" ) { print "\342\224\202"; }           
                if ( $m[$i][$j] eq "x" ) { print "\342\224\274"; }
            }
        }
    }	
}


my @wobbles_seq = ();
for ($i=0; $i<$nt; $i++){
	if ( ($sequence_array[$i] eq "G" && $sequence_array[$map[$i]] eq "U") || ($sequence_array[$i] eq "U" && $sequence_array[$map[$i]] eq "G") || ($sequence_array[$i] eq "K" && $sequence_array[$map[$i]] eq "K") ){ 
		$wobbles_seq[$i] = $sequence_array[$i];
		$wobbles_seq[$map[$i]] = $sequence_array[$map[$i]];
	} else {	$wobbles_seq[$i]="\342\227\246"}
}


print "\n\nHighlighting GU Pairs\n";
    for ($i=0;$i<1000;$i++) {
        for ($j=0;$j<1000;$j++) {
            if ( defined $m[$i][$j] ) { 
                if ( $m[$i][$j] =~ /[35]/ ) { print "$m[$i][$j]"; } 
                if ( $m[$i][$j] =~ /[NXACGURYKMSWVHBD]/ ) { print "$wobbles_seq[$n[$i][$j]-1]"; }                 
                
                if ( $m[$i][$j] eq "T" ) { print "\342\224\200"; } 
                if ( $m[$i][$j] eq "\n" ) { print "$m[$i][$j]"; }
                if ( $m[$i][$j] eq " " ) { print "$m[$i][$j]"; } 
                if ( $m[$i][$j] eq "7" ) { print "\342\225\256"; }
                if ( $m[$i][$j] eq "-" ) { print "\342\224\200"; }
                if ( $m[$i][$j] eq "r" ) { print "\342\225\255"; }
                if ( $m[$i][$j] eq "L" ) { print "\342\225\260"; }
                if ( $m[$i][$j] eq "J" ) { print "\342\225\257"; }
                if ( $m[$i][$j] eq "i" ) { print "\342\224\202"; }
                if ( $m[$i][$j] eq "p" ) { print "\342\224\202"; }
                if ( $m[$i][$j] eq "b" ) { print "\342\224\202"; }
                if ( $m[$i][$j] eq "!" ) { print "\342\224\202"; } 
                if ( $m[$i][$j] eq "\*" ) { print "\342\224\202"; }           
                if ( $m[$i][$j] eq "x" ) { print "\342\224\274"; }
            }
        }
    }


print "\n\nHighlighting Structural Barriers\n\n";

my @struc_array = split(//,$structure_map);

my @barriers = ();

my $topo_count =0;

for ($i=0; $i<$nt; $i++){
	$barriers[$i]="E";  #error code,
	if ($struc_array[$i]eq"."){$barriers[$i]="\342\227\246"; $topo_count=0;} #unpaired are not marked
	if ($struc_array[$i]eq"("){$barriers[$i]="x"; $topo_count=0;} #mark open pairs as potential barriers
	if ($struc_array[$i]eq")"){	
		if ($barriers[$map[$i]] eq "x"){  #closing pairs remove barriers
			$barriers[$i]="\342\227\246"; 
			$barriers[$map[$i]]="\342\227\246"; 
			$topo_count=0;
		} 
		if ($barriers[$map[$i]] eq "X"){ #if the closing pair is already blocked, then we mark the barrier reached and start counting
			if($topo_count>5){
				$barriers[$i]="X"; 
				$barriers[$map[$i]]="~";
			} else { $barriers[$i]="~"; $barriers[$map[$i]]="~";}
			$topo_count++;  #increment the topology counter
		}  
	}
	if ($struc_array[$i]eq"]"){   #closing KLs marks all the intervening sequence is blocked
		$barriers[$i]="\342\227\246";  
		$barriers[$map[$i]]="\342\227\246";
		for($k=$map[$i]; $k<$i; $k++){
			if ($barriers[$k]eq"x"){$barriers[$k]="X"; } #mark all the blocked residues
		}

	}
}

    for ($i=0;$i<1000;$i++) {
        for ($j=0;$j<1000;$j++) {
            if ( defined $m[$i][$j] ) { 
                if ( $m[$i][$j] =~ /[35]/ ) { print "$m[$i][$j]"; } 
                if ( $m[$i][$j] =~ /[NXACGURYKMSTWVHBD]/ ) { print "$barriers[$n[$i][$j]-1]"; }                 

                if ( $m[$i][$j] eq "\n" ) { print "$m[$i][$j]"; }
                if ( $m[$i][$j] eq " " ) { print "$m[$i][$j]"; } 
                if ( $m[$i][$j] eq "7" ) { print "\342\225\256"; }
                if ( $m[$i][$j] eq "-" ) { print "\342\224\200"; }
                if ( $m[$i][$j] eq "r" ) { print "\342\225\255"; }
                if ( $m[$i][$j] eq "L" ) { print "\342\225\260"; }
                if ( $m[$i][$j] eq "J" ) { print "\342\225\257"; }
                if ( $m[$i][$j] eq "i" ) { print "\342\224\202"; }
                if ( $m[$i][$j] eq "p" ) { print "\342\224\202"; }
                if ( $m[$i][$j] eq "b" ) { print "\342\224\202"; }
                if ( $m[$i][$j] eq "!" ) { print "\342\224\202"; } 
                if ( $m[$i][$j] eq "\*" ) { print "\342\224\202"; }           
                if ( $m[$i][$j] eq "x" ) { print "\342\224\274"; }
            }
        }
    }




sub countrepeats {
	my $strand_length = $nt;
	for ($i=0; $i<$strand_length; $i++){$repeat_map[$i]="-";}  ##zero the map
	
    my $trial_sol_seq = join "",@trial_sol;
    my $pattern_sense = "AAAAAAA";
    my $pattern_antisense = "UUUUUUU";
    
    $wobble_zones=0;
	for ( $i=0; $i<($strand_length-$wobble_window); $i++){		
		
		$pattern_sense = substr($trial_sol_seq, $i, $wobble_window);  #grab the substring pattern here
		$pattern_antisense = "";  #reset the pattern here
		for (my $j=0; $j<$wobble_window; $j++){ #generate the complement sequence
			if (substr($pattern_sense,$wobble_window-$j,1) eq "A" ){$pattern_antisense=$pattern_antisense."U";}
			if (substr($pattern_sense,$wobble_window-$j,1) eq "U" ){$pattern_antisense=$pattern_antisense."A";}
			if (substr($pattern_sense,$wobble_window-$j,1) eq "C" ){$pattern_antisense=$pattern_antisense."G";}
			if (substr($pattern_sense,$wobble_window-$j,1) eq "G" ){$pattern_antisense=$pattern_antisense."C";}									
		}
		
		my $test_result = index($trial_sol_seq, $pattern_antisense, 0);
		for (my $j=$i; $j<$strand_length-$wobble_window; $j++){
			$test_result = index($trial_sol_seq, $pattern_antisense, $j);
			if ($test_result != -1){
				for (my $k=1; $k<$wobble_window; $k++){ ##bugfix try $k=1 rather than 0..
					$repeat_map[$i+$k]="P";
					$repeat_map[$test_result-1+$k]="P";
				}
			}
		}
	}
	
	for ($i=0; $i<$strand_length; $i++){
		if ($repeat_map[$i]eq"P"){$wobble_zones+=1;}
	}
	
	$pattern_repeats=0;
	for ( $i=0; $i<$strand_length-7; $i++){	   #count up stretches of A/U or G/C that are 8nts
		if ( ($trial_sol[$i]eq"A" || $trial_sol[$i]eq"U") &&
		 ($trial_sol[$i+1]eq"A" || $trial_sol[$i+1]eq"U") &&
		 ($trial_sol[$i+2]eq"A" || $trial_sol[$i+2]eq"U") &&
		 ($trial_sol[$i+3]eq"A" || $trial_sol[$i+3]eq"U") &&
		 ($trial_sol[$i+4]eq"A" || $trial_sol[$i+4]eq"U") &&
		 ($trial_sol[$i+5]eq"A" || $trial_sol[$i+5]eq"U") &&
		 ($trial_sol[$i+6]eq"A" || $trial_sol[$i+6]eq"U") &&
		 ($trial_sol[$i+7]eq"A" || $trial_sol[$i+7]eq"U") ) {  #all stretch of GC_rich are targeted
			$pattern_repeats+=1;
			$repeat_map[$i]="W"; $repeat_map[$i+1]="W";
			$repeat_map[$i+2]="W"; $repeat_map[$i+3]="W";
			$repeat_map[$i+4]="W"; $repeat_map[$i+5]="W";	   
			$repeat_map[$i+6]="W"; $repeat_map[$i+7]="W";			  
		}
		if ( ($trial_sol[$i]eq"C" || $trial_sol[$i]eq"G") &&
		 ($trial_sol[$i+1]eq"C" || $trial_sol[$i+1]eq"G") &&
		 ($trial_sol[$i+2]eq"C" || $trial_sol[$i+2]eq"G") &&
		 ($trial_sol[$i+3]eq"C" || $trial_sol[$i+3]eq"G") &&
		 ($trial_sol[$i+4]eq"C" || $trial_sol[$i+4]eq"G") &&
		 ($trial_sol[$i+5]eq"C" || $trial_sol[$i+5]eq"G") &&
		 ($trial_sol[$i+6]eq"C" || $trial_sol[$i+6]eq"G") &&
		 ($trial_sol[$i+7]eq"C" || $trial_sol[$i+7]eq"G") ) {  #all stretch of GC_rich are targeted
			$pattern_repeats+=1;
			$repeat_map[$i]="S"; $repeat_map[$i+1]="S";
			$repeat_map[$i+2]="S"; $repeat_map[$i+3]="S";
			$repeat_map[$i+4]="S"; $repeat_map[$i+5]="S";	   
			$repeat_map[$i+6]="S"; $repeat_map[$i+7]="S";			  
		} 
	}

	$poly_repeats=0;
	
		for ( $i=0; $i<$strand_length-6; $i++){	   
		if ( ($trial_sol[$i]eq$trial_sol[$i+1]) &&
			 ($trial_sol[$i]eq$trial_sol[$i+2]) &&
			 ($trial_sol[$i]eq$trial_sol[$i+3]) &&
			 ($trial_sol[$i]eq$trial_sol[$i+4])) {  #all stretch of 5G,5C,%U,5A strand are counted
			$poly_repeats+=1;
			$repeat_map[$i]=$trial_sol[$i]; $repeat_map[$i+1]=$trial_sol[$i+1];
			$repeat_map[$i+2]=$trial_sol[$i+2]; $repeat_map[$i+3]=$trial_sol[$i+3];
			$repeat_map[$i+4]=$trial_sol[$i+4];				
		}
	}


	$restriction_sites =0;
		
	for ( $i=0; $i<$strand_length-6; $i++){	   
		if 	(	( ($trial_sol[$i]eq"G") && ($trial_sol[$i+1]eq"G") && ($trial_sol[$i+2]eq"U") && ($trial_sol[$i+3]eq"C") && ($trial_sol[$i+4]eq"U") && ($trial_sol[$i+5]eq"C") ) ||    ## Bsa1 GGUCUC
				( ($trial_sol[$i]eq"G") && ($trial_sol[$i+1]eq"A") && ($trial_sol[$i+2]eq"G") && ($trial_sol[$i+3]eq"A") && ($trial_sol[$i+4]eq"C") && ($trial_sol[$i+5]eq"C") ) ||    ## Bsa1 GAGACC
				( ($trial_sol[$i]eq"G") && ($trial_sol[$i+1]eq"A") && ($trial_sol[$i+2]eq"A") && ($trial_sol[$i+3]eq"G") && ($trial_sol[$i+4]eq"A") && ($trial_sol[$i+5]eq"C") ) ||    ## BBs1 GAAGAC
				( ($trial_sol[$i]eq"G") && ($trial_sol[$i+1]eq"U") && ($trial_sol[$i+2]eq"C") && ($trial_sol[$i+3]eq"U") && ($trial_sol[$i+4]eq"U") && ($trial_sol[$i+5]eq"C") ) ||    ## BBs1 GUCUUC
				( ($trial_sol[$i]eq"C") && ($trial_sol[$i+1]eq"G") && ($trial_sol[$i+2]eq"U") && ($trial_sol[$i+3]eq"C") && ($trial_sol[$i+4]eq"U") && ($trial_sol[$i+5]eq"C") ) ||    ## BsMB1 CGUCUC
				( ($trial_sol[$i]eq"G") && ($trial_sol[$i+1]eq"A") && ($trial_sol[$i+2]eq"G") && ($trial_sol[$i+3]eq"A") && ($trial_sol[$i+4]eq"C") && ($trial_sol[$i+5]eq"G") ) ||    ## BsMB1 GAGACG
				( ($trial_sol[$i]eq"G") && ($trial_sol[$i+1]eq"C") && ($trial_sol[$i+2]eq"U") && ($trial_sol[$i+3]eq"C") && ($trial_sol[$i+4]eq"U") && ($trial_sol[$i+5]eq"U") && ($trial_sol[$i+6]eq"C")) ||   ## Sap1 GCUCUUC
				( ($trial_sol[$i]eq"G") && ($trial_sol[$i+1]eq"A") && ($trial_sol[$i+2]eq"A") && ($trial_sol[$i+3]eq"G") && ($trial_sol[$i+4]eq"A") && ($trial_sol[$i+5]eq"G") && ($trial_sol[$i+6]eq"C")) ||   ## Sap1 GAAGAGC
				( ($trial_sol[$i]eq"A") && ($trial_sol[$i+1]eq"U") && ($trial_sol[$i+2]eq"C") && ($trial_sol[$i+3]eq"U") && ($trial_sol[$i+4]eq"G") && ($trial_sol[$i+5]eq"U") && ($trial_sol[$i+6]eq"U")) ) {    ## PTH  AUCUGUU

			$restriction_sites+=1;
			$repeat_map[$i]="X"; $repeat_map[$i+1]="X";
			$repeat_map[$i+2]="X"; $repeat_map[$i+3]="X";
			$repeat_map[$i+4]="X"; $repeat_map[$i+5]="X";	   		  
		}
	}
}


####################################

sub map {
	@map = ();
	my @secondary = split('', $_[0]);
	my @last_bracket;	my @last_brace;	my @last_curlybrace;  #markers to remember the 	
	my @last_a_brace;	my @last_b_brace;	my @last_c_brace;
	
	my $counter = 0;
	foreach(@secondary){
		if ($_ eq "."){
			$map[$counter] = $counter;	 ##Single strands and map to themselves.
		}
		elsif ($_ eq "("){
			push @last_bracket, $counter;		   
		}
		elsif ($_ eq ")"){
			my $partner = pop @last_bracket;
			$map[$counter]=$partner;
			$map[$partner]=$counter;
		}
		elsif ($_ eq "["){
			push @last_brace, $counter;	
			#$map[$counter] = $counter;	 ##Single strands and map to themselves.
	 
		}
		elsif ($_ eq "]"){
			my $buddy = pop @last_brace;
			$map[$counter]=$buddy;
			$map[$buddy]=$counter;
			#$map[$counter] = $counter;	 ##Single strands and map to themselves.
		}
		$counter ++;
	}
	
	######## Scrub the sequence and fill all non-AUCG space with a random N
	
	$counter=0;
	foreach(@seq){
		if($seq[$counter] eq 'A' || $seq[$counter] eq 'U' || $seq[$counter] eq 'C' || $seq[$counter] eq 'G'){
			
		} else {
			my $pick = int rand(4);			
			if ($pick == 0){$seq[$counter]= 'A'; $seq[$map[$counter]]='U'}
			if ($pick == 1){$seq[$counter]= 'U'; $seq[$map[$counter]]='A'}
			if ($pick == 2){$seq[$counter]= 'G'; $seq[$map[$counter]]='C'}
			if ($pick == 3){$seq[$counter]= 'C'; $seq[$map[$counter]]='G'}
		}
		$counter ++;
	}
	
}
