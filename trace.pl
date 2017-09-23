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
#  June20, edit: bug fixes in sequence read-in parser
#
#  Fixed support for "!" extra base pairs.
#  Added strand-path diagram as extra output at the end



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

while ( $line = <FILE> ) {
    $l = length $line;
    if ($j<1){$maxlength = $l;}
    $k = 0;
    for ($i=0;$i<$l;$i++) {
        $cols = substr("$line", $i, 1);
        if ( $cols eq "#" ) { $i=$i+$maxlength; }
        if ( $cols eq ">" ) { $i=$i+$maxlength; }
        if ( $cols eq "@" ) { $i=$i+$maxlength; }
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
        if ( $cols eq "!"    ) { $m[$j][$k] = "!"; $k++; }         
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
                   # print "\nThe 5p end is found at row $r, column $c and is running $d.\n"; 
                }
            }
        }
    }
}

# Find total amount of nucleotides in blueprint

my $nt = 0;
for ($i=0;$i<1000;$i++) {
    for ($j=0;$j<1000;$j++) {
        if ( defined $m[$i][$j] ) {
            if ($m[$i][$j] =~ /[NXGACURYKMSWVHBDT]/ ) {
                $nt++;
            }
        }
    } 
}
#print "There are $nt nucleotides in the blueprint file.\n";

# Now trace the structure

my $r2 = scalar ($r);
my $c2 = scalar ($c);
my $d2 = $d;
my $num = 0;
my @seq = ( );
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
   # print "The structure has been successfully traced (3p end found).\n";
} else {
   # print "The trace through the structure failed (3p end not found). Ended at row $r, column $c.\n";
}
# print "The are " . scalar (@seq) . " nts from 5p to 3p.\n";

# Find base pairs

$r = scalar ($r2); # reset row
$c = scalar ($c2); # reset column
$d = $d2; # restore 5p direction
my @a = ( );
my @b = ( );
my @p = ( );
for ($k=0;$k<$nt+10000;$k++) { 
    # TRACE HORIZONTAL STRAND
    if ( $d eq "right" && $m[$r][$c+1] =~ /[xNXGACURYKMSWVHBDT-]/ ) { 
        if ( $m[$r][$c+1] =~ /[NXGACURYKMSWVHBDT]/ ) { 
            if ( $m[$r+1][$c+1] =~ /[!p\*]/ ) {
                push @a, $n[$r][$c+1];
                push @b, $n[$r+2][$c+1]; 
                push @p, $m[$r+1][$c+1]; 
            } elsif ( $m[$r-1][$c+1] =~ /[!p\*]/ ) {
                push @a, $n[$r][$c+1];
                push @b, $n[$r-2][$c+1]; 
                push @p, $m[$r-1][$c+1]; 
            } else {
                push @a, $n[$r][$c+1];
                push @b, 0; 
                push @p, "-";
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
                push @a, $n[$r][$c-1];
                push @b, 0; 
                push @p, "-";
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


# print blueprint

if ( $file2 ne "no" ) {
   # print "\n\nOutput: 2D diagram with sequence\n";
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

