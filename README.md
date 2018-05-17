# Introduction
Trace is a Perl script that can trace a 2D blueprint of a single-stranded RNA origami structure, and output a code that can be used to design sequences in NUPACK (http://www.nupack.org/design/new) that fold the desired RNA origami. The script was developed to help design the single-stranded RNA origami structures published in Science in 2014 (http://science.sciencemag.org/content/345/6198/799.full). A protocol describing the details of the RNA origami design process has also been published (https://link.springer.com/protocol/10.1007%2F978-1-4939-6454-3_5).
# Usage
## Installing
Make sure you have a Perl compiler installed for your platform. Run the command `perl -version` in the command line to check for installed versions. If no installed versions are found, get an installer for you platform from https://www.perl.org/get.html.

Get the source files by running:
```
git pull git@gitlab.au.dk:esa-lab/trace.git
```
## Running
Enter the source folder and trace a 2D blueprint. The script is run with the command:
```
perl trace.pl inputfile.txt > outputfile.txt
```
Example files from the protocol mentioned in the introduction can be found in the "examples" folder. Trace the file "7_2H-AE_raw.txt" with the command:
```
perl trace.pl 7_2H-AE_raw.txt > 8_2H-AE_nupack.txt
```
The script should create the file "8_2H-AE_nupack.txt" that contain the code to be used for design in NUPACK.
