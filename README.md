# Introduction
Tracer is a Perl script that can take RNA motif design files from an ASCII format, and convert it to a format that does stuff..
# Usage
## Installing
Make sure you have a Perl compiler installed for your platform. Run the command `perl -version` in the command line to check for installed versions. If no installed versions are found, get an installer for you platform from https://www.perl.org/get.html.

Get the source files by running:
```
git pull git@gitlab.au.dk:esa-lab/trace.git
```
## Running
Enter the source folder, and trace a design file. The script is run with the command:
```
perl trace.pl inputfile.txt > outputfile.txt
```
Some example design files can be found in the "examples" folder. Trace one of the example designs with the command:
```
perl trace.pl examples/6hbmono.txt > result.txt
```
The file "result.txt" should now be created by the script, and contain the result of the trace.
