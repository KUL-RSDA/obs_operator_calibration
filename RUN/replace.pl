#!/usr/bin/perl
#
($arg1, $arg2, $arg3, $arg4) = @ARGV;
# arg1 is old filename
# arg2 is new filename
# arg3 is old string
# arg4 is new string

$num_args = @ARGV;

if( $num_args < 4 ) {
   print "Please use: replace.pl old_filename new_filename, old_string,
          new_string.  \n ";
   exit(0);
}

$input_file = $arg1;
$output_file = $arg2;
$old_string = $arg3;
$new_string = $arg4;

open(INPUTFILE, "$input_file");
open(OUTPUTFILE, ">$output_file");

while (<INPUTFILE>) {
#      chop($_);  
#      print "This line is: $_. \n";
      s/$old_string/$new_string/g;

      print OUTPUTFILE "$_";

}

close(INPUTFILE);
close(OUTPUTFILE);

