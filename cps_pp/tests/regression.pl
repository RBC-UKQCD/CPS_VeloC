# $Id: regression.pl.in,v 1.20 2004/10/27 14:36:24 zs Exp $
# tests/regression.pl.  Generated from regression.pl.in by configure.
#--------------------------------------------------------------------
 

# directory with good test data in it
$control_dir = "regressions/xmldiff/v6_0_0" ; 

#
#  Different subroutines to convert the CPS output
#  files into XML format
#

$combine_file = sub { 
    system("perl ../combine_files_xml.pl *.dat > $candidate_file"); 
};

$scalar_run = sub { 
    system("$executable  > $output_file") ; 
};


$do_nothing = sub {  };


# assume that the xmldiff utility is in your PATH
$xmldiff = "xmldiff" ; 
$status = system("which $xmldiff > /dev/null") ;
die("ERROR: can NOT find $xmldiff in your path") unless $status == 0 ;


# load in information on which tests to run
require $control_dir."/test_names.pl" ;


# machine parameters
require 'machine_depend.pl' ; 

#------ 

# This makes it run autoconf version:
  $machine = 'x86_64-unknown-linux-gnu';
  $parallel = 'yes';
  $executable = "NOARCH.x";
  $mach = "NOARCH"  ; 




$executable = $machine_names{$mach}{"exec"} ; 

print "Name of program = " . $executable  .  "\n"  ; 



#
# Configurations: These should be relevant always:
#

#The name of the shell script this will create:
$output_dir = "regressions";



#
#  initialize (eg. set up qdeamon for QCDOC)
#

$setup_machine = $machine_names{$mach}{"setup"} ; 
&$setup_machine() ; 


##
##  first compile the applications
##

print "COMPILATION STATUS\n";
foreach $test ( keys %test_names ) 
{
    chdir $test ;


    system("gmake -s -f Makefile_regression clean > log_clean 2>&1") ;
#     on qcdoc the bourne shell and the recursive make get confused 
     system("bash -c \"gmake -d -s -f Makefile_regression ".$executable." \" > log_compile 2>&1") ; 

    print $test ."  " ; 
    if( -f $executable )
    {
	print "SUCCESS\n" ; 
    }
    else
    {
	print "FAILURE\n" ; 
    }

  chdir ".." ; 

# Loop over to the next test:
}


##exit ; 
#
#  now run the code
#

print  "TEST_CASE_STATUS\n";
foreach $test ( keys %test_names ) {
    $do_test = $test_names{$test}{"check_data"} ; 
    if( $do_test == "YES" ) 
    {

#        Construct the names of the files to store the results:
	$output_file = "stdio.".$test.".log";
        $candidate_dir = $machine_names{$mach}{"output"} ; 
	$candidate_file = $candidate_dir."candidate_".$test.".xml";
	$xmldiff_log_file = "xmldiff_".$test.".log";

	$control_file   = "../".$control_dir."/control_".$test.".xml";
	$metric_file    =  "../".$control_dir."/metric_".$test.".xml";

	chdir $test ;
	if( -f $executable )
	{
	    print $test."  " ;
	    $run_exec = $machine_names{$mach}{"run"} ; 
            &$run_exec() ; 

# Grab all of the *.dat output files and sling them into a single file:
	    $xml_the_files = $test_names{$test}{"get_xml"} ; 
	    &$xml_the_files() ; 

# Diff the output and report if the output differs from that expected:

	    $cmp = system("$xmldiff $control_file $candidate_file $metric_file $xmldiff_log_file "); 
	    if( $cmp == 0 ) 
	    {
		print "PASS\n"  ;
	    } else
	    {
		print "FAIL\n"  ;
	    }

	    system("gmake -s -f Makefile_regression clean");
	    chdir ".." ; 

	} ## end if executable exists

    }

# Loop over to the next test:
}

print " ------------------------------\n";
print " DISCLAIMER\n";
print " Please also test this code on \n";
print " a physical system before using in \n";
print " production runs \n";
print " ------------------------------\n";



$shutdown_machine = $machine_names{$mach}{"final"} ; 
&$shutdown_machine() ; 


exit 0;
# the end 

