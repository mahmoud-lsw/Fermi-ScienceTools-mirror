#!/usr/bin/perl -w
#$Id: ST-pulsar-test.pl,v 1.11 2012/11/21 23:40:36 areustle Exp $
#This script is used to run and validate a specific analysis path using
#the FSSC distribution of the ScienceTools.
use strict;
use Sys::Hostname;

#This only imports Test::More if script is called directly
#Needed to use the main part of this in other scripts.
unless(caller(0)){
    require Test::More;
    Test::More->import qw(no_plan);
}

sub test_all3(){
    my @test_list=@_;
    my $return;
    my @scratch;
    foreach (@test_list){
	$return=system("$_ > /dev/null 2>&1");
	@scratch=split(/\//);
	ok ( $return == 0, $scratch[-1]);
    }
    return;
}

sub pulsar_test{
    my $fermiroot=$_[0];
    my $binpath="$fermiroot/bin";
    my $datapath="$fermiroot/refdata/fermi/test-scripts";
    unless (-d "output"){
	mkdir("output") or die "Cannot make directory for output files!";
    }

    my @testcodes=("$binpath/ST-gtbary.pl",
		   "$binpath/ST-gtephem.pl",
		   "$binpath/ST-gtpulsardb.pl",
		   "$binpath/ST-gtpphase.pl",
		   "$binpath/ST-gtophase.pl",
		   "$binpath/ST-gtpsearch.pl",
		   "$binpath/ST-gtpspec.pl",
		   "$binpath/ST-gtptest.pl",
		   );

    &test_all3(@testcodes);

#Check test output (usually using ftdiff)
    my $return;
#Compare files to reference.
    my %testcommands=("testTimeCorrectorApp_par1.fits", 'ftdiff exclude="CREATOR,DATASUM" tolerance=1e-6',
		      "ST-gtephem.log", "$binpath/ephdiff.pl",
		      "testPulsarDbApp_par2.fits", 'ftdiff exclude="CREATOR,FILENAME,DATASUM"',
		      "testPulsePhaseApp_par1.fits", 'ftdiff tolerance=1e-3 exclude="DATASUM"',
		      "testOrbitalPhaseApp_par1.fits", 'ftdiff tolerance=1e-12 exclude="DATASUM"',
		      "testPeriodSearchApp_par1.fits", 'ftdiff exclude="CREATOR,DATASUM" tolerance=1e-9',
		      "testPowerSpectrumApp_par1.fits", 'ftdiff exclude="CREATOR,DATASUM" tolerance=1e-9',
		      "testPeriodicityTestApp_par1.fits", 'ftdiff exclude="CREATOR,DATASUM" tolerance=1e-9',
		      );

    foreach (keys(%testcommands)) {
	$return=system("$testcommands{$_} $datapath/outref/$_ output/$_  > /dev/null 2>&1");
	ok ( $return == 0, $_);
    }

    return;
}

######This is the standalone part to run the test###############
################################################################
unless(caller(0)){
    my $host = hostname;
    my $start=time;

#Just test that ScienceTools and FTools are set up.
    my $fermiroot=$ENV{'FERMI_DIR'};
    ok ($fermiroot, "Is top level environment variable set?") or BAIL_OUT("\$FERMI_DIR must be set!");
    ok (-d $fermiroot, "Does \$FERMI_DIR exist and is it a directory?") or BAIL_OUT("\$FERMI_DIR must point to a directory!");
#But first be sure ftools are present (need this for analysis tests too).
    my $headas=$ENV{'HEADAS'};
    ok ($headas, "Are ftools available?") or BAIL_OUT("ftools are required!");
    my $return=system("ftversion > /dev/null 2>&1");
    ok ( $return == 0, "Do ftools work?") or BAIL_OUT("Must be able to run ftools!");
    &pulsar_test($fermiroot);

    my $stop=time;
    my $runtime=$stop-$start;
    print "Total run time: $runtime s on $host\n";
}
################################################################
1;
