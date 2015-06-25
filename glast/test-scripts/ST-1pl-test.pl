#!/usr/bin/perl -w
#$Id: ST-1pl-test.pl,v 1.17 2012/11/21 23:40:36 areustle Exp $
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

sub test_all2(){
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

sub agn_test{
    my $fermiroot=$_[0];
    my $binpath="$fermiroot/bin";
    my $datapath="$fermiroot/refdata/fermi/test-scripts";
    unless (-d "output"){
	mkdir("output") or die "Cannot make directory for output files!";
    }

    # The first command below ensures that steps after the first step
    # all start with the same data. This is needed (for now) since the
    # first (ST-mkobssim.pl) step involves random number generation,
    # which can vary between 32-bit and 64-bit systems.
    my @testcodes=(
		   "$binpath/ST-mkobssim.pl",
		   "mv ./output/_3c84_sim_events_0000.fits ./output/_3c84_sim_events_0000.fits.test;" .
		   "mv ./output/_3c84_sim_srcIds.txt ./output/_3c84_sim_srcIds.txt.test;" .
		   "cp $fermiroot/refdata/fermi/test-scripts/outref/_3c84_sim_events_0000.fits ./output/;" .
		   "cp $fermiroot/refdata/fermi/test-scripts/outref/_3c84_sim_srcIds.txt ./output/",
		   "$binpath/ST-mkdatacuts.pl",
		   "$binpath/ST-mkcmap.pl",
		   "$binpath/ST-mkexpcube.pl",
		   "$binpath/ST-mkexpmap.pl",
		   "$binpath/ST-mkgtdiffrsp.pl",
		   "$binpath/ST-gtlike.pl",
		   );

    &test_all2(@testcodes);

#Check test output (usually using ftdiff)
    my $return;
#Compare files to reference.
    my %testcommands = (
			'counts_spectra.fits' =>
			'ftdiff exclude="CREATOR,DATASUM" reltol=1e-5',
			'_3c84_ExpCube.fits' =>
			'ftdiff exclude="CREATOR,DATASUM,DATE-END"',
			'_3c84_ExpMap.fits' =>
			'ftdiff exclude="CREATOR,DATASUM" reltol=1e-6',
			'_3c84_sim_cmap.fits' =>
			'ftdiff exclude="CREATOR,DATASUM,DATE-END"',
			'_3c84_sim_events_0000.fits' =>
			'ftdiff exclude="CREATOR,DATASUM"',
			'_3c84_sim_ph02.fits' =>
			'ftdiff exclude="CREATOR,DATASUM,DATE-END" reltol=1e-5',
			'_3c84_sim_srcIds.txt' =>
			'diff',
			'results.dat' =>
			'diff',
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
    &agn_test($fermiroot);

    my $stop=time;
    my $runtime=$stop-$start;
    print "Total run time: $runtime s on $host\n";
}
################################################################
1;
