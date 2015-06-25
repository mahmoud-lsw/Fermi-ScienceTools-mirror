#!/usr/bin/perl -w
#$Id: ST-unit-test.pl,v 1.23 2012/11/21 23:40:36 areustle Exp $
#This script is used to run all the unit tests in the FSSC distribution
#of the ScienceTools and validate the output when possible.

use strict;
use Sys::Hostname;

#This only imports Test::More if script is called directly
#Needed to use the main part of this in other scripts.
unless(caller(0)){
    require Test::More;
    Test::More->import qw(no_plan);
}

sub test_all(){
    my @test_list=@_;
    my $return;
    my @scratch;
    foreach (@test_list){
	#Some tests are not built in some cases (for example when using 
	#--without-root).  These test should be skipped and not count as
	#failures.  If a test has whitespace anywhere in the name, always
	#run it since those confuse my simple which test for existence.
	if($_=~m/\s/){
	    $return=0;
	}else{
	    $return=system("which $_ > /dev/null 2>&1");
	}
      SKIP: {
	  #No idea why I need to use the form:
	  #skip(...) conditional; 
	  #instead of:
	  #skip ... conditional;
	  #like the documentation says.
	  skip("$_ was not built.",1) if ($return);
	  $return=system("$_ > /dev/null 2>&1");
	  ok ( $return == 0, $_);
      }
    }
    return;
}

sub unit_tests{
    my $fermiroot=$_[0];

    unless (-d "output"){
	mkdir("output") or die "Cannot make directory for output files!";
    }
    chdir("output");

    my @testcodes=("test_Likelihood",
		   "test_astro",
		   "test_burstFit",
		   "test_catalogAccess",
		   #"test_GRB",
		   "test_GRBobs",
		   "test_GRBtemplate",
		   "test_eblAtten",
		   "test_genericSources",
		   "test_microQuasar",
		   "test_celestialSources",
		   "test_dataSubselector",
		   "runPy",
		   "testUtil",
		   #"test_env",
		   "test_time",
		   "test_fitsGen",
		   "test_flux",
		   "test_healpix",
		   "test_dc1Response",
		   "test_dc1aResponse",
		   "test_dc2Response",
		   "test_g25Response",
		   "test_handoffResponse",
		   "test_irfInterface",
		   "test_irfLoader",
		   "test_latResponse",
		   "test_testResponse",
		   "test_map_tools",
		   "test_observationSim",
		   "test_optimizers",
		   "test_psearch",
		   "test_pulsarDb",
		   "test_pulsePhase",
		   "test_rspgen",
		   "test_evtbin",
		   "test_st_app",
		   "export STTEST=sttest;test_st_facilities",
		   "test_st_graph",
		   "test_st_stream",
		   "test_timeSystem",
		   "test_tip",
		   "entity_test",
		   #"test_IFile",
		   #"test_altSchema",
		   "test_mem",
		   "test_write",
		   "python -m pyLikelihood",
		   "gtirfs",
		   );

    &test_all(@testcodes);

#Check test output (usually using ftdiff)
    my $return;
#evtbin
    #tolerance set for GBMLC1.lc approved by peachey (package author)
    my $excludes='';
    my @test_output=("CM2.fits",
		     "GBMLC1.lc exclude='DATASUM' tolerance=1e-16",
		     "GBMPHA1.pha",
		     "LC1.lc",
		     "merged_spectrum.pha",
		     "PHA1.pha",
		     "PHA2.pha",
		     "separate_spectrum.pha",
		     "test.ccube",);
    my @scratch;

    foreach (@test_output){
	#if datafile has a special tolerance set individually, use that
	    if ($_=~"tolerance"){
		@scratch=split " ",$_;
		$_=$scratch[0];
		$return=system("ftdiff $scratch[1] $scratch[2] $_ $fermiroot/refdata/fermi/evtbin/outref/$_ > /dev/null 2>&1");
	    }else{
		$return=system("ftdiff $excludes $_ $fermiroot/refdata/fermi/evtbin/outref/$_ > /dev/null 2>&1");
	    }
	    ok ( $return == 0, "evtbin $_");
	}

#rspgen
    #tolerances were approved by peachey (pacakge author)
    $excludes='';
    @test_output=("test_response1.rsp exclude='DATASUM,DRM_NUM,TLMIN4' tolerance=1e-12",
		  "test_response2.rsp exclude='DATASUM,DRM_NUM,TLMIN4' tolerance=1e-10",
		  "test_response3.rsp exclude='DATASUM,DRM_NUM,TLMIN4' tolerance=1e-12",
		  #"test_response4.rsp exclude='DATASUM,DRM_NUM,TLMIN4' tolerance=1e-11",
		  "test_response5.rsp exclude='DATASUM,DRM_NUM,TLMIN4' tolerance=1e-11",
		  #"test_response6.rsp exclude='DATASUM,DRM_NUM,TLMIN4' tolerance=1e-11",
		  #"test_response7.rsp exclude='DATASUM,DRM_NUM,TLMIN4' tolerance=1e-13",
		  #"test_response8.rsp exclude='DATASUM,DRM_NUM,TLMIN4' tolerance=1e-12"
		  );

    foreach (@test_output){
	#if datafile has a special tolerance set individually, use that
	    if ($_=~"tolerance"){
		@scratch=split " ",$_;
		$_=$scratch[0];
		$return=system("ftdiff $scratch[1] $scratch[2] $_ $fermiroot/refdata/fermi/rspgen/outref/$_ > /dev/null 2>&1");
	    }else{
		$return=system("ftdiff $excludes $_ $fermiroot/refdata/fermi/rspgen/outref/$_ > /dev/null 2>&1");
	    }
	    ok ( $return == 0, "rspgen $_");
	}
    #In some cases we can't use ftdiff since changes in columns or things
    #like that will always show up as different.  ftstat should show if 
    #two fits files are almost the same.
    @test_output=("test_response4.rsp",
		     "test_response6.rsp",
		     "test_response7.rsp",
		     "test_response8.rsp",
		  );

    foreach (@test_output){
	system("ftstat $_ > stat1");
	system("ftstat $fermiroot/refdata/fermi/rspgen/outref/$_ > stat2");
	$return=system("diff stat1 stat2 -I '======' > /dev/null 2>&1");
	ok ( $return == 0, "rspgen $_");
	}


#pulsePhase
    #Masaharu.Hirayama (pulsePhase author) said tolerance as high as
    #1e-3 might even be okay but I'm going to try to be more stringent than
    #that if possible.
    $excludes='exclude="DATASUM" tolerance=1e-8';
    @test_output=("testOrbitalPhaseApp_par1.fits",
		  "testOrbitalPhaseApp_par2.log",
		  "testOrbitalPhaseApp_par3.log",
		  "testOrbitalPhaseApp_par4.log",
		  "testPulsePhaseApp_par10.log",
		  "testPulsePhaseApp_par1.fits",
		  "testPulsePhaseApp_par2.fits",
		  "testPulsePhaseApp_par3.fits exclude='DATASUM' tolerance=1e-5",
		  "testPulsePhaseApp_par4.fits",
		  "testPulsePhaseApp_par5.log",
		  "testPulsePhaseApp_par6.log",
		  "testPulsePhaseApp_par7.log",
		  "testPulsePhaseApp_par8.log",
		  "testPulsePhaseApp_par9.log",);

    foreach (@test_output){
	if ($_=~".fits"){
	    #if datafile has a special tolerance set individually, use that
	    if ($_=~"tolerance"){
		@scratch=split " ",$_;
		$_=$scratch[0];
		$return=system("ftdiff $scratch[1] $scratch[2] $_ $fermiroot/refdata/fermi/pulsePhase/outref/$_ > /dev/null 2>&1");
	    }else{
		$return=system("ftdiff $excludes $_ $fermiroot/refdata/fermi/pulsePhase/outref/$_ > /dev/null 2>&1");
	    }
	}else{
	    $return=system("diff $_ $fermiroot/refdata/fermi/pulsePhase/outref/$_ > /dev/null 2>&1");
	}
	ok ( $return == 0, "pulsePhase $_");
    }

#timeSystem
    #tolerance as high as 1e-6 okayed by Masa
    $excludes='exclude="DATASUM,CREATOR" tolerance=1e-7';
    @test_output=("testTimeCorrectorApp_par1.fits",
		  "testTimeCorrectorApp_par2.fits");
    foreach (@test_output){
	if ($_=~".fits"){
	    $return=system("ftdiff $excludes $_ $fermiroot/refdata/fermi/timeSystem/outref/$_ > /dev/null 2>&1");
	}else{
	    $return=system("diff $_ $fermiroot/refdata/fermi/timeSystem/outref/$_ > /dev/null 2>&1");
	}
	ok ( $return == 0, "timeSystem $_");
    }

#st_stream
    $return=system("diff test_st_stream-out $fermiroot/refdata/fermi/st_stream/outref/test_st_stream-out > /dev/null 2>&1");
    ok ( $return == 0, "st_stream");

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
    my $return=system("ftversion >/dev/null 2>&1");
    ok ( $return == 0, "Do ftools work?") or BAIL_OUT("Must be able to run ftools!");
    &unit_tests($fermiroot);

    my $stop=time;
    my $runtime=$stop-$start;
    print "Total run time: $runtime s on $host\n";
}
################################################################
1;
