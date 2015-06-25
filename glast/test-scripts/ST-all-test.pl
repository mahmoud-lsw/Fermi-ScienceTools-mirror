#!/usr/bin/perl -w
#$Id: ST-all-test.pl,v 1.6 2012/11/21 23:40:36 areustle Exp $
#Run all FSSC ScienceTools tests

use strict;
use Sys::Hostname;
use Test::More qw(no_plan);

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

require "$fermiroot/bin/ST-unit-test.pl";
&unit_tests($fermiroot);
require "$fermiroot/bin/ST-1pl-test.pl";
&agn_test($fermiroot);
require "$fermiroot/bin/ST-pulsar-test.pl";
&pulsar_test($fermiroot);

my $stop=time;
my $runtime=$stop-$start;
print "Total run time: $runtime s on $host\n";
