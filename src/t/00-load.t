#!perl -T
use 5.006;
use strict;
use warnings FATAL => 'all';
use Test::More qw( no_plan );
use Data::Dumper;

BEGIN {
    use_ok( 'CandiSNP' ) || print "Bail out!\n";
}

require_ok( 'CandiSNP' );

my $new_r = CandiSNP::R;
isa_ok $new_r, "Statistics::R", "can't create R connection object";

#check file opening
my $fh = CandiSNP::_open_file("sample_data/header.csv");
isa_ok $fh, "Tie::Handle::CSV", "can't create CSV file object";

#check file header
ok(CandiSNP::_header_ok($fh), "CSV file header could not be validated");

#check positions retrieved properly
ok(CandiSNP::_is_snp($fh, -file => "meh", -cutoff => 0.7 ), "");


diag( "Testing CandiSNP $CandiSNP::VERSION, Perl $], $^X" );

