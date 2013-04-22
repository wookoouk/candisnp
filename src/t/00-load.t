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


$new_r->run(qq'x=runif(1)');
my $result = $new_r->get('x');
warn "from test $result";


diag( "Testing CandiSNP $CandiSNP::VERSION, Perl $], $^X" );
