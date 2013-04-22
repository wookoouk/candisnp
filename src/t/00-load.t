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

#check snp filter
my $l = <$fh>;
ok(CandiSNP::_is_snp($l, -cutoff => 0.7 ), "not a SNP");

#check file hash made properly given cutoff
my $hash = {
          'Chr1' => {
                      '1' => {
                               '_in_cds' => "NA",
                               '_ctga' => "FALSE",
                               '_ref' => 'A',
                               '_syn' => "NA",
                               '_allele_freq' => '0.89',
                               '_alt' => 'T'
                             }
                    }
        };

my $data_from_file = CandiSNP::get_positions_from_file(
	-file => "sample_data/header.csv", 
	-cutoff => 0.7,
	-genome => "athalianaTair10"
	 );
	
is_deeply $hash, $data_from_file, "not expected data structure for file";

##check file written correctly...
unlink "sample_data/test_out.csv";
CandiSNP::_data_hash_to_file($data_from_file, "sample_data/test_out.csv");
ok -e "sample_data/test_out.csv";

warn Dumper "Running external snpEff .. may take some time\n";
CandiSNP::get_positions_from_file(
	-file=> "sample_data/snps.csv", 
	-cutoff=> 0.7,
	-genome => "athalianaTair10"
	);

diag( "Testing CandiSNP $CandiSNP::VERSION, Perl $], $^X" );

