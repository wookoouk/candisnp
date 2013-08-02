#!perl -T
use Modern::Perl 2011;
use strict;
use warnings FATAL => 'all';
use Test::More qw( no_plan );
use Data::Dumper;


BEGIN {
    use_ok( 'CandiSNP' ) || print "Bail out!\n";
}
diag( "Testing CandiSNP $CandiSNP::VERSION, Perl $], $^X" );

require_ok( 'CandiSNP' );



#check file opening
#my $fh = CandiSNP::_open_file("sample_data/header.csv");
#isa_ok $fh, "Tie::Handle::CSV", "can't create CSV file object";

#check file header
#ok(CandiSNP::_header_ok($fh), "CSV file header could not be validated");

#check snp filter
#my $l = <$fh>;
#ok(CandiSNP::_is_snp($l, -cutoff => 0.7 ), "not a SNP");

#check file hash made properly given cutoff
#my $hash = {
#          'Chr1' => {
#                      '1' => {
#                               '_in_cds' => "NA",
#                               '_ctga' => "FALSE",
#                               '_ref' => 'A',
#                               '_syn' => "NA",
#                               '_allele_freq' => '0.89',
#                               '_alt' => 'T'
#                             }
#                    }
#        };

#my $data_from_file = CandiSNP::get_positions_from_file(
#	-file => "sample_data/header.csv",
#	-genome => "athalianaTair10"
#	 );
#	
#is_deeply $data_from_file, $hash, "not expected data structure for file";

##check file written correctly...
#unlink "sample_data/test_out.csv";
#CandiSNP::_data_hash_to_file($data_from_file, "sample_data/test_out.csv");
#ok -e "sample_data/test_out.csv", "";

##test snpEff is installed in bin
#ok -e CandiSNP::bin_folder() . "/snpEff.jar", "No snpEff.jar in bin/";
#ok -e CandiSNP::bin_folder() . "/snpEff.config", "No snpEff.config in bin/";
#ok -e CandiSNP::bin_folder() . "/genomes", "No bin/genomes/ folder";

#diag("Running external snpEff .. may take some time");
#my $large_data = CandiSNP::get_positions_from_file(
#	-file=> "sample_data/snps.csv", 
#	-genome => "athalianaTair10"
#	);
#$large_data = CandiSNP::annotate_positions($large_data, -genome => "athalianaTair10");

##check get no snps for unknown chromosome

#my $bad_data = CandiSNP::get_positions_from_file(
#	-file=> "sample_data/no_chr.csv", 
#	-genome => "athalianaTair10"
#	);

#is $bad_data, 0, "got SNPs when shouldnt have";

#my $annot_snp = {
#						'_in_cds' => "FALSE",
#					    '_ctga' => "TRUE",
#					    '_ref' => 'G',
#					    '_syn' => "FALSE",
#					    '_allele_freq' => '1.0',
#					    '_alt' => 'A',
#						'_gene' => 'AT1G29750'
#        };

##pick a snp to check its annotations .. 
#my $snp_selected = $$large_data{'Chr1'}{10417334};

#is_deeply $snp_selected, $annot_snp, "not expected data structure for parsed SNP";

##check filtering removes SNPs
#my $filtered_data = CandiSNP::apply_filter($large_data, "0.99");
#ok( scalar keys %{$$large_data{"Chr1"}} > scalar keys %{$$filtered_data{"Chr1"}}, "filtering didn't remove any SNPs");

#get information on genome length
#my $genome_lengths = CandiSNP::genome_lengths('athalianaTair10');
#ok 34964571 eq $genome_lengths->{Chr1}, "genome length not retrieved";

#returns list of six positions into 
#my ($scale_marks, $scale_labels) = CandiSNP::scale_marks($genome_lengths);

#my @exp_marks = ("10000000", "20000000", "30000000");
#my @exp_labels = ("10Mb", "20Mb", "30Mb");

#is_deeply $scale_marks, \@exp_marks, "didn't get expected divisions";
#is_deeply $scale_labels, \@exp_labels, "didn't get expected labels";


#start R session ..
#my $new_r = CandiSNP::R;
#isa_ok $new_r, "Statistics::R", "can't create R connection object";

##add the filtered data to the R session, along with the scales for these data, plot image
#my $image = CandiSNP::plot_data($new_r, $filtered_data, $scale_marks, $scale_labels);


##how to run the thing...
my $big_data = CandiSNP::get_positions_from_file(
-file => "sample_data/large_snp_no_chr.csv",
-genome => "athalianaTair10"
);

#warn Dumper $big_data;
##once data is loaded, you can get a unique id for it ...
my $filetag = CandiSNP::get_filetag($big_data);


$big_data = CandiSNP::annotate_positions($big_data, -genome => "athalianaTair10", -filetag => $filetag);	
#warn Dumper $big_data;
#my $snp_selected = $$big_data{'Chr1'}{10417334};

#warn Dumper $snp_selected;

my $filtered_big_data = CandiSNP::apply_filter($big_data, "0.9");
#warn Dumper $filtered_big_data;
my $all_genome_lengths = CandiSNP::genome_lengths('athalianaTair10');
my ($big_scale_marks, $big_scale_labels) = CandiSNP::scale_marks($all_genome_lengths);
#warn Dumper $big_scale_marks;
#warn Dumper $big_scale_labels;
my $R = CandiSNP::R;


my $palette = CandiSNP::get_palette('contrast');
#warn Dumper $palette;
#warn Dumper $filtered_big_data;
#warn Dumper $filetag;
#warn Dumper $big_scale_marks;
#warn Dumper $big_scale_labels;
#warn Dumper $all_genome_lengths;
my $big_image = CandiSNP::plot_data($R, $filtered_big_data, $filetag, $big_scale_marks, $big_scale_labels,$all_genome_lengths,$palette);
warn Dumper "$big_image printed";
##dump file of used snps to a csv
CandiSNP::data_hash_to_file($filtered_big_data,$filetag,-format=>'long');
#$R->stop;