#!/usr/bin/perl

use strict;
use warnings;
use lib "../src/lib/";
use CandiSnp;

sub doWork()
{
	
	
	my $snp_file = shift or die "Can't access snp file\n";
	my $allele_freq = shift or die "Can't access allele frequency\n";

	##how to run the thing...
	my $big_data = CandiSNP::get_positions_from_file(
		-file => $snp_file,
		-genome => "athalianaTair10"
		 );
		
	##once data is loaded, you can get a unique id for it ...
	my $filetag = CandiSNP::get_filetag($big_data);
	
	
	$big_data = CandiSNP::annotate_positions($big_data, -genome => "athalianaTair10", -filetag => $filetag);	

	my $filtered_big_data = CandiSNP::apply_filter($big_data, $allele_freq);
	my $all_genome_lengths = CandiSNP::genome_lengths('athalianaTair10');
	my ($big_scale_marks, $big_scale_labels) = CandiSNP::scale_marks($all_genome_lengths);
	my $R = CandiSNP::R;
	
	
	my $palette = CandiSNP::get_palette('contrast');
	my $big_image = CandiSNP::plot_data($R, $filtered_big_data, $filetag, $big_scale_marks, $big_scale_labels,$all_genome_lengths,$palette);

	##dump file of used snps to a csv
	CandiSNP::data_hash_to_file($filtered_big_data,$filetag,-format=>'long');
	#$R->stop;
	return $filetag;
}
1;