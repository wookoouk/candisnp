#!/usr/bin/perl

use strict;
use warnings;
use CGI::Carp qw(fatalsToBrowser);
use File::Spec;
use File::Basename;
use CGI;
require "Main.pl";

my $uploaddir = '../public';
my $IN = new CGI;

my $file;
if ($IN->param('POSTDATA')) {
  $file = $IN->param('POSTDATA');
} else {
  $file = $IN->upload('qqfile');
}

my $allele_freq = $IN->param('allele_freq');
my $genome = $IN->param('genome');
my $palette = $IN->param('palette');
my $filter_cen_snps = $IN->param('filter_snps');

print STDERR "filter snp val = $filter_cen_snps\n";
open(WRITEIT, ">$uploaddir/$file") or die "Cant write to $uploaddir/$file Reason: $!";

binmode(WRITEIT);

if ($IN->param('POSTDATA')) {
    print WRITEIT $file;
} else {
	while (<$file>) {
		print WRITEIT;
    }
}
close(WRITEIT);

my $check_size = -s "$uploaddir/$file";

my $dir = dirname( File::Spec->rel2abs("$uploaddir/$file") );

my $filename = doWork("$dir/$file", $allele_freq, $genome, $palette, $filter_cen_snps);

print $IN->header();
if ($check_size < 1) {
    print STDERR "ooops, its empty - gonna get rid of it!\n";
    print qq|{ "success": false, "error": "File is empty..." }|;
    print STDERR "file has been NOT been uploaded... \n";
}

else
{
  print qq|{ "success": true, "file": "$filename", "allele_frequency" :  $allele_freq}|;
}
