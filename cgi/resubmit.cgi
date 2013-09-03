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

my $file_name = $IN->param('file_name');
my $allele_freq = $IN->param('allele_freq');
my $genome = $IN->param('genome');
my $palette = $IN->param('palette');
my $dir = dirname( File::Spec->rel2abs("$uploaddir/$file_name") );

my $filename = doWork("$dir/$file_name", $allele_freq, $genome, $palette);
print $IN->header();
print qq|{ "success": true, "file": "$filename" }|;
