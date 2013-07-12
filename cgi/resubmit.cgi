#!/usr/bin/perl

use strict;
use warnings;
use CGI::Carp qw(fatalsToBrowser);
use Digest::MD5;
use File::Spec;
use File::Basename;
require "Main.pl";

my $uploaddir = '../public';
open(PARAMS, ">>$uploaddir/RESUBMIT.txt");

#my $maxFileSize = 0.5 * 1024 * 1024; # 1/2mb max file size...

use CGI;
my $IN = new CGI;

my $file;
#if ($IN->param('POSTDATA')) {
#  $file = $IN->param('POSTDATA');
#} else {
#  $file = $IN->upload('qqfile');
#}

my $file_name = $IN->param('file_name');
my $allele_freq = $IN->param('allele_freq');
my $dir = dirname( File::Spec->rel2abs("$uploaddir/$file_name") );

print PARAMS "$dir/$file_name\t$allele_freq\n";


my $filename = doWork("$dir/$file_name", $allele_freq);
print $IN->header();
print qq|{ "success": true, "file": "$filename.svg" }|;
