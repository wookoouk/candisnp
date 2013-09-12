#!/usr/bin/perl

use strict;
use warnings;
use CGI::Carp qw(fatalsToBrowser);
use File::Spec;
use File::Basename;
use CGI;

my $uploaddir = '../public';
my $IN = new CGI;

my $file;
if ($IN->param('POSTDATA')) {
  $file = $IN->param('POSTDATA');
} else {
  $file = $IN->upload('qqfile');
}


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



print $IN->header();
if ($check_size < 1) {
    print STDERR "ooops, its empty - gonna get rid of it!\n";
    print qq|{ "success": false, "error": "File is empty..." }|;
    print STDERR "file has been NOT been uploaded... \n";
}


  print qq|{ "success": true}|;
