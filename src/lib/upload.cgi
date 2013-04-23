#!/usr/bin/perl

use strict;
use CGI::Carp qw(fatalsToBrowser);

use Digest::MD5;

my $uploaddir = '../../uploads';

#my $maxFileSize = 0.5 * 1024 * 1024; # 1/2mb max file size...

use CGI;
my $IN = new CGI;

my $file;
if ($IN->param('POSTDATA')) {
  $file = $IN->param('POSTDATA');
} else {
  $file = $IN->upload('qqfile');
}

binmode(WRITEIT);
open(WRITEIT, ">$uploaddir/$file") or die "Cant write to $uploaddir/$file Reason: $!";
if ($IN->param('POSTDATA')) {
    print WRITEIT $file;
} else {
    while (<$file>) {
      print WRITEIT;
    }
}
close(WRITEIT);

my $check_size = -s "$uploaddir/$file";

#print STDERR qq|Main filesize: $check_size  Max Filesize: $maxFileSize \n\n|;

print $IN->header();
if ($check_size < 1) {
    print STDERR "ooops, its empty - gonna get rid of it!\n";
    print qq|{ "success": false, "error": "File is empty..." }|;
    print STDERR "file has been NOT been uploaded... \n";
} 

else  {
    print qq|{ "success": true }|;

    print STDERR "file has been successfully uploaded... thank you.\n";
}
