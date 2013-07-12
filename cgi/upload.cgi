#!/usr/bin/perl

use strict;
use CGI::Carp qw(fatalsToBrowser);
use Digest::MD5;
use File::Spec;
use File::Basename;
require "Main.pl";

my $uploaddir = '../public';
open(PARAMS, ">>$uploaddir/PARAMS.txt");

#my $maxFileSize = 0.5 * 1024 * 1024; # 1/2mb max file size...

use CGI;
my $IN = new CGI;

my $file;
if ($IN->param('POSTDATA')) {
  $file = $IN->param('POSTDATA');
} else {
  $file = $IN->upload('qqfile');
}

#my $temp_id = $IN->param('temp_id');
my $allele_freq = $IN->param('allele_freq');


# make a random filename, and we guess the file type later on...
#my $name = Digest::MD5::md5_base64( rand );
#$name =~ s/\+/_/g;
#$name =~ s/\//_/g;


#mkdir("$uploaddir/$file");
#my ($ext) = $file =~ /(\.[^.]+)$/;
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

my $dir = dirname( File::Spec->rel2abs("$uploaddir/$file") );


my $filename = doWork("$dir/$file", $allele_freq);
print PARAMS "$dir/$file\t$allele_freq\t$filename.svg\n";
close(PARAMS);
print $IN->header();
if ($check_size < 1) {
    print STDERR "ooops, its empty - gonna get rid of it!\n";
    print qq|{ "success": false, "error": "File is empty..." }|;
    print STDERR "file has been NOT been uploaded... \n";
}
else
{
  print qq|{ "success": true, "file": "$filename.svg", "allele_frequency" :  $allele_freq}|;
}
