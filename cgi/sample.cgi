#!/usr/bin/perl -w
use CGI;
use strict;                             # load CGI routines
use JSON;

my $q = CGI->new;
my $json = JSON->new->allow_nonref;
my $data = $json->encode("hello, im json");

print $q->header('text/json');
print $data;
