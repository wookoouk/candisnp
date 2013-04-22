package CandiSNP;

use 5.006;
use strict;
use warnings FATAL => 'all';
use Statistics::R;
$ENV{PATH} = "/usr/bin/";

=head1 NAME

CandiSNP - The great new CandiSNP!

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Quick summary of what the module does.

Perhaps a little code snippet.

    use CandiSNP;

    my $foo = CandiSNP->new();
    ...

=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 SUBROUTINES/METHODS

=head2 function1

=cut

sub function1 {
}

=head2 function2

=cut

sub function2 {
}


=head2 new_r

Returns a new R interface object

=cut

sub R{
	my $R = Statistics::R->new();
	my $cmd = <<EOF;
	my_func = function(x){ x }
EOF
	$R->run($cmd);
	$R->run(q'x=my_func(99)');
	warn "from module ", $R->get('x');
	return $R;
}



=head1 AUTHOR

Dan MacLean (TSL), C<< <dan.maclean at tsl.ac.uk> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-candisnp at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=CandiSNP>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc CandiSNP


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=CandiSNP>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/CandiSNP>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/CandiSNP>

=item * Search CPAN

L<http://search.cpan.org/dist/CandiSNP/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2013 Dan MacLean (TSL).

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see L<http://www.gnu.org/licenses/>.


=cut

1; # End of CandiSNP
