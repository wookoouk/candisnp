package CandiSNP;

use 5.006;
use strict;
use warnings FATAL => 'all';
use Statistics::R;
use Tie::Handle::CSV;
use Carp;
use Data::Dumper;
use Sort::Key::Natural qw(natsort);
use File::Basename;
use Digest::MD5 qw(md5_hex);
use IPC::Open2;

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
	my $cmd = <<'EOF';
	#get the R functions loaded
	##these are not neccesarily the final R functions.
	library("ggplot2")
	my_plot2 = function(x,colours){
	points = geom_point(position=position_jitter(height=.25,width=2), aes(colour=corrected_type))
	p = ggplot(x, aes(position,chromosome) ) + facet_grid(mutant ~ ., scales="free") + colours + points +  theme(panel.grid.major.y = element_line(size=6, colour= "white")   ) + theme(panel.background=element_rect(colour="gray90")) + theme(panel.grid.major.x = element_line(size=1, colour="gray90") ) + scale_x_continuous(breaks=c(0,10000000,20000000,30000000),labels=c("0Mb", "10Mb", "20Mb", "30Mb" ))
	return(p)
	}
	f = function(data,cutoff=0.7,nmutants=1){
	x = data;
	x = x[x$allele_frequency >= cutoff & x$in_X_mutants == nmutants, ]
	x = x[x$chromosome != "mitochondria" & x$method == "PileUp", ]
	x$chromosome = factor(x$chromosome, levels=rev(levels(x$chromosome)))
	x = x[x$corrected_type != 'indel' & x$corrected_type != "unknown error (likely annotation errors in GFF)",]
	x$chromosome = factor(x$chromosome)
	x$corrected_type = factor(x$corrected_type)
	x$method = factor(x$method)
	return(x)
	}
	get_colour_assignments = function(name_list, colour_list){
	names(colour_list) = levels(name_list)
	colors = scale_colour_manual(name = "SNP Type",values = colour_list)
	return(colors)
	}
	do_picture = function(data,cutoff,colours){
	dat = f(data,cutoff)
	filename = paste(cutoff, ".svg", sep="")
	svg(filename, width=8.3,height=11.7)
	p = my_plot2(dat,colours)
	print(p)
	dev.off()
	}
EOF
	$R->run($cmd);
	return $R;
}


##gets the ruby Bio::Synreport annotations for the positions provided 
##forks a child process and runs SNPeff, parses the result data back into the $data hash
sub annotate_positions{
	my $data = shift;
	my %opts = @_;
	my $bin = bin_folder();
	my $public_folder = public_folder();
	my $md5 = md5_hex(%{$data});
	my $tmpfile = $public_folder . "/" . $md5;
	$opts{-format} = 'short';
	_data_hash_to_file($data, $tmpfile, %opts);
	my($chld_out, $chld_in);
	my $pid = open2($chld_out, $chld_in, "java -Xmx2g -jar $bin/snpEff.jar -c $bin/snpEff.config -i txt -o txt -noLog  -noStats -canon -snp -no-downstream -no-upstream -no-utr $opts{-genome} $tmpfile");
	while (my $line = <$chld_out>){
		next if $line =~ m/^#/;
		chomp $line;
		$data = _parse_snpEff($data,$line);
	}
	unlink $tmpfile;
	return $data;
	
}

sub _parse_snpEff{
	my ($data,$line) = @_;
	my @data = split(/\t/, $line);
	
	
	my ($chr,$pos,$ref,$alt,$gene,$effect,$nucs ) = ($data[0],$data[1],$data[2], $data[3], $data[10],$data[15],$data[16]);
	#warn Dumper join(",",$chr,$pos,$ref,$alt,$gene,$effect,$nucs);
	$chr = 'Chr' . $chr if (grep /Chr$chr/, keys %{$data});

	return $data if defined $$data{$chr}{$pos}{_gene};

	if ($effect ne 'INTERGENIC' || $effect ne 'INTRON'){
		$$data{$chr}{$pos}{_in_cds} = "FALSE";
		$$data{$chr}{$pos}{_syn} = "FALSE";
	}
	else{
		$$data{$chr}{$pos}{_in_cds} = "TRUE";
		$$data{$chr}{$pos}{_syn} = $effect;
	}
	$$data{$chr}{$pos}{_gene} = $gene;
	return $data;
}

sub _data_hash_to_file{
	my $d = shift;
	my $file_name = shift;
	my %data = %{$d};
	my %opts = @_;
	open OUT, ">$file_name";
	print OUT _header() unless defined $opts{-format} and $opts{-format} eq 'short';
	foreach my $chr (natsort keys %data ){
		foreach my $pos (natsort keys %{$data{$chr}}){
			my @line = ($chr, $pos, $data{$chr}{$pos}{_ref},$data{$chr}{$pos}{_alt} ); 
			if (defined $opts{-format} and $opts{-format} eq 'short') {
				print OUT join("\t", @line), "\n";
			}
			else {
				push @line, ($data{$chr}{$pos}{_allele_freq}, $data{$chr}{$pos}{_in_cds}, $data{$chr}{$pos}{_syn}, $data{$chr}{$pos}{_ctga});
				print OUT join(",", @line);
			}
		}
	}
	close OUT;
}

sub _base_folder{
	my $dirname = dirname(__FILE__);
	$dirname =~ s/\/src\/blib\/lib//; 
	return $dirname;
}

sub bin_folder{
	return _base_folder() . "/bin";
}

sub public_folder{
	return _base_folder() . "/public";
}

sub _header{
	return "Chr,Pos,Ref,Alt,Allele_Freq,in_cds,is_synonymous,is_ctga\n";
}

sub _is_ctga{
	my ($ref,$alt) = @_;
	return "TRUE" if ( ($ref =~ /c/i and $alt =~ /t/i) or ($ref =~ /g/i and $alt =~ /a/i) );
	return "FALSE";
}


##from the uploaded text file, gets the positions that pass filters
##adds on whether they are synonymous / non_synonymous ... 
##returns as hash structure
sub get_positions_from_file{

	my %opts = @_;
	my $fh = _open_file($opts{-file}); 
	croak "bad file headers" unless _header_ok($fh);
	my $data = {};
	while (my $l = <$fh>){
		next unless _is_snp($l, -cutoff => $opts{-cutoff});
		$$data{$l->{'chr'}}{$l->{'pos'}}{_alt} = $l->{'alt'};
		$$data{$l->{'chr'}}{$l->{'pos'}}{_ref} = $l->{'ref'};
		$$data{$l->{'chr'}}{$l->{'pos'}}{_allele_freq} = $l->{'allele_freq'};
		$$data{$l->{'chr'}}{$l->{'pos'}}{_syn} = "NA";
		$$data{$l->{'chr'}}{$l->{'pos'}}{_ctga} = _is_ctga($l->{'ref'}, $l->{'alt'});
		$$data{$l->{'chr'}}{$l->{'pos'}}{_in_cds} = "NA";
	}
	return $data;
}

#returns the line if it passes the user supplied threshold
sub _is_snp{
	my $l = shift;
	my %opts = @_;
	if ($l->{'allele_freq'} >= $opts{-cutoff} ){
		return 1;
	}
	return 0;
}


##returns csv file object
sub _open_file{
	my $file = shift;
	my $fh = Tie::Handle::CSV->new($file, 
		header => 1, 
		key_case => 'any', 
		open_mode => '<'
		) || croak "couldn't open file $file\n\n";
	return $fh;	
}

##checks the header
## text file must have headers "Chr" "Pos" "Ref" "Alt" "Allele_Freq"
sub _header_ok{
	my $fh = shift;
	my @headers = split(/,/,$fh->header);
	foreach my $h (qw{Chr Pos Ref Alt Allele_Freq}){
		return 0 unless grep /$h/i, @headers;
	}
	return 1;
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
