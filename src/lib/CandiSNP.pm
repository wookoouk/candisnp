package CandiSNP;

use Modern::Perl 2011;
use autodie;
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
use Storable 'dclone';
use Number::Bytes::Human qw(format_bytes parse_bytes);
use feature qw/switch/; 

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
	options(warn=-1)
	suppressPackageStartupMessages(library("ggplot2"))
	suppressMessages ( library(ggplot2) )
	candi_plot = function(x,colours,marks,labels){
	points = geom_point(position=position_jitter(height=.25,width=2), aes(colour=type))
	facets = facet_grid(chromosome ~ .,scales="free", space="free_x")
	x_axis = theme(axis.title.x = element_blank())
	y_axis = theme(axis.title.y = element_blank())
	p = ggplot(x, aes(position,chromosome) ) + colours + points + scale_x_continuous(breaks=marks,labels=labels) + x_axis + y_axis + facets
	return(p)
	}
	
	get_colours = function(){
		colour_list = structure(c("gray50", "gray50", "red", "gray50"), .Names = c("Synonymous in Coding Region", "Non-Synonymous in Coding Region", "Non-Synonymous in Coding Region C-T or G-A", "Non Coding Region"))
		scale_colour_manual(name = "SNP Type",values = colour_list)
	}
	
	get_height = function(list){
		return(length(levels(list)) * 1.2 )
	}
		
	save_picture = function(p,filename, height){
	svg(filename, height=height, width=8)
	print(p)
	dev.off()
	}
EOF

	$R->run($cmd);
	return $R;
}




#takes the filtered annotated data and writes it into an R data frame with headings:
#chromosome position type
#creates the plot and returns its filename.
sub plot_data{
	my ($R, $data, $scale_marks, $scale_labels) = @_;
	my $md5 = md5_hex(%{$data});
	my $public_folder = public_folder();
	my $tmpfile = $public_folder . "/" . $md5 . ".svg";
	$R->set("marks", $scale_marks);
	$R->set("labels",$scale_labels);
	$R->set("filename",$tmpfile);
	my @chrs;
	my @posns;
	my @types;
	foreach my $chr (natsort keys %{$data}){
		foreach my $pos (keys %{$$data{$chr}}){
			my $type = _get_type($$data{$chr}{$pos});
			if (defined $type){
				push @chrs, $chr;
				push @posns,$pos;
				push @types,$type;
			}
		}
	}
	
	#bung the data into R lists
	$R->set('chrs', \@chrs);
	$R->set('posns', \@posns);
	$R->set('types', \@types);
	#warn Dumper @types;
	my $cmd = <<'EOF';
	data = data.frame(chromosome=chrs,position=posns,type=types)
	colours = get_colours()
	height = get_height(data$chromosome)
	plot = candi_plot(data,colours,marks,labels)
	save_picture(plot,filename,height)
	
EOF
	$R->run($cmd);
	return $tmpfile;
}

sub _get_type{
	my %data = %{$_[0]};
	my $result = undef;
	#warn Dumper %data;
	if ($data{_syn} eq "FALSE" and $data{_ctga} eq "TRUE" and $data{_in_cds} eq "TRUE"){
		$result = "Non-Synonymous in Coding Region C-T or G-A";
	}
	elsif($data{_syn} eq 'FALSE' and $data{_ctga} eq "FALSE" and $data{_in_cds} eq "TRUE" ){
		$result = "Non-Synonymous in Coding Region";
	}
	elsif($data{_syn} eq 'TRUE' and $data{_ctga} eq "FALSE" and $data{_in_cds} eq "TRUE" ){
		$result = "Synonymous in Coding Region";
	}
	elsif( $data{_in_cds} eq "FALSE" ){
		$result = "Non Coding Region";
	}
	else{
		$result = undef; #unclassified SNP
	}
	return $result;
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
		next if $line =~ m/^\n$/;
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
	warn Dumper $effect;
	$chr = 'Chr' . $chr if (grep /Chr$chr/, keys %{$data});
	
	my $syn = "TRUE";
	if ($effect eq "NON_SYNONYMOUS_CODING"){
		$syn = "FALSE";
	}
	return $data if defined $$data{$chr}{$pos}{_gene};

	if ($effect eq 'INTERGENIC' || $effect eq 'INTRON'){
		$$data{$chr}{$pos}{_in_cds} = "FALSE";
		$$data{$chr}{$pos}{_syn} = "FALSE";
	}
	else{
		$$data{$chr}{$pos}{_in_cds} = "TRUE";
		$$data{$chr}{$pos}{_syn} = $syn;
	}
	$$data{$chr}{$pos}{_gene} = $gene;
	return $data;
}

sub _data_hash_to_file{
	my $d = shift;
	my $file_name = shift;
	my %data = %{$d};
	my %opts = @_;
	open my $OUT, ">", $file_name;
	print $OUT _header() unless defined $opts{-format} and $opts{-format} eq 'short';
	foreach my $chr (natsort keys %data ){
		foreach my $pos (natsort keys %{$data{$chr}}){
			my @line = ($chr, $pos, $data{$chr}{$pos}{_ref},$data{$chr}{$pos}{_alt} ); 
			if (defined $opts{-format} and $opts{-format} eq 'short') {
				print $OUT join("\t", @line), "\n";
			}
			else {
				push @line, ($data{$chr}{$pos}{_allele_freq}, $data{$chr}{$pos}{_in_cds}, $data{$chr}{$pos}{_syn}, $data{$chr}{$pos}{_ctga});
				print $OUT join(",", @line);
			}
		}
	}
	close $OUT;
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


##from the uploaded text file, gets the positions that are on chromosomes we can use...
##returns as hash structure
sub get_positions_from_file{

	my %opts = @_;
	my $fh = _open_file($opts{-file}); 
	croak "bad file headers" unless _header_ok($fh);
	my %genome = %{genome_lengths($opts{-genome})};
	my $data = {};
	while (my $l = <$fh>){
		next unless defined $genome{$l->{'chr'}}; #skip any positions on chromosomes not in our genome definition...
		warn Dumper "skipping $l->{'chr'} : $l->{'pos'}" unless defined $genome{$l->{'chr'}};
		$$data{$l->{'chr'}}{$l->{'pos'}}{_alt} = $l->{'alt'};
		$$data{$l->{'chr'}}{$l->{'pos'}}{_ref} = $l->{'ref'};
		$$data{$l->{'chr'}}{$l->{'pos'}}{_allele_freq} = $l->{'allele_freq'};
		$$data{$l->{'chr'}}{$l->{'pos'}}{_syn} = "NA";
		$$data{$l->{'chr'}}{$l->{'pos'}}{_ctga} = _is_ctga($l->{'ref'}, $l->{'alt'});
		$$data{$l->{'chr'}}{$l->{'pos'}}{_in_cds} = "NA";
	}
	if (scalar (keys %{$data}) == 0){
		carp "Couldnt find any chromosomes with those names in this genome";
		return 0;
	}
	return $data;
}

#filter out SNPs with allele freq less than that specified, returns a new hash
sub apply_filter{
	my $data = shift;
	my $cutoff = shift;
	my $new = dclone($data);
	foreach my $chr (keys %{$data} ){
		foreach my $pos (keys %{$$data{$chr}}){
			if ( $$new{$chr}{$pos}{_allele_freq} < $cutoff ){
				 delete $$new{$chr}{$pos};
			}
		}
	}
	return $new;
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

#returns list of contig/chromosome lengths for a given genome file
sub genome_lengths{
	my $genome = shift;
	my $lengths = {
		'athalianaTair10' => {
			"Chr1" => 34964571,
			"Chr2" => 22037565,
			"Chr3" => 25499034,
			"Chr4" => 20862711,
			"Chr5" => 31270811
		}
	};
	return $$lengths{$genome};
}


sub scale_marks{
	my $lengths = shift;
	my @sorted = reverse natsort values %{$lengths} ;
	my $longest = shift @sorted;
	my $interval = scale_units($longest);
	my @marks = ($interval);
	foreach my $gap (@marks){
		my $next = $interval * (scalar @marks + 1);
		last if $next >= $longest;
		push @marks, $next;
	}
	my @labels = ();
	foreach my $mark (@marks){
		push @labels, format_bytes($mark, bs => 1000) . 'b';
	}
	return \@marks, \@labels;
}


sub scale_units{
	my $result = 1000000;
		if (length $_[0] == 3){$result = 100;}
		elsif (length $_[0] == 4){$result = 1000;}
		elsif (length $_[0] == 5){$result = 10000;}
		elsif (length $_[0] == 6){$result = 100000;}
		elsif (length $_[0] == 7){$result = 1000000;}
		elsif (length $_[0] == 8){$result = 10000000;}
	return $result;
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
