package CandiSNP;

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

$ENV{PATH} = "/usr/bin/:/usr/local/bin/";

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
	#options(verbose=TRUE)
	#options(warn=1)
	#suppressPackageStartupMessages(library("ggplot2"))
	#suppressMessages ( library(ggplot2) )
	library(ggplot2)
	library(gridExtra)
	library(plyr)
	require(grid)
	
	candi_plot = function(x,colours,marks,labels,genome_lengths){
	points = geom_point(position=position_jitter(height=.4,width=2), aes(colour=type),alpha = 1 )
	facets = facet_grid(chromosome ~ ., scales="free", space="free")
	x_axis = theme(axis.text.x = element_text(face = "bold", size = 50))
	x_axis = theme(axis.title.x = element_text(size = 20))
	y_axis = theme(axis.text.y =element_text(face = "bold", size = 30))
	y_axis = theme(axis.title.y = element_text(size = 20))
	opts =  opts(strip.background = element_blank(), strip.text.x = element_blank(), strip.text.y = element_blank()) +
	opts(legend.position="top", panel.background = theme_rect(fill='grey99', colour='grey'))
	max_l = max(genome_lengths$length)
	rect = geom_rect(data=genome_lengths, aes(xmin=length, xmax=Inf, ymin=-Inf, ymax=Inf,x=NULL, y=NULL), fill='grey95' )
	
	peakfind <- function(x) 
	{
	  d <- density(x$position, adjust = 1, kernel = "gaussian")
	  peaks <- which(diff(sign(diff(d$y)))==-2)
	  data.frame(x=d$x[peaks],y=d$y[peaks])
	}
	
	#debug_file <- file("/Users/ethering/Desktop/r_output/debug.txt", "w")
	chrs <- unique(unlist(x$chromosome, use.names = TRUE))
	num_chrs <- length(chrs)
	#x_file <- file("/Users/ethering/Desktop/r_output/x_table.txt", "w")
	#write.table(x, file=x_file, sep = "\t", col.names = NA)
	
	dat <- subset(x, type=='Non-Synonymous in Coding Region C-T or G-A' | type=='Non-Synonymous in Coding Region', select=c(position,chromosome, type))
	#dat_file <- file("/Users/ethering/Desktop/r_output/dat_table.txt", "w") 
	#write.table(dat, file=dat_file, sep = "\t", col.names = NA)
	plot_types <- list()
	index <- 0
	plotlist <- list()
	for (i in 1:length(chrs))
	{
		#write(chrs[i], file=debug_file)
	  #all the snps for the current chromosome
	  snpsub <- subset(x, chromosome == chrs[i])
	  #all the non-syn coding snps for the current chromosome
	  sub <- subset(dat, chromosome == chrs[i])
	  
	  if (nrow(snpsub) > 2)
	  {
		index <- index + 1 
		#draw the scatter plot and then add it to the layout
		scatter <- ggplot(snpsub, aes(position, chromosome) ) + points + facets + opts + scale_x_continuous(breaks=marks,labels=labels, limits=c(1, max_l)) + x_axis + y_axis
		plotlist[[index]]=scatter
		plot_types[[index]]="scatter"
		#print(scatter, vp = viewport(layout.pos.row = index, layout.pos.col = 1))
	  }
	  if (nrow(sub) > 2)
	  {
		#move to the next row
		index <- index + 1 
		#find the peaks to plot
		mypeaks <- ddply(sub,.(chromosome),peakfind)
		#write(nrow(sub), file=debug_file)
		#draw the peaks
		#seg <- geom_segment(data = mypeaks, aes(x=x, xend=x, y=y, yend=0))
		denplot <- ggplot(sub, aes(position)) +
		  geom_density(alpha = 0.2) + facet_wrap(~ chromosome, ncol=1) + #draw the density plot line
		  opts + facets + scale_x_continuous(breaks=marks,labels=labels, limits=c(1, max_l)) + x_axis + y_axis
		  plotlist[[index]]=denplot
			plot_types[[index]]="denplot"
		#add it to the layout
		#print(denplot, vp = viewport(layout.pos.row = index, layout.pos.col = 1))
		#get the peak values on the x axis
		#xpeaks <- mypeaks$x
		#find the corresponding bins in the 'sub' dataset
		#d <- density(sub$position)
		#for (p in xpeaks)
		#{
		#  best <- which.min(abs(d$x - p))
		#  #the peakfind fuction seems to struggle with small data and hence gently sloping curves..
		#  #..get the ten flanking bins to make sure we capture the correct area
		#  range <- d$x[(best-10):(best+10)]
		#  #get the highest and lowest values
		#  min <- min(range)
		#  max <- max(range)
		#
		#  #get any snps found within the three bins
		#  bestdat <- subset(x, chromosome == chrs[i] & position >= min & position <= max & (type=='Non-Synonymous in Coding Region C-T or G-A' | type=='Non-Synonymous in Coding Region'))
		#  #print(bestdat)
		#}
	  }
	}
	
	#p = ggplot(x, aes(position,chromosome) )  + colours + points + scale_x_continuous(breaks=marks,labels=labels, limits=c(1, max_l)) + x_axis + y_axis + facets + opts + rect
	return_list <- list("plotlist" = plotlist, "plot_types" = plot_types)
	return(return_list)
	}
	
	get_colours = function(palette){
		colour_list = structure(palette, .Names = c("Non-Synonymous in Coding Region C-T or G-A", "Non-Synonymous in Coding Region", "Synonymous in Coding Region", "Non Coding Region"))
		scale_colour_manual(name = "SNP Type",values = colour_list)
	}
	
	get_height = function(list){
		return(length(levels(list)) * 1.2 )
	}
		
	save_picture = function(p,fname, h){
		ggsave(filename=fname,plot=p, height=h, width=16)
	}
	
	# Multiple plot function
	#
	# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
	# - cols:   Number of columns in layout
	# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
	#
	# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
	# then plot 1 will go in the upper left, 2 will go in the upper right, and
	# 3 will go all the way across the bottom.
	#
	multiplot <- function(..., plotlist=NULL, plot_types = NULL, file, cols=1, rows, layout=NULL) {
		require(grid)
		svg(filename = file,height=length(plotlist)*3, width=15)
		# Make a list from the ... arguments and plotlist
		plots <- c(list(...), plotlist)
		
		numPlots = length(plots)
	    
		scatters <- length(grep("scatter", plot_types))
		dens <- length(grep("denplot", plot_types))
		numRows = (scatters*2) + dens
		

		# If layout is NULL, then use 'cols' to determine layout
		if (is.null(layout)) {
		  # Make the panel
		  # ncol: Number of columns of plots
		  # nrow: Number of rows needed, calculated from # of cols
		  layout <- matrix(seq(1, cols * ceiling(numRows/cols)),
                     ncol = cols, nrow = ceiling(numRows/cols))
		}
		
		if (numPlots==1) {
		  print(plots[[1]])
		  
		} else {
		  # Set up the page
		  grid.newpage()
		  pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
		   rowIndex <- 1
		  # Make each plot, in the correct location
		  for (i in 1:numPlots) {
			# Get the i,j matrix positions of the regions that contain this subplot
			matchidx <- as.data.frame(which(layout == rowIndex, arr.ind = TRUE))
			
			if (plot_types[i] == "scatter")
			{
			  to <- matchidx$row +1
			  print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row:to,
											  layout.pos.col = matchidx$col))
			  rowIndex <- rowIndex +2
			}
			else
			{
			  print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
											  layout.pos.col = matchidx$col))
			  rowIndex <- rowIndex +1
			}
		  }
		}
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
	my ($R, $data,$filetag, $scale_marks, $scale_labels,$genome_lengths, $palette) = @_;
	my $public_folder = public_folder();
	my $tmpfile = $public_folder . "/" . $filetag . ".svg";
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
	my @conts;
	my @lengths;
	
	foreach my $chr (keys %{$data}){
		push @conts, $chr;
		#warn Dumper $$genome_lengths{$chr};
		push @lengths, $$genome_lengths{$chr};
	}
	$R->set('conts', \@conts);
	$R->set('lengths', \@lengths);
	
	$R->set('palette', $palette);
	my $cmd = <<'EOF';

	chrsf = as.factor(chrs)
	data = data.frame(chromosome=chrsf,position=posns,type=types)

	genome_lengths = data.frame(chromosome=conts,length=lengths)

	colours = get_colours(palette)
	height = get_height(data$chromosome)
	#plot = candi_plot(data,colours,marks,labels,genome_lengths)
	return_plots = candi_plot(data,colours,marks,labels,genome_lengths)
	plot <- return_plots$plotlist
	plot_types <- return_plots$plot_types
	#save_picture(plot,filename,(height*2))
	multiplot(plotlist=plot, plot_types=plot_types, cols=1, file=filename)
	
EOF
	$R->run($cmd);
	return $tmpfile;
}

sub get_palette{
	my $type = shift;
	my @pal = ("gray50", "gray50", "red", "gray50");
	if ($type eq 'gradient'){
		 @pal = ("#225EA8", "#41B6C4", "#A1DAB4", "#FFFFCC");
		return \@pal;
	}
	elsif($type eq 'contrast'){
		@pal = ("#D7191C", "#FDAE61", "#2C7BB6", "#ABD9E9");
		return \@pal;
	}
	elsif ($type eq 'diverging')
	{
		@pal = ("#D01C8B", "#F1B6DA", "#B8E186", "#4DAC26");
		return \@pal;
	}
	elsif ($type eq 'sequential')
	{
		@pal = ("#E31A1C", "#FD8D3C", "#FECC5C", "#FFFFB2");
		return \@pal;
	}
	elsif ($type eq 'diverse')
	{
		@pal = ("#CB181D", "#FB6A4A", "#FCAE91", "#FEE5D9");
		return \@pal;
	}
	elsif ($type eq 'CandiSNP')
	{
		my @pal = ("gray50", "gray50", "red", "gray50");
		return \@pal;
	}
	return \@pal;
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

#send it a hashref, returns a big number
sub get_filetag{
	return md5_hex(%{$_[0]});
}

##gets the ruby Bio::Synreport annotations for the positions provided 
##forks a child process and runs SNPeff, parses the result data back into the $data hash
sub annotate_positions{
	my $data = shift;
	my %opts = @_;
	my $bin = bin_folder();
	my $public_folder = public_folder();
	my $filetag = $opts{-filetag};
	my $tmpfile = $public_folder . "/" . $filetag;
	$opts{-format} = 'short';
	data_hash_to_file($data, $tmpfile, %opts);
	my($chld_out, $chld_in);
	my $pid = open2($chld_out, $chld_in, "java -jar -Xmx2g $bin/snpEff.jar -c $bin/snpEff.config -i txt -o txt -noLog  -noStats -canon -snp -no-downstream -no-upstream -no-utr $opts{-genome} $tmpfile");
	while (my $line = <$chld_out>){
		
		next if $line =~ m/^#/;
		next if $line =~ m/^\n$/;
		chomp $line;
		#warn Dumper $line;
		$data = _parse_snpEff($data,$line);
	}
	#unlink $tmpfile;
	return $data;
	
}

sub _parse_snpEff{
	my ($data,$line) = @_;
	my @data = split(/\t/, $line);
	
	
	my ($chr,$pos,$ref,$alt,$gene,$effect,$nucs ) = ($data[0],$data[1],$data[2], $data[3], $data[10],$data[15],$data[16]);
	#warn Dumper join(",",$chr,$pos,$ref,$alt,$gene,$effect,$nucs);
	#warn Dumper $effect;
	#$chr = 'Chr' . $chr if (grep /Chr$chr/, keys %{$data});
	
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

	$gene = "" unless defined $gene;
	$nucs = "" unless defined $nucs;
	$effect = "" unless defined $effect;
	$$data{$chr}{$pos}{_gene} = $gene =~ m/\w+/ ? $gene : "-";
	$$data{$chr}{$pos}{_nucs} = $nucs =~ m/\w+/ ? $nucs : "-";
	$$data{$chr}{$pos}{_effect} = $effect =~ m/\w+/ ? $effect : "-";
	return $data;
}

sub data_hash_to_file{
	my $d = shift;
	my $file_name = shift;
	my %data = %{$d};
	my %opts = @_;
	$file_name = public_folder() . "/" . $file_name . '.csv' if ($opts{-format} eq 'long');
	open my $OUT, ">", $file_name;
	print $OUT _header() unless defined $opts{-format} and $opts{-format} eq 'short';
	foreach my $chr (natsort keys %data ){
		foreach my $pos (natsort keys %{$data{$chr}}){
			my @line = ($chr, $pos, $data{$chr}{$pos}{_ref},$data{$chr}{$pos}{_alt} ); 
			if (defined $opts{-format} and $opts{-format} eq 'short') {
				print $OUT join("\t", @line), "\n";
				#warn Dumper @line;
			}
			else {
				if (not defined ($data{$chr}{$pos}{_gene}) )
				{
					$data{$chr}{$pos}{_gene}= "NA";
				}
				if (not defined ($data{$chr}{$pos}{_nucs}) )
				{
					$data{$chr}{$pos}{_nucs}= "NA";
				}
				if (not defined ($data{$chr}{$pos}{_effect}) )
				{
					$data{$chr}{$pos}{_effect}= "NA";
				}
								
				push @line, ($data{$chr}{$pos}{_allele_freq}, $data{$chr}{$pos}{_in_cds}, $data{$chr}{$pos}{_syn}, $data{$chr}{$pos}{_ctga},$data{$chr}{$pos}{_nucs},$data{$chr}{$pos}{_gene},$data{$chr}{$pos}{_effect});
				#warn Dumper @line;
				print $OUT join(",", @line), "\n";
				
			}
		}
	}
	close $OUT;
}

sub _base_folder{
	my $dirname = dirname(__FILE__);
	$dirname =~ s/\/src\/lib//; 
	return $dirname;
}

sub bin_folder{
	return _base_folder() . "/bin";
}

sub public_folder{
	return _base_folder() . "/public";
}

sub _header{
	return "Chr,Pos,Ref,Alt,Allele_Freq,in_cds,is_synonymous,is_ctga,change,gene,summary\n";
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
	my %centromeres;
	my $filter = $opts{-filter};
	my $windowSize = 500000;#the number of nucleotides up and down stream of the centromere to ignore
	if ($filter eq "yes") {
		%centromeres = %{centromere_positions()};
		print STDERR Dumper %centromeres;
	}
	
	my $data = {};
	while (my $l = <$fh>){
		next unless defined $genome{$l->{'chr'}}; #skip any positions on chromosomes not in our genome definition...
		#skip any positions on chromosomes that's within a centromere...

		next if ($filter eq "yes" and ($l->{'pos'} < ($centromeres{$l->{'chr'}}+$windowSize) and $l->{'pos'} > ($centromeres{$l->{'chr'}}-$windowSize) ));
		#{
		#	print STDERR "****SKIPPING****\n";
		#	print STDERR "chr  = $l->{'chr'}\n";
		#	print STDERR "position  = $l->{'pos'}\n";
		#	next;
		#}
		#next if ($filter eq "yes" and ($l->{'pos'} < ($centromeres{$l->{'chr'}}+50) and $l->{'pos'} > ($centromeres{$l->{'chr'}}-50) ));
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

#these are the centromere centres, from which full centromere may be found +/- 500nt up/down stream
sub centromere_positions{
	my $centromere_centres = {
			"1" => 15086545,
			"2" => 3608429,
			"3" => 14209452,
			"4" => 3956521,
			"5" => 11725524
	};
	return $centromere_centres;
}

#returns list of contig/chromosome lengths for a given genome file
sub genome_lengths{
	my $genome = shift;
	my $lengths = {
		'athalianaTair10' => {
			"1" => 34964571,
			"2" => 22037565,
			"3" => 25499034,
			"4" => 20862711,
			"5" => 31270811
		},
		
		'athalianaTair9' => {
			"1" => 30427671,
			"5" => 26975502,
			"3" => 23459830,
			"2" => 19698289,
			"4" => 18585056
		},
		
		'rice7' => {
			"1" => 43270923,
			"2" => 35937250,
			"3" => 36413819,
			"4" => 35502694,
			"5" => 29958434,
			"6" => 31248787,
			"7" => 29697621,
			"8" => 28443022,
			"9" => 23012720,
			"10" => 23207287,
			"11" => 29021106,
			"12" => 27531856
		},
		
		'SL2.40' =>  {
			"SL2.40ch01" => 90304244,
			"SL2.40ch09" => 67662091,
			"SL2.40ch12" => 65486253,
			"SL2.40ch07" => 65268621,
			"SL2.40ch05" => 65021438,
			"SL2.40ch03" => 64840714,
			"SL2.40ch10" => 64834305,
			"SL2.40ch04" => 64064312,
			"SL2.40ch08" => 63032657,
			"SL2.40ch11" => 53386025,
			"SL2.40ch02" => 49918294,
			"SL2.40ch06" => 46041636,
			"SL2.40ch00" => 21805821
		},
		
		'gmax1.09v8' => {
			"Gm18" => 62308140,
			"Gm01" => 55915595,
			"Gm02" => 51656713,
			"Gm10" => 50969635,
			"Gm15" => 50939160,
			"Gm06" => 50722821,
			"Gm19" => 50589441,
			"Gm14" => 49711204,
			"Gm04" => 49243852,
			"Gm03" => 47781076,
			"Gm08" => 46995532,
			"Gm09" => 46843750,
			"Gm20" => 46773167,
			"Gm07" => 44683157,
			"Gm13" => 44408971,
			"Gm05" => 41936504,
			"Gm17" => 41906774,
			"Gm12" => 40113140,
			"Gm11" => 39172790,
			"Gm16" => 37397385
		},
		
		'grape' => {
			"Un" => 43154196,
			"14" => 30274277,
			"18" => 29360087,
			"5" => 25021643,
			"13" => 24396255,
			"19" => 24021853,
			"4" => 23867706,
			"1" => 23037639,
			"9" => 23006712,
			"12" => 22702307,
			"8" => 22385789,
			"16" => 22053297,
			"6" => 21508407,
			"7" => 21026613,
			"15" => 20304914,
			"11" => 19818926,
			"3" => 19341862,
			"2" => 18779844,
			"10" => 18140952,
			"17" => 17126926
		},
		
		'maizeZmB73' => {
			"1" => 301354135,
			"4" => 241473504,
			"2" => 237068873,
			"3" => 232140174,
			"5" => 217872852,
			"7" => 176764762,
			"8" => 175793759,
			"6" => 169174353,
			"9" => 156750706,
			"10" => 150189435
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
