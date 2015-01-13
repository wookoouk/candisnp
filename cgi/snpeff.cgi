#!/usr/bin/perl -w
use CGI;
use strict;                             # load CGI routines
use JSON;
use File::Basename;
use IPC::Open2;
use Sort::Key::Natural qw(natsort);



#Basename for this script and therefore path to tmp dir
my $dirname = dirname(__FILE__);

my $upload_dir = "$dirname/../tmp";
my $bin = "$dirname/../bin";

my $q = CGI->new;
my $browser_file_name = $q->param('file');
my $species = $q->param('species');

##assume filename checked on client side
my $tmpfile = upload_file_to_tmp($browser_file_name, $upload_dir);

## do SNPeff
my $snp_eff_result = {};
$snp_eff_result = run_snpeff($bin, $tmpfile, $species, $snp_eff_result);

##send back the response
print $q->header('text/json');
my $json = JSON->new->allow_nonref;
print $json->encode($snp_eff_result);


#########################################################################################################################
##
sub upload_file_to_tmp{
	my ($filename,$upload_dir) = @_;
	my $rand = int(rand(10000000));
	my $new_file = $rand . ".csv";
	my $tmpfile = "$upload_dir/$new_file";

	open ( UPLOADFILE, ">$tmpfile" ) or die "$!";
	binmode UPLOADFILE;
	while ( <$filename> ){
		print UPLOADFILE;
	}
	close UPLOADFILE;
	return $tmpfile;
}

sub run_snpeff{
	my ($bin, $tmpfile, $species, $data) = @_;
	my($chld_out, $chld_in);
	my $pid = open2($chld_out, $chld_in, "java -jar -Xmx2g $bin/snpEff.jar -c $bin/snpEff.config -i txt -o txt -noLog  -noStats -canon -snp -no-downstream -no-upstream -no-utr $species $tmpfile");
	while (my $line = <$chld_out>){
		next if $line =~ m/^#/;
		next if $line =~ m/^\n$/;
		chomp $line;
		$data = parse_snpEff($data,$line);
	}
	return data_hash_to_json($data);
}


sub parse_snpEff{
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

sub data_hash_to_json{
	my %data = %{$_[0]};
	my @records = [];
	foreach my $chr (natsort keys %{$data} ){
		foreach my $pos (natsort keys %{$data{$chr}}){

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

							
			push @records, {
				"chromosome" => $chr,
				"position" => $pos,
				"reference_base" => $data{$chr}{$pos}{_ref},
				"alternate_base" => $data{$chr}{$pos}{_alt},
				"allele_freq" => $data{$chr}{$pos}{_allele_freq}, 
				"in_cds" => $data{$chr}{$pos}{_in_cds}, 
				"is_synonymous" => $data{$chr}{$pos}{_syn}, 
				"is_ctga" => $data{$chr}{$pos}{_ctga},
				"change" => $data{$chr}{$pos}{_nucs},
				"gene" => $data{$chr}{$pos}{_gene},
				"effect" => $data{$chr}{$pos}{_effect}
			};


		}
	}
	
}
