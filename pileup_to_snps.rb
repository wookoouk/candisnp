#!/usr/bin/ruby
# encoding: utf-8
#
#  pileup_to_snps.rb
#
#  Created by Dan MacLean (TSL) on 2013-01-09.
#  Copyright (c). All rights reserved.
#
require 'rubygems'
require 'bio-samtools'

class Bio::DB::FastaLengthDB
  require 'bio'
  def initialize(args)
    @file = args[:file]
    @seqs = {}
    file = Bio::FastaFormat.open(@file)
    file.each do |entry|
      @seqs[entry.entry_id] = entry.length
    end
    
    def each
      @seqs.keys.sort.each do |k|
        yield k, @seqs[k]
      end
    end
    
  end
end



sequences = Bio::DB::FastaLengthDB.new(:file => ARGV[0])

bam = Bio::DB::Sam.new(:bam=>ARGV[1], :fasta=>ARGV[0])
bam.open

puts ["Chr", "Pos","Ref","Alt","Allele_Freq"].join(",")

sequences.each do |id,length|
  $stderr.puts "on #{id}:1-#{length}"
  bam.mpileup(
    :r => "#{id}:#{1}-#{length}",
    :Q => 20,
    :q => 20
  ) do |pileup|
    if pileup.is_snp?(:ignore_reference_n => true, :min_depth => 6, :min_non_ref_count => 3) and pileup.consensus != pileup.ref_base
      mut = "FALSE"
      mut = "TRUE" if (pileup.ref_base == 'G' and pileup.consensus == 'A') or (pileup.ref_base == 'C' and pileup.consensus == 'T')
      puts [pileup.ref_name, pileup.pos, pileup.ref_base, pileup.consensus, pileup.non_ref_count/pileup.coverage, mut].join(",")
    end
  end
end
