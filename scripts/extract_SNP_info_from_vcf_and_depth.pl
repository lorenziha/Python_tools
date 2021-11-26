#!/Users/lorenziha/bin/miniconda3/bin/perl

use strict;

my $usage =  "$0 -b bin_size -v <vcf file - single sample> -d <depth file>\n\n";

my %arg = @ARGV;
die $usage unless $arg{-b} and $arg{-v} and $arg{-d};
my $binparam = $arg{-b};

# for binning, want the following output:
#   range (0-10) -> 10
#   range (11-20) -> 20
#          (21-30) -> 30
#                         etc.
#  formula is:
#  x = number
#  binned(x) at seg-y = int ( (x-1+y)/y)

# collect chromosome size data
# chrom_id <tab> size in bp
my %chr;

## Iterate though the depth file to collect depth data and count of bases
## with coverage of 1 or more.
## CM000429	159	5 (chr_id, pos, cov)
my $old_id = "none";
my (%DEPTH, %COV_COUNTS, %DEPTH_MEAN);

open (DEPTH, "<$arg{-d}") || die "ERROR, I cannot open $arg{-v}: $!\n\n";
while(<DEPTH>){
    chomp;
    next if m/^[\s\t]*$/;
    my ($chr_id, $pos, $depth) = split /\t/;
    my $idx = int($pos / $binparam) + 1;
    $COV_COUNTS{$chr_id}->{$idx}++ if $depth > 0;
    $DEPTH{$chr_id}->{$idx} += $depth;
    $DEPTH_MEAN{$chr_id}->{$idx} = $DEPTH{$chr_id}->{$idx} / $binparam;
}
close DEPTH;

# Iterate through chromosome IDs in VCF file
# ./.:0:0:.:.:.:0
my $old_chr = 'none';
my @snps;
my @match;
my @bins_per_chr;

open(VCF, "<$arg{-v}") || die "ERROR, I cannot open $arg{-v}: $!\n\n";

# Print table header
print "CHROM\tBIN_START\tSNP_COUNT\tNORM_SNP_COUNT\tVARIANTS/KB\tDEPTH_MEAN\tWINDOW_COV\n";

while(<VCF>){
	next if m/^[\s\t]*$/;
    if(m/^##contig=<ID=(\S+),length=(\d+)>/){
        $chr{$1} = $2;
    }
	next if m/^#/;
	my @x = split /\t/;
	my ($chrom_id, $pos, $info) = ($x[0], $x[1], $x[9]);
	if ($chrom_id ne $old_chr && $old_chr ne 'none'){
        #print "($chrom_id ne $old_chr && $old_chr ne 'none')\n\n";
        # Send %snps to analysis and print results
        my ($track_ref, $tc, $vc, $vdc_ref) = &generate_bins($old_chr, @snps);  # return(\%tracker, $total_counts, $value_counter, \%value_dist_counter)
        # For matches
        my ($mtrack_ref, $mtc, $mvc, $mvdc_ref) = &generate_bins($old_chr, @match);
        
        &print_results($old_chr, $track_ref, $tc, $vc, $vdc_ref, $mtrack_ref, $mtc);
        # Reset @snps and @match for new chromosome
        @snps = ();
        @match = ();
        
	}
    $old_chr = $chrom_id;
	my ($gt, @z) = split(":", $info);
	
    # skip SNPs = REF allele or missing SNPs.
	if ($gt eq '1/1' || $gt eq '0/1'){
        push @snps, $pos; 
    }
    elsif ($gt eq '0/0'){
        push @match, $pos;
    }
    #print "$chrom_id\t$pos\n";
	
}


#For last chromosome
my ($track_ref, $tc, $vc, $vdc_ref) = &generate_bins($old_chr, @snps);  # return(\%tracker, $total_counts, $value_counter, \%value_dist_counter)
# For matches
my ($mtrack_ref, $mtc, $mvc, $mvdc_ref) = &generate_bins($old_chr, @match);
&print_results($old_chr, $track_ref, $tc, $vc, $vdc_ref, $mtrack_ref, $mtc);


close VCF;

exit(0);
###########################################################################

#analyse data
sub generate_bins {
    my @my_snps = @_; #print "\n\n$my_snps[0] - $my_snps[1] - $my_snps[2] - $my_snps[3] - $my_snps[-1]\n\n";
    my $chr_id = shift(@my_snps);
	my $include_zero_counts = 1;
	my $size = $chr{$chr_id}; ########################### edit this line to include chr size !!!!!!!!!!!!!!!!!!!!
	my $value_counter = 0;
	my %value_dist_counter;

	my $entry_counter = 0;

	my %tracker;
	foreach my $snp (@my_snps) {
    		$snp =~ s/\s+//g;
    		my $orig = $snp;
    		unless ($snp =~ m/^[+-]?[\d\.]+$/) {
		        print STDERR "can't process $snp\n";
		        next;
    		}
    
    		my $val_examine = $orig;
    		if ($val_examine > 0) {
			    $val_examine--;
    		} elsif ($val_examine < 0) {
			    $val_examine++;
    		}
    	
    		my $index = int ( $val_examine/ $binparam); 
    		if ($orig > 0) {
			    $index++;
    		} elsif ($orig < 0) {
			    $index--;
    		}
    
    		$tracker{$index}++;

    		$value_counter += abs($orig);
    		$value_dist_counter{$index} += abs($orig);

    		$entry_counter++;

	}

	my $total_counts = $entry_counter;
    return(\%tracker, $total_counts, $value_counter, \%value_dist_counter);
}

#print results
sub print_results {
    # ($chrom_id, $track_ref, $tc, $vc, $vdc_ref, $mtrack_ref, $mtc)
    my ($chrom_id, $track_ref, $total_counts, $value_counter, $vdc_ref, $mtrack_ref, $m_total_counts) = @_;
    my %value_dist_counter = %{$vdc_ref};
    my %tracker = %{$track_ref}; 
	my @values = sort {$a<=>$b} keys %tracker;
	my $maxvalue = $values[$#values];
	my $minvalue = 1; # $values[0];
    my $size = $chr{$chrom_id};
    my $include_zero_counts = 1;
    #print STDERR "$mtrack_ref\n";
    my %mtracker;
    if ($mtrack_ref){
        %mtracker = %{$mtrack_ref};
    } 
    else {
        
    }

    #print STDERR "max iter = ".int($size / $binparam)."\n";
	
	#for (my $i = $minvalue; $i <= $maxvalue; $i++) {
	for (my $i = $minvalue; $i <= int($size / $binparam); $i++) {
            
    		my $index = ($i)  * ($binparam);
    		my $value = (exists($tracker{$i})) ? sprintf ("%.2f", ($tracker{$i}))  : 
                (exists($mtracker{$i})) ? 0 : -10; # assign -10 to bins with no sequencing data
            my $nor_value = (exists($tracker{$i})) ? sprintf ("%.2f", ($tracker{$i} * $binparam/$COV_COUNTS{$chrom_id}->{$i}))  : (exists($mtracker{$i})) ? 0 : -10; # assign -10 to bins with no sequencing data
            #if ($value != 0 || 1) {

                #my $percentage_counts = sprintf ("%.2f", $value / $total_counts * 100);
                #my $percentage_values = sprintf ("%.2f", $value_dist_counter{$i}/$value_counter * 100);
                my $snp_per_kb = sprintf ("%.2f", ($nor_value * 1000) / $binparam);

                #print "$index\t$value\t$percentage_counts\t$percentage_values\n";
                #print "CHROM\tBIN_START\tSNP_COUNT\tNORM_SNP_COUNT\tVARIANTS/KB\tDEPTH_MEAN\tWINDOW_COV\n";
                print "$chrom_id\t$index\t$value\t$nor_value\t$snp_per_kb\t".$DEPTH_MEAN{$chrom_id}->{$i}."\t".sprintf ("%.2f", (100 * $COV_COUNTS{$chrom_id}->{$i} / $binparam))."\n";
    		#}
	}
}


