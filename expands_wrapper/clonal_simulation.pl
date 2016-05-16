#!/usr/bin/perl
use strict;
use POSIX;
my %clones = ("0.86"=>58,"0.71"=>65,"0.6"=>55,"0.2"=>66);
#my %clones = ("0.1"=>78,"0.9"=>55);
my $clones = join(",", keys %clones);
my $params = "clones_" . $clones;
my $seg = "simu_out_$params.seg";
my $maf = "simu_out_$params.maf";
my $depth = 300;
open SEG, ">$seg" or die "$!\n";
open MAF, ">$maf" or die "$!\n";
print MAF "Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\tReference_Allele\tTumor_Seq_Allele2\tt_ref_count\tt_alt_count\tn_ref_count\tn_alt_count\n";
my %genome_states = (1=>2,2=>2,3=>2,4=>2,5=>2,8=>2,9=>2,10=>2,11=>2,12=>2,13=>3,14=>4,19=>1,22=>2,6=>"LOH",7=>"LOH");

my %clone_genome_states = ("0.86"=>{3=>3,10=>3},"0.71"=>{11=>4,12=>3},"0.6"=>{1=>3},"0.2"=>{2=>3,4=>4});

my %chroms;
for my $clone (keys %clone_genome_states){
	my @chroms = keys(%{$clone_genome_states{$clone}});
	for(@chroms){
	    $chroms{$_}++;
	}
}
my @chroms = keys %chroms;
for my $chr(@chroms){
	for my $clone (keys %clone_genome_states){
		print "$clone $clone_genome_states{$clone}{$chr} ";
	}
	print "\n";
}

sub getRandPloidy{
    my $state = shift;
    my @ploidies;
    if($state == 2){
	@ploidies = (1,1,1,1,1,1,1,2);

    }
    elsif($state == 3){
	@ploidies = (1,1,1,1,2,2,2,2,3);
    }
    elsif($state ==4 ){
	@ploidies = (1,1,1,1,2,2,2,2,3,3,3,3,4);
    }
    my $ploidy = $ploidies[rand @ploidies];
    return($ploidy);
}



my %chrom_signals = (1=>2,2=>2,3=>2,4=>2,5=>2,7=>"LOH",10=>2,11=>2,12=>2,13=>2,14=>2,15=>2,16=>2,17=>2);
#my @chroms    = keys %{$clone_genome_states{"0.84"}};
my $chrom_size = 100000;
my $n = 1;
for my $c_freq (keys %clones){
	#print "$clones{$c_freq}\t$c_freq\n";
	while($clones{$c_freq}>0){
		$clones{$c_freq}--;
		my $mod = $clones{$c_freq}%2;
		my $status = $c_freq;
		my $chrom   = $chroms[rand @chroms];
		#next unless $chrom == 11;
		#print "$clones{$c_freq}\t$c_freq\t$chrom\n";
		my $state = $clone_genome_states{$c_freq}{$chrom};
		my $vaf;
		if($state eq "LOH"){
			#SNV added with high VAF to MAF file and indicated as an LOH region in SEG file
			#$vaf = $c_freq + (1-$c_freq)/2;
			next;
			$status .= "LOH_Ploidy_2";
		}
		elsif($state ==2){
		    my $ploidy = getRandPloidy($state);
			if(!$mod){

			    if($ploidy == 1){
				$vaf = $c_freq * 0.5;
				$status .= "CN_$state";
				$status.="_Ploidy_1";
			    }
			    else{
				$vaf = $c_freq;
                                $status .= "CN_$state";
                                $status.="_Ploidy_2";
			    }
			}
			else{
			    my $mutant_copies = $ploidy * $c_freq;
				#my $mutant_copies = $c_freq; #ploidy 1
				my $cnv_signal = 0;
				my $total_copies;
				for my $clone (keys %clone_genome_states ){
					my $chrom_state = $clone_genome_states{$clone}{$chrom};
					$total_copies+=$chrom_state* $clone;
					next if $chrom_state == 2;
					next unless $chrom_state;
					$cnv_signal = $chrom_state * $clone + 2 - 2*$clone ;
					#print "THIS: $cnv_signal ($clone,$chrom_state,$state)\n";
					#$cnv_signal += $this_cnv_signal * $clone;
					#0.9 * $genome_states{$chr} +2  -2*0.9
					#$cnv_signal +=1;
				}

				#print "$cnv_signal $chrom_signals{$chrom}\n";
				$chrom_signals{$chrom} = $cnv_signal if $cnv_signal > $chrom_signals{$chrom}; #n_clones
				$vaf = $mutant_copies / $total_copies;
				print "$chrom\tS:\t$state\t$ploidy\t$total_copies\t$chrom_signals{$chrom}\t$vaf\n";
			}
		}
		elsif($state == 3){
			if(!$mod){
				$status .= "CN_$state";

				if($clones{$c_freq}%4){
					$vaf = $c_freq * 0.33333;
					$status.="_Ploidy_1";
				}
				else{
					$vaf = $c_freq * 0.66666;
					$status.="_Ploidy_2";
				}
			}
			else{
				#CNV is in a different subclone, VAF will be diluted by extra copies present in other subclone but have normal ploidy
				#CNV signal would be the average of all clones, VAF would be offset by the number of non-mutant copies proportionally added by this subclonal CNV
				my $mutant_copies = $c_freq; #ploidy 1
				my $cnv_signal = 0;
				my $total_copies;
				for my $clone (keys %clone_genome_states ){
					my $chrom_state = $clone_genome_states{$clone}{$chrom};
					next if $chrom_state == 2;
					$cnv_signal = $chrom_state * $clone + 2 - 2*$clone ;
					#print "THIS: $cnv_signal ($clone,$chrom_state,$state)\n";
					#print "THIS: $this_cnv_signal\n";
					#$cnv_signal += $this_cnv_signal * $clone;
					if($clone == $c_freq){
						$total_copies+=$chrom_state * $clone;
						#print "same clone $clone $c_freq $chrom_state\n";
					}
					else{
						$total_copies+=$chrom_state * $clone;
						#print "diff clones $clone $c_freq $chrom_state\n";
					}
				}
				$chrom_signals{$chrom} = $cnv_signal if $cnv_signal > $chrom_signals{$chrom};
#				$chrom_signals{$chrom} = $cnv_signal;
				#$cnv_signal+=1;
				#$vaf = $mutant_copies / $total_copies;
				#print "S:\t$state\t$total_copies\t$cnv_signal\t$vaf\n";
			}
		}
		elsif($state ==4){
			if(!$mod){
				$status .= "CN_$state";
				if($clones{$c_freq}%2){
					$vaf = $c_freq * 0.25;
					$status.="_Ploidy_1";
				}
				else{
					$vaf = $c_freq * 0.75;
					$status.="_Ploidy_3";
				}
			}
			else{
				#CNV is in a different subclone, VAF will be diluted by extra copies present in other subclone but have normal ploidy
				#CNV signal would be the average of all clones, VAF would be offset by the number of non-mutant copies proportionally added by this subclonal CNV
				my $mutant_copies = $c_freq; #ploidy 1
				my $cnv_signal = 0;
				my $total_copies;
				for my $clone (keys %clone_genome_states ){
					my $chrom_state = $clone_genome_states{$clone}{$chrom};
					$cnv_signal = $chrom_state * $clone + 2 - 2*$clone ;
					#print "THIS: $cnv_signal ($clone,$chrom_state)\n";
					#print "THIS: $this_cnv_signal\n";
					#$cnv_signal += $this_cnv_signal * $clone;
					if($clone == $c_freq){
						$total_copies+=$chrom_state * $clone;
						#print "same clone $clone $c_freq $chrom_state\n";
					}
					else{
						$total_copies+=$chrom_state * $clone;
						#print "diff clones $clone $c_freq $chrom_state\n";
					}
				}
				$chrom_signals{$chrom} = $cnv_signal if $cnv_signal > $chrom_signals{$chrom};
#				$chrom_signals{$chrom} = $cnv_signal;
				#$cnv_signal+=1;
				$vaf = $mutant_copies / $total_copies;
				#print "S:\t$state\t$total_copies\t$cnv_signal\t$vaf\n";
			}
		}
		elsif($state == 1){
			$status .= "CN_$state";
			$vaf = $c_freq;
			$status.="_Ploidy_1";
		}
		my $nref;
		if($n%2){
			$nref = ceil($depth*$vaf) + int(rand(3));
		}
		else{
			$nref = ceil($depth*$vaf) - int(rand(3));
		}
		my $ref = $depth-$nref;
		my $m = $n+1;
		$n++;
		print MAF "$vaf\_C_$status\t$chrom\t$n\t$m\tA\tC\t$ref\t$nref\t0\t0\n";

	}
}

use Data::Dumper;
print Dumper %chrom_signals;

print SEG "sample\tchromosome\tstart\tend\tLOH_flag\tBAF\tlogratio\n";
for my $chr (keys %genome_states){
	if($genome_states{$chr} eq "LOH"){
		#just print 10 segments with the same BAF to mimic a series of LOH germline SNVs in the region
		
		for my $c_freq (keys %clones){
			for my $n(1..10){
				my $vaf = 1- ($c_freq + (1-$c_freq)/2);
				my $status = "nLOH";
				my $nref = ceil($depth*$vaf);
				my $ref = $depth-$nref;
				my $m = $n+1;
				print SEG "$vaf\_$c_freq\_$genome_states{$chr}\t$chr\t$n\t$m\t1\t$vaf\t0\n";
			}
		}
		next;
	}
	#this one gives non-exact values diluted by normal DNA
	#my $abs = 0.9 * $genome_states{$chr} +2  -2*0.9;
	#my $logratio;
	#if($abs == 0){
#		$logratio = "-Inf";
#	}
#	else{
#		$logratio = log($abs/2)/log(2);
#	}

	#If instead using exact estimates (rounded)
	my $logratio = log($genome_states{$chr}/2)/log(2);
	print SEG "nosample_$genome_states{$chr}\t$chr\t1\t$chrom_size\t0\t1\t$logratio\n";
}
