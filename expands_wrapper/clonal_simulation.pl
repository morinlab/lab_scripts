#!/usr/bin/perl
use strict;
use POSIX;
my %clones = ("0.89"=>81,"0.64"=>78,"0.39"=>55);

my $clones = join(",", keys %clones);
my $params = "clones_" . $clones;
my $seg = "simu_out_$params.seg";
my $maf = "simu_out_$params.maf";

open SEG, ">$seg" or die "$!\n";
open MAF, ">$maf" or die "$!\n";
print MAF "Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\tReference_Allele\tTumor_Seq_Allele2\tt_ref_count\tt_alt_count\tn_ref_count\tn_alt_count\n";
my %genome_states = (1=>2,2=>2,3=>2,4=>2,5=>2,8=>2,9=>2,10=>2,11=>2,12=>2,13=>3,14=>4,19=>1,22=>2,6=>"LOH",7=>"LOH");
my @chroms    = keys %genome_states;
my $chrom_size = 100000;
my $n = 1;
for my $c_freq (keys %clones){
	#print "$clones{$c_freq}\t$c_freq\n";
	while($clones{$c_freq}>0){
		$clones{$c_freq}--;
		my $status = $c_freq;
		my $chrom   = $chroms[rand @chroms];
		#print "$clones{$c_freq}\t$c_freq\t$chrom\n";
		my $state = $genome_states{$chrom};
		my $vaf;
		if($state eq "LOH"){
			#SNV added with high VAF to MAF file and indicated as an LOH region in SEG file
			#$vaf = $c_freq + (1-$c_freq)/2;
			next;
			$status .= "LOH_Ploidy_2";
		}
		elsif($state ==2){
			$vaf = $c_freq * 0.5;
			$status .= "CN_$state";
			$status.="_Ploidy_1";
		}
		elsif($state == 3){
			$status .= "CN_$state";
			if($clones{$c_freq}%2){
				$vaf = $c_freq * 0.33333;
				$status.="_Ploidy_1";
			}
			else{
				$vaf = $c_freq * 0.66666;
				$status.="_Ploidy_2";
			}
		}
		elsif($state ==4){
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
		elsif($state == 1){
			$status .= "CN_$state";
			$vaf = $c_freq;
			$status.="_Ploidy_1";
		}
		my $nref;
		if($n%2){
			$nref = ceil(100*$vaf) + int(rand(3));
		}
		else{
			$nref = ceil(100*$vaf) - int(rand(3));
		}
		my $ref = 100-$nref;
		my $m = $n+1;
		$n++;
		print MAF "$vaf\_C_$status\t$chrom\t$n\t$m\tA\tC\t$ref\t$nref\t0\t0\n";

	}
}

print SEG "sample\tchromosome\tstart\tend\tLOH_flag\tBAF\tlogratio\n";
for my $chr (keys %genome_states){
	if($genome_states{$chr} eq "LOH"){
		#just print 10 segments with the same BAF to mimic a series of LOH germline SNVs in the region
		
		for my $c_freq (keys %clones){
			for my $n(1..10){
				my $vaf = 1- ($c_freq + (1-$c_freq)/2);
				my $status = "nLOH";
				my $nref = ceil(100*$vaf);
				my $ref = 100-$nref;
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
