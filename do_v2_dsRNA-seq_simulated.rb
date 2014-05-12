#! /usr/bin/ruby

# test dsRNA-seq model with simulated data

# list of RNAs
#loci = ["SNORD36A"]
#loci = ["U1_snRNA", "U3_snRNA", "U5_snRNA", "U15_snoRNA", "U22_snoRNA", "U97_snoRNA", "hsa-let-7a-1", "hsa-mir-17", "5S_rRNA"]
#loci = ["U15_snoRNA", "U22_snoRNA", "U97_snoRNA", "5S_rRNA"]
loci = ["ASH1_validated_region", "HAC1_validated_region", "NM_003234_validated_region", "NM_005080_validated_region", "NM_016332_validated_region"]
enzymes = ["v2_dsRNA-seq_simulated"]
#enzymes = ["A3_dsRNA", "A4_dsRNA", "A5_dsRNA", "A6_dsRNA", "B1_dsRNA", "B2_dsRNA", "C1_dsRNA", "C2_dsRNA", "D1_dsRNA", "D2_dsRNA", "E1_dsRNA", "E2_dsRNA", "F1_dsRNA", "G1_dsRNA"]
uds = "A:0.01,C:0.01,T:0.01,G:0.01"
vds = "A:0.02,C:0.02,T:0.02,G:0.02"
num_reads = 100000
ncuts = "1,2,3"

## New analyses with v2b likelihood model
# estimate structure from real data
enzymes.each { |enzyme|
	# copy base files
	if !File.directory?("/Data03/fli/dsRNA-seq/#{enzyme}/")
		system("mkdir /Data03/fli/dsRNA-seq/#{enzyme}/")
	end
	system("cp /Data03/fli/dsRNA-seq/base_yeast/*1* /Data03/fli/dsRNA-seq/#{enzyme}/")
	system("cp /Data03/fli/dsRNA-seq/base_human/* /Data03/fli/dsRNA-seq/#{enzyme}/")
	#	generate simulated reads
	loci.each { |locus|
		system("./generate_simulated_reads_v2b.R /Data03/fli/dsRNA-seq/#{enzyme}/#{locus}.known.rnafold #{uds} #{vds} 15 40 #{ncuts} #{num_reads} /Data03/fli/dsRNA-seq/#{enzyme}/#{locus}.ds.collapsed.reads")
	}

	# estimate digestion rates based on reads
#	system("ruby /Data02/fli/code/RNA-seq-fold/make_endpoint_files_by_dir.rb /Data03/fli/dsRNA-seq/#{enzyme}/ --log2")
	system("ruby /Data02/fli/code/RNA-seq-fold/make_endpoint_files_by_dir.rb /Data03/fli/dsRNA-seq/#{enzyme}/ --type ds.collapsed --5p")
	system("ruby /Data02/fli/code/RNA-seq-fold/calc_pernucl_digestion_rates.rb /Data03/fli/dsRNA-seq/#{enzyme}/ 0.01,0.01,0.01,0.01 ds.collapsed.endpoints > /Data03/fli/dsRNA-seq/#{enzyme}/#{enzyme}.ds.v2b.rates")

	loci.each { |locus|
		pid = fork do
			exec("perl /Data02/fli/code/RNASeqFold-2.0/RNASeqFold.pl -rmin 15 -rmax 40 -n 100000 -burn_in 10000 -samp 100 -t sd -m 2 -numCPU 1 -outdir /Data03/fli/dsRNA-seq/#{enzyme}/ /Data03/fli/dsRNA-seq/#{enzyme}/#{locus}.ds.collapsed.reads /Data03/fli/dsRNA-seq/#{enzyme}/#{enzyme}.ds.v2b.rates /Data03/fli/dsRNA-seq/#{enzyme}/#{locus}.rnafold")
		end
#		system("renice +11 #{pid}")
	}
	Process.waitall
	
	loci.each { |locus|
#		system("mv /Data03/fli/dsRNA-seq/#{enzyme}/#{locus}.posterior.001.txt /Data03/fli/dsRNA-seq/#{enzyme}/#{locus}.s.posterior.txt")
#		system("mv /Data03/fli/dsRNA-seq/#{enzyme}/#{locus}.mcmc.001.txt /Data03/fli/dsRNA-seq/#{enzyme}/#{locus}.s.mcmc.txt")
#		system("mv /Data03/fli/dsRNA-seq/#{enzyme}/#{locus}.mcmc.001.log /Data03/fli/dsRNA-seq/#{enzyme}/#{locus}.s.mcmc.log")
		
		system("ln -s /Data03/fli/dsRNA-seq/#{enzyme}/#{locus}.posterior.001.txt /Data03/fli/dsRNA-seq/#{enzyme}/#{locus}.s.posterior.txt")
		system("ruby get_posterior_distribution_dinucl.rb /Data03/fli/dsRNA-seq/#{enzyme}/#{locus}.mcmc.001.txt s 10000 100 > /Data03/fli/dsRNA-seq/#{enzyme}/#{locus}.posterior.001.txt")
		system("ruby get_posterior_distribution_dinucl.rb /Data03/fli/dsRNA-seq/#{enzyme}/#{locus}.mcmc.001.txt d 10000 100 > /Data03/fli/dsRNA-seq/#{enzyme}/#{locus}.parameters.001.txt")
		system("ruby annotate_svg_plot.rb /Data03/fli/dsRNA-seq/#{enzyme}/#{locus}.svg /Data03/fli/dsRNA-seq/#{enzyme}/#{locus}.posterior.001.txt -a posterior -s -m 10 -e -z -i #{locus} -colorscheme blue-red > /Data03/fli/dsRNA-seq/#{enzyme}/#{locus}.s.posterior.svg")
##		system("./summ_MCMC_acceptance_rate.R /Data03/fli/dsRNA-seq/#{enzyme}/#{locus}.s.mcmc.log /Data03/fli/dsRNA-seq/#{enzyme}/#{locus}.s.mcmc.pdf")
##		
##		# plot endpoints on MCMC posterior-informed structure
##		system("ruby posterior2constraint.rb /Data03/fli/dsRNA-seq/#{enzyme}/#{locus}.s.posterior.txt /Data03/fli/dsRNA-seq/#{enzyme}/#{locus}.seq 0.45 0.55 > /Data03/fli/dsRNA-seq/#{enzyme}/#{locus}.s.posterior.constraint")
##		system("RNAfold --noPS -C < /Data03/fli/dsRNA-seq/#{enzyme}/#{locus}.s.posterior.constraint > /Data03/fli/dsRNA-seq/#{enzyme}/#{locus}.s.posterior.constrained.rnafold")
##		system("RNAplot -o svg < /Data03/fli/dsRNA-seq/#{enzyme}/#{locus}.s.posterior.constrained.rnafold")
##		system("mv rna.svg /Data03/fli/dsRNA-seq/#{enzyme}/#{locus}.s.posterior.constrained.svg")
##		system("ruby annotate_svg_plot.rb /Data03/fli/dsRNA-seq/#{enzyme}/#{locus}.s.posterior.constrained.svg /Data03/fli/dsRNA-seq/#{enzyme}/#{locus}.ds.collapsed.endpoints -a ds_endpoints -s -m 10 -e -i #{locus} -colorscheme blue-red > /Data03/fli/dsRNA-seq/#{enzyme}/#{locus}.ds.collapsed.endpoints.s.posterior.constrained.svg")
##		
##		# plot endpoints on gold standard structure
		system("ruby annotate_svg_plot.rb /Data03/fli/dsRNA-seq/#{enzyme}/#{locus}.svg /Data03/fli/dsRNA-seq/#{enzyme}/#{locus}.ds.collapsed.endpoints -a ds_endpoints -s -m 10 -e -i #{locus} -colorscheme blue-red > /Data03/fli/dsRNA-seq/#{enzyme}/#{locus}.ds.collapsed.endpoints.svg")
		
		
	}

}


