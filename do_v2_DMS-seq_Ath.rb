#! /usr/bin/ruby


# list of RNAs
#loci = ["SNORD36A"]
#loci = ["U1_snRNA", "U3_snRNA", "U5_snRNA", "U15_snoRNA", "U22_snoRNA", "U97_snoRNA", "hsa-let-7a-1", "hsa-mir-17", "5S_rRNA"]
#loci = ["U15_snoRNA", "U22_snoRNA", "U97_snoRNA", "5S_rRNA"]
loci = ["AT1G29930_validated_region"]
enzymes = ["v2_Col-0-Leaf"]
#enzymes = ["A3_dsRNA", "A4_dsRNA", "A5_dsRNA", "A6_dsRNA", "B1_dsRNA", "B2_dsRNA", "C1_dsRNA", "C2_dsRNA", "D1_dsRNA", "D2_dsRNA", "E1_dsRNA", "E2_dsRNA", "F1_dsRNA", "G1_dsRNA"]
uds = "A:0.01,C:0.01,T:0.01,G:0.01"

## New analyses with v2b likelihood model
# estimate structure from real data
enzymes.each { |enzyme|
	
	# copy base files
	if !File.directory?("/Data03/fli/DMS-seq/#{enzyme}/")
		system("mkdir /Data03/fli/DMS-seq/#{enzyme}/")
	end
	system("cp /Data03/fli/DMS-seq/base_ath/* /Data03/fli/DMS-seq/#{enzyme}/")
	system("cp /Data03/fli/SRP027216/*.reads /Data03/fli/DMS-seq/#{enzyme}/")

	# estimate digestion rates based on reads
	system("ruby /Data02/fli/code/RNA-seq-fold/make_endpoint_files_by_dir.rb /Data03/fli/DMS-seq/#{enzyme}/ --type DMS.collapsed")
	system("ruby /Data02/fli/code/RNA-seq-fold/calc_pernucl_digestion_rates.rb /Data03/fli/DMS-seq/#{enzyme}/ 0.01,0.01,0.01,0.01 DMS.collapsed.endpoints > /Data03/fli/DMS-seq/#{enzyme}/#{enzyme}.DMS.v2b.rates")

	loci.each { |locus|
		pid = fork do
			exec("perl RNASeqFold.pl -rmin 15 -rmax 40 -n 100000 -burn_in 10000 -samp 100 -t s -m 2 -numCPU 1 -outdir /Data03/fli/DMS-seq/#{enzyme}/ /Data03/fli/DMS-seq/#{enzyme}/#{locus}.DMS.collapsed.reads /Data03/fli/DMS-seq/#{enzyme}/#{enzyme}.DMS.v2b.rates /Data03/fli/DMS-seq/#{enzyme}/#{locus}.rnafold")
		end
	}
	Process.waitall
	
##	loci.each { |locus|
###		system("ln -s /Data03/fli/DMS-seq/#{enzyme}/#{locus}.posterior.001.txt /Data03/fli/DMS-seq/#{enzyme}/#{locus}.s.posterior.txt")
###		system("ruby get_posterior_distribution_dinucl.rb /Data03/fli/DMS-seq/#{enzyme}/#{locus}.s.mcmc.txt s 10000 100 > /Data03/fli/DMS-seq/#{enzyme}/#{locus}.s.posterior.txt")
###		system("mv /Data03/fli/DMS-seq/#{enzyme}/#{locus}.posterior.001.txt /Data03/fli/DMS-seq/#{enzyme}/#{locus}.s.posterior.txt")
###		system("mv /Data03/fli/DMS-seq/#{enzyme}/#{locus}.mcmc.001.txt /Data03/fli/DMS-seq/#{enzyme}/#{locus}.s.mcmc.txt")
###		system("mv /Data03/fli/DMS-seq/#{enzyme}/#{locus}.mcmc.001.log /Data03/fli/DMS-seq/#{enzyme}/#{locus}.s.mcmc.log")
##		
##		system("ruby annotate_svg_plot.rb /Data03/fli/DMS-seq/#{enzyme}/#{locus}.svg /Data03/fli/DMS-seq/#{enzyme}/#{locus}.posterior.001.txt -a posterior -s -m 10 -e -z -i #{locus} -colorscheme blue-red > /Data03/fli/DMS-seq/#{enzyme}/#{locus}.s.posterior.svg")
###		system("./summ_MCMC_acceptance_rate.R /Data03/fli/DMS-seq/#{enzyme}/#{locus}.s.mcmc.log /Data03/fli/DMS-seq/#{enzyme}/#{locus}.s.mcmc.pdf")
##		
###		# plot endpoints on MCMC posterior-informed structure
###		system("ruby posterior2constraint.rb /Data03/fli/DMS-seq/#{enzyme}/#{locus}.s.posterior.txt /Data03/fli/DMS-seq/#{enzyme}/#{locus}.seq 0.45 0.55 > /Data03/fli/DMS-seq/#{enzyme}/#{locus}.s.posterior.constraint")
###		system("RNAfold --noPS -C < /Data03/fli/DMS-seq/#{enzyme}/#{locus}.s.posterior.constraint > /Data03/fli/DMS-seq/#{enzyme}/#{locus}.s.posterior.constrained.rnafold")
###		system("RNAplot -o svg < /Data03/fli/DMS-seq/#{enzyme}/#{locus}.s.posterior.constrained.rnafold")
###		system("mv rna.svg /Data03/fli/DMS-seq/#{enzyme}/#{locus}.s.posterior.constrained.svg")
###		system("ruby annotate_svg_plot.rb /Data03/fli/DMS-seq/#{enzyme}/#{locus}.s.posterior.constrained.svg /Data03/fli/DMS-seq/#{enzyme}/#{locus}.DMS.collapsed.endpoints -a DMS_endpoints -s -m 10 -e -i #{locus} -colorscheme blue-red > /Data03/fli/DMS-seq/#{enzyme}/#{locus}.DMS.collapsed.endpoints.s.posterior.constrained.svg")
##		
##		# plot endpoints on gold standard structure
##		system("ruby annotate_svg_plot.rb /Data03/fli/DMS-seq/#{enzyme}/#{locus}.svg /Data03/fli/DMS-seq/#{enzyme}/#{locus}.DMS.collapsed.endpoints -a DMS_endpoints -log -s -m 10 -e -z -i #{locus} -colorscheme blue-red > /Data03/fli/DMS-seq/#{enzyme}/#{locus}.DMS.collapsed.endpoints.svg")
##		
##		
##	}
}
