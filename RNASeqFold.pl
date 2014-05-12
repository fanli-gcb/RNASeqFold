#!/usr/bin/perl -w
# Wrapper to run RNA-seq-fold and extract posterior
use strict;
use Getopt::Long qw( :config posix_default no_ignore_case );
use Pod::Usage;
use IPC::Open3;
use Symbol qw(gensym);
use IO::File;
use Forks::Super;
$Forks::Super::MAX_PROC = 1;
$Forks::Super::ON_BUSY = 'block';

# flags and variables
my $readsFile = "";
my $structureFile = "";
my $ratesFile = "";
my $rmin = 40;
my $rmax = 80;
my $n = 100000;
my $numCPU = 1;
my $numStartingPoints = 0;
my $burnIn = 10000;
my $samp = 100;
my $t = "s";
my $model = 1;
my $mcConstant = 999.0;
my $likelihoodOnly = 0;
my $outDir = "RNASeqFold_output";
my $man = 0;
my $help = 0;
my $id = "";

my $doStructure = 0;
my $doParameters = 0;

# Parse options and print usage if there is a syntax error,
# or if usage was explicitly requested.
GetOptions('rmin=i' => \$rmin, 'rmax=i' => \$rmax, 'n=i' => \$n, 'numCPU=i' => \$numCPU, 'burn_in=i' => \$burnIn, 'samp=i' => \$samp, 't=s' => \&type_handler, 'm=i' => \$model, 'mc=f' => \$mcConstant, 'outdir=s' => \$outDir, 'likelihood_only' => \$likelihoodOnly, 'help|?' => \$help, 'man' => \$man) or pod2usage(-verbose => 1);

sub type_handler {
	$t = $_[1];
#	die("Invalid type $t") if ($t =~ /[^sd]/);
	die("Invalid type $t\n") unless (($t eq "s") || ($t eq "d") || ($t eq "sd") || ($t eq "ds"));
	if ($t =~ /s/) {
		$doStructure = 1;
	}
	if ($t =~ /d/) {
		$doParameters = 1;
	}
}

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

my $remaining = @ARGV;
if ($remaining < 3) {
	printf STDERR "\nERROR: Missing required arguments(s).\n\n";
	pod2usage(-verbose => 1);
}
$readsFile = shift;
$ratesFile = shift;
$structureFile = shift;
if ($structureFile =~ /\/([^\/]+)\.rnafold/) {
	$id = $1;
}

if (@ARGV) {
	printf STDERR "\nWarning: Extra arguments/options detected. Exiting!\n\n";
	pod2usage(-verbose => 1);
}

$outDir =~ s/\/$//;
if (!(-d "$outDir")) {
	mkdir($outDir);
}
$Forks::Super::MAX_PROC = $numCPU;

# separate starting points
open(STRUCT, "$structureFile") || die "Unable to read from $structureFile: $!\n";
while (my $line = <STRUCT>) {
	$numStartingPoints++;
	chomp($line);
	my $line2 = <STRUCT>; chomp($line2);
	open(TOUT, sprintf(">%s/%s\.starting\.%.3d\.rnafold", $outDir, $id, $numStartingPoints));
	print TOUT "$line\n$line2\n";
	close(TOUT);
}
close(STRUCT);

# run MCMC and get posterior if all checks passed
print STDERR sprintf("Starting MCMC with %d CPUs ...\n", $numCPU);
my $i = 1;
my @pids;
my $pids = @pids;
for (my $i=1; $i<=$numStartingPoints; $i++) {
	my $starting_fn = sprintf("%s/%s\.starting\.%.3d\.rnafold", $outDir, $id, $i);
	my $out_mcmc_fn = sprintf("%s/%s\.mcmc\.%.3d\.txt", $outDir, $id, $i);
	my $out_err_fn = sprintf("%s/%s\.mcmc\.%.3d\.err", $outDir, $id, $i);
	my $out_log_fn = sprintf("%s/%s\.mcmc\.%.3d\.log", $outDir, $id, $i);
	
	my $likelihoodOnlyFlag = ($likelihoodOnly) ? "-likelihood_only" : "";
#	print STDERR "start RNASeqFold with $starting_fn $readsFile $ratesFile\n";
	my $cmd = "./RNASeqFold $starting_fn $readsFile $ratesFile -rmin $rmin -rmax $rmax -r 0 -n $n -t $t -m $model -mc $mcConstant -log $out_log_fn $likelihoodOnlyFlag -out $out_mcmc_fn";
	
	open(ERR, ">$out_err_fn") || die "Unable to write $out_err_fn: $!\n";
	my $out_err = "";
	my $pid = Forks::Super::fork { cmd => $cmd, stderr => \$out_err };
	print STDERR "Running MCMC for sequence $i ...\n";
	print ERR "$out_err";
	close(ERR);
}
waitall;

if (!($likelihoodOnly)) {
	for (my $i=1; $i<=$numStartingPoints; $i++) {
		print STDERR "Extracting posterior for sequence $i ...\n";
		my $out_mcmc_fn = sprintf("%s/%s\.mcmc\.%.3d\.txt", $outDir, $id, $i);
		
		if ($doStructure) {
			my $out_err_fn = sprintf("%s/%s\.posterior\.%.3d\.err", $outDir, $id, $i);
			my $out_posterior_fn = sprintf("%s/%s\.posterior\.%.3d\.txt", $outDir, $id, $i);
			my $cmd = "ruby get_posterior_distribution_dinucl.rb $out_mcmc_fn s $burnIn $samp";
	
			open(OUT, ">$out_posterior_fn") || die "Unable to write $out_posterior_fn";
			open(ERR, ">$out_err_fn") || die "Unable to write $out_err_fn";
			my $pid = open3(gensym, ">&OUT", ">&ERR", $cmd);
			waitpid($pid, 0);
			seek $_, 0, 0 for \*OUT, \*ERR;
			close(OUT);
			close(ERR);
		}
		if ($doParameters) {
			my $out_err_fn = sprintf("%s/%s\.parameters\.%.3d\.err", $outDir, $id, $i);
			my $out_posterior_fn = sprintf("%s/%s\.parameters\.%.3d\.txt", $outDir, $id, $i);
			my $cmd = "ruby get_posterior_distribution_dinucl.rb $out_mcmc_fn d $burnIn $samp";
	
			open(OUT, ">$out_posterior_fn") || die "Unable to write $out_posterior_fn";
			open(ERR, ">$out_err_fn") || die "Unable to write $out_err_fn";
			my $pid = open3(gensym, ">&OUT", ">&ERR", $cmd);
			waitpid($pid, 0);
			seek $_, 0, 0 for \*OUT, \*ERR;
			close(OUT);
			close(ERR);
		}
	}
	waitall;
}
## Deprecated open3 version
##print STDERR "Running MCMC in non-batch mode ...\n";
##for (my $i=1; $i<=$numStartingPoints; $i++) {
##	my $starting_fn = sprintf("%s/starting\.%.3d\.rnafold", $outDir, $i);
##	my $out_mcmc_fn = sprintf("%s/mcmc\.%.3d\.txt", $outDir, $i);
##	my $out_err_fn = sprintf("%s/mcmc\.%.3d\.err", $outDir, $i);
##	my $out_log_fn = sprintf("%s/mcmc\.%.3d\.log", $outDir, $i);
##	my $out_posterior_fn = sprintf("%s/posterior\.%.3d\.txt", $outDir, $i);
##	my $out_posterior_err_fn = sprintf("%s/posterior\.%.3d\.err", $outDir, $i);
##	my $cmd = "./RNASeqFold -s $starting_fn $dsArg $ssArg $ratesFile -ncuts $ncuts -rmin $rmin -rmax $rmax -r 0 -n $n -t $t -log $out_log_fn -out $out_mcmc_fn";
##	open(ERR, ">$out_err_fn") || die "Unable to write $out_err_fn";
##	my $pid = open3(gensym, ">&STDOUT", ">&ERR", $cmd);
##	waitpid($pid, 0);
##	seek $_, 0, 0 for \*ERR;
##	close(ERR);

##	print STDERR "Extracting posterior in non-batch mode ...\n";
##	$cmd = "ruby get_posterior_distribution_dinucl.rb $outDir/mcmc.txt $t $burnIn $samp";
##	open(OUT, ">$out_posterior_fn") || die "Unable to write $out_posterior_fn";
##	open(ERR, ">$out_err_fn") || die "Unable to write $out_err_fn";
##	$pid = open3(gensym, ">&OUT", ">&ERR", $cmd);
##	waitpid($pid, 0);
##	seek $_, 0, 0 for \*OUT, \*ERR;
##	close(OUT);
##	close(ERR);
##}


#print "./optimize_structure2b -s $structureFile $dsArg $ssArg $ratesFile -ncuts $ncuts -rmin $rmin -rmax $rmax -r 0 -n $n -t $t -log $outDir/mcmc.log > $outDir/mcmc.txt\n";
#print "ruby get_posterior_distribution_dinucl.rb $outDir/mcmc.txt $t $burnIn $samp > $outDir/posterior.txt\n";

__END__
 
=head1 NAME
 
RNASeqFold - Bayesian Markov chain Monte Carlo inference of RNA secondary structures.
 
=head1 SYNOPSIS
 
RNASeqFold [options] reads_file rates_file structure_file
 
=head1 ARGUMENTS

=over 8
 
=item B<reads_file>

File containing either DMS-seq or dsRNA/ssNRA-seq reads.

=item B<rates_file>

File containing enzyme digestion rates. The first line should contain four comma-delimited rates corresponding to digestion probabilities 3' of paired [A,C,U,G] nucleotides.
The second line should contain rates 3' of unpaired [A,C,U,G] positions.

=item B<structure_file>

File containing the sequence and structure of the RNA of interest. The provided dot-paren structure is used as the initial estimate.

=back

=head1 OPTIONS

=over 8

=item B<-rmin> INT

Minimum read length. Fragments shorter than this will not be included in likelihood calcuations. Defaults to 40.

=item B<-rmax> INT

Maximum read length. Fragments longer than this will not be included in likelihood calcuations. Defaults to 80.

=item B<-n> INT

Number of MCMC iterations to run. Defaults to 100000.

=item B<-numCPU> INT

Maximum number of CPUs to use. Defaults to 1.

=item B<-burn_in> INT

Number of MCMC iterations to discard as burn-in when computing posterior. Defaults to 10000.

=item B<-samp> INT

Sampling frequency of MCMC iterations after burn-in. Defaults to 100.

=item B<-t> [sd]

Selects which parameter to optimize. B<'s'> indicates structure and B<'d'> indicates digestion rates. Defaults to B<'s'> for structure estimation.

=item B<-m> INT

Select likelihood model to use (1=dsRNA/ssRNA-seq, 2=DMS-seq)

=item B<-mc> DOUBLE

Proportionality constant used in Metropolis-Hastings acceptance criterion (Default = 999.0)

=item B<-likelihood_only

Only compute the initial log-likelihood (do not perform MCMC iterations)

=item B<-outdir dir>

Directory to output results. Defaults to 'RNASeqFold_output/'.

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

Estimate base pairing posteriors for a given RNA sequence based on the distribution of read fragments along the locus.

This program requires a file containing the sequence and initial secondary structure of an RNA of interest, a file containing DMS-seq or dsRNA-seq/ssRNA-seq reads, and a file containing the initial estimates of enzyme digestion rates. A Bayesian MCMC algorithm is then used to estimate the base pairing posterior that best fits the observed sequencing reads. The results and intermediate files are written to a directory (RNASeqFold_output/ by default).

The output file posterior.txt contains the base pairing posteriors at each nucleotide position, one entry per line.

=cut

