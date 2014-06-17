#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <tr1/unordered_map>
#include <sstream>
#include <climits>
#include <list>
#include <math.h>
#include <limits.h>
#include "boost/multi_array.hpp"
#include <boost/lexical_cast.hpp>

using namespace std;

// split a std::string -> vector of std::strings
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
  std::stringstream ss(s);
  std::string item;
  while(std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}
std::vector<std::string> split(const std::string &s, char delim) {
  std::vector<std::string> elems;
  return split(s, delim, elems);
}

// split a std::string -> vector of doubles
std::vector<double> &split2double(const std::string &s, char delim, std::vector<double> &elems) {
  std::stringstream ss(s);
  std::string item;
  while(std::getline(ss, item, delim)) {
    elems.push_back(atof(item.c_str()));
  }
  return elems;
}
std::vector<double> split2double(const std::string &s, char delim) {
  std::vector<double> elems;
  return split2double(s, delim, elems);
}

// split a std::string -> vector of ints
std::vector<int> &split2int(const std::string &s, char delim, std::vector<int> &elems) {
  std::stringstream ss(s);
  std::string item;
  while(std::getline(ss, item, delim)) {
    elems.push_back(atoi(item.c_str()));
  }
  return elems;
}
std::vector<int> split2int(const std::string &s, char delim) {
  std::vector<int> elems;
  return split2int(s, delim, elems);
}

// define a Read struct
// these are 0-indexed reads, not 1-indexed!
struct Read {
	std::string id;
	int start;
	int end;
	int nocc;
	int count;
	int len;
	Read() {}
	// read a *.reads row
	Read(std::string &line) {
		std::string junk;
		int extend;
    istringstream is(line);
		is >> junk;
		if (!junk.empty()) {
			std::vector<std::string> elems = split(junk, '@');
			id = elems[0];
			count = atoi(elems[2].c_str());
			len = atoi(elems[4].c_str());
			is >> start;
			is >> end;
			is >> nocc;
		}
	}
};

struct Pair {
	int p1;
	int p2;
	Pair(int a, int b) : p1(a), p2(b) {}
};

struct Move {
	int type;
	int p1;
	int p2;
	Move(int a, int b, int c) : type(a), p1(b), p2(c) {}
};

// ====================================================================================================
// FUNCTION DEFINITIONS
bool canPair(char a, char b);
bool checkStructure(std::string structure, std::string seq);
std::vector<double> updateDigestionParameter(std::vector<double> p);
std::vector<double> get_per_site_digestion_rates(std::vector<double> &per_nuc_digestion_rates, std::string seq);
bool checkDigestionParameters(std::vector<double> &u, std::vector<double> &v, double r);
void writeToLog(ofstream &fp, std::string type, int result);
void writeToLog(ofstream &fp, std::string type, int result, double currentLL, double newLL, double q1, double q2, double mcConstant, std::vector<int> &moveCounts);
double calculateCloneNormalizationFactor(list<Read> &reads, int lcut, int rcut);
std::string vec2string(vector<double> &values, std::string sep);
std::string intvec2string(std::vector<int> &values, std::string sep);

int estimateFragmentSize(list<Read> &reads);
double compute_LL (std::string structure, std::list<Read> &reads, std::vector<double> &u, std::vector<double> &v, double r, int lcut, int rcut, int model);
std::vector<Move> getValidMoves(std::string structure, std::string seq);
int numValidMoves(std::string structure, std::string seq);
std::vector<int> getMoveCounts(std::string structure, std::string seq);
std::string randomMove(std::string structure, std::string seq, bool weightMFE);
std::string makeMove(std::string structure, std::string seq, Move m);
std::vector<double> getMoveWeights(std::vector<Move> moves, std::string seq);

unsigned long mix(unsigned long a, unsigned long b, unsigned long c);
double drand(double min, double max);
int irand(int min, int max);
int weighted_irand(vector<double> &weights);

void doMCMC(std::string structure, std::string seq, list<Read> &reads, vector<double> &u, vector<double> &v, double r, int lcut, int rcut, int model, double mcConstant, int num_iters, bool doStructure, bool doDigestion, bool doFragmentation, bool weightMFE, bool likelihoodOnly, bool doLog, std::string log_fn, std::string out_fn);
// ====================================================================================================

/*
uses Gibbs sampling
NOTE: (06/17/2013) - New v2 likelihood model that is conditioned on N digestion events per molecule, not the "infinite digestion" assumption of previous model with site independence
NOTE: (07/26/2013) - 2 digestion rates [u,v] per nucleotide (total of 2x4=8 digestion rate parameters)
NOTE: (02/12/2014) - Using old likelihood model (independent reads) for speedup
INPUTS:
	structfn - *.rnafold file containing some initial starting structure for the transcript (probably RNAfold prediction); should also contain the sequence
	readsfn - file containing DMS-seq or dsRNA/ssRNA-seq reads
	ratesfn - file containing comma-separated digestion rates for paired and unpaired positions (should be in order [A,C,T,G])
	type - any combination of [sdfa] - structure, digestion, fragmentation, all (overrides other flags, obviously)
	rmin - lower range of size restriction
	rmax - upper range of size restriction
	num_iters - number of iterations to run MCMC
	model - variant of likelihood model (1 = default: 5' and 3' digestion, 2 = yeast PARS: 5' digestion and 3' random fragmentation, 3 = mouse FragSeq: 5' and 3' digestion at 3+ runs of loops)
	r - random fragmentation rate
	weight_mfe - weight move probabilities by their base pair MFE (e.g. GC > AU > GU)
	logfile - file to which the acceptance/rejection decisions are logged
OUTPUTS:
	print to stdout
*/

int main(int argc, char **argv) {
	bool doStructure = false;
	bool doDigestion = false;
	bool doFragmentation = false;
	bool weightMFE = false;
	bool doLog = false;
	bool likelihoodOnly = false;
	std::vector<double> u;
	std::vector<double> v;
	double r = 0.05; // default fragmentation probability at each nucleotide
	int lcut = 20;
	int rcut = 35;
	std::string typestr = "sd";
	int num_iters = 1000000;
	int model = 1;
	double mcConstant = 999.0; // the higher this is, the less likely we are to accept under Metropolis-Hastings
	std::string readsfn, ratesfn, structfn, log_fn, out_fn;
	std::string usage = "USAGE: ./RNASeqFold struct_fn fn rate_fn [-rmin rmin] [-rmax rmax] [-t type] [-n num_iters] [-mc mcConstant] [-m model] [-r frag_param] [-weight_mfe] [-log log_fn] [-out out_fn] [-likelihood only]\n";
	
	for (int i=1; i<argc; i++) {
		std::string arg(argv[i]);
		if (arg.substr(0,1).compare("-")!=0) {
			if (argc < i+5) {
				cerr << usage;
				return(-1);
			}
			structfn = argv[i];
			readsfn = argv[++i];
			ratesfn = argv[++i];
		}
		else if (arg.compare("-rmin")==0) {
			if (argc < i+2) {
				cerr << usage;
				return(-1);
			}
			lcut = atoi(argv[++i]);
		}
		else if (arg.compare("-rmax")==0) {
			if (argc < i+2) {
				cerr << usage;
				return(-1);
			}
			rcut = atoi(argv[++i]);
		}
		else if (arg.compare("-t")==0) {
			if (argc < i+2) {
				cerr << usage;
				return(-1);
			}
			typestr = argv[++i];
		}
		else if (arg.compare("-n")==0) {
			if (argc < i+2) {
				cerr << usage;
				return(-1);
			}
			num_iters = atoi(argv[++i]);
		}
		else if (arg.compare("-m")==0) {
			if (argc < i+2) {
				cerr << usage;
				return(-1);
			}
			model = atoi(argv[++i]);
			if (!((model==1) || (model==2) || (model==3))) {
				cerr << "Invalid model - must be 1 (default dsRNA-seq/ssRNA-seq), 2 (DMS-seq), or 3 (FragSeq)\n";
				return(-1);
			}
		}
		else if (arg.compare("-mc")==0) {
			if (argc < i+2) {
				cerr << usage;
				return(-1);
			}
			mcConstant = atof(argv[++i]);
		}
		else if (arg.compare("-r")==0) {
			if (argc < i+2) {
				cerr << usage;
				return(-1);
			}
			r = atof(argv[++i]);
			if ((r<0) || (r>1)) {
				cerr << "Invalid frag_param - must be between 0 and 1\n";
				return(-1);
			}
		}
		else if (arg.compare("-weight_mfe")==0) {
			weightMFE = true;
		}
		else if (arg.compare("-log")==0) {
			if (argc < i+2) {
				cerr << usage;
				return(-1);
			}
			log_fn = argv[++i];
			doLog = true;
		}
		else if (arg.compare("-out")==0) {
			if (argc < i+2) {
				cerr << usage;
				return(-1);
			}
			out_fn = argv[++i];
		}
		else if (arg.compare("-likelihood_only")==0) {
			likelihoodOnly = true;
		}
		else {
			cerr << "unrecognized option " << arg << endl;
			return(-1);
		}
	}
	// check arguments
	if (structfn.empty()) {
		cerr << "requires struct_fn\n";
		return(-1);
	}
	// parse typestr
	for (int i=0; i<typestr.length(); i++) {
		std::string tmp = typestr.substr(i,1);
		if (tmp.compare("s")==0)
			doStructure = true;
		else if (tmp.compare("d")==0)
			doDigestion = true;
		else if (tmp.compare("f")==0)
			doFragmentation = true;
		else if (tmp.compare("a")==0) {
			doStructure = true;
			doDigestion = true;
			doFragmentation = true;
		}
		else {
			cerr << "Invalid type string " << typestr << endl;
			return(-1);
		}
	}

  // open files for reading
	ifstream reads_is, rates_is;
	reads_is.open(readsfn.c_str());
	if (!reads_is.is_open()) {
	  cerr << "Error: Could not open file " << readsfn << " for reading." << endl;
	  return(1);
	}
	rates_is.open(ratesfn.c_str());
	if (!rates_is.is_open()) {
	  cerr << "Error: Could not open file " << ratesfn << " for reading." << endl;
	  return(1);
	}
  ifstream struct_is(structfn.c_str());
  if (!struct_is.is_open()) {
    cerr << "Error: Could not open file " << structfn << " for reading." << endl;
    return(1);
  }
	
	// variable declaration
	std::string seq;
	std::string structure;
	std::string line;
	list<Read> reads;
	unsigned long seed = mix(clock(), time(NULL), getpid());
	srand(seed);
	
	// read in the sequence, structure, and reads
	getline(struct_is, seq);
	getline(struct_is, structure);
	
	size_t found;
	found = seq.find("T");
	while (found != -1) {
		seq.replace(found,1,"U");
		found = seq.find("T");
	}	
	found = structure.find("\t");
	if (found != -1) // replace tabs with space
		structure.replace(found,1," ");
	found = structure.find(" ");
	structure.erase(found);
	int structure_length = structure.length();
	
  while(!reads_is.eof()) {
		getline(reads_is, line);
		if (line.empty()) {
			continue;
		}
		Read read(line);
		reads.push_back(read);
	}
	
	// read in digestion rates
	getline(rates_is, line);
	u = split2double(line, ',');
	getline(rates_is, line);
	v = split2double(line, ',');
	
//	cout << fixed;
	doMCMC(structure, seq, reads, u, v, r, lcut, rcut, model, mcConstant, num_iters, doStructure, doDigestion, doFragmentation, weightMFE, likelihoodOnly, doLog, log_fn, out_fn);
  
}





// ====================================================================================================
// FUNCTION IMPLEMENTATIONS

// Return true if a can base-pair with b, false otherwise
bool canPair(char a, char b) {
	if ((a == 'A') && (b == 'U'))
		return true;
	else if ((a == 'U') && ((b == 'A') || (b == 'G')))
		return true;
	else if ((a == 'G') && ((b == 'C') || (b == 'U')))
		return true;
	else if ((a == 'C') && (b == 'G'))
		return true;
	else
		return false;
}

// Check that the given structure is properly nested, only contains valid base pairs, and doesn't violate steric constraints
bool checkStructure(std::string structure, std::string seq) {
	std::vector<int> stack;
	int count = 0;
	bool good;
	for (int i=0; i<structure.length(); i++) {
		if (structure[i] == '.') {
			// do nothing
		}
		else if (structure[i] == '(') {
			stack.push_back(i);
			count++;
		}
		else if (structure[i] == ')') {
			// nesting check
			if (count <= 0) {
				return false;
			}
			count--;
			// steric check
			if (abs(i-stack.back()) < 3) {
				return false;
			}
			// base pairing check
			good = canPair(seq[stack.back()], seq[i]);
			if (!good)
				return false;
			stack.pop_back();

		}
		else if (structure[i] == 'x') {
			// do nothing - these are fixed unpaired flanking regions
		}
		else {
			cerr << "Invalid character found in structure " << structure << endl;
			exit(1);
		}
	}
	return true;
	
}

// Check that the digestion parameters are within valid bounds
bool checkDigestionParameters(vector<double> &u, vector<double> &v, double r) {

	for (vector<double>::iterator it=u.begin(); it!=u.end(); ++it) {
		double value = *it;
		if (((value <= 0.0) || (value >= 1.0)) && (value != -1))
			return false;
	}
	for (vector<double>::iterator it=v.begin(); it!=v.end(); ++it) {
		double value = *it;
		if (((value <= 0.0) || (value >= 1.0)) && (value != -1))
			return false;
	}
	if (((r < 0.0) || (r >= 1.0)) && (r != -1))
		return false;
		
	// check that u < v (within some bound)
	for (int i=0; i<u.size(); i++) {
		if (v[i] < u[i]) {
			return false;
		}
	}
	
		
	return true;
}

std::vector<double> updateDigestionParameter(std::vector<double> p) {

	// pick one to modify
	int ind = irand(0, p.size());
	int i = irand(0,1);
	
//	double change;
//	if (i == 0) {
//		change = -0.01;
//	}
//	else {
//		change = 0.01;
//	}
//	p[ind] += change;

	for(int i = 0; i < p.size(); i++) {
		if (p[i] != -1) {
			p[i] += drand(-0.01, 0.01); // random step in [-0.01, 0.01]
//			p[i] += irand(-1, 1) * 0.01; // randomly select from [-0,01, 0, 0.01]
		}
	}
	return p;
}

std::vector<double> get_per_site_digestion_rates(std::vector<double> &per_nuc_digestion_rates, std::string seq) {
	// per_nuc_digestion_rates should be [A,C,T,G]
	std::vector<double> per_site_digestion_rates;
	for (int j=0; j<seq.length(); j++) {
		if (seq[j] == 'A')
			per_site_digestion_rates.push_back(per_nuc_digestion_rates[0]);
		else if (seq[j] == 'C')
			per_site_digestion_rates.push_back(per_nuc_digestion_rates[1]);
		else if (seq[j] == 'U')
			per_site_digestion_rates.push_back(per_nuc_digestion_rates[2]);
		else if (seq[j] == 'G')
			per_site_digestion_rates.push_back(per_nuc_digestion_rates[3]);
		else {
			cerr << "Invalid character in seq " << seq << endl;
			exit(-1);
		}
	}
	return per_site_digestion_rates;
}

void writeToLog (ofstream &fp, std::string type, std::vector<double> &u, std::vector<double> &v, int result) {
	std::string result_str;
	if (type.compare("structure")==0) {
		switch (result) {
			case 1: result_str = "structure:accepted_better_LL"; break;
			case 2: result_str = "structure:accepted_invalid_structure"; break;
			case 3: result_str = "structure:rejected_invalid_structure"; break;
			case 4: result_str = "structure:accepted_inferior_LL"; break;
			case 5: result_str = "structure:rejected_inferior_LL"; break;
		}
	}
	else if (type.compare("digestion")==0) {
		switch (result) {
			case 1: result_str = "digestion:accepted_better_LL:"; break;
			case 2: result_str = "digestion:accepted_inferior_LL:"; break;
			case 3: result_str = "digestion:rejected_inferior_LL:"; break;
		}
		result_str += boost::lexical_cast<std::string>(u[0]) + "," + boost::lexical_cast<std::string>(u[1]) + "," + boost::lexical_cast<std::string>(u[2]) + "," + boost::lexical_cast<std::string>(u[3]) + ":" + boost::lexical_cast<std::string>(v[0]) + "," + boost::lexical_cast<std::string>(v[1]) + "," + boost::lexical_cast<std::string>(v[2]) + "," + boost::lexical_cast<std::string>(v[3]);
	}
	else {
		cout << "Invalid type in writeToLog: " << type << "\n";
		exit(-1);
	}
	fp << result_str << endl;
}


void writeToLog (ofstream &fp, std::string type, int result, double currentLL, double newLL, double q1, double q2, double mcConstant, std::vector<int> &moveCounts) {
	std::string result_str;
	std::string movecount_str = "0-" + boost::lexical_cast<std::string>(moveCounts[0]) + ",1-" + boost::lexical_cast<std::string>(moveCounts[1]) + ",2-" + boost::lexical_cast<std::string>(moveCounts[2]) + ",3-" + boost::lexical_cast<std::string>(moveCounts[3]);
	
	if (type.compare("structure")==0) {
		switch (result) {
			case 1: result_str = "structure:accepted_better_LL:" + boost::lexical_cast<std::string>(currentLL) + ":" + boost::lexical_cast<std::string>(newLL) + ":" + boost::lexical_cast<std::string>(q1) + ":" + boost::lexical_cast<std::string>(q2) + ":" + boost::lexical_cast<std::string>(mcConstant) + ":" + movecount_str; break;
			case 2: result_str = "structure:accepted_invalid_structure:" + boost::lexical_cast<std::string>(currentLL) + ":" + boost::lexical_cast<std::string>(newLL) + ":" + boost::lexical_cast<std::string>(q1) + ":" + boost::lexical_cast<std::string>(q2) + ":" + boost::lexical_cast<std::string>(mcConstant) + ":" + movecount_str; break;
			case 3: result_str = "structure:rejected_invalid_structure:" + boost::lexical_cast<std::string>(currentLL) + ":" + boost::lexical_cast<std::string>(newLL) + ":" + boost::lexical_cast<std::string>(q1) + ":" + boost::lexical_cast<std::string>(q2) + ":" + boost::lexical_cast<std::string>(mcConstant) + ":" + movecount_str; break;
			case 4: result_str = "structure:accepted_inferior_LL:" + boost::lexical_cast<std::string>(currentLL) + ":" + boost::lexical_cast<std::string>(newLL) + ":" + boost::lexical_cast<std::string>(q1) + ":" + boost::lexical_cast<std::string>(q2) + ":" + boost::lexical_cast<std::string>(mcConstant) + ":" + movecount_str; break;
			case 5: result_str = "structure:rejected_inferior_LL:" + boost::lexical_cast<std::string>(currentLL) + ":" + boost::lexical_cast<std::string>(newLL) + ":" + boost::lexical_cast<std::string>(q1) + ":" + boost::lexical_cast<std::string>(q2) + ":" + boost::lexical_cast<std::string>(mcConstant) + ":" + movecount_str; break;
		}
	}
	else if (type.compare("digestion")==0) {
		switch (result) {
			case 1: result_str = "digestion:accepted_better_LL"; break;
			case 2: result_str = "digestion:accepted_inferior_LL"; break;
			case 3: result_str = "digestion:rejected_inferior_LL"; break;
		}
	}
	else {
		cout << "Invalid type in writeToLog: " << type << "\n";
		exit(-1);
	}
	fp << result_str << endl;
}

std::string vec2string(std::vector<double> &values, std::string sep) {
  stringstream ss;
	for(int i = 0; i < values.size(); i++)
	{
		if(i != 0)
		  ss << sep;
		ss << values[i];
	}
	return ss.str();
}

std::string intvec2string(std::vector<int> &values, std::string sep) {
  stringstream ss;
	for(int i = 0; i < values.size(); i++)
	{
		if(i != 0)
		  ss << sep;
		ss << values[i];
	}
	return ss.str();
}

double calculateCloneNormalizationFactor(list<Read> &reads, int lcut, int rcut) {
	double count = 0;
	list<Read>::iterator it;
  for ( it=reads.begin() ; it != reads.end(); it++ ) {
  	count += (it->count * 1.0 / it->nocc);
  }
  return count;
}

// Calculates the average read length (rounded to nearest integer)
int estimateFragmentSize(list<Read> &reads) {
	double top, bottom;
	list<Read>::iterator it;
  for ( it=reads.begin() ; it != reads.end(); it++ ) {
  	top += (it->count * it->len);
  	bottom += (it->count);
  }
  return (int) (top/bottom);
}

// Calculates the value of the log-likelihood function for the given structure, digestion parameters, and read data
double compute_LL (std::string structure, list<Read> &reads, vector<double> &u, vector<double> &v, double r, int lcut, int rcut, int model) {
	// variable declaration
	double totalLL = 0;
	double ll = 0;
	double extend = 0;
	int start, stop;
	
	// pairing state (0 or 1) and smoothed digestion probabilities
	int structure_length = structure.length();	
	int pairing[structure_length];
	double digestion_probs[structure_length];
	double nodigestion_probs[structure_length];
//	digestion_probs[0] = 1;
//	digestion_probs[structure_length] = 1;
//	nodigestion_probs[0] = 0;
//	nodigestion_probs[structure_length] = 0;
	for (int i=0; i<structure_length; i++) {
		if ((structure[i] == '.') || (structure[i] == 'x'))
			pairing[i] = 0;
		else
			pairing[i] = 1;
	}
	// get digestion probability curve (no smoothing) - represents probability of digestion/modification 3' of the given position
	// NOTE: this means the first and last positions are meaningless (since they don't necessarily result from digestion/modification)
	for (int i=0; i<structure_length; i++) {
		int p;
		p = pairing[i];
		digestion_probs[i] = p*u[i] + (1-p)*v[i] + r; // (if paired, then u, else v) + random fragmentation
		nodigestion_probs[i] = 1-digestion_probs[i];
	}
	
	
//	// apply simple smoothing function to get digestion probability curve
//	for (int i=1; i<structure_length; i++) {
//		int p1, p2, p3;
//		p1 = pairing[i-1];
//		p2 = pairing[i];
//		p3 = pairing[i+1];
//		if ((model==1) || (model==2)) {
//			// dsRNA-seq/ssRNA-seq model (cut at unpaired positions 5' and 3')
//			digestion_probs[i] = p2*u + (1-p2)*v + r; // no averaging
//			nodigestion_probs[i] = 1-digestion_probs[i]; // does not cut wp 0.9 at paired, 0.8 at unpaired
//		}
//		else if (model==3) {
//			if ((p1==0) && (p2==0) && (p3==0))
//				digestion_probs[i] = v+r;
//			else
//				digestion_probs[i] = u+r;
//			nodigestion_probs[i] = 1-digestion_probs[i];
//		}
//	}
//	for (int i=0; i<structure_length+1; i++) {
//		printf("i=%d %lf\n", i,digestion_probs[i]);
//	}
/*	exit(0);*/

	// Fill in read fragment probability matrix ll_lookup, where each [i,j] element is the probability of the given fragment [i,j]
	typedef boost::multi_array<double, 2> Double2d;
	Double2d ll_lookup(boost::extents[structure_length][structure_length]);
	double sum = 0;
	for (int i=0; i<structure_length; i++) {
		for (int j=0; j<structure_length; j++) {
			ll_lookup[i][j] = 0;
		}
	}
	for (int i=1; i<structure_length; i++) {
		for (int j=i; j<structure_length; j++) {
			if ((j-i+1 < lcut) || (j-i+1 > rcut)) {
				ll_lookup[i][j] = 0;
			}
			else {
				// passes size restriction
				if (model==1) {
					ll = 0;
					// contribution from digestion 5' of the first position of the read
					ll += log(digestion_probs[i-1]);
					// contribution from non-digestion along the length of the read
					for (int k=i; k<j; k++) {
						ll += log(nodigestion_probs[k]);
					}
					// contribution from digestion 3' of the last position of the read
					ll += log(digestion_probs[j]);
					ll_lookup[i][j] = ll;
				}
				else if (model==2) {
					ll = 0;
					// contribution from modification 5' of the first position of the read
					ll += log(digestion_probs[i-1]);
					// contribution from non-modification along the length of the read
					for (int k=i; k<j; k++) {
						ll += log(nodigestion_probs[k]);
					}
					ll_lookup[i][j] = ll;
				}
			} // end passes size restriction
		} // end j
	} // end i
	
	// convert to likelihood
	for (int i=0; i<structure_length; i++) {
		for (int j=i; j<structure_length; j++) {
			if (ll_lookup[i][j] != 0) {
				// keep 0 values as 0
				ll_lookup[i][j] = exp(ll_lookup[i][j]);
				sum += ll_lookup[i][j];
//				cerr << "converting value L[" << i << "][" << j << "]=" << ll_lookup[i][j] << " to likelihood" << endl;
			}
		}
	}
	// normalize to sum==1 and convert back to log-likelihood
	for (int i=0; i<structure_length; i++) {
		for (int j=i; j<structure_length; j++) {
			if (ll_lookup[i][j] != 0) {
				ll_lookup[i][j] = log(ll_lookup[i][j] / sum);
//				cerr << "normalize value ll_lookup[" << i << "][" << j << "]=" << ll_lookup[i][j] << " and back to LL" << endl;
			}
		}
	}
	
//	for (int i=0; i<structure_length; i++) {
//		for (int j=i; j<structure_length; j++) {
//			cout <<  "ll_lookup[" << i << "][" << j << "] = " << ll_lookup[i][j] << endl;
//		}
//	}
//	exit(-1);
	
	// Go through the reads and keep a running total of the log-likelihood
	// should get start, end, count, nocc
	list<Read>::iterator it;
  for ( it=reads.begin() ; it != reads.end(); it++ ) {
  	// for now, ignore partially mapping reads (use start and stop, instead of checking against the actual read length)
  	// fix this later by extending the matrix ll_lookup to flanking regions beyond the locus
//  	cerr << "looking up read at " << it->start << " " << it->end << endl;
		extend = it->len - (it->end - it->start + 1);
//		cout << "checking len = " << it->len << " against [start, end] " << it->start << " " << it->end << endl;
		if (extend > 0) {
			// this is an ambiguous read - maps in the middle of the locus, but not for the full length of the read
			// we can't tell if the missing part of the full-length read belongs at the front or back of the read, so discard it
			// we only need to check that the stored read length does not match its start/end coordinates because the reads should have already been processed to be the correct length otherwise
//			cout << "skipping read with extend = " << extend << endl;
			continue;
		}
		ll = ll_lookup[it->start][it->end];
//		cerr << "looking up read at " << it->start << " " << it->end << " with count " << it->count << " and LL=" << ll << endl;
		totalLL += (double) it->count * ll;
//  	totalLL += ((double) it->count / (double) it->nocc)*ll; // weight for read count divided by number of genomic alignments
  }
	return(totalLL);
}

// Returns a list of valid move strings
vector<Move> getValidMoves(std::string structure, std::string seq) {
	vector<Move> moves;
	int m = 5; // some arbitrary constant for how far a single paren can move
	vector<int> dots, stack;
	dots.reserve(structure.length());
	stack.reserve(structure.length());
	vector<Pair> open_closes;
	map<int, int> pairs;
	for (int i=0; i<structure.length(); i++) {
		if (structure[i] == '.')
			dots.push_back(i);
		else if (structure[i] == '(')
			stack.push_back(i);
		else if (structure[i] == ')') {
			pairs[stack.back()] = i; // add this pairing to the list to check
			stack.pop_back();
		}
		else if (structure[i] == 'x') {
			// fixed unpaired flanking position - immutable
		}
		else {
			cerr << "Invalid character found in structure " << structure << endl;
			exit(1);
		}
	}
	// get the correctly paired opens and closes
	open_closes.reserve(pairs.size()); // reserves storage for the vector
	map<int, int>::iterator endIter = pairs.end();
	for (map<int, int>::iterator it=pairs.begin(); it!=endIter; it++) {
		Pair p(it->first, it->second);
		open_closes.push_back(p);
	}
	// get valid moves - remove a pair of matching parens
	for (int i=0; i<open_closes.size(); i++) {
		int openpos = open_closes[i].p1;
		int closepos = open_closes[i].p2;
		Move m(2,openpos,closepos);
		moves.push_back(m);
	}
	// get valid moves - insert a pair of matching parens
	if (dots.size() >= 2) {
		for (int i=0; i<dots.size(); i++) {
			for (int j=i+1; j<dots.size(); j++) {
				// for each possible (i,j) combination of dot positions
				std::string newstruct = structure;
				int openpos = dots[i];
				int closepos = dots[j];
				newstruct[openpos] = '(';
				newstruct[closepos] = ')';
				if (checkStructure(newstruct,seq)) {
					Move m(3,openpos,closepos);
					moves.push_back(m);
				}
				// flip order of open and close
				newstruct[openpos] = ')';
				newstruct[closepos] = '(';
				if (checkStructure(newstruct,seq)) {
					Move m(3,closepos,openpos);
					moves.push_back(m);
				}
			}
		}
	}
	return moves;
}

// Returns number of valid moves possible from current structure and sequence
int numValidMoves(std::string structure, std::string seq) {
	vector<Move> moves = getValidMoves(structure, seq);
	return moves.size();
}

// Returns a count of the move types from the current structure and sequence
std::vector<int> getMoveCounts(std::string structure, std::string seq) {
	std::vector<int> counts (4,0);
	vector<Move> moves = getValidMoves(structure, seq);
	for (std::vector<Move>::iterator it = moves.begin(); it != moves.end(); it++) {
		counts[(*it).type]++;
	}
	return counts;
}

// Return a modified structure from a random move in the structure landscape
std::string randomMove(std::string structure, std::string seq, bool weightMFE) {
	// DEFINED MOVE LIST:
	//	1: Shift a single open paren by x positions to either side, where x is random from [1,m]
	//	2: Shift a single close paren by x positions to either side, where x is random from [1,m]
	//	3: Remove a pair of matching parens
	//	4: Add a pair of matching parens
	std::string newstruct = structure;
	vector<Move> moves = getValidMoves(structure, seq);
	int ind;
	int msize = moves.size();
	if (msize == 0) {
		cerr << "Failed to get any valid moves - this is extremely unlikely and therefore probably a bug" << endl;
		exit(-1);
	}

	if (weightMFE) {
		std::vector<double> weights = getMoveWeights(moves, seq);
		ind = weighted_irand(weights);
	}
	else {
	 ind = irand(0,msize-1);
	}
	Move m = moves[ind];
	if (m.type == 0) {
		newstruct[m.p1] = '.';
		newstruct[m.p2] = '(';
	}
	else if (m.type == 1) {
		newstruct[m.p1] = '.';
		newstruct[m.p2] = ')';
	}
	else if (m.type == 2) {
		newstruct[m.p1] = '.';
		newstruct[m.p2] = '.';
	}
	else if (m.type == 3) {
		newstruct[m.p1] = '(';
		newstruct[m.p2] = ')';
	}
	else {
		cerr << "invalid move type " << m.type << endl;
		exit(-1);
	}
	return newstruct;
}

std::string makeMove(std::string structure, std::string seq, Move m) {
// DEFINED MOVE LIST:
	//	1: Shift a single open paren by x positions to either side, where x is random from [1,m]
	//	2: Shift a single close paren by x positions to either side, where x is random from [1,m]
	//	3: Remove a pair of matching parens
	//	4: Add a pair of matching parens
	std::string newstruct = structure;
	if (m.type == 0) {
		newstruct[m.p1] = '.';
		newstruct[m.p2] = '(';
	}
	else if (m.type == 1) {
		newstruct[m.p1] = '.';
		newstruct[m.p2] = ')';
	}
	else if (m.type == 2) {
		newstruct[m.p1] = '.';
		newstruct[m.p2] = '.';
	}
	else if (m.type == 3) {
		newstruct[m.p1] = '(';
		newstruct[m.p2] = ')';
	}
	else {
		cerr << "invalid move type " << m.type << endl;
		exit(-1);
	}
	return newstruct;
}


std::vector<double> getMoveWeights(std::vector<Move> moves, std::string seq) {
	std::vector<double> weights;
	
	for (std::vector<Move>::iterator it=moves.begin(); it!=moves.end(); it++) {
		double weight;
		if ((*it).type == 2) {
			// deleting a pair
			if ( ((seq[(*it).p1] == 'G') && (seq[(*it).p2] == 'C')) || ((seq[(*it).p1] == 'C') && (seq[(*it).p2] == 'G')) ) {
				weight = 1.0;
			}
			else if ( ((seq[(*it).p1] == 'A') && (seq[(*it).p2] == 'U')) || ((seq[(*it).p1] == 'U') && (seq[(*it).p2] == 'A')) ) {
				weight = 2.0;
			}
			else if ( ((seq[(*it).p1] == 'G') && (seq[(*it).p2] == 'U')) || ((seq[(*it).p1] == 'U') && (seq[(*it).p2] == 'G')) ) {
				weight = 3.0;
			}
			else {
				weight = 4.0;
			}
		}
		else if ((*it).type == 3) {
			// adding a pair
			if ( ((seq[(*it).p1] == 'G') && (seq[(*it).p2] == 'C')) || ((seq[(*it).p1] == 'C') && (seq[(*it).p2] == 'G')) ) {
				weight = 4.0;
			}
			else if ( ((seq[(*it).p1] == 'A') && (seq[(*it).p2] == 'U')) || ((seq[(*it).p1] == 'U') && (seq[(*it).p2] == 'A')) ) {
				weight = 3.0;
			}
			else if ( ((seq[(*it).p1] == 'G') && (seq[(*it).p2] == 'U')) || ((seq[(*it).p1] == 'U') && (seq[(*it).p2] == 'G')) ) {
				weight = 2.0;
			}
			else {
				weight = 1.0;
			}
		}
		else {
			cerr << "invalid move type " << (*it).type << endl;
			exit(-1);
		}
		weights.push_back(weight);
	}
	return weights;
}

// function to generate better seeds for srand
// http://www.concentric.net/~Ttwang/tech/inthash.htm
unsigned long mix(unsigned long a, unsigned long b, unsigned long c)
{
    a=a-b;  a=a-c;  a=a^(c >> 13);
    b=b-c;  b=b-a;  b=b^(a << 8);
    c=c-a;  c=c-b;  c=c^(b >> 13);
    a=a-b;  a=a-c;  a=a^(c >> 12);
    b=b-c;  b=b-a;  b=b^(a << 16);
    c=c-a;  c=c-b;  c=c^(b >> 5);
    a=a-b;  a=a-c;  a=a^(c >> 3);
    b=b-c;  b=b-a;  b=b^(a << 10);
    c=c-a;  c=c-b;  c=c^(b >> 15);
    return c;
}

// generate a random double in [min, max]
double drand(double min, double max) {
	double p = (double) rand() / (double) RAND_MAX; // [0,1]
	double range = max - min;
	return p*range + min; // [min,max]
}

// generate a random int in [min, max]
int irand(int min, int max) {
	if (min == max)
		return min;
	else
		return rand() % (max-min+1) + min;
}

// generate a random int in [0, weights.size-1] weighted by &weights
int weighted_irand(vector<double> &weights) {
	double cumsum = 0;
	double r = rand()/double(RAND_MAX); // old school (C++ 11 has real RNG functions)
	for(int i=0; i < weights.size(); ++i) {
		cumsum += weights[i];
		if (r < cumsum)
		  return i;
	}
	return weights.size()-1;
}

void doMCMC(std::string structure, std::string seq, list<Read> &reads, vector<double> &u, vector<double> &v, double r, int lcut, int rcut, int model, double mcConstant, int num_iters, bool doStructure, bool doDigestion, bool doFragmentation, bool weightMFE, bool likelihoodOnly, bool doLog, std::string log_fn, std::string out_fn) {
	double currentLL, ll;
	double newLL;
	bool hadValidStructure;
	bool validStructure;
	int iternum;
	double p, p2;
	
	int structure_length = structure.length();	
	
	int struct_mcmc_counter = 0;
	double q1, q2;
  // open log file for writing
  ofstream logfile, outfile;
  if (doLog) {
  	logfile.open(log_fn.c_str());
	}
	outfile.open(out_fn.c_str());
	
	bool validInitial = checkDigestionParameters(u,v,r);
	if (!validInitial) {
		cerr << "Found invalid starting digestion parameters" << endl;
		exit(-1);
	}
	
	
// calculate initial likelihood
	std::vector<double> per_site_u = get_per_site_digestion_rates(u, seq);
	std::vector<double> per_site_v = get_per_site_digestion_rates(v, seq);
	ll = compute_LL(structure,reads,per_site_u,per_site_v,r,lcut,rcut,model);
	
	// 5/7/13 - calculate effective clone count to scale likelihood ratio into proper Metropolis-Hastings Uniform[0,1] range
	// NOTE: this is necessary because the likelihood table is normalized so that the P(any given read fragment) sums to 1; and we use the product of all these density values as the currentLL/newLL totals
	double norm;
	norm = calculateCloneNormalizationFactor(reads, lcut, rcut);

	if (norm > 0)
		ll /= norm;
	
  currentLL = ll;
  hadValidStructure = checkStructure(structure,seq);
  // print out initial result
	std::string u_string = vec2string(u, ",");
	std::string v_string = vec2string(v, ",");
  outfile << structure << "\t" << currentLL << "\t" << u_string << "\t" << v_string << "\t" << r << "\t0" << endl;
	
	if (likelihoodOnly) {
		outfile.close();
		exit(-1);
	}
	
	// MCMC iterations
	for (iternum=1; iternum<=num_iters; iternum++) {
		// Step A - random move along structure landscape
		int stepresult = -1;
		if (doStructure) {
//			cerr << "getting randomMove" << endl;
			std::string newstruct = randomMove(structure, seq, weightMFE);
//			cerr << "computing numValidMoves" << endl;
			q1 = (double) numValidMoves(structure, seq);
//			cerr << "computing numValidMovies2" << endl;
			q2 = (double) numValidMoves(newstruct, seq);
//			q1 = 5;
//			q2 = 5;
//			cout << "got q1 q2 = " << q1 << " " << q2 << endl;
//			cerr << "getting moveCounts" << endl;
			std::vector<int> moveCounts = getMoveCounts(structure, seq);
//			cout << "got moveCounts " << endl;
			std::vector<double> per_site_u = get_per_site_digestion_rates(u, seq);
			std::vector<double> per_site_v = get_per_site_digestion_rates(v, seq);
//			cerr << "computing likelihood" << endl;
			ll = compute_LL(newstruct,reads,per_site_u,per_site_v,r,lcut,rcut,model);
			
			// added 5/7/13 - normalize ll by number of clones
			if (norm > 0)
				ll /= norm;
			newLL = ll;
			double storedCurrent = currentLL;
//			cerr << "checkStructure" << endl;
//			cout << "computed newLL = " << newLL << endl;
			validStructure = checkStructure(newstruct, seq);
//			cout << "comparing currentLL = " << currentLL << " with newLL = " << newLL << " diff = " << newLL-currentLL << endl;
//			cout << "validStructre = " << validStructure  << " q1 = " << q1 << " q2 = " << q2 << " mcConstant = " << mcConstant << endl;
			if ((newLL > currentLL) && (validStructure)) {
				// found a better valid structure
				currentLL = newLL;
				structure = newstruct;
				hadValidStructure = true;
				stepresult = 1;
			}
			else if (newLL > currentLL) {
				// a better structure, but with invalid base pairing - accept with very low probability
				p2 = drand(0,1);	// changed 6/24/13 - dont care about invalid structure base pairing
				if ((p2 > 0.99) || (!hadValidStructure)) {
					currentLL = newLL;
					structure = newstruct;
					hadValidStructure = false;
					stepresult = 2;
				}
				else {
					stepresult = 3;
				}
			}
			else if ((newLL-currentLL+log(q1)-log(q2)) >= log(drand(0,1)*mcConstant)) {
				// accept under Metropolis-Hastings criterion
				currentLL = newLL;
				structure = newstruct;
				hadValidStructure = validStructure;
				stepresult = 4;
			}
			else {
				// do nothing
				stepresult = 5;
			}

			struct_mcmc_counter++;
			if (doLog) {
				writeToLog(logfile, "structure", stepresult, storedCurrent, newLL, q1, q2, mcConstant, moveCounts);
			}
//			cerr << "finished step A" << endl;
		}
		// end step A

		// step B - random move along digestion rate landscape
		stepresult = -1;
		if (doDigestion) {
			bool goodstep = false;
			vector<double> newu, newv;
			while (!goodstep) {
				newu = updateDigestionParameter(u);
				newv = updateDigestionParameter(v);
				goodstep = checkDigestionParameters(newu,newv,r);
//				cerr << "stuck in parameter random walk loop " << newu << " " << newv << " " << endl;
			}
			std::vector<double> per_site_newu = get_per_site_digestion_rates(newu, seq);
			std::vector<double> per_site_newv = get_per_site_digestion_rates(newv, seq);
			ll = compute_LL(structure,reads,per_site_newu,per_site_newv,r,lcut,rcut,model);
			// added 5/7/13 - normalize ll by number of clones
			if (norm > 0)
				ll /= norm;
			newLL = ll;
			
			if (newLL > currentLL) {
				// better digestion parameters
//				cerr << "accepting new digestion params " << u << " " << v << endl;
				currentLL = newLL;
				u = newu;
				v = newv;
				stepresult = 1;
			}
			else if ((newLL-currentLL) >= log(drand(0,1)*mcConstant)) {
				// accept inferior likelihood under Metropolis-Hastings
				currentLL = newLL;
				u = newu;
				v = newv;
				stepresult = 2;
			}
			else {
				stepresult = 3;
			}
			struct_mcmc_counter = 0;
			if (doLog) {
				writeToLog(logfile, "digestion", u, v, stepresult);
			}
		}
		// end step B

		// print out results for this iteration
		std::string u_string = vec2string(u, ",");
		std::string v_string = vec2string(v, ",");
		outfile << structure << "\t" << currentLL << "\t" << u_string << "\t" << v_string << "\t" << r << "\t" << iternum << endl;
	} // end MCMC iterations
	if ((iternum % 1000) == 0)
		cerr << "completed " << iternum << "/" << num_iters << " iterations" << endl;
	
	logfile.close();
	outfile.close();
}

// END FUNCTION DEFINITIONS
// ====================================================================================================
