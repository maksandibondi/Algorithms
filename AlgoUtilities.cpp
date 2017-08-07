#include "AlgoUtilities.h"

namespace AlgoUtilities {

	Population::Population(){}

	Population::Population(int& size, bool firstIteration) {

		/*
		std::vector<Individual> individuals(0);
		individuals.resize(size);
		if (firstIteration) {

			std::vector<Individual>::iterator i;
			for (i = individuals.begin(); i != individuals.end(); i++) {
				Individual* newIndividual = new Individual();
				*i = *newIndividual;
				delete newIndividual;
			}
			
		}
		*/

		for (int i = 0; i < size; i++) {
			Individual* newIndividual = new Individual();
			individuals.push_back(*newIndividual);
			delete newIndividual;
		}


	}

	Individual& Population::getIndividual(int& index) {
		return individuals[index];
	}

	void Population::setIndividual(int& index, Individual indiv) {
		individuals[index] = indiv;
	}

	Individual Population::getFittest() {
		Individual fittest = individuals[0];
		// Loop through individuals to find fittest
		int sz = individuals.size();
		for (int i = 0; i < sz; i++) {

			if (fittest.getFitness() <= getIndividual(i).getFitness()) {
				fittest = getIndividual(i);
			}
		}
		return fittest;
	}

	Individual Population::getFittestForBS(MarketData md, DealData dd) {
		Individual fittest = individuals[0];
		// Loop through individuals to find fittest
		int sz = individuals.size();
		for (int i = 0; i < sz; i++) {

			if (fittest.getFitnessForBSModel(md,dd) <= getIndividual(i).getFitnessForBSModel(md,dd)) {
				fittest = getIndividual(i);
			}
		}
		return fittest;
	}

	void Population::addAnIndividual(Individual indiv) {
		individuals.push_back(indiv);
	}

	int Population::size() {
		int sz = individuals.size();
		return sz;
	}



	// new individual creates an individual with certain combination of genes
	Individual::Individual() {

		genes.resize(precision, false);
				
		for (int i = 0; i < precision; i++) {
			setGene(i, (bool)std::round((double)rand()/ (double)RAND_MAX));
		}

	}
	
	bool Individual::getGene(int index) {
		bool gene = (bool)genes[index];
		return gene;
	}

	void Individual::setGene(int index, bool value) {
		genes[index] = value;
	}

	int Individual::getFitness() {
		int fit = 0;
		int sz = (this->genes).size();

		for (int i = 0; i < sz; i++) {
			
			if (this->getGene(i) == (solution[i])) {
				fit++;
			}
		}
		this->fitness = fit;
		return fit;
	}

	int Individual::getFitnessForBSModel(MarketData md, DealData dd) {
		int fit = 0;
		int sz = (this->genes).size();
		md.sigma = convert(genes);

			for (int i = 0; i < sz; i++) {

				if (BSSqrDiffBitwise(md, dd)[i] == (solution[i])) { // compare i-th bit of sum of differences with its bit of solution
					fit++;
				}
			}
			this->fitness = fit;
			return fit;
	}

	void Individual::setSolution(boost::dynamic_bitset<> sol) {
		solution = sol;
		precision = sol.size();
	}

	int Individual::getPrecision() {
		return precision;
	}



	Population GeneticAlgo::evolvePopulation(Population pop) {
		int size = pop.size();
		Population* newPopulation = new Population(size, false);
		// Keep the best individual
		if (elitism) {
			newPopulation->addAnIndividual(pop.getFittest());
		}
		
		// add the new individuals to newPopulation individuals vector
		for (int i = 0; i < size-elitism; i++) {
			Individual indiv1 = tournamentSelection(pop);
			Individual indiv2 = tournamentSelection(pop);
			Individual newIndiv = crossover(indiv1, indiv2);
			newPopulation->setIndividual(i, newIndiv);
		}

		// Mutate population
		for (int i = 0; i < size; i++) {
			mutate(newPopulation->getIndividual(i));
		}

		return *newPopulation;
	}
	
	// returns a value of new daughter individual
	Individual GeneticAlgo::crossover(Individual indiv1, Individual indiv2) {
		Individual* newSol = new Individual();
		int precision = Individual::getPrecision();
		// Loop through genes
		for (int i = 0; i < precision; i++) {
			// Crossover
			if ((double)rand()/(double)RAND_MAX <= uniformRate) {
				newSol->setGene(i, indiv1.getGene(i));
			}
			else {
				newSol->setGene(i, indiv2.getGene(i));
			}
		}
		return *newSol;
	}

	// Mutate an individual (changes the genes inside an individual
	void GeneticAlgo::mutate(Individual indiv) {
		int precision = Individual::getPrecision();
		// Loop through genes
		for (int i = 0; i < precision; i++) {
			if ((double)rand()/(double)RAND_MAX <= mutationRate) {
				// Create random gene
				bool gene = (bool)round((double)rand()/(double)RAND_MAX);
				indiv.setGene(i, gene);
			}
		}
	}

	// takes some random individuals from population and chooses the fittest
	Individual GeneticAlgo::tournamentSelection(Population pop) {
		Individual* fittest = new Individual();
		int sz = pop.size();
		// Create a tournament population
		Population* tournament = new Population(tournamentSize, false);
		// For each place in the tournament get a random individual
		for (int i = 0; i < tournamentSize; i++) {
			int randomId = (int)(((double)rand()/(double)RAND_MAX) * (sz-1));
			tournament->setIndividual(i, pop.getIndividual(randomId));
		}
		// Get the fittest
		*fittest = tournament->getFittest();
		return *fittest;
	}

	// initialize static algo input
	void  GeneticAlgo::initializeAlgoInput(double uniformRate, double mutationRate, int tournamentSize, bool elitism) {
		GeneticAlgo::uniformRate = uniformRate;
		GeneticAlgo::mutationRate = mutationRate;
		GeneticAlgo::tournamentSize = tournamentSize;
		GeneticAlgo::elitism = elitism;
	}

	// static variables default values
	boost::dynamic_bitset<> Individual::solution = *(new boost::dynamic_bitset<>(64));
	int Individual::precision = 64;
	double GeneticAlgo::uniformRate = 0.5;
	double GeneticAlgo::mutationRate = 0.05;
	int GeneticAlgo::tournamentSize = 5;
	bool GeneticAlgo::elitism = 0;








	DealData::DealData() {
		K = 50;
		T = { 0.25, 0.5, 0.75, 1 };
	}

	MarketData::MarketData() {
		S = 50;
		r = 0;
		prices = { 1.2, 1.5, 1.7, 1.72 };
	}







	
	boost::dynamic_bitset<> BSSqrDiffBitwise(MarketData md, DealData dd) {
		double S = md.S;
		double r = md.r;
		double K = dd.K;
		double sigma = md.sigma;
		std::vector<double> T = dd.T;
		double sumOfTheSqrDifference = 0;

		int sz = T.size();

		for (int i = 0; i < sz; i++) {

			double d1 = (1 / (sigma * sqrt(T[i])))*(log(S / K) + (r + pow(sigma, 2) / 1)*T[i]);
			//cout << "d1 = " << d1.getValue() << endl;

			double d2 = d1 - sigma * sqrt(T[i]);
			//cout << "d2 = " << d2.getValue() << endl;

			double price = (NormalCDFCody(d1)*S) - (NormalCDFCody(d2)*K*exp(-r*T[i]));

			sumOfTheSqrDifference = sumOfTheSqrDifference + pow((price - md.prices[i]),2);
		}

		boost::dynamic_bitset<> x(64,sumOfTheSqrDifference);

		return x;

		
		
	}

	double NormalCDFCody(double u) {
		double y = abs(u);
		if (y > 35.0) {
			if (u > 0)
				return 1;
			else
				return 0;
		}
		if (y <= 0.662912607) {
			//  evaluate erf() for |u| <= sqrt(2)*0.46875
			double a0 = 1.161110663653770e-2;
			double a1 = 3.951404679838207e-1;
			double a2 = 2.846603853776254e+1;
			double a3 = 1.887426188426510e+2;
			double a4 = 3.209377589138469e+3;

			double b0 = 1.767766952966369e-1;
			double b1 = 8.344316438579620;
			double b2 = 1.725514762600375e+2;
			double b3 = 1.813893686502485e+3;
			double b4 = 8.044716608901563e+3;

			double z = y * y;
			y = u * ((((a0 * z + a1) * z + a2) * z + a3) * z + a4);
			y /= ((((b0 * z + b1) * z + b2) * z + b3) * z + b4);
			return 0.5 + y;
		}
		double zinterm = 0.5 * exp(-y * y / 2);
		if (y <= 4.0) {
			double c0 = 2.15311535474403846e-8;
			double c1 = 5.64188496988670089e-1;
			double c2 = 8.88314979438837594;
			double c3 = 6.61191906371416295e+1;
			double c4 = 2.98635138197400131e+2;
			double c5 = 8.81952221241769090e+2;
			double c6 = 1.71204761263407058e+3;
			double c7 = 2.05107837782607147e+3;
			double c8 = 1.23033935479799725e+3;

			double d0 = 1.0;
			double d1 = 1.57449261107098347e+1;
			double d2 = 1.17693950891312499e+2;
			double d3 = 5.37181101862009858e+2;
			double d4 = 1.62138957456669019e+3;
			double d5 = 3.29079923573345963e+3;
			double d6 = 4.36261909014324716e+3;
			double d7 = 3.43936767414372164e+3;
			double d8 = 1.23033935480374942e+3;

			// evaluate erfc() for sqrt(2)*0.46875 <= |u| <= sqrt(2)*4.0
			y = y / 1.4142135623730950488;
			double num = ((((((((c0 * y + c1) * y + c2) * y + c3) * y + c4) * y + c5) * y + c6) * y + c7) * y + c8);
			double den = ((((((((d0 * y + d1) * y + d2) * y + d3) * y + d4) * y + d5) * y + d6) * y + d7) * y + d8);

			y = num / den;
			y = zinterm * y;
		}
		else {
			double p0 = 1.63153871373020978e-2;
			double p1 = 3.05326634961232344e-1;
			double p2 = 3.60344899949804439e-1;
			double p3 = 1.25781726111229246e-1;
			double p4 = 1.60837851487422766e-2;
			double p5 = 6.58749161529837803e-4;

			double q0 = 1.00000000000000000;
			double q1 = 2.56852019228982242;
			double q2 = 1.87295284992346047;
			double q3 = 5.27905102951428412e-1;
			double q4 = 6.05183413124413191e-2;
			double q5 = 2.33520497626869185e-3;
			// evaluate erfc() for |u| > sqrt(2)*4.0
			double z = zinterm * 1.41421356237309504880 / y;
			y = 2 / (y * y);
			y = y * (((((p0 * y + p1) * y + p2) * y + p3) * y + p4) * y + p5) / (((((q0 * y + q1) * y + q2) * y + q3) * y + q4) * y + q5);
			y = z*(0.564189583547756287 - y);
		}

		if (u < 0)
			return y;
		return 1 - y;
	}

	double convert(boost::dynamic_bitset<> const& bs) {
		std::bitset<64> myBit(0);
		for (boost::dynamic_bitset<>::size_type i = 0; i < bs.size(); i++) {
			myBit[i] = bs[i];
		}		
		static_assert(sizeof(uint64_t) == sizeof(double), "Cannot use this!");

		uint64_t const u = myBit.to_ullong();
		double d;

		// Aliases to `char*` are explicitly allowed in the Standard (and only them)
		char const* cu = reinterpret_cast<char const*>(&u);
		char* cd = reinterpret_cast<char*>(&d);

		// Copy the bitwise representation from u to d
		memcpy(cd, cu, sizeof(u));

		return d;
	}

}