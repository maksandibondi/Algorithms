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

	int Individual::getFitnessForModel(const char* model, MarketData md, DealData dd) {
		int fit = 0;
		int sz = (this->genes).size();

		if (model == "BS") {
			for (int i = 0; i < sz; i++) {

				if (BSPricebitwise(md, dd) == (solution[i])) {
					fit++;
				}
			}
			this->fitness = fit;
			return fit;
		}
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

	std::vector<boost::dynamic_bitset<>> DealData::getMaturityInBits() {
		// here we have to translate each element of T array to bits and create a vector 
	}

	MarketData::MarketData() {
		S = 50;
		r = 0;
		prices = { 1.2, 2.2, 3.2, 4.2};
	}

	std::vector<boost::dynamic_bitset<>> MarketData::getPricesInBits() {
		// here we have to translate each element of prices array to bits and create a vector 
	}

}