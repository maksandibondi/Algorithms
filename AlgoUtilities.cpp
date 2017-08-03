#include "AlgoUtilities.h"

namespace AlgoUtilities {

	Population::Population(){}

	Population::Population(int& size, bool firstIteration) {

		this->individuals = new std::vector<Individual>(size);

		if (firstIteration) {

			std::vector<Individual>::iterator i;
			for (i = (*individuals).begin(); i != (*individuals).end(); i++) {
				Individual* newIndividual = new Individual();
				*i = *newIndividual;
				delete newIndividual;
			}

		}
	}

	Individual& Population::getIndividual(int& index) {
		return (*individuals)[index];
	}

	void Population::setIndividual(int& index, Individual indiv) {
		(*individuals)[index] = indiv;
	}

	Individual& Population::getFittest() {
		Individual fittest = (*individuals)[0];
		// Loop through individuals to find fittest
		int sz = individuals->size();
		for (int i = 1; i < sz; i++) {

			if (fittest.getFitness() <= getIndividual(i).getFitness()) {
				fittest = getIndividual(i);
			}
		}
		return fittest;
	}

	void Population::addAnIndividual(Individual indiv) {
		individuals->push_back(indiv);
	}

	int& Population::size() {
		int sz = individuals->size();
		return sz;
	}




	Individual::Individual() {

				
		for (int i = 0; i < precision; i++) {
			setGene(i, (bool)std::round(rand()));
		}

	}
	
	bool Individual::getGene(int index) {
		return (bool)genes[index];
	}

	void Individual::setGene(int index, bool value) {
		genes[index] = value;
	}

	int Individual::getFitness() {
		int fit = this->fitness;
		int sz = (this->genes).size();

		for (int i = 0; i < sz; i++) {
			
			if (this->getGene(i) == (solution[i])) {
				fit++;
			}
		}
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
		Population* newPopulation = new Population(size, 0);
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
			if (rand() <= uniformRate) {
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
			if (rand() <= mutationRate) {
				// Create random gene
				bool gene = (bool)round(rand());
				indiv.setGene(i, gene);
			}
		}
	}

	// takes some random individuals from population and chooses the fittest
	Individual GeneticAlgo::tournamentSelection(Population pop) {
		int sz = pop.size();
		// Create a tournament population
		Population* tournament = new Population(tournamentSize, 0);
		// For each place in the tournament get a random individual
		for (int i = 0; i < tournamentSize; i++) {
			int randomId = (int)(rand() * sz);
			tournament->setIndividual(i, pop.getIndividual(randomId));
		}
		// Get the fittest
		Individual fittest = tournament->getFittest();
		return fittest;
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


}