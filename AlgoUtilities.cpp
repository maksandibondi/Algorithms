#include "AlgoUtilities.h"

namespace AlgoUtilities {

	Population::Population(){}

	Population::Population(int& size, int& precision) {
		this->individuals = new std::vector<Individual>(size);
		std::vector<Individual>::iterator i;
		for (i = (*individuals).begin(); i != (*individuals).end(); i++) {
			Individual* newIndividual = new Individual(precision);
			*i = *newIndividual;
			delete newIndividual;
		}
	}

	Individual& Population::getIndividual(int& index) {
		return (*individuals)[index];
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

	Individual::Individual() {
	}

	Individual::Individual(int& precision) {

				
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
	}

	boost::dynamic_bitset<> Individual::solution = *(new boost::dynamic_bitset<>(64));
	
}