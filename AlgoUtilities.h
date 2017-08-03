#pragma once

#ifdef ALGOUTILITIES_EXPORTS  
#define ALGOUTILITIES_API __declspec(dllexport)   
#else  
#define ALGOUTILITIES_API __declspec(dllimport)   
#endif  

#include <vector>
#include <algorithm>
#include <cmath>
#include <boost/dynamic_bitset.hpp>


namespace AlgoUtilities {
	
	class Individual {
		static boost::dynamic_bitset<> solution;
		boost::dynamic_bitset<> genes;
		int fitness = 0;
	public:
		Individual();
		Individual(int& precision);
		bool getGene(int index);
		void setGene(int index, bool value);
		int getFitness();
		static void setSolution(boost::dynamic_bitset<> sol);

	};

	class Population {
		int precision;
		std::vector<Individual>* individuals;
	public:
		Population();
		Population(int& size, int& precision);
		Individual& getIndividual(int& index);
		Individual& getFittest();
	};

	class GeneticAlgo {

	};

}