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
		static int precision;
		boost::dynamic_bitset<> genes;
		int fitness = 0;
	public:
		Individual();
		bool getGene(int index);
		void setGene(int index, bool value);
		int getFitness();
		static void setSolution(boost::dynamic_bitset<> sol);
		static int getPrecision();

	};

	class Population {
		std::vector<Individual>* individuals = nullptr;
	public:
		Population();
		Population(int& size, bool firstIteration);
		Individual& getIndividual(int& index);
		void setIndividual(int& index, Individual indiv);
		Individual& getFittest();
		void addAnIndividual(Individual indiv);
		int& size();
	};

	class GeneticAlgo {
		 static double uniformRate;
		 static double mutationRate;
		 static int tournamentSize;
		 static bool elitism;

	public:

		GeneticAlgo();
		GeneticAlgo(double uniformRate, double mutationRate, int tournamentSize, bool elitism);

		static Population evolvePopulation(Population pop); 
		static Individual tournamentSelection(Population pop);
		static Individual crossover(Individual indiv1, Individual indiv2);
		static void mutate(Individual indiv);
		static void initializeAlgoInput(double uniformRate, double mutationRate, int tournamentSize, bool elitism);

	};

}