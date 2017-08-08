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
#include <iostream>
#include <random>
#include <bitset>


namespace AlgoUtilities {


	static class DealData {
	public:
		std::vector<double> T;
		double K;
		DealData();
		//std::vector<boost::dynamic_bitset<>> getMaturityInBits();
	};

	static class MarketData {
	public:
		double S;
		double r;
		double sigma;
		std::vector<double>  prices;
		MarketData();
		//std::vector<boost::dynamic_bitset<>> getPricesInBits();
	};

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
		int getFitnessForBSModel(MarketData md, DealData dd);
		bool acceptGene(bool value, bool stateMin, bool stateMax, int index, bool isSign);
		static void setSolution(boost::dynamic_bitset<> sol);
		static int getPrecision();

	};

	class Population {
		std::vector<Individual> individuals;
	public:
		Population();
		Population(int& size, bool firstIteration);
		Individual& getIndividual(int& index);
		void setIndividual(int& index, Individual indiv);
		Individual getFittest();
		Individual getFittestForBS(MarketData md, DealData dd);
		void addAnIndividual(Individual indiv);
		int size();
	};

	class GeneticAlgo {

		static double uniformRate;
		static double mutationRate;
		static int tournamentSize;
		static bool elitism;

	public:

		//GeneticAlgo();
		//GeneticAlgo(double uniformRate, double mutationRate, int tournamentSize, bool elitism);
		static boost::dynamic_bitset<> minset;
		static boost::dynamic_bitset<> maxset;
		static Population evolvePopulation(Population pop);
		static Individual tournamentSelection(Population pop);
		static Individual crossover(Individual indiv1, Individual indiv2);
		static void mutate(Individual indiv);
		static void initializeAlgoInput(double uniformRate, double mutationRate, int tournamentSize, bool elitism);
		static void setSystemConstraints(boost::dynamic_bitset<> valmin, boost::dynamic_bitset<> valmax);
	};

	static boost::dynamic_bitset<> BSSqrDiffBitwise(MarketData md, DealData dd);

	static double NormalCDFCody(double u);

	static double convert(boost::dynamic_bitset<> const& bs);

}