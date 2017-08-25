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

	static class Matrix {
		size_t d1, d2;
		std::vector<double> data; //linearazation vector for matrix
	public:
		Matrix(size_t d1, size_t d2);

		// rechanrging of operator = is not necessary as next function will return a reference to the member that we wanna set

		double & operator()(size_t i, size_t j);

	};

	static class DealData3D{
	public:
		int discretization_num_T;
		int discretization_num_K;
		DealData3D();
		std::vector<double> T;
		std::vector<double> K;
	};

	static class MarketData3D{
	public:
		MarketData3D();
		double S;
		double r;
		std::vector<std::vector<double>> sigma;
		std::vector<std::vector<double>>  prices;
	};

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
		boost::dynamic_bitset<> genes; // target variable binary expression
		double target;
		int fitness = 0;
		double fitnessDouble = 0;
	public:
		Individual();
		bool getGene(int index);
		double getTarget();
		void setGene(int index, bool value, bool &stateMin, bool &stateMax);
		int getFitness();
		//int getFitnessForBSModel(MarketData md, DealData dd);
		double getFitnessForBSModel(MarketData md, DealData dd);
		bool acceptGene(bool value, bool& stateMin, bool& stateMax, int index, bool isSign);
		static void setSolution(boost::dynamic_bitset<> sol);
		static int getPrecision();

	};

	class Population {
		std::vector<Individual> individuals;
	public:
		Population();
		Population(int& size);
		Individual& getIndividual(int& index);
		void setIndividual(int& index, Individual indiv);
		Individual getFittest();
		Individual getFittestForBS(MarketData md, DealData dd);
		//Individual getFittestForBS2(MarketData md, DealData dd);
		void addAnIndividual(Individual indiv);
		int size();
	};

	class Individual3D {
		static boost::dynamic_bitset<> solution;
		static int precision;
		boost::dynamic_bitset<> genes; // target variable binary expression
		double target;
		int fitness = 0;
		double fitnessDouble = 0;
	public:
		Individual3D();
		bool getGene(int index);
		double getTarget();
		void setGene(int index, bool value, bool &stateMin, bool &stateMax);
		int getFitness();
		//int getFitnessForBSModel(MarketData md, DealData dd);
		double getFitnessForBSModel(MarketData md, DealData dd);
		bool acceptGene(bool value, bool& stateMin, bool& stateMax, int index, bool isSign);
		static void setSolution(boost::dynamic_bitset<> sol);
		static int getPrecision();
	};

	class Population3D {
		std::vector<Matrix> individuals;
	public:
		Population3D();
		Population3D(int& size);
		Individual& getIndividual(int& index);
		void setIndividual(int& index, Individual indiv);
		Individual getFittest();
		Individual getFittestForBS(MarketData3D md, DealData3D dd);
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
		static boost::dynamic_bitset<> minset; // binary constraints for sigma
		static boost::dynamic_bitset<> maxset; // binary constraints for sigma
		static double minval; // double constraint for sigma
		static double maxval; // double constraint for sigma

		static Population evolvePopulation(Population pop);
		static Individual tournamentSelection(Population pop);
		static Individual crossover(Individual indiv1, Individual indiv2);
		static void mutate(Individual indiv);
		static void initializeAlgoInput(double uniformRate, double mutationRate, int tournamentSize, bool elitism);
		static void setSystemBinaryConstraints(boost::dynamic_bitset<> valmin, boost::dynamic_bitset<> valmax);
		static void setSystemDoubleConstraints(double valmin, double valmax);
		
		static boost::dynamic_bitset<> convertDoubleTo64Bit(double value);
		static boost::dynamic_bitset<> convertIntToBit(int value);
		static boost::dynamic_bitset<> convertIntTo11Bit(int value);
		static boost::dynamic_bitset<> convertFractionToBit(double value);
		static boost::dynamic_bitset<> getExponentMantissaByNormalization(boost::dynamic_bitset<>bitIntegralPart, boost::dynamic_bitset<> bitFractionalPart, boost::dynamic_bitset<> &mantissa);
		static int convertBitToInt(boost::dynamic_bitset<> value);
		static double convertBitToFraction(boost::dynamic_bitset<> value);
		static double convertBitToDouble(boost::dynamic_bitset<> value);
		static std::vector<bool> convertBitsetToVector(boost::dynamic_bitset<> array);
	};

	static boost::dynamic_bitset<> BSSqrDiffBitwise(MarketData md, DealData dd);

	std::vector<std::vector<double>> FDMLocalVolpricer(MarketData3D md, DealData3D dd);

	static double NormalCDFCody(double u);

	static double convert(boost::dynamic_bitset<> const& bs);

}