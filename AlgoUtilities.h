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
	
	template <typename T>
	class Matrix {
		size_t d1, d2;
		std::vector<T> data; //linearazation vector for matrix
	public:
		Matrix(size_t d1, size_t d2) {
			this->d1 = d1;
			this->d2 = d2;
			this->data.resize(d1*d2);
		}

		Matrix(size_t d1, size_t d2, T value) {
			this->d1 = d1;
			this->d2 = d2;
			this->data.resize(d1*d2);
			sz = data.size();
			for (int i = 0; i < sz; i++) {
				data[i] = value;
			}
		}

		// rechanrging of operator = is not necessary as next function will return a reference to the member that we wanna set

		T & operator()(size_t i, size_t j) {
			return data[i*d1 + j];
		}

		Matrix<T>* & operator-(Matrix<T>* a) {
			sz = data.size();
			for (i = 0; i < sz; i++) {
				data[i] = data[i] - a.data[i];
			}
			return data;
		}

		Matrix<T>* & operator^(int power) {
			sz = data.size();
			for (i = 0; i < sz; i++) {
				data[i] = pow(data[i],power);
			}
			return data;
		}

		double sumOfElements() {
			double sum = 0;
			sz = d1*d2;
			for (int i = 0; i < sz; i++) {
				sum += data[i];
			}
			return sum;
		}

		size_t size() {
			return d1*d2;
		}

		size_t size(int dim) {
			if (dim == 0) {
				return d1;
			}
			else if (dim == 1) {
				return d2;
			}
			
		}

		Matrix<T> vectorToMatrix(std::vector<std::vector(T)> vec) {
			
			sz1 = vec.size();
			sz2 = vec[0].size();
			Matrix<T> *matrix = new Matrix(sz1, sz2);
			for (int i = 0; i < sz1; i++) {
				for (int j = 0; j < sz2; j++) {
					(*matrix)(i, j) = vec[i][j];
				}
			}
		}

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
		Matrix<double> *sigma;
		Matrix<double> *prices;
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
		friend class Individual3D;
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
		//Matrix <boost::dynamic_bitset<>> genes; // target variable binary expression
		//Matrix<double> target;
		Matrix<Individual> *indivs;
		size_t d1 = 0;
		size_t d2 = 0;
		int fitness = 0;
		double fitnessDouble = 0;
	public:
		Individual3D();
		Individual3D(size_t d1, size_t d2);
		bool getGene(size_t idx1, size_t idx2, int genIndex);
		double getTarget(size_t idx1, size_t idx2);
		Matrix<double>* Individual3D::getTargetMatrix();
		//void setGene(size_t d1, size_t d2, int genIndex, bool value, bool &stateMin, bool &stateMax);
		int getFitness();
		//int getFitnessForBSModel(MarketData md, DealData dd);
		double getFitnessForBSModel(MarketData3D md, DealData3D dd);
		static void setSolution(boost::dynamic_bitset<> sol);
		static int getPrecision();
	};

	class Population3D {
		std::vector<Individual3D> individuals;
	public:
		Population3D();
		Population3D(int& size, size_t d1, size_t d2);
		Individual3D& getIndividual(int& index);
		void setIndividual(int& index, Individual3D indiv);
		Individual3D getFittest();
		Individual3D getFittestForBS(MarketData3D md, DealData3D dd);
		void addAnIndividual(Individual3D indiv);
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

	static boost::dynamic_bitset<> BSSqrDiffBitwise3D(MarketData3D md, DealData3D dd);

	Matrix<double>* FDMLocalVolpricer(MarketData3D md, DealData3D dd);

	static double NormalCDFCody(double u);

	static double convert(boost::dynamic_bitset<> const& bs);

}