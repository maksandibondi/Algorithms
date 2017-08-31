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
			int sz = data.size();
			for (int i = 0; i < sz; i++) {
				data[i] = value;
			}
		}

		Matrix(std::vector<std::vector<T>> vec) {
			this->d1 = vec.size();
			this->d2 = vec[0].size();
			this->data.resize(d1*d2);
			int sz = data.size();
			for (int i = 0; i < d1; i++) {
				for (int j = 0; j < d2; j++) {
					data[i*d1 + j] = vec[i][j];
				}
			}
		}

		// rechanrging of operator = is not necessary as next function will return a reference to the member that we wanna set

		T & operator()(size_t i, size_t j) {
			return data[i*d1 + j];
		}

		Matrix<T> operator+(Matrix<T>& a) {
			size_t sz1 = this->d1;
			size_t sz2 = this->d2;
			Matrix<T> res = Matrix(d1, d2);

			for (int i = 0; i < d1; i++) {
				for (int j = 0; j < d2; j++) {
					res(i, j) = data[i*d1 + j] + a(i, j);
				}
			}
			return res;
		}

		Matrix<T> operator-(Matrix<T>& a) {
			size_t sz1 = this->d1;
			size_t sz2 = this->d2;
			Matrix<T> res = Matrix(d1, d2);

			for (int i = 0; i < d1; i++) {
				for (int j = 0; j < d2; j++) {
					res(i, j) = data[i*d1 + j] - a(i, j);
				}
			}
			return res;
		}

		Matrix<T> operator^(int power) {
			size_t sz1 = this->d1;
			size_t sz2 = this->d2;
			Matrix<T> res = Matrix(d1, d2);

			for (int i = 0; i < d1; i++) {
				for (int j = 0; j < d2; j++) {
					res(i, j) = pow(data[i*d1 + j], 2);
				}
			}
			return res;
		}

		// scalar product
		/*Matrix<T> operator()(Matrix<T>& a) {
			size_t sz1 = this->d1;
			size_t sz2 = this->d2;
			Matrix<T> res = Matrix(d1, d2);

			for (int i = 0; i < d1; i++) {
				for (int j = 0; j < d2; j++) {
					res(i, j) = data[i*d1 + j] - a(i, j);
				}
			}
			return res;
		}*/

		double sumOfElements() {
			double sum = 0;
			int sz = d1*d2;
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
		void setIndividual(int& index, Individual& indiv);
		Individual getFittest();
		Individual getFittestForBS(MarketData md, DealData dd);
		//Individual getFittestForBS2(MarketData md, DealData dd);
		void addAnIndividual(Individual indiv);
		int size();
	};

	class Individual3D {
		static boost::dynamic_bitset<> solution;
		static int precision;
		Matrix<Individual> *indivs;
		static size_t dim1;
		static size_t dim2;
		int fitness = 0;
		double fitnessDouble = 0;
	public:
		Individual3D();
		bool getGene(size_t idx1, size_t idx2, int genIndex);
		double getTarget(size_t idx1, size_t idx2);
		Matrix<Individual>* getIndividualsMatrix();
		Matrix<double>* getTargetMatrix();
		//void setGene(size_t d1, size_t d2, int genIndex, bool value, bool &stateMin, bool &stateMax);
		int getFitness();
		//int getFitnessForBSModel(MarketData md, DealData dd);
		double getFitnessForBSModel(MarketData3D md, DealData3D dd);
		static void setSolution(boost::dynamic_bitset<> sol);
		static int getPrecision();
		static void setDimensionsOf3DIndividuals(size_t d1, size_t d2);

	};

	class Population3D {
		std::vector<Individual3D> individuals;
	public:
		static size_t dim1; // dimension of Individual3D
		static size_t dim2; // dimension of Individual3D

		Population3D();
		Population3D(int& size);
		Individual3D& getIndividual(int& index);
		void setIndividual(int& index, Individual3D indiv);
		Individual3D getFittest();
		Individual3D getFittestForBS(MarketData3D md, DealData3D dd);
		void addAnIndividual(Individual3D indiv);
		int size();
		static void setDimensionsOf3DIndividuals(size_t d1, size_t d2);
	};

	class GeneticAlgo {

		static double uniformRate;
		static double mutationRate;
		static int tournamentSize;
		static bool elitism;

	public:

		static boost::dynamic_bitset<> minset; // binary constraints for sigma
		static boost::dynamic_bitset<> maxset; // binary constraints for sigma
		static double minval; // double constraint for sigma
		static double maxval; // double constraint for sigma

		static Population* evolvePopulation(Population& pop);
		static Individual* tournamentSelection(Population& pop);
		static Individual* crossover(Individual& indiv1, Individual& indiv2);
		static void mutate(Individual& indiv);

		
		static Population3D* evolvePopulation(Population3D& pop);
		static Individual3D* tournamentSelection(Population3D& pop);
		static Individual3D* crossover(Individual3D& indiv1, Individual3D& indiv2);
		static void mutate(Individual3D& indiv);
		

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