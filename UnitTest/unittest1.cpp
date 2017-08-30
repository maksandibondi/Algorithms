#include "stdafx.h"
#include "CppUnitTest.h"
#include "AlgoUtilities.h"
#include <iostream>
#include <vector>
#include <algorithm>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace AlgoUtilities;

namespace UnitTest {		

	TEST_CLASS(UnitTest) {
public:


	TEST_METHOD(TestGeneticAlgoPOC) {
		GeneticAlgo::initializeAlgoInput(0.5, 0.00, 5, 0);
		int precision = 64;
		boost::dynamic_bitset<> sol(precision);
		Individual::setSolution(sol);
		int myPopulationSize = 20;
		GeneticAlgo::setSystemDoubleConstraints(0, 1);

		Population* myPop = new Population(myPopulationSize);
		int generationCount = 0;
		std::vector<int> fit;

		while (myPop->getFittest().getFitness() < precision) {
			generationCount++;
			*myPop = GeneticAlgo::evolvePopulation(*myPop);
		}

		Assert::AreEqual(double(myPop->getFittest().getFitness()), double(precision));

	}

	// Conversion tests
	TEST_METHOD(TestGeneticAlgoSettings) {

		GeneticAlgo::initializeAlgoInput(0.5, 0.05, 5, 0);
		int precision = 64; // precision in bits
		boost::dynamic_bitset<> sol(precision);
		Individual::setSolution(sol);
		int myPopulationSize = 100;
		boost::dynamic_bitset<> valmin = GeneticAlgo::convertDoubleTo64Bit(0); // add the value of valmin, valmax from web or from algo of conversion
		boost::dynamic_bitset<> valmax = GeneticAlgo::convertDoubleTo64Bit(2.5);
		GeneticAlgo::setSystemBinaryConstraints(valmin, valmax);

	}

	TEST_METHOD(TestConversionIntToBit) {

		int i = 10;
		boost::dynamic_bitset<> a = GeneticAlgo::convertIntToBit(i);
		int sz = a.size();
		std::vector<bool> test(sz, 0);
		for (int k = 0; k < sz; k++) {
			test[k] = a[k];
		}
		std::vector<bool> test1 = test;
	}

	TEST_METHOD(TestConversionIntTo11Bit) {
		int i = 1023;
		boost::dynamic_bitset<> a = GeneticAlgo::convertIntTo11Bit(i);
		int sz = a.size();
		std::vector<bool> test(sz, 0);
		for (int k = 0; k < sz; k++) {
			test[k] = a[k];
		}
		std::vector<bool> test1 = test;
	}
	
	TEST_METHOD(TestConversionFractionToBit) {
		double i = 0.15;
		boost::dynamic_bitset<> a = GeneticAlgo::convertFractionToBit(i);
		int sz = a.size();
		std::vector<bool> test(sz, 0);
		for (int k = 0; k < sz; k++) {
			test[k] = a[k];
		}
		std::vector<bool> test1 = test;
	}

	TEST_METHOD(TestConversionBitToInt) {
		std::string str = std::string("01011001110");
		std::reverse(str.begin(), str.end());
		boost::dynamic_bitset<> testarray(str); // string writes values in the opposite order

		int sz = testarray.size();
		std::vector<bool> test(sz, 0);
		for (int k = 0; k < sz; k++) {
			test[k] = testarray[k];
		}
		int value = GeneticAlgo::convertBitToInt(testarray);

	}

	TEST_METHOD(TestConversionBitToFraction) {
		std::string str = std::string ("01101");
		std::reverse(str.begin(), str.end());
		boost::dynamic_bitset<> testarray(str); // string writes values in the opposite order

		int sz = testarray.size();
		std::vector<bool> test(sz, 0);
		for (int k = 0; k < sz; k++) {
			test[k] = testarray[k];
		}

		double value = GeneticAlgo::convertBitToFraction(testarray);
	}

	TEST_METHOD(TestConversionBitToDouble) {
		std::string str = std::string("0011111101010100100000000010100100000000010100100000000010100100");
		std::reverse(str.begin(), str.end());
		boost::dynamic_bitset<> testarray(str); // string writes values in the opposite order

		double result = GeneticAlgo::convertBitToDouble(testarray);
	}
	
	TEST_METHOD(TestConversionDoubleTo64bit) {
		double sigma = 0.56358531449324012;
		std::string str = std::string("0011111101011000101111010110011000100111011111000100010111001100");
		std::reverse(str.begin(), str.end());
		boost::dynamic_bitset<> testarray(str); // string writes values in the opposite order

		boost::dynamic_bitset<> result = GeneticAlgo::convertDoubleTo64Bit(sigma);
		std::vector<bool> test1 = GeneticAlgo::convertBitsetToVector(result);
	}






	// algos to test 
	TEST_METHOD(TestGeneticAlgoBS) {
		GeneticAlgo::initializeAlgoInput(0.5, 0.06, 5, 0);
		int precision = 64; // precision in bits
		boost::dynamic_bitset<> sol(precision);
		Individual::setSolution(sol);
		int myPopulationSize = 20;
		GeneticAlgo::setSystemDoubleConstraints(0, 1);

		Population* myPop = new Population(myPopulationSize);
		int generationCount = 0;
		std::vector<int> fit;

		DealData dd = DealData::DealData();
		MarketData md = MarketData::MarketData();
		//double f = myPop->getFittestForBS(md, dd).getFitnessForBSModel(md, dd);

		while (myPop->getFittestForBS(md, dd).getFitnessForBSModel(md, dd) < -0.0001) {
			generationCount = generationCount + 1;
			*myPop = GeneticAlgo::evolvePopulation(*myPop);
		}
		double gc = generationCount;
		double res = (myPop->getFittestForBS(md, dd)).getTarget();
		Assert::AreEqual(double(myPop->getFittestForBS(md, dd).getFitnessForBSModel(md, dd)), double(-0.0001));

		//Assert::AreEqual(md.sigma, double(precision));
	}

	// algos to test 
	TEST_METHOD(TestGeneticAlgoBSLocalVol) {
		GeneticAlgo::initializeAlgoInput(0.5, 0.06, 5, 0);
		int precision = 64; // precision in bits
		boost::dynamic_bitset<> sol(precision);
		Individual::setSolution(sol);
		int myPopulationSize = 20;
		GeneticAlgo::setSystemDoubleConstraints(0, 1);
		DealData dd = DealData::DealData();
		MarketData md = MarketData::MarketData();

		Population3D* myPop = new Population3D(myPopulationSize, dd.T.size(), dd.K.size());
		int generationCount = 0;
		std::vector<int> fit;

		

		while (myPop->getFittestForBS(md, dd).getFitnessForBSModel(md, dd) < -0.0001) {
			generationCount = generationCount + 1;
			*myPop = GeneticAlgo::evolvePopulation(*myPop);
		}
		double gc = generationCount;
		double res = (myPop->getFittestForBS(md, dd)).getTarget();
		Assert::AreEqual(double(myPop->getFittestForBS(md, dd).getFitnessForBSModel(md, dd)), double(-0.0001));

		//Assert::AreEqual(md.sigma, double(precision));
	}

	};

}