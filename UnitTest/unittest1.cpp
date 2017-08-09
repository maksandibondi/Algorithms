#include "stdafx.h"
#include "CppUnitTest.h"
#include "AlgoUtilities.h"
#include <iostream>
#include <vector>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace AlgoUtilities;

namespace UnitTest {		

	TEST_CLASS(UnitTest) {
public:

	TEST_METHOD(TestMethod1) {
		int precision = 64;
		boost::dynamic_bitset<> sol(precision);
		Individual::setSolution(sol);
		int myPopulationSize = 100;

		Population* myPop = new Population(myPopulationSize, true);
		int generationCount = 0;
		std::vector<int> fit;

		while (myPop->getFittest().getFitness() < precision) {
			//fit.push_back(myPop->getFittest().getFitness());
			generationCount++;
			//std::cout << "Generation: " << generationCount << " Fittest: " << myPop->getFittest().getFitness() << std::endl;
			*myPop = GeneticAlgo::evolvePopulation(*myPop);
		}

		Assert::AreEqual(double(myPop->getFittest().getFitness()), double(precision));

	}

	// bs test
	/*TEST_METHOD(TestMethod2) {
		int precision = 64; // precision in bits
		boost::dynamic_bitset<> sol(precision);
		Individual::setSolution(sol);
		int myPopulationSize = 100;
		boost::dynamic_bitset<> valmin = GeneticAlgo::convertDoubleTo64Bit(0); // add the value of valmin, valmax from web or from algo of conversion
		boost::dynamic_bitset<> valmax = GeneticAlgo::convertDoubleTo64Bit(2);

		GeneticAlgo::setSystemConstraints(valmin, valmax);
		Population* myPop = new Population(myPopulationSize, true);
		int generationCount = 0;
		std::vector<int> fit;

		DealData dd = DealData::DealData();
		MarketData md = MarketData::MarketData();

		while (myPop->getFittestForBS(md,dd).getFitnessForBSModel(md,dd) < precision) {
			//fit.push_back(myPop->getFittest().getFitness());
			generationCount++;
			//std::cout << "Generation: " << generationCount << " Fittest: " << myPop->getFittest().getFitness() << std::endl;
			*myPop = GeneticAlgo::evolvePopulation(*myPop);
		}

		Assert::AreEqual(double(myPop->getFittestForBS(md, dd).getFitnessForBSModel(md, dd)), double(precision));
		//Assert::AreEqual(md.sigma, double(precision));
	}*/

	TEST_METHOD(TestMethod2) {

		GeneticAlgo::initializeAlgoInput(0.5, 0.05, 5, 0);
		int precision = 64; // precision in bits
		boost::dynamic_bitset<> sol(precision);
		Individual::setSolution(sol);
		int myPopulationSize = 100;
		boost::dynamic_bitset<> valmin = GeneticAlgo::convertDoubleTo64Bit(0); // add the value of valmin, valmax from web or from algo of conversion
		boost::dynamic_bitset<> valmax = GeneticAlgo::convertDoubleTo64Bit(2.5);
		GeneticAlgo::setSystemConstraints(valmin, valmax);

	}

	TEST_METHOD(TestMethod3) {
		GeneticAlgo::initializeAlgoInput(0.5, 0.05, 5, 0);
		int precision = 64; // precision in bits
		boost::dynamic_bitset<> sol(precision);
		Individual::setSolution(sol);
		int myPopulationSize = 100;
		boost::dynamic_bitset<> valmin = GeneticAlgo::convertDoubleTo64Bit(0); // add the value of valmin, valmax from web or from algo of conversion
		boost::dynamic_bitset<> valmax = GeneticAlgo::convertDoubleTo64Bit(2.5);
		GeneticAlgo::setSystemConstraints(valmin, valmax);
		Population* myPop = new Population(myPopulationSize, true);
		int generationCount = 0;
		std::vector<int> fit;

		DealData dd = DealData::DealData();
		MarketData md = MarketData::MarketData();

		while (myPop->getFittestForBS(md, dd).getFitnessForBSModel(md, dd) < precision) {
			//fit.push_back(myPop->getFittest().getFitness());
			generationCount++;
			//std::cout << "Generation: " << generationCount << " Fittest: " << myPop->getFittest().getFitness() << std::endl;
			*myPop = GeneticAlgo::evolvePopulation(*myPop);
		}

		Assert::AreEqual(double(myPop->getFittestForBS(md, dd).getFitnessForBSModel(md, dd)), double(precision));
		//Assert::AreEqual(md.sigma, double(precision));
	}


	};

}