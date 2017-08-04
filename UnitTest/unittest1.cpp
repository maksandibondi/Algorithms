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
			std::cout <<"Solution found!" << std::endl;
			std::cout << "Generation: " + generationCount << std::endl;
			std::cout << "Genes: " << myPop->getFittest().getFitness();

			Assert::AreEqual(double(myPop->getFittest().getFitness()), double(precision));

		}

	};

}