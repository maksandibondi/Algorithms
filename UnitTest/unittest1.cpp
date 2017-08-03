#include "stdafx.h"
#include "CppUnitTest.h"
#include "AlgoUtilities.h"
#include <iostream>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace AlgoUtilities;

namespace UnitTest {		

	TEST_CLASS(UnitTest) {
	public:
		
		TEST_METHOD(TestMethod1) {
			int precision = 64;
			boost::dynamic_bitset<> sol(precision);
			Individual::setSolution(sol);
			int myPopulationSize = 10;
			
			Population* myPop = new Population(myPopulationSize, true);
			int generationCount = 0;
			while (myPop->getFittest().getFitness() < precision) {
				generationCount++;
				std::cout << "Generation: " << generationCount << " Fittest: " << myPop->getFittest().getFitness() << std::endl;
				*myPop = GeneticAlgo::evolvePopulation(*myPop);
			}
			std::cout <<"Solution found!" << std::endl;
			std::cout << "Generation: " + generationCount << std::endl;
			std::cout << "Genes: " << myPop->getFittest().getFitness();

			Assert::AreEqual(double(myPop->getFittest().getFitness()), double(0));

		}

	};

}