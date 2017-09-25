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
			myPop = GeneticAlgo::evolvePopulation(*myPop);
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

	TEST_METHOD(TestMatrixOperations) {

		Matrix<double> *prices = new Matrix<double>({ { 7.71436943, 18.12343211, 19.77587219, 21.34355787 },{ 5.740060514, 15.64069661, 17.42982267, 19.0910727 },{ 4.111298367, 13.3226896, 15.2318196, 16.97210398 },{ 2.831019104, 11.19474524, 13.19631397, 14.99605698 } });
		double x = (*prices)(0, 0);
		Matrix<double> *test = new Matrix<double>({ { 6.71436943, 17.12343211, 18.77587219, 20.34355787 },{ 5.740060514, 15.64069661, 17.42982267, 19.0910727 },{ 4.111298367, 13.3226896, 15.2318196, 16.97210398 },{ 2.831019104, 11.19474524, 13.19631397, 14.99605698 } });

		Matrix<double> res_substr = *prices - *test;
		Matrix<double> res_sqr = (*prices) ^ 2;
		double sum = res_sqr.sumOfElements();

		Matrix<double> a(4, 4, char('id'));

		system("pause");
	}

	TEST_METHOD(Test3DEntities) {
		GeneticAlgo::initializeAlgoInput(0.5, 0.05, 5, 0);
		int precision = 64; // precision in bits
		boost::dynamic_bitset<> sol(precision);
		Individual::setSolution(sol);
		int myPopulationSize = 20;
		GeneticAlgo::setSystemDoubleConstraints(0, 1);

		DealData3D* dd = new DealData3D(DealData3D::DealData3D());
		MarketData3D* md = new MarketData3D(MarketData3D::MarketData3D());
		md->prices = BSPriceMatrixCreator(md, dd); // CREATING THE MATRIX OF MARKET PRICES FORM MARKET DATA
		Population3D::setDimensionsOf3DIndividuals(md->prices->size(0), md->prices->size(1));
		Individual3D::setDimensionsOf3DIndividuals(md->prices->size(0), md->prices->size(1));

		Population3D* myPop = new Population3D(myPopulationSize);
		
		int i = 0;
		Matrix<double>* x = myPop->getIndividual(i)->getTargetMatrix();
		bool gene_test = myPop->getIndividual(i)->getGene(0, 0, 63);

		//boost::dynamic_bitset<> x2 = BSSqrDiffBitwise3D(md, dd);
		//double res = GeneticAlgo::convertBitToDouble(x2);


 		//system("pause");
	}

	TEST_METHOD(TestSumOfSqrDif) {
		DealData3D* dd = new DealData3D(DealData3D::DealData3D());
		MarketData3D* md = new MarketData3D(MarketData3D::MarketData3D());
		md->sigma = new Matrix<double>(dd->discretization_num_T, dd->discretization_num_K, 0.2);
		md->prices = BSPriceMatrixCreator(md, dd); // CREATING THE MATRIX OF MARKET PRICES FORM MARKET DATA

		double sumOfSqrDifRel = BSSqrDiff3D(md, dd);

		system("pause");


	}

	TEST_METHOD(TestFDMLocalVolPricer) {
		DealData3D* dd = new DealData3D(DealData3D::DealData3D());
		MarketData3D* md = new MarketData3D(MarketData3D::MarketData3D());
		md->sigma = new Matrix<double>(dd->discretization_num_T, dd->discretization_num_K, 0.2);
		md->prices = BSPriceMatrixCreator(md, dd); // CREATING THE MATRIX OF MARKET PRICES FORM MARKET DATA

		Matrix<double> u = FDMLocalVolpricer(md, dd);
		Matrix<double> dif = (u - *(md->prices))/(*(md->prices));
		Matrix<double> difabs = (u - *(md->prices));
		Matrix<double> res = dif.abs();
		double max = res.max();

		double sumdif = dif.sumOfElements()/u.size();

		system("pause");

	}

	TEST_METHOD(TestFDMLocalVolpricerThetaScheme) {
		DealData3D* dd = new DealData3D(DealData3D::DealData3D());
		MarketData3D* md = new MarketData3D(MarketData3D::MarketData3D());
		md->sigma = new Matrix<double>(dd->discretization_num_T, dd->discretization_num_K, 0.2);
		md->prices = BSPriceMatrixCreator(md, dd); // CREATING THE MATRIX OF MARKET PRICES FORM MARKET DATA

		Matrix<double> u = FDMLocalVolpricerThetaScheme(md, dd, 0); // Implicit scheme pricer
		Matrix<double> u2 = FDMLocalVolpricer(md, dd); // Explicit scheme pricer 
		Matrix<double> dif = (u - *(md->prices)) / (*(md->prices));
		Matrix<double> dif2 = (u2 - *(md->prices)) / (*(md->prices));

		Matrix<double> difabs = (u - *(md->prices));
		Matrix<double> res = dif.abs();
		double max = res.max();

		double sumdif = dif.sumOfElements() / u.size();

		system("pause");

	}

	TEST_METHOD(TestPopRowsColsFromMatrix) {
		Matrix<double> a(3, 4, 0);
		a(0, 0) = 1;
		a(1, 0) = 2;
		a(1, 3) = 1;
		a(2, 2) = 3;

		a.pop({ 0,2 }, { 0,3 });
		//system("pause");
	}

	TEST_METHOD(TestMatrixARamaCont) {
		DealData3D* dd = new DealData3D(DealData3D::DealData3D());
		MarketData3D* md = new MarketData3D(MarketData3D::MarketData3D());
		md->sigma = new Matrix<double>(dd->discretization_num_T, dd->discretization_num_K, 0.2);
		md->prices = BSPriceMatrixCreator(md, dd); // CREATING THE MATRIX OF MARKET PRICES FORM MARKET DATA

		

		system("pause");

	}

	// algos to test 
	TEST_METHOD(TestGeneticAlgoBS) {
		GeneticAlgo::initializeAlgoInput(0.5, 0.06, 5, 0);
		int precision = 64; // precision in bits
		boost::dynamic_bitset<> sol(precision);
		Individual::setSolution(sol);
		int myPopulationSize = 20;
		GeneticAlgo::setSystemDoubleConstraints(0, 1);

		DealData dd = DealData::DealData();
		MarketData md = MarketData::MarketData();
		//Population3D::setDimensionsOf3DIndividuals(md.prices->size(0), md.prices->size(1));
		//Individual3D::setDimensionsOf3DIndividuals(md.prices->size(0), md.prices->size(1));

		Population* myPop = new Population(myPopulationSize);
		int generationCount = 0;
		std::vector<int> fit;


		//double f = myPop->getFittestForBS(md, dd).getFitnessForBSModel(md, dd);

		while (myPop->getFittestForBS(md, dd).getFitnessForBSModel(md, dd) < -0.01) {
			generationCount = generationCount + 1;
			myPop = GeneticAlgo::evolvePopulation(*myPop);
		}
		double gc = generationCount;
		double res = (myPop->getFittestForBS(md, dd)).getTarget();
		Assert::AreEqual(double(myPop->getFittestForBS(md, dd).getFitnessForBSModel(md, dd)), double(-0.01));

		//Assert::AreEqual(md.sigma, double(precision));
	}

	// algos to test 
	TEST_METHOD(TestGeneticAlgoBSLocalVol) {
		GeneticAlgo::initializeAlgoInput(0.5, 0.01, 5, 0);
		int precision = 64; // precision in bits
		boost::dynamic_bitset<> sol(precision);
		Individual::setSolution(sol);
		int myPopulationSize = 30;
		GeneticAlgo::setSystemDoubleConstraints(0, 0.7);

		DealData3D* dd = new DealData3D(DealData3D::DealData3D());
		MarketData3D* md = new MarketData3D(MarketData3D::MarketData3D());
		md->prices = BSPriceMatrixCreator(md, dd); // CREATING THE MATRIX OF MARKET PRICES FORM MARKET DATA
		size_t sz1 = md->prices->size(0);
		size_t sz2 = md->prices->size(1);
		Population3D::setDimensionsOf3DIndividuals(md->prices->size(0), md->prices->size(1));

		Population3D* myPop = new Population3D(myPopulationSize);
		
		int generationCount = 0;

		std::vector<double>* f = new std::vector<double>();
		f->push_back(myPop->getFittestForBS(md, dd)->getFitnessForBSModel(md, dd));
		while ((*f)[generationCount] < -0.01) {
			generationCount = generationCount + 1;
			myPop = GeneticAlgo::evolvePopulation(myPop, md, dd);
			f->push_back(myPop->getFittestForBS(md, dd)->getFitnessForBSModel(md, dd));
		}

		double gc = generationCount;
		Matrix<double>* res = (myPop->getFittestForBS(md, dd))->getTargetMatrix();
		//Assert::AreEqual(md.sigma, double(precision));
		system("pause");
	}
	
	

	};

}