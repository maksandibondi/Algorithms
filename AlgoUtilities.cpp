#include "AlgoUtilities.h"

namespace AlgoUtilities {


	// Data

	DealData::DealData() {
		K = double(100);
		T = { 0.25, 0.5, 0.75, 1 , 1.25, 1.5 };
		//T = { double(30) / double(365) };
	}

	MarketData::MarketData() {
		S = double(110);
		r = 0.05;
		sigma = 0;
		prices = { 11.98832952, 14.07538404, 15.94706732, 17.66295374, 19.26355258, 20.77444988 };
		//prices = { 11.2 };
	}

	DealData3D::DealData3D() {
		discretization_num_T = 20;
		discretization_num_K = 55;
		
		K = *(new std::vector<double>(discretization_num_K, 0));
		for (int j = 1; j < discretization_num_K; j++) {
			K[j] = K[j - 1] + 4;
		}

		T = *(new std::vector<double>(discretization_num_T, 0));
		for (int i = 1; i < discretization_num_T; i++) {
			T[i] = T[i - 1] + 0.1;
		}
		//T = { double(30) / double(365) };
	}

	MarketData3D::MarketData3D() {
		S = double(110);
		r = 0.00;
		sigma = nullptr; // new Matrix<double>(4, 5, 0);
		prices = nullptr; //new Matrix<double>({ {8.133428294, 6.194460012, 4.571019117, 3.265104769, 2.256690667},{ 10.53118471, 8.685045271, 7.06117636, 5.65952825, 4.472227255},{ 12.55583683, 10.75226442, 9.12811083, 7.68294186, 6.412109843}, { 14.37034509, 12.5944325, 10.9718374, 9.502020105, 8.181743876} });
		//prices = { 11.2 };
	}


	// Population and Individuals

	Population::Population() {}

	Population::Population(int& size) {

		for (int i = 0; i < size; i++) {
			Individual* newIndividual = new Individual();
			individuals.push_back(*newIndividual);
			delete newIndividual;
		}


	}

	Individual& Population::getIndividual(int& index) {
		return individuals[index];
	}

	void Population::setIndividual(int& index, Individual& indiv) {
		individuals[index] = indiv;
	}

	Individual Population::getFittest() {
		Individual fittest = individuals[0];
		// Loop through individuals to find fittest
		int sz = individuals.size();
		for (int i = 0; i < sz; i++) {

			if (fittest.getFitness() <= getIndividual(i).getFitness()) {
				fittest = getIndividual(i);
			}
		}
		return fittest;
	}

	Individual Population::getFittestForBS(MarketData md, DealData dd) {
		Individual fittest = individuals[0];
		// Loop through individuals to find fittest
		int sz = individuals.size();
		for (int i = 0; i < sz; i++) {

			if (fittest.getFitnessForBSModel(md, dd) <= getIndividual(i).getFitnessForBSModel(md, dd)) {
				fittest = getIndividual(i);
			}
		}
		return fittest;
	}

	void Population::addAnIndividual(Individual indiv) {
		individuals.push_back(indiv);
	}

	int Population::size() {
		int sz = individuals.size();
		return sz;
	}



	Population3D::Population3D() {}

	Population3D::Population3D(int& size) {
		individuals = new std::vector<Individual3D*>();
		for (int i = 0; i < size; i++) {
			Individual3D* newIndividual = new Individual3D();
			individuals->push_back(newIndividual);
		}
	}

	Individual3D* Population3D::getIndividual(int& index) {
		return (*individuals)[index];
	}

	void Population3D::setIndividual(int& index, Individual3D* indiv) {
		(*individuals)[index] = indiv;
	}

	Individual3D* Population3D::getFittest() {
		Individual3D* fittest = (*individuals)[0];
		// Loop through individuals to find fittest
		int sz = individuals->size();
		for (int i = 0; i < sz; i++) {

			if (fittest->getFitness() <= getIndividual(i)->getFitness()) {
				fittest = getIndividual(i);
			}
		}
		return fittest;
	}

	Individual3D* Population3D::getFittestForBS(MarketData3D* md, DealData3D* dd) {
		Individual3D* fittest = (*individuals)[0];
		// Loop through individuals to find fittest
		int sz = individuals->size();
		for (int i = 0; i < sz; i++) {

			if (fittest->getFitnessForBSModel(md, dd) <= getIndividual(i)->getFitnessForBSModel(md, dd)) {
				fittest = getIndividual(i);
			}
		}
		return fittest;
	}

	void Population3D::addAnIndividual(Individual3D* indiv) {
		individuals->push_back(indiv);
	}

	int Population3D::size() {
		int sz = individuals->size();
		return sz;
	}

	void Population3D::setDimensionsOf3DIndividuals(size_t d1, size_t d2) {
		dim1 = d1;
		dim2 = d2;
		Individual3D::setDimensionsOf3DIndividuals(dim1, dim2);
	}


	// new individual creates an individual with certain combination of genes
	Individual::Individual() {

		bool stateMin;
		bool stateMax;

		genes.resize(precision, false);
		stateMin = false; // initial state for constraints
		stateMax = false; // initial state for constraints
		double sigma = GeneticAlgo::minval + ((double)rand() / (double)RAND_MAX)*(GeneticAlgo::maxval - GeneticAlgo::minval);
		this->target = sigma;
		boost::dynamic_bitset<> binarySigma = GeneticAlgo::convertDoubleTo64Bit(sigma);
		for (int i = 0; i < precision; i++) {
			setGene(i, binarySigma[i], stateMin, stateMax);
		}

		/*std::vector<double> checkbin(64, 0);
		std::vector<double> checkgenes(64, 0);
		for (int i = 0; i < 64; i++) {

			checkbin[i] = binarySigma[i]; // binary sigma is well calculated
			checkgenes[i] = genes[i]; // genes acception algo is not working well
		}*/

		/*double val = GeneticAlgo::convertBitToDouble(genes);
		double x = val;*/
	}

	bool Individual::getGene(int index) {
		bool gene = (bool)genes[index];
		return gene;
	}

	double Individual::getTarget() {
		return target;
	}

	void Individual::setGene(int index, bool value, bool &stateMin, bool &stateMax) {

		bool isSign = (index == 0);
		if (acceptGene(value, stateMin, stateMax, index, isSign)) {
			genes[index] = value;
		}
		else {
			genes[index] = 1 - value;
		}
	}

	bool Individual::acceptGene(bool value, bool& stateMin, bool& stateMax, int index, bool isSign) {
		bool acceptMin = false;
		bool acceptMax = false;
		bool min = GeneticAlgo::minset[index];
		bool max = GeneticAlgo::maxset[index];

		switch (isSign) {

		case (false):
			// for min
			if (stateMin == true) {
				acceptMin = true; // works for negative
			}
			else {
				if (min == 1 && value == 1) {
					acceptMin = true;
				}
				else if (min == 0) {
					acceptMin = true;
					if (value == 1) {
						stateMin = true;
					}
				}
				else if (min == 1) {
					if (value == 1) {
						acceptMin = true;
					}
				}
			}
			// for max
			if (stateMax == true) {
				acceptMax = true; // works for negative
			}
			else {
				if (max == 1 && value == 1) {
					acceptMax = true;
				}
				else if (max == 1) {
					acceptMax = true;
					if (value == 0) {
						stateMax = true;
					}
				}
				else if (max == 0) {
					if (value == 0) {
						acceptMax = true;
					}
				}
			}

			break;

		case (true):

			if (max == 0) { // if max is positive
				acceptMax = true;
			}
			else if (max == 1 && value == 1) { // if max is negative and value is negative
				acceptMax = true;
			}

			if (min == 1) { // if min is negative
				acceptMin = true;
			}
			else if (min == 0 && value == 0) { // if min is positive and value is positive
				acceptMin = true;
			}

			break;
		}

		return acceptMax*acceptMin;

	}

	int Individual::getFitness() {
		int fit = 0;
		int sz = (this->genes).size();

		for (int i = 0; i < sz; i++) {

			if (this->getGene(i) == (solution[i])) {
				fit++;
			}
		}
		this->fitness = fit;
		return fit;
	}

	double Individual::getFitnessForBSModel(MarketData md, DealData dd) {
		double fit = 0;
		md.sigma = this->target;

		fit = -abs(GeneticAlgo::convertBitToDouble(BSSqrDiffBitwise(md, dd)));

		this->fitnessDouble = fit;
		return fit;
	}

	// solution and precision setter
	void Individual::setSolution(boost::dynamic_bitset<> sol) {
		solution = sol;
		precision = sol.size();
	}

	// solution getter
	int Individual::getPrecision() {
		return precision;
	}



	Individual3D::Individual3D() {
		this->indivs = new Matrix<Individual*>(dim1, dim2);
		for (int i = 0; i < dim1; i++) {
			for (int j = 0; j < dim2; j++) {
				(*indivs)(i, j) = new Individual();
			}
		}
	}

	bool Individual3D::getGene(size_t idx1, size_t idx2, int genIndex) {
		bool gene = (*indivs)(idx1, idx2)->getGene(genIndex); //(bool)genes[index];
		return gene;
	}

	Matrix<Individual*>* Individual3D::getIndividualsMatrix() {
		return this->indivs;
	}

	double Individual3D::getTarget(size_t idx1, size_t idx2) {
		return (*indivs)(idx1, idx2)->getTarget();
	}

	Matrix<double>* Individual3D::getTargetMatrix() {
		Matrix<double> *val = new Matrix<double>(dim1, dim2);
		for (size_t i = 0; i < dim1; i++) {
			for (size_t j = 0; j < dim2; j++) {
				(*val)(i, j) = getTarget(i, j);
			}
		}

		return val;
	}

	/*void Individual3D::setGene(size_t d1, size_t d2, int genIndex, bool value, bool &stateMin, bool &stateMax) {

		bool isSign = (genIndex == 0);
		if (acceptGene(value, stateMin, stateMax, genIndex, isSign)) {
			(*indivs)(d1, d2).genes[genIndex] = value;
		}
		else {
			(*indivs)(d1, d2).genes[genIndex] = 1 - value;
		}
	}*/

	int Individual3D::getFitness() {
		int fit = 0;
		int sz = (*indivs).size(0);
		int sz2 = (*indivs).size(1);

		for (int i = 0; i < sz; i++) {
			for (int j = 0; j < sz2; j++) {
				for (int g = 0; g < precision; g++) {
					if (((*indivs)(i, j))->genes[g] == (solution[i])) {
						fit++;
					}
				}
			}
		}
		this->fitness = fit;
		return fit;
	}

	double Individual3D::getFitnessForBSModel(MarketData3D* md, DealData3D* dd) {
		double fit = 0;
		md->sigma = this->getTargetMatrix(); // setting value of sigma to market data to calculate differences

		fit = -abs(BSSqrDiff3D(md, dd));
		//fit = -abs(GeneticAlgo::convertBitToDouble(BSSqrDiffBitwise3D(md, dd))); // we expect fitness to be as close as possible to zero

		this->fitnessDouble = fit;
		return fit;
	}

	// solution and precision setter
	void Individual3D::setSolution(boost::dynamic_bitset<> sol) {
		solution = sol;
		precision = sol.size();
	}

	// solution getter
	int Individual3D::getPrecision() {
		return precision;
	}

	void Individual3D::setDimensionsOf3DIndividuals(size_t d1, size_t d2) {
		dim1 = d1;
		dim2 = d2;
	}







	//default values
	double GeneticAlgo::uniformRate = 0.5;
	double GeneticAlgo::mutationRate = 0.05;
	int GeneticAlgo::tournamentSize = 5;
	bool GeneticAlgo::elitism = 0;
	boost::dynamic_bitset<> Individual::solution (64, 0);
	int Individual::precision = 64;
	boost::dynamic_bitset<> GeneticAlgo::minset(64, 0);
	boost::dynamic_bitset<> GeneticAlgo::maxset(64, 1);
	double GeneticAlgo::minval = 0; // double constraint for sigma
	double GeneticAlgo::maxval = 1; // double constraint for sigma
	boost::dynamic_bitset<> Individual3D::solution(64, 0);
	int Individual3D::precision = 64;
	size_t Population3D::dim1 = 0;
	size_t Population3D::dim2 = 0;
	size_t Individual3D::dim1 = 0;
	size_t Individual3D::dim2 = 0;

	Population* GeneticAlgo::evolvePopulation(Population& pop) {
		int size = pop.size();
		Population* newPopulation = new Population(size);
		// Keep the best individual
		if (elitism) {
			newPopulation->addAnIndividual(pop.getFittest());
		}
		
		// add the new individuals to newPopulation individuals vector
		for (int i = 0; i < size-elitism; i++) {
			Individual* indiv1 = tournamentSelection(pop);
			Individual* indiv2 = tournamentSelection(pop);
			Individual* newIndiv = crossover(*indiv1, *indiv2);
			newPopulation->setIndividual(i, *newIndiv);
		}

		// Mutate population
		for (int i = 0; i < size; i++) {
			mutate(newPopulation->getIndividual(i));
		}

		return newPopulation;
	}
	
	// returns a value of new daughter individual
	Individual* GeneticAlgo::crossover(Individual& indiv1, Individual& indiv2) {
		Individual* newSol = new Individual();
		int precision = Individual::getPrecision();
		// Loop through genes
		bool stateMin = false; // initial state for constraints
		bool stateMax = false; // initial state for constraints
		for (int i = 0; i < precision; i++) {
			// Crossover
			if ((double)rand()/(double)RAND_MAX <= uniformRate) {
				newSol->setGene(i, indiv1.getGene(i), stateMin, stateMax);
			}
			else {
				newSol->setGene(i, indiv2.getGene(i), stateMin, stateMax);
			}
		}
		return newSol;
	}

	// Mutate an individual (changes the genes inside an individual
	void GeneticAlgo::mutate(Individual& indiv) {
		int precision = Individual::getPrecision();
		// Loop through genes
		bool stateMin = false; // initial state for constraints
		bool stateMax = false; // initial state for constraints
		for (int i = 0; i < precision; i++) {
			if ((double)rand()/(double)RAND_MAX <= mutationRate) {
				// Create random gene
				bool gene = (bool)round((double)rand()/(double)RAND_MAX);
				indiv.setGene(i, gene, stateMin, stateMax);
			}
		}
	}

	// takes some random individuals from population and chooses the fittest
	Individual* GeneticAlgo::tournamentSelection(Population& pop) {
		Individual* fittest = new Individual();
		int sz = pop.size();
		// Create a tournament population
		Population* tournament = new Population(tournamentSize);
		// For each place in the tournament get a random individual
		for (int i = 0; i < tournamentSize; i++) {
			int randomId = (int)(((double)rand()/(double)RAND_MAX) * (sz-1));
			tournament->setIndividual(i, pop.getIndividual(randomId));
		}
		// Get the fittest
		*fittest = tournament->getFittest();
		return fittest;
	}

	
	
	Population3D* GeneticAlgo::evolvePopulation(Population3D* pop) {
		int size = pop->size();
		
		Population3D* newPopulation = new Population3D(size);
		// Keep the best individual
		if (elitism) {
			newPopulation->addAnIndividual(pop->getFittest());
		}

		// add the new individuals to newPopulation individuals vector
		for (int i = 0; i < size - elitism; i++) {
			Individual3D* indiv1 = tournamentSelection(pop);
			Individual3D* indiv2 = tournamentSelection(pop);
			Individual3D* newIndiv = crossover(indiv1, indiv2);
			newPopulation->setIndividual(i, newIndiv);
		}

		// Mutate population
		for (int i = 0; i < size; i++) {
			mutate(newPopulation->getIndividual(i));
		}

		return newPopulation;
	}

	// returns a value of new daughter individual
	Individual3D* GeneticAlgo::crossover(Individual3D* indiv1, Individual3D* indiv2) {
		Individual3D* newSol = new Individual3D();
		int precision = Individual3D::getPrecision();
		Matrix<Individual*>* indivs1D_1 = indiv1->getIndividualsMatrix();
		Matrix<Individual*>* indivs1D_2 = indiv2->getIndividualsMatrix();
		Matrix<Individual*>* indivs1D_newSol = newSol->getIndividualsMatrix();
		size_t sz1 = Population3D::dim1;
		size_t sz2 = Population3D::dim2;

		for (int k = 0; k < sz1; k++) {
			
			for (int j = 0; j < sz2; j++) {

				// Loop through genes
				bool stateMin = false; // initial state for constraints
				bool stateMax = false; // initial state for constraints
				Individual* currentIndNewSol = (*indivs1D_newSol)(k, j);
				Individual* currentInd1 = (*indivs1D_1)(k, j);
				Individual* currentInd2 = (*indivs1D_2)(k, j);


				for (int i = 0; i < precision; i++) {
					// Crossover
					if ((double)rand() / (double)RAND_MAX <= uniformRate) {
						currentIndNewSol->setGene(i, currentInd1->getGene(i), stateMin, stateMax);
					}
					else {
						currentIndNewSol->setGene(i, currentInd2->getGene(i), stateMin, stateMax);
					}
				}

			}
		}

		return newSol;
	}

	// Mutate an individual (changes the genes inside an individual
	void GeneticAlgo::mutate(Individual3D* indiv) {
		int precision = Individual3D::getPrecision();
		Matrix<Individual*>* indivs1D = indiv->getIndividualsMatrix();
		size_t sz1 = Population3D::dim1;
		size_t sz2 = Population3D::dim2;

		for (int k = 0; k < sz1; k++) {
			for (int j = 0; j < sz2; j++) {
				// Loop through genes
				bool stateMin = false; // initial state for constraints
				bool stateMax = false; // initial state for constraints
				Individual* currentInd = (*indivs1D)(k, j);

				for (int i = 0; i < precision; i++) {
					if ((double)rand() / (double)RAND_MAX <= mutationRate) {
						// Create random gene
						bool gene = (bool)round((double)rand() / (double)RAND_MAX);
						currentInd->setGene(i, gene, stateMin, stateMax);
					}
				}
			}
		}
		
	}

	// takes some random individuals from population and chooses the fittest
	Individual3D* GeneticAlgo::tournamentSelection(Population3D* pop) {
		Individual3D* fittest = new Individual3D();
		int sz = pop->size();
		// Create a tournament population
		Population3D* tournament = new Population3D(tournamentSize);
		// For each place in the tournament get a random individual
		for (int i = 0; i < tournamentSize; i++) {
			int randomId = (int)(((double)rand() / (double)RAND_MAX) * (sz - 1));
			tournament->setIndividual(i, pop->getIndividual(randomId));
		}
		// Get the fittest
		fittest = tournament->getFittest();
		return fittest;
	}





	// initialize static algo input
	void  GeneticAlgo::initializeAlgoInput(double uniformRate, double mutationRate, int tournamentSize, bool elitism) {
		GeneticAlgo::uniformRate = uniformRate;
		GeneticAlgo::mutationRate = mutationRate;
		GeneticAlgo::tournamentSize = tournamentSize;
		GeneticAlgo::elitism = elitism;

	}

	// setter for constraints
	void GeneticAlgo::setSystemBinaryConstraints(boost::dynamic_bitset<> valmin, boost::dynamic_bitset<> valmax) {
		minset = valmin;
		maxset = valmax;
		minval = convertBitToDouble(valmin);
		maxval = convertBitToDouble(valmax);
	}
	
	void GeneticAlgo::setSystemDoubleConstraints(double valmin, double valmax) {
		minval = valmin;
		maxval = valmax;
		minset = convertDoubleTo64Bit(valmin);
		maxset = convertDoubleTo64Bit(valmax);
	}
	

	// binary functionality for constraints
	boost::dynamic_bitset<> GeneticAlgo::convertDoubleTo64Bit(double value) {

		boost::dynamic_bitset<> result(64, 0);
		bool sign;
		// get sign
		if (value < 0) {
			sign = 1;
		}
		else {
			sign = 0;
		}

		// get the integral and fractional part values
		double integralPart, fractionalPart;
		fractionalPart = std::modf(value, &integralPart);

		// convert integral and fractional parts to bits
		boost::dynamic_bitset<> bitIntegralPart = convertIntToBit(integralPart);
		boost::dynamic_bitset<> bitFractionalPart = convertFractionToBit(fractionalPart);
		std::vector<bool> testInt = convertBitsetToVector(bitIntegralPart);
		std::vector<bool> testFloat = convertBitsetToVector(bitFractionalPart);

		
		// getting binary exponent (11) and mantissa (52) values
		boost::dynamic_bitset<> mantissa (52,0);
		boost::dynamic_bitset<> exponent = getExponentMantissaByNormalization(bitIntegralPart, bitFractionalPart, mantissa);
		std::vector<bool> testExp = convertBitsetToVector(exponent);
		std::vector<bool> testMent = convertBitsetToVector(mantissa);

		// setting resulting values
		result[0] = sign;
		int k = 1;
		for (int i = 0; i < 11; i++) {
			result[k] = exponent[i];
			k++;
		}
		for (int i = 0; i < 52; i++) {
			result[k] = mantissa[i];
			k++;
		}


		return result;
	}

	boost::dynamic_bitset<> GeneticAlgo::convertIntToBit(int value) {
		int res = value;
		boost::dynamic_bitset<> temp;
		if (res == 0) {
			return boost::dynamic_bitset<>(1, 0);
		}

		while (res != 0) {
			int rem = res % 2;
			res = res / 2;
			temp.push_back(rem);
		}

		int sz = temp.size();
		boost::dynamic_bitset<> bitRepresantation(sz,0);
		for (int i = 0; i < sz; i++) {
			bitRepresantation[i] = temp[sz - 1 - i];
		}

		return bitRepresantation;
	}

	boost::dynamic_bitset<> GeneticAlgo::convertIntTo11Bit(int value) {
		int res = value;
		boost::dynamic_bitset<> temp;
		if (res == 0) {
			return boost::dynamic_bitset<>(11, 0);
		}

		while (res != 0) {
			int rem = res % 2;
			res = res / 2;
			temp.push_back(rem);
		}

		int sz = temp.size();


		boost::dynamic_bitset<> bitRepresantation(11, 0);
		int i = 0;
		while (i < 11-sz) {
			bitRepresantation[i] = 0;
			i++;
		}
		while (i < 11) {
			bitRepresantation[i] = temp[11-i-1];
			i++;
		}

		return bitRepresantation;
	}

	boost::dynamic_bitset<> GeneticAlgo::convertFractionToBit(double value) {
		double res = value;
		boost::dynamic_bitset<> bitRepresentation;
		if (res == 0) {
			return boost::dynamic_bitset<>(1, 0);
		}


		while (res != 1) {
			res = res * 2;
			if (res > 1) {
				bitRepresentation.push_back(1);
				res = res - 1;
			}
			else if (res < 1){
				bitRepresentation.push_back(0);
			}
			else {
				bitRepresentation.push_back(1);
				break;
			}
		}

		return bitRepresentation;

	}

	boost::dynamic_bitset<> GeneticAlgo::getExponentMantissaByNormalization(boost::dynamic_bitset<> bitIntegralPart, boost::dynamic_bitset<> bitFractionalPart, boost::dynamic_bitset<> &mantissa) {
		int szInt = bitIntegralPart.size();
		int szFr = bitFractionalPart.size();
		int posFromSeparator = 0; // power of 2
		boost::dynamic_bitset<> exponent(11, 0);

		// check the position from the left
		for (int i = 0; i < szInt; i++) {
			if (bitIntegralPart[i] == 1) {
				posFromSeparator = szInt - 1 - i;
				int k = 0;
				for (int i = szInt-1-posFromSeparator; i < szInt - 1; i++) {
					mantissa[k] = bitIntegralPart[i + 1];
					k++;
				}
				for (int i = 0; i < szFr; i++) {
					mantissa[k] = bitFractionalPart[i];
					k++;
				}
				for (k; k < 52; k++) {
					mantissa[k] = 0;
					k++;
				}
				exponent = convertIntTo11Bit(1023+posFromSeparator); 
				break;
			}
		}

		// if we have 0 on the left check the position from the right
		if (szInt == 1 && bitIntegralPart[0] == 0) {
			for (int i = 0; i < szFr; i++) {
				if (bitFractionalPart[i] == 1) {
					posFromSeparator = -(i+1);
					int k = 0;
					for (int i = -posFromSeparator; i < szFr; i++) {
						mantissa[k] = bitFractionalPart[i];
						k++;
					}
					for (k; k < 52; k++) {
						mantissa[k] = 0;
					}
					exponent = convertIntTo11Bit(1023 + posFromSeparator);
					break;
				}
			}

		}

		return exponent;
		
		
	}

	int GeneticAlgo::convertBitToInt(boost::dynamic_bitset<> value) {
		int sz = value.size();
		
		int sum = 0;
		for (int i = 0; i < sz; i++) {
			sum = sum + value[i]*pow(2, sz - 1 - i);
		}
		return sum;
	}

	double GeneticAlgo::convertBitToFraction(boost::dynamic_bitset<> value) {
		int sz = value.size();
		double res = 1.0;
		for (int i = sz - 1; i > 0; i--) {
			res = res / 2;
			if (value[i-1] == 1) {
				res = 1 + res;
			}
		}
		res = res / 2;

		return res;
	}

	double GeneticAlgo::convertBitToDouble(boost::dynamic_bitset<> value) {
		int sign = 1;
		if (value[0] == 1) { sign = -1; }

		boost::dynamic_bitset<> exponent(11, 0);
		boost::dynamic_bitset<> mantissa(52, 0);
		boost::dynamic_bitset<> binIntegralPart;
		boost::dynamic_bitset<> binFloatingPart;
		int integralPart = 0;
		double floatingPart = 0;

		std::vector<bool>testexp(11,0);
		std::vector<bool>testmant(52, 0);
		
		int k = 0;
		for (int i = 1; i < 12; i++) {
			exponent[k] = value[i];
			k++;
		}

		k = 0;
		for (int i = 12; i < 64; i++) {
			mantissa[k] = value[i];
			k++;
		}
		
		// testing obtained exp and mantissa
		testexp = convertBitsetToVector(exponent);
		testmant = convertBitsetToVector(mantissa);


		int temp = convertBitToInt(exponent);
		if (temp == 0) {
			return 0;
		}
		int powerOfTwo = temp - 1023;

		if (powerOfTwo >= 0) {
			binIntegralPart.push_back(1);
			for (int i = 0; i < powerOfTwo; i++) {
				binIntegralPart.push_back(mantissa[i]);
			}
			for (int i = powerOfTwo; i < 52; i++) {
				binFloatingPart.push_back(mantissa[i]);
			}
		}
		else {
			binIntegralPart.push_back(0);
			for (int i = 0; i < -powerOfTwo-1; i++) {
				binFloatingPart.push_back(0);
			}
			binFloatingPart.push_back(1);
			for (int i = 0; i < 52; i++) {
				binFloatingPart.push_back(mantissa[i]);
			}
		}
		
		integralPart = convertBitToInt(binIntegralPart);
		floatingPart = convertBitToFraction(binFloatingPart);

		return sign*(integralPart + floatingPart);

	}

	std::vector<bool> GeneticAlgo::convertBitsetToVector(boost::dynamic_bitset<> array){
		int sz = array.size();
		std::vector<bool> vec(sz, 0);
		for (int k = 0; k < sz; k++) {
			vec[k] = array[k];
		}
		return vec;
	}



	


	// Financial module
	
	boost::dynamic_bitset<> BSSqrDiffBitwise(MarketData md, DealData dd) {
		double S = md.S;
		double r = md.r;
		double K = dd.K;
		double sigma = md.sigma;
		std::vector<double> T = dd.T;
		double sumOfTheSqrDifference = 0;

		int sz = T.size();

		for (int i = 0; i < sz; i++) {

			double d1 = (1 / (sigma * sqrt(T[i])))*(log(S / K) + (r + pow(sigma, 2) / 2)*T[i]);
			//cout << "d1 = " << d1.getValue() << endl;

			double d2 = d1 - sigma * sqrt(T[i]);
			//cout << "d2 = " << d2.getValue() << endl;

			double price = (NormalCDFCody(d1)*S) - (NormalCDFCody(d2)*K*exp(-r*T[i]));

			sumOfTheSqrDifference = sumOfTheSqrDifference + pow((price - md.prices[i]),2);
		}

		boost::dynamic_bitset<> x = GeneticAlgo::convertDoubleTo64Bit(sumOfTheSqrDifference);
		
		std::vector<bool> test = GeneticAlgo::convertBitsetToVector(x); // testing results

		return x;

		
	}

	boost::dynamic_bitset<> BSSqrDiffBitwise3D(MarketData3D* md, DealData3D* dd) {

		Matrix<double> price = FDMLocalVolpricer(md, dd);

		double sumOfTheSqrDifference = ((price - *(md->prices))^2).sumOfElements();
	
		boost::dynamic_bitset<> x = GeneticAlgo::convertDoubleTo64Bit(sumOfTheSqrDifference);

		std::vector<bool> test = GeneticAlgo::convertBitsetToVector(x); // testing results

		return x;
		
	}

	double BSSqrDiff3D(MarketData3D* md, DealData3D* dd) {

		Matrix<double> price = FDMLocalVolpricer(md, dd);

		double sumOfTheSqrDifference = ((price - *(md->prices)) ^ 2).sumOfElements();

		return sumOfTheSqrDifference;

	}



	// Matrix fill and theta-schema utilities
	std::vector<double> thomasAlgo(Matrix<double>& A, Matrix<double>& B) {
		int szRow = A.size(0);

		std::vector<double> res (szRow, 0);
		std::vector<double> c (szRow, 0);
		std::vector<double> d (szRow, 0);
		//B(0, 0) = 110;
		c[0] = A(0, 1) / A(0, 0);
		d[0] = B(0, 0) / A(0, 0);

		// forward walk: looking for transformed A matrix members
		for (int i = 1; i < szRow-1; i++) {
			c[i] = A(i, i + 1) / (A(i, i) - c[i - 1] * A(i, i - 1));
			d[i] = (B(0, i) - d[i - 1] * A(i, i - 1)) / (A(i, i) - c[i - 1] * A(i, i - 1));
		}
		d[szRow - 1] = (B(0, szRow - 1) - d[szRow - 1 - 1] * A(szRow - 1, szRow - 1 - 1)) / (A(szRow - 1, szRow - 1) - c[szRow - 1 - 1] * A(szRow - 1, szRow - 1 - 1));
		
		// backward walk: finding the solutions
		res[szRow - 1] = d[szRow - 1];

		for (int i = szRow - 1; i > 0; i--) {
			res[i - 1] = d[i - 1] - res[i] * c[i - 1];
		}

		return res;
	}

	Matrix<double> createATriagonalMatrix(Matrix<double>& alpha, Matrix<double>& beta, Matrix<double>& gamma, int timeIndex) {
		int sz = alpha.size(1)-2;
		Matrix<double> triagonal (sz, sz, 0);

		triagonal(0, 0) = beta(timeIndex, 0);
		triagonal(0, 1) = gamma(timeIndex, 0);
		triagonal(sz - 1, sz - 1) = beta(timeIndex, sz - 1);
		triagonal(sz - 1, sz - 2) = alpha(timeIndex, sz - 1);

		for (int i = 1; i < sz - 1; i++) {
			triagonal(i, i) = beta(timeIndex, i);
			triagonal(i, i - 1) = alpha(timeIndex, i);
			triagonal(i, i + 1) = gamma(timeIndex, i);
		}

		return triagonal;

	}

	void fillAlphaBetaGammaFromSigmaDeltaT(Matrix<double>& alpha, Matrix<double>& beta, Matrix<double>& gamma, Matrix<double>& sigma, double r, double delta_K) {
		int d1 = alpha.size(0);
		int d2 = alpha.size(1);


		alpha = ((sigma)^2)*(-1 / (2*pow(delta_K, 2)) - 1 / (4 * delta_K)) - r / (2 * delta_K);
		//alpha.pop({ 0,d1-1 }, { 0,d2-1 });

		beta = ((sigma)^2) / (delta_K*delta_K);
		//beta.pop({ 0,d1-1 }, { 0,d2-1 });

		gamma = ((sigma)^2)*(-1 / (2*pow(delta_K, 2)) + 1 / (4 * delta_K)) + r / (2 * delta_K);
		//gamma.pop({ 0,d1-1 }, { 0,d2-1 });
	}


	// Pricers

	Matrix<double> FDMLocalVolpricer(MarketData3D* md, DealData3D* dd) {
		
		int discretization_num_T = dd->discretization_num_T; //will be the same as size of K,T not to interpolate in the beginning
		int discretization_num_K = dd->discretization_num_K;
		double S = md->S;
		double r = md->r;
		std::vector<double> K = dd->K;
		std::vector<double> T = dd->T;
		Matrix<double>* sigma = md->sigma;
		double sumOfTheSqrDifference = 0;
		double delta_T = (T[discretization_num_T-1] - T[0]) / (discretization_num_T-1);
		double delta_K = (K[discretization_num_K-1] - K[0]) / (discretization_num_K-1);
		
		Matrix<double> u (discretization_num_T, discretization_num_K);
		// setting initial condition at T = t0
		for (int j = 0; j < discretization_num_K; j++) {
			u(0, j) = std::max(S - exp(-r*T[0])*K[j], double(0));
		}

		// setting boundary conditions for K = K0; K = end; 
		//if K>>S0 -> u(:,K) -> 0     if K<<S0 -> u(T(i),K) -> S0-K*exp(-r*T(i)) for all i
		for (int i = 0; i < discretization_num_T; i++) {
			u(i, 0) = S - K[0]; // Dirichlet boundary Call
			u(i, discretization_num_K - 1) = 0; // Dirichlet boundary Call
		}

		// Solution of PDE
		for (int i = 0; i < discretization_num_T-1; i++) {
			for (int j = 1; j < discretization_num_K-1; j++) {
				u(i + 1, j) = u(i, j) + 0.5*pow((*sigma)(i, j),2)*(pow(K[i],2))*(u(i, j + 1) + u(i, j - 1) - 2 * (u(i, j)))*delta_T / (pow(delta_K, 2)) - r*K[i] * (u(i, j + 1) - u(i, j))*delta_T / delta_K;
			}
		}

		return u;

	}

	Matrix<double> FDMLocalVolpricerThetaScheme(MarketData3D* md, DealData3D* dd, double theta) {
		int discretization_num_T = dd->discretization_num_T; //will be the same as size of K,T not to interpolate in the beginning
		int discretization_num_K = dd->discretization_num_K;
		double S = md->S;
		double r = md->r;
		std::vector<double> K = dd->K;
		std::vector<double> T = dd->T;
		Matrix<double>* sigma = md->sigma;
		double sumOfTheSqrDifference = 0;
		double delta_T = (T[discretization_num_T - 1] - T[0]) / (discretization_num_T - 1);
		double delta_K = (K[discretization_num_K - 1] - K[0]) / (discretization_num_K - 1);
		double delta_log_K = (log(K[discretization_num_K - 1]) - log(K[1])) / (discretization_num_K - 1);

		Matrix<double> u (discretization_num_T, discretization_num_K);
		// setting initial condition at T = t0
		for (int j = 0; j < discretization_num_K; j++) {
			u(0, j) = std::max(S - K[j], double(0));
		}
		// setting initial condtion for K = k0
		for (int i = 0; i < discretization_num_T; i++) {
			u(i, 0) = S;
		}


		Matrix<double> B (1, discretization_num_K-2);
		Matrix<double> alpha(discretization_num_T, discretization_num_K);
		Matrix<double> beta (discretization_num_T, discretization_num_K);
		Matrix<double> gamma (discretization_num_T, discretization_num_K);


		// filling the matrix of coefficients alpha, beta, gamma
		fillAlphaBetaGammaFromSigmaDeltaT(alpha, beta, gamma, *sigma, r, delta_log_K);
		
		if (theta == 0) {
			for (int n = 0; n < discretization_num_T - 1; n++) {
				//caclulating inputs for thomas algo to find u(n+1,:)
				Matrix<double> A = createATriagonalMatrix(alpha, beta, gamma, n + 1); // An+1 = 38*38
				Matrix<double> H = A*delta_T + Matrix<double>(A.size(0), A.size(1), char('id')); // Hn
				for (int k = 0; k < discretization_num_K-2; k++) {
					B(0, k) = u(n, k+1) - (alpha)(n + 1, 0)*S*delta_T*std::max(1-k,0); // Bn
				}

				// Run thomas algo
				std::vector<double> res = thomasAlgo(H, B);

				// write results into our solution matrix
				for (int i = 1; i < discretization_num_K-1; i++) {
					u(n + 1, i) = res[i-1];
				}
			}
		}


		return u;

	}

	
	// Additional functionality

	Matrix<double>* BSPriceMatrixCreator(MarketData3D* md, DealData3D* dd) {
		double S = md->S;
		double r = md->r;
		std::vector<double> K = dd->K;
		std::vector<double> T = dd->T;
		

		int sz1 = T.size();
		int sz2 = K.size();
		Matrix<double> sigma (sz1, sz2, 0.2);
		Matrix<double>* prices = new Matrix<double>(sz1, sz2);

		for (int i = 0; i < sz1; i++) {

			for (int j = 0; j < sz2; j++) {

				double test = log(S / K[j]);
				double test2 = (r + pow(sigma(i, j), 2) / 2)*T[i];
				double test3 = sqrt(T[i]);
				
				double d1 = (1 / (sigma(i,j) * sqrt(T[i])))*(log(S / K[j]) + (r + pow(sigma(i,j), 2) / 2)*T[i]);
				//cout << "d1 = " << d1.getValue() << endl;
				if (isnan(d1)){
					d1 = INFINITY;
				}

				double d2 = d1 - sigma(i,j) * sqrt(T[i]);
				//cout << "d2 = " << d2.getValue() << endl;

				double price = (NormalCDFCody(d1)*S) - (NormalCDFCody(d2)*K[j]*exp(-r*T[i]));

				(*prices)(i, j) = price;

			}

		}

		return prices;

	}

	double NormalCDFCody(double u) {
		double y = abs(u);
		if (y > 35.0) {
			if (u > 0)
				return 1;
			else
				return 0;
		}
		if (y <= 0.662912607) {
			//  evaluate erf() for |u| <= sqrt(2)*0.46875
			double a0 = 1.161110663653770e-2;
			double a1 = 3.951404679838207e-1;
			double a2 = 2.846603853776254e+1;
			double a3 = 1.887426188426510e+2;
			double a4 = 3.209377589138469e+3;

			double b0 = 1.767766952966369e-1;
			double b1 = 8.344316438579620;
			double b2 = 1.725514762600375e+2;
			double b3 = 1.813893686502485e+3;
			double b4 = 8.044716608901563e+3;

			double z = y * y;
			y = u * ((((a0 * z + a1) * z + a2) * z + a3) * z + a4);
			y /= ((((b0 * z + b1) * z + b2) * z + b3) * z + b4);
			return 0.5 + y;
		}
		double zinterm = 0.5 * exp(-y * y / 2);
		if (y <= 4.0) {
			double c0 = 2.15311535474403846e-8;
			double c1 = 5.64188496988670089e-1;
			double c2 = 8.88314979438837594;
			double c3 = 6.61191906371416295e+1;
			double c4 = 2.98635138197400131e+2;
			double c5 = 8.81952221241769090e+2;
			double c6 = 1.71204761263407058e+3;
			double c7 = 2.05107837782607147e+3;
			double c8 = 1.23033935479799725e+3;

			double d0 = 1.0;
			double d1 = 1.57449261107098347e+1;
			double d2 = 1.17693950891312499e+2;
			double d3 = 5.37181101862009858e+2;
			double d4 = 1.62138957456669019e+3;
			double d5 = 3.29079923573345963e+3;
			double d6 = 4.36261909014324716e+3;
			double d7 = 3.43936767414372164e+3;
			double d8 = 1.23033935480374942e+3;

			// evaluate erfc() for sqrt(2)*0.46875 <= |u| <= sqrt(2)*4.0
			y = y / 1.4142135623730950488;
			double num = ((((((((c0 * y + c1) * y + c2) * y + c3) * y + c4) * y + c5) * y + c6) * y + c7) * y + c8);
			double den = ((((((((d0 * y + d1) * y + d2) * y + d3) * y + d4) * y + d5) * y + d6) * y + d7) * y + d8);

			y = num / den;
			y = zinterm * y;
		}
		else {
			double p0 = 1.63153871373020978e-2;
			double p1 = 3.05326634961232344e-1;
			double p2 = 3.60344899949804439e-1;
			double p3 = 1.25781726111229246e-1;
			double p4 = 1.60837851487422766e-2;
			double p5 = 6.58749161529837803e-4;

			double q0 = 1.00000000000000000;
			double q1 = 2.56852019228982242;
			double q2 = 1.87295284992346047;
			double q3 = 5.27905102951428412e-1;
			double q4 = 6.05183413124413191e-2;
			double q5 = 2.33520497626869185e-3;
			// evaluate erfc() for |u| > sqrt(2)*4.0
			double z = zinterm * 1.41421356237309504880 / y;
			y = 2 / (y * y);
			y = y * (((((p0 * y + p1) * y + p2) * y + p3) * y + p4) * y + p5) / (((((q0 * y + q1) * y + q2) * y + q3) * y + q4) * y + q5);
			y = z*(0.564189583547756287 - y);
		}

		if (u < 0)
			return y;
		return 1 - y;
	}

	double convert(boost::dynamic_bitset<> const& bs) {
		std::bitset<64> myBit(0);
		for (boost::dynamic_bitset<>::size_type i = 0; i < bs.size(); i++) {
			myBit[i] = bs[i];
		}		
		//static_assert(sizeof(uint64_t) == sizeof(double), "Cannot use this!");
		
		uint64_t const u = myBit.to_ullong(); // returns long long const value
		double d = u;

		// Aliases to `char*` are explicitly allowed in the Standard (and only them)
		char const* cu = reinterpret_cast<char const*>(&u); // long long cons value is rewritten as char
		char* cd = reinterpret_cast<char*>(&d); 

		// Copy the bitwise representation from u to d
		memcpy(cd, cu, sizeof(u)); 
		

		return d;
	}

	

}