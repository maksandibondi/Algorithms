#include "stdafx.h"
#include "CppUnitTest.h"
#include "AlgoUtilities.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace AlgoUtilities;

namespace UnitTest
{		
	TEST_CLASS(UnitTest1)
	{
	public:
		
		TEST_METHOD(TestMethod1)
		{
			boost::dynamic_bitset<> sol(64, 0);
			Individual::setSolution(sol);
			// TODO: Your test code here
		}

	};
}