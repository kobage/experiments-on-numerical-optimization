#include <ctime>
#include "methods.h"
#include "testFunctions.h"

using namespace hb;

int main()
{
	makeTestsVector();
	runHBTest(testsVector);
}