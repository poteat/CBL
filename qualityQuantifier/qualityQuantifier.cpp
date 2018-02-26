#include <assert.h>
#include <string>
#include <functional>

#include "off.h"

using namespace cbl;

int main(int argc, char* argv[])
{
	int num_args = 1;
	std::cout << argc << std::endl;
	assert(argc == num_args + 1 && "Incorrect num arguments, 1 input needed: mrc filename");

	mrc input(argv[1]);

	input.applyDeviationThreshold(2);
	input.setHollowBoundary();

	input.write("foo.mrc");

	return 0;
}