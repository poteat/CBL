#include <assert.h>
#include <string>
#include <functional>

#include "axis.h"

using namespace cbl;

int main(int argc, char* argv[])
{
	int num_args = 2;
	assert(argc == num_args + 1 && "Incorrect num arguments, 2 inputs: mrc, pdb");

	mrc map(argv[1]);
	pdb structure(argv[2]);

	axis line(structure);

	map.applyDeviationThreshold(2);
	map.setHollowBoundary();

	auto avg = [](std::vector<cbl::real> v)
	{
		cbl::real sum;
		auto add = [&sum](auto x) {sum += x; };
		std::for_each(v.begin(), v.end(), add);
		return sum / v.size();
	};

	line.calcError(map, avg);
	// line.removeZeroErrors();
	line.write("scatter.dat", 2);

	system("pause");

	return 0;
}