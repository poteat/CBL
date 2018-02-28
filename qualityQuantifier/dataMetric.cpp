#include <assert.h>
#include <string>
#include <functional>
#include <iostream>
#include <iomanip>

#include "axis.h"
#include "mrc.h"

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

	map.normalize();

	//take regular intervals of the threshold from 1 - 0.5
	//continue running smaller thresholds, if voxels are greater than a certain amount, discount them and stop decreasing
	double threshold;
	threshold = 1.0;

	for (int i = 20; i > -1; i--)
	{
		map.applyThreshold(threshold);
		threshold -= 0.05;

		//take slices of helixes to form circles of density and find av distance from center, and then av distance from true axis
		//use new CBL function for the border mrc points







	}

	//At the end, output best threshold into file along with helix and beta sheet comparisons
	std::string mrc_name = argv[1];
	std::ofstream outFile;
	for (int i = 0; i < 4; i++)
	{
		mrc_name.pop_back;
	}
	mrc_name.append(".txt");


	//---------------------------------------------------------------------------------------------------------------------------
	//Example Use of CBL library for reference later
	//mrc temp = mrc("foo.mrc");
	//temp.print();
	//temp.trim(2);
	//temp.applyThreshold(.019);
	//temp.write("bar.mrc");
	//system("pause");

	outFile.close();

	return 0;
}