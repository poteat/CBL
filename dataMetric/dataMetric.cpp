#include <assert.h>
#include <string>
#include <functional>
#include <sstream>

#include "axis.h"
#include "mrc.h"

using namespace cbl;

void applyLineData(mrc map, pdb structure, cbl::real deviation, std::string& data)
{
	axis line(structure);

	map.applyDeviationThreshold(deviation);
	map.setHollowBoundary();

	auto avg = [](std::vector<cbl::real> v)
	{
		cbl::real sum = 0;
		auto add = [&sum](cbl::real x) {sum += x; };
		std::for_each(v.begin(), v.end(), add);
		return sum / v.size();
	};

	line.calcError(map, avg);

	std::ostringstream out;

	line.write(out, deviation);

	data += out.str();
}

int main(int argc, char* argv[])
{
	int num_args = 2;
	assert(argc == num_args + 1 && "Incorrect num arguments, 2 inputs: mrc, pdb");

	mrc map(argv[1]);
	pdb structure(argv[2]);
	std::string data_result;

	map.normalize();

	for (float t = 0.2f; t <= 2.0f; t += 0.2f)
	{
		applyLineData(map, structure, t, data_result);
	}

	std::ofstream out_file("scatter.dat");

	out_file << data_result << std::endl;

	out_file.close();

	system("display_scatter_plot.bat");

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



	return 0;
}