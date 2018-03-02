#include <assert.h>
#include <string>
#include <functional>
#include <sstream>

#include "axis.h"

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

	return 0;
}