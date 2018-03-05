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
void cylinderCutOut(mrc map, pdb structure)
{
	//Chop out the density around helixes using a cylinder of 5-6 angstroms
	cbl::real cropping_dist = (cbl::real) std::stod();

	mrc near, far;

	std::tie(near, far) = map.crop(structure, cropping_dist);

	map = near;
}
double score(mrc map, pdb structure, float threshold)
{
	//find variance and average the bottom 50%
	double quantification = 0;








	return quantification;
}

int main(int argc, char* argv[])
{
	int num_args = 2;
	assert(argc == num_args + 1 && "Incorrect num arguments, 2 inputs: mrc, pdb");

	mrc map(argv[1]);
	pdb structure(argv[2]);
	std::string data_result;

	cylinderCutOut(map, structure);

	map.write("generatedMRC.mrc");

	map.normalize();

	for (float t = 0.2f; t <= 2.0f; t += 0.2f)
	{
		applyLineData(map, structure, t, data_result);
	}

	std::ofstream out_file("scatter.dat");

	out_file << data_result << std::endl;

	out_file.close();

	system("display_scatter_plot.bat");

	std::ofstream score_out_file("quantification.dat");

	double quantification = 0;
	int loopy = 0;

	for (float threshold = 0; threshold <= 1; threshold+= 0.1)
	{
		double threshold_score = score(map, structure, threshold);
		quantification += threshold_score;
		loopy++;
		score_out_file << threshold << "\t" << threshold_score << std::endl;
	}

	quantification = quantification / loopy;
	score_out_file << "----------------------------------------\n" << "Final Score: " << quantification;

	return 0;
}