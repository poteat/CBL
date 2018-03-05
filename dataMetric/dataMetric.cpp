#include <assert.h>
#include <string>
#include <functional>
#include <sstream>

#include "axis.h"
#include "mrc.h"

using namespace cbl;

cbl::real applyLineData(mrc map, pdb structure, cbl::real deviation, std::string& data)
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

	auto variance = [](std::vector<cbl::real> v, cbl::real avg)
	{
		cbl::real variance_sum = 0;
		auto add = [&](cbl::real x) {variance_sum += pow((x - avg), 2); };
		std::for_each(v.begin(), v.end(), add);
		return (pow((variance_sum / v.size()), .5));
	};

	line.calcError(map, avg);

	std::ostringstream out;

	line.write(out, deviation);

	data += out.str();

	return variance(line.error, avg(line.error));
}

void cylinderCutOut(mrc map, pdb structure)
{
	//Chop out the density around helixes using a cylinder of 5-6 angstroms
	std::string dist = "5";
	cbl::real cropping_dist = (cbl::real) std::stod(dist);

	mrc near, far;

	std::tie(near, far) = map.crop(structure, cropping_dist);

	map = near;
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

	std::ofstream score_out_file("quantification.dat");
	std::vector<double>threshold_score;
	int loopy = 0;
	double final_variance = 0;

	for (float t = 0.2f; t <= 2.0f; t += 0.2f)
	{
		threshold_score.push_back(applyLineData(map, structure, t, data_result));
		score_out_file << t << "\t" << threshold_score[loopy] << std::endl;
		loopy++;
	}

	score_out_file << "\n\nUtilized scores\n";
	std::sort(threshold_score.begin(), threshold_score.end());
	for (int i = 0; i <= (loopy / 2); i++)
	{
		final_variance += threshold_score[i];
		score_out_file <<threshold_score[i] << std::endl;
	}

	final_variance = (final_variance / (loopy / 2));

	score_out_file << "----------------------------------------\n" << "Final Score: " << final_variance;
	score_out_file.close();


	std::ofstream out_file("scatter.dat");

	out_file << data_result << std::endl;

	out_file.close();

	system("display_scatter_plot.bat");

	return 0;
}