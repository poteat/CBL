#include <assert.h>
#include <string>
#include <functional>
#include <sstream>
#include <experimental/filesystem>

#include "axis.h"
#include "mrc.h"

using namespace cbl;

cbl::real applyLineData(mrc map, pdb structure, cbl::real deviation, std::string& data)
{
	axis line(structure);

	cbl::real deviation_value = map.map.standard_deviation();

	cbl::real threshold_used = deviation_value * deviation;

	std::cout << "# Deviations: " << deviation << "     Threshold Used: " << threshold_used << std::endl;

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
		auto add = [&](cbl::real x) {variance_sum += pow(x - avg, 2); };
		std::for_each(v.begin(), v.end(), add);
		return sqrt(variance_sum / v.size());
	};

	line.calcError(map, avg);

	line.mergeError(1.0); // Angstroms

	std::ostringstream out;

	line.write(out, deviation);

	data += out.str();

	return (cbl::real) variance(line.error, avg(line.error));
}

mrc cylinderCutOut(mrc &map, pdb &structure)
{
	//Chop out the density around helixes using a cylinder of 5-6 angstroms
	cbl::real cropping_dist = (cbl::real) 5;

	mrc near, far;

	std::tie(near, far) = map.cylinderCrop(structure, cropping_dist);

	near.minimize();

	// To reduce Chimera blockiness visual effect
	// near.pad(5);

	return near;
}

std::vector<pdb> runAxisComparison(std::string pdb_path)
{
	using path = std::experimental::filesystem::path;

	// Convert relative to absolute path if necessary, using experimental std lib
	// If there are symbolic links, we also convert those to canonical form
	path symbolic_pdb_path = pdb_path;
	symbolic_pdb_path = std::experimental::filesystem::canonical(symbolic_pdb_path);
	pdb_path = symbolic_pdb_path.string();

	// Assert that the extension is ".pdb"
	assert(symbolic_pdb_path.extension() == ".pdb");

	// Make output directory for results to be populated
	path output_path = symbolic_pdb_path.parent_path().string() + "/output";
	std::experimental::filesystem::create_directory(output_path);

	// Get path stripped of file extension for input into axis-comparsion
	path stripped_path = symbolic_pdb_path.parent_path().string() + "/" + 
		symbolic_pdb_path.stem().string();

	std::string command = "leastsquare.exe \"" + stripped_path.string()
		+ "\" \"Empty\" \"Empty\" \"Empty\"";

	// Execute axis comparison
	system(command.c_str());

	// Read all files in this directory that begin with "trueHelix"

	std::vector<cbl::pdb> matched_helix;

	for (auto &p : std::experimental::filesystem::directory_iterator(output_path))
	{
		auto s = p.path().filename().string();

		// If begins with "trueHelix"
		if (s.find("trueHelix") == 0)
		{
			matched_helix.emplace_back(p.path().string());
		}
	}

	// Rename output directory to something we can refer to later

	path new_output_path = symbolic_pdb_path.parent_path().string() + "/"
		+ symbolic_pdb_path.stem().string();

	std::experimental::filesystem::remove_all(new_output_path);

	std::experimental::filesystem::rename(output_path, new_output_path);

	return matched_helix;
}

int main(int argc, char* argv[])
{
	int num_args = 2;
	assert(argc == num_args + 1 && "Incorrect num arguments, 2 inputs: mrc, pdb");

	std::string mrc_file_path_in = argv[1];
	std::string pdb_file_path_in = argv[2];
	
	std::vector<pdb> helices = runAxisComparison(pdb_file_path_in);

	mrc entire_map(mrc_file_path_in);

	// Build stripped mrc file_name for creating out filenames
	size_t period_pos = mrc_file_path_in.find_last_of('.');
	mrc_file_path_in.resize(period_pos);
	mrc_file_path_in += "_helix";

	for (size_t i = 0; i < helices.size(); i++)
	{
		pdb& helix = helices[i];
		std::string data_result;

		std::string quant_file_path_out = mrc_file_path_in + std::to_string(i + 1) + "_quantification.txt";
		std::string cropped_file_path_out = mrc_file_path_in + std::to_string(i + 1) + "_chopped.mrc_";
		std::string scatter_file_path_out = mrc_file_path_in + std::to_string(i + 1) + "_scatter.txt";

		mrc helix_mrc = cylinderCutOut(entire_map, helix);
		helix_mrc.normalize();

		helix_mrc.write(cropped_file_path_out);

		std::ofstream score_out_file(quant_file_path_out);
		std::vector<cbl::real> threshold_score;

		size_t j = 0;

		for (float t = 0.75f; t <= 1.0001f; t += 0.025f)
		{
			threshold_score.push_back(applyLineData(helix_mrc, helix, t, data_result));
			score_out_file << t << "\t" << threshold_score[j] << std::endl;
			j++;
		}

		int num_threshold_samples = j;

		score_out_file << "\n\nUtilized scores\n";

		double final_variance = 0;
		std::sort(threshold_score.begin(), threshold_score.end());
		for (int j = 0; j < (num_threshold_samples); j++)
		{
			final_variance += threshold_score[j];
			score_out_file << threshold_score[j] << std::endl;
		}

		final_variance /= (num_threshold_samples / 2 + 1);

		score_out_file << "----------------------------------------\n" << "Final Score: " << final_variance;
		score_out_file.close();


		// Display scatter plot
		std::ofstream out_file(scatter_file_path_out);
		out_file << data_result << std::endl;
		out_file.close();

		std::string scatter_command = "display_scatter_plot.bat \"" + scatter_file_path_out + "\"";
		// system(scatter_command.c_str());


		// Append total helix result to summary file

		std::ofstream summary("summary.txt", std::ofstream::out | std::ofstream::app);
		summary << mrc_file_path_in << i+1 << ", " << final_variance << std::endl;
		summary.close();
	}

	system("pause");

	return 0;
}