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
		int count = 0;
		auto f = [&](cbl::real x) {sum += x; count += x > 0; };
		std::for_each(v.begin(), v.end(), f);
		sum /= count;
		return sum;
	};

	auto variance = [](std::vector<cbl::real> v, cbl::real avg)
	{
		cbl::real square_sum = 0;
		cbl::real count = 0;
		auto f = [&](cbl::real x) {square_sum += (x > 0) * pow(avg - x, 2); count += x > 0; };
		std::for_each(v.begin(), v.end(), f);
		square_sum /= count;
		if (count > 1)
		{
			return square_sum;
		}
		else
		{
			return (cbl::real) 0.0;
		}
	};

	line.calcError(map, avg);

	line.mergeError(1.0); // Angstroms

	std::ostringstream out;

	line.write(out, threshold_used);

	data += out.str();

	cbl::real alpha = 0.5;
	cbl::real beta = 0.5;

	cbl::real threshold_score = alpha * sqrt(variance(line.error, avg(line.error))) + beta * avg(line.error);


	int count_zero = 0;
	for (size_t i = 0; i < line.error.size(); i++)
	{
		if (line.error[i] == 0)
		{
			count_zero++;
		}
	}

	std::cout << "Percentage zero: " << (cbl::real) count_zero / (cbl::real) line.error.size() << std::endl;

	threshold_score *= ((cbl::real) 1.0 - (cbl::real) count_zero / (cbl::real) line.error.size());

	return threshold_score;
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

std::vector<pdb> runAxisComparisonForHelixGeneration(std::string pdb_path)
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

void cylinderFitting(mrc &map, pdb &structure, std::string pdb_path)
{
	//Fit a cylinder to the true structure and find error to the cylinder

	pdb cylinder;

	for (int i = 0; i < structure.size(); i++)
	{



	}

	//cylinder.emplace_back();

	pdb_path = pdb_path + "_autoCylinder.pdb";

	cylinder.write(pdb_path);
}
void halfPoints(pdb &structure, std::string pdb_path)
{
	//vectors for the increment or decrement of each axis
	std::vector<float> x, y, z;

	for (int i = 0; i < structure.size() - 1 ; i++)
	{
		x.push_back(structure[i].x - structure[i + 1].x);
		y.push_back(structure[i].y - structure[i + 1].y);
		z.push_back(structure[i].z - structure[i + 1].z);
	}

	int switch1, switch2, switch3;
	switch1 = 0;
	switch2 = 0;
	switch3 = 0;

	for (int i = 0; i < structure.size() - 1; i++)
	{
		if ((x[i] > 0 && x[i + 1] < 0) || (x[i] < 0 && x[i + 1] > 0))
			switch1++;
		if ((y[i] > 0 && y[i + 1] < 0) || (y[i] < 0 && y[i + 1] > 0))
			switch2++;
		if ((z[i] > 0 && z[i + 1] < 0) || (z[i] < 0 && z[i + 1] > 0))
			switch3++;
	}

	//the vector with the least amount of direction changes will be the main axis

	int values[3] = { switch1, switch2, switch3 };
	std::sort(values, values + 3);

	char line;
	
	if (values[0] == switch1)
		line = 'x';
	if (values[0] == switch2)
		line = 'y';
	if (values[0] == switch3)
		line = 'z';

	//If a single helix folds in half, list points inbetween the two segments

	if (values[0] > 0)
	{
		pdb halfPoints;
		int l_shift, r_shift;
		cbl::real x1, y1, z1;

		//write any mid-points to file

		if (line == 'x')
		{
			for (int i = 0; i < x.size() - 1; i++)
			{
				if ((x[i] > 0 && x[i + 1] < 0) || (x[i] < 0 && x[i + 1] > 0))
				{
					l_shift = i + 1;
					r_shift = i + 1;
					
					while (l_shift > -1 && r_shift < x.size() + 1)
					{
						l_shift -= 1;
						r_shift += 1;

						x1 = (structure[l_shift].x + structure[r_shift].x) / 2;
						y1 = (structure[l_shift].y + structure[r_shift].y) / 2;
						z1 = (structure[l_shift].z + structure[r_shift].z) / 2;

						halfPoints.emplace_back(x1, y1, z1);
					}
				}
			}
		}

		if (line == 'y')
		{
			for (int i = 0; i < y.size() - 1; i++)
			{
				if ((y[i] > 0 && y[i + 1] < 0) || (y[i] < 0 && y[i + 1] > 0))
				{
					l_shift = i + 1;
					r_shift = i + 1;

					while (l_shift > -1 && r_shift < y.size() + 1)
					{
						l_shift -= 1;
						r_shift += 1;

						x1 = (structure[l_shift].x + structure[r_shift].x) / 2;
						y1 = (structure[l_shift].y + structure[r_shift].y) / 2;
						z1 = (structure[l_shift].z + structure[r_shift].z) / 2;

						halfPoints.emplace_back(x1, y1, z1);
					}
				}
			}
		}

		if (line == 'z')
		{
			for (int i = 0; i < z.size() - 1; i++)
			{
				if ((z[i] > 0 && z[i + 1] < 0) || (z[i] < 0 && z[i + 1] > 0))
				{
					l_shift = i + 1;
					r_shift = i + 1;

					while (l_shift > -1 && r_shift < z.size() + 1)
					{
						l_shift -= 1;
						r_shift += 1;

						x1 = (structure[l_shift].x + structure[r_shift].x) / 2;
						y1 = (structure[l_shift].y + structure[r_shift].y) / 2;
						z1 = (structure[l_shift].z + structure[r_shift].z) / 2;

						halfPoints.emplace_back(x1, y1, z1);
					}
				}
			}
		}

		//pdb_path = pdb_path + "_halfPoints.pdb";

		halfPoints.write("halfPoints.pdb");
	}
}

int main(int argc, char* argv[])
{
	int num_args = 2;
	assert(argc == num_args + 1 && "Incorrect num arguments, 2 inputs: mrc, pdb");

	std::string mrc_file_path_in = argv[1];
	std::string pdb_file_path_in = argv[2];
	
	std::vector<pdb> helices = runAxisComparisonForHelixGeneration(pdb_file_path_in);

	mrc entire_map(mrc_file_path_in);

	// Build stripped mrc file_name for creating out filenames
	size_t period_pos = mrc_file_path_in.find_last_of('.');
	mrc_file_path_in.resize(period_pos);
	mrc_file_path_in += "_helix";

	for (size_t i = 0; i < helices.size(); i++)
	{
		pdb& helix = helices[i];

		if (helix.size() > 1)
		{
			std::string data_result;

			std::string quant_file_path_out = mrc_file_path_in + std::to_string(i + 1) + "_quantification.txt";
			std::string cropped_file_path_out = mrc_file_path_in + std::to_string(i + 1) + "_chopped.mrc_";
			std::string scatter_file_path_out = mrc_file_path_in + std::to_string(i + 1) + "_scatter.txt";

			mrc helix_mrc = cylinderCutOut(entire_map, helix);
			helix_mrc.normalize();

			helix_mrc.write(cropped_file_path_out);

			//cylinderFitting(helix_mrc, helix, pdb_file_path_in);

			halfPoints(helix, pdb_file_path_in);

			std::ofstream score_out_file(quant_file_path_out);
			std::vector<cbl::real> threshold_score;

			size_t j = 0;

			for (float t = 0.75f; t <= 1.501f; t += 0.025f)
			{
				threshold_score.push_back(applyLineData(helix_mrc, helix, t, data_result));
				score_out_file << t << "\t" << threshold_score[j] << std::endl;
				j++;
			}

			int num_threshold_samples = j;

			score_out_file << "\n\nUtilized scores\n";

			int count = 0;
			double avg_variance = 0;
			std::sort(threshold_score.begin(), threshold_score.end());
			for (int j = 0; j < (num_threshold_samples); j++)
			{
				avg_variance += threshold_score[j];
				if (threshold_score[j] != 0)
				{
					count++;
				}
				score_out_file << threshold_score[j] << std::endl;
			}

			avg_variance /= count;

			if (avg_variance != 0)
			{
				score_out_file << "----------------------------------------\n" << "Final Score: " << avg_variance;
				score_out_file.close();


				// Save scatter plot image
				std::ofstream out_file(scatter_file_path_out);
				out_file << data_result << std::endl;
				out_file.close();

				std::string scatter_command = "generate_plot_as_image.bat \"" + scatter_file_path_out + "\"";
				system(scatter_command.c_str());

				// Append total helix result to summary file

				std::ofstream summary("summary.txt", std::ofstream::out | std::ofstream::app);
				summary << mrc_file_path_in << i+1 << ", " << avg_variance << std::endl;
				summary.close();
			}
		}
	}

	return 0;
}