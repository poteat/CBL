#include <assert.h>
#include <string>
#include <functional>
#include <sstream>
#include <experimental/filesystem>

#include "axis.h"
#include "mrc.h"

using namespace cbl;

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
			std::string cropped_file_path_out = mrc_file_path_in + std::to_string(i + 1) + "_chopped.mrc";
			std::string scatter_file_path_out = mrc_file_path_in + std::to_string(i + 1) + "_scatter.txt";

			mrc helix_mrc = cylinderCutOut(entire_map, helix);
			helix_mrc.normalize();

			helix_mrc.write(cropped_file_path_out);

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