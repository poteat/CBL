#include <assert.h>
#include <string>
#include <functional>
#include <sstream>
#include <experimental/filesystem>

#include "axis.h"

using namespace cbl;

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

int main(int argc, char* argv[])
{
	int num_args = 1;
	assert(argc == num_args + 1 && "Incorrect num arguments, 2 inputs: mrc, pdb");

	std::string pdb_file_path_in = argv[1];
	std::vector<pdb> helices = runAxisComparisonForHelixGeneration(pdb_file_path_in);

	for (size_t i = 0; i < helices.size(); i++)
	{
		pdb& helix = helices[i];
		pdb &structure = helix;

		//vectors for the increment or decrement of each axis
		std::vector<float> x, y, z;

		for (int i = 0; i < structure.size() - 1; i++)
		{
			x.push_back(structure[i].x - structure[i + 1].x);
			y.push_back(structure[i].y - structure[i + 1].y);
			z.push_back(structure[i].z - structure[i + 1].z);
		}

		int switch1, switch2, switch3;
		switch1 = 0;
		switch2 = 0;
		switch3 = 0;

		for (int i = 0; i < x.size() - 1; i++)
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

		std::cout << "Selected Line: " << line << std::endl;

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

			std::cout << "Half Points detected" << std::endl << std::endl;

			std::string halfPoints_file_path_out = pdb_file_path_in + std::to_string(i + 1) + "_halfPoints.pdb";

			halfPoints.write(halfPoints_file_path_out);

		}

	}

    return 0;
}

