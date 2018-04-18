#include <assert.h>
#include <string>
#include <functional>
#include <sstream>
#include <experimental/filesystem>

#include "axis.h"

//For every helix in a pdb file, detects nearby structures and outputs a boundary around the helix to isolate it from the rest of the pdb

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
	assert(argc == num_args + 1 && "Incorrect num arguments, 1 input: pdb");

	std::string pdb_file_path_in = argv[1];
	std::vector<pdb> helices = runAxisComparisonForHelixGeneration(pdb_file_path_in);
	size_t period_pos = pdb_file_path_in.find_last_of('.');
	pdb_file_path_in.resize(period_pos);

	int helixFileNum = 1;

	for (size_t i = 0; i < helices.size(); i++)
	{
		pdb& helix = helices[i];
		pdb &structure = helix;

		

	}

	system("PAUSE");

	return 0;
}