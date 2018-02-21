#include <assert.h>
#include <string>
#include <functional>
#include <fstream>

#include "off.h"

using namespace cbl;

int main(int argc, char *argv[])
{
	// Confirm we have correct number of arguments
	int num_arguments = 2;
	assert(argc == num_arguments + 1 && "Needs two arguments (in_list.txt, out.pdb");

	// Verify we could open pdb_list for reading
	std::ifstream pdb_list(argv[1]);
	assert(pdb_list && "Could not open pdb_list.txt file (1st argument)");

	// Open output filename as string
	std::string out_filename(argv[2]);

	std::vector<pdb> input_pdb;
	std::string filename;

	// Read in all pdb files in list
	while (pdb_list >> filename)
	{
		input_pdb.emplace_back(filename);
	}

	pdb merged(input_pdb);

	merged.write(out_filename);

	return 0;
}