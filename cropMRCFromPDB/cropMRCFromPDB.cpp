#include <assert.h>
#include <string>
#include <functional>
#include <fstream>

#include "off.h"

using namespace cbl;

// CMD PARAMETERS: "input.mrc" "input.pdb" dist "out_near.mrc" "out_far.mrc"
int main(int argc, char *argv[])
{
	// Confirm we have correct number of arguments
	int num_arguments = 5;
	assert(argc == num_arguments + 1 && "Usage: input.mrc input.pdb dist out_near.mrc out_far.mrc");

	// Read in each input
	mrc in_mrc(argv[1]);
	pdb in_pdb(argv[2]);

	cbl::real cropping_dist = (cbl::real) std::stod(argv[3]);

	std::string out_near_mrc_filename = argv[4];
	std::string out_far_mrc_filename = argv[5];

	// Perform cropping
	mrc near, far;
	std::tie(near, far) = in_mrc.crop(in_pdb, cropping_dist);

	// Write output to files
	near.write(out_near_mrc_filename);
	far.write(out_far_mrc_filename);

	return 0;
}