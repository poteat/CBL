#include <assert.h>
#include <string>
#include <functional>

#include "mrc.h"

using namespace cbl;

int main(int argc, char *argv[])
{
	// Confirm we have correct number of arguments (The .mrc file)
	assert(argc == 2 && "Program must take in one argument (mrc file path)");
	std::string mrc_file_path(argv[1]);

	auto flip_density = [](cbl::real d) {return -d; };
	auto remove_negative = [](cbl::real d) {return d > 0 ? d : 0; };

	// Perform classification on original MRC file
	mrc mrc_file(mrc_file_path);

	mrc_file.apply(flip_density);
	mrc_file.apply(remove_negative);

	mrc_file.write(mrc_file_path);

	return 0;
}