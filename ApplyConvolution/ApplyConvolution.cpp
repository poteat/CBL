#include <assert.h>
#include <string>
#include <functional>

#include "mrc.h"

using namespace cbl;

int main(int argc, char *argv[])
{
	// Confirm we have correct number of arguments (The two files)
	assert(argc == 2 && "Program must take in one argument, .mrc file");
	std::string mrc_file(argv[1]);

	mrc map(mrc_file);
	cbl::cube<cbl::real> kernel("kernel.txt");

	auto remove_negative = [](cbl::real d) {return d > 0 ? d : 0; };

	map.apply(kernel);
	map.apply(remove_negative);

	map.trim(3);

	size_t period_pos = mrc_file.find_last_of('.');
	mrc_file.resize(period_pos);

	map.write(mrc_file + "_convolved.mrc");

	return 0;
}