#include <assert.h>
#include <string>
#include <functional>

#include "off.h"

using namespace cbl;

int main(int argc, char *argv[])
{
	// Confirm we have correct number of arguments (The .off file path)
	assert(argc == 2 && "Incorrect number of arguments");
	std::string off_file_path(argv[1]);

	// Build mrc file path (Replace last .off with .mrc)
	std::string mrc_file_path = off_file_path;
	size_t period_pos = mrc_file_path.find_last_of('.');
	mrc_file_path.resize(period_pos);
	mrc_file_path += ".mrc";

	// Read in MRC file and OFF file
	mrc original_map(mrc_file_path);
	off skeleton_mesh(off_file_path);

	std::cout << off_file_path << ": ";

	// Perform classification on original MRC file
	mrc helix, sheet;
	std::tie(helix, sheet) = skeleton_mesh.classify(std::ref(original_map));

	pdb pdb_off_file = skeleton_mesh.convertToPDB(std::ref(original_map));

	// Write the segmented MRC files to disk with name appended
	period_pos = mrc_file_path.find_last_of('.');
	mrc_file_path.resize(period_pos);
	helix.write(mrc_file_path + "_helix.mrc");
	sheet.write(mrc_file_path + "_sheet.mrc");
	pdb_off_file.write(mrc_file_path + "_skeleton.pdb");

	return 0;
}