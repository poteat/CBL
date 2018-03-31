#include <assert.h>
#include <string>
#include <functional>

#include "off.h"
#include "plane.h"

using namespace cbl;

int main(int argc, char *argv[])
{
	// Confirm we have correct number of arguments (.off, X.mrc)
	int num_args = 2;
	assert(argc == num_args + 1 && "Incorrect num arguments: .off, X.mrc");

	// Save base path for writing files during body
	std::string base_path = argv[2];
	size_t period_pos = base_path.find_last_of('.');
	base_path.resize(period_pos);


	off skeleton_mesh(argv[1]);
	mrc original_map(argv[2]);


	// Perform classification on original MRC file (CLAS)
	mrc helix, sheet;
	std::tie(helix, sheet) = skeleton_mesh.classify(original_map);

	sheet.applyDeviationThreshold(2.75);

	std::vector<mrc> sheet_list = sheet.cluster(70); // Exclude clusters less than 70 voxels

													 // Write individual clusters to disk
	for (size_t i = 0; i < sheet_list.size(); i++)
	{
		sheet_list[i].write(base_path + "_sheet" + std::to_string(i + 1) + ".mrc");
	}


	// Create n planes based on filtered cluster vector

	std::vector<plane> planes;

	for (size_t i = 0; i < sheet_list.size(); i++)
	{
		planes.emplace_back(sheet_list[i]);
	}



	// Write individual planes to disk

	for (size_t i = 0; i < planes.size(); i++)
	{
		planes[i].write(base_path + "_plane" + std::to_string(i + 1) + ".pdb");
	}




	pdb pdb_off_file = skeleton_mesh.convertToPDB(original_map);

	helix.write(base_path + "_helix.mrc");
	sheet.write(base_path + "_sheet.mrc");
	pdb_off_file.write(base_path + "_skeleton.pdb");

	system("pause");

	return 0;
}