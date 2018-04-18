#include <assert.h>
#include <string>
#include <functional>
#include <sstream>
#include <experimental/filesystem>

#include "axis.h"

//Unfinished code to be worked on after twistIsolation is finished
	//Detects the shape of an mrc file surrounding individual helixes

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

cbl::pdb innerCylinder(pdb &structure)
{
	//Fit a cylinder to the true structure

	pdb cylinder;
	cbl::real x, y, z;

	for (int i = 0; i < structure.size(); i++)
	{
		x = structure[i].x;
		y = structure[i].y;
		z = structure[i].z;


	}

	return cylinder;
}

cbl::pdb outerCylinder(mrc &map, pdb &structure, float threshold)
{
	//Fit a cylinder to the border of the mrc map

	pdb cylinder;
	cbl::real x, y, z;
	cbl::real minx, maxx;
	cbl::real miny, maxy;
	cbl::real minz, maxz;
	minx = INFINITY;
	miny = INFINITY;
	minz = INFINITY;
	maxx = 0;
	maxy = 0;
	maxz = 0;

	for (int i = 0; i < structure.size(); i++)
	{
		x = structure[i].x;
		y = structure[i].y;
		z = structure[i].z;

		//if there is density at this x,y,z

		if (x > maxx)
			maxx = x;
		if (y > maxy)
			maxy = y;
		if (z > maxz)
			maxz = z;
		if (x < maxx)
			minx = x;
		if (x < maxx)
			miny = y;
		if (x < maxx)
			minz = z;
	}

	cbl::real xdenominator, ydenominator, zdenominator;
	xdenominator = .5*((maxx - minx)*(maxx - minx));
	ydenominator = .5*((maxy - miny)*(maxy - miny));
	zdenominator = .5*((maxz - minz)*(maxz - minz));

	for (int i = 0; i < structure.size(); i++)
	{
		x = structure[i].x;
		y = structure[i].y;
		z = structure[i].z;

		//cylinder.emplace_back((x*x) / xdenominator, (y*y) / ydenominator, (z*z) / zdenominator);
		//not correct currently, we only want the points where these three values added together equal 2.5 
		//On the bright side, using this formula we can get the distance from the surface directly
		//i.e. if less than 2.5, the difference is the distance from the surface
	}

	return cylinder;
}

int main(int argc, char* argv[])
{
	if (argc == 1 || argc > 3)
	{
		assert("Incorrect num arguments, 2 or 3 inputs: mrc, pdb, (optional: threshold value)");
		return 1;
	}

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
			std::string inner_cylinder_file_path_out = mrc_file_path_in + std::to_string(i + 1) + "_innerCylinder.pdb";
			std::string outer_cylinder_file_path_out = mrc_file_path_in + std::to_string(i + 1) + "_outerCylinder.pdb";
			std::string cropped_file_path_out = mrc_file_path_in + std::to_string(i + 1) + "_chopped.mrc";

			mrc helix_mrc = cylinderCutOut(entire_map, helix);
			helix_mrc.normalize();

			helix_mrc.write(cropped_file_path_out);

			//this is where the magic happens

			float threshold = 0;

			if (argc == 3)
			{
				threshold = std::stof(argv[3]);
			}
			else
			{
				helix_mrc.applyDeviationThreshold(2);
			}

			pdb inner_cylinder = innerCylinder(helix);
			pdb outer_cylinder = outerCylinder(helix_mrc, helix, threshold);
			inner_cylinder.write(inner_cylinder_file_path_out);
			outer_cylinder.write(outer_cylinder_file_path_out);

			
		}
		else if (helix.size() == 1)
		{
			float threshold = 0;

			if (argc == 3)
			{
				threshold = std::stof(argv[3]);
			}
			else
			{
				entire_map.applyDeviationThreshold(2);
			}

			pdb inner_cylinder = innerCylinder(helix);
			pdb outer_cylinder = outerCylinder(entire_map, helix, threshold);

			std::string inner_cylinder_file_path_out = mrc_file_path_in + "_innerCylinder.pdb";
			inner_cylinder.write(inner_cylinder_file_path_out);

			std::string outer_cylinder_file_path_out = mrc_file_path_in + "_outerCylinder.pdb";
			inner_cylinder.write(outer_cylinder_file_path_out);
		}
	}

	return 0;
}