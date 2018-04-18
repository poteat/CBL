#include <assert.h>
#include <string>
#include <functional>
#include <sstream>
#include <experimental/filesystem>

#include "axis.h"

using namespace cbl;

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