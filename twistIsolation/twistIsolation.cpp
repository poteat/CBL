#include <assert.h>
#include <string>
#include <functional>
#include <sstream>
#include <experimental/filesystem>

#include "axis.h"
#include "core.h"

//This program generates points inbetween each structure in a pdb, to create a boundary for later reference

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

	std::vector<cbl::pdb> structures;

	for (auto &p : std::experimental::filesystem::directory_iterator(output_path))
	{
		auto s = p.path().filename().string();

		// If begins with "trueHelix"
		if (s.find("trueHelix") == 0 || s.find("trueSheet") == 0)
		{
			structures.emplace_back(p.path().string());
		}
	}

	// Rename output directory to something we can refer to later

	path new_output_path = symbolic_pdb_path.parent_path().string() + "/"
		+ symbolic_pdb_path.stem().string();

	std::experimental::filesystem::remove_all(new_output_path);

	std::experimental::filesystem::rename(output_path, new_output_path);

	return structures;
}

int main(int argc, char* argv[])
{
	if (argc > 2)
	{
		assert("Incorrect num arguments, 2 inputs: pdb, (optional: bool value on whether to run testing)");
		return 1;
	}

	bool runTesting;

	if (argc == 2)
	{
		runTesting = argv[2];
	}
	else
	{
		runTesting = false;
	}

	std::string pdb_file_path_in = argv[1];
	std::vector<pdb> allstructures = runAxisComparisonForHelixGeneration(pdb_file_path_in);
	size_t period_pos = pdb_file_path_in.find_last_of('.');
	pdb_file_path_in.resize(period_pos);
	std::vector<std::vector<size_t>> vertices;
	vertices.resize(allstructures.size());

	size_t errors = 0;
	std::vector<size_t> numbers;

	//check for any errors with helixes, if there is an error, then rename files to reflect it
	for (size_t check = 0; check < allstructures.size(); check++)
	{
		if (allstructures[check].size() == 1)
		{
			std::cout << "Structure number " << check + 1 << " is incorrectly detected by Axis Comparison. Deleting... " << std::endl;

			errors++;

			std::vector<pdb>::iterator it = allstructures.begin();
			std::advance(it, check);
			allstructures.erase(it);

			check--;
		}
		else
		{
			numbers.push_back(check + errors);
		}

	}

	std::ofstream pointList;
	std::string pointFile = pdb_file_path_in + "_AllPoints.txt";
	pointList.open(pointFile);
	std::vector<pdb> finalPoints;
	finalPoints.resize(allstructures.size());

	if (allstructures.size() == 1)
	{
		std::cout << "\nNo other structures in PDB detected..." << std::endl;

		system("PAUSE");

		return 0;
	}

	for (size_t n = 0; n < allstructures.size(); n++)
	{
		pdb& helix = allstructures[n];
		pdb &structure = helix;

		vertices[n].resize(2);

		//test #2: detection using net distance

		std::vector<cbl::real> distances;

		cbl::point p1 = structure[0];
		cbl::point p2 = structure[structure.size() - 1];

		//find discrepencies in distance from generated line
		for (size_t i = 0; i < structure.size(); i++)
		{
			cbl::point p3 = structure[i];

			//makes a virtual line from p1 to p2 and finds the distance to the line from p3
			cbl::real distance = distToLine(p1, p2, p3);

			distances.push_back(distance);
		}

		vertices[n].push_back(0);

		for (size_t i = 1; i < distances.size() - 1; i++)
		{
			//try to add index location to a vector if the discrepency is notable

			if (distances[i] > 5)
			{
				//checks if a point is a local maxima

				if ((distances[i] > distances[i - 1]) && (distances[i] > distances[i + 1]))
				{
					vertices[n].resize(vertices.size());
					vertices[n].push_back(i);
				}

				//if two points tie as a local maxima, the first point is chosen and the second is skipped

				if ((distances[i] > distances[i - 1]) && (distances[i] == distances[i + 1]))
				{
					if ((size_t)(i + 2) < distances.size())
					{
						if (distances[i + 1] > distances[i + 2])
						{
							vertices[n].resize(vertices.size());
							vertices[n].push_back(i);
							i++;
						}
					}
				}
			}
		}

		vertices[n].push_back(allstructures.size() - 1);

		//vertices vector of vector now represents different line segments to check for nearby structures
	}

	//The previous code to this point found the different line segments that make up the helixes and strands passed to the program
	//The majority of helixes and strands will only have one line segment (From the first point to the last point) unless they fold in on themselves

	//The number of line segments that will be used for distance calculations / loops
	size_t segments = 0;

	for (size_t i = 0; i < vertices.size(); i++)
	{
		for (size_t j = 0; j < vertices[i].size(); j++)
		{
			segments++;
		}
		segments--;
	}

	std::vector<cbl::real> distanceSort;
	distanceSort.resize(segments);

	std::vector<cbl::point> midpoints;

	//store the midpoint of each line
	for (size_t j = 0; j < allstructures.size(); j++)
	{
		pdb& all = allstructures[j];
		pdb &currentstructure = all;		

		size_t tempindex = (currentstructure.size() - 1) / 2;
		cbl::point p1 = (currentstructure[tempindex]);
		midpoints.push_back(p1);
	}

	std::vector<std::vector<cbl::real>> distances;
	distances.resize(midpoints.size());

	//store the distance from each midpoint to all other midpoints
	for (size_t i = 0; i < midpoints.size(); i++)
	{
		for (size_t j = 0; j < midpoints.size(); j++)
		{
			if (i != j)
			{
				cbl::point p1 = midpoints[i];
				cbl::point p2 = midpoints[j];
				distances[i].push_back(p1.dist(p2));
			}
			if (i == j)
			{
				distances[i].push_back(NULL);
			}
		}
	}

	for (size_t i = 0; i < midpoints.size(); i++)
	{
		size_t strand = NULL;
		size_t index2 = NULL;
		size_t index1 = NULL;
		cbl::real minDist = Infinity;

		//Which strand is closest to line 1?
		for (size_t j = 0; j < midpoints.size(); j++)
		{
			if (distances[j][i] < minDist && i != j)
			{
				strand = j;
				minDist = distances[j][i];
			}
		}

		minDist = INFINITY;
		cbl::real minDist2 = INFINITY;
		pdb& all = allstructures[strand];
		pdb &currentstructure = all;

		//Which point on line 2 is closest to line 1?
		for (size_t k = 0; k < currentstructure.size(); k++)
		{
			cbl::point p1 = currentstructure[k];
			cbl::point p2 = midpoints[i];
			if (p1.dist(p2) < minDist)
			{
				index2 = k;
				minDist = p1.dist(p2);
			}
		}

		//Which point on line 1 is closest to line 2?
		pdb& all1 = allstructures[i];
		pdb &currentstructure1 = all1;
		for (size_t k = 0; k < currentstructure1.size(); k++)
		{
			cbl::point p1 = currentstructure1[k];
			cbl::point p2 = midpoints[strand];
			if (p1.dist(p2) < minDist)
			{
				index1 = k;
				minDist2 = p1.dist(p2);
			}
		}

		//After finding the two closest points to each other, calculate midpoint between them
		cbl::real x1, y1, z1;
		pdb& all2 = allstructures[i];
		pdb &structure1 = all2;
		pdb& all3 = allstructures[strand];
		pdb &structure2 = all3;

		x1 = (structure1[index1].x + structure2[index2].x) / 2;
		y1 = (structure1[index1].y + structure2[index2].y) / 2;
		z1 = (structure1[index1].z + structure2[index2].z) / 2;

		//saves the points based on the structure that they belong to
		pointList << x1 << "\t" << y1 << "\t" << z1 << std::endl;
		finalPoints[i].emplace_back(x1, y1, z1);
		finalPoints[strand].emplace_back(x1, y1, z1);
	}

	pointList.close();

	//Finally, output the generated half-points based on which structure they belong to

	size_t current = 0;

	for (size_t n = 0; n < allstructures.size(); n++)
	{
		if (finalPoints[n].size() > 0)
		{
			//safety net for file numbering, if vector is out of range, throw numbers out the window

			if (numbers.size() != allstructures.size())
			{
				current = numbers[n] + 1;
			}
			else
			{
				current = n;
			}

			//output results

			std::string finalPointsFileName = pdb_file_path_in;
			finalPointsFileName = finalPointsFileName + "/" + finalPointsFileName + "_structure"
				+ std::to_string(current) + "_isolationPoints.pdb";

			finalPoints[n].write(finalPointsFileName);

			std::cout << "\nIsolation Points written to file: " << finalPointsFileName << std::endl << std::endl;
		}

		else
		{
			std::cout << "\nNo other structures in PDB detected..." << std::endl;
		}
	}

	std::cout << "\nAll Isolation points written to file: " << pointFile << std::endl;

	if (runTesting = true)
	{
		//experimental methods for creating 3-d shapes to fit helices






	}

	std::cout << std::endl;

	system("PAUSE");

	return 0;
}