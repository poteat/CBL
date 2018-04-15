#include <assert.h>
#include <string>
#include <functional>
#include <sstream>
#include <experimental/filesystem>

#include "axis.h"

using namespace cbl;

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

		//test #2: detection using net distance

		std::vector<int> vertices = structure.halfPoints();

		//list points inbetween the two segments if vertices are detected

		if (vertices.size() > 0)
		{
			pdb halfPoints;
			int l_shift, r_shift;
			cbl::real x1, y1, z1;
			size_t l_bound, r_bound;


			if (vertices.size() > 1)
			{
				for (int i = 0; i < vertices.size() - 1; i++)
				{
					l_shift = vertices[i];
					r_shift = vertices[i];

					//long if/else statement to find the intervals between two maxima

					if (i == 0)
					{
						//The first interval for a helix with multiple vertices

						l_bound = 1;
						r_bound = vertices[i + 1]-1;
					}
					else if ((i != 0) && (i != vertices.size() - 1))
					{
						//The middle interval(s) for a helix with multiple vertices

						l_bound = vertices[i - 1]+1;
						r_bound = vertices[i + 1]-1;
					}
					else
					{
						//The final interval for a helix with multiple vertices

						l_bound = vertices[i - 1]+1;
						r_bound = structure.size()-2;
					}
				}
			}
			else
			{
				//The interval for a helix with one vertex (i.e. the entire interval for the helix)

				l_shift = vertices[0];
				r_shift = vertices[0];

				l_bound = 1;
				r_bound = structure.size()-2;
			}

			while(l_shift >= l_bound && r_shift < r_bound)
			{
				//calculates points inbetween the two sides

				l_shift -= 1;
				r_shift += 1;

				x1 = (structure[l_shift].x + structure[r_shift].x) / 2;
				y1 = (structure[l_shift].y + structure[r_shift].y) / 2;
				z1 = (structure[l_shift].z + structure[r_shift].z) / 2;

				halfPoints.emplace_back(x1, y1, z1);
			}
			
				
			//double check just in case halfpoints weren't generated

			if (halfPoints.size() > 0)
			{
				//output results

				std::string halfPointsFileName = pdb_file_path_in;
				halfPointsFileName = halfPointsFileName + "/" + halfPointsFileName + "_helix"
					+ std::to_string(helixFileNum) + "_halfPoints.pdb";

				halfPoints.write(halfPointsFileName);

				std::cout << "\nHalf Points written to file: " << halfPointsFileName << std::endl << std::endl;
			}
			else
			{
				std::cout << "\nHalf points generation error!" << std::endl;
			}
		}

		//No notable bends in the helix were detected

		else
		{
			std::cout << "\nNo fold detected!\n" << std::endl;
		}

	}

	std::cout << std::endl;

	system("PAUSE");

	return 0;
}