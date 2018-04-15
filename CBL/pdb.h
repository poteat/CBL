#pragma once

#include <fstream>
#include <vector>
#include <string>
#include <assert.h>
#include <iomanip>
#include <iostream>

#include "core.h"

namespace cbl
{
	class pdb
	{
	private:
		struct pdb_point
		{
			real x = 0, y = 0, z = 0;
			std::string atom;
			pdb_point(real x, real y, real z, std::string atom) : x(x), y(y), z(z), atom(atom) {};
			operator point() { return point(x, y, z); }
		};

		std::vector<pdb_point> points;

	public:
		pdb(){};

		pdb(std::string file_name)
		{
			std::ifstream file(file_name);
			assert(file && "Could not open .pdb file for reading");
			std::string line;
			std::vector<std::string> elements;

			while (std::getline(file, line))
			{
				line = removeConsecutive(line, ' ');
				elements = explode(line, ' ');

				if (elements[0] == "ATOM")
				{
					std::string atom = elements[2];
					real x = (real) std::stod(elements[6]);
					real y = (real) std::stod(elements[7]);
					real z = (real) std::stod(elements[8]);

					points.emplace_back(x, y, z, atom);
				}
			}
		};

		pdb(std::vector<pdb> inputs)
		{
			for (size_t i = 0; i < inputs.size(); i++)
			{
				pdb p = inputs[i];

				for (size_t j = 0; j < p.points.size(); j++)
				{
					points.push_back(p.points[j]);
				}
			}
		};

		void emplace_back(real x, real y, real z, std::string atom = "H")
		{
			points.emplace_back(x, y, z, atom);
		};

		void emplace_back(cbl::point p, std::string atom = "H")
		{
			points.emplace_back(p.x, p.y, p.z, atom);
		};

		void write(std::string file_name)
		{
			std::ofstream file(file_name);
			assert(file && "Could not open .pdb file for writing");
			assert(points.size() && "Attempted to write an empty .pdb file");

			for (size_t i = 0; i < points.size(); i++)
			{
				real x = points[i].x;
				real y = points[i].y;
				real z = points[i].z;
				std::string atom = points[i].atom;

				real epsilon = (real) 0.0001;

				x = abs(x) < epsilon ? 0 : x;
				y = abs(y) < epsilon ? 0 : y;
				z = abs(z) < epsilon ? 0 : z;

				std::string space;
				if (atom.size() == 1)
				{
					space = "  ";
				}
				else
				{
					space = " ";
				}

				file << "ATOM " << std::setw(6) << std::right << i << \
					"  " << atom << space << " HOH A   1    " << std::setprecision(5) \
					<< std::setw(7) << std::right << x << " " \
					<< std::setw(7) << std::right << y << " " \
					<< std::setw(7) << std::right << z << std::endl;
			}
		};

		pdb_point& operator[](size_t i)
		{
			return points[i];
		}

		size_t size()
		{
			return points.size();
		}

		std::vector<int> halfPoints()
		{
			std::vector<int> vertices;
			std::vector<cbl::real> distances;

			cbl::point p1 = points[0];
			cbl::point p2 = points[points.size() - 1];

			//find discrepencies in distance from generated line
			for (int i = 0; i < points.size(); i++)
			{
				cbl::point p3 = points[i];

				//makes a virtual line from p1 to p2 and finds the distance to the line from p3
				cbl::real distance = distToLine(p1, p2, p3);

				distances.push_back(distance);
			}

			std::cout << std::endl;

			for (int i = 1; i < distances.size() - 1; i++)
			{
				//try to add index location to a vector if the discrepency is notable

				if (distances[i] > 5)
				{
					//checks if a point is a local maxima

					if ((distances[i] > distances[i - 1]) && (distances[i] > distances[i + 1]))
					{
						std::cout << "Detected vertex at index: " << i << std::endl;

						vertices.push_back(i);
					}

					//if two points tie as a local maxima, the first point is chosen and the second is skipped

					if ((distances[i] > distances[i - 1]) && (distances[i] == distances[i + 1]))
					{
						if ((i + 2) < distances.size())
						{
							if (distances[i + 1] > distances[i + 2])
							{
								std::cout << "Detected vertex at index: " << i << std::endl;
								vertices.push_back(i);
								i++;
							}
						}
					}
				}
			}
			return vertices;
		}
	};
}
