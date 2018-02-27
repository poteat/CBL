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

				if (elements.size() == 9 && elements[0] == "ATOM")
				{
					std::string atom = elements[2];
					real x = (real) std::stod(elements[6]);
					real y = (real) std::stod(elements[7]);
					real z = (real) std::stod(elements[8]);

					points.emplace_back(x, y, z, atom);
				}
			}

			std::cout << "Read " << points.size() << " atoms from " << file_name << std::endl;
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

		const point& operator[](size_t i)
		{
			return points[i];
		}

		size_t size()
		{
			return points.size();
		}
	};
}
