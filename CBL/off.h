#include <tuple>
#include <limits>
#include <math.h>

#include "mrc.h"
#include "pdb.h"

namespace cbl
{
	class off
	{
	public:
		enum sse_type { Unclassified, Helix, Sheet };

		struct off_point
		{
			real x, y, z;
			int num_helix = 0;
			int num_sheet = 0;
			sse_type type = Unclassified;
			off_point(real x, real y, real z) : x(x), y(y), z(z) {};
			operator point() { return point(x, y, z); }
		};

		int num_vertices, num_faces, num_edges;
		std::vector<off_point> points;

		off(std::string file_name)
		{
			std::ifstream file(file_name);
			assert(file && "Could not open .off file");

			std::string type;

			file >> type;
			assert(type == "OFF" && "File format specifier not correct for .off files");

			file >> num_vertices >> num_faces >> num_edges;

			for (int i = 0; i < num_vertices; i++)
			{
				real x, y, z;
				file >> x >> y >> z;
				points.emplace_back(x, y, z);
			}

			for (int i = 0; i < num_faces; i++)
			{
				int num_points_in_face;
				file >> num_points_in_face;

				for (int i = 0; i < num_points_in_face; i++)
				{
					int index;
					file >> index;

					if (num_points_in_face == 4)
					{
						points[index].num_helix++;
					}
					else if (num_points_in_face == 3)
					{
						points[index].num_sheet++;
					}
				}
			}

			// Classification heuristic

			for (unsigned int i = 0; i < points.size(); i++)
			{
				if (points[i].num_sheet == 0)
				{
					points[i].type = Helix;
				}
				else
				{
					points[i].type = Sheet;
				}
			}

			file.close();
		}

		void print()
		{
			for (unsigned int i = 0; i < points.size(); i++)
			{
				std::string type;

				if (points[i].type == Helix)
				{
					type = "Helix";
				}
				else if (points[i].type == Sheet)
				{
					type = "Sheet";
				}

				std::cout << points[i].x << " " << points[i].y << " " << points[i].z << " "
					<< type << std::endl;
			}
		}

		std::tuple<mrc, mrc> classify(mrc& original_map)
		{
			mrc helix_map(original_map);
			mrc sheet_map(original_map);

			// Identify closest skeleton points in original map, and use these results
			// to mask the copied maps.

			real scale = original_map.scale;
			real x_origin = original_map.header.xorigin;
			real y_origin = original_map.header.yorigin;
			real z_origin = original_map.header.zorigin;

			real avg_x = (real)original_map.header.nx / (real)2;
			real avg_y = (real)original_map.header.ny / (real)2;
			real avg_z = (real)original_map.header.nz / (real)2;

			int num_sheet = 0;
			int num_helix = 0;

			for (size_t i = 0; i < original_map.header.nx; i++)
			{
				for (size_t j = 0; j < original_map.header.ny; j++)
				{
					for (size_t k = 0; k < original_map.header.nz; k++)
					{
						real density = original_map.map(i, j, k);

						if (density > 0)
						{
							// Find closest point from skeleton

							// First we convert int-coordinates to real-coordinates.

							real x = (((real)i - avg_x) + x_origin) * scale;
							real y = (((real)j - avg_y) + y_origin) * scale;
							real z = (((real)k - avg_z) + z_origin) * scale;

							// Loop through skeleton points, record minimum distance

							real min_dist = std::numeric_limits<real>::infinity();
							int min_index = 0;

							for (unsigned int m = 0; m < points.size(); m++)
							{
								real dist = pow(points[m].x - x, 2) +
									pow(points[m].y - y, 2) +
									pow(points[m].z - z, 2);

								if (dist < min_dist)
								{
									min_dist = dist;
									min_index = m;
								}
							}

							auto type = points[min_index].type;

							if (type == Helix)
							{
								num_helix++;
								sheet_map.map(i, j, k) = 0;
							}
							else if (type == Sheet)
							{
								num_sheet++;
								helix_map.map(i, j, k) = 0;
							}
						}
					}
				}
			}

			std::cout << "Classified " << num_helix << " helix points and " << num_sheet << " sheet points" << std::endl;

			return std::make_tuple(helix_map, sheet_map);
		}

		pdb convertToPDB(mrc mrc_reference)
		{
			pdb result;

			real scale = mrc_reference.scale;
			real x_origin = mrc_reference.header.xorigin;
			real y_origin = mrc_reference.header.yorigin;
			real z_origin = mrc_reference.header.zorigin;

			real avg_x = (real)mrc_reference.header.nx / (real)2;
			real avg_y = (real)mrc_reference.header.ny / (real)2;
			real avg_z = (real)mrc_reference.header.nz / (real)2;

			for (size_t i = 0; i < points.size(); i++)
			{
				std::string atom = "H";

				if (points[i].type == Sheet)
				{
					atom = "O"; // b-sheets are red
				}
				else if (points[i].type == Helix)
				{
					atom = "N"; // helices are blue
				}
				else
				{
					atom = "Li"; // uncategorized skeleton points are purple
				}

				real x = points[i].x + x_origin + (avg_x - x_origin) * scale;
				real y = points[i].y + y_origin + (avg_y - y_origin) * scale;
				real z = points[i].z + z_origin + (avg_z - z_origin) * scale;

				result.emplace_back(x, y, z, atom);
			}

			return result;
		}
	};
}