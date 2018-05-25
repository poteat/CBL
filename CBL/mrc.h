#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <assert.h>
#include <tuple>
#include <stack>
#include <random>

#include "core.h"
#include "pdb.h"

#include "Eigen/Dense"

namespace cbl
{
	class mrc
	{
	public:
		struct header           // Header is exactly 1kb (1024 bytes)
		{
			unsigned int nx;             /* # of columns ( fastest changing in the map    */
			unsigned int ny;             /* # of rows                                     */
			unsigned int nz;             /* # of sections (slowest changing in the map    */

			int mode;           // 2 (reals)

			int nxstart;        // 0
			int nystart;        // 0
			int nzstart;        // 0

			int mx;             // Same as nx
			int my;             // Same as ny
			int mz;             // Same as nz

			real xlength;      // Width of entire block in angstroms
			real ylength;      // Height of entire block in angstroms
			real zlength;      // Depth of entire block in angstroms

			real alpha;        // 90.0
			real beta;         // 90.0
			real gamma;        // 90.0

			int mapc;           // 1
			int mapr;           // 2
			int maps;           // 3

			real amin;         /* minimum density value                         */
			real amax;         /* maximum density value                         */
			real amean;        /* mean density value                            */

			int ispg;           // 0
			int nsymbt;         // 0

			int extra[25];

			real xorigin;      // Used to convert between PDB and MRC space
			real yorigin;      // 
			real zorigin;      // 

			char map[4];        // String "MAP"

			int machineStamp;

			real rms;

			int nlabl;
			char label[10][80];
		} header;

		cube<real> map;
		real scale;

		mrc() {};

		mrc(std::string file_name)
		{
			std::ifstream file(file_name, std::ios::binary);
			assert(file && "Could not open .mrc file");

			file.read((char*)&header, 1024);

			// assert(header.map == "MAP" && ".mrc file had incorrect format specifier");

			this->scale = header.xlength / (real)header.nx;

			assert(header.nx > 0);
			assert(header.ny > 0);
			assert(header.nz > 0);

			map = cube<real>(header.nx, header.ny, header.nz);

			// Read in real densities sequentially.

			for (size_t k = 0; k < header.nz; k++)
			{
				for (size_t j = 0; j < header.ny; j++)
				{
					for (size_t i = 0; i < header.nx; i++)
					{
						file.read((char*)&map(i, j, k), 4);
					}
				}
			}
		};

		void write(std::string file_name)
		{
			assert(map.count() > 0 && "Attempted to write MRC file with no positive density");

			updateHeaderMinMax();

			std::ofstream file(file_name, std::ios::binary);

			file.write((char*)&header, 1024);

			for (size_t k = 0; k < header.nz; k++)
			{
				for (size_t j = 0; j < header.ny; j++)
				{
					for (size_t i = 0; i < header.nx; i++)
					{
						file.write((char*)&map(i, j, k), 4);
					}
				}
			}

			file.close();
		}

		void updateHeaderMinMax()
		{
			// Update min, mean, and max density values

			header.amin = map.min();
			header.amean = map.mean();
			header.amax = map.max();
		}

		void print()
		{
			for (size_t i = 0; i < header.nx; i++)
			{
				for (size_t j = 0; j < header.ny; j++)
				{
					for (size_t k = 0; k < header.nz; k++)
					{
						real d = map(i, j, k);
						real x = i * scale + header.xorigin;
						real y = j * scale + header.yorigin;
						real z = k * scale + header.zorigin;

						if (d > 0)
						{
							std::cout << x << " " << y << " " << z << std::endl;
						}
					}
				}
			}

			std::cout << "Origin: " << header.xorigin << " " << header.yorigin << " " << header.zorigin;
		}

		void normalize()
		{
			real min = map.min();
			real max = map.max();

			for (size_t i = 0; i < map.size(); i++)
			{
				if (map[i] > 0)
				{
					map[i] = (map[i] - min) / (max - min);
				}

			}
		}

		void applyDeviationThreshold(real multiplier)
		{
			real threshold = multiplier * map.standard_deviation();

			std::cout << "Applying threshold " << threshold << std::endl;

			applyThreshold(threshold);
		}

		// Given a mrc that has been thresholded, only include voxels which have an empty neighbor in their
		// 26-neighborhood
		void setHollowBoundary()
		{
			auto equal_zero = [](real x) {return x == 0; };

			auto has_zero_value = [&](std::vector<real> v, real d)
			{
				return d * std::any_of(v.begin(), v.end(), equal_zero);
			};

			apply(3, 3, 3, has_zero_value);
		}

		// Given a kernel dimension, and a "folding" function, apply the function to all kernel
		// groups of the image
		void apply(size_t nx, size_t ny, size_t nz, std::function<real(std::vector<real>, real)> f)
		{
			// Confirm that kernel size is odd
			assert(nx % 2);
			assert(ny % 2);
			assert(nz % 2);

			cube<real> temp(header.nx, header.ny, header.nz);

			size_t mid_x = nx / 2;
			size_t mid_y = ny / 2;
			size_t mid_z = nz / 2;

			// Loop through all kernel groups

			for (size_t i = 0; i < header.nx; i++)
			{
				for (size_t j = 0; j < header.ny; j++)
				{
					for (size_t k = 0; k < header.nz; k++)
					{
						std::vector<real> kernel_function_input;

						for (size_t I = 0; I < nx; I++)
						{
							for (size_t J = 0; J < ny; J++)
							{
								for (size_t K = 0; K < nz; K++)
								{
									real d = map(i + I - mid_x, j + J - mid_y, k + K - mid_z);

									kernel_function_input.push_back(d);
								}
							}
						}

						temp(i, j, k) = f(kernel_function_input, map(i, j, k));
					}
				}
			}

			map = temp;
		}

		void apply(std::function<real(real)> f)
		{
			for (size_t i = 0; i < header.nx; i++)
			{
				for (size_t j = 0; j < header.ny; j++)
				{
					for (size_t k = 0; k < header.nz; k++)
					{
						real density = map(i, j, k);
						map(i, j, k) = f(density);
					}
				}
			}
		}

		// If possible, chops the total box to reduce file size
		void minimize()
		{
			// Identify min_x, max_x, min_y, etc...

			int min_x = INT_MAX;
			int min_y = INT_MAX;
			int min_z = INT_MAX;

			int max_x = INT_MIN;
			int max_y = INT_MIN;
			int max_z = INT_MIN;

			for (int i = 0; i < (int)header.nx; i++)
			{
				for (int j = 0; j < (int)header.ny; j++)
				{
					for (int k = 0; k < (int)header.nz; k++)
					{
						real density = map(i, j, k);

						if (density != 0)
						{
							min_x = i < min_x ? i : min_x;
							min_y = j < min_y ? j : min_y;
							min_z = k < min_z ? k : min_z;

							max_x = i > max_x ? i : max_x;
							max_y = j > max_y ? j : max_y;
							max_z = k > max_z ? k : max_z;
						}
					}
				}
			}

			// At this point we have the subset of the original cube.  We need
			// to carefully ensure we do not break the alignment.

			int new_nx = max_x - min_x;
			int new_ny = max_y - min_y;
			int new_nz = max_z - min_z;

			auto new_map = cube<real>(new_nx, new_ny, new_nz);

			for (int i = 0; i < new_nx; i++)
			{
				for (int j = 0; j < new_ny; j++)
				{
					for (int k = 0; k < new_nz; k++)
					{
						new_map(i, j, k) = map(i + min_x, j + min_y, k + min_z);
					}
				}
			}

			header.nx = new_nx;
			header.ny = new_ny;
			header.nz = new_nz;

			header.xorigin += float(min_x) * scale;
			header.yorigin += float(min_y) * scale;
			header.zorigin += float(min_z) * scale;

			map = new_map;
		}

		// Adds zero-padding to box on all six sides
		void pad(size_t num)
		{
			size_t new_nx = header.nx + num * 2;
			size_t new_ny = header.ny + num * 2;
			size_t new_nz = header.nz + num * 2;

			auto new_map = cube<real>(new_nx, new_ny, new_nz);

			for (size_t i = num; i < new_nx - num; i++)
			{
				for (size_t j = num; j < new_ny - num; j++)
				{
					for (size_t k = num; k < new_nz - num; k++)
					{
						new_map(i, j, k) = map(i - num, j - num, k - num);
					}
				}
			}

			header.nx = (unsigned int)new_nx;
			header.ny = (unsigned int)new_ny;
			header.nz = (unsigned int)new_nz;

			header.xorigin -= num * scale;
			header.yorigin -= num * scale;
			header.zorigin -= num * scale;

			map = new_map;
		}

		cube<real>::voxel transformVoxel(cube<real>::voxel vox)
		{
			// Compute "real-space" coordinates from (i,j,k)
			vox.x = vox.x * scale + header.xorigin;
			vox.y = vox.y * scale + header.yorigin;
			vox.z = vox.z * scale + header.zorigin;

			return vox;
		}

		cube<real>::voxel getTransformedVoxel(size_t n)
		{
			cube<real>::voxel vox = map.getVoxel(n);

			return transformVoxel(vox);
		}

		void apply(cube<real> &kernel)
		{
			// Confirm that kernel size is odd
			assert(kernel.nx() % 2);
			assert(kernel.ny() % 2);
			assert(kernel.nz() % 2);

			cube<real> temp(header.nx, header.ny, header.nz);

			size_t mid_x = kernel.nx() / 2;
			size_t mid_y = kernel.ny() / 2;
			size_t mid_z = kernel.nz() / 2;

			for (size_t i = 0; i < header.nx; i++)
			{
				for (size_t j = 0; j < header.ny; j++)
				{
					for (size_t k = 0; k < header.nz; k++)
					{
						real sum_product = 0;

						for (size_t I = 0; I < kernel.nx(); I++)
						{
							for (size_t J = 0; J < kernel.ny(); J++)
							{
								for (size_t K = 0; K < kernel.nz(); K++)
								{
									real d = map(i + I - mid_x, j + J - mid_y, k + K - mid_z);
									real kernel_element = kernel(I, J, K);

									sum_product += d * kernel_element;
								}
							}
						}

						temp(i, j, k) = sum_product;
					}
				}
			}

			map = temp;
		}

		// Sets elements on border to zero.  How many layers to remove is "layers".  0 "layers"
		// results in no effect.
		void trim(size_t layers)
		{
			for (size_t i = 0; i < header.nx; i++)
			{
				for (size_t j = 0; j < header.ny; j++)
				{
					for (size_t k = 0; k < header.nz; k++)
					{
						if (i < layers || j < layers || k < layers ||
							i >= header.nx - layers ||
							j >= header.ny - layers ||
							k >= header.nz - layers)
						{
							map(i, j, k) = 0;
						}
					}
				}
			}
		}

		pdb convertToPDB(real threshold)
		{
			pdb result;

			real x_origin = header.xorigin;
			real y_origin = header.yorigin;
			real z_origin = header.zorigin;

			real avg_x = (real)header.nx / (real)2;
			real avg_y = (real)header.ny / (real)2;
			real avg_z = (real)header.nz / (real)2;

			for (size_t i = 0; i < header.nx; i++)
			{
				for (size_t j = 0; j < header.ny; j++)
				{
					for (size_t k = 0; k < header.nz; k++)
					{
						real d = map(i, j, k);

						if (d > threshold)
						{
							real x = i * scale + header.xorigin;
							real y = j * scale + header.yorigin;
							real z = k * scale + header.zorigin;

							std::string atom = "H";

							result.emplace_back(x, y, z, atom);
						}
					}
				}
			}

			return result;
		}

		void applyThreshold(real threshold)
		{
			for (size_t i = 0; i < header.nx; i++)
			{
				for (size_t j = 0; j < header.ny; j++)
				{
					for (size_t k = 0; k < header.nz; k++)
					{
						if (map(i, j, k) < threshold)
						{
							map(i, j, k) = 0;
						}
					}
				}
			}
		}

		std::tuple<mrc, mrc> crop(pdb& pdb, real cropping_dist)
		{
			assert(pdb.size() > 0 && "Attempted MRC crop on empty PDB");

			mrc near_map(*this);
			mrc far_map(*this);

			for (size_t i = 0; i < header.nx; i++)
			{
				for (size_t j = 0; j < header.ny; j++)
				{
					for (size_t k = 0; k < header.nz; k++)
					{
						// Compute "real-space" coordinates from (i,j,k)
						real x = i * scale + header.xorigin;
						real y = j * scale + header.yorigin;
						real z = k * scale + header.zorigin;

						// For each point in given pdb, compute distance.  Find min_distance
						real min_dist = std::numeric_limits<real>::infinity();
						size_t min_index = 0;

						for (size_t m = 0; m < pdb.size(); m++)
						{
							point p = pdb[m];
							real dist_sq = p.distSq(x, y, z);

							if (dist_sq < min_dist)
							{
								min_dist = dist_sq;
								min_index = m;
							}
							
						}

						min_dist = sqrt(min_dist); // Get actual distance (not squared)

						if (min_dist < cropping_dist)
						{
							far_map.map(i, j, k) = 0;
						}
						else
						{
							near_map.map(i, j, k) = 0;
						}
					}
				}
			}

			return std::make_tuple(near_map, far_map);
		}

		void merge(std::vector<mrc> mrc_stack)
		{
			for (size_t i = 0; i < header.nx; i++)
			{
				for (size_t j = 0; j < header.ny; j++)
				{
					for (size_t k = 0; k < header.nz; k++)
					{
						map(i, j, k) = 0;
					}
				}
			}

			for (size_t Z = 0; Z < mrc_stack.size(); Z++)
			{
				for (size_t i = 0; i < header.nx; i++)
				{
					for (size_t j = 0; j < header.ny; j++)
					{
						for (size_t k = 0; k < header.nz; k++)
						{
							if (mrc_stack[Z].map(i, j, k) > 0)
							{
								map(i, j, k) = mrc_stack[Z].map(i, j, k);
							}
						}
					}
				}
			}
		}

		std::vector<mrc> clusterViaConnectivity(int num_required_voxels = 0)
		{
			// Create a copy of the map, grouping voxels by setting their densities to successive
			// negative integers, representing which cluster they belong to.

			mrc groups = *this;

			// Keep track of voxels yet to be categorized by pushing them onto a stack.  So, we can
			// always get an uncategorized one if one exists.

			struct index
			{
				index(size_t i, size_t j, size_t k) : i(i), j(j), k(k) {};
				size_t i, j, k;
			};

			std::stack<index> uncategorized;

			for (size_t i = 0; i < header.nx; i++)
			{
				for (size_t j = 0; j < header.ny; j++)
				{
					for (size_t k = 0; k < header.nz; k++)
					{
						if (map(i, j, k) > 0)
						{
							index ind(i, j, k);
							uncategorized.push(ind);
						}
					}
				}
			}

			// Create a stack of indices to visit next.  (Recursive to iterative form)

			int iterations = 0;

			real current_group = (real)-1;

			while (uncategorized.size())
			{
				std::stack<index> add_to_group;

				add_to_group.push(uncategorized.top());
				uncategorized.pop();

				bool group_had_elements = false;

				while (add_to_group.size())
				{
					index ind = add_to_group.top();
					add_to_group.pop();

					auto &d = groups.map(ind.i, ind.j, ind.k);

					if (d > 0)
					{
						d = current_group;

						// Iterate through 6-neighborhood recursively

						ind.i += 1;
						add_to_group.push(ind);
						ind.i -= 2;
						add_to_group.push(ind);
						ind.i += 1;

						ind.j += 1;
						add_to_group.push(ind);
						ind.j -= 2;
						add_to_group.push(ind);
						ind.j += 1;

						ind.k += 1;
						add_to_group.push(ind);
						ind.k -= 2;
						add_to_group.push(ind);
						ind.k += 1;

						group_had_elements = true;
					}
				}

				if (group_had_elements)
				{
					current_group -= (real)1;
				}
			}

			int num_of_groups = (int)std::round(-current_group);

			// Now all voxels in groups have a density corresponding to their connected cluster.
			// Loop through and copy each voxel to the corresponding map in the return vector.


			// Create copy of this map with density all zero to initialize cluster vector
			mrc copy = *this;
			copy.map = cube<real>(header.nx, header.ny, header.nz);

			std::vector<mrc> clusters(num_of_groups, copy);

			for (size_t i = 0; i < header.nx; i++)
			{
				for (size_t j = 0; j < header.ny; j++)
				{
					for (size_t k = 0; k < header.nz; k++)
					{
						auto &d = groups.map(i, j, k);

						if (d < 0)
						{
							int index = (int)std::round(-d) - 1;

							clusters[index].map(i, j, k) = map(i, j, k);
						}
					}
				}
			}

			std::vector<mrc> filtered;

			// Loop through clusters, including them in filtered vector only if they have a number
			// of voxels greater than what input specifies

			for (size_t i = 0; i < clusters.size(); i++)
			{
				mrc &m = clusters[i];

				int num_voxels = 0;

				for (size_t i = 0; i < header.nx; i++)
				{
					for (size_t j = 0; j < header.ny; j++)
					{
						for (size_t k = 0; k < header.nz; k++)
						{
							if (m.map(i, j, k) > 0)
							{
								num_voxels++;
							}
						}
					}
				}

				if (num_voxels > num_required_voxels)
				{
					filtered.push_back(m);
				}
			}

			return filtered;
		}

		// This function is similar to crop, but is based on a cylindrical model instead
		// of a spherical one.  The pdb input should be a linear sequence of positions with
		// nearby values.
		// The trick is very easy: we exclude the voxels that project to the endpoints from
		// consideration.
		std::tuple<mrc, mrc> cylinderCrop(pdb& pdb, real cropping_dist)
		{
			mrc near_map(*this);
			mrc far_map(*this);

			for (size_t i = 0; i < header.nx; i++)
			{
				for (size_t j = 0; j < header.ny; j++)
				{
					for (size_t k = 0; k < header.nz; k++)
					{
						// Compute "real-space" coordinates from (i,j,k)
						real x = i * scale + header.xorigin;
						real y = j * scale + header.yorigin;
						real z = k * scale + header.zorigin;

						// For each point in given pdb, compute distance.  Find min_distance
						real min_dist = std::numeric_limits<real>::infinity();
						size_t min_index = 0;

						for (size_t m = 0; m < pdb.size(); m++)
						{
							point p = pdb[m];
							real dist_sq = p.distSq(x, y, z);

							if (dist_sq < min_dist)
							{
								min_dist = dist_sq;
								min_index = m;
							}
						}

						min_dist = sqrt(min_dist); // Get actual distance (not squared)

						bool is_side_voxel = false;

						if (min_index == 0 || min_index == pdb.size() - 1)
						{
							is_side_voxel = true;
						}

						if (min_dist < cropping_dist && !is_side_voxel)
						{
							far_map.map(i, j, k) = 0;
						}
						else
						{
							near_map.map(i, j, k) = 0;
						}
					}
				}
			}

			return std::make_tuple(near_map, far_map);
		}
	};
}
