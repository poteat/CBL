#pragma once

#include <sstream>

#include "core.h"
#include "pdb.h"
#include "mrc.h"

#include "Eigen/Dense"

// the axis object represents a single interpolated helix or strand structure.

namespace cbl
{
	class plane
	{
	public:
		plane(mrc &map)
		{
			// We specify the plane by giving a voxel map, with the resulting plane being a "plane 
			// of best fit", that minimizes the RMSE projection distance.

			// The map should have been thresholded before-hand.

			pdb point_set = map.convertToPDB(0);

			for (size_t i = 0; i < point_set.size(); i++)
			{
				points.push_back(point_set[i]);
			}


			// Use the Eigen library to find the parameters needed to model this plane

			// Count how many voxels in map are nonzero

			size_t count = points.size();

			// Cast point 'std::vector' as actual vector of R^3

			point_matrix.resize(3, count);

			for (size_t i = 0; i < points.size(); i++)
			{
				point_matrix(0, i) = (float)points[i].x;
				point_matrix(1, i) = (float)points[i].y;
				point_matrix(2, i) = (float)points[i].z;
			}
			

			mean_vector = point_matrix.rowwise().mean();
			point_matrix = point_matrix.colwise() - mean_vector;

			int opt = Eigen::ComputeFullU | Eigen::ComputeFullV;
			Eigen::JacobiSVD<Eigen::Matrix3Xf> svd = point_matrix.jacobiSvd(opt);


			normal_vector = svd.matrixU().col(2);

			std::cout << normal_vector << std::endl;

			real error2 = 0;




			t_basis = svd.matrixU().col(0);
			u_basis = svd.matrixU().col(1);
			
			fixBasis();
		}

		void fixBasis()
		{
			Eigen::Vector3f col0 = t_basis;
			Eigen::Vector3f col1 = u_basis;
			Eigen::Vector3f col2 = normal_vector;

			real error0 = 0;
			normal_vector = col0;
			t_basis = col1;
			u_basis = col2;
			for (size_t i = 0; i < points.size(); i++)
			{
				point p = points[i];
				Eigen::Vector3f v(p.x, p.y, p.z);

				real t = (v - mean_vector).dot(t_basis);
				real u = (v - mean_vector).dot(u_basis);

				Eigen::Vector3f sample = (float)t * t_basis + (float)u * u_basis;
				sample += mean_vector;

				real dist = (v - sample).squaredNorm();
				error0 += dist;
			}

			real error1 = 0;
			normal_vector = col1;
			t_basis = col0;
			u_basis = col2;
			for (size_t i = 0; i < points.size(); i++)
			{
				point p = points[i];
				Eigen::Vector3f v(p.x, p.y, p.z);

				real t = (v - mean_vector).dot(t_basis);
				real u = (v - mean_vector).dot(u_basis);

				Eigen::Vector3f sample = (float)t * t_basis + (float)u * u_basis;
				sample += mean_vector;

				real dist = (v - sample).squaredNorm();
				error1 += dist;
			}

			real error2 = 0;
			normal_vector = col2;
			t_basis = col0;
			u_basis = col1;
			for (size_t i = 0; i < points.size(); i++)
			{
				point p = points[i];
				Eigen::Vector3f v(p.x, p.y, p.z);

				real t = (v - mean_vector).dot(t_basis);
				real u = (v - mean_vector).dot(u_basis);

				Eigen::Vector3f sample = (float)t * t_basis + (float)u * u_basis;
				sample += mean_vector;

				real dist = (v - sample).squaredNorm();
				error2 += dist;
			}


			std::cout << error0 << " " << error1 << " " << error2 << std::endl;
			

			if (error0 < error1 && error0 < error2)
			{
				std::cout << "Error0 min" << std::endl;
				normal_vector = col0;
				t_basis = col1;
				u_basis = col2;
			}
			else if (error1 < error0 && error1 < error2)
			{
				std::cout << "Error1 min" << std::endl;
				normal_vector = col1;
				t_basis = col0;
				u_basis = col2;
			}
			else if (error2 < error0 && error2 < error1)
			{
				std::cout << "Error2 min" << std::endl;
				normal_vector = col2;
				t_basis = col0;
				u_basis = col1;
			}

		}

		void write(std::string file_name)
		{
			std::ofstream file(file_name);
			assert(file && "Could not open axis file for writing");
			assert(points.size() && "Attempted to write an empty axis file");

			// Find boundary of point set projected onto plane so we may sample a rectangular patch

			real min_t = Infinity, min_u = Infinity;
			real max_t = -Infinity, max_u = -Infinity;

			for (size_t i = 0; i < points.size(); i++)
			{
				point p = points[i];
				Eigen::Vector3f v(p.x, p.y, p.z);

				v -= mean_vector;

				//Eigen::Vector3f projected = v - normal_vector.dot(v - mean_vector) * normal_vector;

				//Eigen::Vector3f projected = v - (v - mean_vector).dot(normal_vector)*normal_vector;

				// Now we find the t and u multiples that correspond with projected

				Eigen::Vector3f projected = v;

				real t = projected.dot(t_basis);
				real u = projected.dot(u_basis);

				min_t = t < min_t ? t : min_t;
				min_u = u < min_u ? u : min_u;
				max_t = t > max_t ? t : max_t;
				max_u = u > max_u ? u : max_u;
			}

			// Sample points on the plane as a 5x5 angstrom square aligned to basis

			pdb sampling_set;

			real res = (real)1.0;

			for (real t = min_t; t <= max_t; t += res)
			{
				for (real u = min_u; u <= max_u; u += res)
				{
					Eigen::Vector3f sample = (float)t * t_basis + (float)u * u_basis;

					sample += mean_vector;

					sampling_set.emplace_back(sample.x(), sample.y(), sample.z());
				}
			}


			sampling_set.write(file_name);
		}


	public:

		std::vector<point> points;
		Eigen::Matrix3Xf point_matrix;

		// Two vectors necessary to specify a plane
		Eigen::Vector3f mean_vector;
		Eigen::Vector3f normal_vector;

		Eigen::Vector3f t_basis;
		Eigen::Vector3f u_basis;
	};
}
