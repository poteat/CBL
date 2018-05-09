#pragma once

#include <iostream>
#include <map>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef float real;

class strand
{
private:
	std::map<int, Eigen::Vector2f> data;

public:
	void push_back(Eigen::Vector2f v)
	{
		int pos = data.empty() ? 0 : (data.rbegin()->first + 1);
		data[pos] = v;
	}

	void push_front(Eigen::Vector2f v)
	{
		int pos = data.empty() ? 0 : (data.begin()->first - 1);
		data[pos] = v;
	}

	void clear()
	{
		data.clear();
	}

	std::map<int, Eigen::Vector2f>::iterator begin()
	{
		return data.begin();
	}

	std::map<int, Eigen::Vector2f>::iterator end()
	{
		return data.end();
	}

	size_t size()
	{
		return data.size();
	}

	void crop(int min, int max)
	{
		std::map<int, Eigen::Vector2f> new_data;

		for (auto it = data.begin(); it != data.end(); it++)
		{
			if (it->first >= min && it->first < max)
			{
				new_data[it->first] = it->second;
			}
		}

		data = new_data;
	}
};

class strandMap
{
private:
	std::map<int, strand> data;

public:
	void push_back(strand s)
	{
		int pos = data.empty() ? 0 : (data.rbegin()->first + 1);
		data[pos] = s;
	}

	void push_front(strand s)
	{
		int pos = data.empty() ? 0 : (data.begin()->first - 1);
		data[pos] = s;
	}

	void clear()
	{
		data.clear();
	}

	std::map<int, strand>::iterator begin()
	{
		return data.begin();
	}

	std::map<int, strand>::iterator end()
	{
		return data.end();
	}

	size_t size()
	{
		return data.size();
	}

	size_t sizeAll()
	{
		int sum = 0;
		for (auto it = data.begin(); it != data.end(); it++)
		{
			sum += it->second.size();
		}

		return sum;
	}
};

class strands
{
public:
	// Surface which strands will be drawn upon
	surface * surf;

	// Bounding-box of surface patch that voxels project on
	real min_t, max_t, min_u, max_u;

	// Linear segments of 2d positions which represent the current strand result.
	strandMap strand_map;

	// Center point
	Eigen::Vector2f center;

	// Sampling rate
	real sampling_distance = (real).01;

	// Three parameters that control shape of strand.  These represent the current strand model.
	real angle, offset, gap;

	strands(surface *surf, real angle = 0, real offset = 0, real gap = 4.5) :
		surf(surf), angle(angle), offset(offset), gap(gap)
	{
		updateBoundingBox();

		updateStrandMap(angle, offset, gap);
	}

	Eigen::Vector2f euclideanShift(Eigen::Vector2f v, real angle, real dist)
	{
		// Shift a sample coordinate "dist" distance units in 3-space

		real epsilon = (real).001;

		real t_delta = cos(angle) * epsilon;
		real u_delta = sin(angle) * epsilon;

		Eigen::Vector2f sampling;
		sampling << t_delta, u_delta;

		Eigen::Vector2f v2 = v + sampling;

		auto p1 = surf->sample(v(0), v(1));
		auto p2 = surf->sample(v2(0), v2(1));

		Eigen::Vector3f delta = (p2 - p1);

		delta.normalize();

		Eigen::Vector3f shifted = p1 + dist * delta;

		auto shifted_samp = surf->desample(shifted);

		return shifted_samp;
	}

	void updateBoundingBox()
	{
		// Loop through surface projected points and calculate min, max of t,u
		// As well, calculate average

		min_t = Infinity;
		max_t = -Infinity;
		min_u = Infinity;
		max_u = -Infinity;

		real sum_t = 0;
		real sum_u = 0;

		for (int i = 0; i < surf->point_projections.cols(); i++)
		{
			real t = surf->point_projections(0, i);
			real u = surf->point_projections(1, i);

			min_t = t < min_t ? t : min_t;
			max_t = t > max_t ? t : max_t;
			min_u = u < min_u ? u : min_u;
			max_u = u > max_u ? u : max_u;

			sum_t += t;
			sum_u += u;
		}

		sum_t /= surf->point_projections.cols();
		sum_u /= surf->point_projections.cols();

		center << sum_t, sum_u;

		std::cout << "Bounds: min_t: " << min_t << " max_t: " << max_t << std::endl;
		std::cout << "        min_u: " << min_u << " max_u: " << max_u << std::endl;
		std::cout << std::endl;
	}



	// On execution, updates strand modeling result.
	void updateStrandMap(real angle, real offset, real gap)
	{
		angle = angle;
		offset = offset;
		gap = gap;

		// First, we construct the central strand.

		std::cout << "Center Coordinate: " << center(0) << ", " << center(1) << std::endl;
		std::cout << std::endl;

		auto in_box = [&](Eigen::Vector2f &p)
		{
			real t = p(0);
			real u = p(1);
			return t >= min_t && t <= max_t && u >= min_u && u <= max_u;
		};

		real t_delta = cos(angle) * sampling_distance;
		real u_delta = sin(angle) * sampling_distance;

		Eigen::Vector2f sampling;
		sampling << t_delta, u_delta;

		real t_perp_delta = (real)cos(angle + M_PI / 2.0) * gap;
		real u_perp_delta = (real)sin(angle + M_PI / 2.0) * gap;

		Eigen::Vector2f perp_sampling;
		perp_sampling << t_perp_delta, u_perp_delta;

		// Arbitrarily increment delta values until we're outside of the box, adding to initial
		// strand each step.

		auto midseek = center;

		// Apply offset
		midseek = euclideanShift(midseek, angle + (real)M_PI / (real)2.0, offset);

		while (in_box(midseek))
		{
			strand strand;

			auto seek = midseek;

			strand.push_back(seek);
			while (in_box(seek))
			{
				seek += sampling;
				strand.push_back(seek);
			}

			seek = midseek;
			while (in_box(seek))
			{
				seek -= sampling;
				strand.push_front(seek);
			}

			strand_map.push_back(strand);

			midseek = euclideanShift(midseek, angle + (real)M_PI / (real)2.0, gap);
		}

		midseek = center;

		center -= perp_sampling;

		while (in_box(midseek))
		{
			strand strand;

			auto seek = midseek;

			strand.push_back(seek);
			while (in_box(seek))
			{
				seek += sampling;
				strand.push_back(seek);
			}

			seek = midseek;
			while (in_box(seek))
			{
				seek -= sampling;
				strand.push_front(seek);
			}

			strand_map.push_front(strand);

			midseek = euclideanShift(midseek, angle + (real)M_PI / (real)2.0, -gap);
		}

		// Cropping process.  Loop through each strand and save the maximum and minimum indices
		// that correspond to a sample that is close to density-based samples.

		real distance_threshold = (real)0.01;

		for (auto it = strand_map.begin(); it != strand_map.end(); it++)
		{
			auto s = it->second;

			int min_near = INT_MAX;
			int max_near = INT_MIN;

			for (auto jt = s.begin(); jt != s.end(); jt++)
			{
				Eigen::Vector2f v = jt->second;

				real min_dist = Infinity;

				for (int k = 0; k < surf->point_projections.cols(); k++)
				{
					Eigen::Vector2f voxel_proj = surf->point_projections.col(k);

					real dist = (v - voxel_proj).squaredNorm();

					if (dist < min_dist)
					{
						min_dist = dist;
					}
				}

				min_dist = sqrt(min_dist);

				if (min_dist < distance_threshold)
				{
					if (jt->first < min_near)
					{
						min_near = jt->first;
					}

					if (jt->first > max_near)
					{
						max_near = jt->first;
					}
				}
			}

			std::cout << min_near << " " << max_near << " " << it->first << std::endl;
			it->second.crop(min_near, max_near);
		}
	}


	void writePointsToPDB(std::string file_name)
	{
		std::ofstream file(file_name);

		Eigen::Matrix3Xf points;
		std::cout << strand_map.sizeAll() << std::endl;
		points.resize(3, strand_map.sizeAll());

		size_t seek = 0;
		for (auto it = strand_map.begin(); it != strand_map.end(); it++)
		{
			auto s = it->second;
			for (auto jt = s.begin(); jt != s.end(); jt++)
			{
				real t = jt->second(0);
				real u = jt->second(1);
				points.col(seek) = surf->sample(t, u);

				seek++;
			}
		}

		// Unrotate
		Eigen::Vector3f rotation_vector = surf->map->rotation_vector;
		float rotation_angle = -(surf->map->rotation_angle);

		Eigen::AngleAxisf aa(rotation_angle, rotation_vector);
		Eigen::Quaternionf q(aa);
		Eigen::Matrix3f rotation_matrix;
		rotation_matrix = q;

		Eigen::MatrixX3f tmp = points.transpose() * rotation_matrix;
		points = tmp.transpose();

		// Denormalize with respect to mean
		points = points.colwise() + surf->map->mean_vector;

		for (int i = 0; i < points.cols(); i++)
		{
			float x = points(0, i);
			float y = points(1, i);
			float z = points(2, i);

			x += surf->map->header.xorigin;
			y += surf->map->header.yorigin;
			z += surf->map->header.zorigin;

			file << "ATOM " << std::setw(6) << std::right << i << \
				"  N   HOH A   1    " << std::setprecision(5) \
				<< std::setw(7) << std::right << x << " " \
				<< std::setw(7) << std::right << y << " " \
				<< std::setw(7) << std::right << z << std::endl;
		}
	}
};