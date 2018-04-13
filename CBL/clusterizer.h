#pragma once

#include <vector>
#include <functional>
#include <random>
#include <map>

#include "core.h"

namespace cbl
{
	template <typename T>
	class clusterizer
	{
	public:

		clusterizer() {};
		
		clusterizer(std::vector<T> p, std::function<real(T, T)> D) : points(p), dist_func(D)
		{
			std::random_device rd;
			rng = std::mt19937(rd());

			cluster_of_point.resize(points.size());
			distance_of_point.resize(points.size());
		}

		// For all p, assigns membership m that minimizes distance function D.
		void updateClusterMemberships()
		{
			cluster_values.clear();
			cluster_values.resize(clusters.size());

			real total_dist = 0;

			// For all p, assign to closest cluster and save basic distance
			for (size_t i = 0; i < points.size(); i++)
			{
				auto p = points[i];

				std::vector<real> dist;

				real min_dist = Infinity;
				size_t min_index = -1;
				for (size_t j = 0; j < clusters.size(); j++)
				{
					auto c = clusters[j];
					real dist = dist_func(p, c);
					if (dist < min_dist)
					{
						min_dist = dist;
						min_index = j;
					}
				}

				cluster_of_point[i] = min_index;
				distance_of_point[i] = min_dist;

				total_dist += min_dist;

				cluster_values[min_index].push_back(p);
			}

			// Save total distance to object state
			total_distance = total_dist / points.size();
		}

		void setClustersBasedOnMembership()
		{
			// Given a set of memberships (i.e. voxel sets which form clusters), calculate new
			// cluster means

			clusters.clear();

			for (size_t i = 0; i < cluster_values.size(); i++)
			{
				// Calculate mean of point set

				auto v = cluster_values[i];

				T mean;

				for (size_t j = 0; j < v.size(); j++)
				{
					mean += v[j];
				}

				mean /= (real)v.size();

				clusters.push_back(mean);
			}
		}

		// K-means++ Improved Initialization Algorithm
		// https://en.wikipedia.org/wiki/K-means++
		// Returns a vector of points that correspond to initial k clusters
		void kmeans_pp_initialize(size_t k)
		{
			assert(k > 0 && "Passed 0 clusters to k-means++ algorithm");
			assert(points.size() > 0 && "Called k-means++ on zero points");

			clusters.clear();

			// Uniformly choose a point to serve as first cluster position
			std::uniform_int_distribution<int> uni(0, points.size() - 1);
			clusters.push_back(points[uni(rng)]);

			while (clusters.size() < k)
			{
				updateClusterMemberships();

				// Add up all distance, and generate random number in range 0 to sum distance
				// Then iteratively subtract point distances from generated number until number
				// is less than epsilon.  Then current index is new cluster

				// For now, uniformly generate until I know how to fix
				// clusters.push_back(points[uni(rng)]);

				real sum_dist = 0;
				for (auto d : distance_of_point)
				{
					sum_dist += d;
				}

				std::uniform_real_distribution<real> uni_real(0.0, sum_dist);

				real random_number = uni_real(rng);
				size_t d_index = 0;
				while (random_number > 0.00001)
				{
					d_index++;
					random_number -= distance_of_point[d_index];
				}

				clusters.push_back(points[d_index]);
			}

			updateClusterMemberships();
		}

		// Lloyd's Iterative K-means Heuristic
		// https://en.wikipedia.org/wiki/Lloyds_algorithm
		// Iteratively assigns clusters, converging to a centroidal Voronoi tesselation
		void lloyds_algorithm(size_t k, int iterations = 10)
		{
			clusterizer best_so_far;

			for (int i = 0; i < iterations; i++)
			{
				kmeans_pp_initialize(k);
				std::cout << "Iteration " << i + 1 << " Start: " << total_distance << std::endl;

				real epsilon = (real)0.0001;
				real previous_distance;

				do
				{
					previous_distance = total_distance;
					setClustersBasedOnMembership();
					updateClusterMemberships();
				}
				while (previous_distance - total_distance > epsilon);

				std::cout << "Iteration " << i + 1 << " Final: " << total_distance << std::endl;

				// Save best result found so far
				if (total_distance < best_so_far.total_distance)
				{
					best_so_far = *this;
				}
			}

			*this = best_so_far;

			std::cout << std::endl;
			std::cout << "Lloyd's k=" << k << ": " << total_distance << std::endl;
		}

		void elbow_method(size_t max_num_clusters)
		{
			std::vector<real> result_scores;
			std::vector<real> percent_diff;
			std::vector<clusterizer> results;

			for (size_t k = 1; k < max_num_clusters + 1; k++)
			{
				lloyds_algorithm(k);
				result_scores.push_back(total_distance);
				results.push_back(*this);
			}

			for (size_t i = 1; i < result_scores.size() - 1; i++)
			{
				real prev = result_scores[i - 1];
				real next = result_scores[i + 1];
				real cur = result_scores[i];
				percent_diff.push_back((prev - cur) / (cur - next));
			}

			size_t max_index = std::max_element(percent_diff.begin(), percent_diff.end()) - 
				percent_diff.begin();
			size_t max_k = max_index + 2;

			std::cout << "Elbow K: " << max_k << std::endl;

			*this = results[max_index + 1];
		}

		std::vector<std::vector<T>> getClusters()
		{
			return cluster_values;
		}

		std::vector<T> clusterMeans()
		{
			return clusters;
		}

	private:

		std::mt19937 rng;

		std::vector<T> points;
		std::function<real(T, T)> dist_func;

		// Following members updated on updateClusterMembership call.

		std::vector<T> clusters; // Clusters
		std::vector<size_t> cluster_of_point; // Cluster index corresponding to point
		std::vector<real> distance_of_point; // Distance corresponding to point
		std::vector<std::vector<T>> cluster_values; // Copy of all member points of each cluster

		real total_distance = Infinity;
	};
}

