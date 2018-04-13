#include <assert.h>
#include <string>
#include <functional>
#include <iostream>

#include "mrc.h"
#include "off.h"
#include "plane.h"
#include "clusterizer.h"

using namespace cbl;

int main(int argc, char *argv[])
{
	// Confirm we have correct number of arguments (.off, X.mrc)
	int num_args = 1;
	assert(argc == num_args + 1 && "Incorrect num arguments: X.mrc");

	// Save base path for writing files during body
	std::string base_path = argv[1];
	size_t period_pos = base_path.find_last_of('.');
	base_path.resize(period_pos);

	// Read in original map to cluster
	mrc original_map(argv[1]);

	original_map.normalize();
	original_map.applyDeviationThreshold(2.5);
	
	// Build clusterizer instance with euclidean squared distance function
	using voxel = cube<real>::voxel;

	auto dist_func = [](voxel a, voxel b) -> real
	{
		return a.d * b.d * std::pow(a.dist(b), 2);
	};

	clusterizer<voxel> clust(original_map.map.getAllVoxels(), dist_func);
	
	clust.lloyds_algorithm(7);

	auto cluster_memberships = clust.getClusters();
	auto clusters = clust.clusterMeans();

	// Loop through clusters means and add to pdb
	pdb clust_vis;
	for (size_t i = 0; i < clusters.size(); i++)
	{
		auto voxel_trans = original_map.transformVoxel(clusters[i]);

		clust_vis.emplace_back(voxel_trans);
	}
	clust_vis.write(base_path + "_clusterMeans.pdb");

	// Loop through cluster values and write clusters
	for (size_t i = 0; i < clusters.size(); i++)
	{
		mrc cluster = original_map;
		cluster.map.setAllVoxels(cluster_memberships[i]);

		cluster.write(base_path + "_cluster" + std::to_string(i + 1) + ".mrc");
	}

	system("pause");

	return 0;
}