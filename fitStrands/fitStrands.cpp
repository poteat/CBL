#include <iostream>
#include <fstream>
#include <tuple>

#include "_surface.h"
#include "_strands.h"

int main(int argc, char *argv[])
{
	// Initialize map via filename string (will open and create self based on file)
	mrc map(argv[1]);

	// Set density threshold to be 2 standard deviations
	map.setThreshold(2 * map.deviation());

	// Convert to "matrix form", i.e. 3 x N matrix representing N 3d points
	map.convertToPoints();

	// Find plane-of-best-fit and rotate points so they lie against X-Y plane
	map.normalizeOrientation();





	// Initialize surface with the density map
	surface surf(&map);

	// Optimize initial 2x2 surface
	surf.optimizeSurface();

	// Until we reach 4x4 surface, iteratively optimize and elevate degree (2x2) -> (3x3) -> (4x4)
	int max_surface_degree = 4;
	for (int n = 3; n <= max_surface_degree; n++)
	{
		surf.elevateDegree();
		surf.optimizeSurface();
	}




	// Initialize strand model with surface model.  By default, initial pose is ang = 0, offset = 0
	strands model(&surf);




	// Output surface and strand results as two pdb files with same base name as input filename
	std::string s = "_surf.pdb";
	surf.drawRenderGrid(argv[1] + s);

	s = "_strands.pdb";
	model.writePointsToPDB(argv[1] + s);
}