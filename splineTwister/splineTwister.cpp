#include <iostream>
#include <fstream>
#include <tuple>

#include "surface.h"

int main(int argc, char *argv[])
{
	mrc map(argv[1]);

	map.grassfire(4);
	map.convertToPoints();
	map.normalizeOrientation();

	std::string s;
	s = "_fire.mrc";
	map.write(argv[1] + s);

	surface surf(&map);

	surf.optimizeSurface();

	int max_surface_degree = 4;

	for (int n = 3; n <= max_surface_degree; n++)
	{
		surf.elevateDegree();
		surf.optimizeSurface();
	}

	s = "_surf.pdb";
	surf.drawRenderGrid(argv[1] + s);
}