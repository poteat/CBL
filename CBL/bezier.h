#pragma once

#include <vector>
#include <initializer_list>

#include "core.h"

namespace cbl
{
	template <class T>
	class bezier
	{
	public:
		box<T> control_points;

		bezier(std::vector<size_t> dimensionality)
		{
			control_points = box<T>(dimensionality);
		}

		T calc(std::initializer_list<real> position)
		{
			// Confirm position within unit 'hyper-square' [0 - 1]
			for (int i = 0; i < position.size(); i++)
			{
				assert(position[i] >= 0 && "Attempted sample of negative position in bezier");
				assert(position[i] <= 1 && "Attempted sample of pos > 1 in bezier");
			}
		}
	};
}