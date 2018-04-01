#pragma once

#include <vector>
#include <string>
#include <algorithm>
#include <sstream>
#include <math.h>
#include <iterator>

#include "Eigen/Dense"

namespace cbl
{
	using namespace Eigen;

	typedef float real;
	real Infinity = std::numeric_limits<real>::max();

	class mrc;
	class off;
	class pdb;

	struct point
	{
		point(real x, real y, real z) : x(x), y(y), z(z) {};

		real distSq(real _x, real _y, real _z)
		{
			return pow(x - _x, 2) + pow(y - _y, 2) + pow(z - _z, 2);
		};

		real dist(real _x, real _y, real _z)
		{
			return sqrt(distSq(_x, _y, _z));
		};

		real dist(point &p)
		{
			return sqrt(distSq(p.x, p.y, p.z));
		}

		void set(Vector3f &v)
		{
			x = v[0];
			y = v[1];
			z = v[2];
		}

		// Rotate theta degrees around {ux, uy, uz} relative to center {_x, _y, _z}
		void rotate(real _x, real _y, real _z, real ux, real uy, real uz, real theta)
		{
			Vector3f w = { ux, uy, uz }; // rotation axis
			w.normalize(); // Make sure rotation axis is unit vector

			Vector3f c = { _x, _y, _z }; // center of rotation
			Affine3f A = Translation3f(c) * AngleAxisf(theta, w) * Translation3f(-c);

			Vector3f p = { x, y, z };
			p = A * p;

			set(p);
		};

		real x, y, z;
	};

	real distSq(point &a, point &b)
	{
		return pow(a.x - b.x, 2) + pow(a.y - b.y, 2) + pow(a.z - b.z, 2);
	}

	template <class T>
	class cube
	{
		size_t _nx = 0, _ny = 0, _nz = 0;
		std::vector<T> data;

	public:
		cube() {};

		struct voxel
		{
			voxel(size_t i, size_t j, size_t k, T d) : i(i), j(j), k(k), d(d) {};
			operator point() { return point::point((real)i, (real)j, (real)k); };

			size_t i, j, k;
			T d;
		};

		cube(size_t nx, size_t ny, size_t nz) : _nx(nx), _ny(ny), _nz(nz), data(nx*ny*nz) {}

		T &operator()(size_t i, size_t j, size_t k)
		{
			if (i < _nx && j < _ny && k < _nz)
			{
				return data[i + j*_nx + k*_nx*_ny];
			}
			else
			{
				data[0] = 0;
				return data[0];
			}
		}

		T &operator[](size_t i)
		{
			return data[i];
		}

		voxel voxel(size_t n)
		{
			size_t i = n % _nx;
			size_t k = n / _nx / _ny;
			size_t j = (n - _nx * _ny * k - i) / _nx;
			T d = this->operator()(i, j, k);

			return voxel::voxel(i, j, k, d);
		}

		point point(size_t n)
		{
			size_t i = n % _nx;
			size_t k = n / _nx / _ny;
			size_t j = (n - _nx * _ny * k - i) / _nx;

			return point::point((real)i, (real)j, (real)k);
		}

		size_t size()
		{
			return _nx*_ny*_nz;
		}

		cube(std::string file_name)
		{
			std::ifstream file(file_name);
			assert(file && "Could not open cube file for reading");
			std::string line;

			std::vector<std::string> elements;

			file >> _nx >> _ny >> _nz;

			std::getline(file, line);

			data.resize(_nx*_ny*_nz);

			for (size_t k = 0; k < _nz; k++)
			{
				for (size_t j = 0; j < _ny; j++)
				{
					std::getline(file, line);
					line = removeConsecutive(line, ' ');
					elements = explode(line, ' ');

					for (size_t i = 0; i < elements.size(); i++)
					{
						(*this)(i, j, k) = (real) std::stod(elements[i]);
					}
				}
			}
		}

		void print()
		{
			for (size_t k = 0; k < _nz; k++)
			{
				for (size_t j = 0; j < _ny; j++)
				{
					for (size_t i = 0; i < _nx; i++)
					{
						std::cout << (*this)(i, j, k) << " ";
					}

					std::cout << std::endl;
				}

				std::cout << std::endl;
			}
		}

		T min()
		{
			return *std::min_element(data.begin(), data.end());
		}

		T max()
		{
			return *std::max_element(data.begin(), data.end());
		}

		T mean()
		{
			T sum = 0;
			int count = 0;
			auto f = [&](T x) {sum += x; count += x > 0; };
			std::for_each(data.begin(), data.end(), f);
			sum /= count;
			return sum;
		}

		T variance()
		{
			T mean_val = mean();
			T square_sum = 0;
			T count = 0;
			auto f = [&](T x) {square_sum += (x > 0) * pow(mean_val - x, 2); count += x > 0; };
			std::for_each(data.begin(), data.end(), f);
			square_sum /= count;
			return square_sum;
		}

		T standard_deviation()
		{
			return sqrt(variance());
		}

		void rotate(real _x, real _y, real _z, real ux, real uy, real uz, real theta)
		{
			Vector3f w = { ux, uy, uz }; // rotation axis
			w.normalize(); // Make sure rotation axis is unit vector

			Vector3f c = { _x, _y, _z }; // center of rotation
			Affine3f A = Translation3f(c) * AngleAxisf(theta, w) * Translation3f(-c);

			for (int i = 0; i < _nx; i++)
			{
				for (int j = 0; j < _ny; j++)
				{
					for (int k = 0; k < _nz; k++)
					{

					}
				}
			}

			Vector3f p = { x, y, z };
			p = A * p;

			set(p);
		}

		size_t nx() { return _nx; };
		size_t ny() { return _ny; };
		size_t nz() { return _nz; };
	};

	int binomial(int n, int k)
	{
		real mul = 1;
		for (int i = 1; i < k; i++)
		{
			mul *= (real)(n + 1 - i) / (real)i;
		}
		return (int) round(mul);
	}

	std::string removeConsecutive(std::string s, char c)
	{
		auto f = [c](char l, char r) { return (l == c) && (r == c); };
		auto new_end = std::unique(s.begin(), s.end(), f);
		s.erase(new_end, s.end());
		if (s.size() == 1)
		{
			if (s[0] == c)
			{
				s.clear();
			}
		}
		return s;
	}

	std::vector<std::string> explode(const std::string& s, char delim)
	{
		std::istringstream ss{ s, delim };
		using StrIt = std::istream_iterator<std::string>;
		std::vector<std::string> words{ StrIt{ ss }, StrIt{} };
		return words;
	};

	template <typename T>
	std::vector<T> reverse(const std::vector<T> &v)
	{
		std::vector<T> reverse;
		auto identity = [](int x) { return x; };
		std::transform(v.rbegin(), v.rend(), std::back_inserter(reverse), identity);
		return reverse;
	}

	template <class T>
	void removeConsecutive(std::vector<T> &v, T c)
	{
		auto f = [c](T l, T r) { return (l == c) && (r == c); };
		auto new_end = std::unique(v.begin(), v.end(), f);
		v.erase(new_end, v.end());
		if (v.size() == 1)
		{
			if (v[0] == c)
			{
				v.clear();
			}
		}
	}

	template <class T>
	void trimBeginningIf(std::vector<T> &v, std::function<bool(T)> f)
	{
		v = reverse(v);

		for (auto it = v.rbegin(); it < v.rend(); )
		{
			if (f(it))
			{
				v.erase(it++);
			}
			else
			{
				break;
			}
		}

		v = reverse(v);
	}

	template <class T>
	void trimEndingIf(std::vector<T> &v, std::function<bool(T)> f)
	{
		for (auto it = v.rbegin(); it < v.rend(); )
		{
			if (f(it))
			{
				v.erase(it++);
			}
			else
			{
				break;
			}
		}
	}

	std::string fileType(std::string file_name)
	{
		size_t period_pos = file_name.find_last_of('.');
		file_name.erase(file_name.begin(), file_name.begin() + period_pos);
		return file_name;
	}
}