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

		void set(Eigen::Vector3f &v)
		{
			x = v[0];
			y = v[1];
			z = v[2];
		}

		// Rotate theta degrees around {ux, uy, uz} relative to center {_x, _y, _z}
		void rotate(real _x, real _y, real _z, real ux, real uy, real uz, real theta)
		{
			Eigen::Vector3f w = { ux, uy, uz }; // rotation axis
			w.normalize(); // Make sure rotation axis is unit vector

			Eigen::Vector3f c = { _x, _y, _z }; // center of rotation
			Eigen::Affine3f A = Eigen::Translation3f(c) * Eigen::AngleAxisf(theta, w) * 
				Eigen::Translation3f(-c);

			Eigen::Vector3f p = { x, y, z };
			p = A * p;

			set(p);
		};

		real x, y, z;
	};

	real sign(real x)
	{
		return (real)(x > 0.0) - (x < 0.0);
	}

	real distSq(point &a, point &b)
	{
		return pow(a.x - b.x, 2) + pow(a.y - b.y, 2) + pow(a.z - b.z, 2);
	}

	real distToLine(point &a, point &b, point &c)
	{
		real hypotenuse = std::sqrt(distSq(b, c));
		real angle = sin(acos(((b.x - a.x)*(c.x - b.x) + (b.y - a.y)*(c.y - b.y) + (b.z - a.z)*(c.z - b.z))
					/ (std::sqrt(distSq(a, b))*std::sqrt(distSq(b, c)))));

		return hypotenuse * angle;
	}

	real angleBetweenTwoVectors(real x1, real y1, real z1, real x2, real y2, real z2)
	{
		float l1 = sqrt(pow(x1, 2) + pow(y1, 2) + pow(z1, 2));
		float l2 = sqrt(pow(x2, 2) + pow(y2, 2) + pow(z2, 2));

		float dot = x1 * x2 + y1 * y2 + z1 * z2;

		float ang = acos(dot / l1 / l2);

		// Twist direction calculation

		float cross_x = y1 * z2 - z1 * y2;
		float cross_y = z1 * x2 - x1 * z2;
		float cross_z = x1 * y2 - y1 * x2;

		float delta_x = x1 - x2;
		float delta_y = y1 - y2;
		float delta_z = z1 - z2;

		float dot_dir = cross_x * delta_x + cross_y * delta_y + cross_z * delta_z;

		float dir = sign(dot_dir);

		return dir * ang;
	}

	real angleBetweenThreePoints(point &a, point &b, point &c)
	{
		real x1 = b.x - a.x;
		real y1 = b.y - a.y;
		real z1 = b.z - b.z;

		real x2 = c.x - b.x;
		real y2 = c.y - b.y;
		real z2 = c.z - b.z;

		return angleBetweenTwoVectors(x1, y1, z1, x2, y2, z2);
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
			voxel() : x(0), y(0), z(0), d(0) {};
			voxel(real x, real y, real z, T d) : x(x), y(y), z(z), d(d) {};
			operator point() { return point::point(x, y, z); };

			bool operator<(const voxel &rhs) const
			{
				return (x < rhs.x) && (y < rhs.y) && (z < rhs.z);
			};

			voxel& operator+=(const voxel& rhs)
			{
				x += rhs.x;
				y += rhs.y;
				z += rhs.z;
				d += rhs.d;
				
				return *this;
			}

			voxel& operator/=(real div)
			{
				x /= div;
				y /= div;
				z /= div;
				d /= div;

				return *this;
			}

			real dist(const voxel &o) const
			{
				return std::sqrt(std::pow(x - o.x, 2) + std::pow(y - o.y, 2) + 
					std::pow(z - o.z, 2));
			}

			real x, y, z;
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

		voxel getVoxel(size_t n)
		{
			size_t i = n % _nx;
			size_t k = n / _nx / _ny;
			size_t j = (n - _nx * _ny * k - i) / _nx;
			real d = data[n];

			return voxel::voxel((real)i, (real)j, (real)k, d);
		}

		std::vector<voxel> getAllVoxels()
		{
			std::vector<voxel> v;

			for (size_t i = 0; i < data.size(); i++)
			{
				if (data[i] > 0)
				{
					v.push_back(getVoxel(i));
				}
			}

			return v;
		}

		void setAllVoxels(std::vector<voxel> vec)
		{
			std::fill(data.begin(), data.end(), (real)0.0);

			for (auto vox : vec)
			{
				size_t i = (size_t)std::round(vox.x);
				size_t j = (size_t)std::round(vox.y);
				size_t k = (size_t)std::round(vox.z);
				real d = vox.d;

				this->operator()(i, j, k) = d;
			}
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

		T min(bool include_zero = false)
		{
			real min = Infinity;
			for (size_t i = 0; i < data.size(); i++)
			{
				if (include_zero || data[i] > 0)
				{
					if (data[i] < min)
					{
						min = data[i];
					}
				}
			}

			return min;
		}

		T max(bool include_zero = false)
		{
			return *std::max_element(data.begin(), data.end());
		}

		T mean(bool include_zero = false)
		{
			T sum = 0;
			size_t count = 0;

			if (include_zero)
			{
				auto f = [&](T x) {sum += x; count++; };
				std::for_each(data.begin(), data.end(), f);
			}
			else
			{
				auto f = [&](T x) {sum += x * (x > 0); count += x > 0; };
				std::for_each(data.begin(), data.end(), f);
			}

			return sum / count;
		}

		size_t count()
		{
			// Return number of voxels > 0
			return std::count_if(data.begin(), data.end(), [](T x) {return x > 0; });
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
			Eigen::Vector3f w = { ux, uy, uz }; // rotation axis
			w.normalize(); // Make sure rotation axis is unit vector

			Eigen::Vector3f c = { _x, _y, _z }; // center of rotation
			Eigen::Affine3f A = Translation3f(c) * AngleAxisf(theta, w) * Translation3f(-c);

			for (int i = 0; i < _nx; i++)
			{
				for (int j = 0; j < _ny; j++)
				{
					for (int k = 0; k < _nz; k++)
					{

					}
				}
			}

			Eigen::Vector3f p = { x, y, z };
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