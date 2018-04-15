#pragma once

#include <vector>
#include <string>
#include <algorithm>
#include <sstream>
#include <math.h>
#include <iterator>
#include <experimental/filesystem>

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

	std::vector<pdb> runAxisComparisonForHelixGeneration(std::string pdb_path)
	{
		using path = std::experimental::filesystem::path;

		// Convert relative to absolute path if necessary, using experimental std lib
		// If there are symbolic links, we also convert those to canonical form
		path symbolic_pdb_path = pdb_path;
		symbolic_pdb_path = std::experimental::filesystem::canonical(symbolic_pdb_path);
		pdb_path = symbolic_pdb_path.string();

		// Assert that the extension is ".pdb"
		assert(symbolic_pdb_path.extension() == ".pdb");

		// Make output directory for results to be populated
		path output_path = symbolic_pdb_path.parent_path().string() + "/output";
		std::experimental::filesystem::create_directory(output_path);

		// Get path stripped of file extension for input into axis-comparsion
		path stripped_path = symbolic_pdb_path.parent_path().string() + "/" +
			symbolic_pdb_path.stem().string();

		std::string command = "leastsquare.exe \"" + stripped_path.string()
			+ "\" \"Empty\" \"Empty\" \"Empty\"";

		// Execute axis comparison
		system(command.c_str());

		// Read all files in this directory that begin with "trueHelix"

		std::vector<cbl::pdb> matched_helix;

		for (auto &p : std::experimental::filesystem::directory_iterator(output_path))
		{
			auto s = p.path().filename().string();

			// If begins with "trueHelix"
			if (s.find("trueHelix") == 0)
			{
				matched_helix.emplace_back(p.path().string());
			}
		}

		// Rename output directory to something we can refer to later

		path new_output_path = symbolic_pdb_path.parent_path().string() + "/"
			+ symbolic_pdb_path.stem().string();

		std::experimental::filesystem::remove_all(new_output_path);

		std::experimental::filesystem::rename(output_path, new_output_path);

		return matched_helix;
	}

	cbl::real applyLineData(mrc map, pdb structure, cbl::real deviation, std::string& data)
	{
		axis line(structure);

		cbl::real deviation_value = map.map.standard_deviation();

		cbl::real threshold_used = deviation_value * deviation;

		std::cout << "# Deviations: " << deviation << "     Threshold Used: " << threshold_used << std::endl;

		map.applyDeviationThreshold(deviation);
		map.setHollowBoundary();

		auto avg = [](std::vector<cbl::real> v)
		{
			cbl::real sum = 0;
			int count = 0;
			auto f = [&](cbl::real x) {sum += x; count += x > 0; };
			std::for_each(v.begin(), v.end(), f);
			sum /= count;
			return sum;
		};

		auto variance = [](std::vector<cbl::real> v, cbl::real avg)
		{
			cbl::real square_sum = 0;
			cbl::real count = 0;
			auto f = [&](cbl::real x) {square_sum += (x > 0) * pow(avg - x, 2); count += x > 0; };
			std::for_each(v.begin(), v.end(), f);
			square_sum /= count;
			if (count > 1)
			{
				return square_sum;
			}
			else
			{
				return (cbl::real) 0.0;
			}
		};

		line.calcError(map, avg);

		line.mergeError(1.0); // Angstroms

		std::ostringstream out;

		line.write(out, threshold_used);

		data += out.str();

		cbl::real alpha = 0.5;
		cbl::real beta = 0.5;

		cbl::real threshold_score = alpha * sqrt(variance(line.error, avg(line.error))) + beta * avg(line.error);


		int count_zero = 0;
		for (size_t i = 0; i < line.error.size(); i++)
		{
			if (line.error[i] == 0)
			{
				count_zero++;
			}
		}

		std::cout << "Percentage zero: " << (cbl::real) count_zero / (cbl::real) line.error.size() << std::endl;

		threshold_score *= ((cbl::real) 1.0 - (cbl::real) count_zero / (cbl::real) line.error.size());

		return threshold_score;
	}

	mrc cylinderCutOut(mrc &map, pdb &structure)
	{
		//Chop out the density around helixes using a cylinder of 5-6 angstroms
		cbl::real cropping_dist = (cbl::real) 5;

		mrc near, far;

		std::tie(near, far) = map.cylinderCrop(structure, cropping_dist);

		near.minimize();

		// To reduce Chimera blockiness visual effect
		// near.pad(5);

		return near;
	}
}