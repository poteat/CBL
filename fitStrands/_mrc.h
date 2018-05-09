#include <fstream>
#include <algorithm>
#include <vector>
#include <limits>
#include <list>

#include <iostream>
#include <iomanip>

#include "Eigen/Dense"
#include "Eigen/Geometry"

#define Infinity std::numeric_limits<float>::infinity()

struct mrc
{
	// 3-dimensional size of mrc (box)
	int nx;
	int ny;
	int nz;

	// 3d array-based voxel data
	float ***map;
	float scale;

	// Matrix form of voxel information
	Eigen::Matrix3Xf points;
	Eigen::VectorXf densities;

	// Plane rotation data (So we can undo it later)
	Eigen::Vector3f mean_vector;
	Eigen::Vector3f rotation_vector;
	float rotation_angle;

	// Raw header data (directly from mrc file format)
	struct header           // Header is exactly 1kb (1024 bytes)
	{
		int nx;             /* # of columns ( fastest changing in the map    */
		int ny;             /* # of rows                                     */
		int nz;             /* # of sections (slowest changing in the map    */

		int mode;           // 2 (floats)

		int nxstart;        // 0
		int nystart;        // 0
		int nzstart;        // 0

		int mx;             // Same as nx
		int my;             // Same as ny
		int mz;             // Same as nz

		float xlength;      // Width of entire block in angstroms
		float ylength;      // Height of entire block in angstroms
		float zlength;      // Depth of entire block in angstroms

		float alpha;        // 90.0
		float beta;         // 90.0
		float gamma;        // 90.0

		int mapc;           // 1
		int mapr;           // 2
		int maps;           // 3

		float amin;         /* minimum density value                         */
		float amax;         /* maximum density value                         */
		float amean;        /* mean density value                            */

		int ispg;           // 0
		int nsymbt;         // 0

		int extra[25];

		float xorigin;      // 0
		float yorigin;      // 0
		float zorigin;      // 0

		char map[4];        // String "MAP"

		int machineStamp;

		float rms;

		int nlabl;
		char label[10][80];
	} header;

	mrc() {};

	mrc(std::string file_name)
	{
		std::ifstream file(file_name, std::ios::binary);
		assert(file && "Could not open .mrc file");

		file.read((char*)&header, 1024);

		nx = header.nx;
		ny = header.ny;
		nz = header.nz;

		this->scale = header.xlength / (float)nx;

		// Create mrc float cube such that the memory space is contiguous

		createCube<float>(nx, ny, nz, map);

		// Read in float densities sequentially

		for (int k = 0; k < nz; k++)
		{
			for (int j = 0; j < ny; j++)
			{
				for (int i = 0; i < nx; i++)
				{
					file.read((char*)&(map[i][j][k]), 4);
				}
			}
		}

		file.close();
	};

	~mrc()
	{
		deleteCube<float>(map);
	};

	template <class T> void createCube(int nx, int ny, int nz, T***& cube)
	{
		cube = new T**[nx];
		cube[0] = new T*[nx*ny];
		cube[0][0] = new T[nx*ny*nz];

		for (int i = 1; i < nx; i++)
		{
			cube[i] = cube[i - 1] + ny;
		}

		for (int j = 1; j < nx*ny; j++)
		{
			cube[0][j] = cube[0][j - 1] + nz;
		}
	};

	template <class T> void deleteCube(T*** cube)
	{
		delete[] cube[0][0];
		delete[] cube[0];
		delete[] cube;
	};

	void write(std::string file_name)
	{
		std::ofstream file(file_name, std::ios::binary);

		file.write((char*)&header, 1024);

		for (int k = 0; k < nz; k++)
		{
			for (int j = 0; j < ny; j++)
			{
				for (int i = 0; i < nx; i++)
				{
					file.write((char*)&map[i][j][k], 4);
				}
			}
		}

		file.close();
	};

	float meanDensity()
	{
		float mean = 0;
		int count = 0;

		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				for (int k = 0; k < nz; k++)
				{
					float d = map[i][j][k];

					if (d > 0)
					{
						mean += map[i][j][k];
						count++;
					}
				}
			}
		}

		mean = count ? mean / count : 0;

		return mean;
	};

	float deviation()
	{
		float mean = meanDensity();
		int count = 0;

		float dev = 0;

		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				for (int k = 0; k < nz; k++)
				{
					float d = map[i][j][k];

					if (d > 0)
					{
						dev += std::pow(d - mean, 2);
						count++;
					}
				}
			}
		}

		dev = std::sqrt(dev / count);

		return dev;
	};

	void setThreshold(float threshold)
	{
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				for (int k = 0; k < nz; k++)
				{
					float &d = map[i][j][k];

					d = d > threshold ? d : 0;
				}
			}
		}
	};

	void convertToPoints()
	{
		int count = 0;

		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				for (int k = 0; k < nz; k++)
				{
					float d = map[i][j][k];

					if (d > 0)
					{
						count++;
					}
				}
			}
		}

		points.resize(3, count);
		densities.resize(count);

		int col = 0;

		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				for (int k = 0; k < nz; k++)
				{
					float d = map[i][j][k];

					if (d > 0)
					{
						points(0, col) = (float)i;
						points(1, col) = (float)j;
						points(2, col) = (float)k;

						points.col(col) *= scale;

						densities(col) = d;

						col++;
					}
				}
			}
		}
	};

	void normalizeOrientation()
	{
		mean_vector = points.rowwise().mean();
		points = points.colwise() - mean_vector;

		int opt = Eigen::ComputeThinU | Eigen::ComputeThinV;
		Eigen::JacobiSVD<Eigen::Matrix3Xf> svd = points.jacobiSvd(opt);

		Eigen::Vector3f normal_vector = svd.matrixU().col(1);

		float A = normal_vector(0);
		float B = normal_vector(1);
		float C = normal_vector(2);

		float mag = sqrt(A*A + C * C);

		float ux = C / mag;
		float uy = 0.0;
		float uz = -A / mag;
		float rot_theta = acos(B);

		rotation_vector(0) = ux;
		rotation_vector(1) = uy;
		rotation_vector(2) = uz;
		rotation_angle = rot_theta;

		Eigen::AngleAxis<float> aa(rotation_angle, rotation_vector);
		Eigen::Quaternion<float> q(aa);
		Eigen::Matrix3f rotation_matrix;
		rotation_matrix = q;

		auto tmp = points.transpose() * rotation_matrix;
		points = tmp.transpose();
	};

	void writePointsToPDB(std::string file_name)
	{
		std::ofstream file(file_name);

		// Denormalize with respect to mean
		points = points.colwise() + mean_vector;

		for (int i = 0; i < points.cols(); i++)
		{
			float x = points(0, i);
			float y = points(1, i);
			float z = points(2, i);

			x += header.xorigin;
			y += header.yorigin;
			z += header.zorigin;

			file << "ATOM " << std::setw(6) << std::right << i << \
				"  C   HOH A   1    " << std::setprecision(5) \
				<< std::setw(7) << std::right << x << " " \
				<< std::setw(7) << std::right << y << " " \
				<< std::setw(7) << std::right << z << std::endl;
		}

		file.close();
	};

	Eigen::Matrix<float, 3, 4> getXYBoundaryPoints()
	{
		// Return corner points in draw-coord order (Z)
		// x_min x_max x_min x_max
		// y_max y_max y_min y_min
		// 0     0     0     0

		float x_max = -Infinity;
		float x_min = Infinity;
		float y_max = -Infinity;
		float y_min = Infinity;

		for (int i = 0; i < points.cols(); i++)
		{
			float x = points(0, i);
			float y = points(1, i);
			// float z = points(2, i);

			x_max = x > x_max ? x : x_max;
			x_min = x < x_min ? x : x_min;
			y_max = y > y_max ? y : y_max;
			y_min = y < y_min ? y : y_min;
		}

		Eigen::Matrix<float, 3, 4> boundary;

		boundary << x_min, x_max, x_min, x_max,
			y_max, y_max, y_min, y_min,
			0, 0, 0, 0;

		std::cout << "Boundary points" << std::endl;
		std::cout << boundary.transpose();

		return boundary;
	};
};