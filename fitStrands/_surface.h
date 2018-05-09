
#include <assert.h>
#include <thread>
#include "_mrc.h"

#define Infinity std::numeric_limits<float>::infinity()

struct surface
{
	// Size of render_samples (Number of points to use for rendering pdb)
	const int RX = 150;
	const int RY = 150;

	// Size of projection_samples (Number of points to use for brute-force projection search)
	const int PX = 5;
	const int PY = 5;

	// Current size of patch (2x2 -> total of 4 control points)
	int X = 2;
	int Y = 2;

	// Foldedness importance multiplier
	float foldedness_multiplier = (float)100;

	// Thread info used for multithreading speedup
	std::thread::id main_thread;
	int max_hardware_threads;

	// Density map on which to build surface
	mrc *map;

	// Control points (number of rows is X * Y, initially 4)
	Eigen::Matrix3Xf control;

	// Set of points to render (surface and 3-euclidean coordinates)
	Eigen::Matrix2Xf render_samples;
	Eigen::Matrix3Xf render_points;

	// Set of points to initially find projections (brute force -> binary refinement)
	Eigen::Matrix2Xf projection_samples;
	Eigen::Matrix3Xf projection_points;

	// Cached projection results from each voxel (Size equals number of map density voxels)
	Eigen::Matrix2Xf point_projections;
	Eigen::Matrix3Xf projected_points;

	surface(mrc *map)
	{
		this->map = map;

		main_thread = std::this_thread::get_id();
		max_hardware_threads = std::thread::hardware_concurrency();

		Eigen::Matrix<float, 3, 4> corners = map->getXYBoundaryPoints();

		control.resize(3, 4);

		control << corners.col(0), corners.col(1),
			corners.col(2), corners.col(3);

		// Initialize rendering samples
		render_samples.resize(2, RX*RY);
		render_points.resize(3, RX*RY);

		// Interpolate from min to max, in case the patch is smaller than the map.
		float interp_min = 0.0;
		float interp_max = 1.0;
		int n = 0;
		for (int i = 0; i < RX; i++)
		{
			for (int j = 0; j < RY; j++)
			{
				float t = (float)i / (RX - 1);
				float u = (float)j / (RY - 1);

				t = t * (interp_max - interp_min) + interp_min;
				u = u * (interp_max - interp_min) + interp_min;

				render_samples.col(n++) << t, u;
			}
		}

		// Initialize projection samples
		projection_samples.resize(2, PX*PY);
		projection_points.resize(3, PX*PY);

		n = 0;
		for (int i = 0; i < PX; i++)
		{
			for (int j = 0; j < PY; j++)
			{
				float t = (float)i / (PX - 1);
				float u = (float)j / (PY - 1);

				projection_samples.col(n++) << t, u;
			}
		}

		// Initialize point projections
		point_projections.resize(2, map->points.cols());
		updateFullProjection();

		projected_points.resize(3, map->points.cols());
	};

	surface(Eigen::Matrix3Xf control_set)
	{
		main_thread = std::this_thread::get_id();
		max_hardware_threads = std::thread::hardware_concurrency();

		control = control_set;
	};

	void drawRenderGrid(std::string file_name)
	{
		float min_dist_limit = 1.0; // Angstroms

		refineProjection(0);	// Refine voxel projections with setback disabled

		sampleSet(render_samples, render_points);

		projected_points.resize(3, point_projections.cols());

		sampleSet(point_projections, projected_points);

		Eigen::Matrix2Xf filtered_samples;

		std::vector<int> cols_to_include;

		for (int i = 0; i < render_points.cols(); i++)
		{
			Eigen::Vector3f samp = render_points.col(i);

			// Make sure sample is within X of at least one projection

			bool valid = false;

			for (int j = 0; j < projected_points.cols(); j++)
			{
				Eigen::Vector3f proj = projected_points.col(j);

				float dist = sqrt(pow(samp(0) - proj(0), 2) + pow(samp(1) - proj(1), 2) + pow(samp(2) - proj(2), 2));

				if (dist <= min_dist_limit)
				{
					valid = true;
					break;
				}
			}

			if (valid)
			{
				cols_to_include.push_back(i);
			}
		}

		filtered_samples.resize(2, cols_to_include.size());

		for (int i = 0; i < filtered_samples.cols(); i++)
		{
			filtered_samples.col(i) = render_samples.col(cols_to_include[i]);
		}

		writePointsToPDB(file_name, filtered_samples);
	};

	void writePointsToPDB(std::string file_name, Eigen::Matrix2Xf &surface_samples)
	{
		std::ofstream file(file_name);

		Eigen::Matrix3Xf points;
		points.resize(3, surface_samples.cols());

		// Generate RX*RY point matrix
		sampleSet(surface_samples, points);

		// Unrotate
		Eigen::Vector3f rotation_vector = map->rotation_vector;
		float rotation_angle = -(map->rotation_angle);

		Eigen::AngleAxisf aa(rotation_angle, rotation_vector);
		Eigen::Quaternionf q(aa);
		Eigen::Matrix3f rotation_matrix;
		rotation_matrix = q;

		Eigen::MatrixX3f tmp = points.transpose() * rotation_matrix;
		points = tmp.transpose();

		// Denormalize with respect to mean
		points = points.colwise() + map->mean_vector;

		for (int i = 0; i < points.cols(); i++)
		{
			float x = points(0, i);
			float y = points(1, i);
			float z = points(2, i);

			x += map->header.xorigin;
			y += map->header.yorigin;
			z += map->header.zorigin;

			file << "ATOM " << std::setw(6) << std::right << i << \
				"  H   HOH A   1    " << std::setprecision(5) \
				<< std::setw(7) << std::right << x << " " \
				<< std::setw(7) << std::right << y << " " \
				<< std::setw(7) << std::right << z << std::endl;
		}

		file.close();
	};

	float binomial(float n, float k)
	{
		float prod = 1.0;
		for (int i = 1; i <= k; i++)
		{
			prod *= (float)(n + 1.0 - (float)i) / (float)i;
		}

		return prod;
	};

	float basis(float t, float i, float n)
	{
		return (float)binomial(n, i) * std::pow(t, i) * std::pow((float)1.0 - t, n - i);
	};

	float angleBetweenThreePoints(Eigen::Vector3f p1, Eigen::Vector3f p2, Eigen::Vector3f p3)
	{
		enum { X, Y, Z };

		float x1 = p2(X) - p1(X);
		float y1 = p2(Y) - p1(Y);
		float z1 = p2(Z) - p2(Z);

		float x2 = p3(X) - p2(X);
		float y2 = p3(Y) - p2(Y);
		float z2 = p3(Z) - p2(Z);

		return angleBetweenTwoVectors(x1, y1, z1, x2, y2, z2);
	};

	float angleBetweenTwoVectors(float x1, float y1, float z1, float x2, float y2, float z2)
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

		auto sign = [](float x) {return (x > 0.0) - (x < 0.0); };

		float dir = (float)sign(dot_dir);

		return dir * ang;
	};

	Eigen::Vector3f sample(float t, float u)
	{
		return sample(t, u, control);
	};

	Eigen::Vector3f sample(float t, float u, const Eigen::Matrix3Xf &control)
	{
		Eigen::Vector3f p = Eigen::Vector3f::Zero();

		for (int i = 0; i < control.cols(); i++)
		{
			Eigen::Vector3f control_point = control.col(i);

			float x = (float)(i % X);
			float y = (float)(i / Y);

			float product = basis(t, x, (float)(X - 1.0)) * basis(u, y, (float)(Y - 1));

			p += product * control_point;
		}

		return p;
	};

	void sampleSet(const Eigen::Matrix2Xf &input, Eigen::Matrix3Xf &output)
	{
		sampleSet(input, output, control);
	};

	void sampleSet(const Eigen::Matrix2Xf &input, Eigen::Matrix3Xf &output, const Eigen::Matrix3Xf &control)
	{
		for (int i = 0; i < input.cols(); i++)
		{
			float t = input(0, i);
			float u = input(1, i);

			output.col(i) << sample(t, u, control);
		}
	};

	Eigen::Vector2f desample(Eigen::Vector3f p)
	{
		float min_dist = Infinity;
		int min_point = -1;

		for (int j = 0; j < projection_points.cols(); j++)
		{
			Eigen::Vector3f proj = projection_points.col(j);

			float square_dist = std::pow(p(0) - proj(0), 2) + std::pow(p(1) - proj(1), 2) + std::pow(p(2) - proj(2), 2);

			if (square_dist < min_dist)
			{
				min_dist = square_dist;
				min_point = j;
			}
		}

		Eigen::Vector2f initial_projection = projection_samples.col(min_point);

		Eigen::Vector2f initial_delta;
		initial_delta << (float)1.0 / (float)(PX - 1), (float)1.0 / (float)(PY - 1);

		float tolerance = (float).001;

		return localBinaryOptimize(initial_projection, initial_delta, p, tolerance);
	};

	Eigen::Matrix3Xf localBinaryOptimize(const Eigen::Matrix3Xf &initial_values, float initial_delta, const std::function<float(const Eigen::Matrix3Xf&)> &scoring_function, float tolerance)
	{
		Eigen::Matrix3Xf parameter = initial_values;
		Eigen::Matrix3Xf delta = Eigen::Matrix3Xf::Constant(initial_values.rows(), initial_values.cols(), initial_delta);
		Eigen::Matrix3Xf prev_direction = Eigen::Matrix3Xf::Zero(initial_values.rows(), initial_values.cols());

		bool finished;
		float score;
		do
		{
			finished = true;
			score = scoring_function(parameter);

			Eigen::Matrix3Xf old_parameter = parameter;

			for (int i = 0; i < parameter.cols(); i++)
			{
				for (int j = 0; j < parameter.rows(); j++)
				{
					if (delta(j, i) > tolerance)
					{
						float val = parameter(j, i);
						float val_inc = val + delta(j, i);
						float val_dec = val - delta(j, i);

						Eigen::Matrix3Xf parameter_inc = parameter;
						Eigen::Matrix3Xf parameter_dec = parameter;

						parameter_inc(j, i) = val_inc;
						parameter_dec(j, i) = val_dec;

						float score_inc = scoring_function(parameter_inc);
						float score_dec = scoring_function(parameter_dec);

						if (score_inc <= score_dec)
						{
							if (score - score_inc > tolerance)
							{
								parameter = parameter_inc;
								score = score_inc;
							}
							else if (prev_direction(j, i) == 1)
							{
								delta(j, i) /= 2.0;
							}

							prev_direction(j, i) = 1;
						}
						else if (score_dec < score_inc)
						{
							if (score - score_dec > tolerance)
							{
								parameter = parameter_dec;
								score = score_dec;
							}
							else if (prev_direction(j, i) == -1)
							{
								delta(j, i) /= 2.0;
							}

							prev_direction(j, i) = -1;
						}
					}

					finished &= delta(j, i) < tolerance;
				}
			}

			if (old_parameter == parameter)
			{
				finished = true;
			}


			float setback = (float)(X >= 4 ? 0 : 0.2);

			refineProjection(setback);

			std::cout << delta.transpose() << std::endl << std::endl;

		} while (finished == false);

		return parameter;
	};

	Eigen::VectorXf localBinaryOptimize(const Eigen::VectorXf &initial_values, const Eigen::VectorXf &initial_delta, const Eigen::Vector3f &p, float tolerance)
	{
		assert(initial_values.rows() == initial_delta.rows());

		Eigen::VectorXf parameter = initial_values;
		Eigen::VectorXf delta = initial_delta;
		Eigen::VectorXi prev_direction = Eigen::VectorXi::Zero(initial_values.rows());

		bool finished;
		Eigen::Vector3f proj;
		float dx, dy, dz;
		float score;
		do
		{
			finished = true;

			proj = sample(parameter(0), parameter(1));
			dx = p(0) - proj(0);
			dy = p(1) - proj(1);
			dz = p(2) - proj(2);
			score = dx * dx + dy * dy + dz * dz;

			for (int i = 0; i < parameter.rows(); i++)
			{
				float val = parameter(i);
				float val_inc = val + delta(i);
				float val_dec = val - delta(i);

				Eigen::VectorXf parameter_inc = parameter;
				Eigen::VectorXf parameter_dec = parameter;

				parameter_inc(i) = val_inc;
				parameter_dec(i) = val_dec;

				//float score_inc = scoring_function(parameter_inc);
				proj = sample(parameter_inc(0), parameter_inc(1));
				dx = p(0) - proj(0);
				dy = p(1) - proj(1);
				dz = p(2) - proj(2);
				float score_inc = dx * dx + dy * dy + dz * dz;

				//float score_dec = scoring_function(parameter_dec);
				proj = sample(parameter_dec(0), parameter_dec(1));
				dx = p(0) - proj(0);
				dy = p(1) - proj(1);
				dz = p(2) - proj(2);
				float score_dec = dx * dx + dy * dy + dz * dz;

				if (score_inc <= score_dec)
				{
					if (score - score_inc > tolerance)
					{
						parameter = parameter_inc;
					}
					else if (prev_direction(i) == 1)
					{
						delta(i) /= 2.0;
					}

					prev_direction(i) = 1;
				}
				else if (score_dec < score_inc)
				{
					if (score - score_dec > tolerance)
					{
						parameter = parameter_dec;
					}
					else if (prev_direction(i) == -1)
					{
						delta(i) /= 2.0;
					}

					prev_direction(i) = -1;
				}

				finished &= delta(i) < tolerance;
			}
		} while (finished == false);

		return parameter;
	};

	void updateFullProjection()
	{
		sampleSet(projection_samples, projection_points);

		for (int i = 0; i < map->points.cols(); i++)
		{
			Eigen::Vector3f p = map->points.col(i);

			float min_dist = Infinity;
			int min_point = -1;

			for (int j = 0; j < projection_points.cols(); j++)
			{
				Eigen::Vector3f proj = projection_points.col(j);

				float square_dist = std::pow(p(0) - proj(0), 2) + std::pow(p(1) - proj(1), 2) + std::pow(p(2) - proj(2), 2);

				if (square_dist < min_dist)
				{
					min_dist = square_dist;
					min_point = j;
				}
			}

			Eigen::Vector2f initial_projection = projection_samples.col(min_point);

			Eigen::Vector2f initial_delta;
			initial_delta << (float)1.0 / (float)(PX - 1), (float)1.0 / (float)(PY - 1);

			float tolerance = (float).001;

			point_projections.col(i) = localBinaryOptimize(initial_projection, initial_delta, p, tolerance);
		}
	};

	void refineProjection(float setback = 0.2)
	{
		std::vector<std::thread> threads;
		std::function<void(int)> func = [this, setback](int i) {refineProjectionThread(i, setback); };

		for (int i = 1; i < max_hardware_threads; i++)
		{
			std::thread worker(func, i);
			threads.push_back(std::move(worker));
		}

		func(0);

		for (unsigned int i = 0; i < threads.size(); i++)
		{
			threads[i].join();
		}
	};

	void refineProjectionThread(int thread, float setback)
	{
		float min = setback;
		float max = (float) 1.0 - setback;

		for (int i = thread; i < map->points.cols(); i += max_hardware_threads)
		{
			Eigen::Vector3f p = map->points.col(i);
			Eigen::Vector2f initial_projection = point_projections.col(i);

			Eigen::Vector2f initial_delta;
			initial_delta << (float)1.0 / (float)(PX - 1), (float)1.0 / (float)(PY - 1);

			initial_delta /= 100;

			float tolerance = (float).001;

			Eigen::Vector2f new_projection;
			new_projection = localBinaryOptimize(initial_projection, initial_delta, p, tolerance);

			if (new_projection(0) < min)
			{
				new_projection(0) = min;
			}
			else if (new_projection(0) > max)
			{
				new_projection(0) = max;
			}

			if (new_projection(1) < min)
			{
				new_projection(1) = min;
			}
			else if (new_projection(1) > max)
			{
				new_projection(1) = max;
			}

			point_projections.col(i) = new_projection;
		}
	};

	float foldednessScore(const Eigen::Matrix3Xf &control)
	{
		int grid_point_sampling = 5;

		float delta = (float) 1.0 / (float)grid_point_sampling;

		float max_foldedness = 0;

		for (float i = 0.0; i <= 1.0; i += delta)
		{
			for (float j = 0.0; j <= 1.0; j += delta)
			{
				Eigen::Vector3f p1 = sample(i, j, control);
				Eigen::Vector3f p2 = sample(i + delta, j, control); // i+
				Eigen::Vector3f p3 = sample(i, j - delta, control); // j-
				Eigen::Vector3f p4 = sample(i - delta, j, control); // i-
				Eigen::Vector3f p5 = sample(i, j + delta, control); // j+

				float vertical_twist = angleBetweenThreePoints(p5, p2, p3);
				float horizontal_twist = angleBetweenThreePoints(p4, p1, p2);

				float foldedness = std::max(abs(vertical_twist), abs(horizontal_twist)) / 2;

				if (foldedness > max_foldedness)
				{
					max_foldedness = foldedness;
				}
			}
		}

		return (max_foldedness * max_foldedness) * foldedness_multiplier;
	};

	float calcSurfaceScore(const Eigen::Matrix3Xf &control)
	{
		sampleSet(point_projections, projected_points, control);

		float sum = 0;

		for (int i = 0; i < map->points.cols(); i++)
		{
			Eigen::Vector3f p = map->points.col(i);
			Eigen::Vector3f proj = projected_points.col(i);

			float sq_dist = std::pow(p[0] - proj[0], 2) + std::pow(p[1] - proj[1], 2) + std::pow(p[2] - proj[2], 2);
			sq_dist *= map->densities(i);
			sum += sq_dist;
		}

		return sum + foldednessScore(control);
	};

	void optimizeSurface()
	{
		std::function<float(Eigen::Matrix3Xf)> scoring_func = [this](Eigen::Matrix3Xf control) {return calcSurfaceScore(control); };
		control = localBinaryOptimize(control, (float)1.0, scoring_func, (float).1);

		std::cout << control.transpose() << std::endl << std::endl;
	};

	void elevateDegree()
	{
		Eigen::Matrix3Xf new_control;
		new_control.resize(3, (X + 1)*(Y + 1));

		float mul = (X - (float)1.0001) / (float)X;

		for (int i = 0; i < X + 1; i++)
		{
			for (int j = 0; j < Y + 1; j++)
			{
				int patch_x = (int)floor(i * mul);
				int patch_y = (int)floor(j * mul);

				float offset_x = ((float)i / (float)X - (float)patch_x / (float)(X - 1)) * (float)(X - 1);
				float offset_y = ((float)j / (float)Y - (float)patch_y / (float)(Y - 1)) * (float)(Y - 1);

				Eigen::Vector3f p1 = control.col(patch_x + patch_y * X);
				Eigen::Vector3f p2 = control.col(patch_x + 1 + patch_y * X);
				Eigen::Vector3f p3 = control.col(patch_x + (patch_y + 1) * X);
				Eigen::Vector3f p4 = control.col(patch_x + 1 + (patch_y + 1) * X);

				new_control.col(i + j * (X + 1)) = bilinearInterpolate(p1, p2, p3, p4, offset_x, offset_y);
			}
		}

		control = new_control;
		X++;
		Y++;
	};

	Eigen::Vector3f bilinearInterpolate(Eigen::Vector3f p1, Eigen::Vector3f p2, Eigen::Vector3f p3, Eigen::Vector3f p4, float t, float u)
	{
		// Use the given 3x4 matrix to construct a 2x2 bezier surface temporarily.
		Eigen::Matrix3Xf set;
		set.resize(3, 4);

		set << p1, p2, p3, p4;

		surface bilinear(set);

		return bilinear.sample(t, u);
	};
};