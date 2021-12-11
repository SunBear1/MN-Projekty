#include <stdio.h>
#include <math.h>
#include <string>
#include <fstream>
#include <chrono>
#include "Matrix.h"

#define LAGRANGE 0
#define SPLINES 1
#define SAMPLES 480		//number of reference samples, from which interpolation nodes will be selected

using namespace std;


struct sample_t
{
	double x;			//distance of the route given in metres
	double y;			//elevation given in metres above sea level
};


void pivoting(Matrix& U, Matrix& L, Matrix& P, int i)
{
	/* finding pivot */
	double pivot = abs(U.A[i][i]);
	int pivot_index = i;

	for (int j = i + 1; j < U.size; ++j)
	{
		if (abs(U.A[j][i]) > pivot)
		{
			pivot = abs(U.A[j][i]);
			pivot_index = j;
		}
	}

	if (!U.A[pivot_index][i])
	{
		printf_s("ERROR: Matrix is singular.");
		return;
	}

	/* interchanging rows */
	if (pivot_index != i)
	{
		double tmp;

		for (int j = 0; j < U.size; ++j)
		{
			if (j >= i)
			{
				tmp = U.A[i][j];
				U.A[i][j] = U.A[pivot_index][j];
				U.A[pivot_index][j] = tmp;
			}
			else
			{
				tmp = L.A[i][j];
				L.A[i][j] = L.A[pivot_index][j];
				L.A[pivot_index][j] = tmp;
			}

			tmp = P.A[i][j];
			P.A[i][j] = P.A[pivot_index][j];
			P.A[pivot_index][j] = tmp;
		}
	}
}

void LU_decomposition(Matrix* A, double* b, double* x, int N)
{
	Matrix* U = new Matrix(*A);		//upper triangular matrix
	Matrix* L = new Matrix(N);		//lower triangular matrix
	Matrix* P = new Matrix(N);		//permutation matrix

	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j) 
		{
			if (i == j) 
				L->A[i][j] = P->A[i][j] = 1;
		}
	}

	/* creating matrices L and U such that A = L * U (with partial pivoting) */
	for (int i = 0; i < N - 1; ++i)		//columns
	{			
		pivoting(*U, *L, *P, i);
		
		for (int j = i + 1; j < N; ++j)	//rows
		{
			L->A[j][i] = U->A[j][i] / U->A[i][i];

			for (int k = i; k < N; ++k)
				U->A[j][k] = U->A[j][k] - L->A[j][i] * U->A[i][k];
		}
	}
	/* now the main equation to solve may be defined as L * U * x = b  */

	/* solving equation L * y = b for y using forward substitution method */
	b = (*P) * b;
	double* y = new double[N];

	for (int i = 0; i < N; ++i)
	{
		double S = 0;

		for (int j = 0; j < i; ++j) S += L->A[i][j] * y[j];

		y[i] = (b[i] - S) / L->A[i][i];
	}

	/* solving equation U * x = y using backward substitution method */
	for (int i = N - 1; i >= 0; --i)
	{
		double S = 0;

		for (int j = i + 1; j < N; ++j) 
			S += U->A[i][j] * x[j];

		x[i] = (y[i] - S) / U->A[i][i];
	}

	delete U;
	delete L;
	delete P;
}

double Lagrange_polynomial(const sample_t* samples, double distance, const int n)
{	
	double elevation = 0.0;

	for (int i = 0; i < n; ++i)
	{
		double a = 1.0;

		for (int j = 0; j < n; ++j)
		{
			if (i != j) a *= (distance - samples[j].x) / (samples[i].x - samples[j].x);
		}

		elevation += a * samples[i].y;
	}

	return elevation;
}

double splines(const sample_t* samples, double distance, const int nodes_number)
{
	int N = 4 * (nodes_number - 1);		//(n - 1) equations, each having 4 summands
	Matrix* M = new Matrix(N);
	double* b = new double[N];
	double* x = new double[N];			//equations' coefficients � 4 coefficients (a, b, c, d) for one equation

	for (int i = 0; i < N; ++i)
	{
		b[i] = 0;
		x[i] = 1;
	}

	/* building system of equations */
	/* S0(x0) = f(x0) = [a0 * 1] + [b0 * 0] + [c0 * 0] + [d0 * 0] + [a1 * 0] + ...  -->  [a0 * 1] = f(x0) */
	M->A[0][0] = 1;
	b[0] = samples[0].y;

	double h;

	/* S0(x1) = f(x1)  -->  [a0 * 1] + [b0 * h] + [c0 * h^2] + [d0 * h^3] = f(x1) */
	h = samples[1].x - samples[0].x;
	M->A[1][0] = 1;
	M->A[1][1] = h;
	M->A[1][2] = pow(h, 2);
	M->A[1][3] = pow(h, 3);
	b[1] = samples[1].y;

	/* S''0(x0) = 0  -->  [c0 * 1] = 0 */
	M->A[2][2] = 1;
	b[2] = 0;

	/* S''n-1(xn) = 0  -->  [cn-1 * 2] + [dn-1 * 6 * h] = 0 */
	h = samples[nodes_number - 1].x - samples[nodes_number - 2].x;
	M->A[3][4 * (nodes_number - 2) + 2] = 2;
	M->A[3][4 * (nodes_number - 2) + 3] = 6 * h;
	b[3] = 0;

	for (int i = 1; i < nodes_number - 1; ++i)
	{
		h = samples[i].x - samples[i - 1].x;

		/* Si(xi) = f(xi)  -->  [ai * 1] = f(xi) */
		M->A[4 * i][4 * i] = 1;
		b[4 * i] = samples[i].y;

		/* Si(xi+1) = f(x+1)  -->  [ai * 1] + [bi * h] + [ci * h^2] + [di * h^3] = f(xi+1) */
		M->A[4 * i + 1][4 * i] = 1;
		M->A[4 * i + 1][4 * i + 1] = h;
		M->A[4 * i + 1][4 * i + 2] = pow(h, 2);
		M->A[4 * i + 1][4 * i + 3] = pow(h, 3);
		b[4 * i + 1] = samples[i + 1].y;

		/* S'i-1(xi) = S'i(xi)  -->  [bi-1 * 1] + [ci-1 * 2 * h] + [di-1 * 3 * h^2] = [bi]  -->  [bi-1 * 1] + [ci-1 * 2 * h] + [di-1 * 3 * h^2] + [bi * (-1)] = 0 */
		M->A[4 * i + 2][4 * (i - 1) + 1] = 1;
		M->A[4 * i + 2][4 * (i - 1) + 2] = 2 * h;
		M->A[4 * i + 2][4 * (i - 1) + 3] = 3 * pow(h, 2);
		M->A[4 * i + 2][4 * i + 1] = -1;
		b[4 * i + 2] = 0;

		/* S''i-1(xi) = S''i(xi)  -->  [ci-1 * 2] + [di-1 * 6 * h] = [ci * 2]  -->  [ci-1 * 2] + [di-1 * 6 * h] + [ci * (-2)] = 0 */
		M->A[4 * i + 3][4 * (i - 1) + 2] = 2;
		M->A[4 * i + 3][4 * (i - 1) + 3] = 6 * h;
		M->A[4 * i + 3][4 * i + 2] = -2;
		b[4 * i + 3] = 0;
	}

	LU_decomposition(M, b, x, N);

	double elevation;
	for (int i = 0; i < nodes_number - 1; ++i)
	{
		elevation = 0;
		
		if (distance >= samples[i].x && distance <= samples[i + 1].x)
		{
			for (int j = 0; j < 4; ++j)
			{
				double h = distance - samples[i].x;
				
				elevation += x[4 * i + j] * pow(h, j);
			}

			break;
		}
	}

	delete[] x;
	delete[] b;
	delete M;
	
	return elevation;
}

void interpolation(const sample_t* samples, int mode, const int nodes_number, sample_t* nodes, const char* path, double* duration)
{
	string file_path = path;
	if (mode == LAGRANGE) file_path.replace(3, 4, "results\\Lagrange");
	else if (mode == SPLINES) file_path.replace(3, 4, "results\\splines");
	file_path.replace(file_path.end() - 4, file_path.end(), "_");
	file_path.append(to_string(nodes_number));
	file_path.append(".txt");

	ofstream file;
	file.open(file_path, file.out);

	if (!file.is_open())
	{
		printf_s("Failed to save the results.");
		return;
	}

	double result;
	auto start = chrono::high_resolution_clock::now();
	
	for (double i = nodes[0].x; i <= nodes[nodes_number - 1].x; i += 8)
	{
		bool interpolate = true;

		for (int j = 0; j < nodes_number; ++j)
		{
			/* values at given nodes don't need to be interpolated */
			if ((int)nodes[j].x == i)
			{
				interpolate = false;
				result = nodes[j].y;
				break;
			}
		}
		
		if (interpolate)
		{
			if (mode == LAGRANGE) result = Lagrange_polynomial(nodes, i, nodes_number);
			else if (mode == SPLINES) result = splines(nodes, i, nodes_number);
		}

		/* saving results of interpolation */
		file << i << " " << result << endl;
	}

	file.close();

	auto end = chrono::high_resolution_clock::now();
	auto difference = end - start;
	*duration += chrono::duration<double, milli>(difference).count();

	file_path.replace(file_path.end() - 4, file_path.end(), "_nodes.txt");

	file.open(file_path, file.out);

	if (!file.is_open())
	{
		printf_s("Failed to save the results.");
		return;
	}

	/* saving interpolation nodes */
	for (int i = 0; i < nodes_number; ++i) file << nodes[i].x << " " << nodes[i].y << endl;

	file.close();
}

bool read_data(const char* filename, sample_t* samples)
{
	ifstream file;
	file.open(filename, file.in);

	if (!file.is_open())
	{
		printf_s("Failed to open the specified file.");
		return false;
	}

	string sample;

	/* reading data */
	for (int i = 0; i < SAMPLES; ++i)
	{
		getline(file, sample, ' ');
		samples[i].x = stod(sample);
		getline(file, sample);
		samples[i].y = stod(sample);
	}

	file.close();

	return true;
}


int main()
{
	const char* data_sets[4] = { "gdansk.txt", "genoa_rapallo.txt", "grand_canyon.txt", "mount_everest.txt" };
	sample_t* samples = new sample_t[SAMPLES];				//reference samples
	int intervals[4] = { 16, 40, 48, 80 };					//number of reference samples that will be skipped while selecting interpolation nodes
	double average_duration[4] = { 0 };

	for (int i = 0; i < 4; ++i)
	{
		const char* filename = data_sets[i];				//name of file containing data
		const char* directory_path = "..\\data\\";
		string file_path = directory_path;
		file_path.append(filename);
		filename = file_path.c_str();

		if (!read_data(filename, samples)) return -1;

		for (int j = 0; j < 4; ++j)
		{
			int step = intervals[j];

			const int nodes_number = SAMPLES / step;		//number of nodes used to interpolate the elevation profile
			sample_t* nodes = new sample_t[nodes_number];

			for (int k = 0, l = 0; l < nodes_number; k += step, ++l)
			{
				nodes[l].x = samples[k].x;
				nodes[l].y = samples[k].y;
			}

			interpolation(samples, SPLINES, nodes_number, nodes, filename, &average_duration[j]);
		}
	}

	printf_s("\n average duration (splines method):\n");
	for (int i = 0; i < 4; ++i)
	{
		int nodes_number = SAMPLES / intervals[i];
		double duration = average_duration[i] / 4;
		printf_s("\t%d nodes: %f ms\n", nodes_number, duration);
	}


	return 0;
}