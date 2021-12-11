#include <math.h>
#include <stdio.h>
#include "Matrix.h"


Matrix::Matrix(int N)
{
	this->A = new double* [N];
	this->size = N;

	for (int i = 0; i < N; ++i) this->A[i] = new double[N];
	
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			//if (i == j) this->A[i][j] = 1;
			this->A[i][j] = 0;
		}
	}
}

Matrix::Matrix(const Matrix& M)
{
	this->size = M.size;
	this->A = new double* [this->size];

	for (int i = 0; i < this->size; ++i)
	{
		this->A[i] = new double[this->size];
	}

	for (int i = 0; i < this->size; ++i)
	{
		for (int j = 0; j < this->size; ++j)
		{
			this->A[i][j] = M.A[i][j];
		}
	}
}

Matrix::~Matrix()
{
	for (int i = 0; i < this->size; ++i) delete[] this->A[i];

	delete[] this->A;
}

void Matrix::print()
{
	for (int i = 0; i < this->size; ++i)
	{
		for (int j = 0; j < this->size; ++j)
		{
			printf_s("%.3f ", this->A[i][j]);
		}

		printf_s("\n");
	}

	printf_s("\n");
}

Matrix& Matrix::operator=(const Matrix& M)
{	
	for (int i = 0; i < this->size; ++i)
	{
		for (int j = 0; j < this->size; ++j)
		{
			this->A[i][j] = M.A[i][j];
		}
	}

	return *this;
}

double* Matrix::operator*(const double* v)
{
	int N = this->size;
	double* u = new double[N];

	for (int i = 0; i < N; ++i)
	{
		double S = 0;

		for (int j = 0; j < N; ++j) S += this->A[i][j] * v[j];

		u[i] = S;
	}

	return u;
}