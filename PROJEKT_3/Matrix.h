#pragma once


class Matrix
{
public:
	double** A;
	int size;

	Matrix(int N);
	Matrix(const Matrix& M);
	~Matrix();

	void print();

	Matrix& operator=(const Matrix& M);
	double* operator*(const double* v);
};