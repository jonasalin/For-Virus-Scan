#ifndef MATH3D_H
#define MATH3D_H
#define DISTANCE_TO_SCREEN 1024
#include <iostream.h>

enum Axis { x, y, z };

class Vect_3d
{
public:
	Vect_3d() {}; //: x(0), y(0), z(0) {};
	Vect_3d(double x1, double y1, double z1) : x(x1), y(y1), z(z1) {};
	double dot(const Vect_3d &v);
	Vect_3d cross(const Vect_3d &v);
	double length(void) { return sqrt(x*x + y*y + z*z); };
	double lengthSq(void) { return x*x + y*y + z*z; };
	void ins(int pos, double num);
	Vect_3d operator + (const Vect_3d& v) const;
	Vect_3d operator - (const Vect_3d& v) const;
	Vect_3d operator * (const double n) const;
	Vect_3d operator * (const Vect_3d& v) const;
	Vect_3d operator / (const double n) const;
	Vect_3d operator / (const Vect_3d& v) const;
	const Vect_3d& operator += (const Vect_3d& v);
	const Vect_3d& operator -= (const Vect_3d& v);
	double operator() (const int n) const;
private:
	double x, y, z;
};

class Matrix
{
public:
	Matrix(const Matrix &M);
	Matrix(int sizeXY);
	Matrix(int sX, int sY);
	~Matrix() { delete[] data; };
	void insert(double *numbers);
	void ZeroMatrix(void);
	void IdMatrix(void);
	void insert(double nr, int row, int col);
	void incr(double nr, int row, int col);
	double operator() (int row, int col) const;
	Matrix operator * (const Matrix &a) const; 	
	Matrix& operator *= (const Matrix &a);
	Matrix operator * (const double n) const;    		
	Matrix& operator *= (const double n);
	Vect_3d operator * (const Vect_3d &a) const; 	
	const Matrix& operator = (const Matrix &a);
	Matrix operator + (const Matrix &a) const;
	Matrix operator - (const Matrix &a) const;
	int getNRows(void) { return nrows; };
	int getNCols(void) { return ncols; };
	int getSize(void) { return size; };
private:
	double *data;
	int size, ncols, nrows;
};

class point_2d {
public:
	point_2d() {}; //: x(0), y(0) {};
	point_2d(int x, int y) : x(x), y(y) {};
	double operator () (const int p) const { switch (p) { case 0: return x; case 1: return y; default: return NULL;} };
	point_2d operator + (const point_2d& p) const;
	point_2d operator - (const point_2d& p) const;
	const point_2d& operator += (const point_2d& p);
	const point_2d& operator -= (const point_2d& p);
private:
	double x, y;
};

void rot_matr(Matrix &M, Axis ax, int angl);
void rot_matrxyz(Matrix &M, int angl, double x, double y, double z);
void scale_matr(Matrix &M, double x, double y, double z);
void transform(Vect_3d &in, Vect_3d &out, Matrix &m);
void transform(Vect_3d &in, Vect_3d &out, Matrix &m, Vect_3d &zero);
point_2d project(Vect_3d &v, int dist);
Matrix powerm(Matrix &M, int p);
#endif
