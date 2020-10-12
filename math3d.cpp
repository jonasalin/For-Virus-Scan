#include <math.h>
#include ".//math3d.h"
#include ".//sincos_dbl.h"

Matrix::Matrix(const Matrix &M)
{
	int r, c;
	ncols = M.ncols;
	nrows = M.nrows;
	size = ncols*nrows;
	data = new double[ncols*nrows];
	
	for(r = 0; r < nrows; r++) 
		for (c = 0; c < ncols; c++) 
			data[r*ncols + c] = M(r,c);
}

Matrix::Matrix(int sizeXY)
{
	size = sizeXY*sizeXY;
	ncols = sizeXY;
	nrows = sizeXY;
	data = new double[sizeXY*sizeXY];
	
	int r, c;
	for (r = 0; r < nrows; r++) 
		for (c = 0; c < ncols; c++) 
			data[r*ncols + c] = 0;
}

Matrix::Matrix(int sX, int sY)
{
	size = sX*sY;
	nrows = sX;
	ncols = sY;
	data = new double[sX*sY];
	
	int r, c;
	for (r = 0; r < nrows; r++) 
		for (c = 0; c < ncols; c++) 
			data[r*ncols + c] = 0;
}

void Matrix::ZeroMatrix()
{
	int r, c;
	for (r = 0; r < nrows; r++) 
		for (c = 0; c < ncols; c++) 
			data[r*ncols + c] = 0;
}

void Matrix::IdMatrix()
{
	int r, c;
	if (nrows != ncols)
		return;
	
	for (r = 0; r < nrows; r++) 
		for (c = 0; c < ncols; c++) 
			data[r*ncols + c] = 0;
	for(r = 0; r < nrows; r++)
		data[r*ncols + r] = 1;
}

void Matrix::insert(double nr, int row, int col)
{
	if(row >= nrows || row < 0 || col >= ncols || col < 0) 
	{
		cerr << "Insertion of numbers outside of a matrix boundary";
		exit(EXIT_FAILURE);
	}
	data[row*ncols + col] = nr;
}

void Matrix::incr(double nr, int row, int col)
{
	if(row >= nrows || row < 0 || col >= ncols || col < 0) 
	{
		cerr << "Increase of number outside of a matrix boundary";
		exit(EXIT_FAILURE);
	}
	data[row*ncols + col] += nr;
}

double Matrix::operator () (int row, int col) const
{
	if(row >= nrows || row < 0 || col >= ncols || col < 0) 
	{
		cerr << "Access of numbers outside of a matrix boundary";
		exit(EXIT_FAILURE);
	}
	return data[row*ncols + col];
}

void Matrix::insert(double *numbers)
{
	int r, c, p = 0;
	for (r = 0; r < nrows; r++) 
		for (c = 0; c < ncols; c++) 
			data[r*ncols + c] = numbers[p++];
}

Matrix Matrix::operator * (const Matrix &a) const
{
	double sum = 0, g, h;
	int r, c, p;
	
	if(a.nrows != ncols) 
	{
		cerr << "Matrix multiplication with non-matching number of rows of B to number of columns of A.";
		exit(EXIT_FAILURE);
	}
	Matrix M(nrows, a.ncols);
	
	for(r = 0; r < nrows; r++) 
		for(c = 0; c < a.ncols; c++) 
		{
			for(p = 0; p < ncols; p++) 
				sum += (*this)(r,p) * a(p,c);
			M.insert(sum,r,c);
			sum = 0;
		}
		
	return M;
}

Matrix& Matrix::operator *= (const Matrix &a)
{
	(*this) = (*this) * a;
	return *this;
}


Matrix Matrix::operator * (const double n) const
{
	Matrix M(nrows, ncols);
	int r, c;
	for(r = 0; r < nrows; r++) 
		for (c = 0; c < ncols; c++) 
			M.insert(data[r*ncols + c] * n, r, c);
	return M;
}

Matrix& Matrix::operator *= (const double n)
{
	(*this) = (*this) * n;
	return *this;
}

Vect_3d Matrix::operator * (const Vect_3d &a) const
{
	Vect_3d b((*this)(0,0)*a(0) + (*this)(0,1)*a(1) + (*this)(0,2)*a(2), 
			  (*this)(1,0)*a(0) + (*this)(1,1)*a(1) + (*this)(1,2)*a(2), 
			  (*this)(2,0)*a(0) + (*this)(2,1)*a(1) + (*this)(2,2)*a(2));
	return b;
}

Matrix Matrix::operator + (const Matrix &a) const
{
	Matrix b(nrows, ncols);
	if(nrows != a.nrows || ncols != a.ncols)
	{
		cerr << "Non-matching number of rows or columns in matrix addition." << endl;
		exit(EXIT_FAILURE);
	}
	for(int r = 0; r < nrows; r++)
		for(int c = 0; c < ncols; c++)
			b.insert((*this)(r,c) + a(r,c), r,c);

	return b;
}

Matrix Matrix::operator - (const Matrix &a) const
{
	Matrix b(nrows, ncols);
	if(nrows != a.nrows || ncols != a.ncols)
	{
		cerr << "Non-matching number of rows or columns in matrix subtraction." << endl;
		exit(EXIT_FAILURE);
	}
	for(int r = 0; r < nrows; r++)
		for(int c = 0; c < ncols; c++)
			b.insert((*this)(r,c) - a(r,c), r,c);

	return b;
}

const Matrix& Matrix::operator = (const Matrix &a)
{
	if(this != &a)
	{
		delete[] data;
		ncols = a.ncols;
		nrows = a.nrows;
		size = a.size;
		data = new double[size];
		
		for(int r = 0; r < a.nrows; r++) 
			for (int c = 0; c < a.ncols; c++)
				data[r*ncols + c] = a(r,c);
	}
	return *this;
}

void scale_matr(Matrix &M, double x, double y, double z)
{
	M.ZeroMatrix();
	M.insert(x,0,0);
	M.insert(y,1,1);
	M.insert(z,2,2);
}

void rot_matr(Matrix &M, Axis ax, int angl)
{
	double cosv, sinv;

	if(angl < 0)
		angl = 360 - (-angl) % 360;
	else if(angl >= 360)
		angl = angl % 360;
	
	sinv = sintabl[angl];
	cosv = costabl[angl];
	
	switch (ax)
	{
		case x:
			M.insert(1, 0, 0);
			M.insert(cosv, 1, 1);
			M.insert(cosv, 2, 2);
			M.insert(-sinv, 1, 2);
			M.insert(sinv, 2, 1);
			break;
		case y:
			M.insert(1, 1, 1);
			M.insert(cosv, 0, 0);
			M.insert(-sinv, 2, 0);
			M.insert(sinv, 0, 2);
			M.insert(cosv, 2, 2);
			break;
		case z:
			M.insert(1, 2, 2);
			M.insert(cosv, 0, 0);
			M.insert(-sinv, 0, 1);
			M.insert(sinv, 1, 0);
			M.insert(cosv, 1, 1);
			break;
	}
}

void rot_matrxyz(Matrix &M, int angl, double x, double y, double z)
{
	double cosv, sinv;

	if(angl < 0)
		angl = 360 - (-angl) % 360;
	else if(angl >= 360)
		angl = angl % 360;

	sinv = sintabl[angl];
	cosv = costabl[angl];
	
	double 	a00 = cosv + x*x*(1 - cosv), 	a01 = x*y*(1 - cosv) - z*sinv, 	a02 = x*z*(1 - cosv) + y*sinv,
			a10 = x*y*(1 - cosv) + z*sinv, 	a11 = cosv + y*y*(1 - cosv), 	a12 = y*z*(1 - cosv) + x*sinv,
			a20 = x*z*(1 - cosv) - y*sinv, 	a21 = y*z*(1 - cosv) + x*sinv, 	a22 = cosv + z*z*(1 - cosv);
	
	M.ZeroMatrix();
	M.insert(a00, 0, 0);
	M.insert(a10, 1, 0);
	M.insert(a20, 2, 0);

	M.insert(a01, 0, 1);
	M.insert(a11, 1, 1);
	M.insert(a21, 2, 1);
	
	M.insert(a02, 0, 2);
	M.insert(a12, 1, 2);
	M.insert(a22, 2, 2);
}

double Vect_3d::dot(const Vect_3d &v)
{
	return x*v.x + y*v.y + z*v.z;
}

Vect_3d Vect_3d::cross(const Vect_3d &v)
{
	Vect_3d c( y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x );
	return c;
}

Vect_3d Vect_3d::operator + (const Vect_3d& v) const
{
	Vect_3d u( x + v.x, y + v.y, z + v.z );
	return u;
}

Vect_3d Vect_3d::operator - (const Vect_3d& v) const
{
	Vect_3d u( x - v.x, y - v.y, z - v.z );
	return u;
}

Vect_3d Vect_3d::operator * (const double n) const
{
	Vect_3d u(x * n, y * n, z * n);
	return u;
}

Vect_3d Vect_3d::operator * (const Vect_3d& v) const
{
	Vect_3d u(x * v.x, y * v.y, z * v.z);
	return u;
}

Vect_3d Vect_3d::operator / (const Vect_3d& v) const
{
	Vect_3d u(x / v.x, y / v.y, z / v.z);
	return u;
}

Vect_3d Vect_3d::operator / (const double n) const
{
	Vect_3d u(x / n, y / n, z / n);
	return u;
}

const Vect_3d& Vect_3d::operator += (const Vect_3d& v)
{
	x += v.x;
	y += v.y;
	z += v.z;
	return *this;
}

const Vect_3d& Vect_3d::operator -= (const Vect_3d& v)
{
	x -= v.x;
	y -= v.y;
	z -= v.z;
	return *this;
}

double Vect_3d::operator() (const int n) const
{
	switch (n)
	{
		case 0:
			return x;
		case 1:
			return y;
		case 2:
			return z;
		default:
			return 0;
	}
}

void Vect_3d::ins(int pos, double num)
{
	switch (pos)
	{
		case 0:
			x = num;
			break;
		case 1:
			y = num;
			break;
		case 2:
			z = num;
			break;
	}
}

void transform(Vect_3d &in, Vect_3d &out, Matrix &m)
{
	out = m * in;
}

void transform(Vect_3d &in, Vect_3d &out, Matrix &m, Vect_3d &zero)
{
	Vect_3d v = in - zero;
	out = m * v;
	out = out + zero;
}

point_2d project(Vect_3d &v, int dist = DISTANCE_TO_SCREEN)
{
	if(v(2) == 0)
	{
		cerr << "Division by zero in 'point_2d project(Vect_3d &)' function." << endl;
		exit(EXIT_FAILURE);
	}
	point_2d s(v(0) * dist / v(2), v(1) * dist / v(2));
	return s;
}

point_2d point_2d::operator + (const point_2d& p) const
{
	point_2d temp(*this);
	temp += p;
	return temp;
}

const point_2d& point_2d::operator += (const point_2d& p)
{
	x += p.x;
	y += p.y;
	return *this;
}

point_2d point_2d::operator - (const point_2d& p) const
{
	point_2d temp(*this);
	temp -= p;
	return temp;
}

const point_2d& point_2d::operator -= (const point_2d& p)
{
	x -= p.x;
	y -= p.y;
	return *this;
}

Matrix powerm(Matrix &M, int p)
{
	if(M.getNCols() != M.getNRows())
	{
		cerr << "Matrix raised to p must have equal number of rows and columns";
		exit(EXIT_FAILURE);
	}
	if (p < 0)
	{
		cerr << "Exponent p in powerm(M,p) function must be >= 0";
		exit(EXIT_FAILURE);
	}
	if (p == 0)
	{
		Matrix I(M.getNRows(), M.getNRows());
		I.IdMatrix();
		return I;
	}
	if (p == 1)
		return M;
	Matrix M1 = M;
	for (int i = 2; i <= p; i++)
		M1 = M1 * M;
	return M1;		
}