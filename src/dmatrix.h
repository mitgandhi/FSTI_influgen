#pragma once
#include <vector>

//This implements a "dense" style matrix and limited matrix operations using std::vector

template <class T>
class dmatrix
{
	std::vector<T> mat;

public:
	size_t m;
	size_t n;

	//Constructor
	dmatrix(const size_t i, const size_t j)
	{
		m = i;
		n = j;
		
		mat.resize(m*n, 0);
	};

	dmatrix()
	{
	};

	void resize(const size_t i, const size_t j)
	{
		//clear the current matrix
		mat.clear();

		m = i;
		n = j;
		
		mat.resize(m*n, 0);
	};

	//Direct matrix access (i,j)
	inline T & operator() (const size_t i, const size_t j) 
	{
		return mat[i*n+j];
	}

	inline const T & operator() (const size_t i, const size_t j) const
	{
		return mat[i*n+j];
	}

	//Direct matrix access [i]
	inline T & operator[] (const int i) 
	{
		return mat[i];
	}

	//Scalar multiplication
	dmatrix<T> & operator*= (const double s)
	{
		for(int i=0; i<m*n; i++)
		{
			mat[i] *= s;
		}
		return *this;
	}

	//Matrix multiplication
	dmatrix<T> operator* (const dmatrix<T> & M) const
	{
		if(n != M.m)
		{
			//Invalid matrix sizes for multiplication
			return dmatrix<T> (0,0);
		}
		dmatrix<T> r(m,M.n);

		for(int i=0; i<m; i++)
		{
			for(int j=0; j<M.n; j++)
			{
				for(int k=0; k<n; k++)
				{
					//could optimize to use direct indexing?
					r(i,j) += (*this)(i,k)*M(k,j);
				}
			}
		} 

		return r;
	}

	//Matrix addition
	dmatrix<T> & operator+= (const dmatrix<T> M)
	{
		if(m != M.m || n != M.n)
		{
			//Invalid matrix sizes for addition so don't do anything
			return *this;
		}

		for(int i=0; i<m*n; i++)
		{
			mat[i] += M(i);
		}
		return *this;
	}

	//Transpose
	dmatrix<T> t()
	{
		dmatrix<T> trans(n,m);
		for(int i=0; i<m; i++)
		{
			for(int j=0; j<n; j++)
			{
				trans(j,i) = (*this)(i,j);
			}
		}

		return trans;
	}

	//Matrix Determinate
	double det()
	{
		if(m != n)
		{
			return 0;
		}

		if(m == 1)
		{
			return mat[0];
		}

		//We will use Laplace's formula for the determinate
		double d = 0;
		for(int j=0; j<n; j++)
		{
			//Create the 'Minor' matrix
			dmatrix M(m-1, n-1);
			for(int r=1; r<m; r++)
			{
				for(int c=0; c<n; c++)
				{
					if(c < j)
					{
						M(r-1,c) = (*this)(r,c);
					} else if (c > j)
					{
						M(r-1,c-1) = (*this)(r,c);
					}
				}
			}
			
			d += pow(-1.0,j) * (*this)(0,j) * M.det();
		}
		return d;
	}

};
