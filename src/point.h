#pragma once

struct point
{
	double coord[3];

	point(const double x, const double y, const double z)
	{
		coord[0] = x;
		coord[1] = y;
		coord[2] = z;
	};

	point()
	{
		for(int i=0; i<3; i++)
		{
			coord[i] = 0;
		}
	};

	double & operator[] (const int i) 
	{
		return coord[i];
	}

	const double & operator[] (const int i) const
	{
		return coord[i];
	}

	double x() const
	{
		return coord[0];
	};

	double y() const
	{
		return coord[1];
	}

	double z() const
	{
		return coord[2];
	}

	point operator+ (point p) const
	{
	  for(int i=0; i<3; i++)
	  {
		  p.coord[i] += coord[i];
	  }
	  return (p);
	}

	point operator* (double d) const
	{
		point p;
	  for(int i=0; i<3; i++)
	  {
		  p.coord[i] = coord[i]*d;
	  }
	  return (p);
	}

	point & operator+= (point p) 
	{
	  for(int i=0; i<3; i++)
	  {
		  coord[i] += p[i];
	  }
	  return *this;
	}

	point operator- (point p) const
	{
	  for(int i=0; i<3; i++)
	  {
		  p.coord[i] = coord[i] - p.coord[i];
	  }
	  return (p);
	}

	point & operator*= (const double s) 
	{
	  for(int i=0; i<3; i++)
	  {
		  coord[i] *= s;
	  }
	  return *this;
	}

	point & operator/= (const double s) 
	{
	  for(int i=0; i<3; i++)
	  {
		  coord[i] /= s;
	  }
	  return *this;
	}

	double norm()
	{
		double norm2 = 0;
		for(int i=0; i<3; i++)
		{
			norm2 += coord[i]*coord[i];
		}
		return pow(norm2, 0.5);
	}

	double dot(const point p) const
	{
		double d = 0;
		for(int i=0; i<3; i++)
		{
			d += coord[i]*p.coord[i];
		}
		return d;
	}

	static point cross(const point & p1, const point & p2)
	{
		point c;
		c[0] = p1.y()*p2.z()-p1.z()*p2.y();
		c[1] = p1.z()*p2.x()-p1.x()*p2.z();
		c[2] = p1.x()*p2.y()-p1.y()*p2.x();
		return c;
	}
};