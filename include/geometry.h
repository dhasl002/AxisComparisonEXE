
#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <iostream>
#include <cmath>
#include <limits>
#include <math.h>
#include "utilityfunctions.h"
#include "Matrix.h"



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// List Of Data Structures //////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct Coordinate{
	double x,
		   y,
		   z;
	//Initiate values to 0.0
	Coordinate() : x(0.00),
					y(0.00),
					z(0.00) {}

    Coordinate operator= (const Coordinate& c1)
    {
        x = c1.x;
        y = c1.y;
        z = c1.z;
    }

    Coordinate operator+= (const Coordinate& c1)
    {
        x = c1.x+x;
        y = c1.y+y;
        z = c1.z+z;
    }

    //Coordinate operator- (const Coordinate& c1)
};
bool operator== (const Coordinate& c1, const Coordinate& c2)
{

    if ((fabs(c1.x-c2.x) < std::numeric_limits<double>::epsilon()) && (fabs(c1.y-c2.y) < std::numeric_limits<double>::epsilon()) && (fabs(c1.z-c2.z) < std::numeric_limits<double>::epsilon()))
        return true;
    else return false;
}
Coordinate operator- (const Coordinate& c1, const Coordinate& c2)
{
    Coordinate c3;
    c3.x = c1.x-c2.x;  c3.y = c1.y-c2.y;  c3.z = c1.z-c2.z;
    return c3;
}
Coordinate operator+ (const Coordinate& c1, const Coordinate& c2)
{
    Coordinate c3;
    c3.x = c1.x+c2.x;  c3.y = c1.y+c2.y;  c3.z = c1.z+c2.z;
    return c3;
}




// DATA STRUCTURE TO STORE PHI AND PSI ANGLES
struct Torsion
{
	double phi;
	double psi;
	//Initialize torsion angles to 999.00
	Torsion() : phi(999.00), psi(999.00) {}
};


#define PI  3.1415926535897932384626433832  // Pi

#define ABS(x)     ((x) >= 0 ? (x) : -(x))   // absolute value

#define SQR(x)   (x * x)  //square function

// convert from radian to degree
#define toDegree(radian) \
	((double) (radian * 180) / PI)

// converts from degree to radian
#define toRadian(degree) \
	((double) (degree * PI) / 180)

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// End Of Constants and Data Structures /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// List of Functions ////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double ge(Matrix A, Matrix b, Matrix & x);												//Function to perform Gaussian elimination to solve A*x = b
Coordinate pointLineIntersection(Coordinate p1,Coordinate p2,Coordinate p3);			//calculate the point where a line and a point intersect
float linePointIntersection(Coordinate lineP1, Coordinate lineP2, Coordinate p, Coordinate &intersectionPnt);	//same as "pointLineIntersection"..but it looks faster
void buildRotationMatrix (double mtx[4][4],Coordinate p1, Coordinate p, double angle);	//build rotation matrix to rotate a round line
void rotatePoint(double mtrx[4][4],Coordinate &v);										//rotate a point around axis....given is the rotation Matrix
double getDistance(Coordinate p1, Coordinate p2);										//calculate the distance b/w 2 points
Coordinate triangleCenter(Coordinate p1, Coordinate p2, Coordinate p3);					//return the center points for a triangle
Coordinate quadCenter(Coordinate p1, Coordinate p2, Coordinate p3, Coordinate p4);      //return the center points for 4 points
Coordinate quadCenter1(Coordinate p1);
Coordinate quadCenter2(Coordinate p1, Coordinate p2);
Coordinate pointOnLine(Coordinate p1, Coordinate p2, float dist);						//find a point on line P1-P2 that far from P2 dist toward P2
double getAngleDegree(Coordinate p1, Coordinate p2, Coordinate p3);						//compute the angle b/w 3 points in degree
double getAngleRadian(Coordinate p1, Coordinate p2, Coordinate p3);						//compute the angle b/w 3 points in radian
double getTorsionAngle(Coordinate p1, Coordinate p2, Coordinate p3, Coordinate p4);		//calculate the torsion angle b/w given 4 points
double getDistLines(Coordinate l1P1, Coordinate l1P2, Coordinate l2P1, Coordinate l2P2);//get the shortest distance b/w 2 lines segments
double getDistInfLinePoint(Coordinate lineP1, Coordinate lineP2, Coordinate p);			//find the shortest distance between a line (infinite line ) and a point
double getDistLineSegPoint(Coordinate lineP1, Coordinate lineP2, Coordinate p);			//find the shortest distance between a line segment (given by two points) and a point
void torsion2xyz(vector<Torsion> torsions, vector<Coordinate> &XYZlist);				//given a list of torsion angles..returns a list of coordinates for backbone atoms
void extrude (vector<Coordinate> &XYZlist,double DZ, double theta, double TAU, int &l);	//auxiliary function used by torsion2xyz function
Coordinate getPoint (Coordinate A, Coordinate B, double BC, double ABCangle);			//given 2 points in a triangle, the second side length (BC) and the angle (ABC) predict one possible C coordinate
void peakClustering(vector<Coordinate> pnts, vector<vector<Coordinate> > & clusters, Coordinate sPnt, float r);	//given a set of points, cluster radius and a starting point...sample given points into clusters
//void catmullRom(vector<Coordinate>* helixAxis, double stepSize); //given a set of points, interpolate a curved line between them /// NOT TESTED WITH vector LENGTH 2
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// ////////// End Of The List //////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////// Vectors CLASS ///////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Vectors
{
public:

	// Constructors
	Vectors();
	Vectors(Coordinate);
	Vectors(double, double, double);
	Vectors(Coordinate,Coordinate);


	void operator+=(Vectors);		//sum 2 vectors
	void operator+=(Coordinate);	//sum a point to the vector
	void operator+= (double);		//sum a number to the vector
	void operator-=(Vectors);		//subtract a vector from the vector
	void operator-=(Coordinate);	//subtract a point from the vector
	void operator-= (double);		//subtract a number from the vector
	void set(Coordinate,Coordinate);	//given 2 points...set the vector
	void set(Coordinate);				//given  a point make the vector equal this point
	void print();						//print the vector

	Vectors mul(double);				//multiply the vector by a number
	void divide(double);				//divide the vector by a number
	Vectors cross(Vectors);				//calculate the cross product "gives the normal"
	Vectors getDiff(Vectors);			//get the difference b/w 2 vectors
	Vectors unit();						//calculate the unit vector

	double dot(Vectors);				//calculate the dot product
	double getDistance(Vectors);		//calculate the distance be/w 2 vectors
	double getAngleRadian(Vectors);		//calculate the angle b/w 2 vectors in Radian
	double getAngleDegree(Vectors);		//calculate the angle b/w 2 vectors in degree
	double length();					//calculate the length of the vector


	Coordinate getCoordinates();		//convert the from Vectors class to Coordinate data Type
	double getX();						//return the value of x coordinate
	double getY();						//return the value of y coordinate
	double getZ();						//return the value of z coordinate


private:
	double  x_,
			y_,
			z_;

};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////// End of Vectors CLASS ///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////// Implementation of Vectors CLASS /////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Vectors::Vectors()		///////constructors/////
{
	x_ = 0;
	y_ = 0;
	z_ = 0;
}
Vectors::Vectors(Coordinate p)	///////constructors/////
{
	x_	= p.x;
	y_	= p.y;
	z_	= p.z;
}
Vectors::Vectors(double xPrime, double yPrime, double zPrime)		///////constructors/////
{
	x_	= xPrime;
	y_	= yPrime;
	z_	= zPrime;
}
Vectors::Vectors(Coordinate p1,Coordinate p2)			///////constructors/////
{
	x_ = p1.x - p2.x;
	y_ = p1.y - p2.y;
	z_ = p1.z - p2.z;
}
///////////////////////////////////////////////////////////////////////

double Vectors::getDistance(Vectors v)
{
	double xDist = v.x_ - x_;
	double yDist = v.y_ - y_;
	double zDist = v.z_ - z_;
	return sqrt(xDist*xDist + yDist*yDist + zDist*zDist);
}
//////////////////////////////////////////////////////////////////////////
Vectors Vectors::getDiff(Vectors v)
{
// computes the difference b/w the vector with vector v
	Vectors tmpVect;
	tmpVect.x_	= x_ - v.x_;
	tmpVect.y_	= y_ - v.y_;
	tmpVect.z_	= z_ - v.z_;
	return tmpVect;
}
///////////////////////////////////////////////////////////////////////////
void Vectors::operator+=(Vectors v)
{
	x_ += v.x_;
	y_ += v.y_;
	z_ += v.z_;
}
///////////////////////////////////////////////////////////////////////////
void Vectors::operator+=(Coordinate p)
{
	x_ += p.x;
	y_ += p.y;
	z_ += p.z;
}
///////////////////////////////////////////////////////////////////////////
void Vectors::operator+= (double num)
{
	x_ += num;
	y_ += num;
	z_ += num;
}
///////////////////////////////////////////////////////////////////////////
void Vectors::operator-=(Vectors v)
{
	x_ -= v.x_;
	y_ -= v.y_;
	z_ -= v.z_;
}
///////////////////////////////////////////////////////////////////////////
void Vectors::operator-= (Coordinate p)
{
	x_ -= p.x;
	y_ -= p.y;
	z_ -= p.z;
}
///////////////////////////////////////////////////////////////////////////
void Vectors::operator-= (double num)
{
	x_ -= num;
	y_ -= num;
	z_ -= num;
}
///////////////////////////////////////////////////////////////////////////
Vectors Vectors::mul(double num)
{
//Multiply a vector by number
	Vectors tmp;
	tmp.x_	= x_ * num;
	tmp.y_	= y_ * num;
	tmp.z_	= z_ * num;
	return tmp;
}
///////////////////////////////////////////////////////////////////////////
void Vectors::divide (double num)
{
	Vectors tmp;
	x_ /= num;
	y_ /= num;
	z_ /= num;
}
///////////////////////////////////////////////////////////////////////////
double Vectors::dot(Vectors v)
{
// A . B = Ax*Bx + Ay*By + Az*Bz
	return ((x_ * v.x_)+ (y_ * v.y_) + (z_ * v.z_));
}
///////////////////////////////////////////////////////////////////////////
Vectors Vectors::cross(Vectors v)
{
//computes the normal b/w these two vectors
//A x B = <Ay*Bz - Az*By, Az*Bx - Ax*Bz, Ax*By - Ay*Bx>
	Vectors tmp;
	tmp.x_	= y_ * v.z_ - z_ * v.y_;
	tmp.y_	= z_ * v.x_ - x_ * v.z_;
	tmp.z_	= x_ * v.y_ - y_ * v.x_;
	return tmp;
}
///////////////////////////////////////////////////////////////////////////
double Vectors::getAngleRadian(Vectors v)
{
	return acos(((x_ * v.x_ + y_ * v.y_ + z_ * v.z_)/(sqrt(x_ * x_ + y_ * y_ + z_ * z_) * sqrt(v.x_ * v.x_ + v.y_ * v.y_ + v.z_ * v.z_))));
}
///////////////////////////////////////////////////////////////////////////
double Vectors::getAngleDegree(Vectors v)
{
	double angle = acos(((x_ * v.x_ + y_ * v.y_ + z_ * v.z_)/(sqrt(x_ * x_ + y_ * y_ + z_ * z_) * sqrt(v.x_ * v.x_ + v.y_ * v.y_ + v.z_ * v.z_))));
	return ((angle * 180) / PI);
}
///////////////////////////////////////////////////////////////////////////
void Vectors::set(Coordinate p1,Coordinate p2) //p1 - p2
{
	x_ = p1.x - p2.x;
	y_ = p1.y - p2.y;
	z_ = p1.z - p2.z;
}
///////////////////////////////////////////////////////////////////////////
void Vectors::set(Coordinate p)
{
	x_	= p.x;
	y_	= p.y;
	z_	= p.z;
}
///////////////////////////////////////////////////////////////////////////
Coordinate Vectors::getCoordinates()
{
	Coordinate p;
	p.x = x_;
	p.y = y_;
	p.z = z_;
	return p;
}
///////////////////////////////////////////////////////////////////////////
double Vectors::getX()
{
	return x_;
}
///////////////////////////////////////////////////////////////////////////
double Vectors::getY()
{
	return y_;
}
///////////////////////////////////////////////////////////////////////////
double Vectors::getZ()
{
	return z_;
}
///////////////////////////////////////////////////////////////////////////
double Vectors::length()
{
	return sqrt(x_ * x_ + y_ * y_ + z_ * z_);
}
///////////////////////////////////////////////////////////////////////////
Vectors Vectors::unit()
{
	Vectors tmp;
	tmp.x_ = x_/length();
	tmp.y_ = y_/length();
	tmp.z_ = z_/length();

	return tmp;
}
///////////////////////////////////////////////////////////////////////////
void Vectors::print(){
	cout<<"X = "<<x_<<" Y = "<<y_<<" Z = "<<z_<<endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// END OF VECTORS CLASS IMPLEMENTATION///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////// Some Utility Functions /////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// source http://www.algarcia.org/nummeth/Programs2E.html
// ge - Function to perform Gaussian elimination to solve A*x = b
//      using scaled column pivoting
// Inputs
//    A    -    Matrix A (N by N)
//    b    -    Vector b (N by 1)
// Outputs
//    x    -    Vector x (N by 1)
//  determ -    Determinant of matrix A	 (return value)
double ge(Matrix A, Matrix b, Matrix & x) {

  int N = A.nRow();
  assert( N == A.nCol() && N == b.nRow() && N == x.nRow() );

  int i, j, k;
  Matrix scale(N);	// Scale factor
  int *index;  index = new int [N+1];	// Row index list

  //* Set scale factor, scale(i) = max( |A(i,j)| ), for each row
  for( i=1; i<=N; i++ ) {
    index[i] = i;		   // Initialize row index list
    double scaleMax = 0.0;
    for( j=1; j<=N; j++ )
      scaleMax = (scaleMax > fabs(A(i,j))) ? scaleMax : fabs(A(i,j));
    scale(i) = scaleMax;
  }

  //* Loop over rows k = 1, ..., (N-1)
  int signDet = 1;
  for( k=1; k<=(N-1); k++ ) {
	//* Select pivot row from max( |A(j,k)/s(j)| )
    double ratiomax = 0.0;
	int jPivot = k;
    for( i=k; i<=N; i++ ) {
      double ratio = fabs(A(index[i],k))/scale(index[i]);
      if( ratio > ratiomax ) {
        jPivot = i;
        ratiomax = ratio;
      }
    }
	//* Perform pivoting using row index list
	int indexJ = index[k];
	if( jPivot != k ) {	          // Pivot
      indexJ = index[jPivot];
      index[jPivot] = index[k];   // Swap index jPivot and k
      index[k] = indexJ;
	  signDet *= -1;			  // Flip sign of determinant
	}
	//* Perform forward elimination
    for( i=k+1; i<=N; i++ ) {
      double coeff = A(index[i],k)/A(indexJ,k);
      for( j=k+1; j<=N; j++ )
        A(index[i],j) -= coeff*A(indexJ,j);
      A(index[i],k) = coeff;
      b(index[i]) -= A(index[i],k)*b(indexJ);
    }
  }
  //* Compute determinant as product of diagonal elements
  double determ = signDet;	   // Sign of determinant
  for( i=1; i<=N; i++ )
	determ *= A(index[i],i);

  //* Perform backsubstitution
  x(N) = b(index[N])/A(index[N],N);
  for( i=N-1; i>=1; i-- ) {
    double sum = b(index[i]);
    for( j=i+1; j<=N; j++ )
      sum -= A(index[i],j)*x(j);
    x(i) = sum/A(index[i],i);
  }

  delete [] index;	 // Release allocated memory
  return( determ );
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//given 3 points, 2 points represent a line and the last point represents
//a point in anywhere in the 3-D...the function returns the point coordinat
// where the line intersects with the third point...it used ge function (Gauss Elimination)
Coordinate pointLineIntersection(Coordinate p1,Coordinate p2,Coordinate p3){
	Coordinate p;
	vector<double> equations, equationsResult;;

	double b1,b2,b3;
	Matrix coeff(3,3),bs(3,1),resultMtrx(3,1);

	b1 = ((p2.x - p3.x) * p1.x) + ((p2.y - p3.y) * p1.y) + ((p2.z - p3.z) * p1.z);
	b2 = (p2.x * (p2.y - p3.y)) - ((p2.x - p3.x) * p2.y);
	b3 = (p2.x * (p2.z - p3.z)) - ((p2.x - p3.x) * p2.z);


	coeff(1,1) = p2.x - p3.x;
	coeff(1,2) = p2.y - p3.y;
	coeff(1,3) = p2.z - p3.z;
	bs(1,1) = b1;

	coeff(2,1) = p2.y - p3.y;
	coeff(2,2) = -1*(p2.x - p3.x);
	coeff(2,3) = 0;
	bs(2,1) = b2;

	coeff(3,1) = p2.z - p3.z;
	coeff(3,2) = 0;
	coeff(3,3) = -1*(p2.x - p3.x);
	bs(3,1) = b3;

	double determ;

	determ = ge(coeff,bs,resultMtrx);  // Call Gauss Elimination

	p.x = resultMtrx(1,1);
	p.y = resultMtrx(2,1);
	p.z = resultMtrx(3,1);

	return p;

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//given a line and a point, returns the point where these two intersect
//it looks faster than the other one "Coordinate pointLineIntersection(Coordinate p1,Coordinate p2,Coordinate p3)"
float linePointIntersection(Coordinate lineP1, Coordinate lineP2, Coordinate p, Coordinate &intersectionPnt){

	//create the vectors represents the line and a nother vector represents a line from lineP1 and the point p
	Vectors v1(lineP2, lineP1), v2(p,lineP1);

	float dotProd = v1.dot (v2);
	float lineLengthSqr = (lineP2.x - lineP1.x)*(lineP2.x - lineP1.x) + (lineP2.y - lineP1.y)*(lineP2.y - lineP1.y) + (lineP2.z - lineP1.z)*(lineP2.z - lineP1.z);

	float e = dotProd/lineLengthSqr;


	intersectionPnt.x = lineP1.x + e * (lineP2.x - lineP1.x);
	intersectionPnt.y = lineP1.y + e * (lineP2.y - lineP1.y);
	intersectionPnt.z = lineP1.z + e * (lineP2.z - lineP1.z);

	/*
	if e is :
		negative	: the intersection point is before point#1 in the line
		0			: the intersection point is on the first point
		b/w (0,1)	: the intersection point is on the line
		> 1			: the intersection point is after the second point in the line
	*/
	return e;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//builds the rotation matrix for a given two points and an angle (in Radian)
void buildRotationMatrix (double mtx[4][4],Coordinate p1, Coordinate p, double angle)
{
// http://www.mines.edu/~gmurray/ArbitraryAxisRotation/ArbitraryAxisRotation.html


	Coordinate &origin = p1;
	Coordinate director;
	director.x = p.x-p1.x;
	director.y = p.y-p1.y;
	director.z = p.z-p1.z;

	double &a = origin.x;
	double &b = origin.y;
	double &c = origin.z;
	double &u = director.x;
	double &v = director.y;
	double &w = director.z;

	double u2 = u*u;
	double v2 = v*v;
	double w2 = w*w;

	double L2 = u2+v2+w2;
	double L  = sqrt(L2);

	double sin_a = sin(angle);
	double cos_a = cos(angle);


	mtx[0][0] = (u2+(v2+w2)*cos_a) / L2;
	mtx[1][0] = (u*v*(1-cos_a)-w*L*sin_a) / L2;
	mtx[2][0] = (u*w*(1-cos_a)+v*L*sin_a) / L2;
	mtx[3][0] = (a*(v2+w2)-u*(b*v+c*w)+(u*(b*v+c*w)-a*(v2+w2))*cos_a+(b*w-c*v)*L*sin_a) / L2;

	mtx[0][1] = (u*v*(1-cos_a)+w*L*sin_a) / L2;
	mtx[1][1] = (v2+(u2+w2)*cos_a) / L2;
	mtx[2][1] = (v*w*(1-cos_a)-u*L*sin_a) / L2;
	mtx[3][1] = (b*(u2+w2)-v*(a*u+c*w)+(v*(a*u+c*w)-b*(u2+w2))*cos_a+(c*u-a*w)*L*sin_a) / L2;

	mtx[0][2] = (u*w*(1-cos_a)-v*L*sin_a) / L2;
	mtx[1][2] = (v*w*(1-cos_a)+u*L*sin_a) / L2;
	mtx[2][2] = (w2+(u2+v2)*cos_a) / L2;
	mtx[3][2] = (c*(u2+v2)-w*(a*u+b*v)+(w*(a*u+b*v)-c*(u2+v2))*cos_a+(a*v-b*u)*L*sin_a) / L2;

	mtx[0][3] = 0;
	mtx[1][3] = 0;
	mtx[2][3] = 0;
	mtx[3][3] = 1;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// Given a rotation matrix with a point to be rotated
//////// Rotation Matrix should be built before the calling of this function using buildRotationMatrix function
void rotatePoint(double mtrx[4][4],Coordinate &v) {

	double vv[4], r[4];
	vv[0] = v.x;
	vv[1] = v.y;
	vv[2] = v.z;
	vv[3] = 1;
	for(int y=0; y<4; y++)
	{
		r[y] = 0;
		for(int x=0; x<4; x++)
		{
			r[y] += mtrx[x][y] * vv[x];
		}
	}
	v.x = r[0];
	v.y = r[1];
	v.z = r[2];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//compute the distance b/w 2 given points
double getDistance(Coordinate p1, Coordinate p2)
{
	return sqrt(((p1.x - p2.x)*(p1.x - p2.x)) + ((p1.y - p2.y)*(p1.y - p2.y)) + ((p1.z - p2.z)*(p1.z - p2.z)));
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// calculate the center of triangle represented by three points
Coordinate triangleCenter(Coordinate p1, Coordinate p2, Coordinate p3)
{
	Coordinate p;

	p.x = (p1.x + p2.x + p3.x)/3;
	p.y = (p1.y + p2.y + p3.y)/3;
	p.z = (p1.z + p2.z + p3.z)/3;

	return p;

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// calculate the center of triangle represented by three points
Coordinate quadCenter(Coordinate p1, Coordinate p2, Coordinate p3, Coordinate p4)
{
	Coordinate p;

	p.x = (p1.x + p2.x + p3.x + p4.x)/4;
	p.y = (p1.y + p2.y + p3.y + p4.y)/4;
	p.z = (p1.z + p2.z + p3.z + p4.z)/4;

	return p;

}

Coordinate quadCenter1(Coordinate p1)
{
	Coordinate p;

	p.x = (p1.x);
	p.y = (p1.y);
	p.z = (p1.z);

	return p;

}
Coordinate quadCenter2(Coordinate p1, Coordinate p2)
{
	Coordinate p;

	p.x = (p1.x + p2.x)/2;
	p.y = (p1.y + p2.y)/2;
	p.z = (p1.z + p2.z)/2;

	return p;

}
//////////////////////////////////////////////////////////////////////////////////////////
// given a line represent by p1-p2... find a point located on this line but has dist from p1 in the direction of p2
// dist = 0 return p1....dist = dist from p1 to p2 will return p2...
Coordinate pointOnLine(Coordinate p1, Coordinate p2, float dist)
{
	Vectors Pa(p1), Pab(p2, p1),
			Pc;

	float distPercantage = dist/Pab.length();

	Pc = Pab.mul(distPercantage);
	Pc += Pa;

	return Pc.getCoordinates();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//find the angle b/w 3 points...p2 is the angle's point
double getAngleDegree(Coordinate p1, Coordinate p2, Coordinate p3)
{
	Vectors v1(p1,p2), v2(p3,p2);

	return v1.getAngleDegree(v2);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double getAngleRadian(Coordinate p1, Coordinate p2, Coordinate p3)
{
	Vectors v1(p1,p2), v2(p3,p2);

	return v1.getAngleRadian(v2);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double getTorsionAngle(Coordinate p1,Coordinate p2,Coordinate p3,Coordinate p4){
	double A1,A2,B1,B2,C1,C2;
	long double angle;
	int sign = 1;

	// Plane  Ax + By + Cz + D = 0
        //  A = y1 ( z2 - z3 ) + y2 ( z3 - z1 ) + y3 ( z1 - z2 )
        //  B = z1 ( x2 - x3 ) + z2 ( x3 - x1 ) + z3 ( x1 - x2 )
        //  C= x1 ( y2 - y3 ) + x2 ( y3 - y1 ) + x3 ( y1 - y2 )
		//  D = - x1 ( y2z3 - y3z2 ) - x2 ( y3z1 - y1z3 ) - x3 ( y1z2 - y2z1 )
//	cout<<p1.x<<" "<<p2.x<<" "<<p3.x<<" "<<p4.x<<endl;
//	cout<<p1.y<<" "<<p2.y<<" "<<p3.y<<" "<<p4.y<<endl;
//	cout<<p1.z<<" "<<p2.z<<" "<<p3.z<<" "<<p4.z<<endl;

	A1 = p1.y *( p2.z - p3.z ) + p2.y *( p3.z - p1.z ) + p3.y *( p1.z - p2.z );
	B1 = p1.z *( p2.x - p3.x ) + p2.z *( p3.x - p1.x ) + p3.z *( p1.x - p2.x );
	C1 = p1.x *( p2.y - p3.y ) + p2.x *( p3.y - p1.y ) + p3.x *( p1.y - p2.y );
	//cout<<"Plane 1 = "<<A1<<" + "<<B1<<" + "<<C1<<endl;

	A2 = p2.y *( p3.z - p4.z ) + p3.y *( p4.z - p2.z ) + p4.y *( p2.z - p3.z );
	B2 = p2.z *( p3.x - p4.x ) + p3.z *( p4.x - p2.x ) + p4.z *( p2.x - p3.x );
	C2 = p2.x *( p3.y - p4.y ) + p3.x *( p4.y - p2.y ) + p4.x *( p2.y - p3.y );

	// Cos (x) = (A1A2 + B1B2 + C1C2)/ (sqrt(A1A1 + B1B1 + C1C1) * sqrt(A2A2 + B2B2 + C2C2))

//	cout<<"angle is "<<(A1*A2 + B1*B2 + C1*C2)<<" / "<<(sqrt(A1*A1 + B1*B1 + C1*C1) * sqrt(A2*A2 + B2*B2 + C2*C2))<<endl;
	angle = (A1*A2 + B1*B2 + C1*C2)/ (sqrt(A1*A1 + B1*B1 + C1*C1) * sqrt(A2*A2 + B2*B2 + C2*C2));

//	cout<<" ----------- angle "<<angle<<"  acos = "<<acos(angle)<<endl;

	if ((angle>0.999999999999999) && (angle<1.000000000000001))
		return 0;

	//get the sign of the torsion angle
	Vectors v1(p1,p2), v2(p3,p2), v4(p4), vNormal;
	vNormal	= v1.cross(v2);
	v4-=p3;

	if (v4.dot(vNormal)>0)
		sign = -1;				//sign is -

	if ((angle>-1.0000000001) && (angle<-0.9999999999))
		return sign * 180;			//the angle is 180

	if ((angle<-1) || (angle>1))
	{
		cout<<"============== in geometry::getTorsionAngle(Coord, Coord, Coord, Coord) ======="<<endl;
		cout<<"Error in computing Torsion angle..... the angle in radian exceed 1 or -1..it is "<<angle<<endl;
		cout<<"==============================================================================="<<endl;
		exit(1);
	}

	return sign * toDegree(acos(angle));		//return the signed angle in degree
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double getDistLines( Coordinate l1P1 /*line 1...point 1 */, Coordinate l1P2 /* line 1... point 2*/,
					 Coordinate l2P1 /*line 2...point 1 */, Coordinate l2P2 /* line 2... point 2*/){
    Vectors   u(l1P2,l1P1);
    Vectors   v(l2P2,l2P1);
    Vectors   w(l1P1,l2P1);
	double SMALL_NUM = 0.00001;
    double    a = u.dot(u);		//dot(u,u);        // always >= 0
    double    b = u.dot (v);	//dot(u,v);
    double    c = v.dot (v);	//dot(v,v);        // always >= 0
    double    d = u.dot (w);	//dot(u,w);
    double    e = v.dot(w);		//dot(v,w);
    double    D = a*c - b*b;       // always >= 0
    double    sc, sN, sD = D;      // sc = sN / sD, default sD = D >= 0
    double    tc, tN, tD = D;      // tc = tN / tD, default tD = D >= 0

    // compute the line parameters of the two closest points
    if (D < SMALL_NUM) { // the lines are almost parallel
        sN = 0.0;        // force using point P0 on segment S1
        sD = 1.0;        // to prevent possible division by 0.0 later
        tN = e;
        tD = c;
    }
    else {                // get the closest points on the infinite lines
        sN = (b*e - c*d);
        tN = (a*e - b*d);
        if (sN < 0.0) {       // sc < 0 => the s=0 edge is visible
            sN = 0.0;
            tN = e;
            tD = c;
        }
        else if (sN > sD) {  // sc > 1 => the s=1 edge is visible
            sN = sD;
            tN = e + b;
            tD = c;
        }
    }

    if (tN < 0.0) {           // tc < 0 => the t=0 edge is visible
        tN = 0.0;
        // recompute sc for this edge
        if (-d < 0.0)
            sN = 0.0;
        else if (-d > a)
            sN = sD;
        else {
            sN = -d;
            sD = a;
        }
    }
    else if (tN > tD) {      // tc > 1 => the t=1 edge is visible
        tN = tD;
        // recompute sc for this edge
        if ((-d + b) < 0.0)
            sN = 0;
        else if ((-d + b) > a)
            sN = sD;
        else {
            sN = (-d + b);
            sD = a;
        }
    }
    // finally do the division to get sc and tc
    sc = (ABS(sN) < SMALL_NUM ? 0.0 : sN / sD);
    tc = (ABS(tN) < SMALL_NUM ? 0.0 : tN / tD);


    // get the difference of the two closest points
    //Vector   dP = w + (sc * u) - (tc * v);  // = S1(sc) - S2(tc)
		Vectors trm1;
		trm1 = u;
		trm1 = trm1.mul(sc);
		Vectors trm2;
		trm2 = v;
		trm2 = trm2.mul(tc);
		Vectors dP;
		dP = w;
		dP += trm1;
		dP -= trm2;
		Vectors final;
		final = dP;
		return final.length();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//find the shortest distance between a line (infinite line ) and a point
//given any two points on the line
//in the case that you are not interested to find the intersection point then use this function b/s it looks faster than finding the intersection point
//then calculated the distance b/w intersection point and the line, which good in case you want the intersection point
double getDistInfLinePoint(Coordinate lineP1, Coordinate lineP2, Coordinate p){

	//create the vectors represents the line and a nother vector represents a line from lineP1 and the point p
	Vectors v1(lineP2, lineP1), v2(p,lineP1), crossV;
	//find the cross product which finds the area of the parallelogram formed by the two vectors
	crossV = v1.cross (v2);

	return crossV.length()/v1.length ();

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//find the shortest distance between a line segment (given by two points) and a point
double getDistLineSegPoint(Coordinate lineP1, Coordinate lineP2, Coordinate p){

	//create the vectors represents the line and a nother vector represents a line from lineP1 and the point p
	Vectors v1(lineP2, lineP1), v2(p,lineP1);

	float dotProd = v1.dot (v2);
	float lineLengthSqr = (lineP2.x - lineP1.x)*(lineP2.x - lineP1.x) + (lineP2.y - lineP1.y)*(lineP2.y - lineP1.y) + (lineP2.z - lineP1.z)*(lineP2.z - lineP1.z);

	float e = dotProd/lineLengthSqr;

	Coordinate newP;

	if (e < 0){
		newP.x = lineP1.x;
		newP.y = lineP1.y;
		newP.z = lineP1.z;
	}
	else{
		if (e>1){
			newP.x = lineP2.x;
			newP.y = lineP2.y;
			newP.z = lineP2.z;
		}
		else{
			newP.x = lineP1.x + e * (lineP2.x - lineP1.x);
			newP.y = lineP1.y + e * (lineP2.y - lineP1.y);
			newP.z = lineP1.z + e * (lineP2.z - lineP1.z);

		}
	}

	return getDistance(newP, p);

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void torsion2xyz(vector<Torsion> torsions, vector<Coordinate> &XYZlist)
{
	//https://lists.sdsc.edu/pipermail/pdb-l/2005-February/002294.html
	//use keyword dihedral 2 pdb converter
	double theta1 = (180 - C_N_Ca_ANGLE) * PI/180,
		   theta2 = (180 - 110) * PI/180,
		   theta3 = (180 - Ca_C_N_ANGLE) * PI/180;

	double omega = 180*PI/180;

	Coordinate curCoord;

	int i,res,l;


	for (i=0;i<torsions.size ();i++){
		torsions[i].phi = (torsions[i].phi * PI)/180;
		torsions[i].psi = (torsions[i].psi * PI)/180;
	}

	//initiate coordinates
	XYZlist.clear ();
	for (i=0; i<torsions.size ();i++)
	{
		//every residue has 4 atoms (N, CA, C, O);
		XYZlist.push_back (curCoord);
		XYZlist.push_back (curCoord);
		XYZlist.push_back (curCoord);
		XYZlist.push_back (curCoord);
	}

	l=-1;
	for (res=0; res<torsions.size ();res++){

		//set N coordinate
		extrude(XYZlist, C_N_BOND_LENGTH, theta1, torsions[res].phi, l);		//distance b/w Ci-1 and N = 1.32
		//set Ca coordinate
		extrude (XYZlist, N_Ca_BOND_LENGTH, theta2, torsions[res].psi, l);		//distance b/w N and Ca = 1.47
		//set C coordinate
		extrude(XYZlist, Ca_C_BOND_LENGTH, theta3, omega, l);					//distance b/w Ca and C = 1.53

		l++;

		//set O coordinate
		XYZlist[l].x = 0.0;
		XYZlist[l].y = 1.01575;		//1.01575 = 1.24 SIN 125
		XYZlist[l].z = 0.71123;		//0.71123 =-1.24 COS 125
	}
}
////////////////////////////////////////////////////////////////////////////////////////
void extrude (vector<Coordinate> &XYZlist,double DZ, double theta, double TAU, int &l){
	double cos1,sin1,cos2,sin2,y,z;


	cos1	= cos(theta);
	sin1	= sin(theta);
	cos2	= cos(TAU);
	sin2	= sin(TAU);

	for (int i=0;i<=l;i++){
		z	= XYZlist[i].z + DZ;
		y	= XYZlist[i].y * cos1 + z * sin1;
		XYZlist[i].z = z * cos1 - XYZlist[i].y * sin1;
		XYZlist[i].y = y * cos2 + XYZlist[i].x * sin2;
		XYZlist[i].x = -1 * y * sin2 + XYZlist[i].x * cos2;
	}
	l++;

	XYZlist[l].x = 0.0;
	XYZlist[l].y = 0.0;
	XYZlist[l].z = 0.0;

}

////////////////////////////////////////////////////////////////////////////////////////
///// This function generate (calculate) the coordination of the last point (C) in a triangle ABC
///// BC is the distance between point B and the point we are generating (C)
///// Theta is the angle between B in (ABC) in degree
Coordinate getPoint (Coordinate A, Coordinate B, double BC, double ABCangle)
{
	Coordinate Cpoint= B;

	//increment z coordinate of Cpoint by BC length
	Cpoint.z += BC;


	//initiate vectors
	Vectors ab(A,B),
			cb(Cpoint, B), abNormal;

	//get angle b/w two vectors
	double angle = ab.getAngleDegree(cb);

	angle -= ABCangle;

	//find the normal b/w two vectors
	abNormal = ab.cross(cb);

	abNormal += B;		//move the start of the normal from the origina to B point

	double mtx[4][4];

	//build rotation matrix (to rotate around the normal )
	buildRotationMatrix(mtx,B,abNormal.getCoordinates(),-toRadian(angle));

	//rotate around the normal
	rotatePoint(mtx, Cpoint);

	return Cpoint;
}
//////////////////////////////////////////////////////////////////////////////////////////////////
void peakClustering (vector<Coordinate> pnts, vector<vector<Coordinate> > & clusters, Coordinate sPnt, float r){

	int i, j;
	bool *clustered;			//flag for each point...1 if the point is already clustered
	int pClustered = 0;			//number of points clustered so far
	Coordinate cPnt= sPnt;		//cluster centroid


	clustered = new bool [pnts.size ()];
	for (i=0; i<pnts.size ();i++)
		clustered[i] = false;

	while (pClustered<pnts.size ()){

		vector<Coordinate> cluster;				//working on new cluster

		//find the closest point to the last centroid
		//at the beginning sPnt is the last centroid
		float minDist = 9999.0;
		int minIndx = -1;
		for (i=0; i<pnts.size (); i++){
			if (!clustered[i]){
				float pDist = round(getDistance(cPnt, pnts[i]),3);			//find the distance from each point to the starting point
				if (pDist<minDist){
					minDist = pDist;
					minIndx = i;
				}
			}
		}

		if (minIndx != -1){

			//save the point into the new cluster
			cluster.push_back (pnts[minIndx]);

			//set the new centroid
			cPnt = pnts[minIndx];
			clustered[minIndx] = true;
			pClustered++;
			bool cont = true;						//a flag to know if we have added new points in the last iteration or not

			while (pClustered<pnts.size () && cont){
				cont = false;
				i=0;
				while (pClustered<pnts.size () && i<pnts.size ()){
					double pDist = round(getDistance(pnts[i], cPnt),3);
					if (!clustered[i] && pDist<r){
						clustered[i] = true;
						pClustered++;

						cluster.push_back (pnts[i]);		//save the point into the cluster

						//re-calculate the centroid
						cPnt.x = (cPnt.x * (cluster.size ()-1) + pnts[i].x)/cluster.size ();
						cPnt.y = (cPnt.y * (cluster.size ()-1) + pnts[i].y)/cluster.size ();
						cPnt.z = (cPnt.z * (cluster.size ()-1) + pnts[i].z)/cluster.size ();

						cont = true;
					}
					i++;
				}
			}

			//add the centroid of the cluster at the end
			cluster.push_back (cPnt);

			//add the cluster
			clusters.push_back(cluster);
		}
	}
}
//////////////////////////////////////////////////////////////////////////////////////////////////
#endif
