//---------------------------------------------------------
//  Solver v 2.0
//---------------------------------------------------------
//  code produced by UID_0
//  begin at 24.05.2015
//  modified at
//  all rights reserved
//---------------------------------------------------------
//  function f is rhs of the equation y''=z'=f(x,y,z)
//  function g is rhs of y'=z=g(x,y,z)
//  y'' + 4y = cos(3x)
//---------------------------------------------------------
#include <math.h>
#include <fstream>
#include <iostream>

using namespace std;
std::ofstream outfile ("science.dat");

class DSolver
{
public:
					DSolver();
					~DSolver();
double*					X;
double*					Y;
double*					Z;
static const int 			n = 100;
void					solve(double Xo, double Yo, double Zo);
void					write();
double 					g (double x, double y, double z){return z;}
double 					f (double x, double y, double z){return cos(3*x)-4*y;}
};

DSolver::DSolver()
{
	X = new double[n];
	Y = new double[n];
	Z = new double[n];
}

DSolver::~DSolver()
{
	delete[] X;
	delete[] Y;
	delete[] Z;
}

inline void DSolver::solve(double Xo, double Yo, double Zo)
{
	double k1, k2, k3, k4;
	double q1, q2, q3, q4;
	double b = 1.0;
	X[0]=Xo;
	double h = (Xo-b)/n;

for (int i = 0; i < n; ++i)
{
   	X[i+1] = X[i] + h;
   	if (X[i]==b) break;
}

Y[0]=Yo;

for(int i=0; i<n; ++i)
{
	k1 = h * f(Xo, Yo, Zo);
   	q1 = h * g(Xo, Yo, Zo);
 
   	k2 = h * f(Xo + h/2.0, Yo + q1/2.0, Zo + k1/2.0);
   	q2 = h * g(Xo + h/2.0, Yo + q1/2.0, Zo + k1/2.0);
 
   	k3 = h * f(Xo + h/2.0, Yo + q2/2.0, Zo + k2/2.0);
   	q3 = h * g(Xo + h/2.0, Yo + q2/2.0, Zo + k2/2.0);
 
   	k4 = h * f(Xo + h, Yo + q3, Zo + k3);
   	q4 = h * g(Xo + h, Yo + q3, Zo + k3);
 
   	Z[i] = Zo + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
   	Y[i] = Yo + (q1 + 2.0*q2 + 2.0*q3 + q4)/6.0;

   	Zo = Z[i];
   	Yo = Y[i];

}
}

void DSolver::write()
{
	for (int i = 0; i < n; ++i)
   	{
      		outfile << scientific << X[i] <<" "<< Y[i] <<" "<<Z[i]<< endl;
   	}
   	outfile.close();
}

int main(int argc, char const *argv[])
{
	DSolver solver;
	solver.solve(0.0,0.8,2.0);
	solver.write();
	return 0;
}
