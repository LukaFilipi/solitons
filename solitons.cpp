#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <map>
using namespace std;

// Change equation parameters here
double L = 10, tmax = 5, h = 0.2, dt=.001;					// length of space and time intervals, space and time step size
int N;														// number of spatial divisions
double alpha = 1;											// soliton parameter
double alpha1 = 1.1, alpha2 = 1, dist1 = -18, dist2 = -12;	// soliton parameters and displacements of colliding waves
double K = 0.3;												// wavenumber of sinusoidal initial waveform


int _round(double x) { return (int)(x+0.5); }	// function to round numbers to closest integer (only works for positive numbers)


double getx(int i) { return (i-N/2)*h; }		// calculate position x from index i


// normalising index to interval [0,N]
int I(int i) { return (i%N+N)%N; }				// i%N shifts to [-N,N], +N shifts to [0,2N] and second %N to [0,N]


// function to print vector - used to print initial wave
void print_vector(vector<double> v, string name)
{
	ofstream outfile(name.c_str());
	for(int i=0; i<v.size(); i++)
		outfile << (i-((int)v.size()-1)/2) << '\t' << v[i] << endl;
	outfile.close();
}


// Function to print t, x and u
void print_map(map<double,vector<double> > m, string name)
{
	ofstream outfile(name.c_str());

	for(map<double,vector<double> >::iterator j=m.begin(); j!=m.end(); j++)
		for(int i=0; i<=N; i++)
			outfile << j->first << '\t' << getx(i) << '\t' << (j->second)[i] << endl;
	outfile.close();
}


// Function to print the waveform data at 6 regularly spaced time intervals
void print_map_in_rows(map<double,vector<double> > m, string name)
{
	ofstream outfile(name.c_str());

	outfile << "   t ->";
	for(map<double,vector<double> >::iterator j=m.begin(); j!=m.end(); j++)
	{
	  for(int row=0; row<5; row++)
	    if(_round(j->first/dt)==_round(tmax/dt)*row/5)
		  outfile << '\t' << j->first;
	}
	outfile << '\t' << (--m.end())->first << endl;

	for(int i=0; i<=N; i++)
	{
		outfile << getx(i);
		for(map<double,vector<double> >::iterator j=m.begin(); j!=m.end(); j++)
		{
		  for(int row=0; row<5; row++)
		    if(_round(j->first/dt)==_round(tmax/dt)*row/5)
				outfile << '\t' << (j->second)[i];
		}
		outfile << '\t' << ((--m.end())->second)[i] << endl;
	}
	outfile.close();
}


// Function to initialise normal mode wave
vector<double> initialise()
{
	vector<double> u;
	for(int i=0; i<=N; i++)			// from -L to L
		u.push_back( 12*alpha*alpha/(cosh(alpha*getx(i))*cosh(alpha*getx(i))) );  // calculate u at t=0 at all points along the spatial grid using normal mode solution
	return u;
}


// Function to initialise two colliding solitons
vector<double> initialise_collision()
{
	vector<double> u;
	for(int i=0; i<=N; i++)
		u.push_back( 12*alpha1*alpha1/(cosh(alpha1*(getx(i)-dist1))*cosh(alpha1*(getx(i)-dist1))) 
					+ 12*alpha2*alpha2/(cosh(alpha2*(getx(i)-dist2))*cosh(alpha2*(getx(i)-dist2))) );
		return u;
}


// Function to initialise sinusoidal wave of period 2pi/K
vector<double> initialise_sine()
{
	vector<double> u;
	for(int i=0; i<=N; i++)
		u.push_back( 12*alpha*alpha*cos(K*getx(i)) );		// keep same amplitude as normal mode
		return u;
}


// Function to calculate the discretised spatial derivatives for KdeV equation
double f(double plus1, double minus1, double plus2, double minus2, double h)
{
	return -(plus1*plus1-minus1*minus1)/(4*h) - (plus2-2*plus1+2*minus1-minus2)/(2*h*h*h);
}


// Function to calculate discretised spatial derivative for shock waves
double f_shock(double plus1, double minus1, double h)
{
	return -(plus1*plus1-minus1*minus1)/(4*h);
}


// Function to propagate wave in time using RK4 scheme
map<double,vector<double> > propagate(vector<double> u)
{
	map<double,vector<double> > wave;
	vector<double> k1(N+1), k2(N+1), k3(N+1), k4(N+1);  // RK4 coefficients

	tmax += dt/2;		// ensure tmax is included (manually account for errors with adding doubles)
	for(double t=0; t<=tmax; t+=dt)
	{
		if( _round(t/dt)%10==0) { wave[t]=u; }		// only store solution at every 10th time step to save memory

		// the function I is used to normalise the indices i
		for(int i=0; i<=N; i++)
			k1[i] = dt * f( u[I(i+1)], u[I(i-1)], u[I(i+2)], u[I(i-2)], h );
		for(int i=0; i<=N; i++)
			k2[i] = dt * f( u[I(i+1)]+k1[I(i+1)]/2, u[I(i-1)]+k1[I(i-1)]/2, u[I(i+2)]+k1[I(i+2)]/2, u[I(i-2)]+k1[I(i-2)]/2, h );
		for(int i=0; i<=N; i++)
			k3[i] = dt * f( u[I(i+1)]+k2[I(i+1)]/2, u[I(i-1)]+k2[I(i-1)]/2, u[I(i+2)]+k2[I(i+2)]/2, u[I(i-2)]+k2[I(i-2)]/2, h );
		for(int i=0; i<=N; i++)
			k4[i] = dt * f( u[I(i+1)]+k3[I(i+1)], u[I(i-1)]+k3[I(i-1)], u[I(i+2)]+k3[I(i+2)], u[I(i-2)]+k3[I(i-2)], h );

		for(int i=0; i<=N; i++)
			u[i] = u[i] + (k1[i]+2*k2[i]+2*k3[i]+k4[i])/6;
	}
	return wave;
}


// Function to propagate wave in time using RK4 scheme for shock waves
map<double,vector<double> > shock(vector<double> u)
{
	map<double,vector<double> > wave;
	vector<double> k1(N+1), k2(N+1), k3(N+1), k4(N+1);  // RK4 coefficients

	tmax += dt/2;		// ensure tmax is included (manually account for errors with adding doubles)
	for(double t=0; t<=tmax; t+=dt)
	{
		if( _round(t/dt)%10==0) { wave[t]=u; }		// only store solution at every 10th time step to save memory

		// the function I is used to normalise the indices i
		for(int i=0; i<=N; i++)
			k1[i] = dt * f_shock( u[I(i+1)], u[I(i-1)], h );
		for(int i=0; i<=N; i++)
			k2[i] = dt * f_shock( u[I(i+1)]+k1[I(i+1)]/2, u[I(i-1)]+k1[I(i-1)]/2, h );
		for(int i=0; i<=N; i++)
			k3[i] = dt * f_shock( u[I(i+1)]+k2[I(i+1)]/2, u[I(i-1)]+k2[I(i-1)]/2, h );
		for(int i=0; i<=N; i++)
			k4[i] = dt * f_shock( u[I(i+1)]+k3[I(i+1)], u[I(i-1)]+k3[I(i-1)], h );

		for(int i=0; i<=N; i++)
			u[i] = u[i] + (k1[i]+2*k2[i]+2*k3[i]+k4[i])/6;
	}
	return wave;
}


int main(void)
{
	N = _round(L/h) * 2;
	
	vector<double> u;								// initialise normal mode
	map<double,vector<double> > wave;				// u_i,j
	
	vector<double> col;								// initialise colliding solitons
	map<double,vector<double> > collision;			// u_i,j

	vector<double> sine;							// initialise sinusoidal waveform
	map<double,vector<double> > wavebreak;			// u_i,j

	map<double,vector<double> > shockwave;

	// Initialisation, propagation and printing of single soliton solution
	u=initialise();
	wave=propagate(u);
	print_vector(u, "initial_wave.txt");
	print_map_in_rows(wave, "propagation.txt");
	print_map(wave, "wave.txt");
	
	// Initialisation, propagation and printing of two-soliton collision
	col=initialise_collision();
	collision=propagate(col);
	print_vector(col, "initial_collision.txt");
	print_map_in_rows(collision, "collision_propagation.txt");
	print_map(collision, "collision.txt");

	// Initialisation, propagation and printing of sinusoidal solution
	sine=initialise_sine();
	wavebreak=propagate(sine);
	print_vector(sine, "initial_wavebreak.txt");
	print_map_in_rows(wavebreak, "wavebreak_propagation.txt");
	print_map(wavebreak, "wavebreak.txt");
	
	// Propagation and printing of shock waves
	shockwave=shock(u);
	print_map_in_rows(shockwave, "shockwave_propagation.txt");
	print_map(shockwave, "shockwave.txt");

}

