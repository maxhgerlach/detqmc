#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <cstdlib>

using namespace std;

double LinearRegression(const vector<double>& x, const vector<double>& y) {
  const int n = x.size();
  if(x.size() != y.size()) {
    std::cerr << "Error: Invalid vector sizes in linear regression" << std::endl;
    exit(1);
  }

  double ss_xy = 0.0;
  double ss_xx = 0.0;
  double mean_x = 0.0;
  double mean_y = 0.0;

  for(int i=0; i<n; ++i) {
    mean_x += x[i];
    mean_y += y[i];
  }
  mean_x /= n;
  mean_y /= n;
  
  for(int i=0; i<n; ++i) {
    ss_xy += x[i]*y[i];
    ss_xx += x[i]*x[i];
  }
  ss_xy -= n*mean_x*mean_y;
  ss_xx -= n*mean_x*mean_x;

  const double b = ss_xy / ss_xx;       // slope
  const double a = mean_y - b*mean_x;   // offset

  return b;
}   // LinearRegression


int main(int argc, char* argv[]) {
  if(argc==1) {
    std::cerr << "Usage: " << argv[0] << "  filename  [N]" << std::endl;
	std::cerr << "  The filename should indicate the file with the fraction of labeled walkers " << std::endl
			  << "  stored in a table, where the left column is the temperature (low to high) " << std::endl
			  << "  and the right column is the fraction (1 to 0)." << std::endl
			  << "  Optionally, you may indicate the number of temperatures [N] "<< std::endl
			  << "  in the new temperature set." << std::endl;
	exit(1);
  }
  
  // input stream
  std::ifstream in(argv[1]);

  // original temperature set
  std::vector<double> T;
  std::vector<double> deltaT;  
  // optimized temperature set
  std::vector<double> T_prime;
  std::vector<double> deltaT_prime;
  // some more stuff
  std::vector<double> Fraction;
  std::vector<double> Derivative;
  std::vector<double> Diffusivity;
  
  // read temperature and fraction
  while(!in.eof()) {
    double temp, frac;
    in >> temp >> frac >> ws;	
    T.push_back(temp);
	Fraction.push_back(frac);
  }
  
  // determine deltaT
  for(int i=0; i<T.size()-1; ++i) {
    deltaT.push_back(T[i+1]-T[i]);
  }

  // determine derivative of fraction
  for(int i=0; i<T.size()-1; ++i) {
    std::vector<double> x,y;
	x.clear(); y.clear();
	
	// regression with 3 points
    if(i) { 
      x.push_back(T[i-1]); y.push_back(Fraction[i-1]);
    }
	else {
	  // left boundary, add additional point
	  x.push_back(T[i+2]); y.push_back(Fraction[i+2]);
	}
	x.push_back(T[i  ]); y.push_back(Fraction[i  ]);
	x.push_back(T[i+1]); y.push_back(Fraction[i+1]);
	
	double slope = LinearRegression(x,y);
	Derivative.push_back(slope);
  }
  
  std::cerr << "Consistency check on derivatives ..." << std::flush;
  for(int i=0; i<Derivative.size(); ++i) {
    if(Derivative[i]>0) {
	  std::cerr << "\n  Derivative of fraction is positive at temperature " << T[i] << std::endl
	            << "  Your data may not be equilibrated." << std::endl
                    << "  To enforce feedback you may try to determine derivatives using a linear regression with more regression points." << std::endl;
	  exit(1);
	}
  }
  std::cerr << " done." << std::endl;
  
  std::cerr << "Writing calculated local diffusivity to file <diffusivity.dat> ..." << std::flush;
  ofstream out("diffusivity.dat");
  for(int i=0; i<Derivative.size(); ++i) {
    out << T[i] << "\t" << -deltaT[i]/Derivative[i] << std::endl;
  }
  std::cerr << " done." << std::endl;
  
  const int N = T.size()-1;  // number of intervals
  std::cerr << "Your current temperature set has N = " << N+1 << " temperature points." << std::endl;
  
  int NewN;
  if(argc>2) NewN = atoi(argv[2])-1;
  else NewN = N;
  std::cerr << "  The temperature set after feedback will have N = " << NewN+1 << " temperature points." << std::endl;
  

  // feedback of local diffusivity
  for(int t=0; t<N; ++t) {
	deltaT_prime.push_back( 1./sqrt( -1. / deltaT[t] * Derivative[t] ) );
  }
	
  // determine normalization
  double C = .0;
  if(deltaT.size() != N) { std::cerr << "ERROR: Inconsistent size of deltaT" << std::endl; exit(1); }
  for(int t=0; t<N; ++t) {
	C += deltaT[t] / deltaT_prime[t];
  }
  C = NewN/C;

  // determine new temperature set
  T_prime.push_back( T[0] );
  double n = 0.;
  double n_prime = 0.;
  int    t = 0;
  double this_deltaT = 0.;
  do {    
	if(n_prime + C * deltaT[t] / deltaT_prime[t] >= n+1) {
	  double tau = deltaT_prime[t]/C * (n - n_prime + 1);
	  this_deltaT += tau;
	  T_prime.push_back( T_prime.back() + this_deltaT );
	  n++;
	  n_prime += C * tau/deltaT_prime[t];
	  deltaT[t] -= tau;
	  this_deltaT = 0.;
	}
	else {
	  n_prime += C * deltaT[t] / deltaT_prime[t];
	  this_deltaT += deltaT[t];
	  t++;
	}
  } while(n<NewN && t<N);
  if(t==N) {
	T_prime.push_back( T[N] );
  }
  
  // output temperature set
  std::cerr << "\n\nNew temperature set:" << std::endl << std::endl;
  for(int i=0; i<T_prime.size(); ++i) {
    std::cout << i+1 << "\t" << T_prime[i] << std::endl;
  }
}
