#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void init( cmplx* const psi0, const double eta, const double sigma, const double dx,
          const int Nx);

void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin, double t);

void lstep(cmplx* const f1, cmplx* const f0,
          const double dt, const double dx,
          const int Nx);

void nlstep(cmplx* const f1,
          const double dt,
          const int Nx);

//-----------------------------------
int main(){

	const int Nx = 4000;
	const double L = 800;
	const double xmin = 0;
	const double Tend = 50;
	const double dx = L / (Nx - 1);
	const double dt = dx  / 10;
	const int Na = 10;
	int Nk = int(Tend / Na / dt + 0.5);
	double t=0;
	const double eta = 0.2;

	stringstream strm;

	cmplx* psi0 = new cmplx[Nx];
	cmplx* psi1 = new cmplx[Nx];
	cmplx* h;
	
	init(psi0, eta, dx, dt,Nx);

	writeToFile(psi0,"psi_0", dx,Nx,xmin,t);


	for (int i = 1; i <= Na; i++) {
	    lstep(psi1,psi0,dt/2,dx,Nx);
		for (int j = 1; j <= Nk-1; j++) {
		      nlstep(psi1,dt,Nx);
		      h = psi0;
		      psi0 = psi1;
		      psi1 = h;
		      lstep(psi1,psi0,dt,dx,Nx);
		      t += dt;
		}
		    nlstep(psi1,dt,Nx);
		    h = psi0;
		    psi0 = psi1;
		    psi1 = h;
		    lstep(psi1,psi0,dt/2,dx,Nx);
	    
		    h = psi0;
		    psi0 = psi1;
		    psi1 = h;
		strm.str("");
		strm << "psi_" << i;
		writeToFile(psi0,strm.str(), dx,Nx,xmin,t);
	}

	delete[] psi0;
	delete[] psi1;
	
	return 0;
}
//-----------------------------------
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin, double t)
{
	ofstream out(s.c_str());
	for(int i=0; i<Nx; i++){
		double x = xmin + i * dx+5*t; //nebeneinander darstellen
		out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag() << endl;
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double eta,  const double dx, const double dt,
          const int Nx)
{
	const double x0 = dx*Nx * 0.5;
	const double f = sqrt(2) * eta;
	for(int i=0;i<Nx; i++){
		double x = i*dx - x0;
		psi0[i] = 2*f/cosh(eta * x);
	}
}
//-.---------------------------------

void lstep(cmplx* const f1, cmplx* const f0,
          const double dt, const double dx,
          const int Nx)
{
  const cmplx alpha = cmplx (0,-dt/(dx*dx));
  
  cmplx* d=new cmplx[Nx];

  for(int i=0;i<Nx;i++) d[i] = 1.0 + cmplx(2.0,0)*alpha;

  for(int i=1;i<Nx;i++){
    d[i]  -=   alpha*alpha/d[i-1];
    f0[i] -=    -alpha/d[i-1]*f0[i-1];
  }

  f1[Nx-1] = f0[Nx-1]/d[Nx-1];
  for(int i=Nx-2;i>0; i--)
    f1[i] = (f0[i] + alpha*f1[i+1])/d[i];

  delete[] d;
}

//-.---------------------------------

void nlstep(cmplx* const f1,
          const double dt,
          const int Nx)
{
  for(int i=0; i<Nx; i++){
    f1[i] *= exp(cmplx(0.,-(abs(f1[i])*abs(f1[i]))*dt));
  }
}