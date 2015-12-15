#include <cmath>
#include <iostream>
#include <cstdlib>
#include <fstream>

using namespace std;
// f[0] = x, f[1] = x';
// g[0] = y, g[1] = y';

//
void nofun(double* fs, double* gs,double* f, double* g, const double mu, double& r, double& s);
//
int main(){
	
	ofstream dafoog("Meh.txt");
	
	
	//AB für RK45
	const double mu = 0.012277471;
	const double T = 17.065216560157;
	double k1x[2], k2x[2], k3x[2], k4x[2], k5x[2], k6x[2], k7x[2];
	double k1y[2], k2y[2], k3y[2], k4y[2], k5y[2], k6y[2], k7y[2];
	double f[2], g[2], fRK5[2], gRK5[2], xtemp[2], ytemp[2];
	double r = 0;
	double s = 0;
	//
	
	//AB für step size control
	const double tol = 1e-5;
	double chi =10.0;
	double A;
	double dt = 0.0001;
	double dtalt = 0;
	double dtnew = 0;
	//
	
	// Anfangsbedingungen:
	f[0] = 0.994;
	f[1] = 0.0;
	g[0] = 0;
	g[1] = -2.00158510637908;
	
	fRK5[0] = 0.994;
	fRK5[1] = 0.0;
	gRK5[0] = 0;
	gRK5[1] = -2.00158510637908;
	
	r = sqrt((f[0] + mu)*(f[0] + mu) + g[0]*g[0]);
	s = sqrt((f[0] - 1 + mu)*(f[0] - 1 + mu) + g[0]*g[0]);
	//
	
	while (chi > tol) {
		
	nofun(k1x,k1y,f,g,mu,r,s);
	xtemp[0] = f[0] + 1/5.0 * dt * k1x[0];
	xtemp[1] = f[1] + 1/5.0 * dt * k1x[1];
	ytemp[0] = g[0] + 1/5.0 * dt * k1y[0];
	ytemp[1] = g[1] + 1/5.0 * dt * k1y[1];
	nofun(k2x,k2y,xtemp,ytemp,mu,r,s);
	
	xtemp[0] = f[0] + 3/40.0 * dt * k1x[0] + 9/40.0 * dt * k2x[0];
	xtemp[1] = f[1] + 3/40.0 * dt * k1x[1] + 9/40.0 * dt * k2x[1];
	ytemp[0] = g[0] + 3/40.0 * dt * k1y[0] + 9/40.0 * dt * k2y[0];
	ytemp[1] = g[1] + 3/40.0 * dt * k1y[1] + 9/40.0 * dt * k2y[1];	
	nofun(k3x,k3y,xtemp,ytemp,mu,r,s);
	
	xtemp[0] = f[0] + 44/45.0 * dt * k1x[0] - 56/15.0 * dt * k2x[0] + 32/9.0 * dt *k3x[0];
	xtemp[1] = f[1] + 44/45.0 * dt * k1x[1] - 56/15.0 * dt * k2x[1] + 32/9.0 * dt *k3x[1];
	ytemp[0] = g[0] + 44/45.0 * dt * k1y[0] - 56/15.0 * dt * k2y[0] + 32/9.0 * dt *k3y[0];
	ytemp[1] = g[1] + 44/45.0 * dt * k1y[1] - 56/15.0 * dt * k2y[1] + 32/9.0 * dt *k3y[1];
	nofun(k4x,k4y,xtemp,ytemp,mu,r,s);
	
	xtemp[0] = f[0] + 19372/6561.0 * dt * k1x[0] - 25360/2187.0 * dt * k2x[0] + 64448/6561.0 * dt * k3x[0] - 212/729.0 *dt *k4x[0];
	xtemp[1] = f[1] + 19372/6561.0 * dt * k1x[1] - 25360/2187.0 * dt * k2x[1] + 64448/6561.0 * dt * k3x[1] - 212/729.0 *dt *k4x[1];
	ytemp[0] = g[0] + 19372/6561.0 * dt * k1y[0] - 25360/2187.0 * dt * k2y[0] + 64448/6561.0 * dt * k3y[0] - 212/729.0 *dt *k4y[0];
	ytemp[1] = g[1] + 19372/6561.0 * dt * k1y[1] - 25360/2187.0 * dt * k2y[1] + 64448/6561.0 * dt * k3y[1] - 212/729.0 *dt *k4y[1];
	nofun(k5x,k5y,xtemp,ytemp,mu,r,s);
	
	xtemp[0] = f[0] + 9017/3168.0 * dt * k1x[0] - 355/33.0 * dt * k2x[0] + 46732/5247.0 * dt * k3x[0] + 49/176.0 * dt * k4x[0] - 5103/18656.0 * dt * k5x[0];
	xtemp[1] = f[1] + 9017/3168.0 * dt * k1x[1] - 355/33.0 * dt * k2x[1] + 46732/5247.0 * dt * k3x[1] + 49/176.0 * dt * k4x[1] - 5103/18656.0 * dt * k5x[1];
	ytemp[0] = g[0] + 9017/3168.0 * dt * k1y[0] - 355/33.0 * dt * k2y[0] + 46732/5247.0 * dt * k3y[0] + 49/176.0 * dt * k4y[0] - 5103/18656.0 * dt * k5y[0];
	ytemp[1] = g[1] + 9017/3168.0 * dt * k1y[1] - 355/33.0 * dt * k2y[1] + 46732/5247.0 * dt * k3y[1] + 49/176.0 * dt * k4y[1] - 5103/18656.0 * dt * k5y[1];
	nofun(k6x,k6y,xtemp,ytemp,mu,r,s);
	
	xtemp[0] = f[0] + 35/384.0 * dt * k1x[0] + 500/1113.0 * dt * k3x[0] + 125/192.0 * dt * k4x[0] - 2187/6784.0 * dt * k5x[0] + 11/84 * dt * k6x[0];
	xtemp[1] = f[1] + 35/384.0 * dt * k1x[1] + 500/1113.0 * dt * k3x[1] + 125/192.0 * dt * k4x[1] - 2187/6784.0 * dt * k5x[1] + 11/84 * dt * k6x[1];
	ytemp[0] = g[0] + 35/384.0 * dt * k1y[0] + 500/1113.0 * dt * k3y[0] + 125/192.0 * dt * k4y[0] - 2187/6784.0 * dt * k5y[0] + 11/84 * dt * k6y[0];
	ytemp[1] = g[1] + 35/384.0 * dt * k1y[1] + 500/1113.0 * dt * k3y[1] + 125/192.0 * dt * k4y[1] - 2187/6784.0 * dt * k5y[1] + 11/84 * dt * k6y[1];
	nofun(k7x,k7y,xtemp,ytemp,mu,r,s);
	
	// RK4
	f[0] += 35/384.0 * k1x[0] + 500/1113.0 * k3x[0] + 125/192.0 * k4x[0] - 2187/6784.0 * k5x[0] + 11/84 * k6x[0];
	f[1] += 35/384.0 * k1x[1] + 500/1113.0 * k3x[1] + 125/192.0 * k4x[1] - 2187/6784.0 * k5x[1] + 11/84 * k6x[1];
	g[0] += 35/384.0 * k1y[0] + 500/1113.0 * k3y[0] + 125/192.0 * k4y[0] - 2187/6784.0 * k5y[0] + 11/84 * k6y[0];
	g[1] += 35/384.0 * k1y[1] + 500/1113.0 * k3y[1] + 125/192.0 * k4y[1] - 2187/6784.0 * k5y[1] + 11/84 * k6y[1];
	
	//RK5
	fRK5[0] += 5179/57600.0 * k1x[0] + 7571/16695.0 * k3x[0] + 393/640.0 * k4x[0] - 92097/339200.0 * k5x[0] + 187/2100.0 * k6x[0] + 1/40.0 * k7x[0];
	fRK5[1] += 5179/57600.0 * k1x[1] + 7571/16695.0 * k3x[1] + 393/640.0 * k4x[1] - 92097/339200.0 * k5x[1] + 187/2100.0 * k6x[1] + 1/40.0 * k7x[1];
	gRK5[0] += 5179/57600.0 * k1y[0] + 7571/16695.0 * k3y[0] + 393/640.0 * k4y[0] - 92097/339200.0 * k5y[0] + 187/2100.0 * k6y[0] + 1/40.0 * k7y[0];
	gRK5[1] += 5179/57600.0 * k1y[1] + 7571/16695.0 * k3y[1] + 393/640.0 * k4y[1] - 92097/339200.0 * k5y[1] + 187/2100.0 * k6y[1] + 1/40.0 * k7y[1];
	
	
	chi = max(max(fabs(f[0] - fRK5[0]), fabs(f[1] - fRK5[1])),max(fabs(g[0] - gRK5[0]), fabs(g[1] - gRK5[1])));
	A = tol/chi;
	dtnew = dt  * 0.1 * pow(A,1/5.0);
	dtalt = dt;
	dt += dtnew;
	
	
	dtalt = dtalt + dt;
 	
	
	dafoog << dt << "\t" << dtalt << "\t" << f[0] << "\t" << g[0] << endl;
	
	
	}
		dafoog.close();
		return 0;
}
void nofun(double* fs, double* gs, double* f, double* g, const double mu, double& r, double& s){
	
	fs[0] = f[1];
	fs[1] = f[0] + 2 * g[1] - ((1 - mu) * (f[0] + mu))/(r*r*r) - (mu * (f[0] - 1 + mu))/(s*s*s);
	
	gs[0] = g[1];
	gs[1] = g[1] - 2 * f[1] - ((1 - mu) * g[1])/(r*r*r) - (mu * g[1])/(s*s*s);	
	
	r = sqrt((f[0] + mu)*(f[0] + mu) + g[0]*g[0]);
	s = sqrt((f[0] - 1 + mu)*(f[0] - 1 + mu) + g[0]*g[0]);
}






