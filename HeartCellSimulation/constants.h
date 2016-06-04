#ifndef CONSTANTS_H
#define CONSTANTS_H

const double vth = 0.0;
const double s6 = 4.0;
const double taupo = 1.0;			//ms
const double r1 = 0.3;				//0.3
const double r2 = 6;				//6			//3
const double cat = 0.5;				//0.5		//3
const double cpt = 1.5;				//1.5		//6.1
const double xk2 = 0.0005;
const double s1t = 0.00195;
const double xk1t = 0.00413;
const double xk2t = 0.00224;
const double vx = -40.0;			//40
const double sx = 3.0;
const double tau3 = 3.0;			//ms
const double vy = -40.0;			//40
const double sy = 4.0;
const double tca = 114.0;
const double vyr = -40.0;			//40
const double syr = 11.32;

const double e = 2.7182818;
const double pi = 3.1415926;
const double Vp = 0.00126;			//0.00126	//Volume of the proximal space
const double Vjsr = 0.02;						//Volume of the Jsr space
const double Vi = 0.5;							//Volume of the Local cytosolic
const double Vs = 0.025;						//Volume of the Local submembrane space
const double Vnsr = 0.025;						//Volume of the Local Nsr space
const int nryr = 100;						//Number of Ryr channels
const double taup = 0.022;			//0.022		//Diffusion time from the proximal to the submembrane
const double tausi = 0.1;						//Diffusion time from the submembrane to the cytosolic
const double tautr = 5.0;						//Diffusion time from NSR to JSR
const double taunl = 24.0;						//Longitudinal NSR
const double taunt = 7.2;						//Transverse NSR
const double tauil = 2.32;						//Longitudinal cytosolic
const double tauit = 2.93;						//Transverse cytosolic
const double tausl = 3.4;						//Longitudinal submembrane
const double taust = 1.42;						//Transverse submembrane

const int nx = 20;		//25		//65	//Number of Units in the x direction
const int ny = 10;		//10		//27	//Number of Units in the y direction
const int nz = 10;		//4			//11	//Number of Units in the z direction


#endif // CONSTANTS_H


