#ifndef PARAMETERS
#define PARAMETERS


//****************** Variables for all the compartments of each unit*******************

//################## Concentrations of different subunits #################
double  ci[nx*ny*nz];					//local cytosol
double  cjsr[nx*ny*nz];					//Jsr
double  cp[nx*ny*nz];					//proximal space
double  cs[nx*ny*nz];					//submembrance
double  cnsr[nx*ny*nz];					//Nsr
double  dotci[nx*ny*nz];				//For ODEs
double  dotcnsr[nx*ny*nz];				//For ODEs
double  dotcjsr[nx*ny*nz];				//For ODEs
//#####################################################################
double  cati[nx*ny*nz];					//Calcium binding to Troponin in cytosolic space
double  cats[nx*ny*nz];					//Calcium binding to Troponin in submembrane space
//@@@@@@@@@@@@@@@@@ Currents between different subunits @@@@@@@@@@@@@@@
double  xire[nx*ny*nz];					//Ir release current
double  xitr[nx*ny*nz];					//Jsr refilling current Nsr->Jsr
double  xica[nx*ny*nz];					//ica single Lcc current
double  xicat[nx*ny*nz];				//Ica total Lcc current
double  xiup[nx*ny*nz];					//Iup uptake current
double  xisi[nx*ny*nz];					//Isi diffusion from submembrane to myoplasm
double  xips[nx*ny*nz];					//Ips diffusion from proximal to submembrane
double  xitci[nx*ny*nz];				//Itci Troponin buffering current in cytosolic space
double  xitcs[nx*ny*nz];				//Itcs Troponin buffering current in submembrane space
double  xileak[nx*ny*nz];				//Ileak Leak from NSR to cytosolic
double  xinaca[nx*ny*nz];				//Inaca Sodium-calcium exchange current in submembrane space
double  xicoupi[nx*ny*nz];				//Coupling current in Cytosolic
double  xicoups[nx*ny*nz];				//Coupling current in submembrane
double  xicoupn[nx*ny*nz];				//Coupling current in NSR
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
double betalum[nx*ny*nz];				//luminal buffer factor in JSR
double betaci[nx*ny*nz];				//instantaneous buffer factor in cytosolic

//$$$$$$$$$$$$$$$$ Variables for Markov state for single lcc and RyR channels $$$$$$$$$$$$$$$$$
int  lcc[4*nx*ny*nz];					//State of single LCC channel
int  nou[nx*ny*nz];						//Number of the RyR channels in open-unbound state
int  ncu[nx*ny*nz];						//Number of the RyR channels in close-unbound state
int  nob[nx*ny*nz];						//Number of the RyR channels in open-bound state
int  ncb[nx*ny*nz];						//Number of the RyR channels in close-bound state
double po[nx*ny*nz];					//Open probability of RyR channels

//initial values
double  cproxit;						//For calculating the average concentration of the proximal space
double  csubt;							//For calculating the average concentration of the submembrane space
double  cit;							//For calculating the average concentration of the cytosolic space
double  cjsrt;							//For calculating the average concentration of the JSR space
double  cnsrt;							//For calculating the average concentration of the NSR space
double  catsto;							//For calculating the average concentration of the Calcium binding to Troponin in submembrane space
double  catito;							//For calculating the average concentration of the Calcium binding to Troponin in cytosolic space
double  cmax=0.0;						//For constructing the bifurcation curves, cmax would restore the peak value of CTAs
double  ccc0;							//For calculating the nearest-neighbor diffusive currents, Ici, Ics, IcNSR	(nx, ny, nz)
double  ccc1;							//For calculating the nearest-neighbor diffusive currents, Ici, Ics, IcNSR	(nx-1, ny, nz)
double  ccc2;							//For calculating the nearest-neighbor diffusive currents, Ici, Ics, IcNSR	(nx+1, ny. nz)
double  ccc3;							//For calculating the nearest-neighbor diffusive currents, Ici, Ics, IcNSR	(nx, ny-1, nz)
double  ccc4;							//For calculating the nearest-neighbor diffusive currents, Ici, Ics, IcNSR	(nx, ny+1, nz)
double  ccc5;							//For calculating the nearest-neighbor diffusive currents, Ici, Ics, IcNSR	(nx, ny, nz-1)
double  ccc6;							//For calculating the nearest-neighbor diffusive currents, Ici, Ics, IcNSR	(nx, ny, nz+1)
double  xireto;							//For calculating the average SR release current Ire
double  xicatto;						//For calculating the average Lcc calcium current in the proximal space
double  xinacato;						//For calculating the average NCX current in the submembrane space
double  poto;							//For calculating the average open probability of RyR channels
double xiremax;							//For plotting the release curve (bell curve)
double xicatmax;						//For plotting the peak Ica curve (bell curve)
double xicav;							//For plotting the Ica curve//
double csrtotal;						//For plotting the SR depletion
double csr0;							//Record the initial SR load during the first spark
double srmin;							//Record the minimum Sr load for depetion curve
int ndc;
int ndcx;
int jx;									//Loop variable
int jy;									//Loop variable
int jz;									//Loop variable
int j;									//Loop variable
int jj;									//Temporary variable
double uu1;								//Temporary variable. Random number
double uu2;								//Temporary variable. Random number
int nl;									//Record the number of open LCC channels in each unit when doing the LCC Gating dynamic
int m=0;								//For doing the Chudin's action potential clamp model
double vmax;							//For doing the Chudin's action potential clamp model
double vcd;								//For the loop, when xxx=1, vcd=90, so one would get the release curve with vmax varied from -40mV to 80mV.
double vmin=-80.0;						//For doing the Chudin's action potential clamp model
double tperiod = 200.0;							//Pacing period. For doing the Chudin's action potential clamp model.
double nbeat;						    //Number of beats

double xt=2.0/3.0/(2.0/3.0+tperiod/1000.0); //For doing the Chudin's action potential clamp model

double randomi[nx*ny*nz];				//For doing the fluctuations of the volumn of the proximal space
double t = 0.0;
double dt;								//pacing time interval in the unit of ms
double v;								//voltage in the unit of mV
double timeElapsed = 0.0;							//pacing time
double tt1;								//For record the moment that Ire reaches the peak when doing the bell curve
double tt2;								//For record the moment that Ica reaches the peak when doing the bell curve
double ttest;
double duration;						//The whole pacing duration
double xi = 0.39; 							//Coupling strength coefficient

#endif // PARAMETERS


