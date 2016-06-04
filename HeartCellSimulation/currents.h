
#ifndef CURRENTS_H
#define CURRENTS_H

#include "constants.h"

double ncx (double cs, double v,double tperiod)
{

 const double vnaca = 21.0;
 const double Kmcai = 0.00359;
 const double Kmcao = 1.3;
 const double Kmnai = 12.3;
 const double Kmnao = 87.5;
 const double Kda = 0.00011;				//mM
 const double ksat = 0.27;
 const double eta = 0.35;
 const double nao = 136.0;

 const double Farad = 96.5;		//	C/mmol
 const double xR = 8.314;		//	J/mol/K
 const double Temper = 308;		//K
 const double Cext = 1.8;		// mM
 const double a = 78.0;
 const double b = 10.0;

 double Inaca;
 double t1;
 double t2;
 double t3;
 double Ka;
 double nai;
 double za=v*Farad/xR/Temper;
 Ka=1.0/(1.0+pow(Kda/cs,3.0));
 nai=a/(1.0+b*sqrt(tperiod/1000.0));
 t1=Kmcai*pow(nao,3)*(1.0+pow(nai/Kmnai,3.0));
 t2=pow(Kmnao,3.0)*cs*(1.0+cs/Kmcai);						//right expression
// t2=pow(Kmnao,3)*cs+pow(Kmnai,3)*Cext*(1.0+cs/Kmcai);		//wrong expression
 t3=(Kmcao+Cext)*pow(nai,3)+cs*pow(nao,3);
 Inaca=Ka*vnaca*(pow(e,eta*za)*pow(nai,3)*Cext-pow(e,(eta-1.0)*za)*pow(nao,3)*cs)/((t1+t2+t3)*(1.0+ksat*pow(e,(eta-1.0)*za)));
 return (Inaca);

}

double instanbuf(double ci)
{double beta;
 double Kcam=7.0;
 double Bcam=24.0;
 double Ksr=0.6;
 double Bsr=47.0;
 double Kmca=0.033;
 double Bmca=140.0;
 double Kmmg=3.64;
 double Bmmg=140.0;
 double Kslh=0.3;
 double Bslh=13.4;
 double Bt=70.0;
 double kton=0.0327;
 double ktoff=0.0196;

 double temp1;
 double temp2;
 double temp3;
 double temp4;
 double temp5;

 temp1=Kcam*Bcam/pow((ci+Kcam),2.0);
 temp2=Ksr*Bsr/pow((ci+Ksr),2.0);
 temp3=Kmca*Bmca/pow((ci+Kmca),2.0);
 temp4=Kmmg*Bmmg/pow((ci+Kmmg),2.0);
// temp5=Kslh*Bslh/pow(ci+Kslh,2);
 beta=1.0/(1.0+temp1+temp2+temp3+temp4);
 return(beta);
}

double troponin(double CaT, double calciu)
{double Itc;
 double Bt=70.0;
 double kton=0.0327;
 double ktoff=0.0196;
 Itc=kton*calciu*(Bt-CaT)-ktoff*CaT;
 return(Itc);
 }

 double currentdsi( double cs, double ci)
{double Idsi;
 Idsi=(cs-ci)/tausi;
 return(Idsi);
}

double uptake(double ci, double cnsr)
{double Iup;
 double vup=0.3;				//0.3 for T=400ms
 double Ki=0.123;
 double Knsr=1700.0;			//1700 for T=400ms
 double HH=1.787;

 Iup=vup*(pow(ci/Ki,HH)-pow(cnsr/Knsr,HH))/(1.0+pow(ci/Ki,HH)+pow(cnsr/Knsr,HH));
 return(Iup);
}

double leak(double cnsr, double ci)
{double Ileak;
 double gleak=0.00001035;
 double Kjsr=500.0;

 Ileak=gleak*(cnsr-ci)*pow(cnsr,2.0)/(pow(cnsr,2.0)+pow(Kjsr,2.0));
 return(Ileak);
 }



double currenttr(double cnsr, double cjsr)
{

 double Itr;
 Itr=(cnsr-cjsr)/tautr;
 return (Itr);
}



double luminal(double cjsr)
{
const double nM = 15.0;
const double nD = 35.0;
const double ratedimer = 5000.0;
const double kdimer = 850.0;
const double hilldimer = 23.0;
const double CSQbers = 400.0;
const double kbers = 600.0;
double beta;
double roo2;
double mono;
double ene;
double enedot;


     roo2=ratedimer/(1.0+pow(kdimer/cjsr,hilldimer));
     mono=(-1.0+sqrt(1.0+8.0*CSQbers*roo2))/(4.0*roo2*CSQbers);
     ene=mono*nM+(1.0-mono)*nD;
     enedot=(nM-nD)*(-mono + 1.0/(4.0*400.0*mono*roo2 + 1.0))*(hilldimer/cjsr)*(1.0 - roo2/ratedimer);
//	printf("\n%f",roo2);
     beta=1.0/(1.0 + (CSQbers*kbers*ene + enedot*cjsr*(cjsr+kbers))/pow((kbers+cjsr),2.0));

     return(beta);

}

double currentdps(double cp, double cs)
{double Idps;
 Idps=(cp-cs)/taup;
 return (Idps);
}

double couplingI (double ccc0, double ccc1, double ccc2, double ccc3, double ccc4, double ccc5, double ccc6, double taul, double taut)
  { double Icoup;
    Icoup=(ccc2+ccc1-2.0*ccc0)/taul+(ccc4+ccc3-2.0*ccc0)/taut+(ccc6+ccc5-2.0*ccc0)/taut;
    return(Icoup);
  }


double release( double po, double cjsr, double cp)
{
 const double Jmax = 0.0147;					//0			//0.0147
 double Ir;
 Ir=Jmax*po*(cjsr-cp)/Vp;				// /randomi[i] outside
// printf("%f\n",Ir);
 return (Ir);
}

double poisson (double xxx)
{
 int yyy;
 double lll;
 int iii;
 double poi;
 double rrr;
 lll=pow(e,-xxx);
 iii=0;
 poi=1.0;
 while (poi >= lll)
 {
    iii=iii+1;
    rrr=1.0*(double)rand()/(double)RAND_MAX;
    poi=poi*rrr;
//	printf("%d %f\n",iii,lll);
 }
 yyy=iii-1;
// printf("%d %f %f\n ",iii,lll,poi);
 return(yyy);
}



double lcccurrent(double &v, double &cp)
{
 const double Pca = 11.9;		//	umol/C/ms
 const double gammai = 0.341;
 const double gammao = 0.341;
 const double Farad = 96.5;		//	C/mmol
 const double xR = 8.314;		//	J/mol/K
 const double Temper = 308;		//K
 const double Cext = 1.8;		// mM
 double za=v*Farad/xR/Temper;
 double ica;

 if (fabs(za)<0.001)
 {
     ica=2.0*Pca*Farad*gammai*(cp/1000.0*pow(e,2.0*za)-Cext);
 }
 else
 {
     ica=4.0*Pca*za*Farad*gammai*(cp/1000.0*pow(e,2.0*za)-Cext)/(pow(e,2.0*za)-1.0);
 }
     ica=(0.06/Vp)*(ica/Pca)/4.0;		// factor 0.5 is for testing whether cptilde=3 works with this model, because we don't want to change the whole cell Ica
 if (ica > 0.0) ica=0.0;
 return (ica);
}

#endif // CURRENTS_H

