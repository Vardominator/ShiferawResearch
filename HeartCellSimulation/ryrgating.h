#ifndef RYRGATING_H
#define RYRGATING_H

#include "constants.h"

void ryrgating ( double cp, double cjsr, int * ncu, int * nou, int * ncb, int * nob, double *dt)
{

    #define Ku 0.0003				//0.00038
    #define Kb 0.00005
    #define tauu 125.0				//125.0			//270.0
    #define taub 5.0				//5.0
    #define tauc 1.0


    #define ratedimer 5000.0
    #define kdimer 850.0
    #define hilldimer 23.0
    #define CSQbers 400.0
    #define kbers 600.0

    double roo2;
    double aub;
    double abu;

    double pku=Ku*pow(cp,2.0);
    double pkb=Kb*pow(cp,2.0);
    double	pkuminus=tauc;
    double pkbminus=tauc;
    double puu;
    double pkuh;
    double lamplus;
    double pkum;
    double lamminus;
    double pau;
    double lamau;
    double pbu;
    double lambu;
    double pcb;
    double lamcb;
    double pcub;
    double lamcub;
    double pkbh;
    double lamplbs;
    double pkbm;
    double lamminbs;

    double u1;
    double u2;
    double re;

    //int poisson1;
    //int poisson2;
    //int poisson3;
    //int poisson4;

    int n_ou_cu;
    int n_cu_ou;
    int n_ou_ob;
    int n_ob_ou;
    int n_cu_cb;
    int n_cb_cu;
    int n_ob_cb;
    int n_cb_ob;

    int kk;



    if (pku*(*dt) > 1.0) pku = 1.0/(*dt);
    if (pkb*(*dt) > 1.0) pkb = 1.0/(*dt);
    if (pkuminus*(*dt) > 1.0) pkuminus = 1.0/(*dt);
    if (pkbminus*(*dt) > 1.0) pkbminus = 1.0/(*dt);


    roo2=ratedimer/(1.0+pow(kdimer/(cjsr),hilldimer));
    aub=(-1.0+sqrt(1.0+8.0*CSQbers*roo2))/(4.0*roo2*CSQbers)/taub;
    abu=1.0/tauu;

   //__________________________ leaving OU >> CU
    pkuh = pkuminus*(*dt);
    lamplus = pow(e,-(*nou)*pkuh);
    n_ou_cu=-1;
    if ((*nou) <=1 || pkuh < 0.2 || ((*nou) <=5 && pkuh < 0.3))
       {//***********
                   kk = 0;
                   puu = 1.0;
                   while(puu >= lamplus)						//generates poisson number = fraction of closed RyR's that open
                   {
                       kk++;
                       re=1.0*(double)rand()/(double)RAND_MAX;
                       puu=puu*re;
                   }
                   n_ou_cu = kk-1;
               }
               else
               while(n_ou_cu < 0)
               {
               double u1=1.0*(double)rand()/(double)RAND_MAX;
               double u2=1.0*(double)rand()/(double)RAND_MAX;               //next is really a gaussian
               n_ou_cu = floor((*nou)*pkuh +sqrt((*nou)*pkuh*(1.0-pkuh))*sqrt(-2.0*log(1.0-u1))*cos(2.0*pi*u2))+rand()%2;
               }//**********
               if(n_ou_cu > nryr)
               {n_ou_cu = nryr;
   //			 printf("Ay ");
               }
               //____________________
               //_____________________leaving CU --> OU
               pkum = pku*(*dt);
               lamminus = pow(e,-(*ncu)*pkum);
               n_cu_ou = -1;
               if((*ncu) <= 1 ||pkum < 0.2 || ((*ncu) <= 5 && pkum < 0.3))	//checks if we use gaussian or poisson approx
               {//***********
                   kk = 0;
                   puu = 1.0;
                   while(puu >= lamminus)						//generates poisson number = fraction of closed RyR's that open
                       {
                       kk++;
                       re=1.0*(double)rand()/(double)RAND_MAX;
                       puu=puu*re;
                   }
                   n_cu_ou = kk-1;
               }
               else
               while(n_cu_ou < 0)
               {
               u1=1.0*(double)rand()/(double)RAND_MAX;
               u2=1.0*(double)rand()/(double)RAND_MAX;               //next is really a gaussian


               n_cu_ou = floor((*ncu)*pkum +sqrt((*ncu)*pkum*(1.0-pkum))*sqrt(-2.0*log(1.0-u1))*cos(2.0*pi*u2))+rand()%2;
               }//**********
               if(n_cu_ou > nryr) {n_cu_ou = nryr;}
               //____________________


               //_____________________leaving OU --> OB
               pau = aub*(*dt);
               lamau = pow(e,-(*nou)*pau);
               n_ou_ob = -1;
               if((*nou) <= 1 || pau < 0.2 || (*nou <= 5 && pau < 0.3))	//checks if we use gaussian or poisson approx
               {//***********
                   kk = 0;
                   puu = 1.0;
                   while(puu >= lamau)						//generates poisson number = fraction of open RyR's that close
                   {
                       kk++;
                       re=1.0*(double)rand()/(double)RAND_MAX;
                       puu=puu*re;
                   }
                   n_ou_ob = kk-1;
               }
               else
               while(n_ou_ob < 0)
               {
               u1=1.0*rand()/RAND_MAX;
               u2=1.0*(double)rand()/(double)RAND_MAX;					//next is really a gaussian
               n_ou_ob = floor(*nou*pau +sqrt(*nou*pau*(1.0-pau))*sqrt(-2.0*log(1.0-u1))*cos(2.0*pi*u2))+rand()%2;
               }//***********
               if(n_ou_ob > nryr) n_ou_ob = nryr;

               //______________________
               //_____________________leaving OB ---> OU
               pbu = abu*(*dt)*(Ku/Kb);
               lambu = pow(e,-(*nob)*pbu);
               n_ob_ou = -1;

               if((*nob) <= 1 || pbu < 0.2 || (*nob <= 5 && pbu < 0.3))	//checks if we use gaussian or poisson approx
               {//***********
                   kk = 0;
                   puu = 1.0;
                   while(puu >= lambu)						//generates poisson number = fraction of open RyR's that close
                       {
                       kk++;
                       re=1.0*(double)rand()/(double)RAND_MAX;
                       puu=puu*re;
                   }
                   n_ob_ou = kk-1;

               }

               else
               while(n_ob_ou < 0)
               {
               u1=1.0*(double)rand()/(double)RAND_MAX;
               u2=1.0*(double)rand()/(double)RAND_MAX;					//next is really a gaussian
               n_ob_ou = floor(*nob*pbu +sqrt(*nob*pbu*(1.0-pbu))*sqrt(-2.0*log(1.0-u1))*cos(2.0*pi*u2))+rand()%2;
   //			printf("%f\n",pbu);
               }//***********
               if(n_ob_ou > nryr) n_ob_ou = nryr;
               //______________________


               //_____________________leaving CB---->CU
               pcb = abu*(*dt);
               lamcb = pow(e,-(*ncb)*pcb);
               n_cb_cu = -1;
               if((*ncb) <= 1 || pcb < 0.2 	|| (pcb < 0.3 && (*ncb) <= 5))	//checks if we use gaussian or poisson approx
               {//***********
                   kk = 0;
                   puu = 1.0;
                   while(puu >= lamcb)						//generates poisson number = fraction of open RyR's that close
                   {
                       kk++;
                       re=1.0*(double)rand()/(double)RAND_MAX;
                       puu=puu*re;
                   }
                   n_cb_cu = kk-1;
               }
               else
               while(n_cb_cu < 0)
               {
               u1=1.0*(double)rand()/(double)RAND_MAX;
               u2=1.0*(double)rand()/(double)RAND_MAX;					//next is really a gaussian
               n_cb_cu = floor((*ncb)*pcb +sqrt((*ncb)*pcb*(1.0-pcb))*sqrt(-2.0*log(1.0-u1))*cos(2.0*pi*u2))+rand()%2;
               }//***********
               if(n_cb_cu > nryr) {n_cb_cu = nryr;}
               //______________________
               //_____________________leaving CU---->CB
               pcub = aub*(*dt);
               lamcub = pow(e,-*ncu*pcub);
               n_cu_cb = -1;
               if(*ncu <= 1 || pcub < 0.2 || (*ncu <= 5 && pcub < 0.3))	//checks if we use gaussian or poisson approx
               {//***********
                   kk = 0;
                   puu = 1.0;
                   while(puu >= lamcub)						//generates poisson number = fraction of open RyR's that close
                   {
                       kk++;
                       re=1.0*(double)rand()/(double)RAND_MAX;
                       puu=puu*re;
                   }
                   n_cu_cb = kk-1;
               }
               else
               while(n_cu_cb < 0)
               {
               u1=1.0*(double)rand()/(double)RAND_MAX;
               u2=1.0*(double)rand()/(double)RAND_MAX;					//next is really a gaussian
               n_cu_cb = floor(*ncu*pcub +sqrt(*ncu*pcub*(1.0-pcub))*sqrt(-2.0*log(1.0-u1))*cos(2.0*pi*u2))+rand()%2;
               }//***********
               if(n_cu_cb > nryr) {n_cu_cb = nryr;}
               //______________________


               //_____________________leaving OB --> CB
               pkbh = pkbminus*(*dt);
               lamplbs = pow(e,-(*nob)*pkbh);
               n_ob_cb = -1;
               if((*nob) <= 1 || pkbh < 0.2 || ((*nob) <= 5 && pkbh < 0.3))	//checks if we use gaussian or poisson approx
               {//***********
                   kk = 0;
                   puu = 1.0;
                   while(puu >= lamplbs)						//generates poisson number = fraction of closed RyR's that open
                   {
                       kk++;
                       re=1.0*(double)rand()/(double)RAND_MAX;
                       puu=puu*re;
                   }
                   n_ob_cb = kk-1;
               }
               else
               while(n_ob_cb < 0)
               {
               u1=1.0*(double)rand()/(double)RAND_MAX;
               u2=1.0*(double)rand()/(double)RAND_MAX;               //next is really a gaussian
               n_ob_cb = floor((*nob)*pkbh +sqrt((*nob)*pkbh*(1.0-pkbh))*sqrt(-2.0*log(1.0-u1))*cos(2.0*pi*u2))+rand()%2;
               }//**********
               //if(n_ob_cb > nn) {n_ob_cb = nn;printf("Ay ");}
               //____________________
               //_____________________leaving CB --> OB
               pkbm = pkb*(*dt);
               lamminbs = pow(e,-(*ncb)*pkbm);
               n_cb_ob = -1;
               if((*ncb) <= 1 || pkbm < 0.2 || ((*ncb) <= 5 && pkbm < 0.3))	//checks if we use gaussian or poisson approx
               {//***********
                   kk = 0;
                   puu = 1.0;
                   while(puu >= lamminbs)						//generates poisson number = fraction of closed RyR's that open
                   {
                       kk++;
                       re=1.0*(double)rand()/(double)RAND_MAX;
                       puu=puu*re;
                   }
                   n_cb_ob = kk-1;
               }
               else
               while(n_cb_ob < 0)
               {
               u1=1.0*(double)rand()/(double)RAND_MAX;
               u2=1.0*(double)rand()/(double)RAND_MAX;               //next is really a gaussian
               n_cb_ob = floor((*ncb)*pkbm +sqrt((*ncb)*pkbm*(1.0-pkbm))*sqrt(-2.0*log(1.0-u1))*cos(2.0*pi*u2))+rand()%2;
               }//**********
               if(n_cb_ob > nryr) {n_cb_ob = nryr;}
               //____________________

               //printf("\naub = %f    abu=%f",aub,abu);

               if(n_ou_ob  +  n_ou_cu > *nou)
               {
               if(n_ou_cu >= n_ou_ob) n_ou_cu = 0;
               else  n_ou_ob = 0;
               if (n_ou_ob > *nou) n_ou_ob = 0;
               else if(n_ou_cu > *nou) n_ou_cu = 0;
               }

               if(n_ob_ou  +  n_ob_cb > *nob)
               {
               if(n_ob_ou >= n_ob_cb) n_ob_ou = 0;
               else  n_ob_cb = 0;
               if (n_ob_cb > *nob) n_ob_cb = 0;
               else if(n_ob_ou > *nob) n_ob_ou = 0;
               }

               if(n_cu_ou  +  n_cu_cb > *ncu )
               {
               if(n_cu_cb >= n_cu_ou) n_cu_cb = 0;
               else  n_cu_ou = 0;
               if (n_cu_ou > *ncu) n_cu_ou = 0;
               else if(n_cu_cb > *ncu) n_cu_cb = 0;
               }


               *nou += -n_ou_ob  -  n_ou_cu    +n_ob_ou  + n_cu_ou;
               if(*nou<0)				{/*printf("Ay ou = %d %d %d ||",*nou,n_ou_ob,n_ou_cu);*/(*nou)=0;}
               if(*nou>nryr)			{/*printf("Ay ou = %d ",*nou);*/*nou=nryr;}

               *nob += -n_ob_ou  -  n_ob_cb    +n_ou_ob  + n_cb_ob;
               if(*nob<0)				{/*printf("Ay ob = %d %d %d ||",*nob,n_ob_ou,n_ob_cb);*/*nob=0;}
               if(*nob>nryr)			{/*printf("Ay ob = %d ",*nob);*/*nob=nryr;}

               *ncu += -n_cu_ou  -  n_cu_cb    +n_ou_cu  + n_cb_cu;
               if(*ncu<0)			    {/*printf("Ay cu = %d %d %d ||",*ncu,n_cu_ou,n_cu_cb);*/*ncu=0;}
               if(*ncu>nryr)			{/*printf("Ay cu = %d ",*ncu);*/*ncu=nryr;}

               *ncb += -n_cb_cu  -  n_cb_ob    +n_ob_cb  + n_cu_cb;
               if(*ncb<0)			    {/*printf("Ay cb = %d %d %d ||",*ncb,n_cb_cu,n_cb_ob);*/*ncb=0;}
               if(*ncb>nryr)			{/*printf("Ay cb = %d ",*ncb);*/*ncb=nryr;}


   }



#endif // RYRGATING_H


