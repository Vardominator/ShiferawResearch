#ifndef MORKOV_H
#define MORKOV_H

#include "constants.h"

int morkov(int &i, double &cpx,double &v, double &dt)
{


         int jj;
                  double poinf=1.0/(1.0+pow(e,-(v-vth)/s6));

         double alpha=poinf/taupo;
         double beta=(1.0-poinf)/taupo;

        double s2t=s1t*(r1/r2)*(xk2t/xk1t);
        double poi=1.0/(1.0+pow(e,-(v-vx)/sx));

         double xk3=(1.0-poi)/tau3;
         double xk3t=xk3;


         double prv=1.0-1.0/(1.0+pow(e,-(v-vy)/sy));

//       !          !  recovery

         double recov=10.0+4954.0*pow(e,v/15.6);

         double poix=1.0/(1.0+pow(e,-(v-vyr)/syr));

        double ra=1.0*(double)rand()/(double)RAND_MAX;
        double ri=ra/dt;

        double fca=1.0/(1.0+pow(cat/cpx,3.0));

        double s1=0.02*fca;
        double xk1=0.03*fca;

        double s2=s1*(r1/r2)*(xk2/xk1);  //reversibility conditions

        double tau_ca=tca/(1.0+pow(cpx/cpt,4.0));

        double tauca=(recov-tau_ca)*prv+tau_ca;
        double tauba=(recov-450.0)*prv+450.0;



        double xk6=fca*poix/tauca;
        double xk5=(1.0-poix)/tauca;

        double xk6t=poix/tauba;
        double xk5t=(1.0-poix)/tauba;

        double xk4=xk3*(alpha/beta)*(xk1/xk2)*(xk5/xk6);
        double xk4t=xk3t*(alpha/beta)*(xk1t/xk2t)*(xk5t/xk6t);


        double ragg = 1.0*(double)rand()/(double)RAND_MAX;

        double rig = ragg/dt;

                if(i==1)
                {//1
                    if(rig < beta)											jj = 2;
                    else if (beta < rig && rig < beta+r1)					jj = 7;
                    else if(beta+r1 < rig && rig < beta+r1+xk1t)			jj = 4;
                    else if(beta+r1+xk1t < rig && rig < beta+r1+xk1t+xk1)	jj = 3;
                    else													jj = 1;

    /*				if((beta+r1+xk1t+xk1)*dt >=1.0)
                    {
                      dt=1.0/(beta+r1+xk1t+xk1);
                      printf("Warning: Time step changed to dt= %f\n", dt);
                    }
    */
                }//1

                if(i==2)
                {//2
                    if(rig < xk6)											jj = 5;
                    else if(xk6 < rig && rig < xk6+xk6t) 					jj = 6;
                    else if(xk6+xk6t < rig && rig < xk6+xk6t+alpha) 	    jj = 1;
                    else													jj = 2;

    /*				if((xk6+xk6t+alpha)*dt >=1.0)
                    {
                      dt=1.0/(xk6+xk6t+alpha);
                      printf("Warning: Time step changed to dt= %f\n", dt);
                    }
    */
                }//2

                if(i==3)
                {//3
                    if(rig < xk3)											jj = 5;
                    else if(xk3 < rig && rig < xk3+xk2) 					jj = 1;
                    else if(xk3+xk2 < rig && rig <xk3+xk2+s2)				jj = 7;
                    else													jj = 3;

        /*			if((xk3+xk2+s2)*dt >=1.0)
                    { dt=1.0/(xk3+xk2+s2);
                      printf("Warning: Time step changed to dt= %f\n", dt);
                    }
        */
                }//3

                                    //1		c1(j)=0.0d0  ! C1
                                    //2		c2(j)=1.0d0  ! C2
                                    //3		xi1ca(j)=0.0d0 ! I1_Ca
                                    //4		xi1ba(j)=0.0d0 ! I1_Ba
                                    //5		xi2ca(j)=0.0d0 ! I2_Ca
                                    //6		xi2ba(j)=0.0d0 ! I2_Ba
                                    //7     open


                if(i==4)
                {//4
                    if(rig < xk3t)											jj = 6;
                    else if(xk3t < rig && rig < xk3t+xk2t)					jj = 1;
                    else if(xk3t+xk2t < rig && rig < xk3t+xk2t+s2t)			jj = 7;
                    else													jj = 4;

        /*			if((xk3t+xk2t+s2t)*dt >=1.0)
                    { dt=1.0/(xk3t+xk2t+s2t);
                      printf("Warning: Time step changed to dt= %f\n", dt);
                    }
        */
                }//4

                if(i==5)
                {//5
                    if(rig < xk5)											jj = 2;
                    else if(xk5 < rig && rig < xk5+xk4) 					jj = 3;
                    else													jj = 5;

        /*			if((xk5+xk4)*dt >=1.0)
                    { dt=1.0/(xk5+xk4);
                      printf("Warning: Time step changed to dt= %f\n", dt);
                    }
        */
                }//5

                if(i==6)
                {//6
                    if(rig < xk5t)											jj = 2;
                    else if(xk5t < rig && rig < xk5t+xk4t)					jj = 4;
                    else													jj = 6;

        /*			if((xk5t+xk4t)*dt >=1.0)
                    { dt=1.0/(xk5t+xk4t);
                      printf("Warning: Time step changed to dt= %f\n", dt);
                    }
        */
                }//6

                if(i==7)
                {//7
                    if(rig < r2)											jj = 1;
                    else if(r2 < rig && rig < r2+s1)						jj = 3;
                    else if(r2+s1 < rig && rig < r2+s1+s1t)					jj = 4;
                    else													jj = 7;

            /*		if((r2+s1+s1t)*dt >=1.0)
                    { dt=1.0/(r2+s1+s1t);
                      printf("Warning: Time step changed to dt= %f\n", dt);
                    }
            */
                }//7


return(jj);

}


#endif // MORKOV_H

