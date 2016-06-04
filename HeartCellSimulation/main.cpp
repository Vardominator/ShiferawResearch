#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>
#include <cstdlib>
#include <ctime>

#include "constants.h"
#include "morkov.h"
#include "currents.h"
#include "ryrgating.h"
#include "parameters.h"

#include <string>
#include <sstream>


using namespace std;


int main(int argc, char* argv[])

{

    int params[4];
    int nx;
    int ny;
    int nz;

    int thresh;

    ifstream inputParams;
    inputParams.open("params.txt");

	for(int i = 0; i < 4; i++)
	{

		inputParams >> params[i];

	}

	inputParams.close();
	nx = params[0];
	ny = params[1];
	nz = params[2];
	thresh = params[3];

    //nx = atoi(argv[1]);
    //ny = atoi(argv[2]);
    //nz = atoi(argv[3]);

    cout << endl;

    srand(time(0));


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


    double randomi[nx*ny*nz];				//For doing the fluctuations of the volumn of the proximal space



    // Filename variables
    int fileCounter = 1;
    string fileNumToStr = "";
    string newFileName = "";

    // Time rangle variables
    double time1 = 0.0;
    double time2 = time1 + 0.1;

    // Loading screen variables
    int loadingVariable;
    int checkNum = 1;

    // Volume changing variables from Heart Failure simulation
    double height;        //Mean of the normal distribution
    double sigma = (0.10000) * height;  //Standard deviation of the normal distribution
    double vol = 0.0; //Varying proximity volume
    bool isPositive = false;
    double goVol;
    default_random_engine generator;
    height = 5.0;
    sigma = height * 0.1;
    normal_distribution<double> normal(height, sigma);


    // Initiate loading screen
    cout << "Loading Simulation..." << endl << endl << endl;

    // Initial load of cell
    double srLoad = 1750.0;



    // This portion regards the activity mapping
    int effectiveNumUnits = (nx - 2) * (ny - 2) * (nz - 2);   //This is the effective number of units, taking into account the boundaries.
    bool** testArray = new bool*[1202];
    int trueCounter = 0;

    for (int num1 = 0; num1 < 1202; num1++)
    {
	testArray[num1] = new bool[effectiveNumUnits];
    }



    for (int num1 = 0; num1 < 1202; num1++)
    {
  	for (int num2 = 0; num2 < effectiveNumUnits; num2++)
	{

    		testArray[num1][num2] = false;

	}
    }





//$$$$$$$$$$$$$$ Randome volume of proximal space $$$$$$$$$$$$$$$$$$$$$$$$$$$$
    for (jz = 0; jz < nz; jz++)
    for (jy = 0; jy < ny; jy++)
    for (jx = 0; jx < nx; jx++)
    {

        while(!isPositive)
        {

            vol = 0.1 * normal(generator);

            if(vol > 0.0)
            {

                isPositive = true;
                randomi[jx+nx*jy+nx*ny*jz] = vol;

            }

        }

        isPositive = false;

    }


    vmax=0.0;
    vcd=1.0;
    ndc=0;
    ndcx=1;



        for (jz = 0; jz < nz; jz++)
        for (jy = 0; jy < ny; jy++)
        for (jx = 0; jx < nx; jx++)
//		for (i = 0; i < Nsample; i++)						//For each voltage, giving initial values of all the compartments
        {

            ci[jx+nx*jy+nx*ny*jz]=0.122;
            cp[jx+nx*jy+nx*ny*jz]=0.1;
            cs[jx+nx*jy+nx*ny*jz]=0.1;

            if((jx >= 28 && jx <= 31) && (jy >= 8 && jy <= 11) && (jz >= 8 && jz <= 11))
            {
                cjsr[jx+nx*jy+nx*ny*jz]=1600;  //750.0
                cnsr[jx+nx*jy+nx*ny*jz]=1600;  //750.0
                //cout << "Hello";
            }
            else
            {
                cjsr[jx+nx*jy+nx*ny*jz]=1500;  //750.0
                cnsr[jx+nx*jy+nx*ny*jz]=1500;  //750.0
                //cout << "hello";
            }

            cati[jx+nx*jy+nx*ny*jz]=20.0;
            cats[jx+nx*jy+nx*ny*jz]=20.0;

            for(j = 0; j < 4; j++)
            {
                lcc[j+4*(jx+nx*jy+nx*ny*jz)]=2;
            }
            //Initial State of LCC channels

            ncb[jx+nx*jy+nx*ny*jz]=70;					//Initial State of RyR channels
            ncu[jx+nx*jy+nx*ny*jz]=30;
            nob[jx+nx*jy+nx*ny*jz]=0;
            nou[jx+nx*jy+nx*ny*jz]=0;

        }



//********** START TIME PERIOD LOOP *****************		 	 TT

        t=0.0;
        m=0;
        dt=0.1;
        xiremax=0.0;
        xicatmax=0.0;
        xicav=0.0;

        tperiod = 120.0;  //400;
        nbeat = 1;      //40;


        duration=nbeat*tperiod;

        while (m < nbeat)
        {


            ostringstream fullCellOut;
            fullCellOut << "/home/varderes/Desktop/CellSimulation/FullCell/FullCellBig" << fileCounter << ".txt";

            ostringstream spreading;
            spreading << "/home/varderes/Desktop/CellSimulation/Spreading/FullCellSpreading" << fileCounter << ".txt";
            

            ostringstream spreadingActivity;
            spreadingActivity << "/home/varderes/Desktop/CellSimulation/ActivityMapping/SpreadingWithActivity" << fileCounter << ".txt";

            ofstream out(fullCellOut.str().c_str(), ios_base::binary);
            ofstream spreadOut(spreading.str().c_str(), ios_base::binary);
            ofstream spreadActivityOut(spreadingActivity.str().c_str(), ios_base::binary);


            v = -80;

//Boundary conditions//!!!!!!!!!!!!!!!!!!!!!!!!! No flux!!!!!!!!!
            for (jx = 1; jx < nx-1; jx++)
            for (jy = 1; jy < ny-1; jy++)
                {
                    ci[jx+nx*jy+0*nz] = ci[jx+nx*jy+nx*ny*1];					//When the unit is at the border we impose
                    cs[jx+nx*jy+0*nz] = cs[jx+nx*jy+nx*ny*1];					//no-flux boundary conditions
                    cnsr[jx+nx*jy+0*nz] = cnsr[jx+nx*jy+nx*ny*1];

                    ci[jx+nx*jy+nx*ny*(nz-1)] = ci[jx+nx*jy+nx*ny*(nz-2)];
                    cs[jx+nx*jy+nx*ny*(nz-1)] = cs[jx+nx*jy+nx*ny*(nz-2)];
                    cnsr[jx+nx*jy+nx*ny*(nz-1)] = cnsr[jx+nx*jy+nx*ny*(nz-2)];
                }

            for (jx = 1; jx < nx-1; jx++)
            for (jz = 1; jz < nz-1; jz++)
                {
                    ci[jx+nx*0+nx*ny*jz] = ci[jx+nx*1+nx*ny*jz];
                    cs[jx+nx*0+nx*ny*jz] = cs[jx+nx*1+nx*ny*jz];
                    cnsr[jx+nx*0+nx*ny*jz] = cnsr[jx+nx*1+nx*ny*jz];

                    ci[jx+nx*(ny-1)+nx*ny*jz] = ci[jx+nx*(ny-2)+nx*ny*jz];
                    cs[jx+nx*(ny-1)+nx*ny*jz] = cs[jx+nx*(ny-2)+nx*ny*jz];
                    cnsr[jx+nx*(ny-1)+nx*ny*jz] = cnsr[jx+nx*(ny-2)+nx*ny*jz];
                }

            for (jy = 1; jy < ny-1; jy++)
            for (jz = 1; jz < nz-1; jz++)
                {
                    ci[0+nx*jy+nx*ny*jz] = ci[1+nx*jy+nx*ny*jz];
                    cs[0+nx*jy+nx*ny*jz] = cs[1+nx*jy+nx*ny*jz];
                    cnsr[0+nx*jy+nx*ny*jz] = cnsr[1+nx*jy+nx*ny*jz];

                    ci[(nx-1)+nx*jy+nx*ny*jz] = ci[(nx-2)+nx*jy+nx*ny*jz];
                    cs[(nx-1)+nx*jy+nx*ny*jz] = cs[(nx-2)+nx*jy+nx*ny*jz];
                    cnsr[(nx-1)+nx*jy+nx*ny*jz] = cnsr[(nx-2)+nx*jy+nx*ny*jz];
                }

//End up Boundary conditions!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            cproxit=0.0;
            csubt=0.0;
            cit=0.0;
            cjsrt=0.0;
            cnsrt=0.0;
            catsto=0.0;
            catito=0.0;
            xireto=0.0;
            xicatto=0.0;
            xinacato=0.0;
            poto=0.0;
            csrtotal=0.0;

	int counter = 0;

            for (jz = 1; jz < nz-1; jz++)
            for (jy = 1; jy < ny-1; jy++)
            for (jx = 1; jx < nx-1; jx++)
            {

		
		counter++;
		//cout << counter << endl;

                nl = 0;

                for (j = 0; j < 4; j++)
                {

                    jj=morkov(lcc[j+4*(jx+nx*jy+nx*ny*jz)],cp[jx+nx*jy+nx*ny*jz],v,dt);			//"morkov" gives the state of LCC channels
                    lcc[j+4*(jx+nx*jy+nx*ny*jz)]=jj;

                    if (jj == 7)
                    {
                        nl=nl+1;
                    }

                }


                xicat[jx+nx*jy+nx*ny*jz]=nl*(0.6)*lcccurrent( v, cp[jx+nx*jy+nx*ny*jz])/randomi[jx+nx*jy+nx*ny*jz];					//The unit of cp would be converted to [mM] in routine "lcccurrent"
                //cout << xicat[jx+nx*jy+nx*ny*jz] << endl;
                ryrgating (cp[jx+nx*jy+nx*ny*jz], cjsr[jx+nx*jy+nx*ny*jz],  &ncu[jx+nx*jy+nx*ny*jz], &nou[jx+nx*jy+nx*ny*jz], &ncb[jx+nx*jy+nx*ny*jz], &nob[jx+nx*jy+nx*ny*jz], &dt);	//RyR gating dynamic
                po[jx+nx*jy+nx*ny*jz]=(double)(nou[jx+nx*jy+nx*ny*jz]+nob[jx+nx*jy+nx*ny*jz])/nryr;		//the fraction of open RyR
                xire[jx+nx*jy+nx*ny*jz]=release( po[jx+nx*jy+nx*ny*jz], cjsr[jx+nx*jy+nx*ny*jz], cp[jx+nx*jy+nx*ny*jz])/randomi[jx+nx*jy+nx*ny*jz];		//release current

                if (xire[jx+nx*jy+nx*ny*jz] < 0.0) xire[jx+nx*jy+nx*ny*jz]=0.0;

                xinaca[jx+nx*jy+nx*ny*jz]=ncx(cs[jx+nx*jy+nx*ny*jz]/1000.0,v,tperiod);	//cs/1000 [mM]	NCX
                xitci[jx+nx*jy+nx*ny*jz]=troponin( cati[jx+nx*jy+nx*ny*jz], ci[jx+nx*jy+nx*ny*jz]);		//I_TCi
                betaci[jx+nx*jy+nx*ny*jz]=instanbuf(ci[jx+nx*jy+nx*ny*jz]);								//instantaneous buffer
                xitcs[jx+nx*jy+nx*ny*jz]=troponin( cats[jx+nx*jy+nx*ny*jz], cs[jx+nx*jy+nx*ny*jz]);		//I_TCs
                xiup[jx+nx*jy+nx*ny*jz]=uptake(ci[jx+nx*jy+nx*ny*jz], cnsr[jx+nx*jy+nx*ny*jz]);			//I_up
                xileak[jx+nx*jy+nx*ny*jz]=leak(cnsr[jx+nx*jy+nx*ny*jz], ci[jx+nx*jy+nx*ny*jz]);			//I_leak
                xips[jx+nx*jy+nx*ny*jz]=currentdps(cp[jx+nx*jy+nx*ny*jz], cs[jx+nx*jy+nx*ny*jz]);		//I_dps
                xisi[jx+nx*jy+nx*ny*jz]=currentdsi(cs[jx+nx*jy+nx*ny*jz], ci[jx+nx*jy+nx*ny*jz]);		//I_dsi
                xitr[jx+nx*jy+nx*ny*jz]=currenttr(cnsr[jx+nx*jy+nx*ny*jz], cjsr[jx+nx*jy+nx*ny*jz]);	//I_tr
                betalum[jx+nx*jy+nx*ny*jz]=luminal(cjsr[jx+nx*jy+nx*ny*jz]);							//Luminal buffer
                dotcjsr[jx+nx*jy+nx*ny*jz]=betalum[jx+nx*jy+nx*ny*jz]*(xitr[jx+nx*jy+nx*ny*jz]-xire[jx+nx*jy+nx*ny*jz]*Vp*randomi[jx+nx*jy+nx*ny*jz]/Vjsr);





                if((jx >= 28 && jx <= 31) && (jy >= 8 && jy <= 11) && (jz >= 8 && jz <= 11))
                {
                    cp[jx+nx*jy+nx*ny*jz]=cs[jx+nx*jy+nx*ny*jz]+taup*(xire[jx+nx*jy+nx*ny*jz]-xicat[jx+nx*jy+nx*ny*jz]+ 100.0) ;
                    //cout << cp[jx+nx*jy+nx*ny*jz] << endl;
                }
                else
                {
                    cp[jx+nx*jy+nx*ny*jz]=cs[jx+nx*jy+nx*ny*jz]+taup*(xire[jx+nx*jy+nx*ny*jz]-xicat[jx+nx*jy+nx*ny*jz]);
                }








		    if (t >= time1 && t < time2)
		    {


		        // Output involving 3D evolution
		        string everything = to_string(jx) + " " + to_string(jy) + " " + to_string(jz) + " " + to_string(cp[jx+nx*jy+nx*ny*jz]);
		        string enter = "\n";
		        out.write(&everything[0], everything.size());
		        out.write(&enter[0], enter.size());



		        // Output involving Spreading
		    	if (cp[jx+nx*jy+nx*ny*jz] >= thresh)
		        {
				  						
		            string everythingSpreading = to_string(jx) + " " + to_string(jy) + " " + to_string(jz);
		            string enterSpread = "\n";
		            spreadOut.write(&everythingSpreading[0], everythingSpreading.size());
		            spreadOut.write(&enterSpread[0], enterSpread.size());
		    
		        }
		        else
		        {

		            string everythingSpreading = to_string(-1) + " " + to_string(-1) + " " + to_string(-1);
		            string enterSpread = "\n";
		            spreadOut.write(&everythingSpreading[0], everythingSpreading.size());
		            spreadOut.write(&enterSpread[0], enterSpread.size());
			
		        }
	
	

		        // Output involving Activity Mapping
		        if (cp[jx+nx*jy+nx*ny*jz] >= thresh && testArray[fileCounter][jx+nx*jy+nx*ny*jz] == false)
		        {
			
		            for (int alpha = 2; alpha < fileCounter - 2; alpha++)
		            {

		                if (testArray[alpha][jx+nx*jy+nx*ny*jz] == true)
		                {

		                    trueCounter++;		

		                }

		            }       

		            if (trueCounter == 0)
		            {
		                    
		                testArray[fileCounter][jx+nx*jy+nx*ny*jz] = true;
		                string everythingSpreading = to_string(jx) + " " + to_string(jy) + " " + to_string(jz);
		                string enterSpread = "\n";
		                spreadActivityOut.write(&everythingSpreading[0], everythingSpreading.size());
		                spreadActivityOut.write(&enterSpread[0], enterSpread.size());		

		            }
		            else
		            {

		                string everythingSpreading = to_string(-1) + " " + to_string(-1) + " " + to_string(-1);
		                string enterSpread = "\n";
		                spreadActivityOut.write(&everythingSpreading[0], everythingSpreading.size());
		                spreadActivityOut.write(&enterSpread[0], enterSpread.size());

		            }

		            trueCounter = 0;

		        }

		        else 
		        {

		            string everythingSpreading = to_string(-1) + " " + to_string(-1) + " " + to_string(-1);
		            string enterSpread = "\n";
		            spreadActivityOut.write(&everythingSpreading[0], everythingSpreading.size());
		            spreadActivityOut.write(&enterSpread[0], enterSpread.size());

		        }


		    }





                if (cp[jx+nx*jy+nx*ny*jz] < 0.0) cp[jx+nx*jy+nx*ny*jz]=0.0;

                ccc0=cnsr[jx+nx*jy+nx*ny*jz];					//Calculate coupling current in NSR
                ccc1=cnsr[(jx-1)+nx*jy+nx*ny*jz];
                ccc2=cnsr[(jx+1)+nx*jy+nx*ny*jz];
                ccc3=cnsr[jx+nx*(jy-1)+nx*ny*jz];
                ccc4=cnsr[jx+nx*(jy+1)+nx*ny*jz];
                ccc5=cnsr[jx+nx*jy+nx*ny*(jz-1)];
                ccc6=cnsr[jx+nx*jy+nx*ny*(jz+1)];
                xicoupn[jx+nx*jy+nx*ny*jz]=couplingI(ccc0, ccc1, ccc2, ccc3, ccc4, ccc5, ccc6, xi*taunl, xi*taunt);			// xi is the coupling strength, which could be modified as need
                dotcnsr[jx+nx*jy+nx*ny*jz]=(xiup[jx+nx*jy+nx*ny*jz]-xileak[jx+nx*jy+nx*ny*jz])*Vi/Vnsr-xitr[jx+nx*jy+nx*ny*jz]*Vjsr/Vnsr+xicoupn[jx+nx*jy+nx*ny*jz];


                ccc0=ci[jx+nx*jy+nx*ny*jz];						//Calculate coupling current in cytosolic
                ccc1=ci[(jx-1)+nx*jy+nx*ny*jz];
                ccc2=ci[(jx+1)+nx*jy+nx*ny*jz];
                ccc3=ci[jx+nx*(jy-1)+nx*ny*jz];
                ccc4=ci[jx+nx*(jy+1)+nx*ny*jz];
                ccc5=ci[jx+nx*jy+nx*ny*(jz-1)];
                ccc6=ci[jx+nx*jy+nx*ny*(jz+1)];
                xicoupi[jx+nx*jy+nx*ny*jz]=couplingI(ccc0, ccc1, ccc2, ccc3, ccc4, ccc5, ccc6, xi*tauil, xi*tauit);
                dotci[jx+nx*jy+nx*ny*jz]=betaci[jx+nx*jy+nx*ny*jz]*(xisi[jx+nx*jy+nx*ny*jz]*(Vs/Vi)-xiup[jx+nx*jy+nx*ny*jz]+xileak[jx+nx*jy+nx*ny*jz]-xitci[jx+nx*jy+nx*ny*jz]+xicoupi[jx+nx*jy+nx*ny*jz]);

                ccc0=cs[jx+nx*jy+nx*ny*jz];						//Calculate coupling current in submembrane
                ccc1=cs[(jx-1)+nx*jy+nx*ny*jz];
                ccc2=cs[(jx+1)+nx*jy+nx*ny*jz];
                ccc3=cs[jx+nx*(jy-1)+nx*ny*jz];
                ccc4=cs[jx+nx*(jy+1)+nx*ny*jz];
                ccc5=cs[jx+nx*jy+nx*ny*(jz-1)];
                ccc6=cs[jx+nx*jy+nx*ny*(jz+1)];
                xicoups[jx+nx*jy+nx*ny*jz]=couplingI(ccc0, ccc1, ccc2, ccc3, ccc4, ccc5, ccc6, xi*tausl, xi*taust);


                cs[jx+nx*jy+nx*ny*jz]=(Vp*randomi[jx+nx*jy+nx*ny*jz]*cp[jx+nx*jy+nx*ny*jz]/taup+xinaca[jx+nx*jy+nx*ny*jz]*Vs+ci[jx+nx*jy+nx*ny*jz]*Vs/tausi-cs[jx+nx*jy+nx*ny*jz]*Vs+(ccc1+ccc2)*Vs/(xi*tausl)+(ccc3+ccc4)*Vs/(xi*taust)+(ccc5+ccc6)*Vs/(xi*taust))/(Vs/tausi+4.0*Vs/(xi*taust)+2.0*Vs/(xi*tausl)+Vp*randomi[jx+nx*jy+nx*ny*jz]/taup);
                if (cs[jx+nx*jy+nx*ny*jz] < 0.0) cs[jx+nx*jy+nx*ny*jz]=0.0;
                ci[jx+nx*jy+nx*ny*jz]=ci[jx+nx*jy+nx*ny*jz]+dotci[jx+nx*jy+nx*ny*jz]*dt;
                if (ci[jx+nx*jy+nx*ny*jz] < 0.0) ci[jx+nx*jy+nx*ny*jz]=0.0;
                cjsr[jx+nx*jy+nx*ny*jz]=cjsr[jx+nx*jy+nx*ny*jz]+dotcjsr[jx+nx*jy+nx*ny*jz]*dt;
                if (cjsr[jx+nx*jy+nx*ny*jz] < 0.0) cjsr[jx+nx*jy+nx*ny*jz]=0.0;
                cnsr[jx+nx*jy+nx*ny*jz]=cnsr[jx+nx*jy+nx*ny*jz]+dotcnsr[jx+nx*jy+nx*ny*jz]*dt;
                if (cnsr[jx+nx*jy+nx*ny*jz] < 0.0) cnsr[jx+nx*jy+nx*ny*jz]=0.0;
                cats[jx+nx*jy+nx*ny*jz]=cats[jx+nx*jy+nx*ny*jz]+xitcs[jx+nx*jy+nx*ny*jz]*dt;
                if (cats[jx+nx*jy+nx*ny*jz] < 0.0) cats[jx+nx*jy+nx*ny*jz]=0.0;
                cati[jx+nx*jy+nx*ny*jz]=cati[jx+nx*jy+nx*ny*jz]+xitci[jx+nx*jy+nx*ny*jz]*dt;
                if (cati[jx+nx*jy+nx*ny*jz] < 0.0) cati[jx+nx*jy+nx*ny*jz]=0.0;


              cit=cit+ci[jx+nx*jy+nx*ny*jz];
              xicatto=xicatto+randomi[jx+nx*jy+nx*ny*jz]*xicat[jx+nx*jy+nx*ny*jz];

		
            }

	counter = 0;

        t=t+dt;




        // Time range of data for particular iteration
        time1 = time1 + 0.1;
        time2 = time2 + 0.1;

        // Increment number of files.  Based on t and dt.
        fileCounter++;
	//cout << fileCounter << endl;


        // Update loading
        loadingVariable = (int)(t/dt/12);

        // Determine completion of simulation
        if(loadingVariable % checkNum == 0 && loadingVariable != 0)
        {

	  // cout << loadingVariable << "%" << " Complete" << endl;
	   cout << loadingVariable << endl;
            checkNum += 1;
            if(checkNum == 120 + 10)
            {

	      // cout << "Done!" << endl;

            }

        }




        if (t >= tperiod*(m+1))
        {

           m = m + 1;

        }					//Next period

    }												//time period loop ends up

	//cout << "hello" << endl;
	for (int num3 = 0; num3 < 1202; num3++)
	{

	    delete[] testArray[num3];	    

	}
	delete[] testArray;


    return (0);

}


