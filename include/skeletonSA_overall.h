#ifndef SKELETONSA_H
#define SKELETONSA_H

#ifdef _WIN				//if work under Windows
#include <direct.h>
#else					//if work under Linux
#include <sys/stat.h>
#endif

#include "skeleton_overall.h"
#include "utilityfunctions.h"
#include "energyfunction.h"
#include "protein.h"
#include "sideChainPacking_v1.h"
#include "topology_graph.h"
#include "rmsd.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// List Of Data Structures //////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//data structure to store the statistic for each conformation has been tested
struct statistic
{
	vector<vector<short> > topology;
	double intra;
	double RMSDca;		//RMSD of the Ca of the full model
	double RMSDSSE;		//RMSD of the SSE (Ca)
	double RMSDbb;		//RMSD of the back bone for the full model
	vector<short> shift;
	vector<short> translation;
	int pRank;
	int accepted; //1 yes 0 no

	//Initializer
	statistic () : intra (0.0), RMSDca(0), RMSDSSE(0), RMSDbb(0), accepted(1), pRank(0) {}
};

//the data structure to store all conformations tested so far for a given assignemtn (permutation)
struct statsRecord
{
	int topologyIndx;
	vector<statistic> stats;
};

//data structure that stores all statistics for one assignment
//vector<statistic> stats;


#define COOLING_CONSTANT 0.0005
#define MAX_TEMP 2000
#define MIN_TEMP 1
#define MIN_ENERGY -30
#define MAX_ITERATION 100


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// End Of Constants and Data Structures /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// List of Functions ////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void sort(vector <statistic> & A);												//sort the vector of statistics
long double probability(double curEnergy, double nEnergy,double temp);			//get the probability for SA
double getRMSDca(Protein &portion, Protein &nativePDBFile, short sIndx);				//get the RMSD between the two full models
double getRMSDBackBone(Protein &portion, Protein &nativePDBFile, short sIndx);				//get the RMSD between the two full models


/*
vector<Protein> getNextState(vector<Protein> initialSkeleton,
				  Protein nativePDBFile,
				  vector<EdgeStick> & sticks,
				  int currentDeltaTranslation[],
				  int currentDeltaRotation[],										//generate the next conformation
				  int seed);
vector<Protein> SA(vector<Protein> initialPortions,					//initial skeleton
					Protein nativePDBFile,							//native pdb file..... used to calculate the RMSD
					vector<EdgeStick> sticks,						//data structure to store topology-related information
					vector<statistic> & stats,
					int numOfRecKeep = 10);							//number of records (best records) we wanna keep for statistics purposes..stored in stats
					//Simulated Annealing
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// ////////// End Of The List //////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// IMPLEMENTATION ///////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//insertion sort
void sort(vector <statistic> & A)
{
	for(int i=0;i<A.size();i++)
	{
		statistic value;
		value = A[i];
		int j = i-1;
		while ((j >= 0) && (A[j].intra > value.intra))
		{
			A[j + 1] = A[j];
			j = j-1;
		}
		A[j+1] = value;
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
long double probability(double curEnergy, double nEnergy,double temp)
{

//	long double boltzmanConstant = 1.38065 * pow(10,-23);
	long double deltaEnergy = -1 * (nEnergy - curEnergy);
//	cout<<" Deltaenergy = " <<deltaenergy<<" Temp = "<<temp<<endl;
//	cout<<"above e is : "<<deltaenergy/temp<<"  The probability = "<<exp(deltaenergy/temp)<<endl;
	return exp(deltaEnergy/temp);
	//return pow(2.718281828,deltaenergy/(temp*boltzmanconstant));
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double getRMSDca(Protein &portion, Protein &nativePDBFile, short sIndx)
{
	////////////get RMSD with native
	double rmsd = 0;
	Coordinate p1,p2;
	double length=0;

	//nX3 coordinates matrices
	double **native, **model;

	//other data structure needed in calculating RMSD...see rmsd for details
	double mov_com[3], mov_to_ref[3], U[3][3];

    int numOfAA = portion.numOfAA ();


	//resize array of coordinates
	native = AllocateDynamicArray<double> (numOfAA,3);
	model =  AllocateDynamicArray<double> (numOfAA,3);

    for (int j=0; j<numOfAA; j++)
    {
        //cout<<portions[i].AAs[j].chr3<<" "<<portions[i].AAs[j].num<<" pdbfile startAA = "<<nativePDBFile.AAs [sIndx].chr3	\
            <<"  "<<nativePDBFile.AAs [sIndx].num<<" shift= "<<sticks[i].shift<<"  direction= "<<sticks[i].direction<<endl;
        p1 = portion.getAtomCoordinate(j, " CA ");
        p2 = nativePDBFile.getAtomCoordinate(sIndx, " CA ");

        sIndx ++;

		native[j][0] = p2.x;  native[j][1] = p2.y;  native[j][2] = p2.z;
		model[j][0]  = p1.x;  model[j][1]  = p1.y;  model[j][2]  = p1.z;

        //rmsd += ((p1.x - p2.x) * (p1.x - p2.x)) + ((p1.y - p2.y) * (p1.y - p2.y)) + ((p1.z - p2.z) * (p1.z - p2.z));
    }

	//one way of calculating RMSD....finding the rotation matrix first
	calculate_rotation_rmsd (native, model, numOfAA, mov_com, mov_to_ref, U, &rmsd);

	//another way of calculating the RMSD....without finding the rotation matrix....see rmsd.h for more details
	//fast_rmsd(native, model, numOfAA, rmsd);

	//return sqrt(rmsd/numOfAA);
	return rmsd;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double getRMSDBackBone(Protein &portion, Protein &nativePDBFile, short sIndx)
{
	////////////get RMSD with native
	double rmsd = 0;
	Coordinate N,Ca, C, O;
	double length=0;

	//nX3 coordinates matrices
	double **native, **model;

	//other data structure needed in calculating RMSD...see rmsd for details
	double mov_com[3], mov_to_ref[3], U[3][3];

    int numOfAA = portion.numOfAA ();


	//resize array of coordinates
	native = AllocateDynamicArray<double> (numOfAA*4,3);
	model =  AllocateDynamicArray<double> (numOfAA*4,3);

    for (int j=0; j<numOfAA; j++)
    {
        //cout<<portions[i].AAs[j].chr3<<" "<<portions[i].AAs[j].num<<" pdbfile startAA = "<<nativePDBFile.AAs [sIndx].chr3	\
            <<"  "<<nativePDBFile.AAs [sIndx].num<<" shift= "<<sticks[i].shift<<"  direction= "<<sticks[i].direction<<endl;
		N  = portion.getAtomCoordinate(j, " N ");
        Ca = portion.getAtomCoordinate(j, " CA ");
		C  = portion.getAtomCoordinate(j, " C ");
		O  = portion.getAtomCoordinate(j, " O ");

		model[4*j][0]    = N.x;   model[4*j][1]    = N.y;   model[4*j][2]    = N.z;
		model[4*j+1][0]  = Ca.x;  model[4*j+1][1]  = Ca.y;  model[4*j+1][2]  = Ca.z;
		model[4*j+2][0]  = C.x;   model[4*j+2][1]  = C.y;   model[4*j+2][2]  = C.z;
		model[4*j+3][0]  = O.x;   model[4*j+3][1]  = O.y;   model[4*j+3][2]  = O.z;

        N  = nativePDBFile.getAtomCoordinate(sIndx, " N ");
		Ca = nativePDBFile.getAtomCoordinate(sIndx, " CA ");
		C  = nativePDBFile.getAtomCoordinate(sIndx, " C ");
		O  = nativePDBFile.getAtomCoordinate(sIndx, " O ");

		native[4*j][0]   = N.x;   native[4*j][1]   = N.y;   native[4*j][2]   = N.z;
		native[4*j+1][0] = Ca.x;  native[4*j+1][1] = Ca.y;  native[4*j+1][2] = Ca.z;
		native[4*j+2][0] = C.x;   native[4*j+2][1] = C.y;   native[4*j+2][2] = C.z;
		native[4*j+3][0] = O.x;   native[4*j+3][1] = O.y;   native[4*j+3][2] = O.z;

		sIndx ++;
		

        //rmsd += ((p1.x - p2.x) * (p1.x - p2.x)) + ((p1.y - p2.y) * (p1.y - p2.y)) + ((p1.z - p2.z) * (p1.z - p2.z));
    }

	//one way of calculating RMSD....finding the rotation matrix first
	calculate_rotation_rmsd (native, model, numOfAA*4, mov_com, mov_to_ref, U, &rmsd);

	//another way of calculating the RMSD....without finding the rotation matrix....see rmsd.h for more details
	//fast_rmsd(native, model, numOfAA, rmsd);

	//return sqrt(rmsd/numOfAA);
	return rmsd;

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Protein SA(Protein nativePDBFile,
				   vector<vector<cell> > &graph,
				   vector<vector< Coordinate> > &hEdges,
				   Coordinate **traceList,
				   float *traceLength,
				   short *traceNoOfPoints,
				   short **topology,
				   vector<statistic> & stats,
				   string outDir,
				   bool printStructures,
				   int myRank)
{
	Protein neighborState,
            currentState,
            bestState;

	int iteration = 0,
		nAA,
		i,
		nSticks = hEdges.size();

	double currentEnergy = 100,
			neighborEnergy,
			bestEnergy = 10000,
			temp = MAX_TEMP,
			rmsdSSE;

	short *translations, *shifts;
	short startIndx;

	//initiate the pointer of translation to MAX_TRANSLATION_ALLOWED
	translations = new short [nSticks];
	shifts = new short [nSticks];


	for (i=0; i<nSticks; i++)
	{
		translations[i] = 0;
		shifts[i] = 0;
	}

	//currentState = initialPortions;

	statistic tmpStat;		//tmp variable to store statistics

	//tmpStat.topology = new short *[nSticks];
	//tmpStat.shift = new short [nSticks];
	//tmpStat.translation = new short [nSticks];

	tmpStat.topology.resize(nSticks);
	tmpStat.shift.resize(nSticks, 0);
	tmpStat.translation.resize(nSticks, 0);
	for (i=0;i<nSticks;i++)
        //tmpStat.topology[i] = new short [2];
        tmpStat.topology[i].resize(2,0);

	//ensure that stats is empty
	stats.clear();
	if (printStructures){
	    outDir += "/NativeConformations/";
#ifdef _WIN
	mkdir(outDir.c_str ());
#else
	mode_t mode = 0740;
	mkdir(outDir.c_str (), mode);
#endif
	}

    ofstream out;
    string outfile;
    outfile = outDir + "/";
    outfile += "proc_";
    outfile += toString(myRank);
    outfile += ".txt";
    out.open(outfile.c_str());
    if (!out) {
        errMsg("SkeletonSA_overall", "SA", "Unable to open " + outfile, true);
    }
    out<<"Working on topology : "<<endl;
    for (i=0; i<nSticks; i++){
        out<<"["<<topology[i][0]<<" "<<topology[i][1]<<"] ";
    }
    out<<endl<<endl;


	while ((iteration < MAX_ITERATION) && (temp > MIN_TEMP) && (currentEnergy > MIN_ENERGY))
	{
	    vector<vector<Coordinate> > newEdges = hEdges;

	    out<<"iteration # "<<iteration+1<<endl;
	    //get random shift and translation to apply to the original sequence and sticks.
	    getRandomList(-2, 2, shifts, nSticks, (myRank+1) * (iteration+1));
	    getRandomList(-MAX_TRANSLATION_ALLOWED, MAX_TRANSLATION_ALLOWED, translations, nSticks, (myRank+1+iteration)*(iteration+1));

	    //apply the translation on the sticks and save information
	    for (i=0; i<nSticks;i++){
            tmpStat.topology[i][0] = topology[i][0];
            tmpStat.topology[i][1] = topology[i][1];
            tmpStat.shift[i] = shifts[i];
            tmpStat.translation[i] = translations[i];

	        //cout<<"shift applied on "<<i+1<<" is "<<shifts[i]<<" translation applied= "<<translations[i]<<endl;
	        out<<"shift applied on "<<i+1<<" is "<<tmpStat.shift[i]<<" translation applied= "<<tmpStat.translation[i]<<endl;
	        if (translations[i] != 0){
	            translateEdge(newEdges[i], translations[i]);
	        }

	    }
        //build the model
		Protein SSEStructure;
        neighborState = buildFullModel(nativePDBFile, graph, newEdges, topology, traceList, traceLength, traceNoOfPoints, shifts, startIndx, rmsdSSE, SSEStructure);

		if (printStructures)
			SSEStructure.writePDB(outDir + "SSE_"+toString(rmsdSSE) + ".pdb",1, SSEStructure.numOfAA());

        //set torsion angles
        for (i=0; i<neighborState.numOfAA(); i++){
            neighborState.AAs[i].angles.phi = neighborState.getTorsion("phi",i);
            neighborState.AAs[i].angles.psi = neighborState.getTorsion("psi",i);
        }
        // predict Side chains
        predictSideChain(neighborState);

        //copy the topology information
        /*
        for (i=0;i<nSticks;i++){
            tmpStat.topology[i][0] = topology[i][0];
            tmpStat.topology[i][1] = topology[i][1];
            tmpStat.shift[i] = shifts[i];
            tmpStat.translation[i] = translations[i];
        }
        */

		//get inter energy for the portions
		//long double interE = interEnergy(neighborState);

		long double intraE = 0;


		//get the intra energy for each portion
		intraE = intraEnergy(neighborState,nAA);
		intraE /= nAA;

		//the energy of the neighbor state
		neighborEnergy = intraE;

		//store the information about energy function
		tmpStat.intra = intraE;

        //cout<<"finding RMSD.."<<endl;
		tmpStat.RMSDca = getRMSDca(neighborState, nativePDBFile, startIndx);		//get the RMSD between this structure and the native one
		tmpStat.RMSDbb = getRMSDBackBone(neighborState, nativePDBFile, startIndx);	//get the back bone RMSD b/w the two structures
		tmpStat.RMSDSSE = rmsdSSE;

        out<<"   Energy : "<<tmpStat.intra<<"\t";
        out<<"   RMSDca : "<<tmpStat.RMSDca<<"\t";
		out<<"   RMSDbb : "<<tmpStat.RMSDbb<<"\t";
		out<<"   RMSDSSE  :"<<tmpStat.RMSDSSE<<endl;
        if (printStructures){
            neighborState.writePDB(outDir + "RMSD"+toString(tmpStat.RMSDca)+"_E" + toString(tmpStat.intra)+".pdb", 1, neighborState.numOfAA());
        }
		//save the information of neighbor conformation (the one just has been generated)
		stats.push_back(tmpStat);

		if (neighborEnergy < currentEnergy)
		{
			currentState  = neighborState;
			currentEnergy = neighborEnergy;
		}
		else
			if (probability(currentEnergy, neighborEnergy, temp) > (getRandom(1,1000)/1000.00))
			{
				//cout<<"neighbor is accepted..."<<endl;
				currentState = neighborState;
				currentEnergy = neighborEnergy;
			}

		//save best state
		if (neighborEnergy < bestEnergy)
			bestState = neighborState;

		//cooling step
		temp -= temp * COOLING_CONSTANT;

		//iteration progress
		iteration++;
	}


    out.close();

	//sort statistics according to the Multi-Peak Energy Value
	sort(stats);


	//this is a memory expensive way...we need to imporove it when the method is stable

	return bestState;

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// End Of IMPLEMENTATION ////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
