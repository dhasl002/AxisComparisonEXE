#ifndef SKELETONSA_H
#define SKELETONSA_H

#include "skeleton.h"
#include "utilityfunctions.h"
#include "energyfunction.h"
#include "protein.h"
#include "sideChainPacking_v1.h"


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// List Of Data Structures //////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//data structure to store the statistic for each conformation has been tested
struct statistic
{
	vector<EdgeStick> sticks;
	double intra;
	double inter;
	double RMSD;
	int pRank;
	int accepted; //1 yes 0 no

	//Initializer
	statistic () : intra (0.0), inter (0.0), RMSD(0), accepted(1), pRank(0) {} 
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
#define MAX_ITERATION 500


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// End Of Constants and Data Structures /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// List of Functions ////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void sort(vector <statistic> & A);												//sort the vector of statistics
long double probability(double curEnergy, double nEnergy,double temp);			//get the probability for SA
double getRMSD(vector <Protein> portions, Protein nativePDBFile, vector<EdgeStick> sticks);				//get the RMSD
vector<Protein> getNextState(vector<Protein> initialSkeleton,	
				  Protein nativePDBFile,
				  vector<EdgeStick> & sticks,
				  int currentDeltaTranslation[],
				  int currentDeltaRotation[],										//generate the next conformation
				  int seed);
vector<Protein> SA(vector<Protein> initialPortions, 
					Protein nativePDBFile, 
					vector<EdgeStick> sticks,
					vector<statistic> & stats);										//Simulated Annealing
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
		while ((j >= 0) && ((A[j].inter + A[j].intra) > (value.inter + value.intra)))
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

	long double boltzmanConstant = 1.38065 * pow(10,-23);
	long double deltaEnergy = -1 * (nEnergy - curEnergy);
//	cout<<" Deltaenergy = " <<deltaenergy<<" Temp = "<<temp<<endl;
//	cout<<"above e is : "<<deltaenergy/temp<<"  The probability = "<<exp(deltaenergy/temp)<<endl;
	return exp(deltaEnergy/temp);
	//return pow(2.718281828,deltaenergy/(temp*boltzmanconstant));
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double getRMSD(vector <Protein> portions, Protein nativePDBFile, vector<EdgeStick> sticks)
{
	////////////get RMSD with native
	double rmsd = 0;
	Coordinate p1,p2;
	double length=0;


	for (int i=0;i<sticks.size();i++)
	{

		int numOfAA = portions[i].numOfAA ();

		length += numOfAA;

		int walkDirection;
		if (sticks[i].direction == 1)
			walkDirection = 1;
		else
			walkDirection = -1;

		int startIndx = sticks[i].startAAIndx + sticks[i].shift;

		for (int j=0; j<numOfAA; j++)
		{
			//cout<<portions[i].AAs[j].chr3<<" "<<portions[i].AAs[j].num<<" pdbfile startAA = "<<nativePDBFile.AAs [startIndx].chr3	\
				<<"  "<<nativePDBFile.AAs [startIndx].num<<" shift= "<<sticks[i].shift<<"  direction= "<<sticks[i].direction<<endl;
			p1 = portions[i].getAtomCoordinate(j, " CA ");
			p2 = nativePDBFile.getAtomCoordinate(startIndx, " CA ");

			startIndx += walkDirection;
			rmsd += ((p1.x - p2.x) * (p1.x - p2.x)) + ((p1.y - p2.y) * (p1.y - p2.y)) + ((p1.z - p2.z) * (p1.z - p2.z));

		}
		//cout<<"=========="<<endl;
	}
	return sqrt(rmsd/length);

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<Protein> getNextState(vector<Protein> initialSkeleton, 
				Protein nativePDBFile, 
				vector<EdgeStick> & sticks, 
				int currentDeltaTranslation[],
				int currentDeltaRotation[],
				int seed)
{
	int success = 0,
		randShift,
		randTranslation,
		randRotation,
		i;

	Protein allSticks;					//the protein file where the all portions will be connected as one file
	vector<Protein> nextStates;

	srand((seed * 2) + 1);


	nextStates = initialSkeleton;
	
	for (i=0;i <sticks.size(); i++)
	{
		//set the random shift
		randShift = getRandom(-MAX_SHIFT_ALLOWED, MAX_SHIFT_ALLOWED);

		//check if the shift is valid
		while ((randShift > sticks[i].validityRight) || (-randShift > sticks[i].validityLeft))
			randShift = getRandom(-MAX_SHIFT_ALLOWED, MAX_SHIFT_ALLOWED);

		//store the shift value
		sticks[i].shift = randShift;

		//set the random translation
		randTranslation = getRandom(- (MAX_TRANSLATION_ALLOWED + currentDeltaTranslation[i]), MAX_TRANSLATION_ALLOWED - currentDeltaTranslation[i]);
		sticks[i].translation = currentDeltaTranslation[i] + randTranslation;

		//set Random rotation 
		randRotation = getRandom(- MAX_ROTATION_ALLOWED, MAX_ROTATION_ALLOWED);
		sticks[i].rotation = (randRotation + currentDeltaRotation[i]) % 360;

		//wite on the screen some of the random information
		//cout<<"direction = "<<sticks[i].direction <<" shift = "<<randShift<<"  translation = "<<randTranslation<<"  rotation = "<<randRotation<<endl; 

		//assign sequence to sticks hlces
		AssignSeqToSkeleton(nativePDBFile, nextStates[i], randShift, i, sticks);

		//translate sticks hlces
		nextStates[i].translateBy (sticks[i].translation);

		//rotate sticks hlces around axis by rand rotation
		rotateAroundAxis(nextStates[i], sticks[i].rotation);

		//append all sticks in one protein structure for side chain packing purposes
		allSticks.append (nextStates[i],0,nextStates[i].numOfAA()-1);
	}


	// predict Side chains
	predictSideChain(allSticks);

	//split back the structure to initial structures
	double size = 0;
	int start = 1;
	for (i=0;i<initialSkeleton.size();i++)
	{
		size += initialSkeleton[i].numOfAA();
		initialSkeleton[i].AAs.clear();
		initialSkeleton[i].append(allSticks,start-1,size-1);
		start += initialSkeleton[i].numOfAA();
	}

	return initialSkeleton;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<Protein> SA(vector<Protein> initialPortions, 
				   Protein nativePDBFile, 
				   vector<EdgeStick> sticks, 
				   vector<statistic> & stats)
{
	vector<Protein> neighborState,
					currentState,
					bestState;

	int iteration = 0,
		totalnAA,
		nAA,
		i;

	double currentEnergy = 100,
			neighborEnergy,
			bestEnergy = 10000,
			temp = MAX_TEMP;

	int *currentDeltaTranslation;		//the current tranlsation shift from the initial Skeleton (delta translation from initial skeleton)
	int *currentDeltaRotation;			//the current rotaion angle from the initial skeleton orientation	(delta rotation from initial skeleton)

	//initiate the pointer of translation to MAX_TRANSLATION_ALLOWED
	currentDeltaTranslation = new int[initialPortions.size()];

	//initiate the pointer of rotation change
	currentDeltaRotation = new int[initialPortions.size()];

	for (i=0; i<initialPortions.size(); i++)
	{
		currentDeltaTranslation[i]  = 0;
		currentDeltaRotation[i]		= 0;
	}

	currentState = initialPortions;

	statistic tmpStat;		//tmp variable to store statistics
	

	while ((iteration < MAX_ITERATION) && (temp > MIN_TEMP) && (currentEnergy > MIN_ENERGY))
	{
		//get the neighbor state for the current conformation
		neighborState = getNextState(initialPortions,				//initial skeleton vector that contains all hlces
									 nativePDBFile,					//native pdb file
									 sticks,						//data structure that contains information about sticks
									 currentDeltaTranslation,		//the amount of translation allowed from the left
									 currentDeltaRotation,			//the delta of rotation occured on initial skeleton
									 iteration);					//seed

		tmpStat.sticks = sticks;		//store neighbor state (structure) information


		//get inter energy for the portions
		long double interE = interEnergy(neighborState);

		long double intraE = 0;

		totalnAA = 0;	//total number of AA would be included in the energy calculations (intra energy)
		//get the intra energy for each portion
		for (i=0;i<neighborState.size();i++)
		{
			intraE += intraEnergy(neighborState[i],nAA);
			totalnAA += nAA;
		}

		intraE /= totalnAA;

		//the energy of the neighbor state
		neighborEnergy = (intraE + interE) / 2;

		//store the information about energy function
		tmpStat.inter = interE;
		tmpStat.intra = intraE;
		
		tmpStat.RMSD = getRMSD(neighborState, nativePDBFile, sticks);		//get the RMSD between this structure and the native one


		//save the information of neighbor conformation (the one just has been generated)
		stats.push_back(tmpStat);	

		if (neighborEnergy < currentEnergy)
		{
			currentState  = neighborState;
			currentEnergy = neighborEnergy;

			//commit the new configuration
			for(i=0; i<sticks.size(); i++)
			{
				currentDeltaTranslation[i] = sticks[i].translation;		
				currentDeltaRotation[i]	   = sticks[i].rotation;
			}
		}
		else
			if (probability(currentEnergy, neighborEnergy, temp) > (getRandom(1,1000)/1000.00))
			{
				//cout<<"neighbor is accepted..."<<endl;
				currentState = neighborState;
				currentEnergy = neighborEnergy;

				//commit the new configuration
				for(i=0; i<sticks.size(); i++)
				{
					currentDeltaTranslation[i] = sticks[i].translation;
					currentDeltaRotation[i]    = sticks[i].rotation;
				}
			}

		//save best state
		if (neighborEnergy < bestEnergy)
			bestState = neighborState;

		//cooling step
		temp -= temp * COOLING_CONSTANT;

		//iteration progress
		iteration++;
	}



	//sort statistics according to the Multi-Peak Energy Value
	sort(stats);


	return bestState;
	
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// End Of IMPLEMENTATION ////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////				

#endif
