#ifndef SKELETON_H
#define SKELETON_H

#include <sstream>

#include "geometry.h"
#include "protein.h"
#include "constants.h"
#include "utilityfunctions.h"



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// List Of Data Structures //////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct EdgeStick
{
	Coordinate start;			//start point (Coordinate) from the edge file
	Coordinate end;				//end point (Coordinate) from edge file
	int stickNum;				//stick number as in edge file
	int hlxAssigned;			//the indx of hlx assigned to this stick (indx as in pdb.hlces)
	int direction;				//the direction of the sequence (how is the sequence assigned to the stick) 1 forward  0 reverese
	int length;					//the assumed number of AA's
	int startAAIndx;			//the indx of first AA (AAs indx)
	int shift;					//the amound of shift from the original hlx
	int translation;			//the amound of translation
	int validityRight;			//the maximum size of shift could be applied to the sequence to the right...
	int validityLeft;			//the maximum size of shift could be applied to the sequence to left...

	double rotation;			//the amount of rotation from initial skeleton


	//initializer
	EdgeStick(): stickNum(-1), hlxAssigned(-1), direction(1), length(0), startAAIndx(-1), \
		shift(0), translation(0), rotation(0), validityRight(0), validityLeft(0) {}
};

static string initialBB = "initialbb.pdb"; //the template hlx
#define MAX_SHIFT_ALLOWED	2
#define MAX_TRANSLATION_ALLOWED 0
#define MAX_ROTATION_ALLOWED 30

//used to save the sticks given by edge file...Global Variable
//static vector<EdgeStick> sticks;
//data structures that stores all sticks vectors (all topologies)
struct sticksRec
{
	vector<EdgeStick> sticks;
};

vector<sticksRec> sticksVect;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// End Of Constants and Data Structures /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////// list of functions /////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void printSkeleton(vector<Protein> hlces, vector<EdgeStick> sticks, string path, string TAG = "");						//print skeleton as pdb files
vector<EdgeStick> getEdgeSticks(string mrcFile);														//build the sticks vector that contains the first and last points of each stick
void overlapSticks(vector<Protein> & sticksHlces, Protein pdbfile, vector<EdgeStick> sticks);			//translate the helices to overlap the sticks given by Edge file those stored in sticks data structure
vector<Protein> getInitialSkeleton(Protein nativePDBFile, string mrcFile, vector<EdgeStick> & sticks, string outDir);	// builds initial skeleton for the given edge file and native pdb file
void setShiftValidity(Protein nativePDBFile, vector<EdgeStick> & sticks);								//set the amound of shift to left and right of the sticks hlces allowed
void AssignSeqToSkeleton(Protein nativePDBFile, Protein & portion, int shift, int stickIndx, vector<EdgeStick> sticks);	//assign AA seq from native pfb file to the sticks hlces
void rotateAroundAxis(Protein &portion, double angleDegree);											//rotate a portion (hlx mainly) around its axis		
void generateAllTopologies(vector<EdgeStick> sticks, vector<sticksRec> & allTopologies);	//given one permutation stored in sticks, then generate all possible topologies and store it in allTopologies data structure
double getDistancebwSticks(vector<Protein> skeleton, int indx1,int direction1, int indx2, int direction2);	//get the distance between two skeleton hlces
inline
void generateValidTopologies(Protein nativePDBFile, vector<Protein> initialSkeleton, vector<int> permutation, vector<EdgeStick> sticks);
void genrateAllValidTopologies(Protein nativePDBFile, vector<Protein> initialSkeleton, vector<OnePermutation> allPermutations, vector<EdgeStick> sticks); //given the pdb file and initial stick information and all permutations, generate all valid topologies into sticksVect


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// End Of The List //////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// IMPLEMENTATION ///////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////print all skeleton hlces
void printSkeleton(vector<Protein> hlces, vector<EdgeStick> sticks, string path, string TAG)
{

	int i;
	
	for (i=0;i<sticks.size();i++)
	{

		stringstream intStr,intStr2;
		string outFileName = path + "/";
		outFileName += "stick_"  ;//+ "_hlx_" + sticks[i].hlxassigned + "_" + TAG + ".pdb";
		intStr<<sticks[i].stickNum;
		//strcat(outfilename,intstr.str());
		outFileName  += intStr.str();
		outFileName  += "_hlx_";
		intStr2<<sticks[i].hlxAssigned + 1;
		outFileName  += intStr2.str();
		outFileName += "_";
		outFileName += TAG;
		outFileName += ".pdb";

		hlces[i].writePDB (outFileName,1,hlces[i].numOfAA());
		//cout<<"Stick number "<<sticks[i].stickNum<<"  assigned to hlx number  "<<endl;
	}

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//build and return the sticks vector that contains the first and last points of each stick
vector<EdgeStick> getEdgeSticks(string mrcFile)
{
	ifstream inFile;
	string str;
	EdgeStick tmpStick;

	inFile.open(mrcFile.c_str());
	if (!inFile)
	{
		cout<<"=================================== in skeleton::setMRCSticks(string) ============="<<endl;
		cerr << "Unable to open " << mrcFile << endl; 
		cout<<"==================================================================================="<<endl;
		exit(1);
	}
	else
	{
		vector<EdgeStick> sticks;
		tmpStick.stickNum = -1;
		while (!inFile.eof()){
			getline(inFile,str);
			if (str.substr(0,1)==" "){
				int stickNum = atoi(str.substr(42,3).c_str());
				if (stickNum != tmpStick.stickNum){
					tmpStick.stickNum = stickNum;
					tmpStick.start.x  = atof(str.substr(0,13).c_str());
					tmpStick.start.y  = atof(str.substr(14,13).c_str());
					tmpStick.start.z  = atof(str.substr(28,13).c_str());
					sticks.push_back(tmpStick);
				}
				else {
					sticks[sticks.size()-1].end.x  = atof(str.substr(0,13).c_str());
					sticks[sticks.size()-1].end.y  = atof(str.substr(14,13).c_str());
					sticks[sticks.size()-1].end.z  = atof(str.substr(28,13).c_str());
				}
			}
		}

		inFile.close();
		return sticks;
	}
	

//	To check the correctness of coordinates exctracted from mrc file
/*
	for (int i=0;i<sticks.size();i++)
	{
		cout<<"Stick# "<<sticks[i].stickNum<<endl;
		Vectors vStart(sticks[i].start), vEnd(sticks[i].end);
		cout<<"  Start : ";
		vStart.print ();
		cout<<"  End   : ";
		vEnd.print();
	}
*/
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//translate the helices to overlap the sticks given by Edge file those stored in sticks data structure
void overlapSticks(vector<Protein> & sticksHlces, Protein pdbfile, vector<EdgeStick> sticks)
{

	Coordinate pnt1,pnt2,pnt3,pnt4;

	for (int i=0;i<sticksHlces.size();i++)
	{

		pnt1 = sticks[i].start;
		pnt2 = sticks[i].end;

		Coordinate startEndCoord, endEndCoord;
		int tmpNumOfAA = sticksHlces[i].numOfAA();
		int AAsCounter = 0;	//for atoms will be involved in the center axial
		//take at most 3 atoms from each end
		while ((AAsCounter < tmpNumOfAA - 1) && (AAsCounter < 3))
		{
			Coordinate tmpCoord;

			//get coordinate of Ca atoms from the start end
			tmpCoord = sticksHlces[i].getAtomCoordinate(AAsCounter," CA ");	
			startEndCoord.x += tmpCoord.x;
			startEndCoord.y += tmpCoord.y;
			startEndCoord.z += tmpCoord.z;

			//get coordinate of Ca atoms from the end end
			tmpCoord = sticksHlces[i].getAtomCoordinate(tmpNumOfAA - AAsCounter - 1," CA ");
			endEndCoord.x += tmpCoord.x;
			endEndCoord.y += tmpCoord.y;
			endEndCoord.z += tmpCoord.z;

			AAsCounter++;
		}

		//find the average (center)
		startEndCoord.x = startEndCoord.x / AAsCounter;
		startEndCoord.y = startEndCoord.y / AAsCounter;
		startEndCoord.z = startEndCoord.z / AAsCounter;

		endEndCoord.x = endEndCoord.x / AAsCounter;
		endEndCoord.y = endEndCoord.y / AAsCounter;
		endEndCoord.z = endEndCoord.z / AAsCounter;	

/*
		Dr. Weitao's method of finding the center
		to get getCenter refer to skeleton program dated 10_30_2008

		double Radius;
		pnt3 = getCenter(sticksHlces[i].getAtomCoordinate(0," CA "), 
						 sticksHlces[i].getAtomCoordinate(1," CA "),
						 sticksHlces[i].getAtomCoordinate(2," CA "),
						 Radius);

		int hlxNumOfAA = sticksHlces[i].numOfAA();
		pnt4 = getCenter(sticksHlces[i].getAtomCoordinate(hlxNumOfAA-1," CA "),
						 sticksHlces[i].getAtomCoordinate(hlxNumOfAA-2," CA "),
						 sticksHlces[i].getAtomCoordinate(hlxNumOfAA-3," CA "),
						 Radius);
*/		

		//move the hlx to overlap the stick
		sticksHlces[i].overlapLine(startEndCoord,endEndCoord,pnt1,pnt2);

		//find the 2 correspponding points represent the imiginary Ca atoms @ the line ... from any end..here from the end (could be from start end)
		Coordinate imiginaryCaLast = pointLineIntersection(sticksHlces[i].getAtomCoordinate(tmpNumOfAA-1," CA "),pnt1,pnt2);		//for last Ca

		Coordinate imiginaryCaFirst = pointLineIntersection(sticksHlces[i].getAtomCoordinate(0," CA "),pnt1,pnt2);			//for first Ca

		double distance = getDistance(imiginaryCaLast,pnt2);	//find the distance b/w the end end of the hlx and the stick
		distance += getDistance(imiginaryCaFirst,pnt1);			//find the distance b/w the first end of the hlx and the stick

		distance /= ALPHA_RISE;									//determine the number of Ca Rise distance hlx is shifted away from the stick

		sticksHlces[i].translateBy((distance/2 + 0.5));		//translate the hlx to overlap the stick from the center
	}
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// builds initial skeleton for the given edge file and native pdb file
//the hlces vector for the native pdb file should be built before calling this function
vector<Protein> getInitialSkeleton(Protein nativePDBFile, string mrcFile, vector<EdgeStick> & sticks, string outDir)
{
	
	Protein initial(initialBB),		//get initial backbone for hlces
			tmpSkeleton;
	vector<Protein> InitialSkeleton;
	int i;

	//The pdb file should contain hlces and the number of sticks should be greater than or equal the number of hlces
	if ((nativePDBFile.hlces.size()) && (sticks.size() <= nativePDBFile.hlces.size()))
	{
		
		for (i=0;i<sticks.size();i++)
		{
			double length = getDistance(sticks[i].start, sticks[i].end);	//get the length of the stick in the edge file
			
			//The suggested length of the hlx (The number of AA in the hlx)
			length = int ((length/ALPHA_RISE) + 0.5 /* To round it */); // Here the approx. number of aa per sticks..u can also divide by 1.33 which I think more accurate
			
			sticks[i].length = length;		//number of AA's in the initial hlces

			Protein tmp;

			tmp.append(initial,0,length-1);
			
			InitialSkeleton.push_back(tmp);
			//InitialSkeleton[i].writeatoms("tmp.pdb");

			//cout<<"stick : "<<sticks[i].stickNum <<" expected # of AA's : "<<sticks[i].length<<endl;

		}

		
		//Transforms the initial helices to overlap the sticks...the assigning is done randomly
		overlapSticks(InitialSkeleton, nativePDBFile, sticks);

		
		//All helices in one file
		for (i=0;i<InitialSkeleton.size();i++)
			tmpSkeleton.append(InitialSkeleton[i],0,InitialSkeleton[i].numOfAA()-1);
		
		string initialSkeletonFileName = outDir + "/";
		initialSkeletonFileName += "InitialSkeleton_" + nativePDBFile.ID + ".pdb";
		tmpSkeleton.writePDB(initialSkeletonFileName,1, tmpSkeleton.numOfAA());
		
		//write sticks as pdb files...with N and CA atoms to be viewed using chimera
		for (i=0;i<InitialSkeleton.size(); i++)
		{
			Protein tmpHlx;
			string stickName = outDir + "/";
			stickName += nativePDBFile.ID;
			stringstream int2Str;
			int2Str<< sticks[i].stickNum;
			tmpHlx.append(InitialSkeleton[i],0,0);
			tmpHlx.removeHAtoms();
			tmpHlx.AAs [0].atoms[0].coord = sticks[i].start;		//N atom now equal to the start
			tmpHlx.AAs [0].atoms[tmpHlx.getAtomIndx(0," CA ")].coord = sticks[i].end;		//CA atom now equal to end
			tmpHlx.AAs [0].num = sticks[i].stickNum;
			tmpHlx.AAs[0].atoms .erase(tmpHlx.AAs[0].atoms .begin () + 2, tmpHlx.AAs[0].atoms.end());					//remove all other atoms
			stickName += "_stick_";
			stickName += int2Str.str();
			stickName += ".pdb";
			tmpHlx.writePDB (stickName, 1, 1);
		}
		

		return InitialSkeleton;
	}
	else
	{
		cout<<"=============================== In skeleton::SetInitialSkeleton ==============="<<endl;
		cout<<"The native pdb contains no hlces or the number of sticks exceeds the number of hlces.."<<endl;
		cout<<"==============================================================================="<<endl;
		exit(0);
	}
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// hlces should be assigned to sticks first...sticks.hlxassigned should be set before
//this function is supposed to use to determine the maximum of shift could be applied to the stick hlx over the native sequence
//validity right means that the stick hlx could or not be shifted to the right (toward the end of the sequence)...
//validity left means that the stick hlx could be or not be shifted to the left (toward the beginning of the sequence)
void setShiftValidity(Protein nativePDBFile, vector<EdgeStick> & sticks)
{
	if ((sticks.size()) && (sticks[0].hlxAssigned != -1))
	{
		for(int i=0;i<sticks.size();i++)
		{
			sticks[i].validityRight = 0;    //reset the value
			sticks[i].validityLeft  = 0;	 //reset the value
			
			int indx = sticks[i].hlxAssigned;
			int nativePDBNumOfAA = nativePDBFile.numOfAA();


			int numOfAADiff = nativePDBFile.hlces[indx].length - sticks[i].length;

			numOfAADiff /= 2;			//the shift amount left or right to centerally assign the hlx


			if (sticks[i].direction == 1)
			{
				sticks[i].startAAIndx = nativePDBFile.hlces[indx].startIndx + numOfAADiff;
		
				if (sticks[i].startAAIndx < 0)		//make it start from the beginning
					sticks[i].startAAIndx = 0;
		
				if (sticks[i].startAAIndx + sticks[i].length > nativePDBNumOfAA)						//make it start accordingly
					sticks[i].startAAIndx -= (sticks[i].startAAIndx + sticks[i].length) -  nativePDBNumOfAA;

				sticks[i].validityLeft = sticks[i].startAAIndx;	//the number of AA left of the stick if assigned to this hlx
				sticks[i].validityRight = nativePDBNumOfAA - (sticks[i].startAAIndx + sticks[i].length);	//the number of AA right to the stick if assigned to this hlx	
			}
			else
			{
				sticks[i].startAAIndx = nativePDBFile.hlces[indx].endIndx - numOfAADiff;
				//cout<<"startAAindx = "<<sticks[i].startAAIndx<<" "<<nativePDBFile.AAs[indx].chr3<<" "<<nativePDBFile.AAs[indx].num<<endl;
		
				if (sticks[i].startAAIndx >= nativePDBNumOfAA)				//make it start in the end of the pdb
					sticks[i].startAAIndx = nativePDBNumOfAA - 1;

				if (sticks[i].startAAIndx - sticks[i].length < -1)		//make it start accordingly
					sticks[i].startAAIndx += sticks[i].length - sticks[i].startAAIndx - 1;

				sticks[i].validityLeft = sticks[i].startAAIndx - sticks[i].length + 1;
				sticks[i].validityRight = nativePDBNumOfAA - sticks[i].startAAIndx - 1;
			}


		//	cout<<sticks[i].stickNum<<" with hlx"<<sticks[i].hlxAssigned+1<<"  Direction "<<sticks[i].direction		\
				<<"  startAAindx "<<sticks[i].startAAIndx<<"( "<<nativePDBFile.AAs[sticks[i].startAAIndx].chr3				\
				<<" "<<nativePDBFile.AAs[sticks[i].startAAIndx].num<<")  Validity  Right = "<<sticks[i].validityRight<<"  Left = "<<sticks[i].validityLeft <<endl;
		}
	}
	else
	{
		cout<<"========================= in skeleton::setReverseValidity (Protein) ==========="<<endl;
		cout<<"The 'sticks' data structure is empty or the hlces are not assigned to sticks..."<<endl;
		cout<<"==============================================================================="<<endl;
		putchar(BEEP);
		if (SHOW_ERRORS)
		{
			cout<<"Press any key..."<<endl;
			getchar();
		}
	}


}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//this function supposes that the shift value sent is always valid and also supposes that the direction of the stick is set
// portion is the stick hlx u want to assign AA's to it
//it is used to shift and assign directly
void AssignSeqToSkeleton(Protein nativePDBFile, Protein & portion, int shift, int stickIndx, vector<EdgeStick> sticks)
{

	int walkDirection;		//the direction to walk into the sequence...depends on the direction
	if (sticks[stickIndx].direction == 1)
		walkDirection = 1;
	else
		walkDirection = -1;

	int startIndx = sticks[stickIndx].startAAIndx + shift;		//positive to the right and negative to the left

	for (int i=0; i<sticks[stickIndx].length; i++)
	{
		portion.renameAA(i, nativePDBFile.AAs[startIndx].chr3);
		startIndx += walkDirection;
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void rotateAroundAxis(Protein &portion, double angleDegree)
{
	if (angleDegree)		//if there is a need to rotate
	{
		int tmpNumOfAA = portion.numOfAA();
		if (tmpNumOfAA >= 2)
		{

			Coordinate startEndCoord, endEndCoord;
			int AAsCounter = 0;	//for atoms will be involved in the center axial
			//take at most 4 atoms from each end
			while ((AAsCounter < tmpNumOfAA - 1) && (AAsCounter < 4))
			{
				Coordinate tmpCoord;

				//get coordinate of Ca atoms from the start end
				tmpCoord = portion.getAtomCoordinate(AAsCounter," CA ");	
				startEndCoord.x += tmpCoord.x;
				startEndCoord.y += tmpCoord.y;
				startEndCoord.z += tmpCoord.z;

				//get coordinate of Ca atoms from the end end
				tmpCoord = portion.getAtomCoordinate(tmpNumOfAA - AAsCounter - 1," CA ");
				endEndCoord.x += tmpCoord.x;
				endEndCoord.y += tmpCoord.y;
				endEndCoord.z += tmpCoord.z;

				AAsCounter++;
			}

			//find the average (center)
			startEndCoord.x = startEndCoord.x / AAsCounter;
			startEndCoord.y = startEndCoord.y / AAsCounter;
			startEndCoord.z = startEndCoord.z / AAsCounter;

			endEndCoord.x = endEndCoord.x / AAsCounter;
			endEndCoord.y = endEndCoord.y / AAsCounter;
			endEndCoord.z = endEndCoord.z / AAsCounter;

			//rotate around the axis
			portion.rotate(0,tmpNumOfAA-1,0,startEndCoord,endEndCoord, angleDegree);
		}
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////
//given one permutation stored in sticks, then generate all possible topologies and store it in allTopologies data structure
void generateAllTopologies(vector<EdgeStick> sticks, vector<sticksRec> & allTopologies)
{

	int nSticks = sticks.size(),
		i,
		j,
		k;
	vector<int> stickDirections,oneVect;

	sticksRec tmpStickRec;

	int  directionCarry=0;

	//initiate directions
	for (i=0; i<nSticks;i++)
	{
		stickDirections.push_back(0);
		oneVect.push_back(0);
	}

	oneVect[0] = 1;

	/*****
		generate all topologies starts from 001 (in case u have three sticks) and goes to 111 and finally generate 000 (0 means reverese and 1 means forward)
	*****/
	for (i=0; i < pow(2, nSticks); i++)
	{
		//get the next topology by changing the direction
		for (j=0; j<nSticks; j++)
		{
			stickDirections[j] = stickDirections[j] + oneVect[j] + directionCarry;
			if (stickDirections[j] >= 2)
			{
				directionCarry = 1;
				stickDirections[j] = 0;
			}
			else
				directionCarry = 0;
		}

		for (k=0; k< nSticks; k++)
		{
			sticks[k].direction = stickDirections[k];
		}
		tmpStickRec.sticks = sticks;

		allTopologies.push_back (tmpStickRec);
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////
double getDistancebwSticks(vector<Protein> skeleton, int indx1,int direction1, int indx2, int direction2)
{
	Coordinate  firstEnd,
				secondStart;

	if (direction1 == 1)
		//I am not using getAtomIndx here b/s it's guranteed 100% that the atom coordinate will be there
		firstEnd = skeleton[indx1].getAtomCoordinate(skeleton[indx1].numOfAA()-1, " C ");		
	else
		//I am not using getAtomIndx here b/s it's guranteed 100% that the atom coordinate will be there
		firstEnd = skeleton[indx1].getAtomCoordinate(0, " C ");

	if (direction2 == 1)
		//I am not using getAtomIndx here b/s it's guranteed 100% that the atom coordinate will be there
		secondStart = skeleton[indx2].getAtomCoordinate(0, " N ");
	else
		//I am not using getAtomIndx here b/s it's guranteed 100% that the atom coordinate will be there
		secondStart = skeleton[indx2].getAtomCoordinate(skeleton[indx2].numOfAA()-1, " N ");

	return getDistance(firstEnd, secondStart);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////
// generate all valid topologies for a given assignment (order or permutation)
// it pushes valid topologies to stickVect
inline
void generateValidTopologies(Protein nativePDBFile, vector<Protein> initialSkeleton, vector<int> permutation, vector<EdgeStick> sticks) 
{
		bool validPermutation = true;
		double lengthVariationHelix = 0.5;			//the maximum of length variation could occure between stick and sequence segment
		int j, k;
		/*****
				check the first condition, the variation on length should not be less than a predefined threshold
		*****/
		for (j=0; j<sticks.size(); j++)
		{
			if (((double) nativePDBFile.hlces[permutation[j] - 1].length/sticks[j].length < lengthVariationHelix) ||
				((double) sticks[j].length /nativePDBFile.hlces[permutation[j]-1].length < lengthVariationHelix ))
				validPermutation = false;			//not valid permutation
		}
		
		/*****
			if the permutation passes the first condition...then  check the other one
		*****/
		if (validPermutation)
		{
			/*****
				Assign hlces to sticks and according to above permutation and then check all possible topologies by reversing dirrections
			*****/
			sticksRec tmpStickRec;

			//Assign hlces to sticks
			for (j=0; j<sticks.size(); j++)
			{
				sticks[j].hlxAssigned = permutation [j] - 1;
			//	cout<<endl<<"Stick "<<sticks[j].stickNum<<" with hlx "<<sticks[j].hlxAssigned +1<<"   stick length= "<<sticks[j].length \
					<<" hlx actual length= "<<nativePDBFile.hlces[sticks[j].hlxAssigned].length<<" startAAindx= "<<sticks[j].startAAIndx<<endl;
			}
		//	cout<<"=========================================================================="<<endl;
			vector<sticksRec> allTopologies;

			/*****
					generate all topologies for this particular permutation
			*****/
			generateAllTopologies(sticks, allTopologies);
			/*****
					check the second condition for each topology
			*****/
		
			for (j=0; j<allTopologies.size(); j++)
			{
				//set the shift validity right and left... how many AA could the sequence be shifted right and left
				setShiftValidity(nativePDBFile, allTopologies[j].sticks);
		
				bool validTopology = true;

				vector<EdgeStick> tmpSticks = allTopologies[j].sticks;	//data structure to store temporariliy the current topology and then find the loop length b/w hlces
				vector<Protein> tmpInitialSkeleton = initialSkeleton;
				/*****
						sort according to the actual hlces, then find the loop length b/w any two adjacent hlces
				*****/
				int kk;
				for (kk= 0; kk<sticks.size(); kk++)
				{
					for (int m=0; m<sticks.size()-1; m++)
					{
						if (tmpSticks[m+1].hlxAssigned < tmpSticks[m].hlxAssigned)
						{
							EdgeStick tmpEdge = tmpSticks[m];
							tmpSticks[m] = tmpSticks[m+1];	
							tmpSticks[m+1] = tmpEdge;

							//sort the skeletons accordingly
							Protein tmpSkeleton = tmpInitialSkeleton[m];
							tmpInitialSkeleton[m] = tmpInitialSkeleton[m+1];
							tmpInitialSkeleton[m+1] = tmpSkeleton;
						}
					}
				}
				/*****
						I think we should flip tmpInitialSkeleton for those hlces assigned reversely
				*****/

				k=0;
				while ((validTopology) && (k<sticks.size()-1))
				{
					/*****
							find the virtual length of the loop b/w the two hlces 
							(the loop for right now will be considered the maximum loop could be b/w the two hlces, in other words, the maximum
							shift allowed should be considered in this step)
					*****/
					
					int loopLength = (tmpSticks[k].validityLeft  < MAX_SHIFT_ALLOWED) ? (tmpSticks[k].validityLeft) : (MAX_SHIFT_ALLOWED);		//get the maximum from the left
					loopLength += (tmpSticks[k+1].validityRight < MAX_SHIFT_ALLOWED) ? (tmpSticks[k+1].validityRight ) : (MAX_SHIFT_ALLOWED);	//get the maximum from the right

					if (tmpSticks[k].direction == 1)
					{
						if (tmpSticks[k+1].direction == 1)
							loopLength += tmpSticks[k+1].startAAIndx - tmpSticks[k].startAAIndx - tmpSticks[k].length;
						else
							loopLength += tmpSticks[k+1].startAAIndx - tmpSticks[k+1].length - tmpSticks[k].startAAIndx - tmpSticks [k].length + 1;
					}
					else
					{
						if (tmpSticks[k+1].direction == 1)
							loopLength += tmpSticks[k+1].startAAIndx - tmpSticks[k].startAAIndx - 1;
						else
							loopLength += tmpSticks[k+1].startAAIndx - tmpSticks[k+1].length - tmpSticks[k].startAAIndx;
					}

					/*****
							find the distance between the sticks. (the distance between last C terminus and the first N terminus)
					*****/
					double distance = getDistancebwSticks(tmpInitialSkeleton, k, tmpSticks[k].direction, k+1, tmpSticks[k+1].direction);

				//	cout<<"loop length = "<<loopLength<<"  distance = "<<distance<<endl;
					if ((loopLength * 3.8) < distance)
					{
						validTopology = false;
				//		cout<<"  NOT VALID..."<<endl;
					}
					k++;					
				}
			//	cout<<"==========="<<endl;
				/*****
						if the topology passes the second geometry screening test, then it is a valid topology
				*****/
				if (validTopology)
					sticksVect.push_back(allTopologies[j]);
			}
		}
}
//////////////////////////////////////////////////////////////////////////////////////////////////////
void genrateAllValidTopologies(Protein nativePDBFile, vector<Protein> initialSkeleton, vector<OnePermutation> allPermutations, vector<EdgeStick> sticks)
{
	//long int indxJump = 1;//factorial(nativePDBFile.hlces.size() - sticks.size());  was in the previous version
	long int i = allPermutations.size() - 1;
	/*****
		Determine Valid Topologies
			geometry screening used to determine the valid topologies is passing two conditions
			1. The difference in length should be less than 7
			2. The length of the loop in between should not exceed the length of extended loop (3.8 * num. of AA's in the loop)
	*****/
	while (i >= 0)
	{
		generateValidTopologies(nativePDBFile, initialSkeleton, allPermutations[i].permutation, sticks);
		allPermutations.pop_back ();			//delete last permutation...we are done with it...to save memory for big proteins
		i--;
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// End Of IMPLEMENTATION ////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif
