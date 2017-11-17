#ifndef SKELETON_OVERALL_H
#define SKELETON_OVERALL_H

#include <sstream>
#include <string.h>
#include <stdio.h>

#include "geometry.h"
#include "protein.h"
#include "constants.h"
#include "utilityfunctions.h"
#include "FBCCD.h"



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// List Of Data Structures //////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct EdgeStick
{
	Coordinate start;			//start point (Coordinate) from the edge file
	Coordinate end;				//end point (Coordinate) from edge file
	int stickNum;				//stick number as in edge file
	int ssAssigned;				//the indx of the secondary structure assigned to this stick (indx as in SecondaryStruct Data Structure)
	char ssType;				//the type of the secondary structure H Helix or S strand
	int direction;				//the direction of the sequence (how is the sequence assigned to the stick) 1 forward  0 reverese
	int nAA;					//the assumed number of AA's
	int startAAIndx;			//the indx of first AA (AAs indx as in PDB)
	int shift;					//the amound of shift from the original SS (shift in the sequence not on the structure)
	int translation;			//the amound of translation	(move right or left) on the structure
	int validityRight;			//the maximum size of shift could be applied to the sequence to the right...
	int validityLeft;			//the maximum size of shift could be applied to the sequence to left...

	double rotation;			//the amount of rotation from initial skeleton


	//initializer
	EdgeStick(): stickNum(-1), ssAssigned(-1), ssType('H'), direction(1), nAA(0), startAAIndx(-1), \
		shift(0), translation(0), rotation(0), validityRight(0), validityLeft(0) {}
};

//structure used to store all kind of secondary structure and their corresponding assignment in the PDB file
struct SecondaryStruct
{
	int startIndx;		//corresponding indx in the PDB of the first AA in the secondary structure
	int startAAnum;		//seq number of the first AA
	int endIndx;		//corresponding indx in the PDB of the last AA in the secondary structure
	int endAAnum;		//seq number of the last AA
	int nAA;			//the length of secondary structure (how many actual AA is in there...missing AA should be counted)
	char type;			//H Helix or S strand
	short order;		//the order of this SS regards to its type...i.e first , second hlx..ot third strand

	//initializer
	SecondaryStruct(): startIndx(-1), startAAnum(0), endIndx(-1), endAAnum(0), nAA(0), type('H') {}
};


#define MAX_SHIFT_ALLOWED	0
#define MAX_TRANSLATION_ALLOWED 3
#define MAX_ROTATION_ALLOWED 0

//used to save the sticks given by edge file...Global Variable
//static vector<EdgeStick> sticks;
//data structures that stores all sticks vectors (all topologies)
struct sticksRec
{
	vector<EdgeStick> sticks;
};

//vector<sticksRec> sticksVect;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// End Of Constants and Data Structures /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////// list of functions /////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void printSkeleton(vector<Protein> hlces, vector<EdgeStick> sticks, string path, string TAG = "");						//print skeleton as pdb files
void fineStick(vector<Coordinate> &SS, float dist);
// To build the structure of template SS (sticks)
void edge2points(string edgeFile, vector<vector<Coordinate> > & edges);			//read given Edge file (Helix Tracer output) and extract the points and store them in a vector
Protein points2pdb(vector<Coordinate> lPoints, string ="GLY", string = " CA ");							//Convert given vector of points to a viewable pdb file by chimera
Protein points2stick(vector<Coordinate> lPoints);
void setAxis(Protein SS, short axisAA, vector<Coordinate> &axis);			//set the axis of a stright SS by a set of points every axisDist AAs
vector<Coordinate> extractCurves(vector<Coordinate> edge, float rise, float cutOffDist, bool direction);	//extract curved line from given some of points
Protein buildStructure(vector<Coordinate> points, char sType, bool direction);		//given a line, represented by some of points, build SS stick
void combineSticks(vector<EdgeStick> hlces, vector<EdgeStick> strands, vector<EdgeStick> & allSticks, bool prnt=false);
void printSticksInfo(vector<EdgeStick> allSticks, int NhlcesSticks, int NstrandsSticks);
vector<EdgeStick> getHlcesSticks(string mrcFile);														//build the sticks vector that contains the first and last points of each stick
void overlapSticks(vector<Protein> & sticksSS, vector<EdgeStick> sticks);			//translate the helices to overlap the sticks given by Edge file those stored in sticks data structure
void flipHelix(Protein & stickHlx, EdgeStick stick);									//flip a helix
vector<Protein> getInitialSkeleton(Protein nativePDBFile, string mrcFile, vector<EdgeStick> & sticks, string outDir);	// builds initial skeleton for the given edge file and native pdb file
///////
void buildSeqSS(Protein nativePDB, vector<SecondaryStruct> & seqSS, bool prnt=false);				//build the vector contain sequence secondary structures
void printSSlist(Protein nativePDB, vector<SecondaryStruct> SSlist);						//print list of Secondary structre on the original order
void setShiftValidity(Protein nativePDBFile, vector<EdgeStick> & sticks);								//set the amound of shift to left and right of the sticks hlces allowed
void AssignSeqToSkeleton(Protein nativePDBFile, Protein & portion, int shift, int stickIndx, vector<EdgeStick> sticks);	//assign AA seq from native pfb file to the sticks hlces
void rotateAroundAxis(Protein &portion, double angleDegree);											//rotate a portion (hlx mainly) around its axis
void generateAllTopologies(vector<EdgeStick> sticks, vector<sticksRec> & allTopologies);	//given one permutation stored in sticks, then generate all possible topologies and store it in allTopologies data structure
double getDistancebwSticks(vector<Protein> skeleton, int indx1,int direction1, int indx2, int direction2);	//get the distance between two skeleton hlces
inline
void generateValidTopologies(Protein nativePDBFile, vector<Protein> initialSkeleton, vector<SecondaryStruct> seqSS, vector<int> permutation, vector<EdgeStick> sticks);

//void genrateAllValidTopologies(Protein nativePDBFile, vector<Protein> initialSkeleton, vector<OnePermutation> hlcesPermutations, vector<OnePermutation> strandsPermutation, vector<EdgeStick> sticks); //given the pdb file and initial stick information and all permutations, generate all valid topologies into sticksVect


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
		outFileName += "stick_"  ;
		intStr<<sticks[i].stickNum;
		outFileName  += intStr.str();
		if (sticks[i].ssType == 'H')
			outFileName  += "_hlx_";
		else
			outFileName += "_Strand_";
		intStr2<<sticks[i].ssAssigned + 1;
		outFileName  += intStr2.str();
		outFileName += "_";
		outFileName += TAG;
		outFileName += ".pdb";

		hlces[i].writePDB (outFileName,1,hlces[i].numOfAA());
		//cout<<"Stick number "<<sticks[i].stickNum<<"  assigned to hlx number  "<<endl;
	}

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void fineStick(vector<Coordinate> &SS, float dist){

	int j;
	bool cont = true;
	Coordinate center;

	while (cont){
		cont= false;
		j=0;
		while (j<SS.size ()-1){

			if (getDistance(SS[j], SS[j+1]) > dist){
				cont = true;
				center.x = (SS[j].x + SS[j+1].x)/2;
				center.y = (SS[j].y + SS[j+1].y)/2;
				center.z = (SS[j].z + SS[j+1].z)/2;


				SS.insert (SS.begin () + j+1, center);
			}
			j++;

		}
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void reWritePntsOnStick(vector<Coordinate> & stick, float dist){
	int i=0, curP;
	vector<Coordinate> tmpStick;
	Coordinate tmpP;

	//push the first point
	tmpStick.push_back (stick[0]);
	curP=1;

	while (i<tmpStick.size () && curP<stick.size ()){

		if (getDistance(tmpStick[i], stick[curP])>dist){
			tmpP = pointOnLine(tmpStick[i], stick[curP], dist);
			tmpStick.push_back (tmpP);
			i++;
		}
		else
			curP++;
	}

	stick = tmpStick;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// read given Edge file (Helix Tracer output) and extract the points and store them in a vector ///////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void edge2points(string edgeFile, vector<vector<Coordinate> > & edges)
{
	ifstream inFile;
	string line;
	EdgeStick tmpStick;
	Coordinate point;
	vector<vector<Coordinate> > tmpEdges;


	inFile.open(edgeFile.c_str());
	if (!inFile)
	{
		cout<<"=================================== in skeleton::edge2Points(.......) ============="<<endl;
		cerr << "Unable to open " << edgeFile << endl;
		cout<<"==================================================================================="<<endl;
		exit(1);
	}
	else
	{
		int prevEdgeNum = -1, curEdgeNum;
		char *token;
		char delim[]="	, ";
		while (!inFile.eof())
		{
			getline(inFile, line);
			if ((line.length ()) && (line.substr(0,1)==" "))
			{

				char *ss= (char*)line.c_str ();
				//ss>>tmpStick.start.x>>tmpStick.start.y>>tmpStick.start .z>>stickNum;

				//tokanize the line
				token=strtok(ss, delim);
				if (token != NULL)
					point.x = atof(token);

				token = strtok(NULL, delim);
				if (token!= NULL)
					point.y = atof(token);

				token = strtok(NULL, delim);
				if (token!=NULL)
					point.z = atof(token);

				token = strtok(NULL, delim);
				if(token!= NULL)
					curEdgeNum = atoi(token);


				//int curEdgeNum = atoi(line.substr(42,3).c_str());
				if (prevEdgeNum != curEdgeNum)				//new edge
				{
					vector<Coordinate> lPoints;		//one Edge

					prevEdgeNum = curEdgeNum;
					//point.x  = atof(line.substr(0,13).c_str());
					//point.y  = atof(line.substr(14,13).c_str());
					//point.z  = atof(line.substr(28,13).c_str());

					lPoints.push_back (point);

					tmpEdges.push_back (lPoints);
				}
				else
				{
					//working on the same Edge line
					//point.x  = atof(line.substr(0,13).c_str());
					//point.y  = atof(line.substr(14,13).c_str());
					//point.z  = atof(line.substr(28,13).c_str());

					tmpEdges[tmpEdges.size ()-1].push_back (point);
				}
			}
		}

	}
	inFile.close ();


	// remove duplicates in the original edge file...alot of duplicated points
	int i, j, k;
	for (i=0; i<tmpEdges.size (); i++)
	{
		vector<Coordinate> lPoints;
		edges.push_back (lPoints);
		for (j= 0; j<tmpEdges[i].size(); j++)
		{
			bool duplicated = false;
			for (k=0; k<edges[i].size(); k++)
				if ((tmpEdges[i][j].x == edges[i][k].x) &&
					(tmpEdges[i][j].y == edges[i][k].y) &&
					(tmpEdges[i][j].z == edges[i][k].z))
				{
					duplicated = true;
					break;
				}
			if (!duplicated)
				edges[i].push_back (tmpEdges[i][j]);

		}
	}

	/*
	// print out points
	int i, j;
	for (i=0; i<edges.size (); i++)
	{
		cout<<"Edge# "<<i+1<<" #Points= "<<edges[i].size()<<endl;
		for (j=0; j<edges[i].size(); j++)
		{
			cout<<"   "<<edges[i][j].x<<" "<<edges[i][j].y<<" "<<edges[i][j].z<<endl;
		}
		cout<<endl;
	}
	*/
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Convert given vector of points to a viewable pdb file by chimera ///////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Protein points2pdb(vector<Coordinate> lPoints, string AA, string atom)
{
	Protein tmp;

	int i;
		//write sticks as pdb files...with CA atoms to be viewed using chimera
	for (i=0; i<lPoints.size(); i++)
	{
		AminoAcid tmpAA;
		Atom tmpAtom;
		tmpAtom.coord = lPoints[i];
		tmpAtom.name = atom;
		tmpAA.atoms.push_back(tmpAtom);
		tmpAA.num = i+1;

		//cout<<"working on AA# "<<i+1<<" : "<<lPoints[i].x<<" "<<lPoints[i].y<<" "<<lPoints[i].z<<endl;
		tmpAA.chr3 = AA;
		tmpAA.num = i+1;
		tmpAA.chain = "A";
		tmp.AAs.push_back(tmpAA);
	}

	return tmp;

}
////////////////////////////////////
Protein points2stick(vector<Coordinate> lPoints)
{
	Protein tmp;

	int i;
		//write sticks as pdb files...with O atoms to be viewed using chimera
	for (i=0; i<lPoints.size(); i++)
	{
		AminoAcid tmpAA;
		Atom tmpAtom;
		tmpAtom.coord = lPoints[i];
		tmpAtom.name = " CA ";
		tmpAA.atoms.push_back(tmpAtom);

		//cout<<"working on AA# "<<i+1<<" : "<<lPoints[i].x<<" "<<lPoints[i].y<<" "<<lPoints[i].z<<endl;
		tmpAA.chr3 = "GLY";
		tmpAA.num = i+1;
		tmpAA.chain = "A";
		tmp.AAs.push_back(tmpAA);
	}

	return tmp;

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// given a list of points, this function's to set a point every specific number of AAs represented
//// by cutOffDist, for example for helix the cutOffDist would be helix rise * nAA.
//// The direction represents the direction of ur generated curve..follow the same ....
//// direction of given edge 1 (true), or reverse direction 0 (false) (the start of curved line returned is the end of original line sent)
//// rise is the rise used to estimate the last portion of the curved line
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<Coordinate> extractCurves(vector<Coordinate> edge, float rise, float cutOffDist, bool direction)
{
	vector<Coordinate> nEdge;		//new edge points
	if (edge.size())
	{
		short increment,
				sIndx,		//start indx
				eIndx;		//end indx
		if (direction)
		{
			increment	=	1;
			sIndx		=	0;
			eIndx		=	edge.size();
		}
		else
		{
			increment	=	-1;
			sIndx		=	edge.size()-1;
			eIndx		=	-1;
		}

		Coordinate curP = edge[sIndx];	//current point we are working on...first point in curved line

		nEdge.push_back(curP);

		int lastPointIndx = -1;
		int i= sIndx + increment;
		while (i != eIndx)
		{
			if (getDistance(curP, edge[i]) > cutOffDist)
			{
				nEdge.push_back (pointOnLine(curP, edge[i], cutOffDist));
				curP = nEdge[nEdge.size()-1];
				lastPointIndx = i-increment;
			}
			i += increment;
		}

		// for last point...expect the number of remianing AAs by dividing the remining distance by rise
		// then draw the point according this calculations
		if (lastPointIndx != eIndx)
		{
			//calculate the number of remaining number of AAs (expected by rise)
			int nAA = int (getDistance(curP, edge[eIndx-increment])/rise + 0.5);
			if (nAA>0)
				nEdge.push_back(pointOnLine(curP, edge[eIndx-increment], nAA*rise));
		}

		//for (i=0; i<nEdge.size()-1; i++)
		//	cout<<getDistance(nEdge[i], nEdge[i+1])<<" eNAA= "<<getDistance(nEdge[i], nEdge[i+1])/ALPHA_RISE<<endl;

	}
	else
	{
		cout<<"================= in extractCurves (vector<Coordinate> edge) =================="<<endl;
		cout<<"The given edge vector is empty..."<<endl;
		cout<<"==============================================================================="<<endl;
		exit(1);
	}

	return nEdge;

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Given a striaght SS,
//// set the center (Axis) of a SS every specific number of AAs..
//// axisAA represents the number of AAs between any 2 points on the axis
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void setAxis(Protein SS, short axisAA, vector<Coordinate> &axis)
{
	short nAA = SS.numOfAA();	//number of AAs in the structure
	//calculate the axis for the straight structure
	if (nAA > -1)
	{
		if (nAA>3){  //4 so you have 2 triangles which will represent the initial axis of the structure
			int i;
			//clear previous axis
			axis.clear();
			//calculate the center of backbone first and last AAs
			SS.setBBCenter(0, 'H');			//heavy atoms only
			SS.setBBCenter(1, 'H');
			SS.setBBCenter(2, 'H');
			SS.setBBCenter(nAA-3,'H');
			SS.setBBCenter(nAA-2, 'H');
			SS.setBBCenter(nAA-1, 'H');

			// get the center of first three AAs
			Coordinate pOne = triangleCenter(SS.AAs[0].coord, SS.AAs[1].coord, SS.AAs[2].coord),
						pLast = triangleCenter(SS.AAs[nAA-3].coord, SS.AAs[nAA-2].coord, SS.AAs[nAA-1].coord);


			Coordinate iPoint;
			//now for every axisAA, set a point on the axis
			for (i=0; i<nAA; i+= axisAA)
			{
				//the corresponding point of AA (CA) on the axis...the imiginary point of AA on the axis
				iPoint = pointLineIntersection(SS.getAtomCoordinate(i," CA "),pOne,pLast);
				axis.push_back (iPoint);
			}
			i -= axisAA;
			if (i <nAA-1)		//last AA
			{
				iPoint = pointLineIntersection(SS.getAtomCoordinate(nAA-1," CA "),pOne,pLast);
				axis.push_back (iPoint);
			}
		}
		else {
			if (nAA == 2){
				axis.push_back (SS.getAtomCoordinate(0, " N  "));
				axis.push_back (SS.getAtomCoordinate(1, " C  "));
			}
			else
				axis.push_back (SS.getAtomCoordinate(0, " CA "));
		}
	}
	else
	{
		cout<<"================= in setAxis (Protein, axisAA) ================================"<<endl;
		cout<<"The number of AA in the portion sent should be greater than 0..."<<endl;
		cout<<"==============================================================================="<<endl;
		exit(1);
	}

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// points vector represents the list of points extracted from edge file (from HT)..
// direction : true forward, same direction with points list, or false reverse.
// sType is the type of structure you are building
Protein buildStructure(vector<Coordinate> points, char sType, bool direction)
{
	EdgeStick tmpStick;
	vector<Coordinate> nEdge;		//new edge line
	float rise;						// rise is the rise of the structure you are building.
	short nRise=4;					// nRise : helps with where to put the list of points on 3D...
									// for example if it is 2 it means draw a point every 2 rises (3 AAs)....
									// points here represent the imiginary of CA atom of the AA on the central axis of the SS
	float phi, psi;
	int i;

	if (sType == 'H')
	{
		rise	= ALPHA_RISE;
		phi		= ALPHA_HLX_PHI_RIGHT;
		psi		= ALPHA_HLX_PSI_RIGHT;
	}
	else if(sType == 'S')
	{
		rise	= BETA_RISE;
		phi		= BETA_SHEET_PHI;
		psi		= BETA_SHEET_PSI;
	}
	else
	{
		rise	= LOOP_RISE;				//loop rise
		phi		= getRandom(-175, -40);		//random phi
		psi		= getRandom(-60, 175);		//random psi
	}

	// extract a curved line from the given list of points....
	// draw a point every rise * nAAcutOff Angstrom
	nEdge = extractCurves(points, rise, rise * nRise, direction);

//	Protein curvd;
//	curvd = points2pdb(nEdge,"HOH","  O ");
//	curvd.writePDB("curvedhlx.pdb", 1, curvd.numOfAA());

	//calculate the expected length of the stick SS
	double tDist=0;
	for (i=0; i<nEdge.size()-1; i++)
		tDist += getDistance(nEdge[i], nEdge[i+1]);
	int nAA = int ((tDist/rise) + 0.5);		//nAA expected on the structure
	nAA +=1;		//for example : 2 rise on axis = 3 AAs


	// build the structure
	Protein tmpProt;
	AminoAcid tmpAA;

	tmpAA.chr3 = "GLY";
	tmpAA.num = 1;
	tmpAA.angles.phi = phi;
	tmpAA.angles.psi = psi;

	for(i=0; i<nAA; i++)
		tmpProt.AAs.push_back(tmpAA);
	//generate structure from a list of torsion angles
	tmpProt.torsion2coord();

	//set SS type flag
	for (i=0; i<nAA; i++)
		tmpProt.AAs[i].SStype = sType;


	// Now overlap the beginning of the generated structure with the curved line
	//set moving points which represent the axis of the generated structure
	vector<Coordinate> mPoints;
	setAxis(tmpProt, nRise, mPoints);
	//set moving points
	Coordinate p;
	p.x	= nEdge[0].x - mPoints[0].x;
	p.y	= nEdge[0].y - mPoints[0].y;
	p.z	= nEdge[0].z - mPoints[0].z;

	tmpProt.translateBy(p);		//move start to start
	Vectors v1(nEdge[nEdge.size()-1], nEdge[0]),
			v2(tmpProt.AAs[tmpProt.numOfAA()-2].atoms[0].coord, nEdge[0]),
			normal;

	//get angle b/w two vectors...to make them on the same direction (roughly)
	double angle = v1.getAngleDegree(v2);

	//find the normal b/w two vectors
	normal = v1.cross(v2);

	normal += nEdge[0];

	tmpProt.rotate(0, tmpProt.numOfAA()-1, 0, nEdge[0], normal.getCoordinates(), angle);

	//tmpProt.writePDB("builtHlx.pdb",1,tmpProt.numOfAA());
	//reset moving points
	setAxis(tmpProt, nRise, mPoints);		//again after translation and rotation


	if (mPoints.size() == nEdge.size())
		//overlap the moving points (mPoints) with target points nEdge by modifying tmpProt structure
		overlapFBCCD(tmpProt, nRise, mPoints, nEdge, 0.1, 100);
	else
	{
		cout<<"=================================== in buildStructure(..........) ================="<<endl;
		cerr << "Number of moving points ( "<<mPoints.size()<<" ) and target points ( "<<nEdge.size()<<" ) are not equal."<< endl;
		cout<<"==================================================================================="<<endl;
		exit(1);
	}

	return tmpProt;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//this function always builds the structure in the same direction of the given list of points in tracePnts
//if the direction of built structure wanted to be reverse....you should reverse the points in the tracePnts before calling this function
//used mainly to build loop and strand...not helix
Protein buildStructureOnTrace(Protein &portion1, Protein &portion2, int p2Indx, vector<Coordinate> &tracePnts,  vector<vector<Coordinate> > traceToAvoid, vector<char> &aaTypes){

	vector<Coordinate> nEdge,		//edge line represents some points on the trace that we want to overlap the built structure with
						mEdge;		//some points on the axis of the built structure to be overlapped with nEdge

	float distOnTrace = 0.5;		//the distance b/w each pair of points on the trace.
	short loopRise=15;				//the desired rise of the trace of the loop model (where every rise we will have a point)

	int i;
    short nAAportion1 = portion1.numOfAA();				//original number of AA in portion1
    short nLoopAA = aaTypes.size();

	Protein sModel;			//structure model for the given trace

    short loopType=0;           //how many amino acid of type loop

    //build a loop at the end of portion 1 but without colliding with the rest of the structure
    //add one AA at the beginning and one AA at the end
    aaTypes.insert(aaTypes.begin(), 'H');
    aaTypes.push_back('H');
    int nAttempts = 0;
    do{
		//cout<<"trying to generate a random loop.."<<endl;
        //clear the loop if previously connected
        while (nAAportion1 != portion1.numOfAA())
            portion1.deleteAA(portion1.numOfAA()-1);
        portion1.AAsCollide.clear();
        portion1.atomsCollide.clear();

        sModel = generateRandomStructure(aaTypes, aaTypes.size());

        portion1.connect(sModel, 'E');
        nAttempts++;
        //getchar();

	}while (portion1.doesCollide(0.5) && nAttempts<20);

    //update each AA type
    for (i=0; i<nLoopAA; i++){
        if (aaTypes[i+1] == 'L')
            loopType++;
    }

    /*
	portion1.writePDB("initialModel.pdb", 1, portion1.numOfAA());

    Protein curvd;
	//write the moving structures
	curvd = points2pdb(tracePnts,"GLY"," CA ");
	curvd.writePDB("initialTrace.pdb", 1, curvd.numOfAA());
	*/

	//rewrite points on the trace so the distance b/w each two points is distOnTrace
	reWritePntsOnStick(tracePnts, distOnTrace);

	//get the axis of the sModel
	//save the Ca trace of sModel into new vector
	vector<Coordinate> mPoints;
	for (i=0; i<=nLoopAA; i++){
		mPoints.push_back (portion1.getAtomCoordinate(nAAportion1+i, " CA "));
	}

	reWritePntsOnStick(mPoints, distOnTrace);
	reWritePntsOnStick(mPoints, distOnTrace *((float) mPoints.size () / (float) tracePnts.size ()));

	//write the moving structures
	/*
	curvd.initialize();
	curvd = points2pdb(mPoints,"GLY"," CA ");
	curvd.writePDB("mPoints.pdb", 1, curvd.numOfAA());

    cout<<"Number of points on the trace: "<<tracePnts.size()<<"  Number of points on initial Model: "<<mPoints.size()<<endl;
    */

	//set target and moving points
	i = loopRise*2;			//means have a point every loopRise/2 Angstrom....assuming the distance b/w any consecutive points on the trace distOnTrace is always = 0.5
	nEdge.push_back (tracePnts[0]);
	mEdge.push_back (tracePnts[0]);
	while (i<tracePnts.size () && i<mPoints.size ()){
		nEdge.push_back(tracePnts[i]);
		mEdge.push_back(mPoints[i]);

		i += loopRise*2;
	}
	//push last point if is not already pushed
	if (i-loopRise*2 < tracePnts.size()-1 || i-loopRise*2 < mPoints.size()-1){
        nEdge.push_back (tracePnts[tracePnts.size ()-1]);
        mEdge.push_back (mPoints[mPoints.size ()-1]);
    }


	//save the final trace back to tracePnts
	tracePnts = nEdge;

    //write the two traces
    /*
	curvd.initialize();
	curvd = points2pdb(nEdge, "GLY", " CA ");
	curvd.writePDB("tracePnts.pdb", 1, curvd.numOfAA());

	cout<<"#of AA in portion1= "<<nAAportion1<<" loop= "<<nLoopAA<<" portion2= "<<portion2.numOfAA()<<" total after connect= "<<portion1.numOfAA()<<endl;
	cout<<"Number of edge pnts= "<<nEdge.size ()<<"  Number of mEdge= "<<mEdge.size()<<endl;
	*/

	//find the correspondant AA for each point in mEdge
	vector<short> correspondantAA(mEdge.size (), -1);
	float dist, minDist;
	correspondantAA[0] = nAAportion1;	//always the first point is correspondant to the first AA in the loop - 3

	//cout<<"Pnt "<<1<<" corr. to AA num "<<portion1.AAs[correspondantAA[0]].num<<endl;
	for (i=1; i<mEdge.size (); i++){
		minDist = 9999.9;
		short aaIndx = -1;
		for (int j=correspondantAA[i-1]+1; j<portion1.numOfAA(); j++){
			//get the closest distance to the edge point
			dist = getDistance(mEdge[i], portion1.getAtomCoordinate(j, " N  "));
			if ( dist < minDist){
				aaIndx = j;
				minDist = dist;
			}
			dist = getDistance(mEdge[i], portion1.getAtomCoordinate(j, " CA "));
			if ( dist < minDist){
				aaIndx = j;
				minDist = dist;
			}
			dist = getDistance(mEdge[i], portion1.getAtomCoordinate(j, " C  "));
			if ( dist < minDist){
				aaIndx = j;
				minDist = dist;
			}
		}
		if (aaIndx != -1){
            correspondantAA[i] = aaIndx;
            //set the movable point to be at the middle b/w N and CA of that amino acid
            mEdge[i].x = (portion1.getAtomCoordinate(aaIndx, " N  ").x + portion1.getAtomCoordinate(aaIndx, " CA ").x)/2;
            mEdge[i].y = (portion1.getAtomCoordinate(aaIndx, " N  ").y + portion1.getAtomCoordinate(aaIndx, " CA ").y)/2;
            mEdge[i].z = (portion1.getAtomCoordinate(aaIndx, " N  ").z + portion1.getAtomCoordinate(aaIndx, " CA ").z)/2;

            //cout<<"Pnt "<<i+1<<" corr. to AA num "<<portion1.AAs[correspondantAA[i]].num<<" minDist= "<<minDist<<endl;
		}
		else{
            mEdge.pop_back();
            nEdge.pop_back();
            //cout<<" Pnt "<<i+1<<" is deleted becuase no corr. AA could be found."<<endl;
		}
	}
    /*
    curvd.initialize();
	curvd = points2pdb(mEdge,"GLY"," CA ");
	curvd.writePDB("mEdge.pdb", 1, curvd.numOfAA());
	*/


	//start loop Modeling on the top of the target trace nEdge
	if (mEdge.size() == nEdge.size()){
		//cout<<"Loop Modeling..."<<endl;

		//int lastIndxOfPortion1 = portion1.numOfAA()-1;		//the index of last AA (portion1 + loop)

		if (mEdge.size()>2){
			i=1;
			vector<char> segAATypes (aaTypes.begin (), aaTypes.begin () + correspondantAA[1]-nAAportion1+1) ;

			//delete the loop from portion1...we will be modeling the loop segment by segment
			while (nAAportion1 != portion1.numOfAA())
                portion1.deleteAA(nAAportion1);

			while(i<mEdge.size ()-1){

				//cout<<"startAA = "<<correspondantAA[i-1]<<" endAA= "<<correspondantAA[i]<<" portion1 size now = "<<portion1.numOfAA()<<endl;

				randomFBCCD1Point(portion1, nEdge[i], correspondantAA[i], segAATypes, traceToAvoid, 0.5);

				i++;

				//getchar();getchar();
				segAATypes.clear();
				for (int j=correspondantAA[i-1]-nAAportion1+2; j<=correspondantAA[i]-nAAportion1+1; j++){
					//cout<<"copying AA indx "<<j<<endl;
					segAATypes.push_back (aaTypes[j]);
				}

			}
            //remove last AA from the loop (the first AA of the second portion been added before)
            segAATypes.pop_back();
            //cout<<"portion1 + loop #AA= "<<portion1.numOfAA()+segAATypes.size()<<endl;
			//overlap last loop segmen and second portion together
			float rmsdLastPortion = randomFBCCD(portion1, portion2, p2Indx, segAATypes, traceToAvoid, 0.05);
			int kk=0;
			while (rmsdLastPortion > 1.7 && kk+2<correspondantAA.size()){
				//cout<<"failed to overlap...try from the previous segment..."<<endl;
				//delete the second portion
				while (correspondantAA[i-kk-2] != portion1.numOfAA())
					portion1.deleteAA(portion1.numOfAA()-1);
				segAATypes.clear ();
				for (int j=correspondantAA[i-kk-2]-nAAportion1+2; j<=correspondantAA[i]-nAAportion1+1; j++){
					//cout<<"copying A indx "<<j<<endl;
					segAATypes.push_back (aaTypes[j]);
				}
                //cout<<"Portion1 + loop #AA= "<<portion1.numOfAA()+segAATypes.size()<<endl;
				randomFBCCD(portion1, portion2, p2Indx, segAATypes, traceToAvoid, 0.05);
				kk++;
			}
		}
        else{

			//delete the loop from portion1...we will be modeling the loop segment by segment
			while (nAAportion1 != portion1.numOfAA())
                portion1.deleteAA(nAAportion1);

			aaTypes.erase (aaTypes.begin ());
			aaTypes.pop_back();
			if (nLoopAA<6 && loopType<3)
				for (i=0; i<aaTypes.size(); i++)
					aaTypes[i] = 'L';

            //cout<<"portion1 + loop #AA= "<<portion1.numOfAA()+aaTypes.size()<<endl;

			randomFBCCD(portion1, portion2, p2Indx, aaTypes, traceToAvoid, 0.05);
        }

		//add current trave to the trace to avoid
		//traceToAvoid.push_back(tracePnts);
	}
	else
	{
		cout<<"=================================== in buildStructure(..........) ================="<<endl;
		cerr << "Number of moving points ( "<<mEdge.size()<<" ) and target points ( "<<nEdge.size()<<" ) are not equal."<< endl;
		cout<<"==================================================================================="<<endl;
		exit(1);
	}

	return portion1;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//build and return the sticks vector that contains the first and last points of each stick
vector<EdgeStick> getHlcesSticks(string mrcFile)
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

		// fill the following information for each stick
		// stickNum
		// start coordinate
		// end coordinate
		// ssType

		vector<EdgeStick> sticks;
		tmpStick.stickNum = -1;
		float x,y,z;
		int stickNum;
		char delim[]="	, ";
		char *token;

		while (!inFile.eof()){
			getline(inFile,str);
			if (str.substr(0,1)==" "){
				char *ss= (char*)str.c_str ();
				//ss>>tmpStick.start.x>>tmpStick.start.y>>tmpStick.start .z>>stickNum;

				//In the first call to strtok, the first argument is the line to be tokenized
				token=strtok(ss, delim);
				if (token != NULL)
				tmpStick.start.x = atof(token);

				token = strtok(NULL, delim);
				if (token!= NULL)
				tmpStick.start.y = atof(token);

				token = strtok(NULL, delim);
				if (token!=NULL)
				tmpStick.start.z = atof(token);

				token = strtok(NULL, delim);
				if(token!= NULL)
				stickNum = atoi(token);

				//cout<<tmpStick.start.x<<" "<<tmpStick.start.y<<" "<<tmpStick.start.z<<"  "<<stickNum<<endl;



				//int stickNum = atoi(str.substr(42,3).c_str());
				if (stickNum != tmpStick.stickNum){
				//	tmpStick.stickNum = stickNum;
				//	tmpStick.start.x  = atof(str.substr(0,13).c_str());
				//	tmpStick.start.y  = atof(str.substr(14,13).c_str());
				//	tmpStick.start.z  = atof(str.substr(28,13).c_str());
				//	tmpStick.ssType = 'H';			//helix stick
					tmpStick.stickNum = stickNum;
					tmpStick.start.x = tmpStick.start.x;
					tmpStick.start.y = tmpStick.start.y;
					tmpStick.start.z = tmpStick.start.z;

					sticks.push_back(tmpStick);
				}
				else {
				//	sticks[sticks.size()-1].end.x  = atof(str.substr(0,13).c_str());
				//	sticks[sticks.size()-1].end.y  = atof(str.substr(14,13).c_str());
				//	sticks[sticks.size()-1].end.z  = atof(str.substr(28,13).c_str());
					sticks[sticks.size()-1].end.x  = tmpStick.start .x;
					sticks[sticks.size()-1].end.y  = tmpStick.start .y;
					sticks[sticks.size()-1].end.z  = tmpStick.start .z;
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
vector<EdgeStick> getStrandsSticks(Protein nativePDB, vector<EdgeStick> sticks)
{
	vector<EdgeStick> strandSticks;
	EdgeStick tmpStick;
	int i;

	// fill the following information for each stick
	// stickNum
	// start coordinate
	// end coordinate
	// ssType

	//get the maximum stick number for helces...so u can number strands accordingly
	//note that stickNum for hlces is not sequential.... it is not in ascending or descending order
	int maxSticknum = -1;
	for (i=0; i<sticks.size(); i++)
		if (sticks[i].stickNum > maxSticknum)
			maxSticknum = sticks[i].stickNum;

	tmpStick.ssType = 'S';

	//first strand number starts from maximum stick number for hlces plus 2
	tmpStick.stickNum = maxSticknum + 2;

	for (i=0; i<nativePDB.sheets .size (); i++)
	{
		tmpStick.stickNum += 1 ;		//set the number of the stick arbitrarly
		tmpStick.start = nativePDB.getAtomCoordinate (nativePDB.sheets [i].startIndx, " CA ");			//get the coordinate of the first CA atom in the strand
		tmpStick.end   = nativePDB.getAtomCoordinate (nativePDB.sheets [i].endIndx, " CA ");			//get the coordinate of the last CA in the strand

		strandSticks.push_back(tmpStick);
	}

	return strandSticks;

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void combineSticks(vector<EdgeStick> hlces, vector<EdgeStick> strands, vector<EdgeStick> & allSticks, bool prnt)
{
	int i, j;
	double length;
	//set lengths of sticks (The suggested number of AA in each stick)
	for (i=0; i<hlces.size (); i++)
	{
		length = getDistance(hlces[i].start, hlces[i].end);

		//cout<<"length = "<<length<<endl;
		//The suggested length of the hlx (The number of AA in the hlx)
		 hlces[i].nAA = int ((length/ALPHA_RISE) + 0.5 /* To round it */); // Here the approx. number of aa per sticks..u can also divide by 1.33

		allSticks.push_back (hlces[i]);
	}

	//set lengths for strands
	for(i=0; i<strands.size (); i++)
	{
		length = getDistance(strands[i].start, strands[i].end);
		strands[i].nAA = int ((length/BETA_RISE) + 0.5) + 1;	//add another AA
		allSticks.push_back (strands[i]);
	}

	/*	we can ignore this step
	//sort density sticks according to the length
	EdgeStick tmpStick;
	for (i=0; i<allSticks.size (); i++)
		for (j=i+1; j<allSticks.size (); j++)
			if (allSticks[i].nAA > allSticks[j].nAA)
			{
				//swap
				tmpStick = allSticks[i];
				allSticks[i] = allSticks[j];
				allSticks[j] = tmpStick;
			}
	*/

	if (prnt)
	{
		printSticksInfo(allSticks, hlces.size(), strands.size());
	}
}

void printSticksInfo(vector<EdgeStick> allSticks, int NhlcesSticks, int NstrandsSticks)
{
	int i;

	//write information on screen
	cout<<"		======================================"<<endl;
	cout<<"			Density Sticks information"<<endl;
	cout<<"		======================================"<<endl;
	for (i=0; i<allSticks.size (); i++)
		cout<<"		"<<i+1<<". stickNum= "<<allSticks[i].stickNum <<" Type= "<<allSticks[i].ssType <<" vLength= "<<allSticks[i].nAA<<endl;

	cout<<endl<<"		# of Strands Sticks= "<<NstrandsSticks<<endl;
	cout<<"		# of Hlces Sticks= "<<NhlcesSticks<<endl;
	cout<<"		Total # of Sticks= "<<NhlcesSticks + NstrandsSticks<<endl;
	cout<<"		======================================"<<endl;
	cout<<"		======================================"<<endl<<endl;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//translate the secondary structures to overlap the sticks given by Edge files those stored in sticks data structure
void overlapSticks(vector<Protein> & sticksSS, vector<EdgeStick> sticks)
{

	Coordinate pnt1,pnt2,pnt3,pnt4;

	for (int i=0;i<sticksSS.size();i++)
	{
		Coordinate startEndCoord, endEndCoord;
		int tmpNumOfAA = sticksSS[i].numOfAA();

		pnt1 = sticks[i].start;
		pnt2 = sticks[i].end;

		if (sticks[i].ssType == 'H')
		{
			int AAsCounter = 0;	//for atoms will be involved in the center axial
			//take at most 3 atoms from each end
			while ((AAsCounter < tmpNumOfAA - 1) && (AAsCounter < 3))
			{
				Coordinate tmpCoord;

				//get coordinate of Ca atoms from the start end
				tmpCoord = sticksSS[i].getAtomCoordinate(AAsCounter," CA ");
				startEndCoord.x += tmpCoord.x;
				startEndCoord.y += tmpCoord.y;
				startEndCoord.z += tmpCoord.z;

				//get coordinate of Ca atoms from the end end
				tmpCoord = sticksSS[i].getAtomCoordinate(tmpNumOfAA - AAsCounter - 1," CA ");
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
			sticksSS[i].overlapLine(startEndCoord,endEndCoord,pnt1,pnt2);
		}
		else
		{
		// To translate beta strands.... we will try to overlap the center of the stick Strand (center of first AA atoms N-Ca-C) with the stick start and then
		//	rotate the structure so the center of the last AA in stick strand overlap the end of the stick
			startEndCoord.x = (sticksSS[i].getAtomCoordinate(0, " N  ").x + sticksSS[i].getAtomCoordinate(0, " C  ").x + sticksSS[i].getAtomCoordinate(0, " CA ").x) / 3;
			startEndCoord.y = (sticksSS[i].getAtomCoordinate(0, " N  ").y + sticksSS[i].getAtomCoordinate(0, " C  ").y + sticksSS[i].getAtomCoordinate(0, " CA ").y) / 3;
			startEndCoord.z = (sticksSS[i].getAtomCoordinate(0, " N  ").z + sticksSS[i].getAtomCoordinate(0, " C  ").z + sticksSS[i].getAtomCoordinate(0, " CA ").z) / 3;

			endEndCoord.x	= (sticksSS[i].getAtomCoordinate(tmpNumOfAA - 1, " N  ").x + sticksSS[i].getAtomCoordinate(tmpNumOfAA - 1, " C  ").x + sticksSS[i].getAtomCoordinate(tmpNumOfAA - 1, " CA ").x ) / 3;
			endEndCoord.y	= (sticksSS[i].getAtomCoordinate(tmpNumOfAA - 1, " N  ").y + sticksSS[i].getAtomCoordinate(tmpNumOfAA - 1, " C  ").y + sticksSS[i].getAtomCoordinate(tmpNumOfAA - 1, " CA ").y ) / 3;
			endEndCoord.z	= (sticksSS[i].getAtomCoordinate(tmpNumOfAA - 1, " N  ").z + sticksSS[i].getAtomCoordinate(tmpNumOfAA - 1, " C  ").z + sticksSS[i].getAtomCoordinate(tmpNumOfAA - 1, " CA ").z ) / 3;

			//move the strand to overlap the stick
			sticksSS[i].overlapLine(startEndCoord,endEndCoord,pnt1,pnt2);

		}
		Coordinate imiginaryNterminus = pointLineIntersection(sticksSS[i].getAtomCoordinate(0," N  "),pnt1,pnt2);				//for first atom
		Coordinate imiginaryCterminus = pointLineIntersection(sticksSS[i].getAtomCoordinate(tmpNumOfAA-1," C  "),pnt1,pnt2);	//for last atom


		Vectors imiginaryLine(imiginaryNterminus, imiginaryCterminus);
		double  distanceFromNterminus = getDistance(imiginaryNterminus,pnt1),
				imiginaryLength = imiginaryLine.length();


		/*****
				first translate the structure so Nterminus almost above (overlap) pnt1
		*****/
		int sign = 1;   //negative toward pnt2, positive to opposite direction

		//check the position of the Nterminus to the right or to the left of pnt1
		if (getDistance(imiginaryCterminus, pnt1) < imiginaryLength)
			sign = -1;

		double translation = distanceFromNterminus / imiginaryLength;		//how much to translate regards to the size of original vector

		imiginaryLine = imiginaryLine.mul(sign * translation);

		sticksSS[i].translateBy(imiginaryLine.getCoordinates());


		//check the position of the Cterminus to the right or to the left of pnt2
		imiginaryCterminus = pointLineIntersection(sticksSS[i].getAtomCoordinate(tmpNumOfAA-1," C  "),pnt1,pnt2);		//for last Ca
		imiginaryNterminus = pnt1;

		imiginaryLine.set(imiginaryNterminus, imiginaryCterminus);		//re-set the vector

		sign = 1;
		if (getDistance(imiginaryNterminus, pnt2) > imiginaryLength)
			sign = -1;

		translation = getDistance(imiginaryCterminus, pnt2);		//get the distance
		translation /= 2;											//divide over 2 .... the half-way
		translation /= imiginaryLength;								//regards to the original vector

		imiginaryLine = imiginaryLine.mul(sign * translation);

		sticksSS[i].translateBy(imiginaryLine.getCoordinates());
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void flipHelix(Protein & stickHlx, EdgeStick stick)
{

	Coordinate middleP,		//the point in the middle
			   origin;		//the Origin.....arbitrary point used to find a perpendicular line on stick later

	middleP.x = (stick.start .x + stick.end .x) / 2;
	middleP.y = (stick.start .y + stick.end .y) / 2;
	middleP.z = (stick.start .z + stick.end .z) / 2;

	Vectors v1(stick.start, middleP);	//the vector from middle point to the start end of the stick

	Vectors v2(origin, middleP);

	double angle = -3.14159	;	//angle is 180 degree ... to filp

	Vectors vNormal;		//the normal which is perpendicular on both lines
	vNormal = v2.cross(v1);
	vNormal += middleP;

	double mtx[4][4];
	//build rotation matrix to rotate a round the normal which is perpendicular on the stick and pass the middle point
	buildRotationMatrix(mtx, middleP, vNormal.getCoordinates(), -angle);

	//rotate all atoms
	for (int i=0;i<stickHlx.numOfAA();i++)
		for (int j=0;j<stickHlx.numOfAtoms(i);j++)
			rotatePoint(mtx, stickHlx.AAs[i].atoms[j].coord );
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// builds initial skeleton for the given edge file and native pdb file
//the hlces vector for the native pdb file should be built before calling this function
vector<Protein> getInitialSkeleton(Protein nativePDBFile, string mrcFile, vector<EdgeStick> & sticks, string outDir)
{

	Protein tmpSkeleton;
	vector<Protein> InitialSkeleton;
	int i, j;

	//count the number of stick hlces and the number of stick strands
	int nStickHlces = 0,
		nStickStrands = 0;

	for (i=0; i<sticks.size (); i++)
	{
		if (sticks[i].ssType == 'H')
			nStickHlces++;
		if (sticks[i].ssType == 'S')
			nStickStrands++;
	}
	//The pdb file should contain hlces or strands, in addition,  the number of sticks should be greater than or equal the number of SS
	if (((nativePDBFile.hlces.size()) || (nativePDBFile.sheets .size ())) &&
		(nStickHlces <= nativePDBFile.hlces.size()) &&
		(nStickStrands <= nativePDBFile.sheets .size()))
	{

		for (i=0;i<sticks.size();i++)
		{
			Protein tmp;
			AminoAcid tmpAA;

			if (sticks[i].ssType == 'H')
			{
				//generate stick hlx
				for (j=0; j<sticks[i].nAA; j++)
				{
					tmp.AAs .push_back (tmpAA);
					tmp.AAs [j].angles .phi = ALPHA_HLX_PHI_RIGHT;
					tmp.AAs [j].angles .psi = ALPHA_HLX_PSI_RIGHT;
				}
			}
			else
			{

				//generate stick strand
				for (j=0; j<sticks[i].nAA; j++)
				{
					tmp.AAs .push_back (tmpAA);
					tmp.AAs [j].angles .phi = BETA_SHEET_PHI;
					tmp.AAs [j].angles .psi = BETA_SHEET_PSI;
				}

			}

			//build the structure
			tmp.torsion2coord();

			InitialSkeleton.push_back(tmp);
			//InitialSkeleton[i].writeatoms("tmp.pdb");
		}

		//Transform the initial SS to overlap the sticks...
		overlapSticks(InitialSkeleton, sticks);


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
			stickName += "_";
			stickName += sticks[i].ssType;
			stickName += ".pdb";
			tmpHlx.writePDB (stickName, 1, 1);
		}


		return InitialSkeleton;
	}
	else
	{
		cout<<"=============================== In skeleton::SetInitialSkeleton ==============="<<endl;
		cout<<"The native pdb contains more number of sticks (hlces and/or sheets) than actual # of SS.."<<endl;
		cout<<"==============================================================================="<<endl;
		exit(0);
	}
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void buildSeqSS(Protein nativePDB, vector<SecondaryStruct> & seqSS, bool prnt)
{
	SecondaryStruct	tmpSS;

	short	i,
			hPntr = 0,
			sPntr = 0;

	while (hPntr < nativePDB.hlces.size() && sPntr < nativePDB.sheets.size())
	{
		//Theoritically....No overlap should occur
		if ((hPntr < nativePDB.hlces.size()) &&
			(sPntr < nativePDB.sheets.size()) &&
			(nativePDB.hlces[hPntr].startIndx == nativePDB.sheets[sPntr].startIndx))
		{
			cout<<"=================== in Skeleton::buildSeqSS(.........) ======================="<<endl;
			cout<<"Hlx ( "<<nativePDB.hlces[hPntr].serialNum<<" ) starts from the same AA that strand ( "
				<<nativePDB.sheets[sPntr].strandNum<<" ) starts from.."<<endl;
			cout<<"==============================================================================="<<endl;
			if (SHOW_ERRORS)
			{
				putchar(BEEP);
				cout<<"Press any key... The program will be halted"<<endl;
				getchar();
				exit(1);
			}
		}

		if (nativePDB.hlces[hPntr].startIndx < nativePDB.sheets[sPntr].startIndx)
		{
			tmpSS.type			= 'H';
			tmpSS.startIndx		= nativePDB.hlces [hPntr].startIndx;
			tmpSS.startAAnum	= nativePDB.AAs [tmpSS.startIndx ].num;
			tmpSS.endIndx		= nativePDB.hlces [hPntr].endIndx;
			tmpSS.endAAnum		= nativePDB.AAs [tmpSS.endIndx ].num ;
			tmpSS.nAA			= nativePDB.hlces [hPntr].nAA;
			tmpSS.order			= hPntr + 1;
			seqSS.push_back (tmpSS);
			hPntr++;
		}
		else
		{
			tmpSS.type			= 'S';
			tmpSS.startIndx		= nativePDB.sheets [sPntr].startIndx ;
			tmpSS.startAAnum	= nativePDB.AAs [tmpSS.startIndx ].num;
			tmpSS.endIndx		= nativePDB.sheets [sPntr].endIndx ;
			tmpSS.endAAnum		= nativePDB.AAs [tmpSS.endIndx ].num;
			tmpSS.nAA			= nativePDB.sheets [sPntr].nAA;
			tmpSS.order			= sPntr + 1;
			seqSS.push_back (tmpSS);
			sPntr++;
		}
	}

	if (hPntr < nativePDB.hlces.size())
	// strands end before hlces
	{
		//check with last strand
		//Theoritically....No overlap should occur
		if (nativePDB.sheets.size())
			if (nativePDB.hlces[hPntr].startIndx == nativePDB.sheets[sPntr-1].startIndx)
			{
				cout<<"=================== in Skeleton::buildSeqSS(.........) ======================="<<endl;
				cout<<"Hlx ( "<<nativePDB.hlces[hPntr].serialNum<<" ) starts from the same AA that strand ( "
					<<nativePDB.sheets[sPntr-1].strandNum<<" ) starts from.."<<endl;
				cout<<"==============================================================================="<<endl;
				if (SHOW_ERRORS)
				{
					putchar(BEEP);
					cout<<"Press any key... The program will be halted"<<endl;
					getchar();
					exit(1);
				}
			}

		for (i=hPntr; i<nativePDB.hlces.size(); i++)
		{
			tmpSS.type			= 'H';
			tmpSS.startIndx		= nativePDB.hlces [hPntr].startIndx;
			tmpSS.startAAnum	= nativePDB.AAs [tmpSS.startIndx ].num;
			tmpSS.endIndx		= nativePDB.hlces [hPntr].endIndx;
			tmpSS.endAAnum		= nativePDB.AAs [tmpSS.endIndx ].num ;
			tmpSS.nAA			= nativePDB.hlces [hPntr].nAA;
			tmpSS.order			= hPntr + 1;
			seqSS.push_back (tmpSS);
			hPntr++;
		}
	}

	if (sPntr < nativePDB.sheets.size())
	// hlces end before strands
	{
		//check with last hlx
		//Theoritically....No overlap should occur
		if (nativePDB.hlces.size())
			if (nativePDB.hlces[hPntr-1].startIndx == nativePDB.sheets[sPntr].startIndx)
			{
				cout<<"=================== in Skeleton::buildSeqSS(.........) ======================="<<endl;
				cout<<"Hlx ( "<<nativePDB.hlces[hPntr-1].serialNum<<" ) starts from the same AA that strand ( "
					<<nativePDB.sheets[sPntr].strandNum<<" ) starts from.."<<endl;
				cout<<"==============================================================================="<<endl;
				if (SHOW_ERRORS)
				{
					putchar(BEEP);
					cout<<"Press any key... The program will be halted"<<endl;
					getchar();
					exit(1);
				}
			}
		for(i=sPntr; i<nativePDB.sheets.size(); i++)
		{
			tmpSS.type			= 'S';
			tmpSS.startIndx		= nativePDB.sheets [sPntr].startIndx ;
			tmpSS.startAAnum	= nativePDB.AAs [tmpSS.startIndx ].num;
			tmpSS.endIndx		= nativePDB.sheets [sPntr].endIndx ;
			tmpSS.endAAnum		= nativePDB.AAs [tmpSS.endIndx ].num;
			tmpSS.nAA			= nativePDB.sheets [sPntr].nAA;
			tmpSS.order			= sPntr + 1;
			seqSS.push_back (tmpSS);
			sPntr++;
		}
	}
	if (prnt)
		printSSlist(nativePDB, seqSS);
}

//////////////////////////////////////////////////////////////////////////////
void printSSlist(Protein nativePDB, vector<SecondaryStruct> SSlist)
{
	int i;
	//write information on screen
	cout<<endl<<"		===================================================="<<endl;
	cout<<"		Secondary Structures...in Original order"<<endl;
	cout<<"		===================================================="<<endl;

	for (i=0; i<SSlist.size (); i++)
		cout<<"		"<<i+1<<". SS : "<<SSlist[i].type <<" "<<nativePDB.AAs[SSlist[i].startIndx].chr3<<SSlist[i].startAAnum
			<<"-"<<nativePDB.AAs[SSlist[i].endIndx].chr3<<SSlist[i].endAAnum<<" ( AAs indx "<<SSlist[i].startIndx
			<<"-"<<SSlist[i].startIndx + SSlist[i].nAA - 1<<" )  Length= "<<SSlist[i].nAA <<endl;

	cout<<endl<<"		# of Seq. Strands= "<<nativePDB.sheets.size()<<endl;
	cout<<"		# of Seq. Hlces= "<<nativePDB.hlces.size()<<endl;
	cout<<"		Total # of SS= "<<nativePDB.sheets.size() + nativePDB.hlces.size()<<endl;
	cout<<"		===================================================="<<endl;
	cout<<"		===================================================="<<endl<<endl;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// hlces should be assigned to sticks first...sticks.ssAssigned should be set before
//this function is supposed to use to determine the maximum of shift could be applied to the stick hlx over the native sequence
//validity right means that the stick hlx could or not be shifted to the right (toward the end of the sequence)...
//validity left means that the stick hlx could be or not be shifted to the left (toward the beginning of the sequence)
void setShiftValidity(Protein nativePDBFile, vector<EdgeStick> & sticks)
{
	if ((sticks.size()) && (sticks[0].ssAssigned != -1))
	{
		for(int i=0;i<sticks.size();i++)
		{
			sticks[i].validityRight = 0;    //reset the value
			sticks[i].validityLeft  = 0;	 //reset the value

			int indx = sticks[i].ssAssigned;
			int nativePDBNumOfAA = nativePDBFile.numOfAA();


			int numOfAADiff = nativePDBFile.hlces[indx].nAAcur - sticks[i].nAA;

			numOfAADiff /= 2;			//the shift amount left or right to centerally assign the hlx


			if (sticks[i].direction == 1)
			{
				sticks[i].startAAIndx = nativePDBFile.hlces[indx].startIndx + numOfAADiff;

				if (sticks[i].startAAIndx < 0)		//make it start from the beginning
					sticks[i].startAAIndx = 0;

				if (sticks[i].startAAIndx + sticks[i].nAA > nativePDBNumOfAA)						//make it start accordingly
					sticks[i].startAAIndx -= (sticks[i].startAAIndx + sticks[i].nAA) -  nativePDBNumOfAA;

				sticks[i].validityLeft = sticks[i].startAAIndx;	//the number of AA left of the stick if assigned to this hlx
				sticks[i].validityRight = nativePDBNumOfAA - (sticks[i].startAAIndx + sticks[i].nAA);	//the number of AA right to the stick if assigned to this hlx
			}
			else
			{
				sticks[i].startAAIndx = nativePDBFile.hlces[indx].endIndx - numOfAADiff;
				//cout<<"startAAindx = "<<sticks[i].startAAIndx<<" "<<nativePDBFile.AAs[indx].chr3<<" "<<nativePDBFile.AAs[indx].num<<endl;

				if (sticks[i].startAAIndx >= nativePDBNumOfAA)				//make it start in the end of the pdb
					sticks[i].startAAIndx = nativePDBNumOfAA - 1;

				if (sticks[i].startAAIndx - sticks[i].nAA < -1)		//make it start accordingly
					sticks[i].startAAIndx += sticks[i].nAA - sticks[i].startAAIndx - 1;

				sticks[i].validityLeft = sticks[i].startAAIndx - sticks[i].nAA + 1;
				sticks[i].validityRight = nativePDBNumOfAA - sticks[i].startAAIndx - 1;
			}


		//	cout<<sticks[i].stickNum<<" with hlx"<<sticks[i].ssAssigned+1<<"  Direction "<<sticks[i].direction		\
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
/*
	int startIndx = sticks[stickIndx].startAAIndx + shift;		//positive to the right and negative to the left

	if (sticks[stickIndx].direction == 1)
	{
		for (int i=0; i<sticks[stickIndx].nAA; i++)
		{
			portion.renameAA(i, nativePDBFile.AAs[startIndx].chr3);
			startIndx += 1;
		}
	}
	else
	{
		for (int i=sticks[stickIndx].nAA-1; i>=0; i--)
		{
			portion.renameAA(i, nativePDBFile.AAs[startIndx].chr3);
			startIndx -= 1;
		}
	}
*/

	int startIndx = sticks[stickIndx].startAAIndx + shift;
	for (int i=0 ; i<sticks[stickIndx].nAA; i++)
	{
		portion.renameAA(i, nativePDBFile.AAs[startIndx].chr3);
		startIndx++;
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
	for (i=0; i < pow(2.0, nSticks); i++)
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
///////////////////////////////////////////////////////////////////////////////////////////////////////
void translateEdge(vector<Coordinate> &edge, float tAmount){
    int i, nPnts;
    Coordinate c1, c2, pnt;
    if (edge.size()>2){
        //find the center of the first half
        nPnts = (int) (edge.size()/2);
        for (i=0; i<nPnts; i++){
            c1.x += edge[i].x;
            c1.y += edge[i].y;
            c1.z += edge[i].z;
        }
        c1.x /= nPnts;
        c1.y /= nPnts;
        c1.z /= nPnts;
        //find the center of the second half
        for (i=nPnts; i<edge.size(); i++){
            c2.x += edge[i].x;
            c2.y += edge[i].y;
            c2.z += edge[i].z;
        }
        c2.x /= (nPnts + 1);
        c2.y /= (nPnts + 1);
        c2.z /= (nPnts + 1);

        //find the point on the line c1-c2
        if (tAmount>0){      //translate to the right toward c2
            pnt = pointOnLine(c1, c2, getDistance(c1,c2) + tAmount);
            //get translaion on each axis
            pnt.x -= c2.x;
            pnt.y -= c2.y;
            pnt.z -= c2.z;
        }
        else{       //translate to the left toward c1
            pnt = pointOnLine(c2,c1, getDistance(c1,c2) + tAmount);
            //get translaion on each axis
            pnt.x -= c1.x;
            pnt.y -= c1.y;
            pnt.z -= c1.z;
        }
        //apply translation
        for (i=0; i<edge.size(); i++){
            edge[i].x += pnt.x;
            edge[i].y += pnt.y;
            edge[i].z += pnt.z;
        }
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
void generateValidTopologies(Protein nativePDBFile, vector<Protein> initialSkeleton, vector<SecondaryStruct> seqSS, vector<int> permutation, vector<EdgeStick> sticks, vector<sticksRec> & validTopologies)
{
		bool validPermutation = true;
		double lengthVariationHelix = 0.5;			//the maximum of length variation could occure between stick and sequence segment
		int j, k;
		/*****
				check the first condition, the variation on length should not be less than a predefined threshold
		*****/
		for (j=0; j<sticks.size(); j++)
		{
			cout<<seqSS[permutation[j]-1].nAA<<" "<<sticks[j].nAA<<endl;
			if (((double) seqSS[permutation[j] - 1].nAA/sticks[j].nAA < lengthVariationHelix) ||
				((double) sticks[j].nAA /seqSS[permutation[j]-1].nAA < lengthVariationHelix ))
				validPermutation = false;			//not valid permutation
		}
		cout<<"Done.."<<endl;

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
				sticks[j].ssAssigned = permutation [j] - 1;
			//	cout<<endl<<"Stick "<<sticks[j].stickNum<<" with hlx "<<sticks[j].ssAssigned +1<<"   stick length= "<<sticks[j].length \
					<<" hlx actual length= "<<nativePDBFile.hlces[sticks[j].ssAssigned].length<<" startAAindx= "<<sticks[j].startAAIndx<<endl;
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
						if (tmpSticks[m+1].ssAssigned < tmpSticks[m].ssAssigned)
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
							loopLength += tmpSticks[k+1].startAAIndx - tmpSticks[k].startAAIndx - tmpSticks[k].nAA;
						else
							loopLength += tmpSticks[k+1].startAAIndx - tmpSticks[k+1].nAA - tmpSticks[k].startAAIndx - tmpSticks [k].nAA + 1;
					}
					else
					{
						if (tmpSticks[k+1].direction == 1)
							loopLength += tmpSticks[k+1].startAAIndx - tmpSticks[k].startAAIndx - 1;
						else
							loopLength += tmpSticks[k+1].startAAIndx - tmpSticks[k+1].nAA - tmpSticks[k].startAAIndx;
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
					validTopologies.push_back(allTopologies[j]);
			}
		}
}
//////////////////////////////////////////////////////////////////////////////////////////////////////

void genrateAllValidTopologies(Protein nativePDBFile,
							   vector<Protein> initialSkeleton,
							   vector<OnePermutation> hlcesPermutations,
							   vector<OnePermutation> strandsPermutations,
							   vector<EdgeStick> sticks)
{
	//long int indxJump = 1;//factorial(nativePDBFile.hlces.size() - sticks.size());  was in the previous version
//	long int i = allPermutations.size() - 1;
	/*****
		Determine Valid Topologies
			geometry screening used to determine the valid topologies if satisfies two conditions
			1. The difference in length should be less than 7
			2. The length of the loop in between should not exceed the length of extended loop (3.8 * num. of AA's in the loop)
	*****/
//	while (i >= 0)
//	{
//		getNextPermutation(hlcesPermutations, strandsPermutations, permutation);
//		generateValidTopologies(nativePDBFile, initialSkeleton, allPermutations[i].permutation, sticks);
//		allPermutations.pop_back ();			//delete last permutation...we are done with it...to save memory for big proteins
//		i--;
//	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// End Of IMPLEMENTATION ////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif
