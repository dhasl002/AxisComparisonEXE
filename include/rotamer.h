#ifndef ROTAMER_H
#define ROTAMER_H

#include <vector>
#include "constants.h"
#include "protein.h"
#include<fstream>


#include<iostream>
#include<string>
#include <stdio.h>
#include <cmath>


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// List Of Data Structures //////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// The Rotamer Used is Penultimate Rotamer Library
// http://kinemage.biochem.duke.edu/databases/rotamer.php

#define ROTAMER_FILE "rotamer.txt"					// the file that contains the generic information "rotamers" for AA's
#define SIDE_CHAINS_DIRECTORY	"AATemplates/"		//the directory that contains the AA's structures

struct chi
{
	double chi1;
	double chi2;
	double chi3;
	double chi4;

	double prob;
	//initializer
	chi() : chi1(999.00), chi2(999.00), chi3(999.00), chi4(999.00), prob(0) {}
};

struct chiInfo
{
	vector<chi> chiVect;		//all chi that could be in an AA
	int numOfChi;				//number of chi in the AA

	//initializer
	chiInfo() : numOfChi(0) {}
};

// The next two data structures used for BBdep library
struct bins
{
	vector<chi> rotamer;
};
struct AARotamer
{
	vector<bins> bin;
	int numOfChiAngles;

	//initializer
	AARotamer() : numOfChiAngles(0) {}
};

//data structure that contains the template side chains for all standard AA's
static Protein sideChains[STANDARD_AMINO_ACIDS];		

//Global Variable has the rotamer information needed (used for generic rotamer)
static vector <chiInfo> rotVect;

//Global Variable has rotamer information needed (used for BBDep rotamer library)
static vector<AARotamer> BBDepRot;


/*****
		BBDepRot Structure
			
			|		.		|		.			||		.	.
			|		.		|		.			||		.	.
			|		.		|		.			||		.	.
AARotamer 1	|		.	   ||		.			||		.	.
			|______________||___________________||___________________||________		
			|			   ||	............	||	...............	 ||
			|			   ||	...........		||	...............	 ||..........
AARotamer 0	|numOfChiAngles||___________________||___________________||
			|			   ||chi1|chi2|chi3|chi4||chi1|chi2|chi3|chi4||	............. 				
			|			   ||____|____|____|____||____|____|____|____|| .........
			|			   ||chi1|chi2|chi3|chi4||chi1|chi2|chi3|chi4|| ..........
			|______________||____|____|____|____||____|____|____|____||____________________
									bin 0				bin 1				bin 2
								-180	-180		-180	-170		-180	-160	
*****/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// End Of Constants and Data Structures /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// List of Functions ////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//takes the name of the file contains the rotamer info
//rotamer .txt file format
// ABC A #     //the name of the aa (3 letters),space, 1 letter name, space,  and then the number of chi angles
//      _-_- _-_- _-_-        //5 spaces then 4 spaces to the number (the first one for the sign and the last three for the number represents the angle)
//      -140  120  60         //example
//       20  -30  -10         //example

//notice: the indx of rotVect as in maptonum function in protein.h for example (also as in rotamer.txt file
// 0;	//ALA
// 1;	//ARG
// 2;	//ASN
// 3;	//ASP
// 4;	//CYS
void buildRotamerVector();

void buildBBDepRotamer(string folderPath, double threshold = 0.9);			//build the BBDepRotamer ... threshold used to determine how many rotamer to be included for each AA
int getBinIndx(double phi, double psi);			//given the angles phi and psi returns the indx of the bin in the AA
void writeBBDepRot(string outFolder);		//write all rotamers to files....given the out folder path	

void applySideChainRotamer(Protein & portion, int AAIndx, int rotIndx, int * & lastChiUsed); //rotate a given AA side chain according a rotamer for one time

//resolves the collision in the pdb protein
//takes the protein pdb, a vector contains 2 int represent the indx of atoms (atom_vect indx) that collide, and an initialized lastchi array
//returns the protein without collision saved in a new pdb class
//uses a (dynamic array) laschiused that keeps track of last chi (chi angle number) has been used in each aa
//notice: lastchiused and rot_vect should be built (initialized); before calling this function

// 3 cases
// 1: the atoms pair collide are both from side chain...uses lastchiused and rot_vect
// 2: the atoms pair collide are one from a side chain and the other from a backbone .... uses lastchiuse and rot_vect
// 3: the atoms pair collide are both from backbone ..... generates a random phi or psi angle for the backbone within the predefined regaes
int resolveCollision(Protein & portion, int * & lastChiUsed);		//resolves all collisions in a portion of AA's..returns 0 if failed to resolve all collisions

void initializeLastChiUsed(int * &lastChiUsed, int size);	//initiate the array of last chi used that is used in resolveCollision()
void createSideChainsArray();			//builds the array contains the template side chains for all standard AA

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// ////////// End Of The List //////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// IMPLEMENTATION ///////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void buildRotamerVector()
{

	ifstream inFile;
	string line;
	chiInfo tmpChiVectForAA;		//a data structure to store the number of chi and chi information for each AA
	chi  tmpChiRec;					//a record contains the information of chi angles for each AA

	inFile.open(ROTAMER_FILE);		//rotFile : the file contains the info (rotamer ) for each AA
	if (!inFile)
	{
		cout<<"======================== in rotamer::buildRotamerVector(string) ==============="<<endl;
		cerr << "Unable to open " << ROTAMER_FILE << endl; 
		cout<<"==============================================================================="<<endl;
		exit(1);
	}
	else
	{
		tmpChiVectForAA.numOfChi = 0;
		tmpChiVectForAA.chiVect.clear();
		
		//reads first aa info .... get information for ALA
		getline(inFile,line);
		tmpChiVectForAA.numOfChi = atoi(line.substr(6,1).c_str());		// or can be obtained from getNumOfChi(char)
		while (!inFile.eof())
		{
			getline(inFile,line);
			if (line[0]!=' ')
			{
				rotVect.push_back(tmpChiVectForAA);				//commit the previous AA info
				tmpChiVectForAA.chiVect.clear();				//clear the data structure to accept new info for a new AA
				tmpChiVectForAA.numOfChi = atoi(line.substr(6,1).c_str());		// get the number of chi angles for this new AA
			}
			else
			{
				//get the information of chi angles for AA
				if (tmpChiVectForAA.numOfChi == 1)
				{
					tmpChiRec.chi1 = atoi(line.substr(5,4).c_str());	//other chi angles are not exist....initialized to be 999.00
					tmpChiRec.chi2 = 999;		//does not exist
					tmpChiRec.chi3 = 999;		//does not exist
					tmpChiRec.chi4 = 999;		//does not exist					
					tmpChiVectForAA.chiVect.push_back(tmpChiRec);
				}

				if (tmpChiVectForAA.numOfChi == 2)
				{
					tmpChiRec.chi1 = atoi(line.substr(5,4).c_str());
					tmpChiRec.chi2 = atoi(line.substr(10,4).c_str());
					tmpChiRec.chi3 = 999;		//does not exist
					tmpChiRec.chi4 = 999;		//does not exist
					tmpChiVectForAA.chiVect.push_back(tmpChiRec);
				}
				if (tmpChiVectForAA.numOfChi == 3)
				{
					tmpChiRec.chi1 = atoi(line.substr(5,4).c_str());
					tmpChiRec.chi2 = atoi(line.substr(10,4).c_str());
					tmpChiRec.chi3 = atoi(line.substr(15,4).c_str());
					tmpChiRec.chi4 = 999;		//does not exist
					tmpChiVectForAA.chiVect.push_back(tmpChiRec);
				}
				if (tmpChiVectForAA.numOfChi == 4)
				{
					tmpChiRec.chi1 = atoi(line.substr(5,4).c_str());
					tmpChiRec.chi2 = atoi(line.substr(10,4).c_str());
					tmpChiRec.chi3 = atoi(line.substr(15,4).c_str());
					tmpChiRec.chi4 = atoi(line.substr(20,4).c_str());
					tmpChiVectForAA.chiVect.push_back(tmpChiRec);
				}

			}

		}
		rotVect.push_back(tmpChiVectForAA);
		inFile.close();

		//to print rotamer vector
		/*
		for (int i=0;i<rotVect.size();i++){
			cout<<"=============================="<<endl;
			cout<<" i = "<<i<<"  numOfChi= "<<rotVect[i].numOfChi<<endl;
			cout<<"=============================="<<endl;
			for (int j=0;j<rotVect[i].chiVect.size();j++)
				cout<<"j= "<<j<<" "<<rotVect[i].chiVect[j].chi1<<"  "<<rotVect[i].chiVect[j].chi2
					<<"  "<<rotVect[i].chiVect[j].chi3<<"  "<<rotVect[i].chiVect[j].chi4<<endl;
		}
		*/
		
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//folder path is the folder where the files contains rotamers is located
//threshold is used to know how many rotamer for each AA to be included
//it is the accumulative sum of the probablity of the rotamer...the default is 0.9
void buildBBDepRotamer(string folderPath, double threshold)
{

	int i;
	string line;

//	cout<<"======================== in rotamer::buildBBDepRotamer(string) ================"<<endl;
	/*****
			build Rotamers of ALL standard AA's
	*****/
	for (i= 0 ;i<STANDARD_AMINO_ACIDS; i++)		//starts from 1 b/s ALA indexed 0 has no rotamer
	{
		chi  tmpChiRec;					//a record contains the information of chi angles for each AA
		bins tmpBin;
		AARotamer tmpAARot;
		string AAchr3 = num2Chr3(i);
		

//		cout<<"      building Rotamer Library for ( "<<AAchr3<<" ).";
		
		if ((AAchr3 != "ALA") && (AAchr3 != "GLY"))		//b/s ALA and GLY have no rotamers
		{
			/*****
					set the path where the rotamer file is located
			*****/
			string AAFilePath = folderPath + "/";		
			AAFilePath += AAchr3;
			AAFilePath +="_Rotamers.txt";		//now the full path is in AAFilePath
			/*****
					open the file and start parsing it
			*****/
			ifstream inFile;
			inFile.open (AAFilePath.c_str ());
			if (!inFile)
			{
				cout<<"======================== in rotamer::buildBBDepRotamer(string) ================"<<endl;
				cerr << "Unable to open " << AAFilePath << endl; 
				cout<<"==============================================================================="<<endl;
				exit(1);
			}
			else
			{
				int binIndx = 0;				//initial binIndx
				double currProb = 0;			//wht's the Accumulative probability of the rotamers u have parsed so far
				while (!inFile.eof())
				{
					getline(inFile,line);
					if (line.length() > 0)
					{
						//get the current bin (the indx in the BBDepRot)
						int tmpBinIndx = getBinIndx(atoi(line.substr (5,4).c_str ()),atoi(line.substr (10,4).c_str ()));

						//r u still reading in the same bin or not
						if (tmpBinIndx != binIndx)
						{
							/*****
									I have finished working with the bin... commit it's information
							*****/
							//push the information of the current bin
							tmpAARot.bin.push_back (tmpBin);
							
							//clear previous rotamers.. for the previous bin
							tmpBin.rotamer.clear ();

							//ready to work on the next bin
							binIndx++;

							//reset currProb
							currProb = 0;

							//save the first (line) information of the next bin
							tmpChiRec.chi1 = atof(line.substr (43,6).c_str ());
							tmpChiRec.chi2 = atof(line.substr (51,6).c_str ());
							tmpChiRec.chi3 = atof(line.substr (59,6).c_str ());
							tmpChiRec.chi4 = atof(line.substr (67,6).c_str ());

							tmpChiRec.prob = atof(line.substr (32,8).c_str ());

							tmpBin.rotamer.push_back (tmpChiRec);

							currProb +=  tmpChiRec.prob;
						}
						else
						{
							if (currProb < threshold)
							{
								/*****
										I am still reading information for the current bin ... save information
								*****/
								tmpChiRec.chi1 = atof(line.substr (43,6).c_str ());
								tmpChiRec.chi2 = atof(line.substr (51,6).c_str ());
								tmpChiRec.chi3 = atof(line.substr (59,6).c_str ());
								tmpChiRec.chi4 = atof(line.substr (67,6).c_str ());

								tmpChiRec.prob = atof(line.substr (32,8).c_str ());

								tmpBin.rotamer.push_back (tmpChiRec);

								//update the accumulative probability
								currProb +=  tmpChiRec.prob;
							}
						}
					}
				}
				/*****
						Done with current AA ... commit it's all bins and clear rotamer to start working on the next AA
				*****/
				//push the information of the last bin
				tmpAARot.bin.push_back (tmpBin);

				//reset chi angles in the bin
				tmpBin.rotamer.clear ();


				//push the information of all bins
				BBDepRot.push_back (tmpAARot);

				//get The number of chi and store it
				if (BBDepRot[i].bin[0].rotamer[0].chi2 == 0)
					BBDepRot[i].numOfChiAngles = 1;
				else
					if (BBDepRot[i].bin[0].rotamer[0].chi3 == 0)
						BBDepRot[i].numOfChiAngles = 2;
					else
						if (BBDepRot[i].bin[0].rotamer[0].chi4 == 0)
							BBDepRot[i].numOfChiAngles = 3;
						else
							BBDepRot[i].numOfChiAngles = 4;

				inFile.close();

			}
		}
		else		//to keep the indexing method when u get ALA or GLY push an empty AARotamer structure

		{
			//push empty structures for AA's have no side chain rotamers
			tmpAARot.bin .clear ();
			BBDepRot.push_back (tmpAARot);
		}
//		cout<<"       num. of Chi angles = "<<BBDepRot[i].numOfChiAngles <<endl;
	}
//	cout<<"======================== End of rotamer.buildBBDepRot (string) ================"<<endl;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int getBinIndx(double phi, double psi)
{
	phi /= 10;

	if (phi >= 0)
		phi = int (phi + 0.5);
	else
		phi = int (phi - 0.5);

	phi += 18;


	psi /= 10;	
	if (psi >= 0)
		psi = int (psi + 0.5);
	else
		psi = int (psi - 0.5);

	psi += 18;

/*
	phi += 180;
	phi /= 10;		//convert phi to a number between 0 and 18
	psi += 180;
	psi /= 10;		//convert psi to a number between 0 and 18
//	cout<<"phi= "<<phi<<" psi= "<<psi<<" ";
*/

	return (psi + (phi * 37));
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void writeBBDepRot(string outFolder)
{
	for (int i= 0; i<BBDepRot.size(); i++)
	{
		// set the path and the name of the AA rotamers to be written
		string outFilePath = outFolder;
		outFilePath += "\\";
		outFilePath += num2Chr3(i);
		outFilePath += "_Rotamers_out.txt";

		ofstream out;
		out.open(outFilePath.c_str());
		if (!out) 
		{
			cout<<"========================= in Rotamer.writeBBDepRot(string) ===================="<<endl;
			cout<< "Unable to open " << outFilePath << endl; 
			cout<<"==============================================================================="<<endl;
			if (SHOW_ERRORS)
			{
				putchar(BEEP);
				cout<<"Press any key..."<<endl;
				getchar();
			}
		}
		else
		{
			out<<setw(7)<<left<<"//Bin#"<<setw(9)<<"chi1"<<setw(9)<<"chi2"<<setw(9)<<"chi3"<<setw(9)<<"chi4"<<setw(9)<<"Prob. "<<"  number of chi angles = "<<BBDepRot[i].numOfChiAngles<<endl;
			for (int j=0; j<BBDepRot[i].bin.size(); j++)
			{
				for (int k=0; k<BBDepRot[i].bin [j].rotamer.size(); k++)
					out<<setw(7)<<left<<j<<setw(9)<<BBDepRot[i].bin[j].rotamer[k].chi1<<setw(9)<<BBDepRot[i].bin[j].rotamer[k].chi2
						<<setw(9)<<BBDepRot[i].bin[j].rotamer[k].chi3<<setw(9)<<BBDepRot[i].bin[j].rotamer[k].chi4<<setw(9)<<BBDepRot[i].bin [j].rotamer [k].prob<<endl;
			}
			out.close ();
		}
		
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void applySideChainRotamer(Protein & portion, int AAIndx, int rotIndx, int * & lastChiUsed)
{
	lastChiUsed[AAIndx] += 1;

	double chiAngle = portion.getChi (AAIndx, 1);	//get chi 1 for the second AA

	
	if ((chiAngle != 999) && (rotVect[rotIndx].numOfChi))
	{
//		cout<<"working on AA : "<<portion.AAs[AAIndx].chr3<<" "<<portion.AAs [AAIndx].num \
			<<" rotIndx = "<<rotIndx<<"  "<<lastChiUsed[AAIndx]<<" "<<rotVect[rotIndx].chiVect.size()<<endl;
		double deltaChi = chiAngle - rotVect[rotIndx].chiVect[lastChiUsed[AAIndx]].chi1;
		Coordinate p1,
					p2;
		///////////////////
		//  rotate
		///////////////////
//		if (abs(deltaChi) > validError)
//		{
			p1 = portion.getAtomCoordinate (AAIndx," CA"); 
			p2 = portion.getAtomCoordinate (AAIndx," CB");

			//cout<<"Ca.x = "<<p2.x <<" Ca.y = "<<p2.y <<" Ca.z = "<<p2.z<<endl;
			portion.rotateSideChain (AAIndx, portion.getAtomIndx (AAIndx, " CB ")+1, p1, p2, deltaChi);
//		}
		if (rotVect[rotIndx].numOfChi >= 2)
		{
			chiAngle = portion.getChi (AAIndx, 2);	//get the second chi of the second AA
			deltaChi = chiAngle - rotVect[rotIndx].chiVect[lastChiUsed[AAIndx]].chi2 ;
			///////////////
			///rotate
			///////
//			if (abs(deltaChi)>validError)
//			{
				p1 = portion.getAtomCoordinate (AAIndx,"G");   //0 b/s jus 1 (one) aa is there
				portion.rotateSideChain (AAIndx, portion.getAtomIndx (AAIndx, "G")+1, p2, p1, deltaChi);
//			}


			
			if (rotVect[rotIndx].numOfChi >= 3)
			{
				chiAngle = portion.getChi (AAIndx, 3);	//get the third chi of the second AA
				deltaChi = chiAngle - rotVect[rotIndx].chiVect[lastChiUsed[AAIndx]].chi3 ;
				////////////////
				///rotate
				///////
//				if (abs(deltaChi)>validError)
//				{
					p2 = portion.getAtomCoordinate (AAIndx,"D");
					portion.rotateSideChain (AAIndx, portion.getAtomIndx (AAIndx, "D")+1, p1, p2, deltaChi);
//				}

				
				if (rotVect[rotIndx].numOfChi >= 4)
				{
					chiAngle = portion.getChi (AAIndx, 4);	//get the fourth chi angle of the second AA
					deltaChi = chiAngle - rotVect[rotIndx].chiVect[lastChiUsed[AAIndx]].chi4 ;
					//////////////
					///rotate
					/////////
//					if (abs(deltaChi)>validError)
//					{
						p1 = portion.getAtomCoordinate (AAIndx,"E");   
						portion.rotateSideChain (AAIndx, portion.getAtomIndx (AAIndx, "E")+1, p2, p1, deltaChi);
//					}
				}
			}
		}
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//return 0 in the case that no rotamer angles still un-evaluated
int resolveCollision(Protein & portion, int * & lastChiUsed)
{

	if ((!rotVect.size()))
	{
		cout<<"============== in rotamer::resolveCollision(Protein &, int * &) ============"<<endl;
		cout<<"Please build rotVect or lastChiUsed first. . . ."<<endl;
		cout<<"============================================================================"<<endl;
		exit(1);
	}

	while (portion.doesCollide ())
	{
		int aa1indx		= portion.AAsCollide[0],		//get the indx of the first AA collides
			aa2indx		= portion.AAsCollide[1],		//get the indx of the second AA collides
			atom1Indx	= portion.atomsCollide [0],		//get the indx of the first atom collides
			atom2Indx	= portion.atomsCollide [1];		//get the indx of the second atom collides
		
		int rotIndx;		//the indx of the rotamer....which AA r u working on

		if (portion.AAs [aa1indx].atoms [atom1Indx].isSideChain) // is the first atom collides a side chain atom
		{
			if (portion.AAs [aa2indx].atoms [atom2Indx].isSideChain)	//is the second atomcollides a side chain atom
			{
			//Side chain with Side chain .... case 1
				rotIndx = chr2Num(portion.AAs [aa2indx].chr1);	//I gonna work on the second AA
				if ((lastChiUsed[aa2indx] != rotVect[rotIndx].chiVect.size() - 1) && 
					(rotVect[rotIndx].numOfChi))  //not finish....u still have angles to be tested in the rotamer (exclude GLY and ALA)
				{
					
					applySideChainRotamer(portion, aa2indx, rotIndx, lastChiUsed);
				}
				else
				{
					//the rotamer of the second AA is finished but u still have collision....try the first AA now
					lastChiUsed[aa2indx] = -1;		//re-initiate the rotamer of the second AA
					applySideChainRotamer(portion, aa2indx, rotIndx, lastChiUsed);	// return the original side chain for the second AA

					rotIndx = chr2Num(portion.AAs [aa1indx].chr1); //get the rotamer indx of the first AA						
					if ((lastChiUsed[aa1indx] != rotVect[rotIndx].chiVect.size()-1) &&
						(rotVect[rotIndx].numOfChi))//if u still have angles to be tested
					{
						applySideChainRotamer(portion, aa1indx, rotIndx, lastChiUsed);

					}
					else
					{
						cout<<"============= In rotamer:: ResolveCollision(Protein &, int * &) =============="<<endl;
						cout<<"All possible angles for current pair have been used....could not resolve the collision.."<<endl;
						cout<<"=============================================================================="<<endl;
						return 0;		//success (no collision) = 1...fail = 0
					}
				}

			}
			else
			{
			//Side chain with BackBone ... case 2
				rotIndx = chr2Num(portion.AAs[aa1indx].chr1);	//get the indx of the first AA
				
				if ((lastChiUsed[aa1indx] != rotVect[rotIndx].chiVect.size()-1) &&
					(rotVect[rotIndx].numOfChi))
				{
					applySideChainRotamer(portion, aa1indx, rotIndx, lastChiUsed);
				}
				else
				{
					cout<<"============= In rotamer:: ResolveCollision(Protein &, int * &) =============="<<endl;
					cout<<"All possible angles for current pair have been used....could not resolve the collision.."<<endl;
					cout<<"=============================================================================="<<endl;
					return 0;   //success (no collision) = 1 ... fail = 0
				}

			}
		}
		else
		{
			if (portion.AAs[aa2indx].atoms [atom2Indx].isSideChain)
			{
			//BackBone with Side chain ... case 2
				rotIndx = chr2Num(portion.AAs[aa2indx].chr1); // get the rotamer indx of the second AA
				
				if ((lastChiUsed[aa2indx] != rotVect[rotIndx].chiVect.size()-1) &&
					(rotVect[rotIndx].numOfChi))
				{
					applySideChainRotamer(portion, aa2indx, rotIndx, lastChiUsed);
				}
				else
				{
					cout<<"============= In rotamer:: ResolveCollision(Protein &, int * &) =============="<<endl;
					cout<<"All possible angles for current pair have been used....could not resolve the collision....."<<endl;
					cout<<"=============================================================================="<<endl;
					return 0; //success (no collision) = 1 ... fail = 0
				}

			}
			else
			{
			//BackBone with BackBone....case 3
				cout<<".....The collision is between two backbone atoms....( "<<portion.AAs[aa1indx].chr3<<" "<<portion.AAs[aa1indx].num<<" ) with ( " \
					<<portion.AAs [aa2indx].chr3<<" "<<portion.AAs[aa2indx].num<<" )"<<endl;
				//for now
				return 0;
			}
		}

		//clear previous AA's and atoms were collide
		portion.AAsCollide .clear ();
		portion.atomsCollide .clear ();
	}

	return 1;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void applySideChainRotamer_BBDep(Protein & portion, int AAIndx, int rotIndx, int binIndx, int * & lastChiUsed)
{
	lastChiUsed[AAIndx] += 1;

	double chiAngle = portion.getChi (AAIndx, 1);	//get chi 1 for the second AA

	
	if ((chiAngle != 999) && (BBDepRot[rotIndx].numOfChiAngles))
	{
//		cout<<"working on AA : "<<portion.AAs[AAIndx].chr3<<" "<<portion.AAs [AAIndx].num \
			<<" rotIndx = "<<rotIndx<<"  "<<lastChiUsed[AAIndx]<<" "<<rotVect[rotIndx].chiVect.size()<<endl;
		double deltaChi = chiAngle - BBDepRot[rotIndx].bin[binIndx].rotamer[lastChiUsed[AAIndx]].chi1;
		Coordinate p1,
					p2;
		///////////////////
		//  rotate
		///////////////////
//		if (abs(deltaChi) > validError)
//		{
			p1 = portion.getAtomCoordinate (AAIndx," CA"); 
			p2 = portion.getAtomCoordinate (AAIndx," CB");

			//cout<<"Ca.x = "<<p2.x <<" Ca.y = "<<p2.y <<" Ca.z = "<<p2.z<<endl;
			portion.rotateSideChain (AAIndx, portion.getAtomIndx (AAIndx, "CB")+1, p1, p2, deltaChi);
//		}
		if (BBDepRot[rotIndx].numOfChiAngles >= 2)
		{
			chiAngle = portion.getChi (AAIndx, 2);	//get the second chi of the second AA
			deltaChi = chiAngle - BBDepRot[rotIndx].bin[binIndx].rotamer[lastChiUsed[AAIndx]].chi2 ;
			///////////////
			///rotate
			///////
//			if (abs(deltaChi)>validError)
//			{
				p1 = portion.getAtomCoordinate (AAIndx,"G");   //0 b/s jus 1 (one) aa is there
				portion.rotateSideChain (AAIndx, portion.getAtomIndx (AAIndx, "D"), p2, p1, deltaChi);
//			}


			
			if (BBDepRot[rotIndx].numOfChiAngles >= 3)
			{
				chiAngle = portion.getChi (AAIndx, 3);	//get the third chi of the second AA
				deltaChi = chiAngle - BBDepRot[rotIndx].bin [binIndx].rotamer [lastChiUsed[AAIndx]].chi3 ;
				////////////////
				///rotate
				///////
//				if (abs(deltaChi)>validError)
//				{
					p2 = portion.getAtomCoordinate (AAIndx,"D");
					portion.rotateSideChain (AAIndx, portion.getAtomIndx (AAIndx, "E"), p1, p2, deltaChi);
//				}

				
				if (BBDepRot[rotIndx].numOfChiAngles >= 4)
				{
					chiAngle = portion.getChi (AAIndx, 4);	//get the fourth chi angle of the second AA
					deltaChi = chiAngle - BBDepRot[rotIndx].bin [binIndx].rotamer [lastChiUsed[AAIndx]].chi4 ;
					//////////////
					///rotate
					/////////
//					if (abs(deltaChi)>validError)
//					{
						p1 = portion.getAtomCoordinate (AAIndx,"E");   
						portion.rotateSideChain (AAIndx, portion.getAtomIndx (AAIndx, "Z"), p2, p1, deltaChi);
//					}
				}
			}
		}
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//return 0 in the case that no rotamer angles still un-evaluated
int resolveCollision_BBDep(Protein & portion, int * & lastChiUsed)
{

	if ((!BBDepRot.size()))
	{
		cout<<"============== in rotamer::resolveCollision(Protein &, int * &) ============"<<endl;
		cout<<"Please build rotVect or lastChiUsed first. . . ."<<endl;
		cout<<"============================================================================"<<endl;
		exit(1);
	}

	while (portion.doesCollide ())
	{
		int aa1indx		= portion.AAsCollide[0],		//get the indx of the first AA collides
			aa2indx		= portion.AAsCollide[1],		//get the indx of the second AA collides
			atom1Indx	= portion.atomsCollide [0],		//get the indx of the first atom collides
			atom2Indx	= portion.atomsCollide [1];		//get the indx of the second atom collides
		
		int rotIndx;		//the indx of the rotamer....which AA r u working on
		int bin1Indx, bin2Indx;		//the indx of the bin

		if (portion.AAs [aa1indx].atoms [atom1Indx].isSideChain) // is the first atom collides a side chain atom
		{
			if (portion.AAs [aa2indx].atoms [atom2Indx].isSideChain)	//is the second atomcollides a side chain atom
			{
			//Side chain with Side chain .... case 1

				rotIndx = chr2Num(portion.AAs [aa2indx].chr1);	//I gonna work on the second AA
				bin2Indx = getBinIndx(portion.AAs[aa2indx].angles.phi, portion.AAs[aa2indx].angles.psi);
				if ((lastChiUsed[aa2indx] != BBDepRot[rotIndx].bin[bin2Indx].rotamer.size() - 1) && 
					(BBDepRot[rotIndx].numOfChiAngles))  //not finish....u still have angles to be tested in the rotamer (exclude GLY and ALA)
				{
					applySideChainRotamer_BBDep(portion, aa2indx, rotIndx, bin2Indx, lastChiUsed);
				}
				else
				{
					//the rotamer of the second AA is finished but u still have collision....try the first AA now
					lastChiUsed[aa2indx] = -1;		//re-initiate the rotamer of the second AA
					applySideChainRotamer_BBDep(portion, aa2indx, rotIndx, bin2Indx, lastChiUsed);	// return the original side chain for the second AA
					
					rotIndx = chr2Num(portion.AAs [aa1indx].chr1); //get the rotamer indx of the first AA	
					bin1Indx = getBinIndx(portion.AAs[aa1indx].angles.phi, portion.AAs[aa1indx].angles.psi);
					
					if ((lastChiUsed[aa1indx] != BBDepRot[rotIndx].bin[bin1Indx].rotamer.size() - 1) &&
						(BBDepRot[rotIndx].numOfChiAngles))//if u still have angles to be tested
					{
						
						applySideChainRotamer_BBDep(portion, aa1indx, rotIndx, bin1Indx, lastChiUsed);

					}
					else
					{
						cout<<"============= In rotamer:: ResolveCollision(Protein &, int * &) =============="<<endl;
						cout<<"All possible angles for current pair have been used....could not resolve the collision.."<<endl;
						cout<<"=============================================================================="<<endl;
						return 0;		//success (no collision) = 1...fail = 0
					}
				}

			}
			else
			{
			//Side chain with BackBone ... case 2
				rotIndx = chr2Num(portion.AAs[aa1indx].chr1);	//get the indx of the first AA
				
				bin1Indx = getBinIndx(portion.AAs[aa1indx].angles.phi, portion.AAs[aa1indx].angles.psi);
				
				if ((lastChiUsed[aa1indx] != BBDepRot[rotIndx].bin[bin1Indx].rotamer.size() - 1) &&
					(BBDepRot[rotIndx].numOfChiAngles))
				{
					applySideChainRotamer_BBDep(portion, aa1indx, rotIndx, bin1Indx, lastChiUsed);
				}
				else
				{
					cout<<"============= In rotamer:: ResolveCollision(Protein &, int * &) =============="<<endl;
					cout<<"All possible angles for current pair have been used....could not resolve the collision.."<<endl;
					cout<<"=============================================================================="<<endl;
					return 0;   //success (no collision) = 1 ... fail = 0
				}

			}
		}
		else
		{
			if (portion.AAs[aa2indx].atoms [atom2Indx].isSideChain)
			{
			//BackBone with Side chain ... case 2
				rotIndx = chr2Num(portion.AAs[aa2indx].chr1); // get the rotamer indx of the second AA
				bin2Indx = getBinIndx(portion.AAs[aa2indx].angles.phi, portion.AAs[aa2indx].angles.psi);
				
				if ((lastChiUsed[aa2indx] != BBDepRot[rotIndx].bin[bin2Indx].rotamer.size() - 1) &&
					(BBDepRot[rotIndx].numOfChiAngles))
				{
					applySideChainRotamer_BBDep(portion, aa2indx, rotIndx, bin2Indx, lastChiUsed);
				}
				else
				{
					cout<<"============= In rotamer:: ResolveCollision(Protein &, int * &) =============="<<endl;
					cout<<"All possible angles for current pair have been used....could not resolve the collision....."<<endl;
					cout<<"=============================================================================="<<endl;
					return 0; //success (no collision) = 1 ... fail = 0
				}

			}
			else
			{
			//BackBone with BackBone....case 3
				cout<<".....The collision is between two backbone atoms....( "<<portion.AAs[aa1indx].chr3<<" "<<portion.AAs[aa1indx].num<<" ) with ( " \
					<<portion.AAs [aa2indx].chr3<<" "<<portion.AAs[aa2indx].num<<" )"<<endl;
				//for now
				return 0;
			}
		}

		//clear previous AA's and atoms were collide
		portion.AAsCollide .clear ();
		portion.atomsCollide .clear ();
	}

	return 1;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//build and initialize lastchiused array that will be used in resolvecollision function
void initializeLastChiUsed(int * &lastChiUsed, int size)
{

	lastChiUsed = new int [size];

	for (int i=0;i<size;i++)
		lastChiUsed[i] = -1 ;    // will represent the indx in rotVect
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// read side chains info from the directory that contains the AA templates and create a data structure to store each side chain
void createSideChainsArray()
{

	for (int i=0; i<STANDARD_AMINO_ACIDS; i++)
	{
		string AAFileName = SIDE_CHAINS_DIRECTORY + num2Chr3(i) + ".pdb";		//get the file name

		sideChains[i].read(AAFileName);			//store the Protein data stucture for the AA in the array
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// End Of IMPLEMENTATION ////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif