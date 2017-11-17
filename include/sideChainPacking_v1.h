#ifndef SIDECHAINPACKING_H
#define SIDECHAINPACKING_H

//#include <windows.h>
#include "protein.h"
#include "rotamer.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// List Of Data Structures //////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct PEDataStructure
{
	vector<double> PErotamer;		//pair energies for a particular rotamer
	int minIndx;					//the indx of the minimum pair energy entry b/w this rotamer and it's neighbors
	double minPE;
	int maxIndx;					//the indx of the maximal PE entry b/w this rotamer and it's neighbor's rotamer
	double maxPE;
	//initializer
	PEDataStructure() : minIndx(-1), minPE(0), maxIndx(-1), maxPE(0) {};
};
struct PEforAA
{
	int neighborIndx;					//the index of the neighbor in the list vector<AAinformation> nList
	vector<PEDataStructure> neighborPE;
};

struct PEforReducedAA
{
	int neighborIndx;
	int reducedAAIndx;
	vector<PEDataStructure> reducedAAindeces;

	//initializer
	PEforReducedAA() : neighborIndx(-1), reducedAAIndx(-1) {}
};

struct mergedN						//in the case that the reduced AA has just one neighbor
{
	vector<int> indx;			//indeces of neighbor rotamers that has been reduced (the rotamer the gives u the best energy)
	int nIndx;					//the indx of the neighbor
};
struct AAinformation
{
	int nNeighbors;
	double scR;
//	int indx;					//the indx of AA in the protein AAs
	int BBDepRotIndx;			//the indx of the AA in the rotamer library
	int binIndx;				//the indx of the bin used for this AA
	int chosenRotIndx;			//the indx of the rotamer to be chosen finally

	vector<int> neighbor;			//list of neighbors indeces
	vector<double> selfE;			//self-energies of rotamers
	Protein allRotamers;			//all rotamers translated and overlapped in the same position of the backbone
	Protein originalBB;				//original BackBone of the AA
	vector<PEforAA>	PE;				//Pair energy b/w all rotamers .. Matrix

	vector<mergedN> prevN;			//a vector that has the indeces of the rotamers that give best energies from a previous neighbor that has been reduced...the reduced AA has 1 neighbor
	vector<PEforReducedAA> prevNs;	//a vector that has best = ============		in the case that the reduced AA has two neighbors

	bool valid;						//flag used for processing to check if the residue is valid for more reduction or not

	//initializer
	AAinformation() : nNeighbors(0), scR(0), BBDepRotIndx(-1), binIndx(-1), chosenRotIndx(0), valid(false) {}
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// End Of Constants and Data Structures /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// List of Functions ////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void setSCRadius(Protein, vector<double> & scR);							//get the radius of side chain according to R3 paper (residue-rotamer-reduction)
//void predictSideChain(Protein & BBone);				//given a backbone...predict the conformation with GMEC (global minimal energy)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// End Of Function List /////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// IMPLEMENTATION ///////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void getFinalStructure(vector<AAinformation> & nList, Protein & finalStructure)
{
	for (int i=0; i<nList.size (); i++)
	{
		finalStructure.append(nList[i].originalBB ,0,0);			//copy backbone
		if (nList[i].selfE .size () >= 1)
			for (int j=0; j<nList[i].allRotamers.numOfAtoms(0); j++)
				finalStructure.AAs[i].atoms.push_back(nList[i].allRotamers .AAs[nList[i].chosenRotIndx].atoms[j]);
	}

	//finalStructure.writePDB("final.pdb",1,finalStructure.numOfAA());
}
/////////////////////////////////////////////////////////////////////////////////////////////////
void printRotamer(vector<AAinformation> nList, int AAindx)
{
	if ((AAindx>=0) && (AAindx< nList.size ()))
	{
		Protein tmp;
		tmp.append(nList[AAindx].originalBB ,0,0);
		for (int i=0; i<nList[AAindx].selfE .size (); i++)
			tmp.append(nList[AAindx].allRotamers , i, i);
		tmp.writePDB("rotamers.pdb",1,tmp.numOfAA());
	}
	else
		cout<<"error.....out of range indx given ("<<AAindx<<")"<<endl;
}
///////////////////////////////////////////////////////////////////////////////////////////////
void printPEinfo(vector<AAinformation> nList, string outFile)
{
	ofstream out;
	out.open(outFile.c_str());
	if (!out) {
		cout<<"================================ in Protein.writePDBFile() ===================="<<endl;
		cout<< "Unable to open " << outFile << endl;
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
		int i,j,k,l;		//counters
		for (i=0; i<nList.size ();i++)
		{
			out<<nList[i].originalBB .AAs[0].chr3<<" "<<nList[i].originalBB .AAs[0].num<<", numOfRot= "<<nList[i].selfE .size ()<<", neighbors (,";
			for (j=0; j<nList[i].neighbor .size (); j++)
				out<<nList[nList[i].neighbor [j]].originalBB .AAs[0].chr3<<" "<<nList[nList[i].neighbor [j]].originalBB .AAs[0].num<<",";
			out<<",)"<<endl;
		}
		out<<endl;

		for (i=0; i<nList.size (); i++)
		{
			out<<nList[i].originalBB .AAs[0].chr3<<" "<<nList[i].originalBB .AAs[0].num<<", scR= "<<nList[i].scR <<", binIndx= "<<nList[i].binIndx<<", neighbors (";
			for (j=0; j<nList[i].neighbor .size (); j++)
				out<<nList[nList[i].neighbor [j]].originalBB.AAs[0].chr3<<" "<<nList[nList[i].neighbor [j]].originalBB.AAs[0].num<<",";
			out<<")"<<endl;

			for (j=0; j<nList[i].PE .size (); j++)
			{
				out<<","<<j+1<<", With "<<nList[nList[i].PE [j].neighborIndx].originalBB.AAs[0].chr3<<" "<<nList[nList[i].PE [j].neighborIndx].originalBB.AAs[0].num<<endl;
				for (k=0; k<nList[i].PE [j].neighborPE .size (); k++)
				{
					out<<",,, Rotamer# "<<k+1<<", selfE= "<<nList[i].selfE [k]<<", IndxMinRot= "<<nList[i].PE [j].neighborPE [k].minIndx <<", minPE= "<<nList[i].PE [j].neighborPE [k].minPE
						<<", IndxMaxRot = "<<nList[i].PE [j].neighborPE [k].maxIndx <<", maxPE= "<<nList[i].PE [j].neighborPE [k].maxPE <<endl;
					for (l=0; l<nList[i].PE [j].neighborPE [k].PErotamer .size (); l++)
					{
						out<<",,,,,,"<<nList[nList[i].PE [j].neighborIndx].originalBB .AAs[0].chr3<<" "<<nList[nList[i].PE [j].neighborIndx].originalBB .AAs[0].num<<", Rotamer # "<<l+1<<", PE= "<<nList[i].PE [j].neighborPE [k].PErotamer [l]<<endl;
					}
				}
			}
		}

	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//scR vector that contains the radius of sidechain as in R3 (same way they compute it)
void setSCRadius(Protein BBone, vector<AAinformation> & naibors)
{
	int *lastRotIndx;
	int size = BBone.numOfAA();
	int rotIndx;
	double maxDist, dist;
	Coordinate CbCoord;

	initializeLastChiUsed(lastRotIndx, size);

	for (int i=0 ; i<BBone.numOfAA(); i++)
	{
		AAinformation tmpN;

		maxDist = 0;

		rotIndx = chr2Num(BBone.AAs[i].chr1);	//which AA we r working on...the indx in BBDepRot

		tmpN.BBDepRotIndx = rotIndx;					//get the indx of this AA in BBDepRot

		/*****
				copy backbone atoms with CB
		*****/
		vector<Atom> tmpAtoms;
		//AminoAcid tmpAA;
		for (int k=0; k< BBone.numOfAtoms(i); k++)
			if ((!BBone.AAs[i].atoms[k].isSideChain) || (BBone.AAs[i].atoms[k].name == " CB "))
				tmpAtoms.push_back(BBone.AAs[i].atoms[k]);
		//tmpAA.atoms = tmpAtoms;
		tmpN.originalBB .AAs.push_back(BBone.AAs[i]);
		tmpN.originalBB .AAs[0].atoms = tmpAtoms;


		//ALA and GLY has no rotamers
		if (BBDepRot[rotIndx].numOfChiAngles >= 1)
		{
			tmpN.binIndx = getBinIndx(BBone.AAs[i].angles.phi, BBone.AAs[i].angles.psi);			//the indx of the bin in BBDepRot

			CbCoord = BBone.getAtomCoordinate(i, "CB");

			for (int j=0; j<BBDepRot[rotIndx].bin[tmpN.binIndx].rotamer.size(); j++)
			{
				applySideChainRotamer_BBDep(BBone, i, rotIndx, tmpN.binIndx, lastRotIndx);

				/*****
						Add Current Rotamer
				*****/
				vector<Atom> tmpAtoms2;
				for (int cntr=0; cntr< BBone.numOfAtoms(i); cntr++)
					if ((BBone.AAs[i].atoms[cntr].isSideChain) && (BBone.AAs[i].atoms[cntr].name != " CB "))
						tmpAtoms2.push_back(BBone.AAs[i].atoms[cntr]);
				tmpN.allRotamers .AAs.push_back(BBone.AAs[i]);
				tmpN.allRotamers .AAs[j].atoms = tmpAtoms2;

				tmpN.chosenRotIndx = -1;					//no rot is chosen yet

				tmpN.selfE .push_back (0);					//push initial self energy value for this rotamer

				for (int k=0; k<BBone.numOfAtoms(i); k++)
				{
					if ((BBone.AAs[i].atoms[k].isSideChain) && (BBone.AAs[i].atoms[k].name != " CB "))
					{
						dist = getDistance(CbCoord, BBone.AAs[i].atoms[k].coord);
						if (dist > maxDist)
							maxDist = dist;
					}
				}
			}
			tmpN.valid = true;
		}

		tmpN.scR = maxDist;

		naibors.push_back (tmpN);
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double getElib(int BBDepIndx, int binIndx, int rotIndx)
{
	int constantK = 3;

	return (-1 *  constantK * log (BBDepRot[BBDepIndx].bin[binIndx].rotamer[rotIndx].prob/ BBDepRot[BBDepIndx].bin[binIndx].rotamer[0].prob));
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void add_BB_VDW_Energy(vector<AAinformation> & nList, int targetAAindx, int curAAindx)
{
	int k, m, l;
	double dist;

	for (k=0; k< nList[targetAAindx].selfE.size(); k++)		//all rotamers of target AA
	{
		//cout<<"  ROTAMER # "<<k+1<<" ===================================  "<<nList[targetAAindx].originalBB .AAs[0].chr3<<" "<<nList[targetAAindx].originalBB .AAs[0].num<<endl;
		for (m=0; m< nList[targetAAindx].allRotamers.numOfAtoms(k);m++)		//all atoms in the target rotamer
		{
			char targetAtomType = nList[targetAAindx].allRotamers.AAs[k].atoms[m].type;
			Coordinate targetAtomCoord = nList[targetAAindx].allRotamers.AAs[k].atoms[m].coord;
			//cout<<"   working on atom "<<nList[targetAAindx].allRotamers.AAs[k].atoms[m].name<<" with:"<<endl;
			for (l=0; l<nList[curAAindx].originalBB.numOfAtoms(0); l++)			//neighbor AA...
			{
				dist = getDistance(targetAtomCoord, nList[curAAindx].originalBB.AAs[0].atoms[l].coord);
				nList[targetAAindx].selfE[k] += nList[targetAAindx].allRotamers.VDW(getRadius(targetAtomType, 'V'), getRadius(nList[curAAindx].originalBB.AAs[0].atoms[l].type, 'V'), dist);

				//cout<<nList[curAAindx].originalBB.AAs[0].chr3<<" "<<nList[curAAindx].originalBB .AAs[0].num<<" "	\
					<<nList[curAAindx].originalBB .AAs[0].atoms[l].name<<" Dist= "<<dist	\
					<<" VDW= "<<nList[targetAAindx].allRotamers.VDW(getRadius(targetAtomType, 'V'), getRadius(nList[curAAindx].originalBB.AAs[0].atoms[l].type, 'V'), dist)<<endl;
			}
		}
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void setSelfEnergy(vector<AAinformation> & nList)
{
	int i, j;	//counters
	Coordinate CbCoord;
	double dist;

	for (i=0; i<nList.size(); i++)
	{
		//cout<<"working on "<<nList[i].originalBB.AAs[0].chr3<<" "<<nList[i].originalBB .AAs[0].num<<" "<<nList[i].selfE .size()<<" rotamers"<<endl;;
		//ALA and GLY has no rotamers
		if (nList[i].selfE .size () > 0)
		{
			//get Elib energy for each rotamer
			for (j=0; j<nList[i].selfE .size (); j++)
				nList[i].selfE [j] = getElib(nList[i].BBDepRotIndx, nList[i].binIndx, j);

			//get the coordinate of Cb.
			CbCoord = nList[i].originalBB.getAtomCoordinate(0, "CB");

			for (j=0; j<nList.size(); j++)
			{
				if ((j!=i))		//not same AA
				{
					//find the distance b/w Cb and Ca of neighbor AA's
					dist = getDistance(CbCoord, nList[j].originalBB.getAtomCoordinate(0, "CA"));
					if (dist <= (nList[i].scR + 5.0))
					{
						//compute and add VDW energy b/w side chain atoms and backbone atoms of all rotamers of target AA with neighbor AA
						//cout<<"  neighbor "<<nList[j].originalBB .AAs[0].chr3<<" "<<nList[j].originalBB .AAs[0].num<<endl;
						add_BB_VDW_Energy(nList, i, j);
						//getchar();
					}
				}
			}
		}
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void buildPEArray(vector<AAinformation> & nList, int targetAAindx, int curAAindx)
{
	int k,l, m, n;		//counters

	double PE, dist;



	/*****
			build the PE matrix and initiate it to zeros between targetAA (rows) and curAA (cols)
	*****/
	PEforAA tmpPERotamersTarget;
	PEDataStructure tmpRecTarget;
	tmpPERotamersTarget.neighborIndx = curAAindx;			//set the neighbor

	for (k=0; k<nList[targetAAindx].selfE .size (); k++)
	{
		vector<double> tmpPE (nList[curAAindx].selfE .size (), 0);
		tmpRecTarget.PErotamer = tmpPE;
		tmpPERotamersTarget.neighborPE .push_back (tmpRecTarget);
	}

	/*****
			build the PE matrix and initiate it to zeros between curAA (rows) and targetAA (cols)
	*****/
	PEforAA tmpPERotamersCur;
	PEDataStructure tmpRecCur;
	tmpPERotamersCur.neighborIndx = targetAAindx;			//set the neighbor

	for (k=0; k<nList[curAAindx].selfE .size (); k++)
	{
		vector<double> tmpPE(nList[targetAAindx].selfE .size (), 0);
		tmpRecCur.PErotamer = tmpPE;
		tmpPERotamersCur.neighborPE .push_back (tmpRecCur);
	}
	/*****
				find PE b/w two AA's and fill matrices
	*****/
	for (k=0; k<nList[targetAAindx].selfE.size(); k++)						//for all rotamers of the target AA
	{

		int minIndx = -1, maxIndx = -1;
		double minPE = 1000, maxPE = -1000;

		//cout<<"Rotamer # "<<k+1<<"  ==============="<<endl;
		for (l=0; l<nList[curAAindx].selfE.size(); l++)						//for all rotamers of the neighbor AA
		{
			PE = 0;				//initialize PE for each rotamer
			for (m=0; m< nList[targetAAindx].allRotamers.numOfAtoms(k); m++)		//for all atoms of the target rotamer
			{
				//cout<<"  "<<nList[targetAAindx].allRotamers .AAs[k].atoms[m].name<<" with ";
				char targetAtomType = nList[targetAAindx].allRotamers.AAs[k].atoms[m].type;
				Coordinate targetAtomCoord = nList[targetAAindx].allRotamers.AAs[k].atoms[m].coord;

				for (n=0; n<nList[curAAindx].allRotamers.numOfAtoms(l); n++)	//for all atoms of the neighbor rotamer
				{
					dist = getDistance(targetAtomCoord, nList[curAAindx].allRotamers.AAs[l].atoms[n].coord);
					PE += nList[targetAAindx].allRotamers.VDW(getRadius(targetAtomType, 'V'),
															  getRadius(nList[curAAindx].allRotamers.AAs[l].atoms[n].type, 'V'),
															  dist);
					//cout<<nList[curAAindx].allRotamers .AAs[l].atoms[n].name<<" dist= "<<dist \
						<<" VDW= "<<nList[targetAAindx].allRotamers.VDW(getRadius(targetAtomType, 'V'), getRadius(nList[curAAindx].allRotamers.AAs[l].atoms[n].type, 'V'), dist)<<endl;
				}

			}

			//cout<<"PE = "<<PE<<endl;
			if (PE < minPE )
			{
				minIndx = l;			//the indx of the neighbor has minimum PE with target rotamer
				minPE = PE;
			}
			if (PE > maxPE )
			{
				maxIndx = l;
				maxPE = PE;
			}

			tmpPERotamersTarget.neighborPE [k].PErotamer [l] = PE;
			tmpPERotamersCur.neighborPE [l].PErotamer [k] = PE;
		}
		//cout<<"Done..."<<endl;
		//getchar();
		tmpPERotamersTarget.neighborPE [k].minIndx = minIndx;
		tmpPERotamersTarget.neighborPE [k].minPE = minPE;
		tmpPERotamersTarget.neighborPE [k].maxIndx = maxIndx;
		tmpPERotamersTarget.neighborPE [k].maxPE = maxPE;
	}

	/*****
			set min and max for curAA
	*****/
	for (k=0; k<tmpPERotamersCur.neighborPE .size (); k++)
	{
		int minIndx = -1, maxIndx = -1;
		double minPE = 1000, maxPE = -1000;
		double PE;
		for (l=0; l<tmpPERotamersCur.neighborPE [k].PErotamer .size (); l++)
		{
			PE = tmpPERotamersCur.neighborPE [k].PErotamer [l];
			if ( PE < minPE )
			{
				minIndx = l;			//the indx of the neighbor has minimum PE with target rotamer

				minPE = PE;
			}
			if ( PE > maxPE )
			{
				maxIndx = l;
				maxPE = PE;
			}
		}
		tmpPERotamersCur.neighborPE [k].maxIndx = maxIndx;
		tmpPERotamersCur.neighborPE [k].maxPE = maxPE;
		tmpPERotamersCur.neighborPE [k].minIndx = minIndx;
		tmpPERotamersCur.neighborPE [k].minPE = minPE;
	}

	nList[targetAAindx].PE .push_back (tmpPERotamersTarget);
	nList[curAAindx].PE .push_back (tmpPERotamersCur);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void setPairEnergy(vector<AAinformation> & naibors)
{
	int i, j;
	Coordinate CbCoord;
	double dist;

	/*****
			find neigbors of an AA and then build the PairEnergy (PE) data structure
	*****/

	for (i=0; i<naibors.size(); i++)
	{
		//ALA and GLY has No Rotamers
		if ( (naibors[i].selfE.size () > 0))
		{
			//get the coordinate of Cb..
			CbCoord = naibors[i].originalBB.getAtomCoordinate(0, "CB");

			for (j=i+1; j<naibors.size(); j++)
			{
				if (naibors[j].selfE.size () > 0)		//not ALA nor GLY
				{
					dist = getDistance(CbCoord, naibors[j].originalBB.getAtomCoordinate(0, "CB"));

					if (dist <= (naibors[i].scR + naibors[j].scR + 4.0))					//it's a neighbor AA
					{
						//push back the indx of the neighbor
						naibors[i].neighbor.push_back(j);
						naibors[i].nNeighbors ++;				//increment number of neighbors of i by 1

						naibors[j].neighbor .push_back (i);
						naibors[j].nNeighbors ++;				//increment number of neighbors of j by 1

						/*****
								actual build of PE matrices for both neighbors
						*****/
						buildPEArray(naibors, i, j);
					}

				}
			}
		}
	}

/*
	//print the list
	for (i=0; i<naibors.size(); i++)
	{
		if (naibors[i].selfE .size () >0)
		{
			cout<<naibors[i].allRotamers.AAs[0].chr3<<" "<<i<<" nNeighbors= "<<naibors[i].nNeighbors <<" scR= "<<naibors[i].scR	\
				<<" BBDIndx= "<<naibors[i].BBDepRotIndx <<" binIndx= "<<naibors[i].binIndx <<" nRot= "<<naibors[i].selfE.size()<<endl;
			for (j=0; j<naibors[i].neighbor .size (); j++)
				cout<<"   "<<" "<<naibors[i].neighbor [j]<<endl;
		}
	}
*/
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int getPEindx(vector<AAinformation> & nList, int curAAindx, int neighborIndx)
{
	for (int i=0; i<nList[curAAindx].PE .size (); i++)
		if (nList[curAAindx].PE [i].neighborIndx == neighborIndx)
			return i;

	return -1;  //not found
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void filterNeighbors(vector<AAinformation> & nList, double epsilon)
{
	int i=0, j=0, k=0, l=0; //counters
	vector<int> PEindeces;		//indeces of Neighbors to be deleted from PE vector b/s after filtering they are not neighbors any more

	for (i=0; i<nList.size (); i++)		//clear previous neighobrs list to be updated
	{
		nList[i].neighbor .clear ();
		nList[i].nNeighbors = 0;
	}

	for (i=0; i<nList.size (); i++)
	{
		for (j=0; j<nList[i].PE .size () ; j++)
		{
			//while ((k<nList[i].PE [j].neighborPE .size ()) && (!neighbors))
			double max = -1, min= 100000,
					maxPE, minPE;
			/*****
					find max PE value and min PE value in all combinations of rotamers b/w these 2 AA's
			*****/
			for (k=0; k<nList[i].PE [j].neighborPE .size (); k++)
			{
				maxPE = nList[i].PE[j].neighborPE[k].maxPE;					//the maximal PE b/w this rotamer and all other rotamers from other AA
				minPE = nList[i].PE[j].neighborPE[k].minPE;					//the minimal PE b/w this rotamer and all other rotamers from other AA

				if (max < maxPE)
					max = maxPE;

				if (min > minPE)
					min = minPE;
			}

			if (max - min <= epsilon)
				PEindeces.push_back (j);			//the indeces of ex-neighbors to be deleted from the list (PE matrices list)
		}

		/*****
				delete ex-neighbors
		*****/
		int numOfDeletedEntries = 0;
		for (j=0; j<PEindeces.size (); j++)
		{
			int nIndx = nList[i].PE [PEindeces[j] - numOfDeletedEntries].neighborIndx;				//the index of the neighbor to be deleted
			nList[i].PE .erase (nList[i].PE .begin () + PEindeces[j] - numOfDeletedEntries++);		//delete from PairEnergy matrix

			int curIndx = getPEindx(nList, nIndx, i);		//the indx of current AA in the PE list of it's neighbor

			//cout<<nList[i].originalBB .AAs[0].num<<" with "<<nList[nIndx].originalBB .AAs[0].num<<" indx erased= "<<PEindeces[j] - (numOfDeletedEntries - 1)<<" "<<curIndx<<endl;
			nList[nIndx].PE .erase (nList[nIndx].PE .begin () + curIndx);				//delete current AA from it's neighbor PE list
		}

		PEindeces.clear ();				//clear previous indeces for the previous AA

		for (j=0; j<nList[i].PE .size (); j++)				//update neighbors
		{
			nList[i].neighbor .push_back (nList[i].PE [j].neighborIndx );			//for current AA
			nList[i].nNeighbors ++;
		}

	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void solveReducedAAs(vector<AAinformation> & nList, int AAindx)
{
	int i;			//counter

	//cout<<"  Number of related reduced AA for this AA is "<<nList[AAindx].prevN .size ()<<endl;
	for (i = nList[AAindx].prevN .size () - 1; i>=0; i--)
	{
		int nIndx = nList[AAindx].prevN [i].nIndx;					//neighbor indx
		//cout<<"  neighbor = "<<nList[nIndx].originalBB .AAs[0].num<<" ";
		if (nList[nIndx].chosenRotIndx < 0)
		{
			nList[nIndx].chosenRotIndx = nList[AAindx].prevN [i].indx [nList[AAindx].chosenRotIndx];
			//cout<<"solved "<<nList[nIndx].chosenRotIndx<<endl;
		}
		//else
			//cout<<"not solved."<<endl;
	}

	//cout<<" Number of related with two neighbors.. "<<nList[AAindx].prevNs .size ()<<endl;
	for (i = nList[AAindx].prevNs .size () - 1; i>= 0; i--)
	{
		//cout<<" i= "<<i<<" neighborPE size = "<<nList[AAindx].prevNs [i].neighborPE .size ();
		int nIndx = nList[AAindx].prevNs [i].neighborIndx;			//the indx of the neighbor	is stored in this cell
		//cout<<" neighbor is  "<<nList[nIndx].originalBB .AAs[0].num<<"  rotChosenIndx= "<<nList[nIndx].chosenRotIndx<<" ";
		if (nList[nIndx].chosenRotIndx >= 0)
		{
			//cout<<"its rotamer is already chosen.. ";
			if (nList[nList[AAindx].prevNs [i].reducedAAIndx].chosenRotIndx < 0)
			{
				nList[nList[AAindx].prevNs [i].reducedAAIndx].chosenRotIndx = nList[AAindx].prevNs [i].reducedAAindeces [nList[AAindx].chosenRotIndx].PErotamer [nList[nIndx].chosenRotIndx ];
				//cout<<" reducedAA = "<<nList[nList[AAindx].prevNs [i].reducedAAIndx].originalBB .AAs[0].num<<" solved "<<nList[nList[AAindx].prevNs [i].reducedAAIndx].chosenRotIndx<<endl;
			}
			//else
				//cout<<" reducedAA == "<<nList[nList[AAindx].prevNs [i].reducedAAIndx].originalBB .AAs[0].num<<" previously solved."<<endl;
			//cout<<" reduced AA indx "<<nList[AAindx].prevNs [i].neighborIndx<<endl;
		}
		//else
			//cout<<" reducedAA === "<<nList[nList[AAindx].prevNs [i].reducedAAIndx].originalBB .AAs[0].num<<" not solved "<<endl;
	}

	//cout<<endl;


}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//reduced AA is the vector contains Amino Acids that have been reduced from the graph but their candidate rotamer not yet chosen
void reduceResidues(vector<AAinformation> & nList, vector<int> & reducedAAs)
{
	int i, j, k, l;		//counters
	bool cont;
	do
	{
		cont = false;
		for (i=0; i<nList.size (); i++)
		{
			if (nList[i].valid )
			{
				if ((nList[i].nNeighbors == 0) && (nList[i].selfE .size () == 1))
				{
					//cout<<nList[i].originalBB .AAs[0].chr3<<" "<<nList[i].originalBB .AAs[0].num<<" has no neighbors but one Rotamer 00.."<<endl;

					nList[i].valid = false;
					nList[i].chosenRotIndx = 0;

					reducedAAs.push_back (i);

					solveReducedAAs(nList, i);

					continue;
				}

				if ((nList[i].nNeighbors == 0) && (nList[i].selfE.size() > 1))			//has no neighbors and has more than one rotamer
				{
					//cout<<nList[i].originalBB .AAs[0].chr3<<" "<<nList[i].originalBB .AAs[0].num<<" has no neighbors but multiple rotamers.."<<endl;
					//keep the rotamer with minimal self energy and discard other rotamers
					double minSelfE = nList[i].selfE [0];		//initially...first rotamer is the best one
					int minIndx = 0;							//the indx of the rotamer who is the minimal according selfE
					for (j=1; j<nList[i].selfE .size (); j++)
						if (nList[i].selfE [j] < minSelfE)
						{
							minSelfE = nList[i].selfE [j];
							minIndx = j;
						}
					nList[i].chosenRotIndx = minIndx;			//THE INDX OF THE ROTAMER WILL BE CHOSEN FINALLY
					nList[i].valid = false;											//we r done with this AA
					cont = true;

					reducedAAs.push_back (i);

					solveReducedAAs(nList, i);			//for previously reduced AA's ... now u can decide which rotamer to select according to the rotamer chosen here

					continue;
				}

				if ((nList[i].nNeighbors >= 1)	&& (nList[i].selfE .size () == 1))	//with only one rotamer but has neighbors
				{
					//cout<<nList[i].originalBB .AAs[0].chr3<<" "<<nList[i].originalBB .AAs[0].num<<" has one rotamer and multiple neighbors."<<endl;

					for (j=0; j<nList[i].PE .size (); j++)
					{
						int nIndx = nList[i].PE [j].neighborIndx;			//the indx of the neighbor

						for (k=0; k<nList[nIndx].selfE .size (); k++)
							nList[nIndx].selfE [k] += nList[i].PE [j].neighborPE [0].PErotamer [k];		//add PE to self energy of the neighbor's rotamer

						/*****
								Delete PE matrices and the indeces of neighbors from neighbor AA lists
						*****/
						int nPEindx = getPEindx(nList, nIndx, i);
						nList[nIndx].PE.erase (nList[nIndx].PE .begin () + nPEindx);
						nList[nIndx].neighbor .erase (nList[nIndx].neighbor .begin () + nPEindx);
						nList[nIndx].nNeighbors --;
					}

					//Delete all neighbors for current AA
					nList[i].neighbor .clear ();
					nList[i].nNeighbors = 0;
					nList[i].PE .clear ();
					nList[i].chosenRotIndx = 0;
					nList[i].valid = false;
					cont = true;

					reducedAAs.push_back (i);

					solveReducedAAs(nList, i);			//for previously reduced AA's ... now u can decide which rotamer to select according to the rotamer chosen here

					continue;
				}

				if (nList[i].nNeighbors == 1)									//has one neighbor
				{
					//cout<<endl<<nList[i].originalBB .AAs[0].chr3<<" "<<nList[i].originalBB .AAs[0].num<<" has one neighbor..."<<nList[nList[i].neighbor [0]].originalBB .AAs[0].num<<endl;


					int nIndx = nList[i].neighbor [0],				//neighbor indx
						PEindx = getPEindx(nList, nIndx, i);		//the indx of PE matrix in the neighbor AA lists

					/*****
						for each rotamer of neighbor AA, find which rotamer of this AA that contributes with less energy and then
						save it's indx in mergedN for neighbor energy and update the self energy of neighbor's rotamers
					*****/

					mergedN tmpRec;

					tmpRec.nIndx = i;				//the indx of the AA (current) to be reduced

					double curContributionE;
					for (j=0; j<nList[nIndx].PE [PEindx].neighborPE .size (); j++)
					{
						int minIndx = 0;
						double minE = nList[nIndx].PE [PEindx].neighborPE [j].PErotamer [0] + nList[i].selfE [0];
						for (k=1; k<nList[i].selfE .size (); k++)
						{
							curContributionE = nList[nIndx].PE [PEindx].neighborPE [j].PErotamer [k] + nList[i].selfE [k];
							//cout<<" with "<<k<<" self= "<<nList[i].selfE [k]<<" PE= "<<nList[nIndx].PE [PEindx].neighborPE [j].PErotamer [k]<<endl;
							if (curContributionE < minE)				//find the rotamer contributes min energy with this rotamer
							{
								minE = curContributionE;
								minIndx = k;
							}
						}
						//push the indx of the rotamer from currunt AA that contributes with min energy and update self energy of
						// the rotamer of neighbor AA

						tmpRec.indx .push_back (minIndx);
						nList[nIndx].selfE [j] += minE;					//update self energy
						//cout<<endl<<"rot "<<j<<" minE = "<<minE<<" minIndx= "<<minIndx<<" currentSelfE = "<<nList[nIndx].selfE [j]<<endl;
					}
					//getchar();
					// save the indeces in neighbor preN vector
					nList[nIndx].prevN .push_back (tmpRec);

					//clear PE matrices and remove AA's indeces from neighbor lists
					nList[nIndx].PE.erase (nList[nIndx].PE .begin () + PEindx);
					nList[nIndx].neighbor .erase (nList[nIndx].neighbor .begin () + PEindx);
					nList[nIndx].nNeighbors --;					//increment number of neighbors by 1

					nList[i].PE .clear ();
					nList[i].neighbor .clear ();
					nList[i].nNeighbors = 0;
					nList[i].valid = false;

					reducedAAs.push_back (i);				//the AA is reduced but it's rotamer not chosen yet

					cont = true;

					continue;

				}

				if (nList[i].nNeighbors == 2)								//has 2 neighbors
				{
					//cout<<endl<<nList[i].originalBB .AAs[0].chr3<<" "<<nList[i].originalBB .AAs[0].num<<"     has two neighbors.."<<nList[nList[i].neighbor [0]].originalBB .AAs[0].num<<"  and  "<<nList[nList[i].neighbor [1]].originalBB .AAs[0].num<<endl;

					int nIndx1 = nList[i].neighbor [0],				//indx of neighbor 1
						nIndx2 = nList[i].neighbor [1],				//indx of neighbor 2
						PEindx1 = getPEindx(nList, nIndx1, i),		//indx of PE matrix in the first neighbor's lists
						PEindx2 = getPEindx(nList, nIndx2, i),		//indx of PE matrix in the second neigbor's lists
						PEnew1 = getPEindx(nList, nIndx1, nIndx2),	//the indx of PE matrix in neighbor 1 for neighbor 2 (b/w them)
						PEnew2 = getPEindx(nList, nIndx2, nIndx1);

					//check if neighbor 1 and neighbor 2 are already neighbors
					if (PEnew1 == -1)
					{
						//they are not neighbors....then create new PE matrix and store it in nList for both AA's

						//create PE for neighbor 1
						PEforAA tmpPERotamers;
						PEDataStructure tmpRec;
						tmpPERotamers.neighborIndx = nIndx2;

						tmpRec.maxIndx = tmpRec.maxPE = tmpRec.minIndx = tmpRec.minPE = 0;		//initiate max and min indeces

						for (j=0; j<nList[nIndx1].selfE .size (); j++)
						{
							vector<double> tmpPE(nList[nIndx2].selfE .size (), 0);

							tmpRec.PErotamer = tmpPE;

							tmpPERotamers.neighborPE .push_back (tmpRec);
						}

						PEnew1 = nList[nIndx1].PE.size ();			//get the indx of new PE matrix
						nList[nIndx1].PE.push_back  (tmpPERotamers);	//push the new PE matrix

						//create PE for neighbor 2
						tmpPERotamers.neighborPE .clear ();
						tmpPERotamers.neighborIndx = nIndx1;
						for (j=0; j<nList[nIndx2].selfE .size (); j++)
						{
							vector<double> tmpPE(nList[nIndx1].selfE .size (), 0);
							tmpRec.PErotamer = tmpPE;

							tmpPERotamers.neighborPE .push_back (tmpRec);
						}


						PEnew2 = nList[nIndx2].PE .size ();
						nList[nIndx2].PE .push_back (tmpPERotamers);

						//add both AA's to neighbor list of each other
						nList[nIndx1].neighbor .push_back (nIndx2);
						nList[nIndx1].nNeighbors ++;
						nList[nIndx2].neighbor .push_back (nIndx1);
						nList[nIndx2].nNeighbors ++;
					}

					//create an adjacent matrix that contains the indeces of the current AA that contributes minimal E. with
					// the combination of the two neighbors
					int prevIndx1 = nList[nIndx1].prevNs .size (),			//get the indx of prevNs where the adjacent matrix will be stored
						prevIndx2  = nList[nIndx2].prevNs.size ();

					//intialize reducedAA indeces matrix for both neighbors
					PEforReducedAA tmpMtrx;

					tmpMtrx.neighborIndx = nIndx2;
					tmpMtrx.reducedAAIndx = i;
					tmpMtrx.reducedAAindeces = nList[nIndx1].PE [PEnew1].neighborPE;

					nList[nIndx1].prevNs .push_back (tmpMtrx);		//push an empty matrix which contains the indeces of rotamers from current AA

					//and neighbor 2
					tmpMtrx.neighborIndx = nIndx1;
					tmpMtrx.reducedAAindeces = nList[nIndx2].PE [PEnew2].neighborPE ;

					nList[nIndx2].prevNs .push_back (tmpMtrx);		//that contribute with minimal energy to GMEC


					/*****
							store rotamers from current AA that contributes minimal E. and update PE b/w the two neighbors
							check where the current AA in respect to the neighbors first b/s the position (before, after, or in between)
							will let u know where are rotamers in the row or the col. of PE Matrix..b/s always PE matrix is stored with the
							AA with small Indx..and then the rotamers of this AA will be in the row of this matrix and the other rotamers (for the other AA)
							will be in col.
					*****/
					double curContributionE;
					double minE;
					int minIndx;

					for (j=0; j<nList[nIndx1].selfE .size (); j++)
					{
						for (k=0; k<nList[nIndx2].selfE .size (); k++)
						{
							//cout<<"rot "<<j+1<<" "<<nList[smallestIndx].originalBB .AAs[0].num<<" with rot "<<k+1<<" "<<nList[largestIndx].originalBB .AAs[0].num<<endl;

							minIndx = 0;
							minE = nList[i].selfE [0] + nList[i].PE [0].neighborPE [0].PErotamer [j] + nList[i].PE [1].neighborPE [0].PErotamer [k];
							for (l=1; l<nList[i].selfE .size (); l++)
							{
								curContributionE = nList[i].selfE [l] + nList[i].PE [0].neighborPE [l].PErotamer [j] + nList[i].PE [1].neighborPE [l].PErotamer [k];
								//cout<<"    rot "<<l+1<<" "<<nList[i].originalBB .AAs[0].num<<"  contributionE = " <<curContributionE<<endl;

								if (curContributionE < minE)
								{
									minE = curContributionE;
									minIndx = l;
								}
							}

							//update energy term and the indx of the rotamer that gives minimal energy for neighbor 1
							nList[nIndx1].PE [PEnew1].neighborPE [j].PErotamer [k] += minE;
							nList[nIndx1].prevNs[prevIndx1].reducedAAindeces [j].PErotamer [k] = minIndx;

							//and for neighbor 2
							nList[nIndx2].PE [PEnew2].neighborPE [k].PErotamer [j] += minE;
							nList[nIndx2].prevNs [prevIndx2].reducedAAindeces [k].PErotamer [j] = minIndx;
						}
					}


					//delete neighbors from neighbor lists
					nList[i].neighbor .clear ();
					nList[i].nNeighbors = 0;
					nList[i].PE .clear ();
					nList[i].valid = false;
					reducedAAs.push_back (i);			//the AA is reduced but it's rotamer is not chosen yet
					cont = true;

					nList[nIndx1].neighbor .erase (nList[nIndx1].neighbor .begin () + PEindx1);
					nList[nIndx1].PE .erase (nList[nIndx1].PE .begin () + PEindx1);
					nList[nIndx1].nNeighbors--;

					nList[nIndx2].PE .erase (nList[nIndx2].PE .begin () + PEindx2);
					nList[nIndx2].neighbor .erase (nList[nIndx2].neighbor .begin () + PEindx2);
					nList[nIndx2].nNeighbors--;

					continue;
				}
			}
		}
	} while (cont);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void eraseRotamer(vector<AAinformation> & nList, int AAindx, int rotIndx)
{
	int k, l;

	nList[AAindx].selfE .erase (nList[AAindx].selfE .begin () + rotIndx);	//delete from the vector contains self energies first

	nList[AAindx].allRotamers .deleteAA(rotIndx);			//delete from rotamers structures
	//cout<<"Done with SelfEnergies.."<<endl;

	for (k=0; k<nList[AAindx].PE .size (); k++)							//deleted from PE matrices of current AA
		nList[AAindx].PE [k].neighborPE .erase (nList[AAindx].PE [k].neighborPE .begin () + rotIndx);

	for (k=0; k<nList[AAindx].prevNs .size (); k++)
		nList[AAindx].prevNs [k].reducedAAindeces .erase (nList[AAindx].prevNs [k].reducedAAindeces .begin () + rotIndx);

	for (k=0; k<nList[AAindx].prevN .size (); k++)
		nList[AAindx].prevN[k].indx .erase (nList[AAindx].prevN[k].indx .begin () + rotIndx);
	//cout<<"Done with self PE matrices.."<<endl;


	for (k=0; k<nList[AAindx].nNeighbors ; k++)					//delete from it's neighbors PE matrices
	{
		int nIndx = nList[AAindx].neighbor [k];			//the indx of the neighbor
		int PEindx = getPEindx(nList, nIndx, AAindx);
		for (l=0; l<nList[nIndx].PE [PEindx].neighborPE .size () ;l++)
			nList[nIndx].PE [PEindx].neighborPE [l].PErotamer .erase (nList[nIndx].PE [PEindx].neighborPE [l].PErotamer.begin () + rotIndx);

		//cout<<"Done with one neighbor PE matrices.."<<endl;

		for (l=0; l<nList[nIndx].prevN .size (); l++)
			if (nList[nIndx].prevN[l].nIndx == AAindx)
				nList[nIndx].prevN[l].indx .erase (nList[nIndx].prevN[l].indx .begin () + rotIndx);
		//cout<<"Done with prevN..for one neighbor.."<<endl;

		for (l=0; l<nList[nIndx].prevNs .size () ;l++)
			if (nList[nIndx].prevNs [l].neighborIndx  == AAindx)
				for (int m=0; m<nList[nIndx].prevNs [l].reducedAAindeces  .size (); m++)
					nList[nIndx].prevNs [l].reducedAAindeces [m].PErotamer.erase (nList[nIndx].prevNs [l].reducedAAindeces[m].PErotamer.begin () + rotIndx);
		//cout<<"Done with prevNs for one neighbor.."<<endl;
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void reduceRotamers(vector<AAinformation> & nList, bool & hasReducedRot)
{
	int i, j, k, l;		//counters
	double minDiff, PEdiff,
		totalContributionE;				//Cir - Cis + segma(min(PE - PE))	from R3 paper	proposition 4.1
	vector<int> rotsToBeDeleted;		//indeces of rotamers to be deleted

	hasReducedRot = false;

	for (i=0; i<nList.size (); i++)
	{
		if (nList[i].valid )
		{
			rotsToBeDeleted.clear ();
			//cout<<" before k...";
			for (k=0; k<nList[i].selfE.size (); k++)	//for all rotamers
			{
				//cout<<" k = "<<k<<" before l ..."<<endl;
				for (l=0; l<nList[i].selfE .size (); l++)		//for all other rotamers
				{
					if (k != l)
					{
						totalContributionE = nList[i].selfE [k] - nList[i].selfE [l];
						//cout<<nList[i].selfE [k]<<" "<<nList[i].selfE [l];
						for (j=0; j<nList[i].neighbor .size (); j++)		//for all neighbors
						{
							//cout<<" with "<<nList[i].neighbor [j]<<" ";
							minDiff = 10000;
							//for (int m=0; m<nList[i].selfE .size (); m++)
							for (int m=0; m<nList[nList[i].neighbor [j]].selfE .size (); m++)
							{
								PEdiff = nList[i].PE [j].neighborPE [k].PErotamer [m] - nList[i].PE [j].neighborPE [l].PErotamer [m];

								if (PEdiff < minDiff)
									minDiff = PEdiff;
							}
							totalContributionE += minDiff;				//add the minimal difference in energy between  these 2 rotamers
						}
						//cout<<"   with RotIndx "<<l<<" totalCont.E = "<<totalContributionE<<endl;;
						//if satisfies criteria then save the indx of k
						if (totalContributionE > 0)
						{
							rotsToBeDeleted.push_back(k);
							hasReducedRot = true;
							break;							//no need to continue to compare with rest rotamers
						}
					}
				}
			}
			/*****
					Delete Rotamers from the list
			*****/
			//cout<<nList[i].originalBB .AAs[0].chr3<<" "<<nList[i].originalBB .AAs[0].num<<" rotamers to be deleted ";
			//for (l=0;l<rotsToBeDeleted.size (); l++)
			//	cout<<rotsToBeDeleted[l]<<" ";
			//cout<<endl;

			for (j=0; j<rotsToBeDeleted.size ();j++)
				eraseRotamer(nList, i, rotsToBeDeleted[j] - j);
		}
		//cout<<i<<"  Done.."<<endl;
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool doneYet(vector<AAinformation> & nList)
{
	int nDoneAA = 0;
	for (int i=0; i<nList.size (); i++)
		if (!nList[i].valid)
			nDoneAA++;

	if (nDoneAA != nList.size ())
	{
	//	cout<<"doneYet ..... "<<nDoneAA<<" #AAin the nList= "<<nList.size ()<<endl;
		return false;
	}
	else
		return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void solveReducedAAStack(vector<AAinformation> & nList, vector<int> reducedAAs)
{
	int i;

	for (i= reducedAAs.size () - 1; i>= 0; i--)
	{
		//cout<<"working on "<<nList[reducedAAs[i]].originalBB .AAs[0].num<<endl;
		solveReducedAAs(nList, reducedAAs[i]);
	}

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void reduceAARandomly(vector<AAinformation> & nList)
{
	int i, j,
		minNumOfRot = 20,		//initial minimum number of rotamers
		AAindx;					//the indx of AA has minimum numbre of AA

	//search for AA has minimum number of rotamers
	for (i=0; i<nList.size (); i++)
		if (nList[i].valid)
			if (nList[i].selfE .size () < minNumOfRot)
			{
				minNumOfRot = nList[i].selfE .size ();
				AAindx = i;
			}

	//simply choose one rotamer (from the first 2 rotamers) to be selected according to its contribution to the energy
	double rot1E = nList[AAindx].selfE [0],
			rot2E= nList[AAindx].selfE [1];
	for (i=0; i<nList[AAindx].PE .size (); i++)
	{
		for (j=0; j<nList[AAindx].PE [i].neighborPE[0].PErotamer .size (); j++)
		{
			//the summation of PE b/w rotamer and all other rotamers for all neighbors
			rot1E += nList[AAindx].PE [i].neighborPE [0].PErotamer [j];
			rot2E += nList[AAindx].PE [i].neighborPE [1].PErotamer [j];
		}
	}

	//cout<<"AA has minimal number of rotamers is "<<nList[AAindx].originalBB .AAs[0].num <<" has "<<nList[AAindx].selfE .size ()<<" rotamers. "<<endl;

	vector<int> rotsToBeDeleted;
	if (rot1E <= rot2E)
		rotsToBeDeleted.push_back (1);
	else
		rotsToBeDeleted.push_back (0);

	//push all other rotamers
	for (i=2; i<nList[AAindx].selfE .size (); i++)
		rotsToBeDeleted.push_back (i);

	for (i=0; i<nList[AAindx].selfE .size (); i++)
		eraseRotamer(nList, AAindx, rotsToBeDeleted[i] - i);


}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//assumes that the array sideChains from rotamer.h is set...and BBDepRot is built
void predictSideChain(Protein & BBone)
{
	int i;
	vector<double> scR;			//vector that contains the maximal side chain radius of the AA


	/*****
		plug in initial side chains
	*****/
	for (i=0; i<BBone.numOfAA();i++)
		BBone.plugSideChain(i, sideChains[chr2Num(BBone.AAs[i].chr1)]);




	vector<AAinformation> nList;			//list of AA's and their neighbors list
	double eValue = 0;

//    cout<<"HERE.."<<endl;
//	unsigned long time1 = GetTickCount();
	setSCRadius(BBone, nList);		//set  scR	and build all rotamers in nList.allRotamers
//	unsigned long time2 = GetTickCount();
//	cout<<"setSCRadius took "<<time2-time1<<" ms"<<endl;
//    cout<<"radius set"<<endl;

	/*****
			compute Elib + VDW energy (self-energy term Cir) from R3 paper
	*****/
//	time1 = GetTickCount();
	setSelfEnergy(nList);
//	time2 = GetTickCount();
//	cout<<"setSelfEnergy took "<<time2 - time1<<" ms"<<endl;
//    cout<<"finished Elib and VDW calculations.."<<endl;

	/*****
			find neighbors of side chain....if the distance b/w two Cb is less than scR1 + scR2
	*****/
	//setSCNeighbors(nList);

	/*****
			determine neighbors of an AA and compute Pair energy b/w rotamers of all neighbors
	*****/
//	time1 = GetTickCount();
	setPairEnergy(nList);
//	time2 = GetTickCount();
//	cout<<"setPairEnergy took "<<time2-time1<<" ms "<<endl;
//    cout<<"pair energy is set.."<<endl;

	//printPEinfo(nList, "PEout_beforeFiltering.txt");

//	time1 = GetTickCount();
	filterNeighbors(nList, 0);		//Proposition 1 ...... filter the neighbor of an AA, because the previous one setPairEnergy is very general
//	time2 = GetTickCount();
//	cout<<"filterNeighbors took "<<time2 - time1<<" ms "<<endl;

	/*****
			main iteration....reduce residues and rotamers till you solve for all residues
	*****/
	int nAA = nList.size ();
	vector<int> reducedAAs;			//AA's have been reduced so far

	//printPEinfo(nList, "PEout_afterFiltering.txt");

//	Protein finalPortion[nAA];			//final portion that has final rotamers have been choosen


//	printRotamer(nList, 1);

//	getchar();
	bool cont = true;
//	time1 = GetTickCount();

	//cout<<"before loop..."<<endl;
	bool hasReducedRot = false;
	do
	{

		do
		{
			hasReducedRot = false;

			reduceResidues(nList, reducedAAs);

			reduceRotamers(nList, hasReducedRot);

			//cout<<"  =========== One iteration is Done ==================="<<endl;

			//printPEinfo(nList, "PEout_inIteration.txt");

			//getchar();
		} while (hasReducedRot);//(canBeReduced(nList));

		if (doneYet(nList))
			cont = false;
		else
			//residue Unification here;
			reduceAARandomly(nList);


	} while (cont);

	solveReducedAAStack(nList, reducedAAs);

//	time2 = GetTickCount();
//	cout<<"reduceResidues took "<<time2-time1<<" ms"<<endl;

//	getchar();
//	for (i=0; i<nList.size () ; i++)
//		cout<<nList[i].originalBB .AAs[0].num<<" rotamer chosen is "<<nList[i].chosenRotIndx <<endl;

	Protein finalStructure;
	getFinalStructure(nList, finalStructure);

	BBone = finalStructure;

	nList.clear ();
	reducedAAs.clear ();

	//printPEinfo(nList, "PEout_afterReduction.txt");

};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// End Of IMPLEMENTATION ////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
