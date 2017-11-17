#ifndef ALLORIENTATIONS_H
#define ALLORIENTATIONS_H

#include "energyfunction.h"

statistic * allOrientations;

//int nOrientations;

#define intervalSize 60		//rotate each stick evergy interval size

void enumerateAllOrientations(vector<statistic> & allOrientations, Protein nativePDBFile, vector<Protein> initialSkeleton, vector<sticksRec> sticksVect)
{
	//pow(2, nSticks)		total number of skeletons.....forward and backward
	//360/intervalSize the number of valid rotations for each skeleton
	//pow(nSticks, 360/intervalSize)		the number of valid conformations for each skeleton

	int nValidRotations = 360 / intervalSize;
	int nSticks = sticksVect[0].sticks.size();
//	nOrientations = sticksVect.size() * pow(2, nSticks) * pow(nValidRotations, nSticks);

//	allOrientations = new statistic[nOrientations];
	//cout<<"total number of orientations = "<<nOrientations<<endl;
	vector<int> stickDirections,oneVect, stickRotations;

	int i, j, k;
//	int orientationCntr = 0;
	vector<Protein> curConformation;
	long double intraE, interE;

	for (int sticksCntr=0; sticksCntr<sticksVect.size(); sticksCntr++)
	{
		int  directionCarry=0, rotationCarry = 0;
		//initiate directions and rotations vectors
		for (i=0; i<nSticks;i++)
		{
			stickDirections.push_back(0);
			oneVect.push_back(0);
			stickRotations.push_back(0);
		}

		oneVect[0] = 1;

		for (i=0; i<pow(2, nSticks); i++)
		{
			curConformation = initialSkeleton;
			for (k=nSticks-1; k>=0; k--)
			{
				sticksVect[sticksCntr].sticks[k].direction = stickDirections[k];
				AssignSeqToSkeleton(nativePDBFile, curConformation[k], 0, k, sticksVect[sticksCntr].sticks);
				//plug in side chain
				for (j=0; j<sticksVect[sticksCntr].sticks[k].length; j++)
					curConformation[k].plugSideChain (j, sideChains[chr2Num(curConformation[k].AAs[j].chr1)]);
				cout<<stickDirections[k]<<" ";
			}
			cout<<endl;

			for (k=0; k<pow(nValidRotations, nSticks); k++)
			{
				Protein allSticks;
				for (j=0;j<nSticks;j++)
				{
					//cout<<stickRotations[nSticks - j - 1]<<" ";
					sticksVect[sticksCntr].sticks[j].rotation = stickRotations[j] * intervalSize;
					rotateAroundAxis(curConformation[j],stickRotations[j] * intervalSize);	
					//append all sticks in one protein file for collision detection purposes
					allSticks.append (curConformation[j],0,curConformation[j].numOfAA()-1);
				}

				////////////// test and resolve collision
				int * arrOfLastChiUsed;
				initializeLastChiUsed(arrOfLastChiUsed, allSticks.numOfAA());	//initiate arrOfLastChiUsed

				// clear collision lists
				allSticks.AAsCollide .clear ();
				allSticks.atomsCollide .clear ();

				//resolve collision...if success or not
				int success = resolveCollision(allSticks, arrOfLastChiUsed);


				if (success)
				{
					interE = interEnergy(curConformation);
					intraE = 0;
					int totalnAA=0,nAA;
					for (j=0;j<nSticks;j++)
					{
						intraE += intraEnergy(curConformation[j],nAA);
						totalnAA += nAA;
					}

					intraE /= totalnAA;
				}
				else
				{
					interE = 1000;
					intraE = 1000;
				}

				for (j=0;j<nSticks;j++)
				{
					//cout<<"("<<res_vect[j]<<one_vect[j]<<carry<<")"<<endl;
					stickRotations[j] = stickRotations[j] + oneVect[j] + rotationCarry;			
					if (stickRotations[j] >= nValidRotations )
					{
						rotationCarry = 1;
						stickRotations[j] = 0;
					}
					else
						rotationCarry = 0;
				}

				statistic tmpStat;
				tmpStat.sticks = sticksVect[sticksCntr].sticks;
				tmpStat.intra  = intraE;
				tmpStat.inter  = interE;
				tmpStat.RMSD   = getRMSD(curConformation, nativePDBFile, sticksVect[sticksCntr].sticks);
				allOrientations.push_back(tmpStat);
			}

			//determine directions
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
		}
	}
}

void printAllOrientations(vector<statistic> allOrientations, string outFile)
{
	ofstream out;
	outFile += "_";
	outFile += "AllOrientations.txt";
	out.open(outFile.c_str());
	if (!out) {
		cerr << "Unable to open " << outFile << endl; 
		exit(1);
	}

	int i, j;
//	cout<<"nOrientations = "<<nOrientations<<endl;
	for (i=0; i<allOrientations.size(); i++)
	{
		out<<"assignments ( ";
		for (j=0; j<allOrientations[i].sticks.size();j++)
		{
			//cout<<"  j = "<<j<<endl;
			//out<<"stick "<<allOrientations[i].sticks[j].stickNum<<" hlx "<<allOrientations[i].sticks[j].hlxAssigned + 1
			//	<<" direction "<<allOrientations[i].sticks[j].direction<<" rotation "<<allOrientations[i].sticks[j].rotation<<endl;
			out<<allOrientations[i].sticks[j].hlxAssigned + 1<<" ";
		}	
		out<<") directions ( ";
		for (j=0; j<allOrientations[i].sticks.size();j++)
			out<<allOrientations[i].sticks[j].direction<<" ";
		out<<") Rotations ( ";
		for (j=0; j<allOrientations[i].sticks.size();j++)
			out<<setw(4)<<allOrientations[i].sticks[j].rotation<<" ";
		out<<")  ";

		out<<"intraE "<<setw(10)<<allOrientations[i].intra<<" interE "<<setw(10)<<allOrientations[i].inter<<" totalE "<<setw(10)<<allOrientations[i].inter + allOrientations[i].intra
			<<setw(8)<<" RMSD "<<allOrientations[i].RMSD<<endl;
	}
	out.close();
}

void NativeLikeAssignment(Protein nativeProtein, vector<Protein> initialSkeleton, vector<sticksRec> sticksVect)
{
	int i, j, k, nSticks = sticksVect[0].sticks.size();
	long double interE, intraE;
	
	for (i=0; i<sticksVect.size(); i++)
	{
		vector<Protein> curAssignment;
		curAssignment = initialSkeleton;
		Protein allSticks;
		for (j=0; j<nSticks; j++)
		{
			cout<<"stick "<<sticksVect[i].sticks[j].stickNum<<" hlx "<<sticksVect[i].sticks[j].hlxAssigned + 1<<" stick length = "<<sticksVect[i].sticks[j].length<<endl;
			sticksVect[i].sticks[j].direction = 1;
			AssignSeqToSkeleton(nativeProtein, curAssignment[j], 0, j, sticksVect[i].sticks);
			//get the angle b/w the first AA in the skeleton and the native Protein
			//double angle = getAngleDegree(curAssignment[j].getAtomCoordinate(0, " CA "), sticksVect[i].sticks[j].start, nativeProtein.getAtomCoordinate(sticksVect[i].sticks[j].startAAIndx, " CA "));
			double angle = getTorsionAngle(curAssignment[j].getAtomCoordinate(0, " CA "),
											sticksVect[i].sticks[j].start,
											sticksVect[i].sticks[j].end,
											nativeProtein.getAtomCoordinate(sticksVect[i].sticks[j].startAAIndx, " CA "));
			cout<<"Angle = "<<angle<<endl;
			rotateAroundAxis(curAssignment[j], angle);
			for (k=0; k<sticksVect[i].sticks[j].length; k++)
				curAssignment[j].plugSideChain (k, sideChains[chr2Num(curAssignment[j].AAs[k].chr1)]);
			allSticks.append (curAssignment[j],0,curAssignment[j].numOfAA()-1);
		}
		allSticks.writePDB("allSticks.pdb", 1, allSticks.numOfAA());
		cout<<"RMSD = "<<getRMSD(curAssignment, nativeProtein, sticksVect[i].sticks);
		////////////// test and resolve collision
		int * arrOfLastChiUsed;
		initializeLastChiUsed(arrOfLastChiUsed, allSticks.numOfAA());	//initiate arrOfLastChiUsed

		// clear collision lists
		allSticks.AAsCollide .clear ();
		allSticks.atomsCollide .clear ();

		//resolve collision...if success or not
		int success = resolveCollision(allSticks, arrOfLastChiUsed);

		if (success)
		{
			interE = interEnergy(curAssignment);
			intraE = 0;
			int totalnAA=0,nAA;
			for (j=0;j<nSticks;j++)
			{
				intraE += intraEnergy(curAssignment[j],nAA);
				totalnAA += nAA;
			}

			intraE /= totalnAA;
		}
		else
		{
			interE = 1000;
			intraE = 1000;
		}
		cout<<"   intraE = "<<intraE<<"  interE = "<<interE<<"  total = "<<(interE + intraE ) / 2<<endl;
		getchar();
	}

}
#endif