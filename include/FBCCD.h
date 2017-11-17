#ifndef FBCCD_H
#define FBCCD_H

#include "geometry.h"
#include "constants.h"
#include "skeleton_overall.h"



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// List Of Data Structures //////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// End Of Constants and Data Structures /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// List of Functions ////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// ////////// End Of The List //////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////// Functions Implementations //////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//aaTypes is the list of the types of amino acids in the structure
Protein generateRandomStructure(vector<char> aaTypes, int seed){
    Protein sModel;
    AminoAcid tmpAA;

    int i;
    //create a protein with number of aa acids
    for (i=0; i<aaTypes.size(); i++){
        sModel.AAs.push_back(tmpAA);
    }

    float *hListPhi, *sListPhi, *lListPhi, *hListPsi, *sListPsi, *lListPsi;
    hListPhi = new float[aaTypes.size()];
    sListPhi = new float[aaTypes.size()];
    lListPhi = new float[aaTypes.size()];
    hListPsi = new float[aaTypes.size()];
    sListPsi = new float[aaTypes.size()];
    lListPsi = new float[aaTypes.size()];

    getRandomList(MIN_PHI_HLX_ALLOWED, MAX_PHI_HLX_ALLOWED, hListPhi, aaTypes.size(), seed);
    getRandomList(MIN_PHI_BETA_ALLOWED,MAX_PHI_BETA_ALLOWED, sListPhi, aaTypes.size(), seed);
    getRandomList(MIN_PHI_LOOP_ALLOWED, MAX_PHI_LOOP_ALLOWED, lListPhi, aaTypes.size(), seed);

    getRandomList(MIN_PSI_HLX_ALLOWED,MAX_PSI_HLX_ALLOWED, hListPsi, aaTypes.size(), seed);
    getRandomList(MIN_PSI_BETA_ALLOWED,MAX_PSI_BETA_ALLOWED, sListPsi, aaTypes.size(), seed);
    getRandomList(MIN_PSI_LOOP_ALLOWED, MAX_PSI_LOOP_ALLOWED, lListPsi, aaTypes.size(), seed);


    for (i=0; i<aaTypes.size(); i++){
        if (aaTypes[i] == 'H'){
            //sModel.AAs[i].angles.phi = getRandomFloat(MIN_PHI_HLX_ALLOWED, MAX_PHI_HLX_ALLOWED , seed*i*(i+1));	//alpha helix phi;
            //sModel.AAs[i].angles.psi = getRandomFloat(MIN_PSI_HLX_ALLOWED,MAX_PSI_HLX_ALLOWED, seed*i*(i+1));	//alpha helix psi;

            sModel.AAs[i].angles.phi = hListPhi[i];
            sModel.AAs[i].angles.psi = hListPsi[i];

            //cout<<"randomly generated phi= "<<sModel.AAs[i].angles.phi<<"  randomly generated Psi= "<<sModel.AAs[i].angles.psi<<endl;
            continue;
        }
        if (aaTypes[i] == 'S'){
            //sModel.AAs[i].angles.phi = getRandomFloat(MIN_PHI_BETA_ALLOWED,MAX_PHI_BETA_ALLOWED, seed*i*(i+1)); //sheet phi;
            //sModel.AAs[i].angles.psi = getRandomFloat(MIN_PSI_BETA_ALLOWED,MAX_PSI_BETA_ALLOWED, seed*i*(i+1)); //sheet psi;

            sModel.AAs[i].angles.phi = sListPhi[i];
            sModel.AAs[i].angles.psi = sListPsi[i];
            //cout<<"randomly generated phi= "<<sModel.AAs[i].angles.phi<<"  randomly generated Psi= "<<sModel.AAs[i].angles.psi<<endl;
            continue;
        }
        if (aaTypes[i] == 'L'){
            //sModel.AAs[i].angles.phi = getRandomFloat(MIN_PHI_LOOP_ALLOWED, MAX_PHI_LOOP_ALLOWED, seed*i*(i+1));	//random phi;
            //sModel.AAs[i].angles.psi = getRandomFloat(MIN_PSI_LOOP_ALLOWED, MAX_PSI_LOOP_ALLOWED, seed*i*(i+1));		//random psi;

            sModel.AAs[i].angles.phi = lListPhi[i];
            sModel.AAs[i].angles.psi = lListPsi[i];

            //cout<<"randomly generated phi= "<<sModel.AAs[i].angles.phi<<"  randomly generated Psi= "<<sModel.AAs[i].angles.psi<<endl;
        }
    }
    //build the 3-d structure from phi-psi model...number of AA in sModel now is ssType.size()+2
    sModel.torsion2coord();


    for (i=0; i<aaTypes.size(); i++){
        sModel.AAs[i].SStype = aaTypes[i];
    }

    return sModel;
}
//////////////////////////////// Find the Minimum Distance and Angle
// moving is a list of coordinate that you wanna move to overlap target list
// number of points in moving should be equal to number of points in target
// p2 and p3 is the 2 atoms of the phi or psi bond
// for phi --> p2 is Ca and p3 is N
// for psi --> p2 is C and p3 is Ca
// nDist is the RMSD between moving points and target points after you rotate
float getMinAngle (vector<Coordinate> moving,Coordinate p2,Coordinate p3,vector<Coordinate> target,float &nDist){


	Vectors fiv, riv, riu, theta_u, sin, rotationaxis;
	double a=0, b=0, c=0;
	Coordinate oi, Fi, Moi;
	vector<double> ai;
	vector<double> bi;
	vector<double> ci;

	rotationaxis.set(p2,p3);
	theta_u = rotationaxis.unit();


	int i;
	for (i=0; i<moving.size(); i++)
	{
		Fi	=	target[i];
		Moi	=	moving[i];
		oi	=	pointLineIntersection(Moi, p2, p3);
		fiv.set(Fi, oi);
		riv.set(Moi, oi);
		riu	=	riv.unit();
		sin = riu.cross(theta_u);

		ai.push_back (SQR(riv.length()) + SQR(fiv.length()));
		bi.push_back (2 * riv.length() * fiv.dot(riu));
		ci.push_back (2 * riv.length() * fiv.dot(sin));

		a += ai[i];
		b += bi[i];
		c += ci[i];

	}


	float dihedrali,									//dihedral angle (mi, p2, p3, ti)  .. in radian
		  si,											//the expected distance between moving point i and target point i
		  theta = acos(b/sqrt(SQR(b)+SQR(c)));			//in radian



	nDist=0;
	for (i=0; i<moving.size(); i++)
	{
		dihedrali = fabs(getTorsionAngle(moving[i],p2,p3,target[i]));		//find dihedral angle for current conformation
		dihedrali = toRadian(dihedrali);									//convert to radian
		si = ai[i] - sqrt(SQR(bi[i])+SQR(ci[i]))*cos(dihedrali-theta);		//calculate the expected distance
		nDist += si;

		/*
		cout<<" === IN getMinAngle === "<<endl;
		cout<<"  dihedral "<<i+1<<" = "<<dihedrali<<"  theta= "<<theta<<endl;
		cout<<"  ai= "<<ai[i]<<" bi= "<<bi[i]<<" ci= "<<ci[i]<<" cos= "<<cos(dihedrali-theta)<<endl;
		cout<<"  dist "<<i+1<<" = "<<si<<endl<<endl;
		*/
	}

	nDist = sqrt(nDist/moving.size());		//RMSD

	theta = toDegree(theta);

	if (atan2(asin(c/sqrt(SQR(b)+SQR(c))),acos(b/sqrt(SQR(b)+SQR(c))))<0){
		if (theta>0)
			theta *=-1;
	}
	else{
		if (theta<0)
			theta *=-1;
	}

	return theta;
}
//////////////////////////////// Find the Minimum Distance and Angle
// moving is a list of coordinate that you wanna move to overlap target list
// number of points in moving should be equal to number of points in target
// p2 and p3 is the 2 atoms of the phi or psi bond
// for phi --> p2 is Ca and p3 is N
// for psi --> p2 is C and p3 is Ca
// nDist is the RMSD between moving points and target points after you rotate
// sIndx		the indx of the first point to overlap
// eIndx		the indx of the last point to overlap
float getMinAngle (vector<Coordinate> moving,
				   Coordinate p2,
				   Coordinate p3,
				   vector<Coordinate> target,
				   short sIndx,
				   short eIndx,
				   float &nDist){


	Vectors fiv, riv, riu, theta_u, sin, rotationaxis;
	double a=0, b=0, c=0;
	Coordinate oi, Fi, Moi;
	vector<double> ai;
	vector<double> bi;
	vector<double> ci;

	rotationaxis.set(p2,p3);
	theta_u = rotationaxis.unit();


	int i,
		cntr=0;
	for (i=sIndx; i<=eIndx; i++)
	{
		Fi	=	target[i];
		Moi	=	moving[i];
		oi	=	pointLineIntersection(Moi, p2, p3);
		fiv.set(Fi, oi);
		riv.set(Moi, oi);
		riu	=	riv.unit();
		sin = riu.cross(theta_u);

		ai.push_back (SQR(riv.length()) + SQR(fiv.length()));
		bi.push_back (2 * riv.length() * fiv.dot(riu));
		ci.push_back (2 * riv.length() * fiv.dot(sin));

		a += ai[cntr];
		b += bi[cntr];
		c += ci[cntr];
		cntr++;
	}


	float dihedrali,									//dihedral angle (mi, p2, p3, ti)  .. in radian
		  si,											//the expected distance between moving point i and target point i
		  theta = acos(b/sqrt(SQR(b)+SQR(c)));			//in radian



	nDist=0;
	cntr=0;
	for (i=sIndx; i<=eIndx; i++)
	{
		dihedrali = fabs(getTorsionAngle(moving[i],p2,p3,target[i]));		//find dihedral angle for current conformation
		dihedrali = toRadian(dihedrali);									//convert to radian
		si = ai[cntr] - sqrt(SQR(bi[cntr])+SQR(ci[cntr]))*cos(dihedrali-theta);		//calculate the expected distance
		nDist += si;
		cntr++;

		/*
		cout<<" === IN getMinAngle === "<<endl;
		cout<<"  dihedral "<<i+1<<" = "<<dihedrali<<"  theta= "<<theta<<endl;
		cout<<"  ai= "<<ai[i]<<" bi= "<<bi[i]<<" ci= "<<ci[i]<<" cos= "<<cos(dihedrali-theta)<<endl;
		cout<<"  dist "<<i+1<<" = "<<si<<endl<<endl;
		*/

	}

	nDist = sqrt(nDist/(eIndx-sIndx+1));		//RMSD

	theta = toDegree(theta);

	if (atan2(asin(c/sqrt(SQR(b)+SQR(c))),acos(b/sqrt(SQR(b)+SQR(c))))<0){
		if (theta>0)
			theta *=-1;
	}
	else{
		if (theta<0)
			theta *=-1;
	}

	return theta;
}
///////////////////////////////////////////////////
// update moving points coordinate after rotation
// calculation of rotation matrix would be computationaly expensive for more iterations
// this could be avoided if we instead of calculate the matrix, just get it from Protein.rotate function
// which is already built, and we will use the same rotation matrix b/s we are rotating around same line
void updatePoints(vector<Coordinate> &mPoints,
				  Coordinate p2,
				  Coordinate p3,
				  float angle)
{
	int i;

	double m[4][4];
	buildRotationMatrix(m, p2, p3, -toRadian(angle));

	for (i=0; i<mPoints.size(); i++)
		rotatePoint(m, mPoints[i]);

}
///////////////////////////////////////////////////
// update moving points coordinate after rotation
// calculation of rotation matrix would be computationaly expensive for more iterations
// this could be avoided if we instead of calculate the matrix, just get it from Protein.rotate function
// which is already built, and we will use the same rotation matrix b/s we are rotating around same line
void updatePoints(vector<Coordinate> &mPoints,
				  Coordinate p2,
				  Coordinate p3,
				  short sIndx,
				  short eIndx,
				  float angle)
{
	int i;

	double m[4][4];
	buildRotationMatrix(m, p2, p3, -toRadian(angle));

	for (i=sIndx; i<=eIndx; i++)
		rotatePoint(m, mPoints[i]);

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
float getRMSD(vector<Coordinate> &pnts1, vector<Coordinate> &pnts2){
	int i;
	float xDiff, yDiff, zDiff, total=0;
	for (i=0; i<pnts1.size (); i++){
		xDiff = pnts1[i].x - pnts2[i].x;
		yDiff = pnts1[i].y - pnts2[i].y;
		zDiff = pnts1[i].z - pnts2[i].z;
		total = total + SQR(xDiff) + SQR(yDiff) + SQR(zDiff);
	}
	total = sqrt(total/pnts1.size ());

	return total;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//check of the applying of an angle to a toirsion angle will cause a collision with the other part in the structure
//startAA is the start AA to start check from
bool collisionFree(Protein &mStructure, vector<vector<Coordinate> > &avoidTrace, short startAA, short endAA, short cAA, string torsion, float deltaAngle){

		Protein tmp = mStructure;
		Coordinate p2, p3;
		int i, j;

		//apply the angle on the structure
		if (torsion == "psi"){
			//Working on Psi
			p2 = tmp.getAtomCoordinate(cAA, " CA ");
			p3 = tmp.getAtomCoordinate(cAA, " C  ");
			//rotate
			tmp.rotate(cAA, tmp.numOfAA()-1, tmp.getAtomIndx(cAA, " C  ")+1, p2, p3, deltaAngle);
		}
		else{
			//Working on phi .... phi of the next AA...
			p2 = tmp.getAtomCoordinate(cAA, " N  ");
			p3 = tmp.getAtomCoordinate(cAA, " CA ");

			//rotate
			tmp.rotate(cAA, tmp.numOfAA()-1, tmp.getAtomIndx(cAA, " CA ")+1, p2, p3, deltaAngle);
		}

		//cout<<"inCollisionFree start AA is "<<mStructure.AAs[startAA].num<<" endAA is "<<mStructure.AAs[endAA].num<<" cAA is "<<tmp.AAs[cAA].num<<endl;
		//check for collision
		Coordinate curCA;
		//for (i=cAA; i<=endAA; i++){
		for (i=cAA+1; i<=endAA; i++){
			//get CA atom
			curCA = tmp.getAtomCoordinate(i, " CA ");
			//should not collide with itself
			for (j=startAA; j<i; j++){
				if (getDistance(curCA, tmp.getAtomCoordinate(j, " CA "))<2.5){
					//cout<<"AA num "<<tmp.AAs[i].num<<" collides with AA num "<<tmp.AAs[j].num<<" dist= "<<getDistance(curCA, tmp.getAtomCoordinate(j, " CA "))<<endl;
					//tmp.writePDB("prot_collision.pdb", 1, tmp.numOfAA());getchar();
					return false;
				}
			}
			//should not collide with the other part of protein represented by a trace
			for (j=0; j<avoidTrace.size(); j++){
				for (int k =0; k<avoidTrace[j].size()-1; k++){
					//get the distance b/w the curCA point and the segment of the line
					if (getDistLineSegPoint(avoidTrace[j][k],avoidTrace[j][k+1], curCA)<4){
						//cout<<"AA num "<<tmp.AAs[i].num<<" collides with line segment "<<j+1<<" - "<<k+1<<" dist= "<<getDistLineSegPoint(avoidTrace[j][k],avoidTrace[j][k+1], curCA)<<endl;
						//cout<<avoidTrace[j][k].x<<" "<<avoidTrace[j][k].y<<" "<<avoidTrace[j][k].z<<" and "<<avoidTrace[j][k+1].x<<" "<<avoidTrace[j][k+1].y<<" "<<avoidTrace[j][k+1].z<<endl;
						//tmp.writePDB("prot_collision.pdb", 1, tmp.numOfAA());getchar();
						return false;
					}
				}
			}

		}

	return true;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// check if the new torsion is valid
// type is the type of the structure phi, or psi
bool validTorsion(float torsion, char sType, string type)
{
	if (sType == 'H')				//for Helix
	{
		if (type == "phi")
			if ((torsion >= MIN_PHI_HLX_ALLOWED) && (torsion <= MAX_PHI_HLX_ALLOWED))
				return true;
			else
				return false;
		if (type == "psi")
			if ((torsion >= MIN_PSI_HLX_ALLOWED) && (torsion<= MAX_PSI_HLX_ALLOWED))
				return true;
			else
				return false;
	}
	else if (sType == 'S')			//for Strand
	{
		if (type == "phi")
			if ((torsion >= MIN_PHI_BETA_ALLOWED) && (torsion <= MAX_PHI_BETA_ALLOWED))
				return true;
			else
				return false;
		if (type == "psi")
			if ((torsion >= MIN_PSI_BETA_ALLOWED) && (torsion<= MAX_PSI_BETA_ALLOWED))
				return true;
			else
				return false;
	}
	else							// for Loop
	{
		if (type == "phi")
			if ((torsion >= MIN_PHI_LOOP_ALLOWED) && (torsion <= MAX_PHI_LOOP_ALLOWED))
				return true;
			else
				return false;
		if (type == "psi")
			if ((torsion >= MIN_PSI_LOOP_ALLOWED) && (torsion<= MAX_PSI_LOOP_ALLOWED))
				return true;
			else
				return false;
	}
}

//////////
// Forward Backward Cyclic Coordinate Descent algorithm
float FBCCD(Protein &mStructure,					//the structure we are going to change so mPoints (which are some points in this mStructure) move to be as close as possible to tPoints (target points)
		  // char sType,							//the type of the structure H, S, or L
		   short sAAindx,						//the indx of the first AA to update...be carefull of the direction
		   short eAAindx,						//the last AA in mStructure to update.... it could be the beginning or the end of the structure according to the direction
		   short direction,						//the firection of the walk (1 forward, -1 backward);
		   vector<Coordinate> &mPoints,			//moving points either on mStructure or some other reference points (such as center of AA)
		   vector<Coordinate> tPoints,			//target points (Coordinates)
		   float RMSDthr,						//RMSD cut off threshhold
		   int tIteration)						//number of iteration cut off
{

	Coordinate p2, p3;
	int itCntr = 0,
		i,
		cAA = sAAindx;		//current AA we are working on
	float nRMSD = RMSDthr + 1,	//new expected RMSD
		  torsionDelta,			//the angle we should rotate
		  curTorsion,			//current torsion angle
		  nTorsion;				//new (the expected) torsion after rotation
	char sType;					//the type of the current AA we are working on
    //variable to get the range of phi and psi angles for a particular AA in the structure
	float	min_phi_allowed,
			max_phi_allowed,
			min_psi_allowed,
			max_psi_allowed;

	tIteration *= (abs(eAAindx-sAAindx)+1);

	while ((nRMSD > RMSDthr) && (itCntr < tIteration))
	{
		//cout<<itCntr+1<<" : Working on AA indx = "<<cAA<<endl;

		//Working on Psi
		p2 = mStructure.getAtomCoordinate(cAA, " CA ");
		p3 = mStructure.getAtomCoordinate(cAA, " C  ");

		torsionDelta = getMinAngle(mPoints, p3, p2, tPoints, nRMSD);
		curTorsion = mStructure.getTorsion("psi", cAA);
		nTorsion = curTorsion - torsionDelta;		//new torsion after you rotate

		if (nTorsion<-180)
			nTorsion = 360+nTorsion;
		else if (nTorsion>180)
			nTorsion = 360-nTorsion;
		//cout<<"================================="<<endl;
		//cout<<"  curTorsion= "<<curTorsion<<"  deltaAngle= "<<torsionDelta<<"  expected nRMSD= "<<nRMSD<<endl;

		//get the type of the amino acid...H, L, or S
		sType = mStructure.AAs[cAA].SStype;
        if (sType=='H')
        {
            min_phi_allowed = MIN_PHI_HLX_ALLOWED;
            max_phi_allowed = MAX_PHI_HLX_ALLOWED;
            min_psi_allowed = MIN_PSI_HLX_ALLOWED;
            max_psi_allowed = MAX_PSI_HLX_ALLOWED;
        }
        else if (sType=='S')
        {
            min_phi_allowed = MIN_PHI_BETA_ALLOWED;
            max_phi_allowed = MAX_PHI_BETA_ALLOWED;
            min_psi_allowed = MIN_PSI_BETA_ALLOWED;
            max_psi_allowed = MAX_PSI_BETA_ALLOWED;
        }
        else
        {
            min_phi_allowed = MIN_PHI_LOOP_ALLOWED;
            max_phi_allowed = MAX_PHI_LOOP_ALLOWED;
            min_psi_allowed = MIN_PSI_LOOP_ALLOWED;
            max_psi_allowed = MAX_PSI_LOOP_ALLOWED;
        }

		if (validTorsion(nTorsion, sType, "psi"))		//mStructure.getTorsion("psi", cAA) - torsionDelta : is the torsion angle after we rotate
		{
			//rotate
			//mStructure.rotate(cAA, eAAindx, mStructure.getAtomIndx(cAA, " C  ")+1, p2, p3, torsionDelta);
			mStructure.rotate(cAA, mStructure.numOfAA()-1, mStructure.getAtomIndx(cAA, " C  ")+1, p2, p3, torsionDelta);


			//udate moving points....because of rotation..we need to update these points everytime
			updatePoints(mPoints, p2,p3, torsionDelta);
		}
		// you can use else section in case the angle is not valid then you rotate the maximum valid angle or a specific angle
		// for example if the new angle is out of valid range, suppose valid range is [-40, 30], then you can rotate so the
		// new torsion is -40 or 30. or you can rotate certain angle, for example 10 degrees or 20 degrees so you still in the range
		// if you don't like to do so...just comment else section
		/*
		else
		{
			//cout<<"  NOT VALID PSI.."<<endl;
			//case #1 : you want to rotate so you are on the boundary of the valid range (you rotate the maximum you can)
			// first determine which boundary...which boundary is closer
			float delta1, delta2;

			delta1 = abs(nTorsion - MIN_PSI_ALLOWED);
			if (delta1 > 180)
				delta1 = 360 - delta1;
			delta2 = abs(nTorsion - MAX_PSI_ALLOWED);
			if (delta2 > 180)
				delta2 = 360 - delta2;

			//cout<<"  nTorsion= "<<nTorsion<<" delta1= "<<delta1<<" delta2= "<<delta2<<endl;

			if (delta1 <= delta2)
			{
				//rotate to MIN boundary
				mStructure.rotate(cAA, eAAindx, mStructure.getAtomIndx(cAA, " C  ")+1, p2, p3, curTorsion-MIN_PSI_ALLOWED);
				//udate moving points....because of rotation..we need to update these points everytime
				updatePoints(mPoints, p2,p3, curTorsion-MIN_PSI_ALLOWED);
			}
			else
			{
				//rotate to Max boundary
				mStructure.rotate(cAA, eAAindx, mStructure.getAtomIndx(cAA, " C  ")+1, p2, p3, curTorsion-MAX_PSI_ALLOWED);
				//udate moving points....because of rotation..we need to update these points everytime
				updatePoints(mPoints, p2,p3, curTorsion-MAX_PSI_ALLOWED);
			}

			//cout<<"  PSI now= "<<mStructure.getTorsion("psi", cAA)<<endl;
		}
		*/

		//check the RMSD...
		//because FBCCD is not greedy..so it is not guranteed that any rotation will lead to better RMSD
		//so, don't rotate any more if you reach the desired RMSD
		if (nRMSD <= RMSDthr)		//in case you reach the RMSD you desire...
			break;

		//Working on phi .... phi of the next AA...
		p2 = mStructure.getAtomCoordinate(cAA+1, " N  ");
		p3 = mStructure.getAtomCoordinate(cAA+1, " CA ");

		torsionDelta = getMinAngle(mPoints, p3, p2, tPoints, nRMSD);
		curTorsion = mStructure.getTorsion("phi", cAA+1);
		nTorsion = curTorsion - torsionDelta;				//new torsion after you rotate

		if (nTorsion<-180)
			nTorsion = 360+nTorsion;
		else if (nTorsion>180)
			nTorsion = 360-nTorsion;
		//cout<<"  curTorsion= "<<curTorsion<<"  deltaAngle= "<<torsionDelta<<"  expected nRMSD= "<<nRMSD<<endl;

		//get the type of the amino acid...H, L, or S
		sType = mStructure.AAs[cAA+1].SStype;
        if (sType=='H')
        {
            min_phi_allowed = MIN_PHI_HLX_ALLOWED;
            max_phi_allowed = MAX_PHI_HLX_ALLOWED;
            min_psi_allowed = MIN_PSI_HLX_ALLOWED;
            max_psi_allowed = MAX_PSI_HLX_ALLOWED;
        }
        else if (sType=='S')
        {
            min_phi_allowed = MIN_PHI_BETA_ALLOWED;
            max_phi_allowed = MAX_PHI_BETA_ALLOWED;
            min_psi_allowed = MIN_PSI_BETA_ALLOWED;
            max_psi_allowed = MAX_PSI_BETA_ALLOWED;
        }
        else
        {
            min_phi_allowed = MIN_PHI_LOOP_ALLOWED;
            max_phi_allowed = MAX_PHI_LOOP_ALLOWED;
            min_psi_allowed = MIN_PSI_LOOP_ALLOWED;
            max_psi_allowed = MAX_PSI_LOOP_ALLOWED;
        }

		if (validTorsion(nTorsion, sType, "phi"))	//mStructure.getTorsion("phi", cAA) - torsionDelta : is the torsion angle after we rotate
		{
			//rotate
			//mStructure.rotate(cAA+1, eAAindx, mStructure.getAtomIndx(cAA+1, " CA ")+1, p2, p3, torsionDelta);
			mStructure.rotate(cAA+1, mStructure.numOfAA()-1, mStructure.getAtomIndx(cAA+1, " CA ")+1, p2, p3, torsionDelta);

			//udate moving points....because of rotation..we need to update these points everytime
			updatePoints(mPoints, p2,p3, torsionDelta);
		}
		// you can use else section in case the angle is not valid then you rotate the maximum valid angle or a specific angle
		// for example if the new angle is out of valid range, suppose valid range is [-40, 30], then you can rotate so the
		// new torsion is -40 or 30. or you can rotate certain angle, for example 10 degrees or 20 degrees so you still in the range
		// if you don't like to do so...just comment else section
		/*
		else
		{
			//cout<<"  NOT VALID PHI.."<<endl;
			//case #1 : you want to rotate so you are on the boundary of the valid range (you rotate the maximum you can)
			// first determine which boundary...which boundary is closer
			float delta1, delta2;

			delta1 = abs(nTorsion - MIN_PHI_ALLOWED);
			if (delta1 > 180)
				delta1 = 360 - delta1;
			delta2 = abs(nTorsion - MAX_PHI_ALLOWED);
			if (delta2 > 180)
				delta2 = 360 - delta2;

			//cout<<"  nTorsion= "<<nTorsion<<" delta1= "<<delta1<<" delta2= "<<delta2<<endl;

			if (delta1 <= delta2)
			{
				//rotate
				mStructure.rotate(cAA+1, eAAindx, mStructure.getAtomIndx(cAA+1, " CA ")+1, p2, p3, curTorsion-MIN_PHI_ALLOWED);
				//udate moving points....because of rotation..we need to update these points everytime
				updatePoints(mPoints, p2,p3, curTorsion-MIN_PHI_ALLOWED);
			}
			else
			{
				//rotate
				mStructure.rotate(cAA+1, eAAindx, mStructure.getAtomIndx(cAA+1, " CA ")+1, p2, p3, curTorsion-MAX_PHI_ALLOWED);
				//udate moving points....because of rotation..we need to update these points everytime
				updatePoints(mPoints, p2,p3, curTorsion-MAX_PHI_ALLOWED);
			}

			//cout<<"  PHI now= "<<mStructure.getTorsion("phi", cAA+1)<<endl;
		}
		*/


		cAA += direction;

		itCntr++;

		if (cAA == eAAindx)
			cAA = sAAindx;
	}

	return nRMSD;
}
///////////////////////////////////////////////////////////////////////////////////////
// modify a structure by FBCCD so mPoints overlap tPoints
// this function was mainly implemented to generate helices and strands for a given axis (which represents the curved hlx or strand)
void overlapFBCCD(Protein &mStructure,			//the structure we are going to change so mPoints (which are some points in this mStructure) move to be as close as possible to tPoints (target points)
				  //char sType,					//structure type H, S, or L
				  short nAAperSegment,			//number of Amino Acids between two consecutive points in mPoints and tPoints
				  vector<Coordinate> mPoints,	//moving points either on mStructure or some other reference points (such as center of AA)
				  vector<Coordinate> tPoints,	//target points (Coordinates)
				  float RMSDthr,				//RMSD cut off threshhold
				  int tIteration)				//number of iteration cut off
{

	Coordinate p2, p3;
	int itCntr = 0,
		i,
		pCntr,
		cAA;					//current AA we are working on
	float nRMSD = RMSDthr + 1,	//new expected RMSD
		  torsionDelta,			//the angle we should rotate
		  curTorsion,			//current torsion angle
		  nTorsion;				//new (the expected) torsion after rotation

	float	min_phi_allowed,
			max_phi_allowed,
			min_psi_allowed,
			max_psi_allowed;

	char sType;			//the type of amino acid we are working on (H, S, L)

	vector<vector<Coordinate> > emptyTrace;


	//cout<<"#mPoints= "<<mPoints.size()<<" #tPoints= "<<tPoints.size()<<endl;
	//start from point number two...overlap two consecutive points everytime
	for (pCntr=1; pCntr<mPoints.size()-1; pCntr++)
	{
		nRMSD = RMSDthr + 1;
		itCntr=0;
		cAA = nAAperSegment*(pCntr-1)+1;		//start from the following AA after the last point u worked on

		while ((nRMSD > RMSDthr) && (itCntr < tIteration))
		{

			sType = mStructure.AAs[cAA].SStype;
			if (sType=='H')
			{
				min_phi_allowed = MIN_PHI_HLX_ALLOWED;
				max_phi_allowed = MAX_PHI_HLX_ALLOWED;
				min_psi_allowed = MIN_PSI_HLX_ALLOWED;
				max_psi_allowed = MAX_PSI_HLX_ALLOWED;
			}
			else if (sType=='S')
			{
				min_phi_allowed = MIN_PHI_BETA_ALLOWED;
				max_phi_allowed = MAX_PHI_BETA_ALLOWED;
				min_psi_allowed = MIN_PSI_BETA_ALLOWED;
				max_psi_allowed = MAX_PSI_BETA_ALLOWED;
			}
			else
			{
				min_phi_allowed = MIN_PHI_LOOP_ALLOWED;
				max_phi_allowed = MAX_PHI_LOOP_ALLOWED;
				min_psi_allowed = MIN_PSI_LOOP_ALLOWED;
				max_psi_allowed = MAX_PSI_LOOP_ALLOWED;
			}
			//cout<<itCntr+1<<" : Working on AA indx = "<<cAA<<" pCntr= "<<pCntr<<endl;

			//Working on Psi
			p2 = mStructure.getAtomCoordinate(cAA, " CA ");
			p3 = mStructure.getAtomCoordinate(cAA, " C  ");

			torsionDelta = getMinAngle(mPoints, p3, p2, tPoints, pCntr, pCntr+1, nRMSD);
			curTorsion = mStructure.getTorsion("psi", cAA);
			nTorsion = curTorsion - torsionDelta;		//new torsion after you rotate


			//cout<<"  curTorsion psi= "<<curTorsion<<"  deltaAngle= "<<torsionDelta<<" nTorsion= "<<nTorsion<<" nRMSD= "<<nRMSD<<endl;
			//cout<<"before psi: mPoints["<<pCntr<<"].x= "<<mPoints[pCntr].x<<" y= "<<mPoints[pCntr].y<<" z= "<<mPoints[pCntr].z<<endl;

			if (validTorsion(nTorsion, sType, "psi")){// &&
				//collisionFree(mStructure, emptyTrace, nAAperSegment*(pCntr-1)+1, MIN(nAAperSegment*(pCntr+1) - 2, mStructure.numOfAA()-3), cAA, "psi", torsionDelta)){
				//rotate
				mStructure.rotate(cAA, mStructure.numOfAA()-1, mStructure.getAtomIndx(cAA, " C  ")+1, p2, p3, torsionDelta);

				//udate moving points....because of rotation..we need to update these points everytime
				updatePoints(mPoints, p2,p3, pCntr, mPoints.size()-1, torsionDelta);
			}
			// you can use else section in case the angle is not valid then you rotate the maximum valid angle or a specific angle
			// for example if the new angle is out of valid range, suppose valid range is [-40, 30], then you can rotate so the
			// new torsion is -40 or 30. or you can rotate certain angle, for example 10 degrees or 20 degrees so you still in the range
			// if you don't like to do so...just comment else section
			else
			{
				//cout<<"  NOT VALID PSI.."<<endl;
				//case #1 : you want to rotate so you are on the boundary of the valid range (you rotate the maximum you can)
				// first determine which boundary...which boundary is closer
				float delta1, delta2;

				delta1 = fabs(nTorsion - min_psi_allowed);
				if (delta1 > 180)
					delta1 = 360 - delta1;
				delta2 = fabs(nTorsion - max_psi_allowed);
				if (delta2 > 180)
					delta2 = 360 - delta2;

				//cout<<"  nTorsion= "<<nTorsion<<" delta1= "<<delta1<<" delta2= "<<delta2<<endl;

				if (delta1 <= delta2){// &&
					//collisionFree(mStructure, emptyTrace, nAAperSegment*(pCntr-1)+1, MIN(nAAperSegment*(pCntr+1) - 2, mStructure.numOfAA()-3), cAA, "psi", torsionDelta)){
					//rotate to MIN boundary
					mStructure.rotate(cAA, mStructure.numOfAA()-1, mStructure.getAtomIndx(cAA, " C  ")+1, p2, p3, curTorsion - min_psi_allowed);
					//udate moving points....because of rotation..we need to update these points everytime
					updatePoints(mPoints, p2,p3, pCntr, mPoints.size()-1, curTorsion-min_psi_allowed);
				}
				else
				{
					//if (collisionFree(mStructure, emptyTrace, nAAperSegment*(pCntr-1)+1, MIN(nAAperSegment*(pCntr+1) - 2, mStructure.numOfAA()-3), cAA, "psi", torsionDelta)){
						//rotate to Max boundary
						mStructure.rotate(cAA, mStructure.numOfAA()-1, mStructure.getAtomIndx(cAA, " C  ")+1, p2, p3, curTorsion - max_psi_allowed);
						//udate moving points....because of rotation..we need to update these points everytime
						updatePoints(mPoints, p2,p3, pCntr, mPoints.size()-1, curTorsion - max_psi_allowed);
					//}
				}

				//cout<<"  PSI now= "<<mStructure.getTorsion("psi", cAA)<<endl;
			}

			//cout<<"after psi: mPoints["<<pCntr<<"].x= "<<mPoints[pCntr].x<<" y= "<<mPoints[pCntr].y<<" z= "<<mPoints[pCntr].z<<endl;


			//check the RMSD...
			//because FBCCD is not greedy..so it is not guranteed that any rotation will lead to better RMSD
			//so, don't rotate any more if you reach the desired RMSD
			if (nRMSD <= RMSDthr)		//in case you reach the RMSD you desire...
				break;

			sType = mStructure.AAs[cAA+1].SStype;
			if (sType=='H')
			{
				min_phi_allowed = MIN_PHI_HLX_ALLOWED;
				max_phi_allowed = MAX_PHI_HLX_ALLOWED;
				min_psi_allowed = MIN_PSI_HLX_ALLOWED;
				max_psi_allowed = MAX_PSI_HLX_ALLOWED;
			}
			else if (sType=='S')
			{
				min_phi_allowed = MIN_PHI_BETA_ALLOWED;
				max_phi_allowed = MAX_PHI_BETA_ALLOWED;
				min_psi_allowed = MIN_PSI_BETA_ALLOWED;
				max_psi_allowed = MAX_PSI_BETA_ALLOWED;
			}
			else
			{
				min_phi_allowed = MIN_PHI_LOOP_ALLOWED;
				max_phi_allowed = MAX_PHI_LOOP_ALLOWED;
				min_psi_allowed = MIN_PSI_LOOP_ALLOWED;
				max_psi_allowed = MAX_PSI_LOOP_ALLOWED;
			}

			//Working on phi .... phi of the next AA...
			p2 = mStructure.getAtomCoordinate(cAA+1, " N  ");
			p3 = mStructure.getAtomCoordinate(cAA+1, " CA ");

			torsionDelta = getMinAngle(mPoints, p3, p2, tPoints, pCntr, pCntr+1, nRMSD);
			curTorsion = mStructure.getTorsion("phi", cAA+1);
			nTorsion = curTorsion - torsionDelta;				//new torsion after you rotate

			//cout<<"  curTorsion phi= "<<curTorsion<<"  deltaAngle= "<<torsionDelta<<" nTorsion= "<<nTorsion<<" nRMSD= "<<nRMSD<<endl;
			if (validTorsion(nTorsion, sType, "phi")){// &&
				//collisionFree(mStructure, emptyTrace, nAAperSegment*(pCntr-1)+1, MIN(nAAperSegment*(pCntr+1) - 2, mStructure.numOfAA()-3), cAA+1, "phi", torsionDelta)){
				//rotate
				mStructure.rotate(cAA+1, mStructure.numOfAA()-1, mStructure.getAtomIndx(cAA+1, " CA ")+1, p2, p3, torsionDelta);
				//udate moving points....because of rotation..we need to update these points everytime
				updatePoints(mPoints, p2,p3, pCntr, mPoints.size()-1, torsionDelta);

			}
			// you can use else section in case the angle is not valid then you rotate the maximum valid angle or a specific angle
			// for example if the new angle is out of valid range, suppose valid range is [-40, 30], then you can rotate so the
			// new torsion is -40 or 30. or you can rotate certain angle, for example 10 degrees or 20 degrees so you still in the range
			// if you don't like to do so...just comment else section

			else
			{
				//cout<<"  NOT VALID PHI.."<<endl;
				//case #1 : you want to rotate so you are on the boundary of the valid range (you rotate the maximum you can)
				// first determine which boundary...which boundary is closer
				float delta1, delta2;

				delta1 = fabs(nTorsion - min_phi_allowed);
				if (delta1 > 180)
					delta1 = 360 - delta1;
				delta2 = fabs(nTorsion - max_phi_allowed);
				if (delta2 > 180)
					delta2 = 360 - delta2;

				//cout<<"  nTorsion= "<<nTorsion<<" delta1= "<<delta1<<" delta2= "<<delta2<<endl;

				if (delta1 <= delta2){ //&&
					//collisionFree(mStructure, emptyTrace, nAAperSegment*(pCntr-1)+1, MIN(nAAperSegment*(pCntr+1) - 2, mStructure.numOfAA()-3), cAA+1, "phi", torsionDelta)){
					//rotate
					mStructure.rotate(cAA+1, mStructure.numOfAA()-1, mStructure.getAtomIndx(cAA+1, " CA ")+1, p2, p3, curTorsion - min_phi_allowed);
					//udate moving points....because of rotation..we need to update these points everytime
					updatePoints(mPoints, p2,p3, pCntr, mPoints.size()-1, curTorsion - min_phi_allowed);
				}
				else
				{
					//if (collisionFree(mStructure, emptyTrace, nAAperSegment*(pCntr-1)+1, MIN(nAAperSegment*(pCntr+1) - 2, mStructure.numOfAA()-3), cAA+1, "phi", torsionDelta)){
						//rotate
						mStructure.rotate(cAA+1, mStructure.numOfAA()-1, mStructure.getAtomIndx(cAA+1, " CA ")+1, p2, p3, curTorsion - max_phi_allowed);
						//udate moving points....because of rotation..we need to update these points everytime
						updatePoints(mPoints, p2,p3, pCntr, mPoints.size()-1, curTorsion - max_phi_allowed);
					//}
				}

				//cout<<"  PHI now= "<<mStructure.getTorsion("phi", cAA+1)<<endl;
			}
			//cout<<"before phi: mPoints["<<pCntr<<"].x= "<<mPoints[pCntr].x<<" y= "<<mPoints[pCntr].y<<" z= "<<mPoints[pCntr].z<<endl;

/*
			Protein tmp;

			for (int k=0; k<mPoints.size(); k++)
			{
				AminoAcid tmpAA;
				Atom tmpAtom;
				tmpAtom.coord = mPoints[k];
				tmpAtom.name = " CA ";
				tmpAA.atoms.push_back(tmpAtom);
				tmpAA.num = k+1;

				//cout<<"working on AA# "<<i+1<<" : "<<lPoints[i].x<<" "<<lPoints[i].y<<" "<<lPoints[i].z<<endl;
				tmpAA.chr3 = "GLY";
				tmpAA.num = k+1;
				tmpAA.chain = "A";
				tmp.AAs.push_back(tmpAA);
			}

			tmp.writePDB("axisFBCCD.pdb", 1, tmp.numOfAA());
			mStructure.writePDB("inFBCCD.pdb", 1, mStructure.numOfAA());
			getchar();
*/



			itCntr++;

			if ((cAA >= nAAperSegment * (pCntr+1) - 2) || (cAA >= mStructure.numOfAA()-3))
				cAA = nAAperSegment*(pCntr-1);

			cAA += 1;
		}
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//this  function was mainly implemented to model a loop....
//the target and movable points should be determined in advance and saved in tPoints and mPoints.
void loopFBCCD(Protein &mStructure,			    //the structure we are going to change so mPoints (which are some points in this mStructure) move to be as close as possible to tPoints (target points)
				vector<Coordinate> &mPoints,	//moving points either on mStructure or some other reference points (such as center of AA)
				vector<Coordinate> tPoints,	    //target points (Coordinates)
				vector<short> corrAA,			//correspondant AA for each point in mPoint
				vector<float> corrRMSD,			//the desired correspondant RMSD for each point
				//float RMSDthr,				    //RMSD cut off threshhold
				int tIteration)				    //number of iteration cut off
{


    Coordinate p2, p3;
	int itCntr = 0,
		i,
		pCntr,
		cAA;					//current AA we are working on
	float RMSDthr,			    //RMSD cut off threshhold
		  nRMSD,				//new expected RMSD
		  torsionDelta,			//the angle we should rotate
		  curTorsion,			//current torsion angle
		  nTorsion;				//new (the expected) torsion after rotation

    //variable to get the range of phi and psi angles for a particular AA in the structure
	float	min_phi_allowed,
			max_phi_allowed,
			min_psi_allowed,
			max_psi_allowed;
    char sType;                 //the SS type of an AA

	//for greedy
	int maxSteps = 0.1*tIteration;
	float **rmsd = new float *[mStructure.numOfAA()];			//save the minimum rmsd for each angle..0 psi and 1 phi
	float **torsion = new float *[mStructure.numOfAA()];
	float **newTorsion = new float *[mStructure.numOfAA()];
	//initialize arrays
	for (i=0; i<mStructure.numOfAA();i++){
		rmsd[i] = new float [2];
		torsion[i] = new float [2];
		newTorsion[i] = new float [2];
		rmsd[i][0] = rmsd[i][1] = torsion[i][0] = torsion[i][1] = newTorsion[i][0] = newTorsion[i][1] = 9999.9;
	}

	//cout<<"#mPoints= "<<mPoints.size()<<" #tPoints= "<<tPoints.size()<<endl;
	//start from point number two...overlap two consecutive points everytime
	for (pCntr=1; pCntr<mPoints.size()-1; pCntr++)
	{
		RMSDthr = (corrRMSD[pCntr] + corrRMSD[pCntr+1])/2;

		nRMSD = RMSDthr + 1;
		//cRMSD = nRMSD+1;
		itCntr=0;
		//cAA = nAAperSegment*(pCntr-1)+1;		//start from the following AA after the last point u worked on
		cAA = corrAA[pCntr-1]-1;        //starting from first AA

		if (cAA>=mStructure.numOfAA()-2)
			continue;

		//run some greedy iterations at the beginning
		while (itCntr<maxSteps && nRMSD>RMSDthr){

			//cout<<itCntr+1<<" : Working on AA num = "<<mStructure.AAs[cAA].num<<" type= "<<mStructure.AAs[cAA].SStype<<" pCntr= "<<pCntr<<endl;

			//get the type of the amino acid...H, L, or S
			sType = mStructure.AAs[cAA].SStype;
            if (sType=='H')
            {
                min_phi_allowed = MIN_PHI_HLX_ALLOWED;
                max_phi_allowed = MAX_PHI_HLX_ALLOWED;
                min_psi_allowed = MIN_PSI_HLX_ALLOWED;
                max_psi_allowed = MAX_PSI_HLX_ALLOWED;
            }
            else if (sType=='S')
            {
                min_phi_allowed = MIN_PHI_BETA_ALLOWED;
                max_phi_allowed = MAX_PHI_BETA_ALLOWED;
                min_psi_allowed = MIN_PSI_BETA_ALLOWED;
                max_psi_allowed = MAX_PSI_BETA_ALLOWED;
            }
            else
            {
                min_phi_allowed = MIN_PHI_LOOP_ALLOWED;
                max_phi_allowed = MAX_PHI_LOOP_ALLOWED;
                min_psi_allowed = MIN_PSI_LOOP_ALLOWED;
                max_psi_allowed = MAX_PSI_LOOP_ALLOWED;
            }
			//Working on Psi
			p2 = mStructure.getAtomCoordinate(cAA, " CA ");
			p3 = mStructure.getAtomCoordinate(cAA, " C  ");

			if (cAA < corrAA[pCntr])
				torsionDelta = getMinAngle(mPoints, p3, p2, tPoints, pCntr, pCntr+1, nRMSD);
			else{
				torsionDelta = getMinAngle(mPoints, p3, p2, tPoints, pCntr+1, pCntr+1, nRMSD);
				//add the RMSD of the current point pCntr
				nRMSD += getDistance(mPoints[pCntr], tPoints[pCntr]);
				nRMSD /= 2;
			}

			curTorsion = mStructure.getTorsion("psi", cAA);
			nTorsion = curTorsion - torsionDelta;		//new torsion after you rotate
			if (nTorsion<-180)
				nTorsion = 360+nTorsion;
			else if (nTorsion>180)
				nTorsion = 360-nTorsion;

			//cout<<"  curTorsion psi= "<<curTorsion<<"  deltaAngle= "<<torsionDelta<<" nTorsion= "<<nTorsion<<" nRMSD= "<<nRMSD<<endl;
			//save the expected rmsd and the angle delta
			rmsd[cAA][0] = nRMSD;
			torsion[cAA][0] = torsionDelta;
			newTorsion[cAA][0] = nTorsion;

            //check the SS type of the next AA...cAA+1
			sType = mStructure.AAs[cAA+1].SStype;
            if (sType=='H')
            {
                min_phi_allowed = MIN_PHI_HLX_ALLOWED;
                max_phi_allowed = MAX_PHI_HLX_ALLOWED;
                min_psi_allowed = MIN_PSI_HLX_ALLOWED;
                max_psi_allowed = MAX_PSI_HLX_ALLOWED;
            }
            else if (sType=='S')
            {
                min_phi_allowed = MIN_PHI_BETA_ALLOWED;
                max_phi_allowed = MAX_PHI_BETA_ALLOWED;
                min_psi_allowed = MIN_PSI_BETA_ALLOWED;
                max_psi_allowed = MAX_PSI_BETA_ALLOWED;
            }
            else
            {
                min_phi_allowed = MIN_PHI_LOOP_ALLOWED;
                max_phi_allowed = MAX_PHI_LOOP_ALLOWED;
                min_psi_allowed = MIN_PSI_LOOP_ALLOWED;
                max_psi_allowed = MAX_PSI_LOOP_ALLOWED;
            }
			//Working on phi .... phi of the next AA...
			p2 = mStructure.getAtomCoordinate(cAA+1, " N  ");
			p3 = mStructure.getAtomCoordinate(cAA+1, " CA ");

			if (cAA+1<corrAA[pCntr])
				torsionDelta = getMinAngle(mPoints, p3, p2, tPoints, pCntr, pCntr+1, nRMSD);
			else{
				torsionDelta = getMinAngle(mPoints, p3, p2, tPoints, pCntr+1, pCntr+1, nRMSD);

				//add the RMSD of the first point
				nRMSD += getDistance(mPoints[pCntr], tPoints[pCntr]);
				nRMSD /=2;
			}

			curTorsion = mStructure.getTorsion("phi", cAA+1);
			nTorsion = curTorsion - torsionDelta;				//new torsion after you rotate

			if (nTorsion<-180)
				nTorsion = 360+nTorsion;
			else if (nTorsion>180)
				nTorsion = 360-nTorsion;

			//cout<<"  curTorsion phi= "<<curTorsion<<"  deltaAngle= "<<torsionDelta<<" nTorsion= "<<nTorsion<<" nRMSD= "<<nRMSD<<endl;
			//save the expected rmsd and the angle delta
			rmsd[cAA+1][1] = nRMSD;
			torsion[cAA+1][1] = torsionDelta;
			newTorsion[cAA+1][1] = nTorsion;

			if (cAA >= corrAA[pCntr+1]-2 || cAA>=mStructure.numOfAA()-3){
                cAA = corrAA[pCntr-1]-2;        //starting from first corresponding AA for this point

				//cout<<" =============== End of the cycle ==========="<<endl;
				//rotate the angle gives the minimum RMSD
				short aaIndx=-1;
				float minRMSD = 9999.9;
				string torsionType;
				for (i=0; i<mStructure.numOfAA(); i++){
					//cout<<mStructure.AAs[i].num<<" phi "<<validTorsion(newTorsion[i][1], mStructure.AAs[i].SStype, "phi")<<" psi "<<validTorsion(newTorsion[i][0], mStructure.AAs[i].SStype, "psi")<<endl;
					if (rmsd[i][0]<minRMSD && validTorsion(newTorsion[i][0], mStructure.AAs[i].SStype, "psi")){
						aaIndx = i;
						minRMSD = rmsd[i][0];
						torsionType = "psi";
					}
					if (rmsd[i][1]<minRMSD && validTorsion(newTorsion[i][1], mStructure.AAs[i].SStype, "phi")){
						aaIndx = i;
						minRMSD = rmsd[i][1];
						torsionType = "phi";
					}
				}
				if (aaIndx != -1)		//mStructure.getTorsion("psi", cAA) - torsionDelta : is the torsion angle after we rotate
				{
					//cout<<"minimum RMSD was applied on AA num "<<mStructure.AAs[aaIndx].num<<" RMSD= "<<minRMSD<<endl;

					if (torsionType == "psi"){
						//Working on Psi
						p2 = mStructure.getAtomCoordinate(aaIndx, " CA ");
						p3 = mStructure.getAtomCoordinate(aaIndx, " C  ");
						//rotate
						mStructure.rotate(aaIndx, mStructure.numOfAA()-1, mStructure.getAtomIndx(aaIndx, " C  ")+1, p2, p3, torsion[aaIndx][0]);
						//udate moving points....because of rotation..we need to update these points everytime
						if (aaIndx < corrAA[pCntr]){
							updatePoints(mPoints, p2,p3, pCntr, mPoints.size()-1, torsion[aaIndx][0]);
						}
						else{
							updatePoints(mPoints, p2,p3, pCntr+1, mPoints.size()-1, torsion[aaIndx][0]);
						}

					}
					else{
						//Working on phi .... phi of the next AA...
						p2 = mStructure.getAtomCoordinate(aaIndx, " N  ");
						p3 = mStructure.getAtomCoordinate(aaIndx, " CA ");
						mStructure.rotate(aaIndx, mStructure.numOfAA()-1, mStructure.getAtomIndx(aaIndx, " CA ")+1, p2, p3, torsion[aaIndx][1]);
						//udate moving points....because of rotation..we need to update these points everytime
						if (aaIndx < corrAA[pCntr]){
							updatePoints(mPoints, p2,p3, pCntr, mPoints.size()-1, torsion[aaIndx][1]);
						}
						else{
							updatePoints(mPoints, p2,p3, pCntr+1, mPoints.size()-1, torsion[aaIndx][1]);
						}
					}
					//updatePoints(mPoints, p2,p3, torsionDelta);
				}
				for (i=0; i<mStructure.numOfAA(); i++){
					rmsd[i][0] = rmsd[i][1] = torsion[i][0] = torsion[i][1] = newTorsion[i][0] = newTorsion[i][1] = 9999.9;
				}
				//getchar();
				nRMSD = minRMSD;
			}
			cAA += 1;
			itCntr++;
		}

		cout<<"minRMSD= "<<nRMSD<<endl;
		while ((nRMSD > RMSDthr) && (itCntr < tIteration))
		{
			//cout<<itCntr+1<<" : Working on AA indx = "<<cAA<<" pCntr= "<<pCntr<<endl;

		//get the type of the amino acid...H, L, or S
			sType = mStructure.AAs[cAA].SStype;
            if (sType=='H')
            {
                min_phi_allowed = MIN_PHI_HLX_ALLOWED;
                max_phi_allowed = MAX_PHI_HLX_ALLOWED;
                min_psi_allowed = MIN_PSI_HLX_ALLOWED;
                max_psi_allowed = MAX_PSI_HLX_ALLOWED;
            }
            else if (sType=='S')
            {
                min_phi_allowed = MIN_PHI_BETA_ALLOWED;
                max_phi_allowed = MAX_PHI_BETA_ALLOWED;
                min_psi_allowed = MIN_PSI_BETA_ALLOWED;
                max_psi_allowed = MAX_PSI_BETA_ALLOWED;
            }
            else
            {
                min_phi_allowed = MIN_PHI_LOOP_ALLOWED;
                max_phi_allowed = MAX_PHI_LOOP_ALLOWED;
                min_psi_allowed = MIN_PSI_LOOP_ALLOWED;
                max_psi_allowed = MAX_PSI_LOOP_ALLOWED;
            }

			//Working on Psi
			p2 = mStructure.getAtomCoordinate(cAA, " CA ");
			p3 = mStructure.getAtomCoordinate(cAA, " C  ");

			if (cAA < corrAA[pCntr])
				torsionDelta = getMinAngle(mPoints, p3, p2, tPoints, pCntr, pCntr+1, nRMSD);
			else{
				torsionDelta = getMinAngle(mPoints, p3, p2, tPoints, pCntr+1, pCntr+1, nRMSD);

				//add the RMSD of the first point
				nRMSD += getDistance(mPoints[pCntr], tPoints[pCntr]);
				nRMSD /=2;
			}

			//torsionDelta = getMinAngle(mPoints, p3, p2, tPoints, pCntr, mPoints.size ()-1, nRMSD);
			//torsionDelta = getMinAngle(mPoints, p3, p2, tPoints, nRMSD);
			curTorsion = mStructure.getTorsion("psi", cAA);
			nTorsion = curTorsion - torsionDelta;		//new torsion after you rotate
			if (nTorsion<-180)
				nTorsion = 360+nTorsion;
			else if (nTorsion>180)
				nTorsion = 360-nTorsion;

			//cout<<"  curTorsion psi= "<<curTorsion<<"  deltaAngle= "<<torsionDelta<<" nTorsion= "<<nTorsion<<" nRMSD= "<<nRMSD<<endl;
			//cout<<"before psi: mPoints["<<pCntr<<"].x= "<<mPoints[pCntr].x<<" y= "<<mPoints[pCntr].y<<" z= "<<mPoints[pCntr].z<<endl;

			if (validTorsion(nTorsion, sType, "psi"))		//mStructure.getTorsion("psi", cAA) - torsionDelta : is the torsion angle after we rotate
			{
				//rotate
				mStructure.rotate(cAA, mStructure.numOfAA()-1, mStructure.getAtomIndx(cAA, " C  ")+1, p2, p3, torsionDelta);

				//udate moving points....because of rotation..we need to update these points everytime
				if (cAA < corrAA[pCntr]){
					updatePoints(mPoints, p2,p3, pCntr, mPoints.size()-1, torsionDelta);

				}
				else{
					updatePoints(mPoints, p2,p3, pCntr+1, mPoints.size()-1, torsionDelta);
				}

				//updatePoints(mPoints, p2,p3, torsionDelta);
			}
			/*
			// you can use else section in case the angle is not valid then you rotate the maximum valid angle or a specific angle
			// for example if the new angle is out of valid range, suppose valid range is [-40, 30], then you can rotate so the
			// new torsion is -40 or 30. or you can rotate certain angle, for example 10 degrees or 20 degrees so you still in the range
			// if you don't like to do so...just comment else section
			else
			{
				//cout<<"  NOT VALID PSI.."<<endl;
				//case #1 : you want to rotate so you are on the boundary of the valid range (you rotate the maximum you can)
				// first determine which boundary...which boundary is closer
				float delta1, delta2;

				delta1 = fabs(nTorsion - min_psi_allowed);
				if (delta1 > 180)
					delta1 = 360 - delta1;
				delta2 = fabs(nTorsion - max_psi_allowed);
				if (delta2 > 180)
					delta2 = 360 - delta2;

				//cout<<"  nTorsion= "<<nTorsion<<" delta1= "<<delta1<<" delta2= "<<delta2<<endl;

				if (delta1 <= delta2)
				{
					//rotate to MIN boundary
					mStructure.rotate(cAA, mStructure.numOfAA()-1, mStructure.getAtomIndx(cAA, " C  ")+1, p2, p3, curTorsion - min_psi_allowed);
					//udate moving points....because of rotation..we need to update these points everytime
					updatePoints(mPoints, p2,p3, pCntr, mPoints.size()-1, curTorsion-min_psi_allowed);
					//updatePoints(mPoints, p2,p3, curTorsion-min_psi_allowed);
				}
				else
				{
					//rotate to Max boundary
					mStructure.rotate(cAA, mStructure.numOfAA()-1, mStructure.getAtomIndx(cAA, " C  ")+1, p2, p3, curTorsion - max_psi_allowed);
					//udate moving points....because of rotation..we need to update these points everytime
					updatePoints(mPoints, p2,p3, pCntr, mPoints.size()-1, curTorsion - max_psi_allowed);
					//updatePoints(mPoints, p2,p3, curTorsion - max_psi_allowed);
				}

				//cout<<"  PSI now= "<<mStructure.getTorsion("psi", cAA)<<endl;
			}
			*/
			//cout<<"after psi: mPoints["<<pCntr<<"].x= "<<mPoints[pCntr].x<<" y= "<<mPoints[pCntr].y<<" z= "<<mPoints[pCntr].z<<endl;



			//check the RMSD...
			//because FBCCD is not greedy..so it is not guranteed that any rotation will lead to better RMSD
			//so, don't rotate any more if you reach the desired RMSD
			if (nRMSD <= RMSDthr)		//in case you reach the RMSD you desire...
				break;

            //check the SS type of the next AA...cAA+1
			sType = mStructure.AAs[cAA+1].SStype;
            if (sType=='H')
            {
                min_phi_allowed = MIN_PHI_HLX_ALLOWED;
                max_phi_allowed = MAX_PHI_HLX_ALLOWED;
                min_psi_allowed = MIN_PSI_HLX_ALLOWED;
                max_psi_allowed = MAX_PSI_HLX_ALLOWED;
            }
            else if (sType=='S')
            {
                min_phi_allowed = MIN_PHI_BETA_ALLOWED;
                max_phi_allowed = MAX_PHI_BETA_ALLOWED;
                min_psi_allowed = MIN_PSI_BETA_ALLOWED;
                max_psi_allowed = MAX_PSI_BETA_ALLOWED;
            }
            else
            {
                min_phi_allowed = MIN_PHI_LOOP_ALLOWED;
                max_phi_allowed = MAX_PHI_LOOP_ALLOWED;
                min_psi_allowed = MIN_PSI_LOOP_ALLOWED;
                max_psi_allowed = MAX_PSI_LOOP_ALLOWED;
            }
			//Working on phi .... phi of the next AA...
			p2 = mStructure.getAtomCoordinate(cAA+1, " N  ");
			p3 = mStructure.getAtomCoordinate(cAA+1, " CA ");

			if (cAA+1<corrAA[pCntr])
				torsionDelta = getMinAngle(mPoints, p3, p2, tPoints, pCntr, pCntr+1, nRMSD);
			else{
				torsionDelta = getMinAngle(mPoints, p3, p2, tPoints, pCntr+1, pCntr+1, nRMSD);
				//add the RMSD of the first point
				nRMSD += getDistance(mPoints[pCntr], tPoints[pCntr]);
				nRMSD /=2;
			}

			//torsionDelta = getMinAngle(mPoints, p3, p2, tPoints, pCntr, mPoints.size ()-1, nRMSD);
			//torsionDelta = getMinAngle(mPoints, p3, p2, tPoints, nRMSD);
			curTorsion = mStructure.getTorsion("phi", cAA+1);
			nTorsion = curTorsion - torsionDelta;				//new torsion after you rotate
			if (nTorsion<-180)
				nTorsion = 360+nTorsion;
			else if (nTorsion>180)
				nTorsion = 360-nTorsion;

			//cout<<"  curTorsion phi= "<<curTorsion<<"  deltaAngle= "<<torsionDelta<<" nTorsion= "<<nTorsion<<" nRMSD= "<<nRMSD<<endl;
			if (validTorsion(nTorsion, sType, "phi"))	//mStructure.getTorsion("phi", cAA) - torsionDelta : is the torsion angle after we rotate
			{

				//rotate
				mStructure.rotate(cAA+1, mStructure.numOfAA()-1, mStructure.getAtomIndx(cAA+1, " CA ")+1, p2, p3, torsionDelta);

				//udate moving points....because of rotation..we need to update these points everytime
				if (cAA+1 < corrAA[pCntr]){
					updatePoints(mPoints, p2,p3, pCntr, mPoints.size()-1, torsionDelta);
				}
				else{
					updatePoints(mPoints, p2,p3, pCntr+1, mPoints.size()-1, torsionDelta);
				}
				//updatePoints(mPoints, p2,p3, torsionDelta);
			}


			// you can use else section in case the angle is not valid then you rotate the maximum valid angle or a specific angle
			// for example if the new angle is out of valid range, suppose valid range is [-40, 30], then you can rotate so the
			// new torsion is -40 or 30. or you can rotate certain angle, for example 10 degrees or 20 degrees so you still in the range
			// if you don't like to do so...just comment else section
			/*
			else
			{
				//cout<<"  NOT VALID PHI.."<<endl;
				//case #1 : you want to rotate so you are on the boundary of the valid range (you rotate the maximum you can)
				// first determine which boundary...which boundary is closer
				float delta1, delta2;

				delta1 = fabs(nTorsion - min_phi_allowed);
				if (delta1 > 180)
					delta1 = 360 - delta1;
				delta2 = fabs(nTorsion - max_phi_allowed);
				if (delta2 > 180)
					delta2 = 360 - delta2;

				//cout<<"  nTorsion= "<<nTorsion<<" delta1= "<<delta1<<" delta2= "<<delta2<<endl;

				if (delta1 <= delta2)
				{
					//rotate
					mStructure.rotate(cAA+1, mStructure.numOfAA()-1, mStructure.getAtomIndx(cAA+1, " CA ")+1, p2, p3, curTorsion - min_phi_allowed);
					//udate moving points....because of rotation..we need to update these points everytime
					updatePoints(mPoints, p2,p3, pCntr, mPoints.size()-1, curTorsion - min_phi_allowed);
					//updatePoints(mPoints, p2,p3, curTorsion - min_phi_allowed);
				}
				else
				{
					//rotate
					mStructure.rotate(cAA+1, mStructure.numOfAA()-1, mStructure.getAtomIndx(cAA+1, " CA ")+1, p2, p3, curTorsion - max_phi_allowed);
					//udate moving points....because of rotation..we need to update these points everytime
					updatePoints(mPoints, p2,p3, pCntr, mPoints.size()-1, curTorsion - max_phi_allowed);
					//updatePoints(mPoints, p2,p3, curTorsion - max_phi_allowed);
				}

				//cout<<"  PHI now= "<<mStructure.getTorsion("phi", cAA+1)<<endl;
			}
			//cout<<"before phi: mPoints["<<pCntr<<"].x= "<<mPoints[pCntr].x<<" y= "<<mPoints[pCntr].y<<" z= "<<mPoints[pCntr].z<<endl;
			*/

/*
			cout<<"dist b/w Pnts "<<pCntr<<" and "<<pCntr+1<<" "<<getDistance(mPoints[pCntr-1], mPoints[pCntr])<<endl;
			cout<<"dist b/w Pnts "<<pCntr+1<<" and "<<pCntr+2<<" "<<getDistance(mPoints[pCntr], mPoints[pCntr+1])<<endl;

			Protein tmp;

			for (int k=0; k<mPoints.size(); k++)
			{
				AminoAcid tmpAA;
				Atom tmpAtom;
				tmpAtom.coord = mPoints[k];
				tmpAtom.name = " CA ";
				tmpAA.atoms.push_back(tmpAtom);
				tmpAA.num = k+1;

				//cout<<"working on AA# "<<i+1<<" : "<<lPoints[i].x<<" "<<lPoints[i].y<<" "<<lPoints[i].z<<endl;
				tmpAA.chr3 = "GLY";
				tmpAA.num = k+1;
				tmpAA.chain = "A";
				tmp.AAs.push_back(tmpAA);
			}

			tmp.writePDB("axisFBCCD.pdb", 1, tmp.numOfAA());
			mStructure.writePDB("inFBCCD.pdb", 1, mStructure.numOfAA());
			getchar();
*/



			itCntr++;

			//if ((cAA >= nAAperSegment * (pCntr+1) - 2) || (cAA >= mStructure.numOfAA()-3))          //reached the end
				//cAA = nAAperSegment*(pCntr-1);
			if (cAA >= corrAA[pCntr+1]-2 || cAA>=mStructure.numOfAA()-3)
                cAA = corrAA[pCntr-1]+1;        //starting from first corresponding AA for this point
			cAA += 1;
			//getchar();
		}
		cout<<"pnt "<<pCntr+1<<" and "<<pCntr+2<<"  minRMSD= "<<nRMSD<<" nIterations= "<<itCntr<<endl;
	}

	delete [] *rmsd;
	delete [] rmsd;
	delete [] *torsion;
	delete [] torsion;
	delete [] *newTorsion;
	delete [] newTorsion;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//greedy FBCD...it assumes that the set of amino acids you wanna overlap are prior the target points, if the target points on the mStructure or not as well
//the algorithm works on the torsion angles as the following
//between [psi-startAA-1 ... phi-endAA-1] in case of forward
//between [psi-endAA-1 .... phi-startAA-1] in case of reverse
float greedyFBCCD(Protein &mStructure,	//the structure we are going to change so mPoints (which are some points in this mStructure) move to be as close as possible to tPoints (target points)
					vector<Coordinate> &mPoints,		//moving points either on mStructure or some other reference points (such as center of AA)
					vector<Coordinate> tPoints,	    //target points (Coordinates)
					vector<vector<Coordinate> > avoidTrace,			//the trace we need to avoid when connect the two portions
					short startAA,			//the index of aa to start rotating angles
					short endAA,			//the index of aa to end
					short collStartAA,		//the index of start AA to check collision from
					short collEndAA,        //the end index of the last AA to check for collision
					short direction,		//the direction of walk....-1 reverse...1 forward
					float RMSDthr,			//RMSD cut off threshhold
					int gIterations,			//number of cycles cut off
					bool dontCheckColl= true){	//a flag to chekc for collision or not


	if (direction != 1 && direction != -1){
		errMsg("FBCD.h", "greedyFBCD", "The direction should be either 1 (Forward walk) or -1 (backward walk).");
		return 9999.9;
	}

	int sAA, eAA;		//the index of start and end AA

	if (direction == 1){
		if (startAA<1 || endAA>=mStructure.numOfAA()){
			errMsg("FBCD.h", "greedyFBCD", "The given two indeces ("+toString(startAA)+" , "+toString(endAA)+ "are not allowed...should be b/w [1 ,"+toString(mStructure.numOfAA()-2)+"].");
			return 9999.9;
		}
		sAA = startAA-1;
		eAA = endAA-2;
	}
	else{
		if (endAA<1 || startAA>=mStructure.numOfAA()){
			errMsg("FBCD.h", "greedyFBCD", "The given two indeces ("+toString(startAA)+" , "+toString(endAA)+ "are not allowed...should be b/w ["+toString(mStructure.numOfAA()-2)+",1].");
			return 9999.9;
		}
		sAA = startAA-2;
		eAA = endAA-1;
	}

	Coordinate p2,p3;				//temp 3d Coordinate points to hold the atoms of psi or phi
	float nRMSD = 9999.9,		//the expected RMSD for each angle
		  torsionDelta,			//the angle we should rotate
		  curTorsion,			//current torsion angle
		  nTorsion;				//new (the expected) torsion after rotation
	int itCntr=0,				//iteration counter
		cAA,				//counter over AAs
		mPointsLastIndx=mPoints.size ()-1,		//last index in mPoints
		i;
	char sType;						//the type of current AA (H, S, L)
    //variable to get the range of phi and psi angles for a particular AA in the structure
	float	min_phi_allowed,
			max_phi_allowed,
			min_psi_allowed,
			max_psi_allowed;

	//for greedy
	float **rmsd = new float *[mStructure.numOfAA()];			//save the minimum rmsd for each angle..0 psi and 1 phi
	float **torsion = new float *[mStructure.numOfAA()];
	float **newTorsion = new float *[mStructure.numOfAA()];
	//initialize arrays
	for (i=0; i<mStructure.numOfAA();i++){
		rmsd[i] = new float [2];
		torsion[i] = new float [2];
		newTorsion[i] = new float [2];
		rmsd[i][0] = rmsd[i][1] = torsion[i][0] = torsion[i][1] = newTorsion[i][0] = newTorsion[i][1] = 9999.9;
	}

	//run some greedy iterations at the beginning
	while (nRMSD>RMSDthr && itCntr<gIterations){
		cAA=sAA;
		while (cAA!=eAA+direction){

			//cout<<itCntr+1<<" : Working on AA num = "<<mStructure.AAs[cAA].num<<" type= "<<mStructure.AAs[cAA].SStype<<endl;

			//get the type of the amino acid...H, L, or S
			sType = mStructure.AAs[cAA].SStype;
			if (sType=='H')
			{
				min_phi_allowed = MIN_PHI_HLX_ALLOWED;
				max_phi_allowed = MAX_PHI_HLX_ALLOWED;
				min_psi_allowed = MIN_PSI_HLX_ALLOWED;
				max_psi_allowed = MAX_PSI_HLX_ALLOWED;
			}
			else if (sType=='S')
			{
				min_phi_allowed = MIN_PHI_BETA_ALLOWED;
				max_phi_allowed = MAX_PHI_BETA_ALLOWED;
				min_psi_allowed = MIN_PSI_BETA_ALLOWED;
				max_psi_allowed = MAX_PSI_BETA_ALLOWED;
			}
			else
			{
				min_phi_allowed = MIN_PHI_LOOP_ALLOWED;
				max_phi_allowed = MAX_PHI_LOOP_ALLOWED;
				min_psi_allowed = MIN_PSI_LOOP_ALLOWED;
				max_psi_allowed = MAX_PSI_LOOP_ALLOWED;
			}
			//Working on Psi
			p2 = mStructure.getAtomCoordinate(cAA, " CA ");
			p3 = mStructure.getAtomCoordinate(cAA, " C  ");

			torsionDelta = getMinAngle(mPoints, p3, p2, tPoints, 0, mPointsLastIndx, nRMSD);

			curTorsion = mStructure.getTorsion("psi", cAA);
			nTorsion = curTorsion - torsionDelta;		//new torsion after you rotate
			if (nTorsion<-180)
				nTorsion = 360+nTorsion;
			else if (nTorsion>180)
				nTorsion = 360-nTorsion;

			//cout<<"  curTorsion psi= "<<curTorsion<<"  deltaAngle= "<<torsionDelta<<" nTorsion= "<<nTorsion<<" nRMSD= "<<nRMSD<<endl;
			//save the expected rmsd and the angle delta
			rmsd[cAA][0] = nRMSD;
			torsion[cAA][0] = torsionDelta;
			newTorsion[cAA][0] = nTorsion;

			//check the SS type of the next AA...cAA+1
			sType = mStructure.AAs[cAA+1].SStype;
			if (sType=='H')
			{
				min_phi_allowed = MIN_PHI_HLX_ALLOWED;
				max_phi_allowed = MAX_PHI_HLX_ALLOWED;
				min_psi_allowed = MIN_PSI_HLX_ALLOWED;
				max_psi_allowed = MAX_PSI_HLX_ALLOWED;
			}
			else if (sType=='S')
			{
				min_phi_allowed = MIN_PHI_BETA_ALLOWED;
				max_phi_allowed = MAX_PHI_BETA_ALLOWED;
				min_psi_allowed = MIN_PSI_BETA_ALLOWED;
				max_psi_allowed = MAX_PSI_BETA_ALLOWED;
			}
			else
			{
				min_phi_allowed = MIN_PHI_LOOP_ALLOWED;
				max_phi_allowed = MAX_PHI_LOOP_ALLOWED;
				min_psi_allowed = MIN_PSI_LOOP_ALLOWED;
				max_psi_allowed = MAX_PSI_LOOP_ALLOWED;
			}
			//Working on phi .... phi of the next AA...
			p2 = mStructure.getAtomCoordinate(cAA+1, " N  ");
			p3 = mStructure.getAtomCoordinate(cAA+1, " CA ");

			torsionDelta = getMinAngle(mPoints, p3, p2, tPoints, 0, mPointsLastIndx, nRMSD);

			curTorsion = mStructure.getTorsion("phi", cAA+1);
			nTorsion = curTorsion - torsionDelta;				//new torsion after you rotate

			if (nTorsion<-180)
				nTorsion = 360+nTorsion;
			else if (nTorsion>180)
				nTorsion = 360-nTorsion;

			//cout<<"  curTorsion phi= "<<curTorsion<<"  deltaAngle= "<<torsionDelta<<" nTorsion= "<<nTorsion<<" nRMSD= "<<nRMSD<<endl;
			//save the expected rmsd and the angle delta
			rmsd[cAA+1][1] = nRMSD;
			torsion[cAA+1][1] = torsionDelta;
			newTorsion[cAA+1][1] = nTorsion;

			cAA += direction;
		}

		//cout<<" =============== End of the cycle ==========="<<endl;
		//rotate the angle gives the minimum RMSD
		short aaIndx=-1;
		float minRMSD = 9999.9;
		string torsionType;
		for (i=0; i<mStructure.numOfAA(); i++){
			//cout<<mStructure.AAs[i].num<<" phi "<<validTorsion(newTorsion[i][1], mStructure.AAs[i].SStype, "phi")<<" psi "<<validTorsion(newTorsion[i][0], mStructure.AAs[i].SStype, "psi")<<endl;
			if (rmsd[i][0]<minRMSD && validTorsion(newTorsion[i][0], mStructure.AAs[i].SStype, "psi")){
				aaIndx = i;
				minRMSD = rmsd[i][0];
				torsionType = "psi";
			}
			if (rmsd[i][1]<minRMSD && validTorsion(newTorsion[i][1], mStructure.AAs[i].SStype, "phi")){
				aaIndx = i;
				minRMSD = rmsd[i][1];
				torsionType = "phi";
			}
		}
		if (aaIndx != -1)		//mStructure.getTorsion("psi", cAA) - torsionDelta : is the torsion angle after we rotate
		{
			//cout<<"AA num "<<mStructure.AAs[aaIndx].num<<" nRMSD= "<<minRMSD;

			if (torsionType == "psi"){

				if (dontCheckColl || collisionFree(mStructure, avoidTrace, collStartAA, collEndAA, aaIndx, "psi", torsion[aaIndx][0])){
					//cout<<" delta= "<<torsion[aaIndx][0]<<" psi nTorsion= "<<newTorsion[aaIndx][0]<<endl;
					//Working on Psi
					p2 = mStructure.getAtomCoordinate(aaIndx, " CA ");
					p3 = mStructure.getAtomCoordinate(aaIndx, " C  ");
					//rotate
					mStructure.rotate(aaIndx, mStructure.numOfAA()-1, mStructure.getAtomIndx(aaIndx, " C  ")+1, p2, p3, torsion[aaIndx][0]);
					//udate moving points....because of rotation..we need to update these points everytime
					updatePoints(mPoints, p2,p3, 0, mPointsLastIndx, torsion[aaIndx][0]);
				}
			}
			else{
				if (dontCheckColl || collisionFree(mStructure, avoidTrace, collStartAA, collEndAA, aaIndx, "phi", torsion[aaIndx][1])){
					//cout<<" delta= "<<torsion[aaIndx][1]<<" phi nTorsion= "<<newTorsion[aaIndx][1]<<endl;
					//Working on phi .... phi of the next AA...
					p2 = mStructure.getAtomCoordinate(aaIndx, " N  ");
					p3 = mStructure.getAtomCoordinate(aaIndx, " CA ");
					mStructure.rotate(aaIndx, mStructure.numOfAA()-1, mStructure.getAtomIndx(aaIndx, " CA ")+1, p2, p3, torsion[aaIndx][1]);
					//udate moving points....because of rotation..we need to update these points everytime
					updatePoints(mPoints, p2,p3, 0, mPointsLastIndx, torsion[aaIndx][1]);
				}
			}
			//mStructure.writePDB("inFBCD_"+toString(endAA)+"_"+toString(itCntr+1)+".pdb", 1, mStructure.numOfAA());
			nRMSD = minRMSD;
		}
		for (i=0; i<mStructure.numOfAA(); i++){
			rmsd[i][0] = rmsd[i][1] = torsion[i][0] = torsion[i][1] = newTorsion[i][0] = newTorsion[i][1] = 9999.9;
		}
		//getchar();
		itCntr++;
	}
	delete [] *rmsd;
	delete [] rmsd;
	delete [] *torsion;
	delete [] torsion;
	delete [] *newTorsion;
	delete [] newTorsion;

	return nRMSD;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//this  function was mainly implemented to model a loop....
//the target and movable points should be determined in advance and saved in tPoint and mPoint.
//the algorithm assumes that no target point should be between startAA and endAA
//greedy FBCD...it assumes that the set of amino acids you wanna overlap are prior the target points, if the target points on the mStructure or not as well
//the algorithm works on the torsion angles as the following
//between [psi-startAA-1 ... phi-endAA-1] in case of forward
//between [psi-endAA-1 .... phi-startAA-1] in case of reverse
float FBCCDnPoints(Protein &mStructure,	//the structure we are going to change so mPoints (which are some points in this mStructure) move to be as close as possible to tPoints (target points)
					vector<Coordinate> &mPoints,	//moving points either on mStructure or some other reference points (such as center of AA)
					vector<Coordinate> tPoints,    //target points (Coordinates)
					vector<vector<Coordinate> > avoidTrace,		//the trace of already built structure to avoid collision
					short startAA,			//the index of aa to start rotating angles
					short endAA,			//the index of the last aa
					short collStartAA,		//the indexs of start AA to check collision
					short collEndAA,        //the index of the end AA to check the collision...always start AA should be less than start
					float RMSDthr,			//RMSD cut off threshhold
					int tIteration,			//number of cycles cut off
					short direction,		//the direction of walk
					bool relax=false,		//flag to allow relaxing the angle in case it is not valid
					bool greedy=false,		//flag to allow run greedy steps at the beginning
					float greedyPercentage=0.2,	//the percentage of greedy steps out of total number of iterations
					bool dontCheckColl= true)	//a flag to chekc for collision or not
{

	if (direction != 1 && direction != -1){
		errMsg("FBCD.h", "FBCDnPoints", "The direction should be either 1 (Forward walk) or -1 (backward walk).");
		return 9999.9;
	}

	int sAA, eAA;		//the index of start and end AA

	if (direction == 1){
		if (startAA<1 || endAA>=mStructure.numOfAA()){
			errMsg("FBCD.h", "FBCDnPoints", "The given two indeces ("+toString(startAA)+" , "+toString(endAA)+ "are not allowed...should be b/w [1 ,"+toString(mStructure.numOfAA()-2)+"].");
			return 9999.9;
		}
		sAA = startAA-1;
		eAA = endAA-2;
	}
	else{
		if (endAA<1 || startAA>=mStructure.numOfAA()){
			errMsg("FBCD.h", "FBCDnPoints", "The given two indeces ("+toString(startAA)+" , "+toString(endAA)+ "are not allowed...should be b/w ["+toString(mStructure.numOfAA()-2)+",1].");
			return 9999.9;
		}
		sAA = startAA-2;
		eAA = endAA-1;
	}

    Coordinate p2, p3;
	int itCntr = 0,
		i,
		mPointsLastIndx=mPoints.size ()-1,		//last index in mPoints
		cAA;					//current AA we are working on
	float nRMSD=9999.9,				//new expected RMSD
			cRMSD=nRMSD+1,				//the last updated RMSD
		  torsionDelta,			//the angle we should rotate
		  curTorsion,			//current torsion angle
		  nTorsion;				//new (the expected) torsion after rotation

	Protein bestStructure = mStructure;
	float bestRMSD=cRMSD;
    //variable to get the range of phi and psi angles for a particular AA in the structure
	float	min_phi_allowed,
			max_phi_allowed,
			min_psi_allowed,
			max_psi_allowed;
    char sType;                 //the SS type of an AA
	bool rotated;				//a flag to indicate wether the structure has been rotated or not

	//cout<<"#mPoints= "<<mPoints.size()<<" #tPoints= "<<tPoints.size()<<endl;
	if (greedy){
		itCntr = greedyPercentage*tIteration;
		nRMSD = greedyFBCCD(mStructure,mPoints,tPoints, avoidTrace, startAA,endAA, collStartAA, collEndAA, direction, RMSDthr,itCntr, dontCheckColl);
		if (nRMSD < RMSDthr)
			return nRMSD;
		bestStructure = mStructure;
		bestRMSD = nRMSD;
	}

    cRMSD = nRMSD;
	while (itCntr<tIteration){
		cAA=sAA;
		while (cAA!=eAA+direction){

			rotated = false;
			//cout<<itCntr+1<<" : Working on AA num = "<<mStructure.AAs[cAA].num<<" type= "<<mStructure.AAs[cAA].SStype<<endl;

			//mStructure.writePDB("mStructureAtBeginning.pdb", 1, mStructure.numOfAA());

			//get the type of the amino acid...H, L, or S
			sType = mStructure.AAs[cAA].SStype;
            if (sType=='H')
            {
                min_phi_allowed = MIN_PHI_HLX_ALLOWED;
                max_phi_allowed = MAX_PHI_HLX_ALLOWED;
                min_psi_allowed = MIN_PSI_HLX_ALLOWED;
                max_psi_allowed = MAX_PSI_HLX_ALLOWED;
            }
            else if (sType=='S')
            {
                min_phi_allowed = MIN_PHI_BETA_ALLOWED;
                max_phi_allowed = MAX_PHI_BETA_ALLOWED;
                min_psi_allowed = MIN_PSI_BETA_ALLOWED;
                max_psi_allowed = MAX_PSI_BETA_ALLOWED;
            }
            else
            {
                min_phi_allowed = MIN_PHI_LOOP_ALLOWED;
                max_phi_allowed = MAX_PHI_LOOP_ALLOWED;
                min_psi_allowed = MIN_PSI_LOOP_ALLOWED;
                max_psi_allowed = MAX_PSI_LOOP_ALLOWED;
            }

			//Working on Psi
			p2 = mStructure.getAtomCoordinate(cAA, " CA ");
			p3 = mStructure.getAtomCoordinate(cAA, " C  ");


			//get the angle gives the minimum distance
			torsionDelta = getMinAngle(mPoints, p3, p2, tPoints, 0, mPointsLastIndx, nRMSD);

			curTorsion = mStructure.getTorsion("psi", cAA);
			nTorsion = curTorsion - torsionDelta;		//new torsion after you rotate
			if (nTorsion<-180)
				nTorsion = 360+nTorsion;
			else if (nTorsion>180)
				nTorsion = 360-nTorsion;

			//cout<<"  curPsi= "<<curTorsion<<"  delAngle= "<<torsionDelta<<" nTorsion= "<<nTorsion<<" nRMSD= "<<nRMSD<<" cRMSD= "<<cRMSD<<endl;
			if (validTorsion(nTorsion, sType, "psi") &&
				(dontCheckColl || collisionFree(mStructure, avoidTrace, collStartAA, collEndAA, cAA, "psi", torsionDelta))){
				//rotate will be rotated in collisionFree
				mStructure.rotate(cAA, mStructure.numOfAA()-1, mStructure.getAtomIndx(cAA, " C  ")+1, p2, p3, torsionDelta);

				//udate moving points....because of rotation..we need to update these points everytime
				updatePoints(mPoints, p2,p3, 0, mPointsLastIndx, torsionDelta);

				cRMSD = nRMSD;
				//cout<<"roatatd with the first option."<<endl;

				rotated = true;
			}
			// you can use else section in case the angle is not valid then you rotate the maximum valid angle or a specific angle
			// for example if the new angle is out of valid range, suppose valid range is [-40, 30], then you can rotate so the
			// new torsion is -40 or 30. or you can rotate certain angle, for example 10 degrees or 20 degrees so you still in the range
			// if you don't like to do so...just comment else section
			else
			{
				if (relax){
					//cout<<"  NOT VALID PSI.."<<endl;
					//case #1 : you want to rotate so you are on the boundary of the valid range (you rotate the maximum you can)
					// first determine which boundary...which boundary is closer
					float delta1, delta2;

					delta1 = fabs(nTorsion - min_psi_allowed);
					if (delta1 > 180)
						delta1 = 360 - delta1;
					delta2 = fabs(nTorsion - max_psi_allowed);
					if (delta2 > 180)
						delta2 = 360 - delta2;

					//cout<<" delta1= "<<delta1<<" delta2= "<<delta2<<endl;

					if (delta1 <= delta2){

						if (dontCheckColl || collisionFree(mStructure, avoidTrace, collStartAA, collEndAA, cAA, "psi", curTorsion - min_psi_allowed)){

							//rotate to MIN boundary
							mStructure.rotate(cAA, mStructure.numOfAA()-1, mStructure.getAtomIndx(cAA, " C  ")+1, p2, p3, curTorsion - min_psi_allowed);

							//udate moving points....because of rotation..we need to update these points everytime
							updatePoints(mPoints, p2,p3, 0, mPointsLastIndx, curTorsion-min_psi_allowed);

							cRMSD = nRMSD;

							//mStructure.writePDB("portionAfterPsi.pdb",1,mStructure.numOfAA());
							//cout<<"min PSI applied "<<min_psi_allowed<<" cRMSD = "<<cRMSD<<endl;

							rotated = true;
						}
					}
					else
					{

                        if (dontCheckColl || collisionFree(mStructure, avoidTrace, collStartAA, collEndAA, cAA, "psi", curTorsion - max_psi_allowed)){
                            //rotate to Max boundary
                            mStructure.rotate(cAA, mStructure.numOfAA()-1, mStructure.getAtomIndx(cAA, " C  ")+1, p2, p3, curTorsion - max_psi_allowed);
                            //udate moving points....because of rotation..we need to update these points everytime
                            updatePoints(mPoints, p2,p3, 0, mPointsLastIndx, curTorsion - max_psi_allowed);

                            //get new RMSD
                            cRMSD = nRMSD;

                            //mStructure.writePDB("portionAfterPsi.pdb",1,mStructure.numOfAA());
                            //cout<<"max PSI applied "<<max_psi_allowed<<" cRMSD = "<<cRMSD<<endl;

                            rotated = true;
                        }
					}
				}
				//cout<<"  PSI now= "<<mStructure.getTorsion("psi", cAA)<<endl;
			}

			if (cRMSD < bestRMSD){
				bestRMSD = cRMSD;
				bestStructure = mStructure;
			}

			//check the RMSD
			if (cRMSD<=RMSDthr && rotated){
				itCntr = tIteration+1;
				break;
			}

            //check the SS type of the next AA...cAA+1
			sType = mStructure.AAs[cAA+1].SStype;
            if (sType=='H')
            {
                min_phi_allowed = MIN_PHI_HLX_ALLOWED;
                max_phi_allowed = MAX_PHI_HLX_ALLOWED;
                min_psi_allowed = MIN_PSI_HLX_ALLOWED;
                max_psi_allowed = MAX_PSI_HLX_ALLOWED;
            }
            else if (sType=='S')
            {
                min_phi_allowed = MIN_PHI_BETA_ALLOWED;
                max_phi_allowed = MAX_PHI_BETA_ALLOWED;
                min_psi_allowed = MIN_PSI_BETA_ALLOWED;
                max_psi_allowed = MAX_PSI_BETA_ALLOWED;
            }
            else
            {
                min_phi_allowed = MIN_PHI_LOOP_ALLOWED;
                max_phi_allowed = MAX_PHI_LOOP_ALLOWED;
                min_psi_allowed = MIN_PSI_LOOP_ALLOWED;
                max_psi_allowed = MAX_PSI_LOOP_ALLOWED;
            }
			//Working on phi .... phi of the next AA...
			p2 = mStructure.getAtomCoordinate(cAA+1, " N  ");
			p3 = mStructure.getAtomCoordinate(cAA+1, " CA ");

			//get the torsion angle gives the minimum distance
			torsionDelta = getMinAngle(mPoints, p3, p2, tPoints, 0, mPointsLastIndx, nRMSD);

			curTorsion = mStructure.getTorsion("phi", cAA+1);
			nTorsion = curTorsion - torsionDelta;				//new torsion after you rotate

			if (nTorsion<-180)
				nTorsion = 360+nTorsion;
			else if (nTorsion>180)
				nTorsion = 360-nTorsion;

			//cout<<"  curPhi= "<<curTorsion<<"  delAngle= "<<torsionDelta<<" nTorsion= "<<nTorsion<<" nRMSD= "<<nRMSD<<" cRMSD= "<<cRMSD<<endl;
			if (validTorsion(nTorsion, sType, "phi") &&
				(dontCheckColl || collisionFree(mStructure, avoidTrace, collStartAA, collEndAA, cAA+1, "phi", torsionDelta)))
			{

				//rotate
				mStructure.rotate(cAA+1, mStructure.numOfAA()-1, mStructure.getAtomIndx(cAA+1, " CA ")+1, p2, p3, torsionDelta);

				//udate moving points....because of rotation..we need to update these points everytime
				updatePoints(mPoints, p2,p3, 0, mPointsLastIndx, torsionDelta);

				cRMSD = nRMSD;

				//cout<<"rotated with the first option."<<endl;


				rotated = true;
			}
			// you can use else section in case the angle is not valid then you rotate the maximum valid angle or a specific angle
			// for example if the new angle is out of valid range, suppose valid range is [-40, 30], then you can rotate so the
			// new torsion is -40 or 30. or you can rotate certain angle, for example 10 degrees or 20 degrees so you still in the range
			// if you don't like to do so...just comment else section
			else
			{
				if (relax){
					//cout<<"  NOT VALID PHI.."<<endl;
					//case #1 : you want to rotate so you are on the boundary of the valid range (you rotate the maximum you can)
					// first determine which boundary...which boundary is closer
					float delta1, delta2;

					delta1 = fabs(nTorsion - min_phi_allowed);
					if (delta1 > 180)
						delta1 = 360 - delta1;
					delta2 = fabs(nTorsion - max_phi_allowed);
					if (delta2 > 180)
						delta2 = 360 - delta2;

					//cout<<" delta1= "<<delta1<<" delta2= "<<delta2<<endl;

					if (delta1 <= delta2){

						if (dontCheckColl || collisionFree(mStructure, avoidTrace, collStartAA, collEndAA, cAA+1, "phi", curTorsion - min_phi_allowed)){

							//rotate
							mStructure.rotate(cAA+1, mStructure.numOfAA()-1, mStructure.getAtomIndx(cAA+1, " CA ")+1, p2, p3, curTorsion - min_phi_allowed);

							//udate moving points....because of rotation..we need to update these points everytime
							updatePoints(mPoints, p2,p3, 0, mPointsLastIndx, curTorsion - min_phi_allowed);

							//get the current RMSD
							cRMSD = nRMSD;

							//mStructure.writePDB("portionAfterPhi.pdb",1,mStructure.numOfAA());
							//cout<<"min PHI applied "<<min_phi_allowed<<" cRMSD = "<<cRMSD<<endl;

							rotated = true;
						}
					}
					else
					{
                        if (dontCheckColl || collisionFree(mStructure, avoidTrace, collStartAA, collEndAA, cAA+1, "phi", curTorsion - max_phi_allowed)){
                            //rotate
                            mStructure.rotate(cAA+1, mStructure.numOfAA()-1, mStructure.getAtomIndx(cAA+1, " CA ")+1, p2, p3, curTorsion - max_phi_allowed);
                            //udate moving points....because of rotation..we need to update these points everytime
                            updatePoints(mPoints, p2,p3, 0, mPointsLastIndx, curTorsion - max_phi_allowed);

                            cRMSD = nRMSD;

                            //mStructure.writePDB("portionAfterPhi.pdb",1,mStructure.numOfAA());
                            //cout<<"max PHI applied "<<max_phi_allowed<<" cRMSD = "<<cRMSD<<endl;

                            rotated = true;
                        }
					}
				}

				//cout<<"  PHI now= "<<mStructure.getTorsion("phi", cAA+1)<<endl;
			}

			if (cRMSD < bestRMSD){
				bestRMSD = cRMSD;
				bestStructure = mStructure;
			}

			//check the RMSD
			if (cRMSD<=RMSDthr && rotated){
				itCntr = tIteration+1;
				break;
			}

			//cout<<"before phi: mPoints["<<pCntr<<"].x= "<<mPoints[pCntr].x<<" y= "<<mPoints[pCntr].y<<" z= "<<mPoints[pCntr].z<<endl;

			//getchar();
			cAA += direction;
		}

		itCntr++;
		//cout<<"  bestRMSD= "<<bestRMSD<<" nIterations= "<<itCntr<<endl;
		//getchar();
		//getchar();
	}

	mStructure = bestStructure;

	//cout<<" RMSD = "<<nRMSD<<endl;
	return bestRMSD;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
float randomFBCCD(Protein & portion1, Protein &portion2, int p2Indx, vector<char> aaTypes, vector<vector<Coordinate> > &AvoidTrace, float distTHR){

	float nRMSD = distTHR+1, minRMSD=9999.9;
	vector<Coordinate> tPnts, tPnts2,			//target points
						mPnts;			//movable points

	Protein bestStructure, originalPortion1 = portion1;
	vector<vector<Coordinate> > tmpAvoidTrace = AvoidTrace;
	AminoAcid tmpAA;
	int i;
	int nAAportion1 = portion1.numOfAA();
	int cntr = 0;
	int lastIndxOfPortion1 = nAAportion1+aaTypes.size();
	int nLoopAA = aaTypes.size();
	int MAXit = 10;						//max number of iterations

	//cout<<"portion1 size = "<<portion1.numOfAA()<<" loop= "<<nLoopAA<<" portion2= "<<portion2.numOfAA()<<" startIndx= "<<nAAportion1<<endl;

	//set first target points to be the first AA in the second portion
	tPnts.push_back (portion2.getAtomCoordinate(0, " N  "));
	tPnts.push_back (portion2.getAtomCoordinate(0, " CA "));
	tPnts.push_back (portion2.getAtomCoordinate(0, " C  "));


	//assumes that the length of the second portion is at least 2 AAs
	//set the new target points (the axis of second portion)
	tPnts2.push_back (triangleCenter(portion2.getAtomCoordinate(0, " N  "), portion2.getAtomCoordinate(0, " CA "), portion2.getAtomCoordinate(0, " C  ")));
	tPnts2.push_back(triangleCenter(portion2.getAtomCoordinate(portion2.numOfAA()-1, " N  "), portion2.getAtomCoordinate(portion2.numOfAA()-1, " CA "), portion2.getAtomCoordinate(portion2.numOfAA()-1, " C  ")));
	//a point in the middle
	short portion2middle = portion2.numOfAA()/2;
	if (portion2.numOfAA()>2){
		tPnts2.push_back (triangleCenter(portion2.getAtomCoordinate(portion2middle, " N  "),portion2.getAtomCoordinate(portion2middle, " CA "),portion2.getAtomCoordinate(portion2middle, " C  ")));
	}

	//for short loops...set the type to be loop
	if (aaTypes.size () < 6)
		for (i=0; i<aaTypes.size (); i++)
			aaTypes[i] = 'L';

    //add one aa at the beginning and another one at the end
    aaTypes.insert(aaTypes.begin(), portion1.AAs[nAAportion1-1].SStype);
    aaTypes.push_back(portion2.AAs[0].SStype);

    Protein loop;

	while (nRMSD>distTHR && cntr<MAXit){

		tmpAvoidTrace = AvoidTrace;

		portion1 = originalPortion1;

		//srand(nRMSD*2+5);

		int nAttempts=0;

		do{
			//cout<<"Trying to generate a random loop..."<<endl;
            //clear the loop if previously connected
            while (nAAportion1 != portion1.numOfAA())
                portion1.deleteAA(portion1.numOfAA()-1);

            portion1.AAsCollide.clear();
            portion1.atomsCollide.clear();

            loop = generateRandomStructure(aaTypes, (aaTypes.size() + cntr) * (nAttempts+1) );

            portion1.connect(loop, 'E');			//connect the loop with portion1

            //getchar();

        }while (portion1.doesCollide(0.5) && nAttempts++<20);

		//set the types
		//for (i=0;i<aaTypes.size(); i++){
			//portion1.AAs[i+nAAportion1].SStype = aaTypes[i];
		//}


		//get the last AA close to the first AA in portion2
		//set target and moving points
		mPnts.clear ();
		mPnts.push_back (portion1.getAtomCoordinate(lastIndxOfPortion1, " N  "));
		mPnts.push_back (portion1.getAtomCoordinate(lastIndxOfPortion1, " CA "));
		mPnts.push_back (portion1.getAtomCoordinate(lastIndxOfPortion1, " C  "));

		//portion1.writePDB("randomFBCCDbefore.pdb",1,portion1.numOfAA());

		//more cycles for short loops
		if (nLoopAA < 10)
			nRMSD = FBCCDnPoints(portion1,mPnts,tPnts,tmpAvoidTrace, MAX(nAAportion1-4, 1), lastIndxOfPortion1, 0, lastIndxOfPortion1, distTHR, 500, 1, true, true, 0.1, false);
		else
			nRMSD = FBCCDnPoints(portion1,mPnts,tPnts,tmpAvoidTrace, MAX(nAAportion1-2,1 ), lastIndxOfPortion1, 0, lastIndxOfPortion1, distTHR, 300, 1, true, true, 0.1, false);

		//cout<<"first RMSD = "<<nRMSD<<endl;
		if (nRMSD < 0.1){

			//portion1.AAsCollide.clear();
            //portion1.atomsCollide.clear();
			//check for collision
			//if (!portion1.doesCollide(0.5)){
				//if the aa RMSD less than 0.1 then simple append and return
				portion1.append(portion2, 1, portion2.numOfAA()-1);
				//cout<<"  Done..."<<endl;

				//set the trace for the new structure built
				/*
				vector<Coordinate> tmpTrace;
				for (i=nAAportion1-1; i<lastIndxOfPortion1+1; i++)
					tmpTrace.push_back (portion1.getAtomCoordinate(i, " CA "));
				AvoidTrace.push_back (tmpTrace);
				*/

				return nRMSD;
			//}
		}

		portion1.connect(portion2, 'E');		//connect the second portion

		//portion1.writePDB("randomFBCCDconnected.pdb",1, portion1.numOfAA());


		//set the movable points
		mPnts.clear  ();
		mPnts.push_back (triangleCenter(portion1.getAtomCoordinate(lastIndxOfPortion1, " N  "), portion1.getAtomCoordinate(lastIndxOfPortion1, " CA "), portion1.getAtomCoordinate(lastIndxOfPortion1, " C  ")));
		mPnts.push_back (triangleCenter(portion1.getAtomCoordinate(portion1.numOfAA()-1, " N  "), portion1.getAtomCoordinate(portion1.numOfAA()-1, " CA "), portion1.getAtomCoordinate(portion1.numOfAA()-1, " C  ")));

		if (portion2.numOfAA()>2){
			mPnts.push_back (triangleCenter(portion1.getAtomCoordinate(lastIndxOfPortion1 + portion2middle," N  "),portion1.getAtomCoordinate(lastIndxOfPortion1 + portion2middle," CA "),portion1.getAtomCoordinate(lastIndxOfPortion1 + portion2middle," C  ")));
		}


		//erase the helix from the traces to avoid
		tmpAvoidTrace.erase (tmpAvoidTrace.begin ()+p2Indx);

		if (nLoopAA < 10)
			nRMSD = FBCCDnPoints(portion1,mPnts,tPnts2,tmpAvoidTrace, lastIndxOfPortion1, MAX(nAAportion1-4,1), 0, MIN(lastIndxOfPortion1+4,portion1.numOfAA()-1), distTHR, 500, -1, true, true, 0.3, false);
		else
			nRMSD = FBCCDnPoints(portion1,mPnts,tPnts2,tmpAvoidTrace, lastIndxOfPortion1, MAX(nAAportion1-2,1), 0, MIN(lastIndxOfPortion1+4,portion1.numOfAA()-1), distTHR, 300, -1, true, true, 0.3, false);

		//portion1.writePDB("randomFBCCDafter.pdb",1,portion1.numOfAA());

        if (nRMSD<0.5){
            //cout<<"second RMSD = "<<nRMSD<<" .... will return"<<endl;
            return nRMSD;
        }

		if (nRMSD<minRMSD){
			//portion1.AAsCollide.clear();
            //portion1.atomsCollide.clear();
			//if (!portion1.doesCollide(0.5)){
				minRMSD = nRMSD;
				bestStructure = portion1;
			//}
		}

		cntr++;
		//cout<<portion1.getTorsion("phi", startIndx+1);

		//cout<<"iteration... best RMSD so far "<<minRMSD<<" current RMSD = "<<nRMSD<<endl<<"Press any key.."<<endl;getchar();
		//cout<<"nRMSD now is : "<<nRMSD<<endl;getchar();
	}

	//bestStructure.connect(portion2, 'E');
	portion1 = bestStructure;

	//set the trace for the new structure built
	/*
	vector<Coordinate> tmpTrace;
	for (i=nAAportion1-1; i<lastIndxOfPortion1+1; i++)
		tmpTrace.push_back (portion1.getAtomCoordinate(i, " CA "));
	AvoidTrace.push_back (tmpTrace);
	*/

	//cout<<"Done..."<<endl;
	return minRMSD;

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
float randomFBCCD1Point(Protein & portion1, Coordinate trgtPnt, short indxMov, vector<char> aaTypes, vector<vector<Coordinate> > &AvoidTrace, float distTHR){

	float nRMSD = distTHR+1, minRMSD=9999.9;
	vector<Coordinate> tPnts(1),		//target points
						mPnts(1);			//movable points

	Protein bestStructure,
			originalPortion1 = portion1;
	int MAXit = 10;			//max number of iterations

	AminoAcid tmpAA;
	int i;
	int nAAportion1 = portion1.numOfAA();
	int cntr = 0;
	int lastIndxOfPortion1 = nAAportion1+aaTypes.size()-1;
	int nLoopAA = aaTypes.size();

	//cout<<"portion1 size = "<<portion1.numOfAA()<<" indxMov= "<<indxMov<<" loop= "<<nLoopAA<<" startIndx= "<<nAAportion1<<endl;

	//set first target points to be the first AA in the second portion
	tPnts[0] = trgtPnt;

    //add one aa at the beginning
    aaTypes.insert(aaTypes.begin(), portion1.AAs[nAAportion1-1].SStype);

    Protein loop;

	while (nRMSD>distTHR && cntr<MAXit){

		portion1 = originalPortion1;

		srand(nRMSD*2+5);

		short nAttempts=0;
		do{
			//cout<<"Trying to generate a random loop...for 1 point function "<<aaTypes.size()<<" nAttempts= "<<nAttempts<<endl;
            //clear the loop if previously connected
            while (nAAportion1 != portion1.numOfAA())
                portion1.deleteAA(portion1.numOfAA()-1);

            portion1.AAsCollide.clear();
            portion1.atomsCollide.clear();

            loop = generateRandomStructure(aaTypes, indxMov*(nAttempts+1));

            portion1.connect(loop, 'E');			//connect the loop with portion1

            //getchar();
            //portion1.writePDB("randomFBCCDafter.pdb",1,portion1.numOfAA());

        }while (portion1.doesCollide(0.5) && nAttempts++ < 20);

		//get the last AA close to the first AA in portion2
		//set target and moving points
        mPnts[0].x = (portion1.getAtomCoordinate(indxMov, " N  ").x + portion1.getAtomCoordinate(indxMov, " CA ").x)/2;
        mPnts[0].y = (portion1.getAtomCoordinate(indxMov, " N  ").y + portion1.getAtomCoordinate(indxMov, " CA ").y)/2;
        mPnts[0].z = (portion1.getAtomCoordinate(indxMov, " N  ").z + portion1.getAtomCoordinate(indxMov, " CA ").z)/2;

		//portion1.writePDB("randomFBCCDbefore.pdb",1,portion1.numOfAA());

		nRMSD = FBCCDnPoints(portion1,mPnts,tPnts,AvoidTrace, MAX(nAAportion1-2,1 ), lastIndxOfPortion1, 0, lastIndxOfPortion1, distTHR, 300, 1, true, true, 0.1, false);


		//portion1.writePDB("randomFBCCDafter.pdb",1,portion1.numOfAA());

		if (nRMSD<minRMSD){
			minRMSD = nRMSD;
			bestStructure = portion1;
		}

		cntr++;

		//cout<<"iteration... best RMSD so far "<<minRMSD<<" current RMSD = "<<nRMSD<<endl<<"Press any key.."<<endl;
		//getchar();
	}

	portion1 = bestStructure;

	//cout<<"Done..."<<endl;getchar();

	return minRMSD;

}
#endif
