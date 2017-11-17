#ifndef PROTEIN_H
#define PROTEIN_H

#include<iostream>
#include<iomanip>
#include<string>
#include<fstream>
#include<vector>
#include <stdio.h>
#include <sstream>
#include <windows.h>



#include "constants.h"
#include "geometry.h"
#include "utilityfunctions.h"


using namespace std;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////// Some Constants and Data Structures ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// DATA STRUCTURE TO STORE PHI AND PSI ANGLES
//struct Torsion
//{
//	double phi;
//	double psi;
//	//Initialize torsion angles to 999.00
//	Torsion() : phi(999.00), psi(999.00) {}
//};

// list of missing AAs if any.
struct missingAA
{
	string chr3;				//AA in 3 letters representation
	char chr1;					//AA in one letter representation
	int num;					//sequence number
	char resInsertion;			//Code for insertions
	char SStype;				//the type of the AA...H helix, S strand, or L loop

	missingAA(): chr3(""), chr1(' '), num(0), resInsertion(' '), SStype('L') {}

};

								//////////////////////////////////////////////////////////////////////////
								//////////////////// Helices Secondary Structure Struct //////////////////
								//////////////////////////////////////////////////////////////////////////
struct HelicesSecondaryStructure
{
	int startIndx;		//The index (AA vector indx) of the first AA in the helix
	int endIndx;		//The index (AA vector indx) of the last AA in the helix
	int serialNum;		//Helix serial Number
	int hlxID;			//Helix identifier
	int nAAcur;			//The length (the missing AA's are not counted) (# of AA)  of the helix
	int nAA;			//The Actual length of the helix...missing AA's are counted
	string comment;		//any comment

	int type;			//Contains the type of the helix

						//The type as follows:
						// 1  Right-handed alpha (default)		6  Left-handed alpha
						// 2  Right-handed omega				7  Left-handed omega
						// 3  Right-handed pi					8  Left-handed gamma
						// 4  Right-handed gamma				9  2/7 ribbon/helix
						// 5  Right-handed 3/10					10  Polyproline

	//Initializer
	HelicesSecondaryStructure() : startIndx(-1), endIndx(-1), serialNum(0), hlxID(0), nAAcur(0), nAA(0), comment(""), type(1) {}
};
								//////////////// END OF Helices SECONDARY STRUCTURES /////////////////////

								//////////////////////////////////////////////////////////////////////////
								//////////////////// Sheets Secondary Structure Struct ///////////////////
								//////////////////////////////////////////////////////////////////////////
struct SheetsSecondaryStructure
{
	//store the information of the strand
	int startIndx;		//The index (AA vector indx) of the first AA in the strand
	int endIndx;		//The index (AA vector indx) of the last AA in the strand
	int strandNum;		//strand number
	string sheetID;		//sheet identifier
	int nStrand;		//number of strands in current sheet

	int nAAcur;			//The length (the missing AA's are not counted) (# of AA)  of the beta strand
	int nAA;			//The Actual length of the strand...missing AA's are counted
	int sense;			//Contains strand sense with respect to previous
						//0 first strand			1 parrallel		-1 antiparallel
	//The following fields identify two atoms involved in a hydrogen bond, the first in the current strand and the second in the previous strand.  \
    //These fields should be blank for strand 1 (the first strand in a sheet).
	string curAtomName;		//the name of the atom in the current strand involved in a hydrogen bond
	string prevAtomName;	//the name of the atom in the previous strand involved in a hydrogen bond
	int AACurIndx;			//the index of the AA in the current strand involved in a hydrogen bond
	int	AAPrevIndx;			//the index of the AA in the previous strand involved in the hydrogen bond


	//Initializer
	SheetsSecondaryStructure() : startIndx(-1), endIndx(-1), strandNum(0), sheetID(""), nStrand(0), nAAcur(0), nAA(0), sense(1), AACurIndx(-1), AAPrevIndx(-1) {}
};
								//////////////// END OF Sheets SECONDARY STRUCTURES /////////////////////

								/////////////////////////////////////////////////////////////////////////
								////////////////////// ATOM STRUCTURE ///////////////////////////////////
								/////////////////////////////////////////////////////////////////////////
struct Atom
{
	string name;								//Atom name ...considered as 4 characters
	char locIndicator;							//Alternate_location_indicator			position (16,1)
	string occupancy;							//Occupancy							//	position (54,6)		right alignment
	string tempFactor;							//Temperature factor				//	position (60,6)		right alignment
	string charge;								//Charge							//	position (78,2)

	Coordinate coord;							// X,Y,Z position of the atom

	char type;									//The type of the atom..(N,C,O....)
	bool isSideChain;							//side chain flag....true if it is a side chain atom...false otherwise

	//Initializer
	Atom() : name("NONE"), locIndicator(' '), occupancy(""), tempFactor(""), charge(""), type(' '), isSideChain(false) {}
};
								////////////////// END OF ATOM STRUCTURE ////////////////////////////////


								/////////////////////////////////////////////////////////////////////////
								/////////////////// AMINO ACID STRUCTURE ////////////////////////////////
								/////////////////////////////////////////////////////////////////////////
struct AminoAcid
{
	vector<Atom> atoms;			//list of atoms from atom struct
	/*
					atoms structure will look like
					N
					N-H atoms if any
					Ca
					Ca-H atoms if any
					CB if any
					CB-H atoms if any
					SideChain atoms followed by their H atoms
					...
					...
					C
					O
					OXT
	*/
	char chr1;					//AA name (1-character)
	string chr3;				//AA name (3-character)
	string chain;				//the current chain
	int num;					//AA serial number...same as the number in PDB file
	char resInsertion;			//Code for insertions of residues	//	position (26,1)
	Coordinate coord;			//could be the center of side-chain or the all-atoms center
	char SStype;				//stores the type of secondary structure this AA represent
								//  H : Helix
								//	S : Sheet
								//  L : Loop or other;
	Coordinate ScEndPoint;		//The end point that represents the end of the side chain.... a vector from CA to this point represent the size and
								// direction of the side chain
	char whtInCoord;			//what is stored in coord variable...could be
								//	S : side-chain center
								//	A : All-atoms center
								//  B : backbone-atoms center
								//	N : Nothing Interesting (initial)
	char atomsIncluded;			//H heavy atoms, A all, N nothing (initial)..included in calculation stored in coord
	double gyration;			//stores the gyration radius of amino acid (for side chain)...if it is equal to 0 then gyCoord is invalid
	Coordinate gyCoord;		    // cotains the Coordinate of the radius of gyration
	Torsion angles;				//data structure to store phi psi angles

	//Initializer
	AminoAcid() : chr1(' '), chr3(""), chain(""), num(0), resInsertion(' '), whtInCoord('N'), atomsIncluded('N'), gyration(0), SStype('L') {}
};
								////////////////// END OF AMINO ACID STRUCTURE /////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// End Of Constants and Data Structures /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// List of Functions ////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Protein copyBBone(Protein originalPortion, int startIndx, int endIndx);		//copy the backbone from a particular portion of protein

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// ////////// End Of The List //////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////// PROTEIN CLASS PROTOTYPE ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Protein
{
	friend void copyBBone(Protein, int, int, Protein &);	//copy backbone from one protein to another for specified length
public:
	Protein();									//The default constructor
	Protein(string, string = "");				//Another Constructor..given file name (path) and chain specifier
	void initialize();							//initialize variables


	string ID;									//The ID of the protein that the portion came from
	Coordinate centOfCharge;					//center of charge
	vector<AminoAcid> AAs;						//Contains the list of AA's that form the Portion
	vector<missingAA> missingAAs;				// a list of amino acids are missing from 3D structure
	vector<string> header;						//The header information from the PDB file
	vector<HelicesSecondaryStructure> hlces;			//The list of hlces in the portion
	vector<SheetsSecondaryStructure> sheets;			//The list of sheets in the portion
	vector<int> AAsCollide;						//stores the Indeces of first 2 AA's that collide with each other ...
	vector<int> atomsCollide;					//stores the Indeces of atoms that collide with each other within the 2 AA's collide from above variable


	void read(string, bool &moreChains, string &nextChain, string = "");					//given a pdb file path and a target chain....reads the PDB file
	void writePDB(string, int, int, bool = false);	//write a specified range of the portion to a PDB file
	void writeSCModel(string, int, int);			//write a specified range of AA's using the simplified model of Side Chain
	void writeAAInfo(string);						//write information of AA's in the portion to a file


	inline int numOfAA();						//returns the number of AA in the Protein
	inline int numOfAtoms();					//returns the number of atoms in the whole portion
	inline int numOfAtoms(int);					//given the indx of AA....returns the number of atoms in that AA
	int numOfSCAtoms(int, char = 'H');			//get the number of side chain atoms, H: heavy or A: all
	int getAAIndx(int);							//given AA num (as in PDB) returns the indx of this AA in the AA vector .. .-1 if not found
	int getAtomIndx(int, string);				//given AA indx and the name of atom (or substring)...returns the indx of this atom or -1 if not exist
	int getAtomIndx2(int, string);
	string getSequence(int, int);				//return the corresponding sequence with a specified range
	string getSequenceWmissing(int, int);		//return the corresponding sequence with missing AAs
	Coordinate getAtomCoordinate(int, string);	//return the coordinate of a given atom (atom name or a substring of it) in a particular AA
	AminoAcid operator() (int);					//overload of () operator...reprsents the rank...for example Protein(1) = Protein.AAs[0]...but it is slow
	double getTorsion(string, int);				//return the torsion angle phi or psi...given AAIndx
	double getChi(int, int);					//return the chi angle for a particular AA
	int doesCollide(double = 0.0);				//check if there is a collision within AA's of the portion....stores those AA's and atoms in AAsCollide and atomsCollide
	double VDW(double, double , double);		//find VDW energy between two atoms
	bool isSSAA(int);							//return the kind of secondary structure the AA lies in. H helix, S strand, N not a SS Amino acid
	int inHlces(int);							//return the index of helix the AA is located... -1 if it is not a helix AA


	void fillMissingAAs();						//fill missing AAs if any in missingAAs data structure
	void removeHAtoms();						//remove H atoms from the portion
	void renameAA(int, string);					//rename the AA (AAindx given).... nothing will happen to coordinates or atoms
	void deleteAA(int);							//delete a specific AA from the list....just drop it off
	void deleteAtom(int, string);				//delete a particualr atom froma specific AA
	void append(Protein, int, int, int = 0);	//Append a number of AA to the end of the portion...given a range and the gap u wanna leave b/w current portion and the appended one
	void setSCCenter(int, char = 'H');			//given an index of an AA....this stores the center of sidechain in AAs.coord and set whtInCoord to S
	void setAACenter(int, char = 'H');			//given AA index...compute the center (mass center of all atoms) of this AA.and set whtInCoord to A
	void setBBCenter(int, char = 'H');			//given AA indx ... copmute the center(mass center) of all backbone atoms
	void setRgyration(int, char ='H');			//set the radius of gyration and the gyration coordinate...H heavy atoms included or A all atoms
	void setCentOfCharge();						//set the center of charge for all atoms in the portion
	void setAxisSegments(int, int, float, vector<Coordinate> &, vector<bool> &shortHelix, int);		//stor the axis of a segment every "given distance" as a set of points
	void setAxisSegments2 (int startIndx, int endIndx, float dist, vector<Coordinate> &axis, Coordinate last, vector<bool> &shortHelix, int);

	void buildSS();								//build secondary structure vector.....
	void translateBy(Coordinate);				//given a translation value stored in a point....do transformation
	void translateBy(int);						//translate (move) the portion by number of Ca atoms
	void overlapLine(Coordinate, Coordinate, Coordinate, Coordinate);	//overlap the portion with a line....Do translation and rotation
	void connect(Protein, char = 'E');			//Connect a list of AAs with the portion from start (S) or end (E)
	void concat(Protein, char = 'E', double = -60, double = 60);	//concat a list of AAs with the portion from start or end according to the phi and psi given
	void plugSideChain(int, Protein);			//plug a side chain for a specific AA
	void setScEndPoint(int, char = 'H');//set the side chain end point for a given AA, H: heavy atoms to be included A: all atoms
	void torsion2coord();						//generate xyz coordinate system according to torsion (phi, psi) angles....will delete the previous coordinate system
												//BE CAREFULL...you will loose side chains if there is any
	void rotateSideChain(int,int,Coordinate, Coordinate, double);		//given the indx of AA and the indx of atom and the bond  coordinates...rotate the side change by an angle
	void rotate(int,int, int, Coordinate, Coordinate, double);			//rotate a specific number of AAs around line

private:
	string path;								//the path from where the protein came
	char chr3ToChr1(string);					//converts from 3-letters format to 1-letter format
	char getAtomType(string);					//returns the type of given Atom...N, C, O, or S...
	bool isSideChainAtom(string);				//true if the atom is a sidechain atom
	void getCPosition(Protein &, char = 'E');	//determine the position of C atom when connecting or concatenating to portions (it moves the connected portion)
	AminoAcid reOrderAtoms(int);
	int getNumOfChi(char);						//returns the number of chi angles for a particular AA
	bool isHeavyAtom(char);						//check if the atom is heavy or not
	void sortHlces(vector<HelicesSecondaryStructure> &, const unsigned int, unsigned int);	//sort hlces according to the indeces of the first AA
	void sortStrands(vector<SheetsSecondaryStructure> &, const unsigned int, unsigned int);	//sort sheets according to the indeces of the first AA
	void setShortAxis(int, int, float, vector<Coordinate> &);

};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////// END OF PROTEIN CLASS PROTOTYPE /////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// PROTEIN Class Implementation /////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Protein::Protein()		//constructor
{
	//initialize variables
	initialize();
}
////////////////////////////////////////////////////////////////////////////////////////
Protein::Protein(string filePath, string chain)		//constructor
{
    bool temp = false;
	read(filePath, temp, chain);
}
////////////////////////////////////////////////////////////////////////////////////////////
/*
Protein::~Protein ()
{
	AAs.~vector();
	header.~vector();
	hlces.~vector();
	sheets.~vector ();
	AAsCollide.~vector ();
	atomsCollide.~vector();

}
*/
////////////////////////////////////////////////////////////////////////////////////////////
void Protein::initialize()
{
	path = "Unknown";
	ID = "Unknown";
	AAs.clear();
	missingAAs.clear ();
	header.clear();
	hlces.clear();
	sheets.clear();
	AAsCollide.clear ();
	atomsCollide.clear();

}
/////////////////////////////////////////////////////////////////////////////////////////////

// read PDB file
///modified to return number of chains
void Protein::read(string fileName, bool &moreChains, string &nextChain, string chain)
{
	string line,tmpChain, tempLine;

	ifstream infile;
	ifstream tempInfile;

	bool exitReading = false;
	bool reachedChain = false;

	//initialize variables .... every time u wanna read a protein
	initialize();

	path = fileName;		//the path is where the portion came from

	//get portion ID
	if (fileName.length() > 8)
		if ((atoi(fileName.substr (fileName.length () - 9,1).c_str ()) != 0) &&		// "\" symbol in windows
			(fileName.substr (fileName.length () - 9,1) != "\\") &&					// "/" symbol
			(atoi(fileName.substr (fileName.length ()- 9,1).c_str()) != 92 ))		// "\" symbol in linux
			ID = fileName.substr(fileName.length()-9,4);
		else
			//the path does not contain "\" or "/" but it may contain the chain ID
			ID = fileName.substr (fileName.length () - 8,4);
	else
		if (fileName.length() == 8)
			ID = fileName.substr (fileName.length ()-8,4);
		else
			//the ID in this case could be wrong since the name is irrelative to the ID
			ID = fileName.substr (0, fileName.length() - 4);
	infile.open(fileName.c_str());
	tempInfile.open(fileName.c_str());
	if (!infile)
	{
		errMsg("Protein", "read", "Unable to open "+ fileName);

		exit(1);
	}
	else
	{
		AminoAcid aa;
		Atom tmpAtom;

		while((!infile.eof()) && (!exitReading))
		{

			getline(infile, line);
			getline(tempInfile, tempLine);                                                //get a line from pdb file

			if(line.size() > 23 && (line.substr(21,1).c_str() == chain || chain == "") && (line.substr(0,4) == "ATOM"))
                reachedChain = true;

			if(line.substr (0,3) == "TER" && reachedChain == true)
			{
			    getline(tempInfile, tempLine);
			    if(tempLine.substr(0,4) == "ATOM")
                {
                    nextChain = (tempLine.substr(21,1).c_str());
                    cout << nextChain << endl;
                    moreChains = true;
                }

			}

			if ((line.substr(0,4) == "ATOM") &&										//ATOM word token
				(line[17] != ' ') &&												//not a DNA
				((line.substr(21,1).c_str() == chain) || (chain == "")) &&			//same target chain
				(!exitReading))														//remain atoms to be read
			{
				/*****
						Samples:
						ATOM      1  N   VAL A   3      18.481  17.489  41.966  1.00 33.08           N
						ATOM      6  CG1 VAL A   3      14.652  17.181  41.934  1.00 33.32           C
				*****/

				/*****
						 begins with first chain if no target chain was sent
				*****/
				if ((numOfAA() == 0))
				{
					tmpChain = line.substr(21,1).c_str();
					aa.num = -1000; //atoi(line.substr(22,4).c_str());

				}
				/*****
						the current version of the system does not work with portions have code of residue inserting....
						the portion will be all the part before the first AA has a code of inserting
				*****/
				if (line[26] != ' ')
				{
					cout<<ID<<" has code of residue inserting"<<endl;
					exitReading = true;
				}

				/*****
						saves the first chain, whatever this chain is A,B,or C...etc		terminate if u reach the conditions below
				*****/
				if (tmpChain!=line.substr(21,1).c_str() ||								//another chain or
					(atoi(line.substr(22,4).c_str()) < aa.num) ||						//AA with less seq num
					(line[26] != ' ') ||												//Code for inserting residue
					(line.substr (0,3) == "TER"))           //TER reserved word
					{
					    exitReading = true;
					}
				else
				{
					/*****
							check if we reach a new AA and it is not a dublicate AA (by locIndicator)
					*****/
					if (atoi(line.substr(22,4).c_str()) > aa.num )
					{
						aa.chr3			= line.substr(17,3);
						aa.num			= atoi(line.substr(22,4).c_str());
						aa.chain		= line.substr(21,1).c_str();
						aa.chr1			= chr3ToChr1(line.substr(17,3).c_str ());
						aa.resInsertion = line[26];

						/*****
								if it is the only alternative take it...or it has more than one alternative...take alternative A
						*****/
						if ((line[16] == ' ') || (line[16] == 'A'))
							AAs.push_back(aa);							//Push the information of the AA
					}
					/*****
							if the atom is not for a duplicated AA then take it or it is alternative A
					*****/
					if ((line[16] == ' ') || (line[16] == 'A'))
					{
						tmpAtom.name         = line.substr(12,4);
						tmpAtom.locIndicator = line[16];
						tmpAtom.coord.x		 = atof(line.substr(30,8).c_str());
						tmpAtom.coord.y		 = atof(line.substr(38,8).c_str());
						tmpAtom.coord.z		 = atof(line.substr(46,8).c_str());
					//	tmpAtom.type		 = line[77];
						/*****
								if one of atoms has no type then get it
						*****/
					//	if (int (tmpAtom.type) == 0)
							tmpAtom.type = getAtomType(tmpAtom.name);
						/*****
								Extra Information for atom
						*****/
						if (line.length()>54)
							tmpAtom.occupancy    = line.substr(54,6).c_str();
						if (line.length()>60)
							tmpAtom.tempFactor	 = line.substr(60,6).c_str();
						if (line.length()>78)
							tmpAtom.charge       = line.substr(78,2).c_str();
						/*****
								check if it is a side chain atom..set the flag
						*****/
						tmpAtom.isSideChain  = isSideChainAtom(tmpAtom.name);
						/*****
								push the atom to the AAs DataStructure
						*****/
						AAs[AAs.size()-1].atoms.push_back(tmpAtom);

					}
				}
			}
			else
			{
				/*****
						save header information from the beginning of the Pdb file..exclude any chain before the specified one (if any has been specified)
				*****/
				if ((!exitReading) && (numOfAA()==0) && (line.substr(0,4) != "ATOM"))
				{
					//cout<<line<<endl;
					header.push_back(line);
				}
			}
		}

		/*****
				close the file
		*****/
		infile.close();

		/*****
				reOrder Atoms...so every atom is followed by its H atoms in the AAs vector..last three atoms are C, O, and OXT
				and set the Side chain end point (where the line from Ca to this end point represents the direction and the length of the side chin
		*****/
		int i;
		for (i=0;i<numOfAA();i++)
		{
			AAs[i] = reOrderAtoms(i);
			setScEndPoint(i);
		}

		/*****
				build the secondary structure information
		*****/
		buildSS();

		/*****
				set the information of secondary structure for each AA
		*****/
		for (i=0;i<numOfAA(); i++)
		{
			for (int m=0; m<sheets.size (); m++)
				if ((i>= sheets[m].startIndx) && (i <= sheets[m].endIndx ))
					AAs[i].SStype = 'S';
			for (int k=0; k<hlces.size (); k++)
				if ((i>=hlces[k].startIndx ) && (i<= hlces[k].endIndx  ))
					AAs[i].SStype = 'H';
		}

		/*****
					fill missingAAs data structure for any missing residue from the 3D structure
		*****/
		fillMissingAAs();
	}
}
////////////////////////////////////////////////////////////////////////////////
void Protein::writePDB(string outFile,int startAARank, int endAARank, bool wHeader)
// The whole Protein could be written to a file or a portion of it
//startAARank is the rank of the first AA to be written to outfile...for example 4 is the 4th AA
//endAARank is the rank of the last AA would be written to the outfile PDB file
//if you want to print with the header info in the original pdb file..then send wHeader = true;
{
	int tmpNumOfAA = numOfAA();

	//assure from the ranks given
	if ((startAARank>=1) && (startAARank <= tmpNumOfAA) && (endAARank <= tmpNumOfAA))
	{
		int i,
			j;
		ofstream out;
		out.open(outFile.c_str());
		if (!out) {
			errMsg("Protein", "writePDBFile", "Unable to open " + outFile, true);
		}

		/*****
				write header information...if the header information is chosen to be written then no need to
				do the next operation which writes the secondary structures in the specified range if the SS
				has been built using buildSS() before;
		*****/
		if (wHeader)
			for (i=0;i<header.size();i++)
				out<<header[i]<<endl;
		else
		{
			/*****
					In the case of buildSS() has been called, SS existance, and wHeader flag = false....write the list of SS first
			*****/
			int hCounter = 1;	//counter for helices
			int sCounter = 1;	//counter for sheets

			/*****
					write hlces if any
			*****/
			for (i=0;i<hlces.size();i++)
			{
				//The whole hlx is within the range
				if ((startAARank - 1 <= hlces[i].startIndx) && (endAARank - 1 >= hlces[i].endIndx))
				{
					out<<"HELIX  "<<setw(3)<<right<<hlces[i].serialNum<<" "<<setw(3)<<right<<hlces[i].hlxID<<" "<<AAs[hlces[i].startIndx].chr3<<" "<<AAs[hlces[i].startIndx].chain<<" "
					   <<setw(4)<<right<<AAs[hlces[i].startIndx].num<<AAs[hlces[i].startIndx].resInsertion<<" "<<AAs[hlces[i].endIndx].chr3<<" "<<AAs[hlces[i].endIndx].chain<<" "<<setw(4)<<right<<AAs[hlces[i].endIndx].num
					  <<AAs[hlces[i].endIndx].resInsertion<<setw(2)<<right<<hlces[i].type<<setw(30)<<left<<hlces[i].comment<<" "<<setw(5)<<right<<hlces[i].nAAcur<<endl;
						++hCounter;
						continue;
				}

				//both ends within the range
				if ((startAARank > hlces [i].startIndx ) && (endAARank < hlces [i].endIndx))
				{
					out<<"HELIX  "<<setw(3)<<right<<hCounter<<" "<<setw(3)<<right<<hCounter<<" "<<
					//new AA1
					AAs[startAARank-1].chr3<<" "<<AAs[hlces[i].startIndx].chain<<" "<<setw(4)<<right<<AAs[startAARank-1].num<<AAs[startAARank-1].resInsertion<<" "<<
					//new AA2
					AAs[endAARank-1].chr3<<" "<<AAs[hlces[i].endIndx].chain<<" "<<setw(4)<<right<<AAs[endAARank-1].num<<AAs[endAARank-1].resInsertion<<setw(2)<<right<<
					hlces[i].type <<setw(30)<<left<<hlces[i].comment<<" "<<setw(5)<<right<<endAARank - startAARank + 1<<endl;
					++hCounter;
					continue;
				}
				//The lower end is within the range
				if ((startAARank - 1<= hlces[i].startIndx) && (endAARank > hlces[i].startIndx) && (endAARank - 1 < hlces[i].endIndx))
				{
					out<<"HELIX  "<<setw(3)<<right<<hCounter<<" "<<setw(3)<<right<<hCounter<<" "<<AAs[hlces[i].startIndx].chr3<<" "<<
						AAs[hlces[i].startIndx ].chain<<" "<<setw(4)<<right<<AAs[hlces[i].startIndx].num<<AAs[hlces[i].startIndx].resInsertion<<" "
					   //new AA
					   <<AAs[endAARank - 1].chr3<<" "<<AAs[hlces[i].endIndx].chain<<" "<<setw(4)<<right
		 			   //endAARank - 1 is the new end right now
					   <<AAs[endAARank - 1].num<<AAs[endAARank-1].resInsertion<<setw(2)<<right<<hlces[i].type<<setw(30)<<left<<hlces[i].comment<<" "<<setw(5)<<right
					   //The new length
					   <<endAARank - hlces[i].startIndx<<endl;
					   ++hCounter;
					   continue;
				}

				//The upper end is within the range
				if ((startAARank -1 > hlces[i].startIndx) && (startAARank -1 <= hlces[i].endIndx) && (endAARank - 1 >= hlces[i].endIndx))
				{
					out<<"HELIX  "<<setw(3)<<right<<hCounter<<" "<<setw(3)<<right<<hCounter<<" "
						//New AA
						<<AAs[startAARank - 1].chr3<<" "<<AAs[hlces[i].startIndx].chain<<" "<<setw(4)<<right
						//startAARank - 1 is the new start right now
						<<AAs[startAARank -1 ].num<<AAs[startAARank-1].resInsertion<<" "<<AAs[hlces[i].endIndx].chr3<<" "<<AAs[hlces[i].endIndx].chain
						<<" "<<setw(4)<<right<<AAs[hlces[i].endIndx].num<<AAs[hlces[i].endIndx ].resInsertion<<setw(2)<<right<<hlces[i].type
						<<setw(30)<<left<<hlces[i].comment<<" "<<setw(5)<<right
					  //new length
					  <<hlces[i].endIndx - startAARank + 2<<endl;
						++hCounter;
						continue;
				}
			}

			/*****
					write sheets if any...the operation of hlces above should be done for sheets as well
			*****/
			for (i=0;i<sheets.size();i++)
			{
				//if ((startAARank - 1 <= sheets[i].startIndx) && (endAARank - 1 >= sheets[i].endIndx))
				//{
					out<<"SHEET  "<<setw(3)<<right<<sheets[i].strandNum<<" "<<sheets[i].sheetID<<setw(2)<<sheets[i].nStrand<<" "
						<<AAs[sheets[i].startIndx ].chr3<<" "<<AAs[sheets[i].startIndx].chain<<setw(4)<<AAs[sheets[i].startIndx].num
						<<AAs[sheets[i].startIndx ].resInsertion<<" "<<AAs[sheets[i].endIndx].chr3<<" "<<AAs[sheets[i].endIndx].chain
						<<setw(4)<<AAs[sheets[i].endIndx].num<<AAs[sheets[i].endIndx].resInsertion<<setw(2)<<sheets[i].sense<<" ";
					//write information for AA's involved in Hydrogen Bond starting from second strand
					if (sheets[i].strandNum != 1)
					{
						out<<sheets[i].curAtomName<<AAs[sheets[i].AACurIndx].chr3<<" "<<AAs[sheets[i].AACurIndx].chain<<setw(4)<<AAs[sheets[i].AACurIndx].num
							<<AAs[sheets[i].AACurIndx ].resInsertion <<" "<<sheets[i].prevAtomName <<AAs[sheets[i].AAPrevIndx ].chr3 <<" "
							<<AAs[sheets[i].AAPrevIndx ].chain <<setw(4)<<AAs[sheets[i].AAPrevIndx ].num<<AAs[sheets[i].AAPrevIndx ].resInsertion;

					}
					out<<endl;

					++sCounter;
				//}
				//the 2 cases left r to be written later
			}
		}

		/*****
				write Atoms
		*****/
		int atomsCounter = 1;
		for(i=startAARank-1;i<endAARank;i++)
		{
			for (j=0;j<numOfAtoms(i);j++)
			{

				out<<"ATOM"<<setw(7)<<atomsCounter<<" "<<setw(4)<<AAs[i].atoms[j].name<<AAs[i].atoms [j].locIndicator<<setw(3)<<right<<AAs[i].chr3
					<<setw(2)<<right<<AAs[i].chain<<setw(4)<<right<<AAs[i].num<<AAs[i].resInsertion<<setw(3)
					<<" "<<setw(8)<<right<<setiosflags(ios::fixed) << setprecision(3)<<AAs[i].atoms[j].coord.x<<setw(8)<<right<<AAs[i].atoms[j].coord.y
					<<setw(8)<<right<<AAs[i].atoms[j].coord.z<<setw(6)<<right<<AAs[i].atoms[j].occupancy<<setw(6)<<right<<AAs[i].atoms[j].tempFactor
					<<setw(10)<<" "<<setw(2)<<right<<AAs[i].atoms[j].type <<setw(2)<<AAs[i].atoms[j].charge<<endl;
				atomsCounter++;
			}
		}

		out.close();
	}
	else
	{
		string eMsg = "The given range (";
		eMsg += toString(startAARank);
		eMsg += ", ";
		eMsg += toString(endAARank);
		eMsg += ") is incorrect..or the portion os empty (no.AA= ";
		eMsg += toString(numOfAA());
		eMsg += ")";

		errMsg("Protein", "writePDB", eMsg, true);

	}
}
////////////////////////////////////////////////////////////////////////////////
void Protein::writeSCModel(string outFile,int startAARank, int endAARank)
//the written file will be written as a simplified model...the side chain as a line
// The whole Protein could be written to a file or a portion of it
//startAARank is the rank of the first AA to be written to outfile...for example 4 is the 4th AA
//endAARank is the rank of the last AA would be written to the outfile PDB file
{
	int tmpNumOfAA = numOfAA();

	//assure from the ranks given
	if ((startAARank>=1) && (startAARank <= tmpNumOfAA) && (endAARank <= tmpNumOfAA))
	{
		ofstream out;
		out.open(outFile.c_str());
		if (!out)
		{
			errMsg("Protein", "writeSCModel", "Unable to open "+outFile);
		}
		else
		{
			int atomsCounter = 1;
			for (int i=startAARank - 1; i<endAARank; i++)
			{
				if (AAs[i].whtInCoord != 'S')
					setScEndPoint(i);				//set the side chain end point using heavy atoms
				/*****
						write just backbone atoms
				*****/
				for (int j=0; j<numOfAtoms(i); j++)
				{

					if (!AAs[i].atoms [j].isSideChain)
					{
						out<<"ATOM"<<setw(7)<<atomsCounter<<" "<<setw(4)<<AAs[i].atoms[j].name<<AAs[i].atoms [j].locIndicator<<setw(3)<<right<<AAs[i].chr3
							<<setw(2)<<right<<AAs[i].chain<<setw(4)<<right<<AAs[i].num<<AAs[i].resInsertion<<setw(3)
							<<" "<<setw(8)<<right<<setiosflags(ios::fixed) << setprecision(3)<<AAs[i].atoms[j].coord.x<<setw(8)<<right<<AAs[i].atoms[j].coord.y
							<<setw(8)<<right<<AAs[i].atoms[j].coord.z<<setw(6)<<right<<AAs[i].atoms[j].occupancy<<setw(6)<<right<<AAs[i].atoms[j].tempFactor
							<<setw(10)<<" "<<setw(2)<<right<<AAs[i].atoms[j].type <<setw(2)<<AAs[i].atoms[j].charge<<endl;
						atomsCounter++;
					}
				}
				/*****
						write the line represents the side chain
				*****/
				out<<"ATOM"<<setw(7)<<atomsCounter<<" "<<setw(4)<<" CB "<<' '<<setw(3)<<right<<AAs[i].chr3
					<<setw(2)<<right<<AAs[i].chain<<setw(4)<<right<<AAs[i].num<<AAs[i].resInsertion<<setw(3)
					<<" "<<setw(8)<<right<<setiosflags(ios::fixed) << setprecision(3)<<AAs[i].ScEndPoint .x<<setw(8)<<right<<AAs[i].ScEndPoint .y
					<<setw(8)<<right<<AAs[i].ScEndPoint .z<<setw(6)<<right<<AAs[i].atoms[0].occupancy<<setw(6)<<right<<AAs[i].atoms[0].tempFactor
					<<setw(10)<<" "<<setw(2)<<right<<AAs[i].atoms[0].type <<setw(2)<<AAs[i].atoms[0].charge<<endl;
				atomsCounter++;

			}
		}

	}
	else
	{
		string eMsg = "The given range (";
		eMsg += toString(startAARank);
		eMsg += ", ";
		eMsg += toString(endAARank);
		eMsg += "is incorrect..or the portion is empty (no.AA= ";
		eMsg += toString(numOfAA());
		eMsg += ")";

		errMsg("Protein", "writeSCModel", eMsg, true);
	}
}
////////////////////////////////////////////////////////////////////////////
void Protein::writeAAInfo(string outFile)
{
	ofstream out;

	out.open (outFile.c_str ());

	if (!out)
	{
		errMsg("Protein", "writeAAInfo", "Unable to open " +outFile, true);
	}
	else
	{
		out<<"AA#  AA  chain SeqNum   RG    RG.x    RG.y    RG.z   WhtInCoord   wht.x   wht.y   wht.z   phi     psi    #ofAtoms"<<endl;
		out<<"==============================================================================================================="<<endl;
		for (int i=0;i<numOfAA();i++)
			out<<setw(4)<<i+1<<setw(4)<<AAs[i].chr3<<" "<<setw(5)<<AAs[i].chain<<" "<<setw(6)<<AAs[i].num<<" "<<setprecision(3)<<setw(5)<<AAs[i].gyration
				<<" "<<setw(6)<<AAs[i].gyCoord.x<<"  "<<setw(6)<<AAs[i].gyCoord.y<<"  "<<setw(6)<<AAs[i].gyCoord.z<<" "<<setw(7)<<AAs[i].whtInCoord
				<<setw(7)<<" "<<setw(6)<<AAs[i].coord .x <<"  "<<setw(6)<<AAs[i].coord .y <<"  "<<setw(6)<<AAs[i].coord .z <<" "<<setw(6)<<AAs[i].angles.phi
				<<" "<<setw(6)<<AAs[i].angles.psi<<setw(3)<<" "<<setw(8)<<numOfAtoms(i)<<endl;
	}
}
////////////////////////////////////////////////////////////////////////////
//number of AA in the whole portion
inline int Protein::numOfAA()
{
	return AAs.size();
}
/////////////////////////////////////////////////////////////////////////////
//number of atoms in the whole portion
inline int Protein::numOfAtoms()
{
	int numOfAtoms = 0;
	for (int i=0;i<AAs.size();i++)
		numOfAtoms += AAs[i].atoms.size();

	return numOfAtoms;
}
///////////////////////////////////////////////////////////////////////////////
//number of atoms for a particular AA
inline int Protein::numOfAtoms(int AAIndx)
{
	if ((AAIndx>=0) && (AAIndx<numOfAA()))
		return AAs[AAIndx].atoms.size();
	else
	{
		string eMsg = "The given index (";
		eMsg += toString (AAIndx);
		eMsg += ") is out of range..or the portion is empty (no.AAs= ";
		eMsg += toString(numOfAA());
		eMsg += ")";

		errMsg("Protein", "numOfAtoms", eMsg);
	}

	return 0;
}
//////////////////////////////////////////////////////////////////////////////////////////
//atomsincluded: H heavy atoms, A all
int Protein::numOfSCAtoms(int AAIndx, char atomsIncluded)
{
	if ((AAIndx >= 0) && (AAIndx < numOfAA()))
	{
		int cntr=0, i;
		/*****
				look for Heavy atoms
		*****/
		if (atomsIncluded == 'H')
			for (i=0; i<numOfAtoms(AAIndx); i++)
				if ((AAs[AAIndx].atoms[i].isSideChain) && (isHeavyAtom(AAs[AAIndx].atoms[i].type)))
					cntr++;
		/*****
				look for All atoms in the side chain (H atom will be included)
		*****/
		if (atomsIncluded == 'A')
			for (i = 0; i<numOfAtoms(AAIndx); i++)
				if (AAs[AAIndx].atoms[i].isSideChain)
					cntr++;
		return cntr;
	}
	else
	{
		string eMsg = "AAindx given (";
		eMsg += toString(AAIndx);
		eMsg += "is out of range..";

		errMsg("Protein", "numOfSCAtoms", eMsg);
	}

	return 0;
}
////////////////////////////////////////////////////////////////////////////
//given AA sequence number...return the indx in AAs
int Protein::getAAIndx(int AANum)
{
	for (int i=0;i<numOfAA();i++)
		if (AAs[i].num == AANum)
			return i;

	//Not Found
	return -1;
}
///////////////////////////////////////////////////////////////////////////////
//given AA indx in AAs and atom name (or a substring of atom name)...return atom indx in AAs.atoms
int Protein::getAtomIndx(int AAIndx, string atomName)
{

	if ((AAIndx>=0) && (AAIndx < numOfAA()))
	{
		for (int i=0;i<numOfAtoms(AAIndx);i++)
		{
			if (AAs[AAIndx].atoms[i].name.find(atomName) != AAs[AAIndx].atoms[i].name.npos)
				return i;		//The indx of the atom
		}
	}
	else
	{
		string eMsg = "The given index (";
		eMsg += toString(AAIndx);
		eMsg += ") is out of range..or the portion is empty (No.AAs= ";
		eMsg += toString(numOfAA());
		eMsg += ")";

		errMsg("Protein", "getAtomIndx", eMsg);
	}

	//Not found
	return -1;

}

int Protein::getAtomIndx2(int AAIndx, string atomName)
{

	if ((AAIndx>=0) && (AAIndx < numOfAA()))
	{
		return numOfAtoms(AAIndx)-1;		//The indx of the atom
    }
	else
	{
		string eMsg = "The given index (";
		eMsg += toString(AAIndx);
		eMsg += ") is out of range..or the portion is empty (No.AAs= ";
		eMsg += toString(numOfAA());
		eMsg += ")";

		errMsg("Protein", "getAtomIndx", eMsg);
	}

	//Not found
	return -1;

}
////////////////////////////////////////////////////////////////////////////////////
//() overloading
AminoAcid Protein::operator () (int AARank)
{
	return AAs[AARank-1];
}
/////////////////////////////////////////////////////////////////////////////////////
//given the name of the torsion u gonna compute (phi or psi) and the indx of AA in the AAs...return the torsion angle
//if the angle is phi and the AAIndx is 0 (first one) then the computed torsion is with origin (0, 0, 0) instead of Ci-1
//if the angle is psi and the AAIndx is numOfAA()-1 (last one) then the computed torsion is with origin instead of Ni+1
double Protein::getTorsion(string angleName, int AAIndx)
{
	//Coordinate p1,p2,p3,p4;
	Coordinate origin;
	int p1, p2, p3, p4; //the indeces of the points

	int tmpNumOfAA = numOfAA();

	if ((AAIndx>=0) && (AAIndx < tmpNumOfAA))
	{
		if (angleName == "phi")
		{
			//else
			//p1 is 0 0 0 (origin)
			//p2 = getAtomCoordinate(AAIndx," N  ");
			p2 = getAtomIndx(AAIndx, " N  ");
			//p3 = getAtomCoordinate(AAIndx," CA ");
			p3 = getAtomIndx(AAIndx, " CA ");
			//p4 = getAtomCoordinate(AAIndx," C  ");
			p4 = getAtomIndx(AAIndx, " C  ");

			if (AAIndx>0)		//it is not the first AA
			{
				p1 = getAtomIndx(AAIndx - 1, " C  ");				//from previoud AA
				if ((p1 != -1) && (p2 != -1) && (p3 != -1) && (p4 != -1))
					return getTorsionAngle(AAs[AAIndx-1].atoms [p1].coord,
											AAs[AAIndx].atoms [p2].coord,
											AAs[AAIndx].atoms [p3].coord,
											AAs[AAIndx].atoms [p4].coord);
				else
					return 999.0;
			}
			else
			{
				if ((p2 != -1) && (p3 != -1) && (p4 != -1))
					return getTorsionAngle(origin,
											AAs[AAIndx].atoms [p2].coord,
											AAs[AAIndx].atoms [p3].coord,
											AAs[AAIndx].atoms [p4].coord);
				else
					return 999.0;
			}

		}
		else
			if (angleName == "psi") {

				//p1 = getAtomCoordinate(AAIndx," N  ");
				p1 = getAtomIndx(AAIndx, " N  ");
				//p2 = getAtomCoordinate(AAIndx," CA ");
				p2 = getAtomIndx(AAIndx, " CA ");
				//p3 = getAtomCoordinate(AAIndx," C  ");
				p3 = getAtomIndx(AAIndx, " C  ");

				if (AAIndx < tmpNumOfAA - 1)		//it is not the last AA
				{
					//p4 = getAtomCoordinate(AAIndx + 1," N  ");		//next AA
					p4 = getAtomIndx(AAIndx+1, " N  ");			//next AA
					if ((p1 != -1) && (p2 != -1) && (p3 != -1) && (p4 != -1))
						return getTorsionAngle(AAs[AAIndx].atoms [p1].coord,
												AAs[AAIndx].atoms [p2].coord,
												AAs[AAIndx].atoms [p3].coord,
												AAs[AAIndx+1].atoms [p4].coord);
					else
						return 999.0;
				}
				else
				{
					if ((p1 != -1) && (p2 != -1) && (p3 != -1))
						return getTorsionAngle(AAs[AAIndx].atoms [p1].coord,
												AAs[AAIndx].atoms [p2].coord,
												AAs[AAIndx].atoms [p3].coord,
												origin);
					else
						return 999.0;

				}
			}
			else
			{
				string eMsg = "The name of the torsion angle sent (" + angleName;
				eMsg += ") is wrong...plz choose phi or psi.";

				errMsg("Protein", "getTorsion", eMsg, true);

			}
	}
	else
	{
		string eMsg = "The given index (";
		eMsg += toString(AAIndx);
		eMsg += ") is out of range..or the portion is empty (No.AAs= ";
		eMsg += toString(numOfAA());
		eMsg += ")";

		errMsg("Protein", "getTorsion", eMsg, true);
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//takes the indx of AA (AAs indx) and the chi number...e.x. 1 means chi1..2 means chi2 ....11 means chi11..12 means chi12
//returns the chi angle
/*////////////////////////////
phi main chain torsion angle for atoms C,N,CA,C.
psi main chain torsion angle for atoms N,CA,C,N.
omega main chain torsion angle for atoms CA,C,N,CA.
alpha virtual torsion angle between consecutive C-alpha atoms.

chi1 side chain torsion angle for atoms N,CA,CB,*G.
	chi11 side chain torsion angle for atoms N,CA,CB,*G1.
	chi12 side chain torsion angle for atoms N,CA,CB,*G2.
chi2 side chain torsion angle for atoms CA,CB,*G,*D.
	chi21 side chain torsion angle for atoms CA,CB,*G,*D1.
	chi22 side chain torsion angle for atoms CA,CB,*G,*D2.
chi3 side chain torsion angle for atoms CB,*G,*D,*E.
	chi31 side chain torsion angle for atoms CB,*G,*D,*E1.
	chi32 side chain torsion angle for atoms CB,*G,*D,*E2.
chi4 side chain torsion angle for atoms *G,*D,*E,*Z.
	chi51 side chain torsion angle for atoms *D,*E,*Z, NH1.
	chi52 side chain torsion angle for atoms *D,*E,*Z, NH2.
*///
// in the case u send for example chi 1 and the function could not find *G it will computes Chi11 instead, if could not then chi12
double Protein::getChi(int AAIndx, int chiNum){

	Coordinate N,CA,CB;
	int _G, _D, _E, _Z, NH1, NH2;			//the indeces of such atoms

	int Nindx = getAtomIndx(AAIndx, " N  ");
	int CAindx = getAtomIndx (AAIndx, " CA ");
	int CBindx = getAtomIndx(AAIndx, " CB ");
	if ((Nindx != -1) && (CBindx != -1) && (CAindx != -1))
	{
		N  = AAs[AAIndx].atoms [Nindx].coord ;
		CA = AAs[AAIndx].atoms [CAindx].coord ;
		CB = AAs[AAIndx].atoms [CBindx].coord ;

		if ((AAIndx>=0) && (AAIndx<numOfAA()))
		{

			/*****
					To check that the chi angle you requested is exist within the amino acid
					for example u may request chi 21 but the aa contains just chi 1
					I think that number of chi's in AA means the number of rotatable chi
					for example PHE: the number of chi is 2 (the only rotatable angles) but it still has CD, and CZ atoms that form chi3 and chi 4
			*****/
			int chiCategory;
			if (chiNum>10)
			{
				chiCategory = chiNum/10;
				if (chiCategory > getNumOfChi(AAs[AAIndx].chr1))
					return 999.00;
			}
			else
				if (getNumOfChi(AAs[AAIndx].chr1)<chiNum)
					return 999.0;

			//chi1 or chi11 or chi12
			if ((chiNum == 1) || (chiCategory == 1))
			{
				switch (chiNum)
				{
					case 1:
							_G = getAtomIndx(AAIndx, "G");
							if (_G != -1)
								return getTorsionAngle(N,CA,CB,AAs[AAIndx].atoms [_G].coord);
							else
								return 999;  //in the case of error
					case 11: _G = getAtomIndx(AAIndx, "G1");
							if (_G != -1)
								return getTorsionAngle(N,CA,CB,AAs[AAIndx].atoms [_G].coord);
							else
								return 999;  //in the case of error
					case 12: _G = getAtomIndx(AAIndx, "G2");
							if (_G != -1)
								return getTorsionAngle(N,CA,CB,AAs[AAIndx].atoms[_G].coord);
							else
								return 999;  //in the case of error
					default :	return 999.0;
				}
			}

			//chi 2 or chi21 or chi22
			if ((chiNum == 2) || (chiCategory == 2))
			{
				switch (chiNum)
				{
					case 2: _G = getAtomIndx(AAIndx, "G");
							_D = getAtomIndx(AAIndx, "D");
							if ((_G != -1) && (_D != -1))
								return getTorsionAngle(CA,CB,AAs[AAIndx].atoms [_G].coord ,AAs[AAIndx].atoms [_D].coord );
							else
								return 999;  //in the case of error
					case 21: _G = getAtomIndx(AAIndx, "G");
							_D = getAtomIndx(AAIndx, "D1");
							if ((_G != -1) && (_D != -1))
								return getTorsionAngle(CA,CB,AAs[AAIndx].atoms [_G].coord ,AAs[AAIndx].atoms [_D].coord);
							else
								return 999;  //in the case of error
					case 22: _G = getAtomIndx(AAIndx, "G");
							_D = getAtomIndx(AAIndx, "D2");
							if ((_G != -1) && (_D != -1))
								return getTorsionAngle(CA,CB,AAs[AAIndx].atoms [_G].coord ,AAs[AAIndx].atoms [_D].coord);
							else
								return 999;  //in the case of error
					default :	return 999.0;
				}
			}

			//chi3 or chi31 or chi32
			if ((chiNum == 3) || (chiCategory == 3))
			{
				switch (chiNum)
				{
					case 3: _G = getAtomIndx(AAIndx, "G");
							_D = getAtomIndx(AAIndx, "D");
							_E = getAtomIndx(AAIndx, "E");
							if ((_G != -1) && (_D != -1) && (_E != -1))
								return getTorsionAngle(CB,AAs[AAIndx].atoms [_G].coord ,AAs[AAIndx].atoms [_D].coord,AAs[AAIndx].atoms [_E].coord );
							else
								return 999;  //in the case of error
					case 31: _G = getAtomIndx(AAIndx, "G");
							_D = getAtomIndx(AAIndx, "D");
							_E = getAtomIndx(AAIndx, "E1");
							if ((_G != -1) && (_D != -1) && (_E != -1))
								return getTorsionAngle(CB,AAs[AAIndx].atoms [_G].coord ,AAs[AAIndx].atoms [_D].coord,AAs[AAIndx].atoms [_E].coord );
							else
								return 999;  //in the case of error
					case 32: _G = getAtomIndx(AAIndx, "G");
							_D = getAtomIndx(AAIndx, "D");
							_E = getAtomIndx(AAIndx, "E2");
							if ((_G != -1) && (_D != -1) && (_E != -1))
								return getTorsionAngle(CB,AAs[AAIndx].atoms [_G].coord ,AAs[AAIndx].atoms [_D].coord,AAs[AAIndx].atoms [_E].coord);
							else
								return 999;  //in the case of error
					default :	return 999.0;
				}
			}

			//chi4
			if (chiNum == 4)
			{
				 _G = getAtomIndx(AAIndx, "G");
				_D = getAtomIndx(AAIndx, "D");
				_E = getAtomIndx(AAIndx, "E");
				_Z = getAtomIndx(AAIndx, "Z");
				if ((_G != -1) && (_D != -1) && (_E != -1) && (_Z != -1))
					return getTorsionAngle(AAs[AAIndx].atoms [_G].coord ,AAs[AAIndx].atoms [_D].coord,AAs[AAIndx].atoms [_E].coord,AAs[AAIndx].atoms [_Z].coord);
				else
					return 999;  //in the case of error
			}

			//chi51 or chi52
			if (chiCategory == 5)
			{
				switch (chiNum)
				{
					case 51: _G = getAtomIndx(AAIndx, "G");
							_D = getAtomIndx(AAIndx, "D");
							_E = getAtomIndx(AAIndx, "E");
							NH1 = getAtomIndx(AAIndx, "NH!");
							if ((_G != -1) && (_D != -1) && (_E != -1) && (NH1 != -1))
								return getTorsionAngle(AAs[AAIndx].atoms [_G].coord ,AAs[AAIndx].atoms [_D].coord,AAs[AAIndx].atoms [_E].coord,AAs[AAIndx].atoms [NH1].coord);
							else
								return 999;  //in the case of error
					case 52: _G = getAtomIndx(AAIndx, "G");
							_D = getAtomIndx(AAIndx, "D");
							_E = getAtomIndx(AAIndx, "E");
							NH2 = getAtomIndx(AAIndx, "NH2");
							if ((_G != -1) && (_D != -1) && (_E != -1) && (NH2 != -1))
								return getTorsionAngle(AAs[AAIndx].atoms [_G].coord ,AAs[AAIndx].atoms [_D].coord,AAs[AAIndx].atoms [_E].coord,AAs[AAIndx].atoms [NH2].coord);
							else
								return 999;  //in the case of error
					default :	return 999.00;
				}
			}

		}
	}
	else
		return 999.0;


}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//method V = VDW;  C= calculated; O= Covalent
//find the self collision
int Protein::doesCollide (double cutOff)
{

	int tmpNumOfAA = numOfAA();
	if ((tmpNumOfAA) && (AAsCollide.size() == 0))
	{
		double  distance   = 0.0,
				frstRadius = 0.0,
				scndRadius = 0.0;


		for (int i=0;i<tmpNumOfAA; i++)
		{
			int tmpNumOfAtoms = numOfAtoms(i);
			for (int j=0;j<tmpNumOfAtoms; j++)
			{
				Coordinate jCoordinate = AAs[i].atoms [j].coord;

				frstRadius = getRadius(AAs[i].atoms[j].type, 'C');
				/*****
						the rest AA's
				*****/
				for (int m=i+1; m<tmpNumOfAA; m++)
				{
					/*****
							The rest of atoms from the next AA after atom N
					*****/
					for (int l=1;l<numOfAtoms(m);l++)
					{
						distance = getDistance(AAs[m].atoms[l].coord, jCoordinate);
						scndRadius = getRadius(AAs[m].atoms[l].type, 'C');
						if ( frstRadius + scndRadius - distance >= cutOff)
						{
							AAsCollide.push_back(i);		//the indx of first AA collide
							AAsCollide.push_back (m);		//the Indx of the second AA collide ..
							atomsCollide.push_back(j);		//the indx of the first atom collide
							atomsCollide.push_back (l);		//the indx of the second atom collide

							//cout<<"Self Collision found b/w "<<AAs[i].chr3<<" "<<AAs[i].num<<" : "<<AAs[i].atoms[j].name<<" and "<<AAs[m].chr3<<" "<<AAs[m].num<<" : "<<AAs[m].atoms[l].name<<"  distance = "<<distance<<endl;
							return 1;
						}

					}
				}


			}

		}

	}
	else
	{
		string eMsg = "The portion is empty (No.AAs= ";
		eMsg += toString(numOfAA());
		eMsg += ") or a list of collide AAs already computed in AAsCollide and atomsCollide.";

		errMsg("Protein", "doesCollide", eMsg);

	}

	return 0;		//no collision was found

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Protein::VDW(double VDW_radius1, double VDW_radius2, double dist)
{
	double Rij = VDW_radius1 + VDW_radius2;

	if (dist >= Rij)
		return 0;
	if (dist < 0.8254 * Rij)
		return 10;
	if ((0.8254 * Rij < dist) && (dist < Rij))
		return 57.272 * (1 - dist/Rij);

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Protein::isSSAA(int AAIndx)
{
	if ((AAIndx >= 0) && (AAIndx < numOfAA()))
	{
		if ((AAs[AAIndx].SStype == 'H') || (AAs[AAIndx].SStype == 'S'))
			return true;
	}
	else
	{
		string eMsg = "The given index (";
		eMsg += toString(AAIndx);
		eMsg += ") is out of range..or the portion is empty (No.AAs= ";
		eMsg += toString(numOfAA());
		eMsg += ")";

		errMsg("Protein", "isSSAA", eMsg);
	}

	return false; //not a SS Amino Acid
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Protein::inHlces (int aaIndx){
	if ((aaIndx >= 0) && (aaIndx < numOfAA()))
	{
		if (AAs[aaIndx].SStype == 'H'){
			for (int i=0; i<hlces.size (); i++)
				if (hlces[i].startIndx <= aaIndx && hlces[i].endIndx >= aaIndx)
					return i;
		}
		else return -1;
	}
	else
	{
		string eMsg = "The given index (";
		eMsg += toString(aaIndx);
		eMsg += ") is out of range..or the portion is empty (No.AAs= ";
		eMsg += toString(numOfAA());
		eMsg += ")";

		errMsg("Protein", "inHlces", eMsg);
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//given a range of the portion...return the sequence
string Protein::getSequence(int startIndx, int endIndx)
{
	string proteinSeq = "";
	if ((startIndx>=0) && (endIndx < numOfAA()))
	{
		for (int i= startIndx;i<=endIndx;i++)
			proteinSeq = proteinSeq + AAs[i].chr1;
	}
	else
	{
		string eMsg = "The given range (";
		eMsg += toString(startIndx);
		eMsg += ", ";
		eMsg += toString(endIndx);
		eMsg += ") is out of range..or the portion is empty (No.AAs= ";
		eMsg += toString(numOfAA());
		eMsg += ")";

		errMsg("Protein", "getSequence", eMsg);
	}

	return proteinSeq;

}
//////////////////////////////////////////////////////////////////////////////////////////////
//given a range of the portion...return the sequence with missing AAs included
string Protein::getSequenceWmissing (int startIndx, int endIndx)
{
	string proteinSeq = "";
	if ((startIndx>=0) && (endIndx < numOfAA()))
	{
		if (missingAAs.size () == 0)
			return getSequence(startIndx, endIndx);
		else
		{

			int missingCntr = 0, i = startIndx;

			while (i<=endIndx)
			{
				if ((missingCntr < missingAAs.size ()) && (AAs[i].num > missingAAs[missingCntr].num))
					proteinSeq += missingAAs[missingCntr++].chr1;
				else

					proteinSeq += AAs[i++].chr1;
			}

			//list the rest of missing AAs
			while (missingCntr < missingAAs.size ())
				proteinSeq += missingAAs[missingCntr++].chr1;

		}

		return proteinSeq;
	}
	else
	{
		string eMsg = "The given range (";
		eMsg += toString(startIndx);
		eMsg += ", ";
		eMsg += toString(endIndx);
		eMsg += ") is out of range..or the portion is empty (No.AAs= ";
		eMsg += toString(numOfAA());
		eMsg += ")";

		errMsg("Protein", "getSequenceWmissing", eMsg);

	}

}
///////////////////////////////////////////////////////////////////////////////////////
//given the indx of AA in AAs and the name of the atom in that AA...return the coordinate of that atom
Coordinate Protein::getAtomCoordinate(int AAIndx, string atomName)
{
	Coordinate coord;

	int indx = getAtomIndx(AAIndx,atomName);
	if ( indx != -1 ) //-1 means not found
		return AAs[AAIndx].atoms[indx].coord;
	else
	{
		string eMsg = "No such atom ( " + atomName;
		eMsg += ") found ";
		if (AAIndx >=0 && AAIndx < numOfAA()){

			eMsg += "in AA (";
			eMsg += AAs[AAIndx].chr3 ;
			eMsg += "-";
			eMsg += toString(AAs[AAIndx].num);
			eMsg += " )";

			errMsg("Protein", "getAtomCoordinate", eMsg);

		}
		else{

			eMsg += "...pad AA indx (";
			eMsg += toString(AAIndx);
			eMsg += ").";

			errMsg("Protein", "getAtomCoordinate", eMsg);
		}

		coord.x = -999.0;
		coord.y = -999.0;
		coord.z = -999.0;

	}

	//Not Found....return the initial values 0 0 0
	return coord;
}
///////////////////////////////////////////////////////////////////////////////////////
void Protein::fillMissingAAs ()
{


	if ((missingAAs.size () == 0) && (numOfAA()))		//fill in the case it is not filled before
	{
		bool stop = false;
		int i = 0;
		size_t found;
		while ((i<header.size()) && (!stop))
		{
			found = header[i].find("MISSING RESIDUES");		//REMARK 465 referes to missing residues
			if ( found != header[i].npos)
				stop = true;
			i++;
		}

		if (stop)		//found some missing residues
		{
			found = header[i].find("RES C SSSEQI");
			while (header[i].find("RES C SSSEQI") == header[i].npos)
				i++;
			i++;

			missingAA tmpmissAA;
			string model;
			for ( int cntr = i; (cntr<header.size() && header[cntr].substr (0, 10) == "REMARK 465"); cntr++)
			{
				if (header[cntr].substr (19,1) == AAs[0].chain)		//should be in the same chain
				{
					if (missingAAs.size () == 0)
						model = header[cntr].substr (13, 1).c_str ();

					if (header[cntr].substr (13,1).c_str () == model)		//if in the same model
					{
						tmpmissAA.chr3 = header[cntr].substr(15, 3).c_str ();
						tmpmissAA.chr1 = chr3ToChr1(tmpmissAA.chr3);
						tmpmissAA.num  = atoi(header[cntr].substr (22, 4).c_str ());
						tmpmissAA.resInsertion = header[cntr][26];


						bool checkStrands = true;

						//check helices
						int SScntr = 0;
						while(SScntr < hlces.size())
						{
							if ((tmpmissAA.num >= AAs[hlces[SScntr].startIndx].num) &&
								(tmpmissAA.num <= AAs[hlces[SScntr].endIndx].num))
							{
								checkStrands = false;
								tmpmissAA.SStype = 'H';
								break;

							}
							SScntr++;
						}
						if (checkStrands)
						{
							SScntr = 0;
							while (SScntr < sheets.size ())
							{
								if ((tmpmissAA.num >= AAs[sheets[SScntr].startIndx].num) &&
									(tmpmissAA.num <= AAs[sheets[SScntr].endIndx].num))
								{
									tmpmissAA.SStype = 'S';
									break;
								}
								SScntr++;
							}
						}
						missingAAs.push_back (tmpmissAA);
					}
				}
			}
		}
	}
}
//////////////////////////////////////////////////////////////////////////////////////
//remove all H atoms from the portion
void Protein::removeHAtoms()
{

	int i,
		j,
		numOfRemovedAtoms;

	for (i=0;i<numOfAA();i++)
	{
		numOfRemovedAtoms = 0;		//number of H atoms removed so far
		j = 0;
		while (j<numOfAtoms(i) + numOfRemovedAtoms)
		{
			/*****
					The type of atom is considered valid all the time in AAs[].atoms.type
			*****/
			if (AAs[i].atoms[j - numOfRemovedAtoms].type == 'H')
			{
				AAs[i].atoms.erase(AAs[i].atoms.begin() + j - numOfRemovedAtoms);
				numOfRemovedAtoms++;
			}
			j++;
		}

		/*****
				if some of H atoms removed, in general wheather the removed atom is from SC or from backbone...reset some variables
		*****/
		if (numOfRemovedAtoms)
		{
			//reset the radius of gyration
			AAs[i].gyration = 0;

			//reset whtInCoord flag
			AAs[i].whtInCoord = 'N';

			//reset the two atoms collide if one of them at least was 'H' atom
			if (atomsCollide.size ())
				if ((AAs[AAsCollide[0]].atoms [atomsCollide[0]].type == 'H') || (AAs[AAsCollide[1]].atoms [atomsCollide[1]].type == 'H'))
				{
					//remove the two atoms collide
					atomsCollide.clear ();
					//remove the two AA's collide ...
					AAsCollide.clear ();
				}


		}
	}

}
/////////////////////////////////////////////////////////////////////////////////////
void Protein::renameAA(int AAIndx, string newName)
{
	/*****
			To be done later.......
			the newName should be checked for validity first
	*****/
	if ((AAIndx>=0) && (AAIndx < numOfAA()))
	{
		AAs[AAIndx].chr3 = newName;
		AAs[AAIndx].chr1 = chr3ToChr1(newName);
	}
	else
	{
		string eMsg = "The given index (";
		eMsg += toString(AAIndx);
		eMsg += ") is out of range..or the portion is empty (No.AAs= ";
		eMsg += toString(numOfAA());
		eMsg += ")";

		errMsg("Protein", "renameAA", eMsg);
	}
}
/////////////////////////////////////////////////////////////////////////////////////
void Protein::deleteAA(int AAIndx)
{
	if ((AAIndx>=0) && (AAIndx < numOfAA()))
	{
		AAs.erase(AAs.begin() + AAIndx);
		/*****
				check if the current AA deleted is a member of secondary structure and update accordingly
		*****/
		int i;
		for (i=0;i<hlces.size(); i++)
		{
			if ((AAIndx <= hlces[i].endIndx) && (AAIndx > hlces[i].startIndx))		//the AA deleted is in the middle or the last one in a hlx
			{
				hlces[i].endIndx -= 1;
				continue;
			}

			if (AAIndx == hlces[i].startIndx )		//the AA deleted is the first AA in the hlx
			{
				hlces[i].startIndx += 1;
				continue;
			}

			cout<<i<<endl;
		}

		for (i=0; i<sheets.size ();i++)
		{
			if ((AAIndx <= sheets[i].endIndx ) && (AAIndx > sheets[i].startIndx ))
			{
				sheets[i].endIndx -= 1;
				continue;
			}

			if (AAIndx == sheets[i].startIndx )
			{
				sheets[i].startIndx += 1;
				continue;
			}

		}

	}
	else
	{
		string eMsg = "The given index (";
		eMsg += toString(AAIndx);
		eMsg += ") is out of range..or the portion is empty (No.AAs= ";
		eMsg += toString(numOfAA());
		eMsg += ")";

		errMsg("Protein", "deleteAA", eMsg);
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////
void Protein::deleteAtom (int AAIndx, string atomName)
{
	if ((AAIndx>=0) && (AAIndx < numOfAA()))
	{
		int atomIndx = getAtomIndx(AAIndx, atomName);
		if ( atomIndx!= -1)
		{
			AAs[AAIndx].atoms.erase (AAs[AAIndx].atoms .begin () + atomIndx);
			//if it is in side chain and the center of side chain is already computed....reset it
			if ((AAs[AAIndx].atoms [atomIndx].isSideChain)  && (AAs[AAIndx].whtInCoord == 'S'))
			{
				AAs[AAIndx].whtInCoord = 'N';
				//reset radius of gyration also
				AAs[AAIndx].gyration = 0;
			}
			//if it is in the back bone and the center of whole AA was computed . . . reset it
			if ((!AAs[AAIndx].atoms [atomIndx].isSideChain ) && (AAs[AAIndx].whtInCoord == 'A'))
				AAs[AAIndx].whtInCoord = 'N';

		}
	}
	else
	{
		string eMsg = "The given index (";
		eMsg += toString(AAIndx);
		eMsg += ") is out of range..or the portion is empty (No.AAs= ";
		eMsg += toString(numOfAA());
		eMsg += ")";

		errMsg("Protein", "deleteAtom", eMsg);

	}
}
/////////////////////////////////////////////////////////////////////////////////////////
//Add a list of AAs (specified range of AAs from another portion) to the end of the portion
//this function does not append secondary structure information (hlces and sheets data structures)...for example if the added AA's are secondary
//structures AA's then that will be mentioned just in SStype attribute.
//skip is the size of the gap in sequence number u wanna leave b/w the current portion and the appended one, defualt is 0
void Protein::append(Protein AAList,int startIndx, int endIndx, int skip)
{
	int listSize = AAList.numOfAA();

	if ((startIndx >= 0) && (endIndx < listSize))	//if the AA list is not empty and the range sent is valid
	{


/*
		//if less than 0...make it starts from the beginning
		if (startIndx < 0)
			startIndx = 0;

		//if exceeds the size...make it ends @ the end
		if (endIndx >= listSize)
			endIndx = listSize - 1;
*/
		string newPath = AAList.path;		//take the source path if the current portion is empty
		string newID = AAList.ID;			//take the source ID if the current portion is empty

		int currentListSize = numOfAA();
		int AASeqNum = 0;		//the last AA Seq number
		string newChain = AAList.AAs[0].chain;		//the new chain is the same chain of the comming AA's if the current portion is empty

		if (currentListSize)
		{
			newPath  = path;		//keep the original path
			newID	 = ID;			//keep the original ID
			AASeqNum = AAs[currentListSize - 1].num;	//if the list already has AA's ...then remember the last seq number
			newChain = AAs[0].chain;					//if the list already has AA's...then the new chain is the same chian
		}

		int numOfAAAdded = 0;		//number of AA's have been added so far

		//add skip value to sequence number
		AASeqNum += skip;

		path = newPath;		//assign path
		ID   = newID;		//assign ID

		for (int i=startIndx;i<=endIndx;i++)
		{
			AAs.push_back(AAList.AAs[i]);								//insert the AA at the end of current portion
			AAs[currentListSize + numOfAAAdded].chain = newChain;		//set the chain
			AAs[currentListSize + numOfAAAdded++].num = ++AASeqNum;		//set the sequence number
		}
	}
	else
	{
		string eMsg = "The given portion to be appended is empty or the given range (";
		eMsg += toString(startIndx);
		eMsg += "-";
		eMsg += toString(endIndx);
		eMsg += " is incorrect.";

		errMsg("Protein", "append", eMsg);


	}
}
////////////////////////////////////////////////////////////////////////////////////////
//calculate the center of sidechain and store it in AAs.coord
void Protein::setSCCenter(int AAIndx, char atomsIncluded)
{
	if ((AAIndx >= 0) && (AAIndx < numOfAA()))
	{
		if ((AAs[AAIndx].whtInCoord != 'S') || (AAs[AAIndx].atomsIncluded != atomsIncluded))
		{
			double segmaX  = 0,
					segmaY = 0,
					segmaZ = 0;
			int numSC = 0;  //counter for number of side chain atoms included...heavy atoms or all

			int i;
			if (atomsIncluded == 'H')
			{
				for (i=0;i<numOfAtoms(AAIndx);i++)
					if ((AAs[AAIndx].atoms[i].isSideChain) && (isHeavyAtom(AAs[AAIndx].atoms[i].type)))
					{
						segmaX += AAs[AAIndx].atoms[i].coord.x;
						segmaY += AAs[AAIndx].atoms[i].coord.y;
						segmaZ += AAs[AAIndx].atoms[i].coord.z;
						numSC++;
					}
			}
			if (atomsIncluded == 'A')
			{
				for (i=0;i<numOfAtoms(AAIndx);i++)
					if (AAs[AAIndx].atoms[i].isSideChain)
					{
						segmaX += AAs[AAIndx].atoms[i].coord.x;
						segmaY += AAs[AAIndx].atoms[i].coord.y;
						segmaZ += AAs[AAIndx].atoms[i].coord.z;
						numSC++;
					}
			}
			//if sidechain atoms are exist
			if (numSC)
			{
				AAs[AAIndx].coord.x = segmaX/numSC;
				AAs[AAIndx].coord.y = segmaY/numSC;
				AAs[AAIndx].coord.z = segmaZ/numSC;
			}
			//if no side chain atoms were found
			else
				//The center of SC is considered to be Ca
				AAs[AAIndx].coord = getAtomCoordinate(AAIndx," CA ");

			AAs[AAIndx].whtInCoord = 'S';		//set the flag to indicate that side chain mass center is stored in coord variable
			AAs[AAIndx].atomsIncluded = atomsIncluded;	//set the type of atoms were included in the calculations
		}
	}
	else
	{
		string eMsg = "The given index (";
		eMsg += toString(AAIndx);
		eMsg += ") is out of range..or the portion is empty (No.AAs= ";
		eMsg += toString(numOfAA());
		eMsg += ")";

		errMsg("Protein", "setSCCenter", eMsg);

	}
}
////////////////////////////////////////////////////////////////////////////////////////
//calculate the center of all atoms and store it in AAs.coord
void Protein::setAACenter(int AAIndx, char atomsIncluded)
{

	if ((AAIndx >= 0) && (AAIndx < numOfAA()))
	{
		if ((AAs[AAIndx].whtInCoord != 'A') || (AAs[AAIndx].atomsIncluded != atomsIncluded))
		{
			double segmaX  = 0,
					segmaY = 0,
					segmaZ = 0;

			int numOfAtomsIncluded = 0;
			int i;

			//heavy atoms in calculation
			if (atomsIncluded == 'H')
			{
				for (i=0;i<numOfAtoms(AAIndx);i++)
				{
					if (isHeavyAtom(AAs[AAIndx].atoms[i].type))
					{
						segmaX += AAs[AAIndx].atoms[i].coord.x;
						segmaY += AAs[AAIndx].atoms[i].coord.y;
						segmaZ += AAs[AAIndx].atoms[i].coord.z;
						numOfAtomsIncluded++;
					}
				}
			}
			//all atoms in calculation
			if (atomsIncluded == 'A')
			{
				for (i=0;i<numOfAtoms(AAIndx);i++)
				{
					segmaX += AAs[AAIndx].atoms[i].coord.x;
					segmaY += AAs[AAIndx].atoms[i].coord.y;
					segmaZ += AAs[AAIndx].atoms[i].coord.z;
					numOfAtomsIncluded++;
				}
			}

			//if AA has atoms
			if (numOfAtomsIncluded)
			{
				AAs[AAIndx].coord.x = segmaX/numOfAtomsIncluded;
				AAs[AAIndx].coord.y = segmaY/numOfAtomsIncluded;
				AAs[AAIndx].coord.z = segmaZ/numOfAtomsIncluded;
			}

			AAs[AAIndx].whtInCoord = 'A';	//set the falg to indicate that all-atoms mass center is stored in the coord variable
			AAs[AAIndx].atomsIncluded = atomsIncluded; //set the type of atoms included in the calculations
		}
	}
	else
	{
		string eMsg = "The given index (";
		eMsg += toString(AAIndx);
		eMsg += ") is out of range..or the portion is empty (No.AAs= ";
		eMsg += toString(numOfAA());
		eMsg += ")";

		errMsg("Protein", "setAACenter", eMsg);

	}

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
void Protein::setBBCenter (int AAIndx, char atomsIncluded)
{
	if ((AAIndx >= 0) && (AAIndx < numOfAA()))
	{
		if ((AAs[AAIndx].whtInCoord != 'B') || (AAs[AAIndx].atomsIncluded != atomsIncluded))
		{
			double segmaX  = 0,
					segmaY = 0,
					segmaZ = 0;

			int numOfAtomsIncluded = 0;
			int i;

			//heavy atoms in calculation
			if (atomsIncluded == 'H')
			{
				for (i=0;i<numOfAtoms(AAIndx);i++)
				{
					//heavy and not a side chain atom
					if ((isHeavyAtom(AAs[AAIndx].atoms[i].type)) && (!AAs[AAIndx].atoms [i].isSideChain))
					{
						segmaX += AAs[AAIndx].atoms[i].coord.x;
						segmaY += AAs[AAIndx].atoms[i].coord.y;
						segmaZ += AAs[AAIndx].atoms[i].coord.z;
						numOfAtomsIncluded++;
					}
				}
			}
			//all atoms in calculation
			if (atomsIncluded == 'A')
			{
				for (i=0;i<numOfAtoms(AAIndx);i++)
				{
					//not a side chain atom
					if (!AAs[AAIndx].atoms [i].isSideChain)
					{
						segmaX += AAs[AAIndx].atoms[i].coord.x;
						segmaY += AAs[AAIndx].atoms[i].coord.y;
						segmaZ += AAs[AAIndx].atoms[i].coord.z;
						numOfAtomsIncluded++;
					}
				}
			}

			//if AA has atoms
			if (numOfAtomsIncluded)
			{
				AAs[AAIndx].coord.x = segmaX/numOfAtomsIncluded;
				AAs[AAIndx].coord.y = segmaY/numOfAtomsIncluded;
				AAs[AAIndx].coord.z = segmaZ/numOfAtomsIncluded;
			}

			AAs[AAIndx].whtInCoord = 'B';				//set the flag to indicate that Backbone mass center is stored in the coord variable
			AAs[AAIndx].atomsIncluded = atomsIncluded;  //set the type of atoms included in the calculations
		}
	}
	else
	{
		string eMsg = "The given index (";
		eMsg += toString(AAIndx);
		eMsg += ") is out of range..or the portion is empty (No.AAs= ";
		eMsg += toString(numOfAA());
		eMsg += ")";

		errMsg("Protein", "setBBCenter", eMsg);

	}
}
//////////////////////////////////////////////////////////////////////////////////////
//atomsincluded H heavy A all
void Protein::setRgyration(int AAIndx, char atomsIncluded)
{
	if ((AAIndx >= 0) && (AAIndx < numOfAA()))
	{
		if (AAs[AAIndx].gyration == 0)
		{
			int numSC = numOfSCAtoms(AAIndx, atomsIncluded);
			if (numSC == 0)
			{
				AAs[AAIndx].gyration = sqrt(getRadius('C', 'C')); // gyration radius is the same as the radius of Ca	.... calculated radius is used
				AAs[AAIndx].gyCoord .x = AAs[AAIndx].gyCoord .y = AAs[AAIndx].gyCoord .z = AAs[AAIndx].gyration ;
			}
			else
			{
				int numAtoms = numOfAtoms(AAIndx);
				//calculate mass center if it is not calculated before
				if ((AAs[AAIndx].whtInCoord != 'S') || (AAs[AAIndx].atomsIncluded != atomsIncluded))
					setSCCenter(AAIndx);

				Coordinate massCenter = AAs[AAIndx].coord;
				//cout<<"mass center = "<<massCenter.x <<" "<<massCenter.y <<" "<<massCenter.z<<endl;

				int i;

				double gyrationR = 0;
				//initialize gyration coordinate
				AAs[AAIndx].gyCoord .x = 0;
				AAs[AAIndx].gyCoord .y = 0;
				AAs[AAIndx].gyCoord .z = 0;
				double atomR,
						rx,
						ry,
						rz,
						rr;
				if (atomsIncluded == 'A')
				{
					for (i=0; i<numAtoms; i++)
					{
						if (AAs[AAIndx].atoms[i].isSideChain)
						{
							atomR = getRadius(AAs[AAIndx].atoms [i].type , 'C');
							rx = AAs[AAIndx].atoms[i].coord .x - massCenter.x;
							ry = AAs[AAIndx].atoms[i].coord .y - massCenter.y;
							rz = AAs[AAIndx].atoms[i].coord .z - massCenter.z;

							rr = sqrt(rx*rx + ry*ry + rz*rz) + atomR;
							gyrationR += rr * rr;

							AAs[AAIndx].gyCoord .x += (rx + atomR) * (rx + atomR);
							AAs[AAIndx].gyCoord .y += (ry + atomR) * (ry + atomR);
							AAs[AAIndx].gyCoord .z += (rz + atomR) * (rz + atomR);
						}
					}
				}
				if (atomsIncluded == 'H')
				{
					for (i=0; i<numAtoms; i++)
					{
						if ((AAs[AAIndx].atoms[i].isSideChain) && (isHeavyAtom(AAs[AAIndx].atoms[i].type)))
						{
							atomR = getRadius(AAs[AAIndx].atoms [i].type , 'C');
							rx = AAs[AAIndx].atoms[i].coord .x - massCenter.x;
							ry = AAs[AAIndx].atoms[i].coord .y - massCenter.y;
							rz = AAs[AAIndx].atoms[i].coord .z - massCenter.z;

							rr = sqrt(rx*rx + ry*ry + rz*rz) + atomR;
							gyrationR += rr * rr;

							AAs[AAIndx].gyCoord .x += (rx + atomR) * (rx + atomR);
							AAs[AAIndx].gyCoord .y += (ry + atomR) * (ry + atomR);
							AAs[AAIndx].gyCoord .z += (rz + atomR) * (rz + atomR);


							//cout<<rx<<" "<<ry<<" "<<rz<<" "<<gyrationR<<endl;
						}
					}
				}

				AAs[AAIndx].gyCoord .x = sqrt(AAs[AAIndx].gyCoord .x / numSC);
				AAs[AAIndx].gyCoord .y = sqrt(AAs[AAIndx].gyCoord .y / numSC);
				AAs[AAIndx].gyCoord .z = sqrt(AAs[AAIndx].gyCoord .z / numSC);

				gyrationR = gyrationR/numSC;
				AAs[AAIndx].gyration = sqrt(gyrationR);
			}
		}
	}
	else
	{
		string eMsg = "The given index (";
		eMsg += toString(AAIndx);
		eMsg += ") is out of range..or the portion is empty (No.AAs= ";
		eMsg += toString(numOfAA());
		eMsg += ")";

		errMsg("Protein", "setRgyration", eMsg);

	}
}
/////////////////////////////////////////////////////////////////////////////////////
void Protein::setCentOfCharge (){

	float	totmass=0,
			totCharge = 0;

	int i, j;

	for (i=0; i<numOfAA(); i++){
		for (j=0; j<numOfAtoms(i); j++){
			float e=0;

			switch (AAs[i].atoms[j].type)
			{
				case 'H': e=1.0; totmass+=1.00794; break;
				case 'C': e=6.0; totmass+=12.0107; break;
				case 'A':       // treat 'A'mbiguous atoms as N, not perfect, but good enough
				case 'N': e=7.0; totmass+=14.00674; break;
				case 'O': e=8.0; totmass+=15.9994; break;
				case 'P': e=15.0; totmass+=30.973761; break;
				case 'S': e=16.0; totmass+=32.066; break;
				case 'W': e=18.0; totmass+=1.00794*2+15.9994; break;    // ficticious water atom
				/*
				case 'U':
						if (buf[12]=='A')
						{
							e=79.0; totmass+=196.96655;
						}
						break;  // gold
				*/
				default:
						printf("Unknown atom %c\n",AAs[i].atoms[j].type);
						e=0;
			}


			centOfCharge.x += AAs[i].atoms[j].coord.x * e;
			centOfCharge.y += AAs[i].atoms[j].coord.y * e;
			centOfCharge.z += AAs[i].atoms[j].coord.z * e;

			totCharge += e;
		}
	}
	centOfCharge.x /= totCharge;
	centOfCharge.y /= totCharge;
	centOfCharge.z /= totCharge;
}
///////////////////////////////////////////////////////////////////////////////////////
//the axis will be represented by set of points, the distance b/w the points will be given by "dist"
//except the last portion of the axis which could be less than this dist
void Protein::setAxisSegments (int startIndx, int endIndx, float dist, vector<Coordinate> &axis, vector<bool> &shortHelix, int num){

	string eMsg;
	//find the length of protein segment (in terms of # of AAs)
	short segLength = endIndx - startIndx + 1;
	axis.clear();		//clear the data structure

	int i;

	if (segLength>4){

		//set Backbone center for all AAs
		for (i=startIndx; i<=endIndx; i++)
			setBBCenter(i, 'H');		//work on heavy atoms only

		vector<Coordinate> triangles;
		vector<Coordinate> triangles2;
		vector<Coordinate> axis2;
		int stop = 0;
        if(segLength>8){
            //find center of triangles for each 3 consecutive AAs
            for (i=startIndx; i<=endIndx-3; i++) {
                //triangles.push_back (triangleCenter(AAs[i].coord, AAs[i+1].coord , AAs[i+2].coord));
                triangles.push_back (quadCenter(AAs[i].coord, AAs[i+1].coord , AAs[i+2].coord, AAs[i+3].coord));
                axis2.push_back(triangles[triangles.size ()-1]);
                stop = i;
            }
            for (i=0; i<stop-startIndx-2; i++) {
                //triangles.push_back (triangleCenter(AAs[i].coord, AAs[i+1].coord , AAs[i+2].coord));
                triangles2.push_back (quadCenter(axis2[i], axis2[i+1], axis2[i+2], axis2[i+3]));
                axis.push_back(triangles2[triangles2.size ()-1]);
            }
        }
        else{
            for (i=startIndx; i<=endIndx-3; i++) {
                //triangles.push_back (triangleCenter(AAs[i].coord, AAs[i+1].coord , AAs[i+2].coord));
                triangles.push_back (quadCenter(AAs[i].coord, AAs[i+1].coord , AAs[i+2].coord, AAs[i+3].coord));
                axis.push_back(triangles[triangles.size ()-1]);
            }
        }
/**

		//find N terminal...if missing Ca ..error otherwise
		Coordinate terminalP,			//terminal point
					newP;				//the point where the line intersects the point
		short indx = getAtomIndx(startIndx, " N  ");
		if (indx != -1)
			terminalP = AAs[startIndx].atoms [indx].coord;
		else{
			indx = getAtomIndx(startIndx, " CA ");
			if (indx != -1)
				terminalP = AAs[startIndx].atoms [indx].coord;
			else{

				eMsg  = "First AA (";
				eMsg += AAs[startIndx].chr3 ;
				eMsg += "-";
				eMsg += toString(AAs[startIndx].num);
				eMsg += ") has no N or Ca atoms..";

				errMsg("Protein", "setAxisSegments", eMsg);

			}
		}
		//find the projection of this point with the first two triangles
		linePointIntersection(triangles[0], triangles[1], terminalP, newP);
		triangles.insert (triangles.begin (), newP);

		//find C terminal..if missing Ca....error otherwise
		indx = getAtomIndx(endIndx, " C  ");
		if (indx != -1)
			terminalP = AAs[endIndx].atoms [indx].coord;
		else{
			indx = getAtomIndx(endIndx, " CA ");
			if (indx != -1)
				terminalP = AAs[endIndx].atoms [indx].coord;
			else{

				eMsg  = "Last AA ("+AAs[endIndx].chr3 ;
				eMsg += "-";
				eMsg += toString(AAs[endIndx].num) ;
				eMsg += ") has no C or Ca atoms..";

				errMsg("Protein", "setAxisSegments", eMsg);
			}
		}
		//find the projection of this point with last two triangles
		linePointIntersection(triangles[triangles.size ()-2], triangles[triangles.size ()-1], terminalP, newP);
		triangles.push_back (newP);



		short 	sIndx	= 0,		//start indx
				eIndx	= triangles.size ()-1;		//end indx


		Coordinate curP = triangles[sIndx];	//current point we are working on...first point in the set of triangles

		axis.push_back(curP);


		int lastPointIndx = -1;
		i= sIndx + 1;
		while (i <= eIndx){
			if (getDistance(curP, triangles[i]) > dist){
				axis.push_back (pointOnLine(curP, triangles[i], dist));
				curP = axis[axis.size()-1];
				lastPointIndx = i;
			}
			else
				i +=1;
		}

		// add the last point if needed
		if (lastPointIndx <= eIndx)
			if (getDistance(curP, triangles[triangles.size ()-1]) > 0.0)
				axis.push_back (triangles[triangles.size ()-1]);


		/*
		cout<<"# of points in triangles = "<<triangles.size ()<<endl;
		float totalDist=0;
		for (i=0; i<triangles.size()-1; i++){
			cout<<i+1<<" distance = "<<getDistance(triangles[i], triangles[i+1])<<endl;
			totalDist += getDistance(triangles[i], triangles[i+1]);
		}

		cout<<totalDist<<endl;
		cout<<"=============="<<endl;
		for (i=0; i<axis.size()-1; i++)
			cout<<i+1<<"  distance = "<<getDistance(axis[i], axis[i+1])<<endl;
		*/

/**

	}
	else{

		errMsg("Protein", "setAxisSegments", "Be sure the segment sent is longer than 3 AAs..no segments ");
	*/ }
else
{
    shortHelix[num] = true;
   setShortAxis(startIndx, endIndx, dist, axis);
}


}

void Protein::setAxisSegments2 (int startIndx, int endIndx, float dist, vector<Coordinate> &axis, Coordinate last, vector<bool> &shortHelix, int num){

	string eMsg;
	//find the length of protein segment (in terms of # of AAs)
	short segLength = endIndx - startIndx + 1;
	axis.clear();		//clear the data structure

	int i;
	Coordinate p;

	if (true){

		//set Backbone center for all AAs
		for (i=startIndx-1; i<=endIndx+1; i++)
			setBBCenter(i, 'H');		//work on heavy atoms only

		vector<Coordinate> triangles;

		//find center of triangles for each 3 consecutive AAs
		for (i=startIndx-1; i<=endIndx-1; i++)
        {
			//triangles.push_back (triangleCenter(AAs[i].coord, AAs[i+1].coord , AAs[i+2].coord));
			triangles.push_back (quadCenter2(AAs[i].coord, AAs[i+1].coord));
			axis.push_back(triangles[triangles.size ()-1]);
			//triangles.push_back (quadCenter1(AAs[i].coord));
			//axis.push_back(triangles[triangles.size ()-1]);
        }
        triangles.push_back (last);
        axis.push_back(triangles[triangles.size ()-1]);

}
else
{
    shortHelix[num] = true;
    setShortAxis(startIndx, endIndx, dist, axis);
}


}

void Protein::setShortAxis(int startIndx, int endIndx, float dist, vector<Coordinate> &axis)
{
	short segLength = endIndx - startIndx + 1;

	if (segLength>=3){
	    Coordinate c1, c2, c3, c4;

        c1 = AAs[startIndx].atoms[getAtomIndx(startIndx, "N")].coord;
        c2 = AAs[startIndx+1].atoms[getAtomIndx(startIndx+1, "N")].coord;
        c3 = AAs[startIndx+2].atoms[getAtomIndx(startIndx+2, "N")].coord;
    if(segLength==4)
        {
            c4 = AAs[endIndx].atoms[getAtomIndx(endIndx, "N")].coord;
            axis.push_back(quadCenter(c1, c2, c3, c4));
        }
        else axis.push_back(triangleCenter(c1, c2, c3));
/*
        c1 = AAs[startIndx].atoms[getAtomIndx(startIndx, "CA")].coord;
        c2 = AAs[startIndx+1].atoms[getAtomIndx(startIndx+1, "CA")].coord;
        c3 = AAs[startIndx+2].atoms[getAtomIndx(startIndx+2, "CA")].coord;
    if(segLength==4)
        {
            c4 = AAs[endIndx].atoms[getAtomIndx(endIndx, "N")].coord;
            axis.push_back(quadCenter(c1, c2, c3, c4));
        }
        else axis.push_back(triangleCenter(c1, c2, c3));

        c1 = AAs[startIndx].atoms[getAtomIndx(startIndx, "C ")].coord;
        c2 = AAs[startIndx+1].atoms[getAtomIndx(startIndx+1, "C ")].coord;
        c3 = AAs[startIndx+2].atoms[getAtomIndx(startIndx+2, "C ")].coord;
    if(segLength==4)
        {
            c4 = AAs[endIndx].atoms[getAtomIndx(endIndx, "N")].coord;
            axis.push_back(quadCenter(c1, c2, c3, c4));
        }
        else axis.push_back(triangleCenter(c1, c2, c3));
    */
    }
}
/////////////////////////////////////////////////////////////////////////////////////
void Protein::buildSS()
{

	if (numOfAA())
	{
		//string line;
		int i, j;

		//clear both vectors
		hlces.clear();
		sheets.clear();

		for (i=0; i<header.size(); i++)
		{

			HelicesSecondaryStructure tmpHlx;
			SheetsSecondaryStructure  tmpStrand;

			// if encounter HELIX
			if ((header[i].substr(0,5) == "HELIX") && (header[i].substr(19,1).c_str() == AAs[0].chain))
			{
				tmpHlx.serialNum		= atoi(header[i].substr(7, 3).c_str ());
				tmpHlx.hlxID			= atoi(header[i].substr(11, 3).c_str ());
				tmpHlx.startIndx		= getAAIndx(atoi(header[i].substr(21,4).c_str()));
				tmpHlx.endIndx			= getAAIndx(atoi(header[i].substr(33,4).c_str()));
				tmpHlx.type				= atoi(header[i].substr(38,2).c_str());
				tmpHlx.comment			= header[i].substr(40, 30);

				if ((tmpHlx.startIndx != -1) &&									//does the AA exist
					(tmpHlx.endIndx != -1) &&									//does the 2nd AA exist
					(header[i].substr(31,1) == header[i].substr(19,1)) )	//check the other chain ... normally the hlx should be in the same chain
				{

					tmpHlx.nAA				= abs(AAs[tmpHlx.endIndx ].num - AAs[tmpHlx.startIndx].num) + 1;
					tmpHlx.nAAcur		= tmpHlx.endIndx - tmpHlx.startIndx + 1;
					hlces.push_back(tmpHlx);
				}

			}

			//build sheets vector
			if ((header[i].substr(0,5) == "SHEET") && (header[i].substr(21,1).c_str() == AAs[0].chain))
			{

				tmpStrand.strandNum = atoi(header[i].substr (7, 3).c_str ());
				tmpStrand.sheetID	= header[i].substr (11, 3).c_str ();
				tmpStrand.nStrand	= atoi(header[i].substr (14, 2).c_str ());
				tmpStrand.startIndx = getAAIndx(atoi(header[i].substr(22,4).c_str()));
				tmpStrand.endIndx	= getAAIndx(atoi(header[i].substr(33,4).c_str()));

				tmpStrand.sense		= atoi(header[i].substr(38,2).c_str());

				if (header[i].length () > 44){

					tmpStrand.curAtomName = header[i].substr (41, 4);
					tmpStrand.AACurIndx		= getAAIndx(atoi(header[i].substr(50, 4).c_str()));
					tmpStrand.prevAtomName = header[i].substr (56, 4);
					tmpStrand.AAPrevIndx	= getAAIndx(atoi(header[i].substr(65, 4).c_str()));
				}

				if ((tmpStrand.startIndx != -1) &&
					(tmpStrand.endIndx != -1))
				{
					tmpStrand.nAAcur		=  tmpStrand.endIndx - tmpStrand.startIndx + 1;
					tmpStrand.nAA			= abs(AAs[tmpStrand.endIndx ].num - AAs[tmpStrand.startIndx].num) + 1;
					sheets.push_back(tmpStrand);
				}
			}
		}

		//sort hlces and strands according to the index of the first AA assigned to them
		if (hlces.size () > 0)
		{
			sortHlces(hlces, 0, hlces.size()-1);
			//remove duplicates
			i=0;
			while (i <hlces.size())
			{
				j = i+1;
				while (j< hlces.size())
				{
					if (hlces[i].startIndx == hlces[j].startIndx ){
						hlces.erase (hlces.begin () + j);		//delete and then decrement counter
						j--;
					}
					j++;
				}
				i++;
			}

		}
		if (sheets.size ())
		{
			sortStrands(sheets, 0, sheets.size ()-1);

			//remove duplicates
			i=0;

			while (i < sheets.size())
			{
				j = i+1;
				while (j < sheets.size())
				{
					if (sheets[i].startIndx  == sheets[j].startIndx){
						sheets.erase (sheets.begin () + j);		//delete and then decrement counter
						j--;
					}

					j++;
				}
				i++;
			}

		}
	}
	else
	{
		string eMsg = "Be sure the portion (" + path;
		eMsg += ") contains AAs or you have read the pdb file..";

		errMsg("Protein", "buildSS", eMsg);

	}


}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//transform the portion by a given x, y, and z
void Protein::translateBy(Coordinate transAmount)
{
	for (int i=0; i<numOfAA();i++){
		for (int j=0;j<numOfAtoms(i);j++){
			AAs[i] .atoms[j] .coord .x += transAmount.x;
			AAs[i] .atoms[j] .coord .y += transAmount.y;
			AAs[i] .atoms[j] .coord .z += transAmount.z;
		}
		//initiate whtInCoord variable
		AAs[i].whtInCoord = 'N';
		//initiate the gyration
		AAs[i].gyration = 0;

		/*****
				the coordinate of Rgyration is incorrect right now
		*****/
		AAs[i].gyCoord .x = 0;
		AAs[i].gyCoord .y = 0;
		AAs[i].gyCoord .z = 0;
	}
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
void Protein::translateBy(int numOfCaAtoms)
{
	int tmpNumOfAA = numOfAA();
	if (tmpNumOfAA >= 2)
	{
		if (numOfCaAtoms)
		{
			Coordinate startEndCoord, endEndCoord;
			int AAsCounter = 0;	//for atoms will be involved in the center axial
			//take at most 4 atoms from each end
			while ((AAsCounter < tmpNumOfAA - 1) && (AAsCounter < 4))
			{
				Coordinate tmpCoord;

				//get coordinate of Ca atoms from the start end
				tmpCoord = getAtomCoordinate(AAsCounter," CA ");
				startEndCoord.x += tmpCoord.x;
				startEndCoord.y += tmpCoord.y;
				startEndCoord.z += tmpCoord.z;

				//get coordinate of Ca atoms from the end end
				tmpCoord = getAtomCoordinate(tmpNumOfAA - AAsCounter - 1," CA ");
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

			//find the 2 correspponding points represent the imiginary Ca atoms @ the line ... from any end..here from the end (could be from start end)
			Coordinate imiginaryCaLast = pointLineIntersection(getAtomCoordinate(tmpNumOfAA-1," CA "),startEndCoord,endEndCoord);
			Coordinate imiginaryCaBeforeLast = pointLineIntersection(getAtomCoordinate(tmpNumOfAA-2," CA "),startEndCoord,endEndCoord);


			//move right or left
			Vectors tmpV(imiginaryCaLast,imiginaryCaBeforeLast);
			tmpV = tmpV.mul(numOfCaAtoms);
			translateBy(tmpV.getCoordinates());
		}
	}
	else
	{
		string eMsg = "The size of the portion (";
		eMsg += toString(tmpNumOfAA);
		eMsg += ") should be greater or equal 2 AA's...";

		errMsg("Protein", "translateBy", eMsg);
	}
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//overlap the portion with the second line (represented by p3 and p4).....
//atom1Coord and atom2Coord represent the line (bond) in the portion to be overlaped with line 2
//the coordinate of points will be changed to become equal to p3 and p4 respectively
// so atom1 coordinate becomes as p3 and atom2 coordinate becomes as p4
void Protein::overlapLine(Coordinate atom1Coord, Coordinate atom2Coord, Coordinate p3, Coordinate p4)
{

	//transform the portion so atom1Coord overlaps p3
	Vectors vectorL2(p3);

	vectorL2 -= atom1Coord;
	translateBy(vectorL2.getCoordinates());

	//transform atom2Coord also...no need to translate atom1Coord b/s it is the same value with p3 right now
	atom2Coord.x += vectorL2.getX();
	atom2Coord.y += vectorL2.getY();
	atom2Coord.z += vectorL2.getZ();

	//find the angle between points...it is better than calling getAngleRadian b/s I need v1 and v2 here...
	Vectors v1(atom2Coord,p3);
	Vectors v2(p4,p3);

	double angle = v1.getAngleRadian(v2);

	//if there is a need to rotate
	if (angle)
	{
		Vectors vNormal;
		vNormal = v2.cross(v1);
		vNormal += p3;

		double mtx[4][4];
		//build rotation matrix
		buildRotationMatrix(mtx,p3,vNormal.getCoordinates(),-angle);

		//rotate all atoms
		for (int i=0;i<numOfAA();i++)
		{
			for (int j=0;j<numOfAtoms(i);j++)
				rotatePoint(mtx,AAs[i].atoms[j].coord );

			//initiate the whtInCoord variable
			AAs[i].whtInCoord = 'N';
			//initiate the gyration variable
			AAs[i].gyration = 0;
		}

	}


}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
// connect a portion with the portion from the begining or the end....
//if connecting from the end....portionToBeConnected will be connected from the second AA
//if connecting from the beginning...the last AA of portionToBeConnected will be deleted
// by overlapping the three atoms (N,CA, and C) of the last residue of the portionToBeConnected with the corresponding atoms of the first residue of the portion (S)
//or overlapping the three atoms (N, CA, and C) of the first residue of the portionToBeConnected with the corresponding atoms of the last residue of the portion (E)
void Protein::connect(Protein portionToBeConnected,char method)
{
	int portionNumOfAA = portionToBeConnected.numOfAA();
	int tmpNumOfAA = numOfAA();
	if ((tmpNumOfAA >= 1) && (portionNumOfAA >= 2))
	{

		if (method == 'E')	//connect portionToBeConnected with the portion from the end
		{
			//Move the portionToBeConnected to the Appropriate position
			getCPosition(portionToBeConnected,'E');

			//change O atom of the last AA to avoid collision
			//O atom is considered the last atom in the atom vector
			//int OIndx = getAtomIndx(tmpNumOfAA-1," O  ");							//O atom of the last AA
			//int portionOIndx = portionToBeConnected.getAtomIndx (0," O  ");			//O atom of the first AA
			AAs[tmpNumOfAA-1].atoms[getAtomIndx(tmpNumOfAA-1," O  ")].coord = portionToBeConnected.AAs[0].atoms[portionToBeConnected.getAtomIndx (0," O  ")].coord ;

			//if there is an OXT atom in the terminal Carbon...delete it
			int oxtIndx = getAtomIndx(tmpNumOfAA-1," OXT");
			if ( oxtIndx != -1)
				AAs[tmpNumOfAA-1].atoms .erase (AAs[tmpNumOfAA-1].atoms .begin () + oxtIndx );

			//reset the gyration for this AA
			AAs[tmpNumOfAA - 1].gyration = 0;
			//reset the wht in cooordinate
			AAs[tmpNumOfAA-1].whtInCoord = 'N';

			//Add AAs from portionToBeConnected to the portion starting from the second AA
			append(portionToBeConnected,1,portionNumOfAA-1);


		}
		else
			if (method == 'S')	//connect portionToBeConnected with the portion from the beginning
			{
				//Move the portionToBeConnected to the appropriate position before the beginning of the portion
				getCPosition(portionToBeConnected,'S');

				//Add AA's from portionToBeConnected...the last AA will be ignored
				int seqNum = AAs[0].num;		//get the seq number that portion starts from
				string chain = AAs[0].chain;	//save the chain name

				//connect the portionToBeConnected from the start
				for (int i = portionNumOfAA - 2;i>=0;i--)
				{
					AAs.insert(AAs.begin(),portionToBeConnected.AAs[i]);		//insert reversely to the portion
					AAs[0].chain = chain;										//change the chain name
					AAs[0].num = --seqNum;										//change sequence number
				}
				//in the case that seq number is below 1..reset seq numbers
				if (AAs[0].num < 1)
					for (int i=0;i<numOfAA();i++)	//not tmpNumOfAA b/s it has been changed
						AAs[i].num = i+1;

/*				Another way to add AA's from portionToBeConnected

				int seqNum = AAs[0].num - portionToBeConnected.numOfAA();		//get the seq number that portion starts from
				string chain = AAs[0].chain;	//save the chain name
				portionToBeConnected.append(*this,0,numOfAA()-1);
				for (int i=0;i<portionToBeConnected.numOfAA();i++)
				{
					portionToBeConnected.AAs[i].chain = chain;
					portionToBeConnected.AAs[i].num = seqNum++;
				}
				*this = portionToBeConnected;
*/
			}
			else
			{
				string eMsg = "The method sent (" + method;
				eMsg += ") is incorrect..plz choose S (from the beginning) or E (from the end)";

				errMsg("Protein", "connect", eMsg);

			}
	}
	else
	{
		string eMsg = "The length of connected portion should be 2 AAs or more and the number of AAs in the portion (the one u r connecting with) should be 1 AA at least.";

		errMsg("Protein", "connect", eMsg);
	}
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
// concatening a portionToBeConnected with the current portion from the beginning or the end
// the length of the current portion should be at least 1 and for the portion to be connected 2
// phi and psi given as the phi and the psi for the bonds connecting both portions
// the portion will not lose any AA like wht happend in connect(Protein, method)
void Protein::concat(Protein portionToBeConnected, char method, double desiredPhi, double desiredPsi)
{
	int portionNumOfAA = portionToBeConnected.numOfAA();
	int oldNumOfAA = numOfAA();
	if ((oldNumOfAA >= 1) && (portionNumOfAA >= 2))
	{

		if (method == 'E')	//connect portionToBeConnected with the portion from the end
		{

			//Move the portionToBeConnected to the Appropriate position
			getCPosition(portionToBeConnected,'E');

			//get the new coordinate for the last O atom...O atom is the last atom in the atoms list
			//int OIndx = getAtomIndx(oldNumOfAA-1," O  ");
			//int portionOIndx = portionToBeConnected.getAtomIndx(0," O  ");
			AAs[oldNumOfAA-1].atoms[getAtomIndx(oldNumOfAA-1," O  ")].coord = portionToBeConnected.AAs[0].atoms[portionToBeConnected.getAtomIndx(0," O  ")].coord;

			//if there is an OXT atom in the terminal Carbon...delete it
			int oxtIndx = getAtomIndx(oldNumOfAA-1," OXT");
			if ( oxtIndx != -1)
				AAs[oldNumOfAA-1].atoms .erase (AAs[oldNumOfAA-1].atoms .begin () + oxtIndx );

			//reset the gyration for this AA
			AAs[oldNumOfAA-1].gyration = 0;
			//reset the wht in cooordinate
			AAs[oldNumOfAA-1].whtInCoord = 'N';

			//translate the portionToBeConnected to N position of the second AA
			Coordinate portionNAtomCoord = portionToBeConnected.getAtomCoordinate(1," N  ");
			Coordinate deltaCoord = portionToBeConnected.getAtomCoordinate(0, " N  ");

			deltaCoord.x = portionNAtomCoord.x - deltaCoord.x;
			deltaCoord.y = portionNAtomCoord.y - deltaCoord.y;
			deltaCoord.z = portionNAtomCoord.z - deltaCoord.z;

			//portionToBeConnected.writePDB("portionBFTranslation.pdb",1,portionNumOfAA);

			//translate the portion
			portionToBeConnected.translateBy(deltaCoord);

			//portionToBeConnected.writePDB("portionafterTranslation.pdb",1,portionNumOfAA);

			//Add AAs from portionToBeConnected to the portion starting from the second AA
			append(portionToBeConnected,0,portionNumOfAA-1);

			//make Omega equal to 178
			Coordinate CaAtomCoord    = getAtomCoordinate(oldNumOfAA-1," CA "),
				       CAtomCoord     = getAtomCoordinate(oldNumOfAA-1," C  "),
			           newCaAtomCoord = getAtomCoordinate(oldNumOfAA," CA ");		//for the first AA appended


			//make Omega always equal to 178 degree
			double deltaOmega = getTorsionAngle(CaAtomCoord, CAtomCoord, portionNAtomCoord, newCaAtomCoord) - 178;

			int curNumOfAA = numOfAA();

			//rotate atoms around Omega
			rotate(oldNumOfAA,				//first AA appended
					curNumOfAA-1,			//last AA
					1,						//from the second atom...the one after N
					CAtomCoord,				//
					portionNAtomCoord,		//Around C-Ni+1
					deltaOmega);			//angle in degree


			//writePDB("portionafterRotation.pdb",1,curNumOfAA-1);

			//find the amount of degrees psi should be rotated
			double deltaPsi = getTorsion("psi",oldNumOfAA-1) - desiredPsi;	//find the psi of the last (before appending) AA of current portion...


			//rotate atoms around Psi
			rotate(oldNumOfAA-1,								//start AA
					curNumOfAA-1,								//end AA
					getAtomIndx(oldNumOfAA-1," C  "),			//start rotation from C directly
					CaAtomCoord,								//
					CAtomCoord,									//around bond CA-C
					deltaPsi);									//angle in degree

			//find the amound of degrees phi should be rotated
			double deltaPhi = getTorsion("phi",oldNumOfAA) - desiredPhi;		//one AA after

			//rotate atoms aroud Phi
			rotate(oldNumOfAA,									//start AA
					curNumOfAA-1,								//end AA
					getAtomIndx(oldNumOfAA," CA "),				//start atom within start AA
					getAtomCoordinate(oldNumOfAA," N  "),		//	has been rotated
					getAtomCoordinate(oldNumOfAA," CA "),		// (has been rotated)  ... .Around N-CA
					deltaPhi);									//angle in degree


		}
		else
			if (method == 'S')	//connect portionToBeConnected with the portion from the beginning
			{
				//Move the portionToBeConnected to the appropriate position before the beginning of the portion
				getCPosition(portionToBeConnected,'S');

				//translate the portionToBeConnected to new C position of AA before that last
				Coordinate deltaCoord = portionToBeConnected.getAtomCoordinate(portionNumOfAA - 1, " C  ");			//last C
				Coordinate portionCAtomCoord = portionToBeConnected.getAtomCoordinate(portionNumOfAA - 2," C  ");	//C atom of the AA before the last

				deltaCoord.x = portionCAtomCoord.x - deltaCoord.x;
				deltaCoord.y = portionCAtomCoord.y - deltaCoord.y;
				deltaCoord.z = portionCAtomCoord.z - deltaCoord.z;

				//portionToBeConnected.writePDB("portionBFTranslation.pdb",1,portionNumOfAA);


				portionToBeConnected.translateBy(deltaCoord);

				//portionToBeConnected.writePDB("portionafterTranslation.pdb",1,portionNumOfAA);


				//Add AA's from portionToBeConnected...
				int seqNum = AAs[0].num;		//get the seq number that portion starts from
				string chain = AAs[0].chain;	//save the chain name

				//connect the portionToBeConnected from the start
				for (int i = portionNumOfAA-1;i>=0;i--)
				{
					AAs.insert(AAs.begin(),portionToBeConnected.AAs[i]);		//insert reversely to the portion
					AAs[0].chain = chain;										//change the chain name
					AAs[0].num = --seqNum;										//change sequence number
				}
				//in the case that seq number is below 1..reset seq numbers
				if (AAs[0].num < 1)
					for (int i=0;i<oldNumOfAA;i++)
						AAs[i].num = i+1;

				//make Omega equal to 178
				Coordinate NAtomCoord      = getAtomCoordinate(portionNumOfAA," N  "),
						   CaAtomCoord     = getAtomCoordinate(portionNumOfAA," CA "),
						   newCaAtomCoord  = getAtomCoordinate(portionNumOfAA-1," CA ");


				//make Omega always equal to 178 degree
				double deltaOmega = getTorsionAngle(newCaAtomCoord, portionCAtomCoord, NAtomCoord, CaAtomCoord) - 178;

				int curNumOfAA = numOfAA();

				//rotate atoms around Omega
				rotate(portionNumOfAA,				//first AA appended
						curNumOfAA-1,			//last AA
						1,						//from the second atom...the one after N
						portionCAtomCoord,		//
						NAtomCoord,				//Around C-Ni+1
						deltaOmega);			//angle in degree


				//find the amount of degrees psi should be rotated
				double deltaPsi = getTorsionAngle(getAtomCoordinate(portionNumOfAA-1," N  "),
												newCaAtomCoord,
												portionCAtomCoord,
												NAtomCoord) - desiredPsi;

				rotate(portionNumOfAA-1,
					   curNumOfAA-1,
					   getAtomIndx(portionNumOfAA-1," C  "),
					   newCaAtomCoord,
					   portionCAtomCoord,
					   deltaPsi);


				//find the amount of degrees phi should be rotated
				NAtomCoord = getAtomCoordinate(portionNumOfAA," N  ");
				CaAtomCoord = getAtomCoordinate(portionNumOfAA," CA ");
				Coordinate CAtomForPhi = getAtomCoordinate(portionNumOfAA," C  ");

				double deltaPhi = getTorsionAngle ( portionCAtomCoord,
													NAtomCoord,
													CaAtomCoord,
													CAtomForPhi) - desiredPhi;

				rotate(portionNumOfAA,
						curNumOfAA-1,
						getAtomIndx(portionNumOfAA," CA "),
						NAtomCoord,
						CaAtomCoord,
						deltaPhi);

				//writePDB("portion.pdb",1,curNumOfAA);
			}
			else
			{
				string eMsg = "The method sent (" + method;
				eMsg += ") is incorrect..plz choose S (from the beginning) or E (from the end)";

				errMsg("Protein", "concat", eMsg);

			}
	}
	else
	{
		string eMsg = "The length of connected portion should be 2 AAs or more and the number of AAs in the portion (the one u r connecting with) should be 1 AA at least.";

		errMsg("Protein", "concat", eMsg);

	}
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//AAWithSideChain is the AA to replace the old AA side chain specefied by AAIndx
void Protein::plugSideChain(int AAIndx, Protein AAWithSideChain)
{
	if ((AAIndx >= 0) && (AAIndx < numOfAA()))
	{
		int sideChainNumOfAA = AAWithSideChain.numOfAA();
		if ((sideChainNumOfAA == 1) && (AAWithSideChain.numOfAtoms(0)))
		{


			Coordinate NCoord   = getAtomCoordinate(AAIndx," N "),
						CaCoord = getAtomCoordinate(AAIndx," CA "),
						CCoord  = getAtomCoordinate(AAIndx," C "),
						sideChainNCoord  = AAWithSideChain.getAtomCoordinate(0," N "),
						sideChainCaCoord = AAWithSideChain.getAtomCoordinate (0," CA "),
						sideChainCCoord  = AAWithSideChain.getAtomCoordinate (0," C ");

			//overlap CA-C bond
			AAWithSideChain.overlapLine(sideChainCaCoord, sideChainCCoord, CaCoord, CCoord);

			//overlap N-CA bond
			Coordinate tmpCoord;

			sideChainNCoord = AAWithSideChain.getAtomCoordinate(0," N "); //new N atom after overlapping

			//originalTorison - currentTorsion is the amount of degree to be rotated
			double originalTorsion = getTorsionAngle(NCoord, CaCoord, CCoord, tmpCoord),
					currentTorsion = getTorsionAngle(sideChainNCoord, CaCoord, CCoord, tmpCoord);


			//rotate side chain
			//here the start atomIndx is not important b/s I will discard the backbone and take just the sidechain..b/s that I've started from atom 0
			AAWithSideChain.rotate(0, sideChainNumOfAA - 1, 0, CaCoord, CCoord, originalTorsion - currentTorsion);

			//erase old side chain atoms
			int numOfAtomsErased = 0;
			int i = 0;
			while ((i - numOfAtomsErased) < numOfAtoms(AAIndx))
			{
				if (AAs[AAIndx].atoms[i - numOfAtomsErased].isSideChain)
				//if ((AAs[AAIndx].atoms[i - numOfAtomsErased].name != " N  ") || (AAs[AAIndx].atoms[i - numOfAtomsErased].name != " CA ") || (AAs[AAIndx].atoms[i - numOfAtomsErased].name != " C  "))
					AAs[AAIndx].atoms.erase(AAs[AAIndx].atoms.begin () + i - numOfAtomsErased++);
				i++;
			}

			//take all atoms after C atoms temporarly
			int CAtomIndx = getAtomIndx(AAIndx, " C  ");
			vector<Atom> afterCAtom (AAs[AAIndx].atoms.begin () + CAtomIndx, AAs[AAIndx].atoms.end ());

			//clear all atoms after C atom....C atom is included
			AAs[AAIndx].atoms.erase (AAs[AAIndx].atoms.begin () + CAtomIndx, AAs[AAIndx].atoms .end ());

			//push back the new side chain atoms
			for (i=0;i<AAWithSideChain.numOfAtoms(0);i++)
				if (AAWithSideChain.AAs[0].atoms[i].isSideChain)
					AAs[AAIndx].atoms.push_back(AAWithSideChain.AAs[0].atoms[i]);

			//re-push atom C and after
			for (i=0;i<afterCAtom.size();i++)
				AAs[AAIndx].atoms.push_back(afterCAtom[i]);

			//reset wht in coordinate flag
			AAs[AAIndx].whtInCoord = 'N';
			//reset the radius of gyration
			AAs[AAIndx].gyration = 0;		//if the gyration is 0...the gyCoord in invalid as well

		}
		else
		{
			errMsg("Protein", "plugSideChain", "AAWithSideChain given is empty, or greater than 1 AA...");
		}
	}
	else
	{
		string eMsg = "The given index (";
		eMsg += toString(AAIndx);
		eMsg += ") is out of range..or the portion is empty (No.AAs= ";
		eMsg += toString(numOfAA());
		eMsg += ")";

		errMsg("Protein", "setRgyration", eMsg);
	}
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
void Protein::setScEndPoint(int AAIndx, char atomsIncluded)
{
	if ((AAIndx >= 0) && (AAIndx < numOfAA()))
	{
		if ((AAs[AAIndx].whtInCoord != 'S') || (AAs[AAIndx].atomsIncluded != atomsIncluded))
		{
			int CaIndx = getAtomIndx(AAIndx, "CA");

			if (CaIndx != -1)
			{
				setSCCenter(AAIndx, atomsIncluded);		//set the center of side chain according to the atom included heavy or all
				setRgyration(AAIndx, atomsIncluded);	//set the radius of gyration for the side chain

				Coordinate CaCoordinate = AAs[AAIndx].atoms [CaIndx].coord;

				Vectors v(AAs[AAIndx].coord, CaCoordinate);				//the vector from Ca atom to the massCenter
				double factorNum = v.length() / AAs[AAIndx].gyration;	// the number is used to divide the vector v by to get the exact vector to be added to the vector v
																		// to add the portion of the vector after the massCenter toward the end of the sidechain. this
																		//addition will be equal to the Rg
				Vectors vToBeAdded = v;
				vToBeAdded.divide (factorNum);							//get the vector of length equal to Rg to be added to the current vector (from CA to massCenter)
				v += vToBeAdded;

				AAs[AAIndx].ScEndPoint.x = v.getX() + CaCoordinate.x;
				AAs[AAIndx].ScEndPoint.y = v.getY() + CaCoordinate.y;
				AAs[AAIndx].ScEndPoint.z = v.getZ() + CaCoordinate.z;
			}
			else
			{
				//string eMsg = "Ca atom is not found (" + AAs[AAIndx].chr3;
				//eMsg += "-";
				//eMsg += toString(AAs[AAIndx].num);
				//eMsg += ") ....no EndPoint has been computed";

				//errMsg("Protein", "setScEndPoint", eMsg);
			}


		}
	}
	else
	{
		string eMsg = "The given index (";
		eMsg += toString(AAIndx);
		eMsg += ") is out of range..or the portion is empty (No.AAs= ";
		eMsg += toString(numOfAA());
		eMsg += ")";

		errMsg("Protein", "setScEndPoint", eMsg);
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////
void Protein::torsion2coord()
{
	int i;
	bool convert = true;
	vector <Torsion>   torsions;		//list of torsion angles of current protein

	for (i=0; i<numOfAA(); i++)
	{
		if ((AAs[i].angles.phi == 999.00) || (AAs[i].angles.psi == 999.00))
		{
			string eMsg = "Torsion angles of the protein could not be converted to xyz if they are not set for AA (" + AAs[i].chr3;
			eMsg += toString(AAs[i].num);
			eMsg += ") they are not set...";

			errMsg("Protein", "torsion2coord", eMsg);

			convert = false;

			break;
		}

		torsions.push_back (AAs[i].angles);
	}
	if (convert)
	{

		vector<Coordinate> XYZlist;			//where the coordinate values will be stored
		torsion2xyz(torsions, XYZlist);

		AAs.clear();		//delete all previous atoms
		Atom tmpAtom;

		int seqNum = 1;

		for (i=0; i<XYZlist.size (); i+=4)
		{
			AminoAcid tmpAA;
			tmpAA.chr3 = "GLY";
			tmpAA.chr1 = 'G';
			tmpAA.chain = "A";

			tmpAtom.name = " N  ";
			tmpAtom.coord = XYZlist[i];
			tmpAtom.type = 'N';
			tmpAA.atoms .push_back (tmpAtom);

			tmpAtom.name = " CA ";
			tmpAtom.coord = XYZlist[i+1];
			tmpAtom.type = 'C';
			tmpAA.atoms .push_back (tmpAtom);

			tmpAtom.name = " C  ";
			tmpAtom.coord = XYZlist[i+2];
			tmpAtom.type = 'C';
			tmpAA.atoms .push_back (tmpAtom);

			tmpAtom.name = " O  ";
			tmpAtom.coord = XYZlist[i+3];
			tmpAtom.type = 'O';
			tmpAA.atoms .push_back (tmpAtom);

			tmpAA.num = seqNum++;
			AAs.push_back(tmpAA);
			AAs[seqNum-2].whtInCoord = 'N';
			AAs[seqNum-2].atomsIncluded = 'N';
			AAs[seqNum-2].SStype = 'L';
			AAs[seqNum-2].gyration = 0;

			hlces.clear();
			sheets.clear();

		}


	}

}
///////////////////////////////////////////////////////////////////////////////////////////////
// rotate the side chain of a given AA (by indx) around a given bond (by p1-p2)..starting from specific atom (by atom indx)
//this function supposes that the structure of the atoms vector in AAs would be side chain atoms are before C and O atoms
// so the AA is reordered so every big atom is followed by its H atom (if any) and C and O atoms are moved to the last
void Protein::rotateSideChain (int AAIndx, int atomIndx, Coordinate p1, Coordinate p2, double degreeAngle)
{
	if ((AAIndx >= 0) && (AAIndx < numOfAA()) && (atomIndx >= 0) & (atomIndx<numOfAtoms(AAIndx)))
	{
		if (degreeAngle)
		{
			double m[4][4];
			buildRotationMatrix(m, p1, p2, -toRadian(degreeAngle));

			int cAtomIndx = getAtomIndx(AAIndx, " C  ");
			for (int i=atomIndx;i<cAtomIndx;i++)		//rotate till u reach C... where after C atom all atoms remain are BB atom
				rotatePoint(m,AAs[AAIndx].atoms[i].coord);

			//initiate whtInCoord variable
			if (AAs[AAIndx].whtInCoord == 'S')
				AAs[AAIndx].whtInCoord = 'N';

			//initiate the gyration
			AAs[AAIndx].gyration = 0;
		}
	}
	else
	{
		string eMsg = "AAindx given (";
		eMsg += toString(AAIndx);
		eMsg += ") is out of range..or atomIndx (";
		eMsg += toString(atomIndx);
		eMsg += ") is out of range";

		errMsg("Protein", "rotateSideChain", eMsg);

	}
}
///////////////////////////////////////////////////////////////////////////////////////////////
// % Rotate number of AAs around a line
// % Given the start indx (AAs vector indx) and end indx...the line and the
//   angle in Degree, the function rotates All AAs
// % startAtomIndx is the indx of the atom in the start AA u want to start rotation from
// % for example if u want to rotate from third to 6th amino acid around phi angle of the third then u call something like
//	 rotate(2,5,AAs[2].indx(CA),N coord, CA coord, angleInDegree)
//////////////////////////////////////////////////////////////////////////////
void Protein::rotate(int startAAIndx,int endAAIndx,int startAtomIndx, Coordinate p1,Coordinate p2, double degreeAngle){


/*
	int tmpNumOfAA = numOfAA();

	//make it starts from 0 if it is lower than 0
	if (startAAIndx<0)
	{
		startAAIndx = 0;
		cout<<"============================ In Protein::rotate =========================="<<endl;
		cout<<"The startAAIndx is modified to be 0"<<endl;
		cout<<"=========================================================================="<<endl;
		putchar(BEEP);
		if (SHOW_ERRORS)
		{
			cout<<"Press any key..."<<endl;
			getchar();
		}
	}

	//make it ends at last AA in the AA list
	if (endAAIndx >= tmpNumOfAA)
	{
		endAAIndx = tmpNumOfAA-1;
		cout<<"================================ In Protein::rotate ======================"<<endl;
		cout<<"The endAAIndx is modified to be "<<endAAIndx<<endl;
		cout<<"=========================================================================="<<endl;
		putchar(BEEP);
		if (SHOW_ERRORS)
		{
			cout<<"Press any key..."<<endl;
			getchar();
		}
	}
*/
	int tmpStartAANumOfAtoms = numOfAtoms(startAAIndx);

	if ((startAAIndx>=0) && (endAAIndx<numOfAA()) && (startAtomIndx >= 0) && (startAtomIndx < tmpStartAANumOfAtoms))
	{
		if (degreeAngle)
		{
			double m[4][4];
			buildRotationMatrix(m, p1, p2, -toRadian(degreeAngle));

			int i;
			//rotate the rest of atoms in the first AA
			for (i=startAtomIndx;i<tmpStartAANumOfAtoms;i++)
				rotatePoint(m,AAs[startAAIndx].atoms[i].coord);

			//rotate atoms of the next AA till endAAindx
			for(i=startAAIndx+1;i<=endAAIndx;i++)
			{
				for (int j=0;j<numOfAtoms(i);j++)
					// Rotate an atom
					rotatePoint(m,AAs[i].atoms[j].coord);

				//initiate the whtInCoord variable
				AAs[i].whtInCoord = 'N';
				//initiate the gyration variable
				AAs[i].gyration = 0;
			}
		}
	}
	else
	{
		string eMsg = "The given index (";
		eMsg += toString(startAtomIndx);
		eMsg += ") is out of range..or the portion is empty (No.AAs= ";
		eMsg += toString(numOfAA());
		eMsg += ")";

		errMsg("Protein", "rotate", eMsg);
	}
}
/////////////////////////////////////////////////////////////////////////////////////////
char Protein::chr3ToChr1(string chr3)
{

	if ( !chr3.compare("ALA"))
	return 'A';
	else if ( !chr3.compare("CYS"))
	return 'C';
	else if ( !chr3.compare("ASP"))
	return 'D';
	else if ( !chr3.compare("GLU") )
	return 'E';
	else if ( !chr3.compare("PHE") )
	return 'F';
	else if ( !chr3.compare("GLY") )
	return 'G';
	else if ( !chr3.compare("HIS") )
	return 'H';
	else if ( !chr3.compare("ILE") )
	return 'I';
	else if ( !chr3.compare("LYS") )
	return 'K';
	else if ( !chr3.compare("LEU") )
	return 'L';
	else if ( !chr3.compare("MET") )
	return 'M';
	else if ( !chr3.compare("ASN") )
	return 'N';
	else if ( !chr3.compare("PRO") )
	return 'P';
	else if ( !chr3.compare("GLN") )
	return 'Q';
	else if ( !chr3.compare("ARG") )
	return 'R';
	else if ( !chr3.compare("SER") )
	return 'S';
	else if ( !chr3.compare("THR") )
	return 'T';
	else if ( !chr3.compare("VAL") )
	return 'V';
	else if ( !chr3.compare("TRP") )
	return 'W';
	else if ( !chr3.compare("TYR") )
	return 'Y';
	else
	{
		//string eMsg = "Unknown character (X) has been returned for AA (" + chr3;
		//eMsg += " )..";

		//errMsg("Protein", "Chr3toChr1", eMsg);

		return 'X';		// Unknown Characters;
	}

}

/////////////////////////////////////////////////////////////////////////////////////////////
//returns the type of the atom
char Protein::getAtomType(string atomName)
{
	size_t found;

	found = atomName.find("N");
	if (found!= atomName.npos)
		return 'N';

	found = atomName.find("C");
	if (found !=atomName.npos)
		return 'C';

	found = atomName.find("O");
	if (found != atomName.npos)
		return 'O';

	found = atomName.find("P");
	if (found != atomName.npos)
		return 'P';

	found = atomName.find("S");
	if (found != atomName.npos)
		return 'S';

	found = atomName.find("H");
	if (found != atomName.npos)
		return 'H';

	errMsg("Protein", "getAtomType", "Unknown atom type.....");

	exit(1);

}
//////////////////////////////////////////////////////////////////////////////////////////
bool Protein::isSideChainAtom(string atomName)
{
	if ((atomName == " N  ") || (atomName == " O  ") || (atomName == " C  ") || (atomName == " CA "))
		return false;

	if ((atomName == " H  ") || (atomName == " H1 ") || (atomName == " H2 ") || (atomName == " H3 "))
		return false;

	if ((atomName == " HA ") || (atomName == " HA1") || (atomName == " HA2") || (atomName == " HA3"))
		return false;
	if (atomName == " OXT")
		return false;

	return true;
}
////////////////////////////////////////////////////////////////////////////////////////////
//method is E (connect or concat from the end) or S (connect or concat from the beginning
void Protein::getCPosition (Protein & movablePortion, char method)
{
	//No need to check for lengths because it should be checked from connect function or concat function

	Coordinate  CaAtomCoord,		//Ca Atom Coordinate of the portion
				NAtomCoord,			//N atom coordinate of the portion
				CaAtomMovable,		//Ca atom coordinate of the movable portion
				NAtomMovable,		//N atom coordinate of the movable portion
				CAtomCoord,			//C atom coordinate of the portion
				CAtomMovableCoord;	//C atom coordinate of the movable portion
	int targetCAtomMovableIndx;		//the indx of C atom in the movable Portion..if E then first C atom ..if S then last C atom

	int tmpNumOfAA = numOfAA();
	int tmpPortionNumOfAA = movablePortion.numOfAA();

	if (method == 'E')
	{
		//get coordinates of the last N-CA bond
		CaAtomCoord = getAtomCoordinate(tmpNumOfAA - 1, " CA ");
		NAtomCoord = getAtomCoordinate(tmpNumOfAA - 1, " N  ");
		//get coordinates of the first N-CA bond of the movable portion
		CaAtomMovable = movablePortion.getAtomCoordinate(0," CA ");
		NAtomMovable = movablePortion.getAtomCoordinate(0," N  ");
		//get the coordinate of the last C atom
		CAtomCoord = getAtomCoordinate(tmpNumOfAA - 1," C  ");
		//get the indx of the first C atom of the movable portion
		targetCAtomMovableIndx = 0;
	}
	else
	{
		//get coordinates of first N-CA bond
		CaAtomCoord = getAtomCoordinate(0, " CA ");
		NAtomCoord = getAtomCoordinate(0, " N  ");
		//get coordinates of the first N-CA bond of the movable portion
		CaAtomMovable = movablePortion.getAtomCoordinate(tmpPortionNumOfAA -1 ," CA ");
		NAtomMovable = movablePortion.getAtomCoordinate(tmpPortionNumOfAA -1 ," N  ");
		//get the coordinate of the first C atom
		CAtomCoord = getAtomCoordinate(0," C  ");
		//get the indx of the last C atom of the movable portion
		targetCAtomMovableIndx = tmpPortionNumOfAA -1;
	}

	//overlap N-CA bonds
	movablePortion.overlapLine(CaAtomMovable,		//CA atom of the first residue	(Movable)
								NAtomMovable,		//N atom of the first residue	(Movable)
								CaAtomCoord,		//CA Atom of the last residue (Fixed)
								NAtomCoord);		//N atom of the last residue (Fixed)

	//Now...N-CA bond of the movablePortion is moved and overlapped with N-CA bond of the portion

	//find the torsion angle value to overlap CA-C bonds
	//first..find the coordinate of C atom of the movable portion
	CAtomMovableCoord = movablePortion.getAtomCoordinate(targetCAtomMovableIndx," C  ");


	Coordinate tmpCoordinate;		// extra point to calculate the torsion angle...good in the case that the portion length is 1 AA

	//find phi for the last AA of the portion
	double targetTorsion = getTorsionAngle (tmpCoordinate,		//Random Point	(Fixed)
											NAtomCoord,			//N Atom of the last residue (Fixed)
											CaAtomCoord,		//CA atom of the last residue (Fixed)
											CAtomCoord);		//C atom of the last residue of the portion (Fixed)



	//find the amount of degrees to be rotated
	double currentTorsion = getTorsionAngle(tmpCoordinate,		//Random Point	(Fixed)
											NAtomCoord,			//N Atom of the last residue (Fixed)
											CaAtomCoord,		//CA atom of the last residue (Fixed)
											CAtomMovableCoord);//C atom of the first residue of movablePortion (Movable)

	//rotate the portionToBeConnected
	movablePortion.rotate(0,tmpPortionNumOfAA -1,0,CaAtomCoord,NAtomCoord, targetTorsion - currentTorsion);
}
////////////////////////////////////////////////////////////////////////////////////////////
AminoAcid Protein::reOrderAtoms (int AAIndx)
{
	if ((AAIndx >= 0) && (AAIndx < numOfAA()))
	{
		int	tmpNumOfAtoms = numOfAtoms(AAIndx);

		AminoAcid tmpAA;
		tmpAA = AAs[AAIndx];
		tmpAA.atoms .clear();
		for (int i=0;i<tmpNumOfAtoms;i++)
		{
			if ((AAs[AAIndx].atoms [i].name != " C  ") && (AAs[AAIndx].atoms [i].name != " O  ") && (AAs[AAIndx].atoms [i].name != " OXT"))
			{
				if (AAs[AAIndx].atoms [i].type != 'H')
				{
					size_t found;


					tmpAA.atoms .push_back (AAs[AAIndx].atoms [i]);		//push the atom
					string suffix = AAs[AAIndx].atoms [i].name .substr (2,2).c_str();
					//if the fourth letter of the atom name is not empty
					if (suffix.substr(suffix.length()-1, 1) != " ")
					{
						for (int j=i+1; j<tmpNumOfAtoms; j++)
						{
							found = AAs[AAIndx].atoms [j].name.find(suffix);

							if ((found != AAs[AAIndx].atoms[j].name.npos) &&
								(AAs[AAIndx].atoms [j].name != " C  ") &&
								(AAs[AAIndx].atoms [j].name != " O  ") &&
								(AAs[AAIndx].atoms [j].name != " OXT"))
							{
								tmpAA.atoms .push_back (AAs[AAIndx].atoms [j]);
							}
						}
					}
					else
						//if the third letter of the atom name is not empty ...for Ca and Cb or other 2 letters atom
						if (suffix.substr(suffix.length()-2,1) != " ")
						{
							for (int j=i+1; j<tmpNumOfAtoms; j++)
							{
								found = AAs[AAIndx].atoms [j].name.substr(2,2).find(suffix.substr(0,1).c_str());

								if ((found != AAs[AAIndx].atoms[j].name.npos) &&
									(AAs[AAIndx].atoms [j].name != " C  ") &&
									(AAs[AAIndx].atoms [j].name != " O  ") &&
									(AAs[AAIndx].atoms [j].name != " OXT"))
								{
									tmpAA.atoms .push_back (AAs[AAIndx].atoms [j]);
								}
							}
						}
						else
							//for N...push all H atoms
							if (suffix == "  ")
							{
								for (int j=i+1; j<tmpNumOfAtoms; j++)
								{
									if ((AAs[AAIndx].atoms [j].name == " H  ") || (AAs[AAIndx].atoms [j].name == " H1 ")
										|| (AAs[AAIndx].atoms [j].name == " H2 ") || (AAs[AAIndx].atoms [j].name == " H3 "))
										tmpAA.atoms .push_back (AAs[AAIndx].atoms [j]);

								}
							}
							else
							{
								errMsg("Protein", "reOrderAtoms", "Error occured in the name of the AA " + AAs[AAIndx].chr3);
								exit(1);

							}

				}
			}
		}

		//Always...the last 3 atoms are C, O, and OXT
		int cAtomIndx = getAtomIndx(AAIndx," C  ");
		if (cAtomIndx != -1)
			tmpAA.atoms.push_back(AAs[AAIndx].atoms[cAtomIndx]);		//push C atom
		int oAtomIndx = getAtomIndx(AAIndx," O  ");
		if (oAtomIndx != -1)
			tmpAA.atoms.push_back(AAs[AAIndx].atoms[oAtomIndx]);		//push O atom
		int oxtIndx = getAtomIndx(AAIndx," OXT");
		if ( oxtIndx!= -1)
			tmpAA.atoms.push_back(AAs[AAIndx].atoms[oxtIndx]);						//push OXT if exist

		return tmpAA;
	}
	else
	{
		string eMsg = "The given index (";
		eMsg += toString(AAIndx);
		eMsg += ") is out of range..or the portion is empty (No.AAs= ";
		eMsg += toString(numOfAA());
		eMsg += ")";

		errMsg("Protein", "reOrderAtoms", eMsg);

		exit(1);

	}
}
///////////////////////////////////////////////////////////////////////////////////////////
int Protein::getNumOfChi(char aa1chr){
//takes the aa in one character mode

		if ( aa1chr == 'A' ) return 0;	//ALA
		if ( aa1chr == 'R' ) return 4;	//ARG
		if ( aa1chr == 'N' ) return 2;	//ASN
		if ( aa1chr == 'D' ) return 2;	//ASP
		if ( aa1chr == 'C' ) return 1;	//CYS
		if ( aa1chr == 'Q' ) return 3;	//GLN
		if ( aa1chr == 'E' ) return 3;	//GLU
		if ( aa1chr == 'G' ) return 0;	//GLY
		if ( aa1chr == 'H' ) return 2;	//HIS
		if ( aa1chr == 'I' ) return 2;	//ILE
		if ( aa1chr == 'L' ) return 2;	//LEU
		if ( aa1chr == 'K' ) return 4;	//LYS
		if ( aa1chr == 'M' ) return 3;	//MET
		if ( aa1chr == 'F' ) return 2;	//PHE
		if ( aa1chr == 'P' ) return 2;	//PRO
		if ( aa1chr == 'S' ) return 1;	//SER
		if ( aa1chr == 'T' ) return 1;	//THR
		if ( aa1chr == 'W' ) return 2;	//TRP
		if ( aa1chr == 'Y' ) return 2;	//TYR
		if ( aa1chr == 'V' ) return 1;	//VAL

		if ( aa1chr == 'X' )
		{
			errMsg("Protein", "getNumOfChi", "Unknown character (X) has been reached..no chi angles will be returned");
		}

		return 0; //Unknown

}
///////////////////////////////////////////////////////////////////////////////////////////////////////
bool Protein::isHeavyAtom(char type)
{
	if ((type == 'C') || (type == 'N') || (type == 'O') || (type == 'S'))
		return true;

	return false;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
void Protein::sortHlces(vector<HelicesSecondaryStructure> & a, const unsigned int leftArg, unsigned int rightArg)
{
  if (leftArg < rightArg)
  {

    int pivotvalue = a[leftArg].startIndx;
    int left = leftArg - 1;
    int right = rightArg + 1;

  for(;;)
  {

    while (a[--right].startIndx > pivotvalue);
    while (a[++left].startIndx < pivotvalue);

    if (left >= right) break;

    HelicesSecondaryStructure temp = a[right];
    a[right] = a[left];
    a[left] = temp;
  }

  int pivot = right;
  sortHlces(a, leftArg, pivot);
  sortHlces(a, pivot + 1, rightArg);
  }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////
void Protein::sortStrands(vector<SheetsSecondaryStructure> & a, const unsigned int leftArg, unsigned int rightArg)
{
  if (leftArg < rightArg)
  {

    int pivotvalue = a[leftArg].startIndx;
    int left = leftArg - 1;
    int right = rightArg + 1;

  for(;;)
  {

    while (a[--right].startIndx > pivotvalue);
    while (a[++left].startIndx < pivotvalue);

    if (left >= right) break;

    SheetsSecondaryStructure temp = a[right];
    a[right] = a[left];
    a[left] = temp;
  }

  int pivot = right;
  sortStrands(a, leftArg, pivot);
  sortStrands(a, pivot + 1, rightArg);
  }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////// END OF PROTEIN CLASS IMPLEMENTATION ////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// IMPLEMENTATION ///////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// copy a specified range of backbone AAs from a portion (originalPortion) to another portion
void copyBBone(Protein originalPortion, int startIndx, int endIndx, Protein & newPortion)
{

	if ((startIndx >= 0 ) && (endIndx < originalPortion.numOfAA()))
	{
		int AACntr = 0, i;
		newPortion.path  = originalPortion.path;		//take the source file path
		newPortion.ID = originalPortion.ID;				//take the source file ID
		newPortion.header = originalPortion.header;		//copy the original header lines
		newPortion.AAsCollide.clear ();					//reset the list of two AA's collide
		newPortion.atomsCollide .clear ();				//reset the list of two atoms collide

		//copy hlces in the range specified...if any
		//similar operation should be made for sheets
		for (i=0;i<originalPortion.hlces.size();i++)
		{
			//The whole hlx is within the range
			if ((startIndx <= originalPortion.hlces[i].startIndx) && (endIndx >= originalPortion.hlces[i].endIndx))
			{
				//update indeces
				originalPortion.hlces [i].startIndx -=  startIndx;
				originalPortion.hlces [i].endIndx	-= startIndx;

				//update the length
				originalPortion.hlces [i].nAAcur = originalPortion.hlces [i].endIndx - originalPortion.hlces [i].startIndx + 1;

				newPortion.hlces .push_back (originalPortion.hlces[i]);
				continue;
			}
			//both end are within the range
			if ((startIndx > originalPortion.hlces [i].startIndx ) && (endIndx < originalPortion.hlces [i].endIndx))
			{
				//new indeces
				originalPortion.hlces [i].startIndx = startIndx - originalPortion.hlces [i].startIndx;
				originalPortion.hlces [i].endIndx	= endIndx - startIndx;

				//update the length
				originalPortion.hlces [i].nAAcur = originalPortion.hlces [i].endIndx - originalPortion.hlces [i].startIndx + 1;

				newPortion.hlces .push_back (originalPortion.hlces [i]);
				continue;


			}
			//The lower end is within the range
			if ((startIndx<= originalPortion.hlces[i].startIndx) && (endIndx >= originalPortion.hlces[i].startIndx) && (endIndx < originalPortion.hlces[i].endIndx))
			{
				//new indeces
				originalPortion.hlces [i].startIndx -=  startIndx ;
				originalPortion.hlces [i].endIndx    =	endIndx - startIndx;

				//update the length
				originalPortion.hlces [i].nAAcur = originalPortion.hlces [i].endIndx - originalPortion.hlces [i].startIndx + 1;

				//push the new hlx
				newPortion.hlces .push_back (originalPortion.hlces [i]);
				continue;
			}

			//The upper end is within the range
			if ((startIndx > originalPortion.hlces[i].startIndx) && (startIndx <= originalPortion.hlces[i].endIndx) && (endIndx >= originalPortion.hlces[i].endIndx))
			{
				//new indeces
				originalPortion.hlces [i].endIndx	 -= startIndx;
				originalPortion.hlces [i].startIndx  = startIndx - originalPortion.hlces [i].startIndx;

				//update the length
				originalPortion.hlces [i].nAAcur = originalPortion.hlces [i].endIndx - originalPortion.hlces [i].startIndx + 1;

				//push new hlx
				newPortion.hlces .push_back (originalPortion.hlces [i]);
				continue;
			}
		}

		newPortion.sheets = originalPortion.sheets;		//take the sheets...in the future..should be done like hlces

		for (i= startIndx; i<=endIndx; i++)
		{
			//copy the AA
			newPortion.AAs.push_back(originalPortion.AAs[i]);

			//clear All old atoms
			newPortion.AAs[AACntr].atoms.clear();

			//push all backbone atoms for the i's AA
			for (int j=0;j<originalPortion.numOfAtoms(i);j++)
				if (!originalPortion.AAs[i].atoms[j].isSideChain)
					newPortion.AAs[AACntr].atoms.push_back(originalPortion.AAs[i].atoms[j]);

			//reset radius of gyration
			newPortion.AAs[AACntr].gyration = 0;
			//reset wht in Coord
			newPortion.AAs [AACntr].whtInCoord = 'N';

			AACntr++;
		}
	}
	else
	{
		string eMsg = "The given range (";
		eMsg += toString(startIndx);
		eMsg += "-";
		eMsg += toString(endIndx);
		eMsg += ") is out of range..or the portion is empty (No.AAs= ";
		eMsg += toString(originalPortion.numOfAA());
		eMsg += ")";

		errMsg("Protein", "copyBBone", eMsg);

	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// End of IMPLEMENTATION ///////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif
