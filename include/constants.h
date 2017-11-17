#ifndef CONSTANTS_H
#define CONSTANTS_H

//#define _WIN              //if work on windows

#define _MPI_RUN			//uncomment it if u wanna run the program using MPI

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////// list Of Constants and Data Structures /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//All possible atom names could be found in a pdb file

/*
ALL_ATOMS_NAMES
		" N  "," CA "," C  "," O  "," CB "," CG "," CG1"," CG2",
		" CD "," CD1"," CD2"," CE "," CE1"," CE2"," CE3"," CZ "," CZ1"," CZ2"," CZ3"," CH2",		//20 atoms
		" ND "," ND1"," ND2"," NZ "," NZ1"," NZ2"," NE "," NE1"," NE2"," NH "," NH1"," NH2",		//12 atoms
		" OD "," OD1"," OD2"," OG "," OG1"," OG2"," OE "," OE1"," OE2"," OH "," OH1"," OH2"," OXT",	//13 atoms
		" SG "," SD ",																				//2 atoms
		//H atoms"
		" H  ","1H  " ,"2H  ","3H  "," HA ","1HA " ,"2HA " ," HB ","1HB ","2HB ", "3HB ",
		" HG ","1HG ", "2HG ","3HG "," HG1","1HG1", "2HG1", "3HG1"," HG2","1HG2", "2HG2", "3HG2",	//23 atoms
		" HD ","1HD ", "2HD ","3HD "," HD1","1HD1", "2HD1", "3HD1"," HD2","1HD2", "2HD2", "3HD2",
		" HE ","1HE ", "2HE ","3HE "," HE1","1HE1", "2HE1", "3HE1"," HE2","1HE2", "2HE2", "3HE2",	//24 atoms
		" HZ ","1HZ ", "2HZ ","3HZ "," HZ1","1HZ1", "2HZ1", "3HZ1"," HZ2","1HZ2", "2HZ2", "3HZ2",
		" HH ","1HH ", "2HH ","3HH "," HH1","1HH1", "2HH1", "3HH1"," HH2","1HH2", "2HH2", "3HH2"	//24 atoms
*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//define bond angles
#define C_N_Ca_ANGLE	123
#define O_C_N_ANGLE		125
#define Ca_C_O_ANGLE	121
#define Ca_C_N_ANGLE	114
#define H_N_Ca_ANGLE	114

//define bond lengths
#define Ca_C_BOND_LENGTH	1.53
#define N_Ca_BOND_LENGTH	1.47
#define C_N_BOND_LENGTH		1.32
#define C_O_BOND_LENGTH		1.24
//////////////////////////////////////////////////////////////////////////////////////
// The ranges of allowed torsion angles for Hlx...
#define MIN_PHI_HLX_ALLOWED		-80					//was -110
#define MAX_PHI_HLX_ALLOWED		-40					//was -60
#define MIN_PSI_HLX_ALLOWED		-60					//was -70
#define MAX_PSI_HLX_ALLOWED		-20					//was 30

// The ranges of allowed torsion angles for Strand.
#define MIN_PHI_BETA_ALLOWED	-170
#define MAX_PHI_BETA_ALLOWED	-60
#define MIN_PSI_BETA_ALLOWED	 90
#define MAX_PSI_BETA_ALLOWED	 175

// The ranges of allowed torsion angles for Loop.
#define MIN_PHI_LOOP_ALLOWED	-170 		//was -175
#define MAX_PHI_LOOP_ALLOWED	 -40//170		//was -40
#define MIN_PSI_LOOP_ALLOWED	-60//-170		//was -60
#define MAX_PSI_LOOP_ALLOWED	 170		//was 175
//////////////////////////////////////////////////////////////////////////////////////
//										phi		psi
//ALPHA HLX (right-handed)				-57		-47			or Approx. -60 and -50
#define ALPHA_HLX_PHI_RIGHT -57
#define ALPHA_HLX_PSI_RIGHT -47

#define ALPHA_RISE	1.5

//Alpha  hlx (left-handed )				 57		 47
#define ALPHA_HLX_PHI_LEFT 57
#define ALPHA_HLX_PSI_LEFT 47


//3-10   hlx (right-handed)				-49     -26
#define HLX_3_10_PHI	-49
#define HLX_3_10_PSI	-26

#define RISE_3_10	2.0

//Pi     hlx (right-handed)				-57		-80
#define PI_HLX_PHI	-57
#define PI_HLX_PSI	-80

#define PI_RISE		1.15


//TypeII hlx (left-handed)				-79		 150
#define TYPEii_HLX_PHI	-79
#define TYPEii_HLX_PSI	150

//Collagen   (right-handed)				-51		 153
#define COLLAGEN_HLX_PHI	-51
#define COLLAGEN_HLX_PSI	153


//betasheet phi and psi angles
#define BETA_SHEET_PHI	-139		//or -130 -->	as from book "Protein structure and function" By Gregory A. Petsko, Dagmar Ringe
#define BETA_SHEET_PSI	135			//or 125	-->
#define BETA_RISE		3			//as from the same book

#define LOOP_RISE	3.8

#define STANDARD_AMINO_ACIDS	20			//the number of standard AA

#define SHOW_ERRORS 0		//1 : if you want the program to show you errors. 0: if you won't


#define BEEP 007			//the code of the beep sound of the internal speaker of the computer

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// End Of Constants and Data Structures /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
