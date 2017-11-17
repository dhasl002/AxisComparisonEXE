#ifndef UTILITYFUNCTIONS_H
#define UTILITYFUNCTIONS_H

#include <sstream>
#include <algorithm>		//used mainly to generate permutations
#include <vector>
#include <iostream>

#include "constants.h"

///////////////////////////////////////////////////////////////////////////////// for random number generator //////////////
#include "randomGenerator/randomc.h"
// This template class combines two different random number generators
// for improved randomness. R1 and R2 are any two different random number
// generator classes.
template <class RG1, class RG2>
class TRandomCombined : private RG1, private RG2 {
public:
   TRandomCombined(int seed) : RG1(seed), RG2(seed+1) {};

   void RandomInit(int seed) {        // re-seed
      RG1::RandomInit(seed);
      RG2::RandomInit(seed+1);
   }

   double Random() {
      double r = RG1::Random() + RG2::Random();
      if (r >= 1.) r -= 1.;
      return r;
   }

   int IRandom(int min, int max){       // output random integer
      // get integer random number in desired interval
      int iinterval = max - min + 1;
      if (iinterval <= 0) return 0x80000000; // error
      int r = int(iinterval * Random()); // truncate
      if (r >= iinterval) r = iinterval-1;
      return min + r;
   }
};
#include "randomGenerator/userintf.cpp"
#include "randomGenerator/mersenne.cpp"
#include "randomGenerator/mother.cpp"
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef _WIN
#include <time.h>					//used to find elapsed time
#else
#include <sys/time.h>
#endif

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// List Of Data Structures //////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct OnePermutation
{//used to store one permutation
	vector<int> permutation;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// End Of Constants and Data Structures /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// List of Functions ////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

vector<OnePermutation> getAllPermutations (vector<OnePermutation> & allCombinations);	//generate all permutations for a combinations....
vector<OnePermutation> getAllCombinations (int nItems, int nChoices);			//choose R from N items
double nCombination(short n, short k);						//returns the number of combinations
double getCombIndx(short V, short K, short* Num);	//given a combination...return the rank (indx) among n choose k entries
void printPermutations(vector<OnePermutation> allPermutations);	//print the permutation
double getRadius (char atomType,char method = 'O');				//given the type of atom and the method ... return the radius of this atom
string num2Chr3 (int AARefNum);			//convert from integer representation of an AA to character representation ( 3 letters)
int chr2Num (char chr1);				//convert from character (1 letter) representation to integer representation of an AA
void errMsg(string className, string methodName, string msg, bool important = false);		//errMsg
// Generate a random integer with a specified range
//#define getRandom(min, max) \
 //   ((rand()%(int)(((max) + 1)-(min)))+ (min))
/*
int getRandom(int min, int max)//  linear congruential pseudo random number generator based on D. Lehmer
{
	static unsigned long a = 2007, b = 4194301, c = 2147483647, z = b;
	if ( max < 0 )
	{
		max = -max;
		a = max;
	}
	z = (a + b * z) % c;
	return z % max + min;
}
*/
int getRandom (int minpossible, int maxpossible){
    //time_t seconds;

    //time(&seconds);
    //Generate random number in a arange
    //srand((unsigned int) ( seconds * seed));
	return minpossible + rand() % (maxpossible - minpossible + 1);
}
float getRandomFloat(float min, float max, int seed = 1){

    TRandomCombined<CRandomMersenne,CRandomMother> RG(seed);

    return (RG.IRandom((int) min, (int) max-1) + RG.Random());
}
void getRandomList(float min, float max, float *rList, int nList, int seed = 1){

    seed *= seed;
    TRandomCombined<CRandomMersenne,CRandomMother> RG(seed);

    for (int i=0; i<nList; i++)
        rList[i] = RG.IRandom((int) min, (int) max-1) + RG.Random();
}
void getRandomList(short min, short max, short *rList, short nList, int seed = 1){

    //seed *= (int)time(0);
    seed *= seed;
    TRandomCombined<CRandomMersenne,CRandomMother> RG(seed);

    for (int i=0; i<nList; i++)
        rList[i] = RG.IRandom(min, max);
}

double factorial(short num);
template<class T> inline void swap(T& v1,T& v2);							//perform template swap
template <class T> inline std::string toString (const T& t);				//convert a data type to string
template <class T> inline T round(const T& t, short places);					//round the number to ith tenth - to decimal whole number, 1 to the tenth,..so on
template <typename T> T **AllocateDynamicArray(int nRows, int nCols);		//Allocate 2 dimensional array
template <typename T> void FreeDynamicArray(T** dArray);					//dispose memory of 2 dimenstional array
#define MAX(x, y)   ((x) >= y ? (x) : (y))				// return max
#define MIN(x, y)	((x) < y ? (x) : (y))				//return min
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// ////////// End Of The List //////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// IMPLEMENTATION ///////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<OnePermutation> getAllPermutations(vector<OnePermutation> & allCombinations)
{

//	vector<int> v;
	vector<OnePermutation> All_permutations;

	int i;

	for (i=0; i<allCombinations.size (); i++)			//for each combination...find it's all permutations
	{

		//OnePermutation tmp;
		//tmp.permutation = allCombinations[i].permutation;		//initialize first permutations
		All_permutations.push_back(allCombinations[i]);

		while (next_permutation(allCombinations[i].permutation.begin(), allCombinations[i].permutation.end() ) )
		{
		// Loop until all permutations are generated.
			//tmp.permutation = v;
			All_permutations.push_back(allCombinations[i]);
		}

	}
	return All_permutations;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<OnePermutation> getAllCombinations (int nItems, int nChoices)
{

	OnePermutation c;
	vector<OnePermutation> allCombinations;
	int i, k;

	// create initial combination (0..nChoices - 1)
    for (i=0; i<nChoices; i++)
		c.permutation.push_back(i+1);
    do
    {
		i = 0;
		k = 1;

		allCombinations.push_back (c);		//commit combination

		while ( i <nChoices-1 && c.permutation[i+1] == c.permutation[i]+1) // for each bump
		{
			c.permutation[i] = k++;                 // fall back
			i++;
		}

		c.permutation[i]++;		//push forward
    }
    while (nItems - c.permutation[i] + 1);		//verify

	return allCombinations;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double nCombination(short n, short k){
	short jC;
	double va=1;

	for (jC=n-k+1;jC<=n; jC++)
		va = va*jC;

	return va/factorial(k);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////
double getCombIndx(short V, short K, short* Num){
//	int V=5;
//	int K=2;
	vector<short> C(K,0);
//	int Num[2] = {2,3};
	double R;

	double LI=0, P1 = K-1;

	if (K == 1)
		return Num[0];

	for (int i=0;i<P1;i++){
		C[i] = 0;
		if (i!=0)
			C[i] = C[i-1];

		C[i]++;
		R = nCombination(V-C[i], K-i-1);
		LI = LI + R;
		//cout<<"R= "<<R<<endl;

		while (C[i] < Num[i]){
			C[i]++;
			R = nCombination(V-C[i], K-i-1);
			LI = LI + R;
			//cout<<"while : "<<C[i]<<" "<<Num[i]<<" V-C[i]= "<<V-C[i]<<" R= "<<R<<endl;
		}
		LI = LI - R;
		//cout<<"   LI= "<<LI<<endl;
	}
	return (LI + Num[K-1] - Num[(int) (P1-1)]);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void printPermutations(vector<OnePermutation> allPermutations)
{
	for (int i=0; i<allPermutations.size();i++)
	{
		cout<<i+1<<" - ";
		for (int j=0; j<allPermutations[i].permutation .size();j++)
			cout<<allPermutations[i].permutation [j]<<" ";
		cout<<endl;
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int chr2Num (char chr1)
{
	if ( chr1 == 'A' ) return 0;	//ALA
	if ( chr1 == 'R' ) return 1;	//ARG
	if ( chr1 == 'N' ) return 2;	//ASN
	if ( chr1 == 'D' ) return 3;	//ASP
	if ( chr1 == 'C' ) return 4;	//CYS
	if ( chr1 == 'Q' ) return 5;	//GLN
	if ( chr1 == 'E' ) return 6;	//GLU
	if ( chr1 == 'G' ) return 7;	//GLY
	if ( chr1 == 'H' ) return 8;	//HIS
	if ( chr1 == 'I' ) return 9;	//ILE
	if ( chr1 == 'L' ) return 10;	//LEU
	if ( chr1 == 'K' ) return 11;	//LYS
	if ( chr1 == 'M' ) return 12;	//MET
	if ( chr1 == 'F' ) return 13;	//PHE
	if ( chr1 == 'P' ) return 14;	//PRO
	if ( chr1 == 'S' ) return 15;	//SER
	if ( chr1 == 'T' ) return 16;	//THR
	if ( chr1 == 'W' ) return 17;	//TRP
	if ( chr1 == 'Y' ) return 18;	//TYR
	if ( chr1 == 'V' ) return 19;	//VAL

	if ( chr1 == 'X' )
	{
		cout<<"========================= in chr2Num(char) ==========================="<<endl;
		cout<<"Unknown character (X) has been reached..int number 20 will be returned"<<endl;
		cout<<"======================================================================"<<endl;
	}

	if (SHOW_ERRORS)
	{
		putchar(BEEP);
		cout<<"Press any key..."<<endl;
		getchar();
	}

	return 20;	//Unknown Character

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//error msg
// given the name of the class, method, and the error msg, print out the msg on the screen
//if the error msg is important then the program should stop working untill the user hit any key
void errMsg(string className, string methodName, string msg, bool important){

	cout<<endl;
	cout<<"================================ in "<<className<<"."<<methodName<<"() ===================="<<endl;
	cout<< msg<< endl;
	cout<<"==============================================================================="<<endl;
	if (SHOW_ERRORS || important)
	{
		putchar(BEEP);
		cout<<"Press any key..."<<endl;
		getchar();
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
string num2Chr3 (int AARefNum)
{
	if ( AARefNum == 0 ) return "ALA";	//ALA
	if ( AARefNum == 1 ) return "ARG";	//ARG
	if ( AARefNum == 2 ) return "ASN";	//ASN
	if ( AARefNum == 3 ) return "ASP";	//ASP
	if ( AARefNum == 4 ) return "CYS";	//CYS
	if ( AARefNum == 5 ) return "GLN";	//GLN
	if ( AARefNum == 6 ) return "GLU";	//GLU
	if ( AARefNum == 7 ) return "GLY";	//GLY
	if ( AARefNum == 8 ) return "HIS";	//HIS
	if ( AARefNum == 9 ) return "ILE";	//ILE
	if ( AARefNum == 10 ) return "LEU";	//LEU
	if ( AARefNum == 11 ) return "LYS";	//LYS
	if ( AARefNum == 12 ) return "MET";	//MET
	if ( AARefNum == 13 ) return "PHE";	//PHE
	if ( AARefNum == 14 ) return "PRO";	//PRO
	if ( AARefNum == 15 ) return "SER";	//SER
	if ( AARefNum == 16 ) return "THR";	//THR
	if ( AARefNum == 17 ) return "TRP";	//TRP
	if ( AARefNum == 18 ) return "TYR";	//TYR
	if ( AARefNum == 19 ) return "VAL";	//VAL

	cout<<"============================= in num2Chr3(char) ======================"<<endl;
	cout<<"Unknown number has given.... XXX AA will be returned"<<endl;
	cout<<"======================================================================"<<endl;

	if (SHOW_ERRORS)
	{
		putchar(BEEP);
		cout<<"Press any key..."<<endl;
		getchar();
	}
	return "XXX";  //Unknown Character
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//method V = VDW;  C= calculated; O= Covalent (default)
double getRadius(char atomType,char method)
{
	/*
	// covalent bond radius (2008 values) /www.webelements.com
	#define RADIUS_C        0.76
	#define RADIUS_N        0.71
	#define RADIUS_O        0.66
	#define RADIUS_H        0.31
	#define RADIUS_S        1.03
	#define RADIUS_P        1.07

	  // calculated radius  //www.webelements.com
	#define RADIUS_C        0.67
	#define RADIUS_N        0.56
	#define RADIUS_O        0.48
	#define RADIUS_H        0.53
	#define RADIUS_S        0.83
	#define RADIUS_P        0.98


	// van der waals force radius
	#define RADIUS_C        1.7
	#define RADIUS_N        1.55
	#define RADIUS_O        1.52
	#define RADIUS_H        1.2
	#define RADIUS_S        1.85
	#define RADIUS_P        1.9

	*/
	switch (method){
		case 'V': switch (atomType){
						/*
						case 'N': return 1.55;
						case 'C': return 1.7;
						case 'O': return 1.52;
						case 'P': return 1.9;
						case 'S': return 1.85;
						case 'H': return 1.2;
						*/
						/*
						case 'N': return 1.5;
						case 'C': return 1.7;
						case 'O': return 1.4;
						case 'P': return 1.9;
						case 'S': return 1.8;
						case 'H': return 1.0;
						*/
						///
						///			as in SCWRL 3.0 Dunbrack paper
						///
						case 'N': return 1.3;
						case 'C': return 1.6;
						case 'O': return 1.3;
						case 'P': return 1.8;
						case 'S': return 1.7;
						case 'H': return 1.0;

					}
		case 'C': switch (atomType){
						/*
						case 'N': return 0.56;
						case 'C': return 0.67;
						case 'O': return 0.48;
						case 'P': return 0.98;
						case 'S': return 0.83;
						case 'H': return 0.53;
						*/
						//Dr. Weitao
						case 'N': return 0.70;
						case 'C': return 0.77;
						case 'O': return 0.66;
						case 'P': return 1.1;
						case 'S': return 1.04;
						case 'H': return 0.32;

				  }
		case 'O': switch (atomType){

						case 'N': return 0.71;
						case 'C': return 0.76;
						case 'O': return 0.66;
						case 'P': return 1.07;
						case 'S': return 1.03;
						case 'H': return 0.31;

			/*
						// Dr. Weitao's values
						case 'N': return 0.73;
						case 'C': return 0.77;
						case 'O': return 0.74;
						case 'P': return 1.1;
						case 'S': return 1.03;
						case 'H': return 0.3;
			*/
				  }
	}

	cout<<"============================= in getRadius(char, char) ==============="<<endl;
	cout<<"The atomtype or method class sent are incorrect"<<endl;
	cout<<"======================================================================"<<endl;
	exit(1);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
double factorial(short num)
{
	 double result=1;
	 for (short i=1; i<=num; i++)
		result *= i;
	 return result;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class T>
inline void swap(T& v1,T& v2)
{
  T temp=v2;
  v2=v1;
  v1=temp;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <class T>
inline std::string toString (const T& t)
{
std::stringstream ss;
ss << t;
return ss.str();
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <class T> inline T round(const T& t, short places){
	T off= pow(10, places);
	return int(t*off + (t<0? -0.5 : 0.5))/off;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Allocate 2 dimentional array
template <typename T>
T **AllocateDynamicArray( int nRows, int nCols)
{
      T **dynamicArray;

      dynamicArray = new T*[nRows];
      for( int i = 0 ; i < nRows ; i++ )
      dynamicArray[i] = new T [nCols];

      return dynamicArray;
}

//Dispose Memeory
template <typename T>
void FreeDynamicArray(T** dArray)
{
      delete [] *dArray;
      delete [] dArray;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// End Of IMPLEMENTATION ////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif


