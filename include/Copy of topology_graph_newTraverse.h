#ifndef	TOPOLOGY_GRAPH_H
#define TOPOLOGY_GRAPH_H

#include <vector>
#include <queue>
#include "skeleton_overall.h"

using namespace std;


/////////// DATA STRUCTURE /////////////////////////////////////
struct node
{
	short row;								//represent the index of the SS 
	short col;								//represent the index of the stick

	vector<bool> allowed;					//to check that this SS has assigned before to this stick
	vector<short> path;						//The path (list of SS chosen) ... represent the topology
	short curIndx;							//the link to the next level
	float w;								//toTal weight

	short nStickSolved;						//number of sticks solved so far...used to terminate

	//Initializer
	node() : row(-1), col(-1), curIndx(-1), w(0), nStickSolved(0) {}
};

struct loopLink			//represent the correspondance between the loop and the link b/w two sticks
{
	short rowIndx;
	short colIndx;
	float w;

	loopLink () : rowIndx(-1), colIndx(-1), w(999.9) {}
};

struct cell
{
	node	parent;			//you can instead...use 3 variables...col, row and parent indx
							//if no parent then row = -1 and col = -1 parent Indx = -1

	Protein skeleton;		//the structure of the SS corresponding to this stick

//	short pRowIndx;			//row indx of parent node
//	short pColIndx;			//col. indx of the parent
//	short pLinkIndx;		//the indx of link inside the cell

	int validityRight;			//the maximum size of shift could be applied to the sequence to the right...
	int validityLeft;			//the maximum size of shift could be applied to the sequence to left...
	int startAAindx;			//the start indx of AAs to be assigned to this stick

	Coordinate	cTerminal,		//c terminal of the initial structure
				nTerminal;		//n terminal of the initial structure

	float corrW;			//corresponding weight of assigning this SS to this stick

	vector<loopLink> links;		//links 
//	vector<short> links;	//nodes from the next level connected with this node...the stored value represent the indx of the coloum in G
							//so the index of the linked node is (i+1, links)
//	vector<float> w;		//list of links weights... corresponding weights of links data structure


	//Initializer
	cell() : /* pRowIndx(-1), pColIndx(-1), pLinkIndx(-1), */ validityRight(0), validityLeft(0), startAAindx (-1), corrW(999.0) {}
};

struct paths
{
	vector<short> col;
	float w;				//total weight

	paths(): w(0) {}
};

// Determine priority.
bool operator<(const paths &a, const paths &b)
{
  return a.w   > b.w ;
}

vector<vector<cell> > graph;				//The graph represents topology problem

/////////// DATA STRUCTURE ///////////////////////////////////////////

//////////////////////////// LIST OF FUNCTIONS ///////////////////////
//////////////////////////////////////////////////////////////////////

inline
void getChild(vector<vector<cell> > & g, node current, node & child)
{

	//I will not check if it has children or not...because it should have if this has been called
	child = current;

	//child.col = g[current.row][current.col].links [current.curIndx];
	child.col = g[current.row ][current.col ].links [current.curIndx].colIndx;
	//child.row += 1;//current.row + 1;
	child.row = g[current.row ][current.col ].links [current.curIndx].rowIndx;

	//child.w += g[current.row ][current.col ].w [current.curIndx];
	child.w += g[current.row ][current.col].links [current.curIndx ].w;

	//add corresponding weight
	child.w += g[current.row][current.col].corrW;


	child.allowed [child.col] = 1;
	child.path[child.col] = child.row;

	//disable the twin SS
	if (child.col % 2 == 0)
	{
		child.allowed [child.col + 1] = 1;
		child.path [child.col+1] = 0;
	}
	else
	{
		child.allowed [child.col - 1] = 1;
		child.path [child.col - 1] = 0;
	}

	
	child.curIndx = 0;
	child.nStickSolved += 1;

	g[child.row][child.col].parent = current;
}

inline
void getParent(vector<vector<cell> > & g, node current, node & parent)
{
	parent = g[current.row][current.col].parent;
	parent.curIndx += 1;

	//cout<<"parent.row= "<<parent.row<<" row= "<<parent.col<<" subindx= "<<parent.curIndx<<" w= "<<parent.w<<endl;
}
inline
node traverse(vector<vector<cell> > & g, priority_queue<paths> & topK, node current, bool & done)
{
	node child, parent;

	//cout<<"cur.row= "<<current.row<<"  cur.Col= "<<current.col<<" cur.indx= "<<current.curIndx<<"  cur.w= "<<current.w<<endl;
	//getchar();
	if (current.curIndx >= g[current.row][current.col].links .size())
	{
		//cout<<"no more links.."<<endl;
		//this happens when no node remains in this subtree and need to look to the parent
		if (g[current.row][current.col].parent.row != -1)		//it is not the root
		{

			//cout<<"get parent..no more links"<<endl;
			getParent(g, current, parent);

			return parent;
		}
		else
		{
			//stop traversing...u r done
			done = true;
			return current;
		}

	}
	else
	{
		//check if the link goes to end node
		if (g[current.row][current.col].links [current.curIndx].rowIndx == g.size ())
		{
			//cout<<"Goes to end.. "<<g.size ()<<endl;
			//if the link goes to end node...then check the number of sticks assigned so far \
			  if all sticks assigned you are done...other wise select sibling of current node
			if (current.nStickSolved > current.allowed .size () / 2 - 1)
			{
				//cout<<"    Done...found a path.."<<current.path .size()<<" path : ";
				paths tmp;

				tmp.col = current.path;
				tmp.w = current.w;

				//tmp.w += g[current.row][current.col ].w [current.curIndx];
				tmp.w += g[current.row ][current.col].links [current.curIndx ].w;
				//tmp.col [current.row+1] = g[current.row][current.col ].links [current.curIndx];
				tmp.col [current.col] = current.row;

				if (current.col % 2 == 0)
					tmp.col [current.col + 1] = 0;
				else
					tmp.col [current.col-1] = 0;

				topK.push(tmp);	

				//for(int i=0; i<tmp.col .size (); i++)
				//	cout<<tmp.col [i]<<" ";
				//cout<<endl;

				//getchar();
				//cout<<tmp.w<<endl<<endl;

				current.curIndx ++;
				
				return current;
			}
			else
			{
				//cout<<"not Done yet.."<<endl;
				//get sibling
				if(current.curIndx < g[current.row][current.col].links.size() - 1)
				{
					//cout<<"get sibling.. not Done yet"<<endl;
					current.curIndx += 1;
					return current;
				}
				else
				{
					//cout<<"get parent...not done yet.."<<endl;
					getParent(g, current, parent);
					return parent;
				}
			}
		}
		else
		{
			
			//check if the link goes to valid assignment
			if (current.allowed [g[current.row][current.col].links [current.curIndx].colIndx] == 0)
			{
				//cout<<"getchild...."<<endl;
				getChild(g, current, child);
				return child;
			}
			else
			{
				//cout<<"the assignment to go to is not valid.."<<endl;
				//get sibling
				if(current.curIndx < g[current.row][current.col].links.size() - 1)
				{
					//cout<<"get sibling..NV Assignment"<<endl;
					current.curIndx += 1;
					return current;
				}
				else
				{
					//cout<<"get parent... NV assignment.."<<endl;
					getParent(g, current, parent);
					return parent;
				}
			}
		}
	}

	//getchar();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// hlces should be assigned to sticks first...sticks.ssAssigned should be set before
//this function is supposed to use to determine the maximum of shift could be applied to the stick hlx over the native sequence
//validity right means that the stick hlx could or not be shifted to the right (toward the end of the sequence)...
//validity left means that the stick hlx could be or not be shifted to the left (toward the beginning of the sequence)
void setShifts(vector<SecondaryStruct> SSs, vector<vector<cell> > & g, int nativePDBNumOfAA)
{
	int i, 
		j,
		nSticks = g[0].size()/2;		//cols are sticks (each two cols represent 1 stick (2 directions))

	if (nSticks)				
	{
		for (i=0; i<SSs.size(); i++)		//SSs
		{
			for (j=0; j<nSticks; j++)		//sticks
			{
				
				short colIndx = j*2;			//col indx in the graph

				

				if (g[i][colIndx].corrW == 999.0)
					continue;		//the assignment is not valid...no need to check the shift validity

				short stickNAA = g[i][colIndx].skeleton.numOfAA();		//number of amino acids in the stick
				int numOfAADiff =  SSs[i].nAA - stickNAA;

				//the shift amount left or right to centerally assign the SS
				numOfAADiff /= 2;

				// The new way of calculating the startIndx and validity...
				g[i][colIndx].startAAindx = SSs[i].startIndx + numOfAADiff;
				
				//make it start from the beginning
				if (g[i][colIndx].startAAindx < 0) 
					g[i][colIndx].startAAindx = 0;

				//make it start accordingly	
				if (g[i][colIndx].startAAindx + stickNAA > nativePDBNumOfAA)
					g[i][colIndx].startAAindx -= (g[i][colIndx].startAAindx + stickNAA) -  nativePDBNumOfAA;

				//the number of AA left of the stick if assigned to this SS
				g[i][colIndx].validityLeft = g[i][colIndx].startAAindx;

				//the number of AA right to the stick if assigned to this SS	
				g[i][colIndx].validityRight = nativePDBNumOfAA - (g[i][colIndx].startAAindx + stickNAA);


				//here the number AA's in reversed stick can differ by 1 or 2 amino acids b/s we directly build the structure reverely 
				//and do not flip the stick, due to that, the number of AA's not neccessary equal to number of AA's in forward direction
				//set the start Indx of the reverse order		
				stickNAA = g[i][colIndx+1].skeleton.numOfAA();		//number of amino acids in the reversed stick
				numOfAADiff =  SSs[i].nAA - stickNAA;
				//the shift amount left or right to centerally assign the SS
				numOfAADiff /= 2;
				// The new way of calculating the startIndx and validity...
				g[i][colIndx+1].startAAindx = SSs[i].startIndx + numOfAADiff;
				//make it start from the beginning
				if (g[i][colIndx+1].startAAindx < 0) 
					g[i][colIndx+1].startAAindx = 0;
				//make it start accordingly	
				if (g[i][colIndx+1].startAAindx + stickNAA > nativePDBNumOfAA)
					g[i][colIndx+1].startAAindx -= (g[i][colIndx+1].startAAindx + stickNAA) -  nativePDBNumOfAA;

				//g[i][colIndx+1].startAAindx = g[i][colIndx].startAAindx;		not valid any more...#AA's not necessary the same

				//the number of AA left of the stick if assigned to this SS
				g[i][colIndx+1].validityLeft = g[i][colIndx+1].startAAindx;

				//the number of AA right to the stick if assigned to this SS	
				g[i][colIndx+1].validityRight = nativePDBNumOfAA - (g[i][colIndx+1].startAAindx + stickNAA);	
				////////////////// end of new the new way of calculation
			}
		}

		//print out the corresponding results
		/*
		for (i=0; i<g.size(); i++)
		{
			cout<<"stick# "<<sticks[i].stickNum<<"  Type "<<sticks[i].ssType<<endl;
			for (j=0; j<g[i].size(); j++)
				cout<<"  with "<<SSs[j/2].type<<" "<<j<<"  startAAindx "<<g[i][j].startAAindx<<"  Validity  Right = "<<g[i][j].validityRight<<"  Left = "<<g[i][j].validityLeft<<" corrW= "<<g[i][j].corrW<<endl;
		}
		*/

	}
	else
	{
		cout<<"========================= in Topology_Graph::setShifts (Protein) =============="<<endl;
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

inline
float getSSScore(int nAASeq, int nAAstick, float parcent)
{
	/*
	//get the score according to the distance not # of AAs
	double seqLength ;
	double denLength = getDistance(stick.start, stick.end);

	if (stick.ssType == 'H')
		//comparing two hlces
		seqLength = ALPHA_RISE * nAASeq;
	else
		//comparing two strands
		seqLength = BETA_RISE * nAASeq;

	if ((seqLength * parcent > denLength) ||		\
		(denLength * parcent > seqLength))
		return 999.00;

	return fabs(seqLength - denLength);
	/////////////////////////////////////////////////////////
	*/

	
	//or get the score according to the difference in # of AAs
	if ((nAASeq * parcent > nAAstick) ||		\
		(nAAstick * parcent > nAASeq))
		return 999.0;

	return fabs(nAASeq - nAAstick);
	////////////////////////////////////////////////////////////
	

	
	
}
//I think we sould consider the maximum possible AA could be considered to the left 
//or to the right of the loop instead of just add MAX_SHIFT_ALLOWED of AA.
//b/s the number of actual AA could be shifted to the left or to the right could 
//be less..... so do something like determining the MAX shift ALLOWED to the left or to the right
inline
void setLinkScore(vector<vector<cell> > & g,				//the graph that will represent the topo. problem
				   vector<SecondaryStruct> SSs,				//list of secondary structures in their original order
				   //vector<short> shortestAssignments,		//shortest assignment for each SS
				   short fromRowIndx,						//the index of stick we are working on...we will deal also with the next stick
				   short fromColIndx,						//the indx of the SS assigned to the current stick
				   short toRowIndx,							//the indx of the next SS
				   short toColIndx)							//the indx of the next stick stick
{

	/*****
			find the virtual length of the loop b/w the two SS 
			(the loop for right now will be considered the maximum loop could be b/w the two SSs,
			in other words, the maximum shift allowed should be considered in this step)
	*****/	

	//get the maximum from the left
	int loopLength = (g[fromRowIndx][fromColIndx].validityLeft  < MAX_SHIFT_ALLOWED) ? (g[fromRowIndx][fromColIndx].validityLeft) : (MAX_SHIFT_ALLOWED);

	//get the maximum from the right
	loopLength += (g[toRowIndx][toColIndx].validityRight < MAX_SHIFT_ALLOWED) ? (g[toRowIndx][toColIndx].validityRight ) : (MAX_SHIFT_ALLOWED);	



	//this is optional step
	//calculate the number of AAs and the length of SSs in between (if any)
	float tSSlengths = 0;
	int nSSAA = 0;
	
	/*
	for (int i=fromRowIndx+1; i<toRowIndx; i++)
	{
		nSSAA += shortestAssignments[i];
		if (SSs[i].type == 'H')
			tSSlengths += shortestAssignments[i] * ALPHA_RISE;
		else
			tSSlengths += shortestAssignments[i] * BETA_RISE;
		//cout<<"nSSAA= "<<nSSAA<<endl;

	}
	*/
	
	loopLength += g[toRowIndx][toColIndx].startAAindx - g[fromRowIndx][fromColIndx].startAAindx - g[fromRowIndx][fromColIndx].skeleton .numOfAA() - nSSAA;

	//forward (from) stick with forward (to) stick
	float mapDistFF = getDistance( g[fromRowIndx][fromColIndx].cTerminal, g[toRowIndx][toColIndx].nTerminal);
	//forward (from) stick with revere (to) stick
	float mapDistFR = getDistance( g[fromRowIndx][fromColIndx].cTerminal, g[toRowIndx][toColIndx+1].nTerminal);
	//reverse (from) stick with forward (to) stick
	float mapDistRF = getDistance( g[fromRowIndx][fromColIndx+1].cTerminal, g[toRowIndx][toColIndx].nTerminal);
	//reverse (from) stick with reverse (to) stick
	float mapDistRR = getDistance( g[fromRowIndx][fromColIndx+1].cTerminal, g[toRowIndx][toColIndx+1].nTerminal);
	
	float loopDist = loopLength*3.8 + tSSlengths;	//the length of the loop with SSs in between in Angstrom

	loopLink tmplink;
	tmplink.rowIndx = toRowIndx;
	tmplink.colIndx = toColIndx;

	//FF
	if (loopDist >= mapDistFF)
	{

		tmplink.w = (int) (loopLength + nSSAA - (mapDistFF/3.8));

		g[fromRowIndx][fromColIndx].links.push_back(tmplink);
	}

	//RF
	if (loopDist >= mapDistRF)
	{
		tmplink.w = (int) (loopLength + nSSAA - (mapDistRF/3.8));

		g[fromRowIndx][fromColIndx+1].links .push_back (tmplink);
	}

	//update col indx
	tmplink.colIndx = toColIndx+1;

	//FR
	if (loopDist >= mapDistFR)
	{
		tmplink.w = (int) (loopLength + nSSAA - (mapDistFR/3.8));

		g[fromRowIndx][fromColIndx].links .push_back (tmplink);
	}

	//RR
	if (loopDist >= mapDistRR)
	{
		tmplink.w = (int) (loopLength + nSSAA - (mapDistRR/3.8));

		g[fromRowIndx][fromColIndx+1].links .push_back (tmplink);
	}


	/*
	//print info
	cout<<"from  strtIndx= "<<g[toRowIndx][toColIndx].startAAindx
		<<" to strtIndx= "<<g[fromRowIndx][fromColIndx].startAAindx
		<<" stk nAA= "<<sticks[fromColIndx/2].nAA<<" nSSAA= "<<nSSAA<<endl;
	cout<<" ["<<fromRowIndx<<"] ["<<fromColIndx<<"] loopNaa= "<<loopLength<<"  loopLength= "<<loopDist<<endl;
	cout<<"    with ["<<toRowIndx<<"] ["<<toColIndx<<"] mapDist= "<<mapDistFF<<" wLink= "<<((int) (loopLength + nSSAA - (mapDistFF/3.8)))<<endl;
	cout<<"         ["<<toRowIndx<<"] ["<<toColIndx+1<<"] mapDist= "<<mapDistFR<<" wLink= "<<((int) (loopLength + nSSAA - (mapDistFR/3.8)))<<endl;
	cout<<" ["<<fromRowIndx<<"] ["<<fromColIndx+1<<"] loopNaa= "<<loopLength<<"  loopLength= "<<loopDist<<endl;
	cout<<"    with ["<<toRowIndx<<"] ["<<toColIndx<<"] mapDist= "<<mapDistRF<<" wLink= "<<(int) (loopLength + nSSAA - (mapDistRF/3.8))<<endl;
	cout<<"         ["<<toRowIndx<<"] ["<<toColIndx+1<<"] mapDist= "<<mapDistRR<<" wLink= "<<(int) (loopLength + nSSAA - (mapDistRR/3.8))<<endl;
	*/

//	Protein twoSS;
//	twoSS.append(g[fromRowIndx][fromColIndx].skeleton, 0, g[fromRowIndx][fromColIndx].skeleton.numOfAA()-1);
//	twoSS.append(g[toRowIndx][toColIndx].skeleton, 0, g[toRowIndx][toColIndx].skeleton.numOfAA()-1);

//	twoSS.writePDB("twoSS.pdb", 1, twoSS.numOfAA());

//	getchar();
//	getchar();


	/*
	//get score in Dist
	if (loopLength * 3.8 > mapDist)
		return fabs(loopLength*3.8 - mapDist);
	*/
		
}

void buildGraph(vector< vector<Coordinate> > edges,		//edges for SS..... should be consistent with sticks datastructure..it means that first edge in edges should be also first stick in sticks
				vector<vector<cell> > & graph,			//graph data structure which will represent correspondance b/w sequence and density map
				vector<SecondaryStruct> SSs,			//sequence SSs
				vector<EdgeStick> sticks,				//density map sticks 
				int proteinNAA,							//number of AAs in the protein....missing AAs are not counted
				short nMissH,							//number of missing map hlces
				short nMissS)							//number of missing map strands
{

	int i,
		j,
		k,
		nSticks = sticks.size(),
		nSS     = SSs.size();

	//SSlengthVariation : This means that the shortest among (stick or SS) should be at least  \
	  (SSlengthVariation) of the length of the longest one
	float SSlengthVariation = 0.4;		

	//rows of the graph represent SSs and Cols represent sticks in two direction (forward and reverse)
	vector< vector<cell> > tmpG (nSS, vector<cell> (nSticks * 2));


	//first initialize the corrW between sticks and SS.... 
	for (i=0; i<nSS; i++)
	{
		for (j=0; j<nSticks; j++)
		{
			short colIndx = 2*j;
			if (sticks[j].ssType == SSs[i].type)
			{
				//set the corresponding weight b/w SS and stick
				// using initial value (forward direction) of #AA's in the stick is not very accurate, b/s forward may differ by 1 or 2 from reverse one
				//but I used them b/s I don't want a case where I correspond to  reverse while not correspond to forward
				tmpG[i][colIndx].corrW = getSSScore(SSs[i].nAA, sticks[j].nAA, SSlengthVariation);
				tmpG[i][colIndx+1].corrW = tmpG[i][colIndx].corrW;

				//save the initial Skeleton (structure) for each valid correspondance
				if (tmpG[i][colIndx].corrW < 999.0)
				{
					//tmpG[i][colIndx].skeleton = initialSkeleton[j];		//the direction is forward
					//tmpG[i][colIndx+1].skeleton = initialSkeleton[j];		//reverse direction
					//flip the structure
					//flipHelix(tmpG[i][colIndx+1].skeleton, sticks[j]);		no need now b/s we directly build it in reverse order


					tmpG[i][colIndx].skeleton = buildStructure(edges[j], sticks[j].ssType, 1);		//send ss type and the direction (forward)

					tmpG[i][colIndx+1].skeleton  = buildStructure(edges[j], sticks[j].ssType, 0);	//send ss type and the direction (reverse)

					//calculate c and n terminals for each stick
					tmpG[i][colIndx].cTerminal = tmpG[i][colIndx].skeleton.getAtomCoordinate(tmpG[i][colIndx].skeleton.numOfAA()-1, " C ");
					tmpG[i][colIndx].nTerminal = tmpG[i][colIndx].skeleton.getAtomCoordinate(0, " N ");

					//for reversed stick
					tmpG[i][colIndx+1].cTerminal = tmpG[i][colIndx+1].skeleton.getAtomCoordinate(tmpG[i][colIndx+1].skeleton.numOfAA()-1, " C ");
					tmpG[i][colIndx+1].nTerminal = tmpG[i][colIndx+1].skeleton.getAtomCoordinate(0, " N ");

				}
			}
		}
	}

	//initiate maximum shift allowed from left and right for each assignment and  \
	  initiate start indeces
	setShifts(SSs, tmpG, proteinNAA);

	//print out the corresponding results
	for (i=0; i<nSS; i++)
	{
		cout<<"SS# "<<i+1<<"  Type "<<SSs[i].type<<" order= "<<SSs[i].order<<" nAA= "<<SSs[i].nAA<<endl;
		for (j=0; j<2*nSticks; j++)
			cout<<"  with "<<sticks[j/2].ssType<<" "<<setw(2)<<j<<" nAA= "<<setw(2)<<tmpG[i][j].skeleton.numOfAA()<<" startAAindx= "<<setw(3)<<tmpG[i][j].startAAindx
				<<"  Valid Right= "<<setw(3)<<tmpG[i][j].validityRight<<"  Left= "<<setw(3)<<tmpG[i][j].validityLeft<<" corrW= "<<tmpG[i][j].corrW<<endl;
	}


	/*
	///////////////////////////////////////////////////////////////////////////////////////////
	//This is an optional step
	//find worst assignment for each SS .... for a given SS, find the shortest stick ..

	vector<short> shortestAssignments;			// # of AAs 
	for (i=0; i<nSS; i++)
	{

		short minLength = 1000;
		short nOfAA = 0;
		//short indx = -1;
		for (j=0; j<nSticks; j++)
			if ((sticks[j].nAA - SSs[i].nAA < minLength) &&		
				(tmpG[i][j*2].corrW != 999.0))				//should be valid assignment also
			{

				minLength = sticks[j].nAA - SSs[i].nAA;
				nOfAA = sticks[j].nAA;
			}

		if(minLength == 1000)
			nOfAA = SSs[i].nAA;

		shortestAssignments.push_back(nOfAA);
	}
	//print shortest Assignments
	cout<<"Shortest Assignments..."<<endl;
	for (i=0; i<nSS; i++)
		cout<<"nAA in shortest Assignments= "<<shortestAssignments[i]<<endl;
	////////////////////////////////////////////////////////////////////////////////////////////
	*/


	short cRow, cCol,			//current row and col
		  nRow, nCol;			//new row and col
	//build links between nodes...
	for (cRow=0; cRow<nSS-1; cRow++)			//current Row is the current SS we are working on
	{
		for (cCol=0; cCol<2*nSticks; cCol += 2)				//current col is the current stick (forward and reverse direction)
		{
			if (tmpG[cRow][cCol].corrW != 999.0)	//if the assignment of this stick with this SS is valid
			{
				short	missH = nMissH,
						missS = nMissS;

				nRow = cRow+1;		//with next Row

				//to consider missing hlces and strands in between
				do
				{
					for (nCol=0; nCol<2*nSticks; nCol += 2)	//next stick for next SSs
					{
						if (tmpG[nRow][nCol].corrW != 999.00)		//if the corresponding Assignment is valid also 
							//don't calculate any link b/w sticks assigned to the same SS
							if (cCol/2 != nCol/2)
								setLinkScore(tmpG, SSs, cRow, cCol, nRow, nCol);
					}
					
					//check the type of next SS to check if it should be skipped
					if (nRow < nSS)
					{
						if (SSs[nRow].type == 'H')
							missH--;
						else
							missS--;
						
					}

					nRow++;
				} while ((missH>=0) && (missS>=0) && (nRow< nSS));
			}
		}
		//getchar();
	}

	//add two rows (start and end)
	vector< vector<cell> > tmpG2 (2, vector<cell> (nSticks * 2));
	tmpG2.insert (tmpG2.begin()+1, tmpG.begin (), tmpG.end());
		
	graph = tmpG2;
	tmpG.clear();
	tmpG2.clear();

	//update links
	for (i=0; i<graph.size(); i++)
		for (j=0; j<graph[i].size(); j++)
			for (k=0; k<graph[i][j].links .size (); k++)
				graph[i][j].links [k].rowIndx += 1;

	//build links from start node
	cRow = 0;
	cCol=0;

	short	missH = nMissH,
			missS = nMissS;
	graph[cRow][cCol].corrW = 0;

	nRow = cRow+1;
	do
	{
		for (nCol=0; nCol < graph[nRow].size (); nCol++)
		{
			if (graph[nRow][nCol].corrW != 999.0)
			{
				loopLink tmplink;
				tmplink.rowIndx = nRow;
				tmplink.colIndx = nCol;
				tmplink.w = 0;
				
				graph[cRow][cCol].links .push_back (tmplink);
			}
		}
		//check the type of next SS to check if it should be skipped
		if (nRow < nSS)
		{
			if (SSs[nRow+1].type == 'H')
				missH--;
			else
				missS--;
			
		}

		nRow++;
	} while ((missH>=0) && (missS>=0) && (nRow< nSS));



	//build links to end node
	cRow = graph.size ()-1;
	cCol = 0;

	missH = nMissH,
	missS = nMissS;
	graph[cRow][cCol].corrW = 0;

	nRow = cRow-1;
	do
	{
		for (nCol=0; nCol < graph[nRow].size (); nCol++)
		{
			if (graph[nRow][nCol].corrW != 999.0)
			{

				loopLink tmplink;
				tmplink.rowIndx = cRow;
				tmplink.colIndx = cCol;
				tmplink.w = 0;
				
				graph[nRow][nCol].links .push_back (tmplink);

			}
		}
		//check the type of next SS to check if it should be skipped
		if (nRow >= 0)
		{
			if (SSs[nRow-1].type == 'H')
				missH--;
			else
				missS--;
			
		}

		nRow--;
	} while ((missH>=0) && (missS>=0) && (nRow >= 0));


	cout<<"Graph is built....Press any key.."<<endl;
	//getchar();
	//print main Information of the graph
	cout<<"///////////////// Graph Info /////////// "<<endl;
	for (i=0; i<graph.size(); i++)
	{
		if ((i>0) && (i<=nSS))
			cout<<i<<" "<<SSs[i-1].type<<endl;
		else
			cout<<"Extra Row..."<<endl;

		for (j=0; j<graph[i].size(); j++)
		{
			cout<<" with "<<sticks[j/2].ssType<<"  corrW= "<<graph[i][j].corrW<<endl;
			for (int k=0; k<graph[i][j].links.size(); k++)
			{
				cout<<"   ["<<graph[i][j].links[k].rowIndx<<", "<<graph[i][j].links[k].colIndx<<"] link weight= "<<graph[i][j].links[k].w<<endl;
			}
		}

		//getchar();
	}
}


void setTopoCandidates(vector<vector<cell> > graph, priority_queue<paths> & topoCandidates)
{
	int kShortest=3000;
	int i,j;

	//graph.erase(graph.begin() + graph.size()-1, 1);
	graph.pop_back ();
	int nCol = graph[0].size(),
		nRow = graph.size();
	//sticks represented by Cols and SS represented by row...foreach stick two cols (forward and reverse directions)
	for(i=0; i< nCol; i++)
	{
		if (graph[0][i].links .size ())
		{
			cout<<"Find paths start from SS indx= "<<i<<endl;
			bool done = false;
			node current;
			current.row = 0;
			current.col = i;
			current.curIndx = 0;

			for (j=0; j<nCol; j++)
			{
				current.allowed.push_back(0);
				current.path .push_back (-1);
			}

			while (!done)
			{
				current = traverse(graph, topoCandidates, current, done);
			}

			//getchar();

			/*
			// keep top kShortest paths //////////////////////
			//////////////////////////////////////////////////
			vector<paths> tmpTopKPaths;

			for (k=0; k<kShortest; k++)
			{
				if (topKPaths.empty())
					break;
				tmpTopKPaths.push_back (topKPaths.top());
				topKPaths.pop();
			}

			while (!topKPaths.empty())
				topKPaths.pop ();

			for (k=0; k<tmpTopKPaths.size (); k++)
				topKPaths.push(tmpTopKPaths[k]);
			//////////////////////////////////////////////////
			*/
		}
	}

	/*
	//print Candidate topologies
	i=1;
	while (!topoCandidates.empty ())
	{
		cout<<i<<"  ";
		for (j=0; j< nCol; j++)
			cout<<topoCandidates.top().col [j]<<" ";
		cout<<"  Total weight = "<<topoCandidates.top ().w <<endl;
		i++;
		topoCandidates.pop ();

		//if (i % 20 == 0)
		//	getchar();
	}
	*/
	
}

/*
vector<vector<float> > corrSS;				//correspondant scores between density sticks and sequence SS
vector<vector<float> > linksScores;			//correspondant scores between links in sequence and links in the density
short	maxMissing		= 0,				//the maximum number of missing SS ... MAX(missing hlces, missing strands)
		nHlcesSeq		= 0,				//number of hlces in the sequence
		nStrandsSeq		= 0,				//number of strands in the sequence
		nHlcesDensity	= 0,				//number of hlces in the density
		nStrandsDensity	= 0;				//number of strands in the density


/////////// DATA STRUCTURE ////////////////
struct node
{
	short SSnum;			//current seqSS num or loop indx
	short level;			//keep track of the level
	float g;				// g(n) function
	float h;				// h(n) function
	float f;				// f(n) function
	priority_queue<pair<short, short> >	order;			//priority queue of the corresponding order b/w seqSS and sticks
	vector<short> loops;				//list of loops for current order.... represents rows in corrLinks matrix
	vector<short> links;				//list of links on the density map for current loops... represents cols in corrLinks matrix

	//initilizer
	node(): SSnum(-1), level(-1), g(0.0), h(0.0), f(0.0) {}
};

// Determine priority.
bool operator<(const node &a, const node &b)
{
  return a.f  < b.f;
}
////
//bool operator<(const short &a, const short &b)
//{
 // return a < b;
//}


//////// PRIORITY QUEUE which will represent the search tree
priority_queue<node> sTree;			//search tree
///////////////////////////////////////////


inline
float getSSScore(int nAASeq, EdgeStick stick, float parcent, char ss)
{
	double seqLength ;
	double denLength = getDistance(stick.start, stick.end);

	//cout<<nAASeq<<" "<<stick.nAA <<"  "<<denLength<<endl;

	if (ss == 'H')
		//comparing two hlces
		seqLength = ALPHA_RISE * nAASeq;
	else
		//comparing two strands
		seqLength = BETA_RISE * nAASeq;

	if ((seqLength * parcent > denLength) ||		\
		(denLength * parcent > seqLength))
		return 9999.00;

	return fabs(seqLength - denLength);
	
}
void buildCorrSS(vector<vector<float> > & corrSS, vector<SecondaryStruct> seqSS, vector<EdgeStick> sticks)
{
	//create and initiate the two dimention array that contains the the corresponding scores b/w 
	// sequence SS and density Sticks
	//rows represent sticks
	//col's represent seq SS

	short nSeq = seqSS.size (),
		  nSticks = sticks.size ();

	vector<vector<float> > tmpCorr (nSticks, vector<float>(nSeq));

	corrSS = tmpCorr;

	short i, j;

	float lengthVariationHlx = 0.5,
		  lengthVariationStrand = 0.4;	

	for (i=0; i<nSticks; i++)
		for (j=0; j<nSeq; j++)
		{
			//strand with hlx or vice versa
			if (seqSS[j].type == sticks[i].ssType)
			{
				if (seqSS[j].type == 'H')
					corrSS[i][j] = getSSScore(seqSS[j].nAA, sticks[i], lengthVariationHlx, 'H');
				else	
					corrSS[i][j] = getSSScore(seqSS[j].nAA, sticks[i], lengthVariationStrand, 'S');	
			}
		}
}

//I think we sould consider the maximum possible AA could be considered to the left 
//or to the right of the loop instead of just add MAx_SHIFT_ALLOWED of AA.
//b/s the number of actual AA could loop be shifted to the left or to the right could 
//be less..... so do something like determining the MAX shift ALLOWED to the left or to the right
inline
float getLinkScore(short seqLength, double dist)
{

	
	//short densityLength = dist/3.8;
	float seqDist = (seqLength + MAX_SHIFT_ALLOWED + MAX_SHIFT_ALLOWED) * 3.8;
	//cout<<setprecision(5)<<seqLength<<"  "<<seqDist<<" "<<dist<<endl;
	if ( seqDist >= dist)
		return fabs(seqDist - dist);

	return 9999.00;
}

void buildLinksScores(vector<vector<float> > & linksScores, vector<SecondaryStruct> seqSS, vector<EdgeStick> sticks, vector<short> & loops)
{
	short i, j, k, m, nSeq = seqSS.size (), nSticks = sticks.size ();

	//get the number of hlces and strands in the sequence
	for (i=0; i<nSeq; i++)
	{
		if (seqSS[i].type == 'H')
			nHlcesSeq++;
		else
			nStrandsSeq++;
	}

	//get the number of hlces sticks and strand sticks in the density
	for (i=0; i<nSticks; i++)
	{
		if (sticks[i].ssType == 'H')
			nHlcesDensity++;
		else
			nStrandsDensity++;
	}

	//get the maximum number of missing SS b/w sequence and sticks
	if (nHlcesSeq - nHlcesDensity >= nStrandsSeq - nStrandsDensity)
		maxMissing = nHlcesSeq - nHlcesDensity;
	else
		maxMissing = nStrandsSeq - nStrandsDensity;

	//calculate the number of links could be exist in the sequence
	short nSeqLinks = nSeq - 1;
	for (i=1; i<=maxMissing; i++)
		nSeqLinks += nSeq - i - 1;

	cout<<"number of links in the sequence is "<<nSeqLinks<<endl;

	short nDensityLinks = 0;

	//for (i=1; i<nSticks; i++)
	//	nDensityLinks += (nSticks-i) * 4;

	nDensityLinks= nSticks * (nSticks-1) * 2;

	cout<<"number of links in the density is "<<nDensityLinks<<endl;

	vector< vector<float> > tmpLinks (nSeqLinks, vector<float>(nDensityLinks));
	linksScores = tmpLinks;
	tmpLinks.clear ();

	int rowCntr=0, colCntr;
	for (i=0; i<=maxMissing; i++)
	{
		for (j=0; j<nSeq-i-1; j++)
		{
			short nAA = seqSS[j+i+1].startAAnum - seqSS[j].endAAnum - 1;

			loops.push_back (nAA);
	
			colCntr = 0;
			for (k=0; k<nSticks-1; k++)
			{
				for (m=k+1; m<nSticks; m++)
				{
					
					double dist = getDistance(sticks[k].start, sticks[m].start);

					linksScores[rowCntr][colCntr] = getLinkScore(nAA, dist);
					colCntr++;
					dist = getDistance(sticks[k].start, sticks[m].end);
					linksScores[rowCntr][colCntr] = getLinkScore(nAA, dist);
					colCntr++;
					dist = getDistance(sticks[k].end, sticks[m].start);
					linksScores[rowCntr][colCntr] = getLinkScore(nAA, dist);
					colCntr++;
					dist = getDistance(sticks[k].end, sticks[m].end);
					linksScores[rowCntr][colCntr] = getLinkScore(nAA, dist);
					colCntr++;
				}
			}
			rowCntr++;
		}
	}
}

short getSSLinkIndx(short SS1num, short SS2num)
{
	int i;
	int offset=0;
	for (i=0; i<SS2num - SS1num-1; i++)
		offset = offset + nHlcesSeq + nStrandsSeq - 1 - i;
	return offset + SS1num-1;
}

short getDensityLinkIndx(short stick1indx, short stick2indx)
{
	short nTotal = nHlcesDensity + nStrandsDensity - 1;
	short offset = 0, i;
	if (stick1indx < stick2indx)
	{
		for (i=0; i<stick1indx; i++)
			offset += nTotal- i;
		return offset + stick2indx - i - 1;
	}
	else
	{
		for (i=0; i<stick2indx; i++)
			offset += nTotal - i;
		return offset + stick1indx - i - 1;
	}


}
*/

#endif
