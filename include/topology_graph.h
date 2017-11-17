#ifndef	TOPOLOGY_GRAPH_H
#define TOPOLOGY_GRAPH_H

#include <vector>
#include <bitset>
#include <queue>

//#define _WIN

#ifdef _WIN
#include <time.h>					//used to find elapsed time
#else
#include <sys/time.h>
#endif

#include "skeleton_overall.h"
#include "rmsd.h"

#include "MRC.h"

using namespace std;


#define MAX_NUM_STICKS	32				//maximum number of sticks...used by the recurrence


/////////// DATA STRUCTURE /////////////////////////////////////
struct node					//this data structure to store the information of current node we are working on the graph
{
	short row;									//represent the index of the SS
	short col;									//represent the index of the stick
	short curLink;								//the indx of the link we are working on

	vector<bool> visited;						//to check that this stick has been assigned to any SS before
	vector<short> topology;						//The path (list of SS chosen) ... represent the topology

	float w;									//toTal weight from start node to this node

	short nStickSolved;							//number of sticks solved (assigned) so far...used to terminate

	//Initializer
	node() : row(-1), col(-1), curLink(-1), w(0), nStickSolved(0) {}
};


struct loopLink			//represent the correspondance between the loop and the link b/w two sticks
{
	short rowIndx;		//row indx of the sink node
	short colIndx;		//col indx of the sink node
	float w;			//weight of the link

	//Initializer
	loopLink () : rowIndx(-1), colIndx(-1), w(999.9) {}
};

struct cell
{
	node	parent;			//you can instead...use 3 variables...col, row and parent indx
							//if no parent then row = -1 and col = -1 parent Indx = -1

	/*
	 *	Information of coresponding generated structure
	 */
	Protein skeleton;					//the structure of the SS corresponding to this stick
	int validityRight;					//the maximum size of shift could be applied to the sequence to the right...
	int validityLeft;					//the maximum size of shift could be applied to the sequence to left...
	int startAAindx;					//the start indx of AAs to be assigned to this stick
	Coordinate	cTerminal,				//c terminal of the initial structure
				nTerminal;				//n terminal of the initial structure

	double nPaths;						//estimated number of paths from this node to end node


	float	corrW;						//corresponding weight of assigning this SS to this stick


	vector<loopLink>	outLinks,		//outgoing links (edges)
						inLinks;		//ingoing links  (edges)
	vector<short>		nAAloop;		//number of AAs in the loop corresponding to each outLink
	vector<string>		traceNum;		//the best trace associated with the loop (trace from density map);

	short curLink;						//the link (edge) that we are working on currently


	//Initializer
	cell() : validityRight(0), validityLeft(0), startAAindx (-1), nPaths(0), corrW(999.0), curLink(-1) {}
};

struct paths
{
	vector<short> ssAssigned;				//list of SS assigned to sticks
	float w;								//total weight from start to end ... the weight of such a topology

	paths(): w(0) {}
};

// Determine priority.
bool operator<(const paths &a, const paths &b)
{
 return a.w   > b.w ;
}

struct cellSets{

	char	**nxtCol;				//the index of the col of the node next to this node in the path...nxtCol[0] for forward and nxtCol[1] for backward col
	char	**nxtRow;				//the indx of the row of the node loacated after this node in the path
	unsigned int	**nxtSetIndx;	//the indx of the set in the next node. equal the set in the next node plus my col
	float	**weight;				//weight[0][i] is the weight for forward node...weight[1][i] is the weight for backward node
	char	*length;				//the length of the path...it is a parallel data structure with sets

	unsigned int *sets;				//sets represented as bit arrays
};

// data structure to store sequence of nodes (by rowIndx and colIndx) represents a path with its weight
//used to generate K-shortest paths
struct iPath {
	float w;
	vector<pair<short, short> >  pCandidate;
	//short *pCandidateRow;
	//short *pCandidateCol;
	int parentPath;										//the index of parent path
	int deviationNode;									//the node this path deviate from the parent path

	iPath(): w(9999.0), parentPath(-1), deviationNode(-1) {}
};
double nPaths = 0;

/////////// DATA STRUCTURE ///////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
//////////////////////////// LIST OF FUNCTIONS ///////////////////////
//////////////////////////////////////////////////////////////////////
inline void getChild(vector<vector<cell> > & g, node & current, node & child);
inline node getParent(vector<vector<cell> >  & g, node & current);
inline node traverse(vector<vector<cell> > & g, priority_queue<paths> & topK, node & current, bool & done);
void writeGraph(vector<vector<cell> > &graph, string gOUT, bool wIN=false);
void printTopos(priority_queue<paths> topK);
void writeTopology(vector<vector<cell> > &graph, Protein nativePDB, paths topology, string outFile);				//given a topology "path" print out a sample structure for it
void setShifts(vector<SecondaryStruct> SSs, vector<vector<cell> > & g, int nativePDBNumOfAA);
inline float getSSScore(int nAASeq, int nAAstick, float parcent);
inline void setLinkScore(vector<vector<cell> > & g,			//the graph that will represent the topo. problem
				   vector<SecondaryStruct> SSs,				//list of secondary structures in their original order
				   //vector<short> shortestAssignments,		//shortest assignment for each SS
				   short fromRowIndx,						//the index of stick we are working on...we will deal also with the next stick
				   short fromColIndx,						//the indx of the SS assigned to the current stick
				   short toRowIndx,	 						//the indx of the next SS
				   short toColIndx);						//the indx of the next stick stick
inline void cleanUp(vector<vector<cell> > & graph);
void buildGraph(vector< vector<Coordinate> > edges,		//edges for SS..... should be consistent with sticks datastructure..it means that first edge in edges should be also first stick in sticks
				vector<vector<cell> > & graph,			//graph data structure which will represent correspondance b/w sequence and density map
				vector<SecondaryStruct> SSs,			//sequence SSs
				vector<EdgeStick> sticks,				//density map sticks
				int proteinNAA,							//number of AAs in the protein....missing AAs are not counted
				short nMissH,							//number of missing map hlces
				short nMissS);							//number of missing map strands
void setTopoCandidates(vector<vector<cell> > graph, priority_queue<paths> & topoCandidates);
iPath getShortestPath(cellSets** graphSets,
				vector<vector<cell> > & graph,
				long** cntr,
				short sRow,							//the node to find the shortest path from (given the indx of row and col)
				short sCol,
				short &nRslvd,
				//vector<pair<short, short> > skipLinks,
				short *skipLinksRow,
				short *skipLinksCol,
				short skipLinksSize,
				int visitedCols);
void gShortestRec(vector<vector<cell> > graph, short ***topoLis, int K);	//find shortest path using a recurrence
inline
void getShortestTraces(short *nodes, vector<vector<short> >	&adjMtrx, float **linksW, float *maxLength, int source, short *(&prev), float *(&D));

double getRMSDSSE(vector<Protein> &portions, Protein &nativePDBFile, int *sIndeces);	//get the RMSD b/w two SSEs
Protein buildFullModel(Protein &pdb1, vector<vector<cell> > graph, vector<vector<Coordinate> > hEdges, short **topology, Coordinate **traceList, float *traceLength, short *traceNoOfPoints, short *shifts, short &sIndx, double *rmsd);
//////////////////////////////////////////////////////////////////////

inline
void getChild(vector<vector<cell> > & g, node & current, node & child)
{

	//I will not check if it has children or not...because it should have if this has been called
	child = current;


	child.col = g[current.row ][current.col ].outLinks [current.curLink].colIndx;

	child.row = g[current.row ][current.col ].outLinks [current.curLink].rowIndx;


	child.w += g[current.row ][current.col].outLinks [current.curLink].w;

	//add corresponding weight
	child.w += g[current.row][current.col].corrW;


	child.visited [child.col] = 1;
	child.topology[child.col] = child.row;

	//disable the twin SS
	if (child.col % 2 == 0)
	{
		child.visited [child.col + 1] = 1;
		child.topology[child.col+1] = 0;
	}
	else
	{
		child.visited [child.col - 1] = 1;
		child.topology [child.col - 1] = 0;
	}


	child.curLink = 0;
	child.nStickSolved += 1;

	g[child.row][child.col].parent = current;
}

inline
node getParent(vector<vector<cell> >  & g, node & current)
{
	node parent;

	parent = g[current.row][current.col].parent;
	parent.curLink += 1;

	return parent;
}

inline
node traverse(vector<vector<cell> > & g, priority_queue<paths> & topK, node & current, bool & done)
{
	node child;

	//cout<<"cur.row= "<<current.row<<"  cur.Col= "<<current.col<<" cur.indx= "<<current.curIndx<<"  cur.w= "<<current.w<<endl;
	//getchar();
	if (current.curLink >= g[current.row][current.col].outLinks .size())
	{
		//cout<<"no more links.."<<endl;
		//this happens when no node remains in this subtree and need to look to the parent
		if (g[current.row][current.col].parent.row != -1)		//it is not the root
		{

			//cout<<"get parent..no more links"<<endl;
			return getParent(g, current);

			//return parent;
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
		if (g[current.row][current.col].outLinks [current.curLink].rowIndx == g.size ())
		{
			//cout<<"Goes to end.. "<<g.size ()<<endl;
			//if the link goes to end node...then check the number of sticks assigned so far \
			  if all sticks assigned you are done...other wise select sibling of current node
			if (current.nStickSolved > current.visited .size () / 2 - 1)
			{
				//cout<<"    Done...found a path.."<<current.path .size()<<" path : ";

				/////////////////////////////////////////////////////////////////////////////
				//	This change is made just to know the number of Paths without saving them....I need to save paths If I wanna
				//	work on them later....to do so...just uncomment next statements
				///////////////////////////////////////////////////////////////////////////////
				paths tmp;

				tmp.ssAssigned = current.topology;
				tmp.w = current.w;

				tmp.w += g[current.row ][current.col].outLinks [current.curLink ].w;
				tmp.w += g[current.row ][current.col].corrW;

				tmp.ssAssigned [current.col] = current.row;

				if (current.col % 2 == 0)
					tmp.ssAssigned [current.col + 1] = 0;
				else
					tmp.ssAssigned [current.col-1] = 0;

				topK.push(tmp);

				//////////////////////////////////////////////////////////////////////////////


				current.curLink ++;

				return current;
			}
			else
			{
				//cout<<"not Done yet.."<<endl;
				//get sibling
				if(current.curLink < g[current.row][current.col].outLinks.size() - 1)
				{
					//cout<<"get sibling.. not Done yet"<<endl;
					current.curLink += 1;
					return current;
				}
				else
				{
					//cout<<"get parent...not done yet.."<<endl;
					return getParent(g, current);
					//return parent;
				}
			}
		}
		else
		{

			//check if the link goes to valid assignment
			if (current.visited [g[current.row][current.col].outLinks [current.curLink].colIndx] == 0)
			{
				//cout<<"getchild...."<<endl;
				getChild(g, current, child);
				return child;
			}
			else
			{
				//cout<<"the assignment to go to is not valid.."<<endl;
				//get sibling
				if(current.curLink < g[current.row][current.col].outLinks.size() - 1)
				{
					//cout<<"get sibling..NV Assignment"<<endl;
					current.curLink += 1;
					return current;
				}
				else
				{
					//cout<<"get parent... NV assignment.."<<endl;
					return getParent(g, current);
					//return parent;
				}
			}
		}
	}

	//getchar();
}
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
inline
bool getChild1(vector<vector<cell> > & g, node & current, node & child)
{

	//search for the first valid child
	while (current.curLink  < g[current.row ][current.col ].outLinks .size ())
	{
		child.col = g[current.row ][current.col ].outLinks [current.curLink].colIndx;
		child.row = g[current.row ][current.col ].outLinks [current.curLink].rowIndx;

		//check if it is not refering to end node
		if (child.row < g.size()-1)
		{

			if (current.visited [child.col] == 0)		//has not been assigned before
			{
				child = current;
				child.col = g[current.row ][current.col ].outLinks [current.curLink].colIndx;
				child.row = g[current.row ][current.col ].outLinks [current.curLink].rowIndx;
				child.w += g[current.row ][current.col].outLinks [current.curLink ].w;
				//child.w += g[current.row][current.col].corrW;        //uncomment if u want to add node weight
				child.visited [child.col] = 1;
				child.topology[child.col] = child.row;

				//disable the twin SS
				if (child.col % 2 == 0)
				{
					child.visited [child.col + 1] = 1;
					child.topology [child.col+1] = 0;
				}
				else
				{
					child.visited [child.col - 1] = 1;
					child.topology [child.col - 1] = 0;
				}

				child.curLink = 0;
				child.nStickSolved += 1;

				g[child.row][child.col].parent = current;

				return true;
			}
		}
		else
			return false;

		current.curLink++;
	}

	return false;

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///// different method of traversing....in this method you should keep the end node
/// qRow : quit row is the level you want to stop in the case you come back to
inline
node traverse1(vector<vector<cell> > & g, priority_queue<paths> & topK, node & current, short qRow, bool & done, double & nPaths)
{
	node parent, child;

	if (getChild1(g, current, child))
		return child;						//still has a valid child
	else
	{
		if (child.row+1 >= g.size())			//check if it reaches the end node (refers to end node)
		{
			if (current.nStickSolved > current.visited .size () / 2 - 1)		//all sticks have been assigned?
			{
				//save the topology

				paths tmp;

				tmp.ssAssigned = current.topology;
				tmp.w = current.w;
				//tmp.w += g[current.row ][current.col ].corrW;		//uncomment if u want to add node weight

				topK.push(tmp);


				nPaths++;
			}

			//cout<<"checking current.row = "<<current.row<<endl;
			//if (current .row != qRow)
			//	return getParent(g, current);
			//else
			//{
			//	done = true;
			//	return current;
			//}

		}
		//else
		//{
			if (g[current.row ][current.col ].parent .row == qRow )			//start node...we are done
			{
				done = true;
				return current;
			}
			else
				return getParent(g, current);
		//}
	}

}
//////////////////////////////////////////////////////////////////////////////////////////
//wIN is a boolean flag true if you want to print the graph with inlinks
void writeGraph(vector<vector<cell> > &graph, string gOUT, bool wIN)
{
	ofstream out;
	out.open(gOUT.c_str());
	if (!out) {
		errMsg("TplgyGraph", "writeGraph", "Unable to open " + gOUT, true);
	}

	int i, j, k;
	//print main Information of the graph

	out<<setw((graph[0].size() * 30)/2)<<" "<<"////////////////////////////////////////////////////////////////////////////// "<<endl;
	out<<setw((graph[0].size() * 30)/2)<<" "<<"////////////////////////////////////// Graph Info //////////////////////////// "<<endl;
	out<<setw((graph[0].size() * 30)/2)<<" "<<"////////////////////////////////////////////////////////////////////////////// "<<endl<<endl;

	//write start node information
	out<<setw((graph[0].size() * 34)/2)<<"        "<<"START"<<endl;
	out<<setw((graph[0].size() * 34)/2)<<"        "<<"=OUT="<<endl;

	for (j=0; j<graph[0].size (); j++)
	{
		for (k=0; k<graph[0][j].outLinks .size (); k++)
			out<<setw((graph[0].size() * 34)/2-4)<<" "<<"(nAAloop)(trace#)["<<graph[0][j].outLinks[k].rowIndx<<", "<<graph[0][j].outLinks[k].colIndx<<"] lnkW= "<<graph[0][j].outLinks[k].w<<endl;
	}
	out<<endl<<endl;

	short prcsion = 3;
	short adjustment = 0;

	if (!wIN)
		adjustment = 5;
	for (i=1; i<graph.size()-1; i++)
	{
		out<<right;
		for (j=0; j<graph[i].size (); j++){
			if (wIN)
				out<<setw(10)<<" "<<"row_"<<setw(2)<<i<<" "<<"Col_"<<setw(2)<<j<<setw(24)<<" ";
			else
				out<<setw(6)<<" "<<"row_"<<setw(2)<<i<<" "<<"Col_"<<setw(2)<<j<<setw(14)<<" ";
		}
		out<<endl;

		for (j=0; j<graph[i].size (); j++){
			if (wIN)
				out<<setw(12)<<" "<<"corrW= "<<setw(prcsion)<<graph[i][j].corrW<<setw(25)<<" ";
			else
				out<<setw(8)<<" "<<"corrW= "<<setw(prcsion)<<graph[i][j].corrW<<setw(15)<<" ";
		}
		out<<endl;
		for (j=0; j<graph[i].size (); j++){
			if (wIN)
				out<<setw(12)<<" "<<"sIndx= "<<setw(prcsion)<<graph[i][j].startAAindx<<setw(25)<<" ";
			else
				out<<setw(8)<<" "<<"sIndx= "<<setw(prcsion)<<graph[i][j].startAAindx<<setw(15)<<" ";
		}
		out<<endl;
		for (j=0; j<graph[i].size (); j++){
			if (wIN)
				out<<setw(3)<<" "<<"=IN="<<setw(17)<<" "<<"=OUT="<<setw(18)<<" ";
			else
				out<<"lpAA    Trace      =OUT=         ";
		}
		out<<endl;

		out<<setprecision(prcsion);
		//find the max iterator (the col with max number of inlinks and outlinks)
		short maxItOUT = 0, maxItIN=0;
		for (j=0; j<graph[i].size (); j++){
			if (graph[i][j].outLinks .size ()>maxItOUT)
				maxItOUT = graph[i][j].outLinks .size ();
			if (graph[i][j].inLinks .size ()>maxItIN)
				maxItIN = graph[i][j].inLinks .size ();
		}
		out<<left;
		if (wIN){
			for (k=0; k<MAX(maxItIN, maxItOUT); k++){			//assume the maximum possible number of in or out links

				for (j=0; j<graph[i].size (); j++)
				{
					//int cntr = MAX(graph[i][j].outLinks .size (), graph[i][j].inLinks .size ());

					if (k<graph[i][j].inLinks .size ())
						out<<"["<<setw(2)<<graph[i][j].inLinks [k].rowIndx <<","<<setw(2)<<graph[i][j].inLinks [k].colIndx <<"]="<<setw(6)<<setfill(' ')<<floor(graph[i][j].inLinks [k].w * 100)/100;
					else
						out<<setw(14)<<" ";

					if (k<graph[i][j].outLinks .size ())
						out<<"("<<setw(3)<<graph[i][j].nAAloop [k]<<")("<<setw(11)<<graph[i][j].traceNum[k]<<")["<<setw(2)<<graph[i][j].outLinks[k].rowIndx<<","<<setw(2)<<graph[i][j].outLinks[k].colIndx<<"]="<<setw(6)<<setfill(' ')<<floor(graph[i][j].outLinks[k].w * 100)/100<<" ";
					else
						out<<setw(33)<<" ";

				}
				out<<endl;
			}
		}
		else{
			for (k=0; k<maxItOUT; k++){
				for (j=0; j<graph[i].size (); j++){
					if (k<graph[i][j].outLinks .size ())
						out<<"("<<setw(3)<<graph[i][j].nAAloop [k]<<")("<<setw(11)<<graph[i][j].traceNum[k]<<")["<<setw(2)<<graph[i][j].outLinks[k].rowIndx<<","<<setw(2)<<graph[i][j].outLinks[k].colIndx<<"]="<<setw(6)<<setfill(' ')<<floor(graph[i][j].outLinks[k].w * 100)/100<<" ";
					else
						out<<setw(33)<<" ";
				}
				out<<endl;
			}
		}

		out<<endl;
	}
	//write End node information
	out<<endl<<setw((graph[0].size() * 34)/2)<<" "<<"END"<<endl;
	out<<setw((graph[0].size() * 34)/2)<<" "<<"=IN="<<endl;

	for (j=0; j<graph[graph.size ()-1].size (); j++)
	{
		for (k=0; k<graph[graph.size ()-1][j].inLinks.size (); k++)
			out<<setw((graph[0].size() * 34)/2 - 4)<<" "<<"["<<graph[graph.size ()-1][j].inLinks[k].rowIndx<<", "<<graph[graph.size ()-1][j].inLinks[k].colIndx<<"] lnkW= "<<graph[graph.size ()-1][j].inLinks[k].w<<endl;
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
void printTopos(priority_queue<paths> topoCandidates)
{
	int j;
	while (!topoCandidates.empty ())
	{
		for (j=0; j< topoCandidates.top().ssAssigned .size(); j++)
			cout<<topoCandidates.top ().ssAssigned [j]<<" ";
		cout<<"  Total weight = "<<topoCandidates.top ().w <<endl;
		topoCandidates.pop ();
	}
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void writeTopology(vector<vector<cell> > &graph, Protein nativePDB, paths topology, string outFile){
	vector<Protein> ssStructures;			//list of SS structures
	vector<EdgeStick> sticks;
	EdgeStick tmpStick;

	int i=0, j=0;
	while (i<topology.ssAssigned .size ()/2){

		//save sticks according to sequence order
		j=0;
		short min=9999.0;
		short minIndx=-1;
		//find next seq. SS
		while (j<topology.ssAssigned .size ()){
			if (topology.ssAssigned [j] != 0 && topology.ssAssigned [j] < min){
				min = topology.ssAssigned [j];
				minIndx = j;
			}
			j++;
		}

		if (minIndx != -1){
			topology.ssAssigned [minIndx] = 9999.0;

			//save information of the SS
			tmpStick.ssAssigned = min-1;
			tmpStick.startAAIndx = graph[min][minIndx].startAAindx;
			tmpStick.direction = minIndx % 2;
			tmpStick.validityRight = graph[min][minIndx].validityRight;
			tmpStick.validityLeft = graph[min][minIndx].validityLeft;
			tmpStick.nAA = graph[min][minIndx].skeleton.numOfAA();

			//set shift randomly
			short minShift = MIN(tmpStick.validityRight, tmpStick.validityLeft);
			tmpStick.shift = getRandom(-MIN(2,minShift), MIN(2,minShift));

			sticks.push_back (tmpStick);
			ssStructures.push_back (graph[min][minIndx].skeleton);
		}
		i++;
	}

	Protein fStructure;

	//write all structures into one file
	for (i=0; i<ssStructures.size (); i++){
		AssignSeqToSkeleton(nativePDB, ssStructures[i], sticks[i].shift , i, sticks);
		fStructure.append(ssStructures[i], 0, ssStructures[i].numOfAA()-1);
	}

	fStructure.writePDB(outFile, 1, fStructure.numOfAA());
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
//////////////////////////////////////////////////////////////////////////////////////////////////
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
	if ((nAASeq * parcent >= nAAstick) ||		\
		(nAAstick * parcent >= nAASeq))
		return 999.0;

	return fabs(nAASeq - nAAstick);
	////////////////////////////////////////////////////////////




}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
	short loopLength = (g[fromRowIndx][fromColIndx].validityLeft  < MAX_SHIFT_ALLOWED) ? (g[fromRowIndx][fromColIndx].validityLeft) : (MAX_SHIFT_ALLOWED);

	//get the maximum from the right
	loopLength += (g[toRowIndx][toColIndx].validityRight < MAX_SHIFT_ALLOWED) ? (g[toRowIndx][toColIndx].validityRight ) : (MAX_SHIFT_ALLOWED);

	//updated 3/27/2011
	// the old way is modified
	//loopLength += g[toRowIndx][toColIndx].startAAindx - g[fromRowIndx][fromColIndx].startAAindx - g[fromRowIndx][fromColIndx].skeleton .numOfAA();

	//The new way of calculating the number of AA b/w two SS segments is now simply the number of AA in b/w on sequence plus max Shift allowed...
	//and I add 2 AA for prediction error on the skeleton

	loopLength += SSs[toRowIndx].startIndx  - SSs[fromRowIndx].startIndx - SSs[fromRowIndx].nAA;

	//add 2 AA to tolerate prediction error on skeleton for short loops
	if (loopLength < 3)
		loopLength += 2;
	if (loopLength < 7)
		loopLength += 2;



	//forward (from) stick with forward (to) stick
	float mapDistFF = getDistance( g[fromRowIndx][fromColIndx].cTerminal, g[toRowIndx][toColIndx].nTerminal);
	//forward (from) stick with revere (to) stick
	float mapDistFR = getDistance( g[fromRowIndx][fromColIndx].cTerminal, g[toRowIndx][toColIndx+1].nTerminal);
	//reverse (from) stick with forward (to) stick
	float mapDistRF = getDistance( g[fromRowIndx][fromColIndx+1].cTerminal, g[toRowIndx][toColIndx].nTerminal);
	//reverse (from) stick with reverse (to) stick
	float mapDistRR = getDistance( g[fromRowIndx][fromColIndx+1].cTerminal, g[toRowIndx][toColIndx+1].nTerminal);

	float loopDist = loopLength*3.8;	//the length of the loop with SSs in between in Angstrom



//	cout<<"loopLngth= "<<loopLength<<"mapFF= "<<mapDistFF<<" mapFR= "<<mapDistFR<<" mapRF= "<<mapDistRF<<" mapRR= "<<mapDistRR<<endl;


	loopLink tmpOutLink,
			 tmpInLink;
	tmpOutLink.rowIndx = toRowIndx;
	tmpOutLink.colIndx = toColIndx;


	//FF
	if (loopDist >= mapDistFF)
	{

		tmpInLink.w  = tmpOutLink.w = (int) (loopLength - (mapDistFF/3.8));

		g[fromRowIndx][fromColIndx].outLinks.push_back(tmpOutLink);
		g[fromRowIndx][fromColIndx].nAAloop .push_back (loopLength);
		g[fromRowIndx][fromColIndx].traceNum .push_back ("NONE");				//initiate trace number


		tmpInLink.rowIndx = fromRowIndx;
		tmpInLink.colIndx = fromColIndx;
		g[toRowIndx][toColIndx].inLinks .push_back (tmpInLink);
	}

	//RF
	if (loopDist >= mapDistRF)
	{
		tmpInLink.w  = tmpOutLink.w = (int) (loopLength - (mapDistRF/3.8));

		g[fromRowIndx][fromColIndx+1].outLinks .push_back (tmpOutLink);
		g[fromRowIndx][fromColIndx+1].nAAloop .push_back (loopLength);
		g[fromRowIndx][fromColIndx+1].traceNum .push_back ("NONE");

		tmpInLink.rowIndx = fromRowIndx;
		tmpInLink.colIndx = fromColIndx+1;
		g[toRowIndx][toColIndx].inLinks .push_back (tmpInLink);

	}

	//update col indx
	tmpOutLink.colIndx = toColIndx+1;

	//FR
	if (loopDist >= mapDistFR)
	{
		tmpInLink.w = tmpOutLink.w = (int) (loopLength - (mapDistFR/3.8));

		g[fromRowIndx][fromColIndx].outLinks .push_back (tmpOutLink);
		g[fromRowIndx][fromColIndx].nAAloop.push_back  (loopLength);
		g[fromRowIndx][fromColIndx].traceNum .push_back ("NONE");

		tmpInLink.rowIndx = fromRowIndx;
		tmpInLink.colIndx = fromColIndx;
		g[toRowIndx][toColIndx+1].inLinks .push_back (tmpInLink);

	}

	//RR
	if (loopDist >= mapDistRR)
	{
		tmpInLink.w = tmpOutLink.w = (int) (loopLength - (mapDistRR/3.8));

		g[fromRowIndx][fromColIndx+1].outLinks .push_back (tmpOutLink);
		g[fromRowIndx][fromColIndx+1].nAAloop .push_back (loopLength);
		g[fromRowIndx][fromColIndx+1].traceNum .push_back ("NONE");

		tmpInLink.rowIndx = fromRowIndx;
		tmpInLink.colIndx = fromColIndx+1;
		g[toRowIndx][toColIndx+1].inLinks .push_back (tmpInLink);

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
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
inline
void cleanUp(vector<vector<cell> > & graph)
{
	int i, j, k, m;

	for (i=1; i<graph.size()-1; i++)		//skip start and end node
	{
		for (j=0; j<graph[i].size(); j++)
		{
			if (graph[i][j].outLinks .size () == 0)		//no outgoing links
			{
				//clean up
				for (k=0; k<graph[i][j].inLinks.size(); k++)
				{
					short	inRow = graph[i][j].inLinks[k].rowIndx,
							inCol = graph[i][j].inLinks[k].colIndx;

					for (m=0; m<graph[inRow][inCol].outLinks.size(); m++)
					{
						if ((graph[inRow][inCol].outLinks[m].rowIndx == i) &&
							(graph[inRow][inCol].outLinks[m].colIndx == j))
						{
							graph[inRow][inCol].outLinks.erase(graph[inRow][inCol].outLinks.begin() + m);
							graph[inRow][inCol].nAAloop.erase(graph[inRow][inCol].nAAloop.begin() + m);
							break;
						}
					}
				}

				graph[i][j].inLinks .clear();

			}

			//this is not that much important.....if no inlinks then delete all outlinks
			if (graph[i][j].inLinks .size () == 0)
			{
				for (k=0; k<graph[i][j].outLinks .size (); k++)
				{
					short outRow = graph[i][j].outLinks [k].rowIndx,
							outCol = graph[i][j].outLinks [k].colIndx;
					for (m=0; m<graph[outRow][outCol].inLinks .size (); m++)
					{
						if ((graph[outRow][outCol].inLinks [m].rowIndx == i) &&
							(graph[outRow][outCol].inLinks [m].colIndx == j))
						{
							graph[outRow][outCol].inLinks .erase (graph[outRow][outCol].inLinks .begin () + m);
							break;
						}
					}
				}

				graph[i][j].outLinks .clear();

			}

		}
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
	
	/*
	//print out the corresponding results
	for (i=0; i<nSS; i++)
	{
		cout<<"SS# "<<i+1<<"  Type "<<SSs[i].type<<" order= "<<SSs[i].order<<" nAA= "<<SSs[i].nAA<<endl;
		for (j=0; j<2*nSticks; j++)
			cout<<"  with "<<sticks[j/2].ssType<<" "<<setw(2)<<j<<" nAA= "<<setw(2)<<tmpG[i][j].skeleton.numOfAA()<<" startAAindx= "<<setw(3)<<tmpG[i][j].startAAindx
				<<"  Valid Right= "<<setw(3)<<tmpG[i][j].validityRight<<"  Left= "<<setw(3)<<tmpG[i][j].validityLeft<<" corrW= "<<tmpG[i][j].corrW<<endl;
	}
	*/


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
							if (cCol/2 != nCol/2){
								//cout<<endl<<"["<<cRow+1<<" "<<cCol<<"]  and ["<<nRow+1<<" "<<nCol<<"]"<<endl;
								setLinkScore(tmpG, SSs, cRow, cCol, nRow, nCol);
								//if (cRow==6 && cCol==18) getchar();
							}
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
		{
			for (k=0; k<graph[i][j].outLinks .size (); k++)
				graph[i][j].outLinks[k].rowIndx += 1;
			for (k=0; k<graph[i][j].inLinks .size (); k++)
				graph[i][j].inLinks [k].rowIndx += 1;
		}

	//build links from start node
	cRow = 0;
	cCol=0;

	short	missH = nMissH,
			missS = nMissS;
	graph[cRow][cCol].corrW = 0;

	nRow = cRow+1;		//new row
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

				graph[cRow][cCol].outLinks .push_back (tmplink);
				graph[cRow][cCol].nAAloop .push_back (0);
				graph[cRow][cCol].traceNum .push_back ("NONE");


				tmplink.rowIndx = 0;
				tmplink.colIndx = 0;

				graph[nRow][nCol].inLinks .push_back (tmplink);
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
	cRow = graph.size ()-1;					//last row ..  the new row
	cCol = 0;			//last coloumn ...the new coloumn..the one exceeds the original

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

				graph[nRow][nCol].outLinks .push_back (tmplink);
				graph[nRow][nCol].nAAloop .push_back (0);
				graph[nRow][nCol].traceNum .push_back ("NONE");

				tmplink.rowIndx = nRow;
				tmplink.colIndx = nCol;

				graph[cRow][cCol].inLinks .push_back (tmplink);

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


	//cout<<"Graph is built...information before clean up"<<endl;

	
	//printGraph(graph);

	cleanUp(graph);

	//printGraph(graph);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
void setTopoCandidates(vector<vector<cell> > graph, priority_queue<paths> & topoCandidates)
{
	int kShortest=3000;
	int j;

	// remove end node.... if you use traverse , uncomment this statement
	//graph.pop_back ();


	int nCol = graph[0].size(),
		nRow = graph.size();

	double nPaths=0;			//number of paths

	//sticks represented by Cols and SS represented by row...foreach stick two cols (forward and reverse directions)

	bool done = false;
	node current;
	current.row = 0;
	current.col = 0;
	current.curLink = 0;

	for (j=0; j<nCol; j++)
	{
		current.visited.push_back(0);
		current.topology .push_back (-1);
	}

	while (!done)
		current = traverse1(graph, topoCandidates, current, -1, done, nPaths);


	/*
	//print Candidate topologies
	int i=1;
	while (!topoCandidates.empty ())
	{
		cout<<i<<"  ";
		for (j=0; j< nCol; j++)
			cout<<topoCandidates.top().ssAssigned [j]<<" ";
		cout<<"  Total weight = "<<topoCandidates.top ().w <<endl;
		i++;
		topoCandidates.pop ();

		//if (i % 20 == 0)
		//	getchar();
	}
	*/

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 *		Read mesh file produced by Gorgon software for density skeleton and store it in a vector of points
 *		each point can be represent late in one atom (or one CA atom)
 */
void readMeshSkeleton(string meshFile, vector<Coordinate> &points){

	ifstream inFile;
	string line;
	Coordinate point;

	inFile.open(meshFile.c_str());
	if (!inFile)
	{
		cout<<"================================ in readMeshSkeleton(.......) ====================="<<endl;
		cerr << "Unable to open " << meshFile << endl;
		cout<<"==================================================================================="<<endl;
		exit(1);
	}
	else
	{
		//skip first two lines
		if (!inFile.eof ())
			getline(inFile, line);
		if (!inFile.eof ())
			getline(inFile, line);

		string token;

		while (!inFile.eof())
		{
			getline(inFile, line);
			stringstream lineStream(line);
			lineStream >> token;
			if (token != "3" && token != "4"){
				point.x = atof(token.c_str ());
				lineStream >> token;
				point.y = atof(token.c_str ());
				lineStream>> token;
				point.z = atof(token.c_str ());

				points.push_back (point);
				//cout<<point.x<<" "<<point.y<<" "<<point.z<<endl;
			}
			else
				break;
		}

	}
	inFile.close ();
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Bron-Kerbosc algorithm to find all cliques in a graph
void BronKerbosch(vector<vector<short> > &adjMtrx, vector<vector<short> > &clx, vector<short> R,vector<short> P,vector<short> X){

	if (P.empty() && X.empty ()){
		clx.push_back(R);
	}
	//choose a pivot node from P U X
	//I'll choose the node with maximum number of neighbors
	short u= -1;		//the index of the pivot node
	short maxN = 0;		//maximum number of neighbors
	short i,j,k;

	for (i=0; i<P.size(); i++){
		if (adjMtrx[P[i]].size() > maxN){
			u = P[i];
			maxN = adjMtrx[P[i]].size ();
		}
	}
	for (i=0; i<X.size(); i++){
		if (adjMtrx[X[i]].size() > maxN){
			u=X[i];
			maxN = adjMtrx[X[i]].size();
		}
	}

	//cout<<setw(cntr)<<" "<<"U choosen to be cluster# : "<<u+1<<endl;
	bool found = false;
	i=0;
	while (i<P.size ()){
		vector<short> tmpP, tmpX, tmpR;
		short v = P[i];
		found = false;

		if (!found){
			tmpR = R;
			tmpR.push_back(v);				//R U v
			//P intesect. N(v)
			for (j=0; j<P.size(); j++){
				for (k=0; k<adjMtrx[v].size(); k++){
					if (P[j] == adjMtrx[v][k]){
						tmpP.push_back (P[j]);
						break;
					}
				}
			}
			//X intesect. N(v)
			for (j=0; j<X.size(); j++){
				for (k=0; k<adjMtrx[v].size(); k++){
					if (X[j] == adjMtrx[v][k]){
						tmpX.push_back (X[j]);
						break;
					}
				}
			}

			BronKerbosch(adjMtrx, clx, tmpR, tmpP, tmpX);

			//update P and X
			X.push_back (v);

			P.erase (P.begin () + i);
			i--;

		}
		i++;
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Bron-Kerbosc algorithm to find all cliques in a graph
void BronKerbosch(vector<vector<short> > &adjMtrx,
				  short **(&clx),
				  short *nElements,				//number of leements in each clique
				  short *R,
				  short nR,					//number of elements in R
				  short *P,
				  short nP,					//number of elements in P
				  short *X,
				  short nX,					//number of elements in X
				  short &clxSze){				//number of cliques

	if (nP == 0 && nX == 0){
		//cout<<setw(1)<<" "<<"Clique is found.. {";
		//for (int i=0; i<nR; i++)
		//	cout<<R[i]+1<<" ";
		//cout<<"}"<<endl;

		nElements [clxSze] = nR;
		clx[clxSze] = new short [nR];

		for (int i=0; i<nR; i++){
			clx[clxSze][i] = R[i];
		//	R[i] = -1;
		}
		//nR = 0;
		clxSze++;
	}
	//choose a pivot node from P U X
	//I'll choose the node with maximum number of neighbors
	short u= -1;		//the index of the pivot node
	short maxN = 0;		//maximum number of neighbors
	short i,j,k;

	/*
	//print P
	cout<<endl<<setw(3)<<" "<<"P = {";
	for (i=0; i<nP; i++)
		cout<<P[i]+1<<" ";
	cout<<"}"<<endl;
	//Print R
	cout<<setw(3)<<" "<<"R = {";
	for (i=0;i<nR;i++)
		cout<<R[i]+1<<" ";
	cout<<"}"<<endl;
	//print X
	cout<<setw(3)<<" "<<"X = {";
	for (i=0; i<nX; i++)
		cout<<X[i]+1<<" ";
	cout<<"}"<<endl;
	*/



	for (i=0; i<nP; i++){
		if (adjMtrx[P[i]].size() > maxN){
			u = P[i];
			maxN = adjMtrx[P[i]].size ();
		}
	}
	for (i=0; i<nX; i++){
		if (adjMtrx[X[i]].size() > maxN){
			u=X[i];
			maxN = adjMtrx[X[i]].size();
		}
	}

	//cout<<setw(cntr)<<" "<<"U choosen to be cluster# : "<<u+1<<endl;
	bool found = false;
	i=0;
	while (i<nP){
		short *tmpP, *tmpX, *tmpR, nTmpR=0, nTmpP=0, nTmpX=0;

		tmpR = new short [adjMtrx.size ()];
		tmpP = new short [adjMtrx.size ()];
		tmpX = new short [adjMtrx.size ()];

		for (k=0; k<adjMtrx.size (); k++){
			tmpR[k] = R[k];
			tmpP[k] = -1;
			tmpX[k] = -1;
		}

		if (P[i] != -1){
			short v = P[i];
			//cout<<setw(2)<<" "<<i<<" : working on cluster : "<<P[i]+1<<endl;
			found = false;

			if (!found){
				nTmpR = nR;
				tmpR[nTmpR++] = v;			//R U v

				//cout<<setw(3)<<"set tmpR is  : ";
				//for (k=0; k<nTmpR; k++){
				//	cout<<" "<<tmpR[k];
				//}
				//cout<<endl;


				//P intesect. N(v)
				for (j=0; j<nP; j++){
					for (k=0; k<adjMtrx[v].size(); k++){
						if (P[j] == adjMtrx[v][k]){
							//tmpP.push_back (P[j]);
							tmpP[nTmpP] = P[j];
							nTmpP++;
							break;
						}
					}
				}
				//X intesect. N(v)
				for (j=0; j<nX; j++){
					for (k=0; k<adjMtrx[v].size(); k++){
						if (X[j] == adjMtrx[v][k]){
							//tmpX.push_back (X[j]);
							tmpX[nTmpX] = X[j];
							nTmpX++;
							break;
						}
					}
				}
				/*
				cout<<setw(3)<<" "<<"P intersect N("<<v+1<<")= {";
				for (j=0; j<nTmpP; j++)
					cout<<tmpP[j]+1<<" ";
				cout<<"}"<<endl;
				cout<<setw(3)<<" "<<"X intersect N("<<v+1<<")= {";
				for (j=0; j<nTmpX;j++)
					cout<<tmpX[j]+1<<" ";
				cout<<"}"<<endl;
				*/
				//getchar();


				BronKerbosch(adjMtrx, clx, nElements, tmpR, nTmpR, tmpP, nTmpP, tmpX, nTmpX, clxSze);

				//update P and X
				//X.push_back (v);
				X[nX++] = v;


				/*
				cout<<setw(cntr)<<" "<<"P is : {";
				for (int m=0; m<P.size (); m++)
					cout<<P[m]+1<<" ";
				cout<<"}"<<endl;
				cout<<"  deleting... "<<P[i]+1;
				*/
				//P.erase (P.begin () + i);
				P[i] = -1;
				i--;

				//cout<<" Done.. i="<<i<<endl;
			}
		}
		i++;

		delete [] tmpR;
		delete [] tmpP;
		delete [] tmpX;
		//cout<<setw(cntr)<<" "<<"i= "<<i<<" P.size= "<<P.size ()<<endl;
	}
	//cout<<setw(cntr)<<" "<<"Finish"<<endl;
	//getchar();
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//finding cliques in a graph...given adjacency mtrx
//we use Bron and Kerbosch algorithm
//we keep iterately form cliques untill we reach a point that no more cliques with 3 elements or more can be formed
//nClusters is the number of initial clusters without clusters for SS ends
//adjMtrx is the initial adjacency Matrix b/w initial clusters
void findallCliques(vector<vector<Coordinate> > inClusters,					//the coordinate points of initial clusters
					vector<vector<short> > adjMtrx,						//adjacency matrix of initial clusters
					vector<vector<short> > &clx,						//final cliques
					vector<vector<short> > &clxAdjMtrx,					//final adjacenvy matrix for cliques
					vector<vector<Coordinate> > &clxClusters){				//final coordinate for points form cliques

	int i,j,k,l,m, maxClq;


	//clear previous contents
	clxAdjMtrx.clear ();
	clx.clear ();
	clxClusters.clear();

	//clock_t start, finish;
	//start = clock();

	//initiate sets
	vector<short>	P,			//the set contains vertices
					X,			//the NOT set
					R;			//the set will contain the maximum clique
	for (i=0; i<adjMtrx.size(); i++)
		P.push_back(i);


	//find cliques and save the indices of cluster form cliques in clx
	BronKerbosch(adjMtrx, clx, R,P,X);

	//finish = clock();
	//cout<<"Time taken is : "<<finish-start<<endl;

	/*
	//print all cliques
	cout<<endl<<"Number of cliques found = "<<clx.size ()<<endl;
	for (i=0; i<clx.size (); i++){
		cout<<"Clique# "<<i+1<<" : ";
		for (int j=0; j<clx[i].size (); j++)
			cout<<clx[i][j]+1<<" ";
		cout<<endl;
	}
	*/

	/*
	start = clock();
	short	**clxAr,
			*nElements,				//number of elements in each clique
			 *P1, *R1, *X1,
			 nP1=adjMtrx.size (), nR1=0, nX1=0,
			 clxSize=0;

	clxAr = new short *[adjMtrx.size ()];			//max number of cliques is the number of original clusters
	nElements = new short [adjMtrx.size ()];
	P1 = new short [adjMtrx.size ()];
	R1 = new short [adjMtrx.size ()];
	X1 = new short [adjMtrx.size ()];

	for (i=0; i<adjMtrx.size(); i++){
		P1[i] = i;
		R1[i] = -1;
		X1[i] = -1;
	}
	BronKerbosch(adjMtrx, clxAr, nElements, R1, nR1, P1, nP1, X1, nX1, clxSize);

	finish = clock();
	cout<<"Time taken is : "<<finish-start<<endl;

	//print all cliques
	cout<<endl<<"Number of cliques found = "<<clxSize<<endl;
	for (i=0; i<clxSize; i++){
		cout<<"Clique# "<<i+1<<" : ";
		for (int j=0; j<nElements[i]; j++)
			cout<<clxAr[i][j]+1<<" ";
		cout<<endl;
	}
	*/

	//delete cliques with only 1 or two nodes...
	i=0;
	maxClq=0;
	while (i<clx.size ()){
		if (clx[i].size () > maxClq)
			maxClq=clx[i].size ();			//get the size of largest clique
		//delete small cliques (of size 1 or 2)
		if (clx[i].size() < 3){
			clx.erase (clx.begin () + i);
			i--;
		}
		i++;
	}

	//print all cliques
    /*
	cout<<"Number of cliques found = "<<clx.size ()<<endl;
	for (i=0; i<clx.size (); i++){
		cout<<"Clique# "<<i+1<<" : ";
		for (int j=0; j<clx[i].size (); j++)
			cout<<clx[i][j]+1<<" ";
		cout<<endl;
	}
	*/

	if (clx.size ()){
		//merge cliques
		//merge any two cliques of size 3...if they share 2 nodes
		//merge any two cliques of size 4...if they share 3 nodes
		//merge any two cliques of size 5 ..if they share 4 nodes and so on
		short sharedNodes,
				diff,							//the different cluster
				clqSze;

		short ** clxArr;						//array version of cliques (faster to deal with)
		short *nEntries;						//number of entries (clusters) in each clique;
		int nClxArr = clx.size ();				//number of cliques

		bool *isValidClq;						//a flag is true if the clique is valid
		bool *processedBefore;					//did this clique merged with another clique before

		short	**tmpClq;						//new produced clique from merging two prev. cliques
		short *tmpClqEntries;					//number of entries in the tmp cliq vector
		bool foundSimilar,						//flag is set if any two cliques found to be very similar
				contMerging = true;				//keep mergning if u still have cliques to merge

		//copy cliques data into clxArr
		clxArr = new short *[nClxArr];
		nEntries = new short[nClxArr];

		for (i=0; i<nClxArr; i++){
			nEntries[i] = clx[i].size ();
			clxArr[i] = new short [nEntries[i]];
			for (j=0; j<nEntries[i]; j++)
				clxArr[i][j] = clx[i][j];
		}

		clqSze=3;
		while (contMerging && clqSze < 9){

			//cout<<"clqSze= "<<clqSze<<endl;

			//initiate flags
			isValidClq = new bool [nClxArr];
			processedBefore = new bool [nClxArr];
			for (i=0; i<nClxArr; i++){
				isValidClq[i] = true;
				processedBefore[i] = false;
			}
			contMerging=false;

			//expected maximum number of merged cliques would be double the original number of cliques
			tmpClq = new short *[10*nClxArr];
			tmpClqEntries = new short [10*nClxArr];

			int tmpClqCntr = 0;		//initiate number of new cliques
			i=0;
			while (i<nClxArr){
				foundSimilar=false;
				bool pushedBefore = false;
				if (nEntries[i] == clqSze && isValidClq[i]){
					j=i+1;
					while (j<nClxArr){
						if (nEntries[j] == clqSze && isValidClq[j]){

							//find the different cluster b.w two cliques if any
							sharedNodes = 0;
							diff = -1;
							bool missing;
							for (l=0; l<nEntries[j];l++){
								missing = true;
								for (k=0; k<nEntries[i]; k++){
									if (clxArr[i][k] == clxArr[j][l]){
										sharedNodes++;
										missing = false;
										break;
									}
								}
								if (missing)
									diff = clxArr[j][l];
							}
							//if two cliques share n-1 clusters
							if (sharedNodes > clqSze-2){
								foundSimilar=true;
								if (diff == -1){			//to avoid any two identical cliques after many merging steps
									if (!pushedBefore){
										//copy the clique into tmpClq vector and update the counter
										tmpClq[tmpClqCntr] = new short [nEntries[i]];
										tmpClqEntries[tmpClqCntr] = nEntries[i];
										for (k=0; k<nEntries[i]; k++)
											tmpClq[tmpClqCntr][k] = clxArr[i][k];
										tmpClqCntr++;

										pushedBefore = true;
										//cout<<"clq# "<<i+1<<" and clq# "<<j+1<<" ARE IDENTICALS.. pushed to new list to clq# "<<tmpClq.size ()<<endl;
									}

									//delete j clique
									isValidClq[j] = false;

								}
								else{
									//copy the clique into tmpClique vector with the different cluster
									tmpClqEntries[tmpClqCntr] = nEntries[i]+1;
									tmpClq[tmpClqCntr] = new short [nEntries[i]+1];
									for (k=0; k<nEntries[i]; k++)
										tmpClq[tmpClqCntr][k] = clxArr[i][k];
									tmpClq[tmpClqCntr][k] = diff;					//add different cluster
									tmpClqCntr++;

									contMerging = true;
									processedBefore[j] = true;
									//cout<<"clq "<<i+1<<" and clq "<<j+1<<" have "<<sharedNodes<<" shared nodes.. diff is "<<diff+1<<".. merged to clq#"<<tmpClq.size ()<<endl;
								}
							}
						}
						j++;
					}
				}
				if (!foundSimilar && isValidClq[i]){
					if (!processedBefore[i]){
						//add the clique as is
						tmpClqEntries[tmpClqCntr] = nEntries[i];
						tmpClq[tmpClqCntr] = new short [nEntries[i]];
						for (k=0; k<nEntries[i]; k++)
							tmpClq[tmpClqCntr][k] = clxArr[i][k];
						tmpClqCntr++;
						//cout<<"clq "<<i+1<<" has no similar clq.... pushed to clq# "<<tmpClq.size ()<<endl;
					}
				}

				i++;

			}
			clqSze++;
			//cout<<"  Done.."<<endl;

			//re-copy everything to clxArr
			delete [] *clxArr;
			delete [] clxArr;
			delete [] nEntries;

			clxArr = new short *[tmpClqCntr];
			nEntries = new short [tmpClqCntr];
			for (i=0; i<tmpClqCntr; i++){
				clxArr[i] = new short [tmpClqEntries[i]];
				nEntries[i] = tmpClqEntries[i];
				for (j=0; j<tmpClqEntries[i]; j++)
					clxArr[i][j] = tmpClq[i][j];
			}
			nClxArr = tmpClqCntr;

			delete [] *tmpClq;
			delete [] tmpClq;
			delete [] tmpClqEntries;

			delete [] isValidClq;
			delete [] processedBefore;

		}


		/*
		cout<<"Number of cliques found = "<<clx.size ()<<endl; getchar();
		for (i=0; i<clx.size (); i++){
			cout<<"Clique# "<<i+1<<" : ";
			for (int j=0; j<clx[i].size (); j++)
				cout<<clx[i][j]+1<<" ";
			cout<<endl;
		}
		getchar();getchar();
		*/


		//delete any clique completely contained in another clique
		//cout<<"# of cliques = "<<nClxArr<<"  Deleting duplicate cliques... ";

		bool *isValid;;
		isValid = new bool [nClxArr];
		for (i=0; i<nClxArr; i++){
			isValid[i] = true;
		}

		i=0;
		while (i<nClxArr){
			if (isValid[i]){
				j=i+1;
				while (j<nClxArr){
					if (isValid[j]){
						short sharedNodes=0;
						for (k=0;k<nEntries[i]; k++){
							for (l=0; l<nEntries[j]; l++){
								if (clxArr[i][k] == clxArr[j][l])
									sharedNodes++;
							}
						}
						if (sharedNodes == nEntries[i]){
							isValid[i] = false;
							break;
						}
						if (sharedNodes == nEntries[j]){
							isValid[j] = false;
						}
					}
					j++;
				}
			}
			i++;
		}

		vector<vector< short> > tmpCliques;
		for (i=0; i<nClxArr; i++){
			vector<short> tmpEntry;
			if (isValid[i]){
				for (j=0; j<nEntries[i]; j++)
					tmpEntry.push_back(clxArr[i][j]);
				tmpCliques.push_back (tmpEntry);
			}

		}
		clx = tmpCliques;
		tmpCliques.clear ();
		//cout<<"Done."<<endl;
	}

	//print valid cliques
	/*
	cout<<"After DELETING SUB CLIQUES...."<<endl;
	for (i=0; i<clx.size (); i++){
		cout<<"Clq# "<<i+1<<" : ";
		for (j=0; j<clx[i].size (); j++){
			cout<<clx[i][j]+1<<" ";
		}
		cout<<endl;
	}
	*/

	//resolve shared clusters....if a cluster appear in two cliques....keep it in only one of them according to the closest distances
	//cout<<"# of cliques = "<<clx.size()<<" Deleting completely shared cliques... ";
	float dist;
	for (i=0; i<clx.size (); i++){
		for (j=i+1; j<clx.size (); j++){
			k=0;
			while (k<clx[i].size ()){
				l=0;
				while (l<clx[j].size ()){
					bool deleted=false;
					//cout<<i+1<<"-"<<k+1<<" with "<<j+1<<"-"<<l+1<<endl;
					if (clx[i][k] == clx[j][l]){
						//cout<<"   found a shared cluster "<<clx[i][k]+1<<endl;
						//find to which clique this cluster is closest
						float minDist = 99999.0;
						short x, p;
						for (m=0; m<clx[i].size (); m++){
							if (m!= k){
								for (p=0; p<inClusters[clx[i][m]].size (); p++){
									for (x=0; x<inClusters[clx[i][k]].size (); x++){
										dist = getDistance(inClusters[clx[i][m]][p], inClusters[clx[i][k]][x]);
										if ( dist < minDist)
											minDist = dist;
									}
								}
							}
						}
						deleted = false;
						for (m=0; m<clx[j].size (); m++){
							if (m!=l){
								for (p=0; p<inClusters[clx[j][m]].size (); p++){
									for (x=0; x<inClusters[clx[j][l]].size (); x++){
										dist = getDistance(inClusters[clx[j][m]][p], inClusters[clx[j][l]][x]);
										if ( dist < minDist){
											//cout<<" minDist now = "<<minDist<<endl;
											//keep it in clique j..delete it from clique i
											clx[i].erase (clx[i].begin () + k);
											k--;
											deleted = true;
											//cout<<"  clq# "<<i+1<<" size now = "<<clx[i].size ()<<" k= "<<k<<endl;
											break;
										}
									}
									if (deleted)
										break;
								}
							}
							if (deleted)
								break;
						}

						if (!deleted){
							clx[j].erase (clx[j].begin () + l);
							l--;
							//cout<<"clq# "<<j+1<<" size now = "<<clx[j].size ()<<endl;
						}
					}
					if (deleted)
						break;
					l++;
				}
				k++;
			}
		//	getchar();
		}
	}
	//cout<<"Done."<<endl;

	//delete empty clusters
	//cout<<"Deleting empty cliques... ";
	i=0;
	while (i<clx.size ()){
		if (clx[i].size () == 0){
			clx.erase (clx.begin () + i);
			i--;
		}
		i++;
	}
	//cout<<"Done."<<endl;

	//print valid cliques
    /*
	cout<<"AFTER DELETING SHARED CLUSTERS.."<<endl;
	for (i=0; i<clx.size (); i++){
		cout<<"Clq# "<<i+1<<" : ";
		for (j=0; j<clx[i].size (); j++){
			cout<<clx[i][j]+1<<" ";
		}
		cout<<endl;
	}
	*/


	//copy all other clusters that were not part of cliques
	vector<short> tmpClx(1);
	bool found;
	for (k=0; k<adjMtrx.size (); k++){
		found = false;
		for (i=0; i<clx.size (); i++){
			for (j=0;j<clx[i].size (); j++){
				if (clx[i][j] == k){
					found = true;
					break;
				}
			}
			if (found)
				break;
		}
		if (!found){
			tmpClx[0] = k;
			clx.push_back (tmpClx);
		}
	}

	//merge clusters (points form clusters) to form cliques
	vector<Coordinate> nClq;		//new clique
	Coordinate nCentroid;
	for (i=0; i<clx.size(); i++){
		nClq.clear ();
		nCentroid.x = 0.0;
		nCentroid.y = 0.0;
		nCentroid.z = 0.0;

		for (j=0; j<clx[i].size(); j++){
			for (k=0; k<inClusters[clx[i][j]].size ()-1; k++)			//don't save the centroid of the old cluster
				nClq.push_back(inClusters[clx[i][j]][k]);

			//add the centroid of this cluster to the new centroid to be calculated
			nCentroid.x += inClusters[clx[i][j]][k].x;
			nCentroid.y += inClusters[clx[i][j]][k].y;
			nCentroid.z += inClusters[clx[i][j]][k].z;
		}
		//new centroid
		nCentroid.x /= clx[i].size ();
		nCentroid.y /= clx[i].size ();
		nCentroid.z /= clx[i].size ();

		//save new centroid
		nClq.push_back (nCentroid);

		clxClusters.push_back (nClq);
	}

	//print final cliques after update
	/*
	cout<<"Final cliques mapping.."<<endl;
	for (i=0; i<clx.size (); i++){
		cout<<"Clique# "<<i+1<<" : ";
		for (j=0; j<clx[i].size (); j++){
			cout<<clx[i][j]+1<<" ";
		}
		cout<<endl;
	}
	*/

	clxAdjMtrx.resize (clx.size ());		//adjacency matrix for clx...neighbors are the indices of cliques not the initial clusters
	// build adjacency matrix for cliques
	for (i=0; i<clx.size (); i++){
		for (j=0; j<clx[i].size (); j++){
			for (l=0; l<adjMtrx[clx[i][j]].size (); l++){
				for (k=i+1; k<clx.size (); k++){
					for (m=0; m<clx[k].size (); m++){
						if ( clx[i][j] != clx[k][m] && adjMtrx[clx[i][j]][l] == clx[k][m]){
							clxAdjMtrx[i].push_back (k);
							clxAdjMtrx[k].push_back (i);
							break;
						}
					}
				}
			}
		}
	}
	//remove duplicates in clxAdjMtrx
	for (i=0; i<clxAdjMtrx.size (); i++){
		for (j=0; j<clxAdjMtrx[i].size (); j++){
			k=j+1;
			while (k<clxAdjMtrx[i].size ()){
				if (clxAdjMtrx[i][j] == clxAdjMtrx[i][k]){
					clxAdjMtrx[i].erase (clxAdjMtrx[i].begin () + k);
					k--;
				}
				k++;
			}
		}
	}

	//print clx Adjacency matrix
	/*
	for (i=0; i<clxAdjMtrx.size (); i++){
		cout<<"N (clq# "<<i+1<<" )  : ";
		for (j=0; j<clxAdjMtrx[i].size (); j++){
			cout<<clxAdjMtrx[i][j]+1<<" ";
		}
		cout<<endl;
	}
	*/
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 *		find continueous traces in the density map starting from a particular point....starting point,
 *		2 dimension matrix have distances b/w each pair of clusters
 *		sDist is a matrix has the shortest distance b/w each pair of clusters
 */
void findTracePaths(vector<vector<short> > &adjMtrx,			//adjacency matrix for clusters graph
					float **sDist,								//matrix has centroid-centroid distances b/w each pair of clusters
					short startPnt,								//the starting point where you want to start ur tracing
					short nClusters,							//number of clusters without clusters for SS ends
					vector<vector<short> > &paths,				//trace paths found for density map
					vector<float> &pathsCosts,					//parallel vector with paths contains the cost for each trace path
					float maxLngth = 500.0){					//maximum length of path allowed


	pair<short, short>	*open;			//open set represetns clusters that are visited already and the number of children in the stack
	short nOpen;						//number of active elements in open
	pair<short, short>	*stack;			//clusters being processed and the parent for each cluster
	short nStack;						//number of active elements in the stack (also the index of the top of the stack)
	float				*openLength;	//the length of the trace so far for each node in open

	bool *visitFlag;
	pair<short, short> tmp, top;
	vector<short> onePath;
	float pathLength = 0.0;
	short i;
	short nChildMax = 0;			//max number of chilren for a node


	//get maximum number of children for a node
	for (i=0; i<adjMtrx.size (); i++)
		if (adjMtrx[i].size () > nChildMax)
			nChildMax = adjMtrx[i].size ();

	//initialize data structures
	stack = new pair<short, short> [adjMtrx.size ()];
	open = new pair<short, short> [adjMtrx.size () * nChildMax];
	openLength = new float [adjMtrx.size() * nChildMax];
	visitFlag = new bool [adjMtrx.size ()];
	nStack = 0;
	nOpen = 0;

	for (i=0; i<adjMtrx.size (); i++)
		visitFlag[i] = false;


	//push the starting point (the indx of starting point) to the stack
	tmp.first = startPnt;
	tmp.second = -1;				//no parent

	//stack.push_back (tmp);
	stack[nStack++] = tmp;

	while (nStack!=0){
		//pop the top of the stack and put its children (those clusters are satisfying threshold distance)
		top = stack[nStack-1];

		//cout<<top.first +1;			//cluster being processed

		//decrement number of neighbors in the stack for its parent
		if (top.second != -1){		//if it is not the start pnt
			open[top.second].second--;
			//cout<<"  PathLength= "<<openLength[nOpen-1];
			//cout<<"     ( parent: "<<open[top.second].first+1<<"  p.Indx= "<<top.second<<"      #children: "<<open[top.second].second<<" ).";
		}

		//cout<<endl;

		//stack.pop_back ();
		nStack--;

		//initialize number of children and push to open
		top.second = 0;
		//open.push_back (top);
		open[nOpen++] = top;

		visitFlag[top.first] = true;		//set visit flag

		//calculate the length of the path so far
		if (nOpen > 1){
			//get distance b/w last two clusters in open list
			short	openLastIndx = nOpen-1,
					addedClusterIndx = open[openLastIndx].first,
					clusterIndx = open[openLastIndx-1].first;

			//distance b/w centroids of the last two clusters in open list
			pathLength =	openLength[openLastIndx-1] + sDist[clusterIndx][addedClusterIndx];

			//openLength.push_back (pathLength);
			openLength [nOpen-1] = pathLength;
		}
		else
			//openLength.push_back (0.0);				//path length for first node
			openLength [nOpen-1] = 0.0;

		//don't continue with this node if you already exceeded the max length allowed
		if (pathLength > maxLngth && top.first < nClusters){
			//cout<<" The Path exceeds the maximum length allowd which is : "<<maxLngth<<endl;
			visitFlag[open[nOpen-1].first] = false;
			//open.pop_back ();
			nOpen--;
			//openLength.pop_back ();
			continue;
		}

		//if the cluster we r working on is one of SS ends. then do not generate its children b/s this would be a trace end
		if ( top.first < nClusters || top.first == startPnt){
			//generate children (neighbors) for the cluster
			for (i=0 ; i<adjMtrx[top.first].size(); i++){
				//cout<<"  "<<adjMtrx[top.first][i]+1<<" ";

				if (!visitFlag[adjMtrx[top.first][i]]){
					tmp.first  = adjMtrx[top.first][i];
					tmp.second = nOpen-1;			//save the index of the parent in open set
					//stack.push_back (tmp);
					stack[nStack++] = tmp;					//push neighbor to stack
					open[tmp.second].second++;				//increment number of neighbors in stack
					//cout<<"pushed to stack.   parentIndx= "<<tmp.second<<" ( "<<open[tmp.second].first+1<<" )."<<endl;
				}
				//else
				//	cout<<" found in open..."<<endl;
			}
		}


		//if we reach a leaf node...with no neighbor in the stack or SS end
		if (open[nOpen-1].second<1){
			//cout<<"  a leaf node has been reached "<<open[nOpen-1].first+1<<" .... "<<endl;
			onePath.clear();
			for (i=0; i<nOpen-1; i++)
				onePath.push_back(open[i].first);
			//push last cluster
			onePath.push_back (open[i].first );

			//save the path and its cost
			paths.push_back (onePath);
			pathsCosts.push_back (openLength[i]);

			//delete all nodes in open set with 0 neighbors
			i=nOpen-1;
			while (i>=0 && open[i].second <=0){
				//cout<<"       "<<open[i].first +1<<" deleted from open"<<endl;
				visitFlag[open[i].first] = false;
				//open.pop_back ();
				nOpen--;
				//openLength.pop_back ();
				i--;
			}
			//cout<<"    Open contents"<<endl;
			//for (i=0; i<nOpen;i++)
			//	cout<<"("<<open[i].first+1 <<","<<open[i].second <<") ";
			//cout<<endl;
		}
	}

	delete [] open;
	delete [] openLength;
	delete [] stack;
	delete [] visitFlag;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void buildAdjMtrx(vector<vector<Coordinate> > &clusters, float **(&cDist), vector<vector<short> > &adjMtrx, float apix, short nClustersOriginal, float continuatyTHR){

	int i,k;
	int nClusters = clusters.size ();

	Coordinate cPnt, nPnt;			//the two points we are measuring distance b/w...in MAP indexing system

	//initiate cluster distance matrices
	cDist = AllocateDynamicArray<float> (nClusters, nClusters);

	for (i=0; i<nClusters; i++){
		for (k=0; k<nClusters; k++)
			cDist[i][k] = 9999.0;
	}

	//
	//		find minimum distance b/w each pair of clusters
	//
	adjMtrx.clear ();
	adjMtrx.resize (nClusters);			//the size is : the initial clusters and the clusters for SSs' ends
	float dist, distTHR;
	bool cont;
	for (i=0; i<nClusters-1; i++){
		for (k= i+1; k<nClusters; k++){
			if (k >= nClustersOriginal)
				distTHR = 3.0;
			else
				distTHR = continuatyTHR;
			cont = true;
			for (int j=0; j<clusters[i].size (); j++){
				cPnt = clusters[i][j];
				for (int l=0; l<clusters[k].size (); l++){
					nPnt = clusters[k][l];
					dist  = round(getDistance(cPnt, nPnt)/apix,3);
					if (dist <= distTHR){
						adjMtrx[i].push_back (k);
						adjMtrx[k].push_back (i);
						cont = false;
						break;
					}
				}
				if (!cont)
					break;
			}
			//find the distance b/w centroids
			cDist[i][k] = cDist[k][i] = round(getDistance(clusters[i][clusters[i].size ()-1], clusters[k][clusters[k].size ()-1]),3);
		}
	}

	//cout<<"Building adjacency matrix is done.."<<endl;

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//select the best trace path among all paths b/w two SS ends
// this function assumes that all paths of a particular end are follow each other
float getBestPathTrace(vector<vector<Coordinate> > &clusters, short **tracePaths, int nTraces, short *tracesSize, float *traceLengths, short nClusters, float loopLength, int firstEnd, int secondEnd, int &bestTraceIndx1, int &bestTraceIndx2){
	int i,j,k;
	//get the starting and ending indeces of the paths belong to first End
	int firstIndx1 = -1, firstIndx2=-1, secondIndx1=-1, secondIndx2=-1;

	//always tracePaths[i][0] is the starting end
	i=0;
	while (i<nTraces && tracePaths[i][0] != firstEnd) i++;
	firstIndx1 = i;
	while (i<nTraces && tracePaths[i][0] == firstEnd) i++;
	firstIndx2 = i-1;

	i=0;
	while (i<nTraces && tracePaths[i][0] != secondEnd) i++;
	secondIndx1 = i;
	while (i<nTraces && tracePaths[i][0] == secondEnd) i++;
	secondIndx2 = i-1;


	float mapDist = getDistance(clusters[firstEnd][clusters[firstEnd].size ()-1], clusters[secondEnd][clusters[secondEnd].size ()-1]);
	float bestWeight = fabs( mapDist - loopLength);

	if (mapDist > 7.0)
		bestWeight += 50.0;			//add some weights for missing trace ... for far away sticks

	int lastCluster1Indx, lastCluster2Indx, centroid1Indx, centroid2Indx;
	Coordinate trace1LastPoint, trace2LastPoint;			//the coordinate of last point in both traces

	//if (firstIndx1 == 932 && secondIndx1 == 1082)
	//	cout<<"1stIndx1= "<<firstIndx1<<"  2ndIndx1= "<<secondIndx1<<"  1stIndx2= "<<firstIndx2<<"  2ndIndx2= "<<secondIndx2<<"  initial bestWeight= "<<bestWeight<<endl;

	for (i=firstIndx1; i<=firstIndx2; i++){
		float matchingLength = 9999.0;
		if (tracePaths[i][tracesSize[i]] == secondEnd){
			//cout<<"Found complete path  trace Length= "<<traceLengths[i]<<"   and loop length= "<<loopLength<<endl;
			if (traceLengths[i] <= (loopLength + 5.0)){
				if (fabs(loopLength-traceLengths[i]) < bestWeight){
					bestWeight = fabs(loopLength - traceLengths[i]);
					bestTraceIndx1 = i;
					bestTraceIndx2 = i;
				}
				//bestWeight -= 10.0;					//reward a connected full path
			}
		}
		else{
			float THR = 10.0;			//discontinuty threshold distance (for real maps use 15. simulated maps use 8.0)
			lastCluster1Indx	= tracePaths[i][tracesSize[i]];
			centroid1Indx		= clusters[lastCluster1Indx].size ()-1;
			trace1LastPoint		= clusters[lastCluster1Indx][centroid1Indx];

			for (j=secondIndx1; j<=secondIndx2; j++){

				//if (firstIndx1 == 932 && secondIndx1 == 1082)
				//	cout<<"checking path "<<i+1<<" with Path "<<j+1<<endl;

				matchingLength = traceLengths[i];
				if ((tracePaths[j][tracesSize[j]] < nClusters) &&		//  || tracePaths[j].size () == 1) &&
					(tracePaths[i][tracesSize[i]] < nClusters)){		//  || tracePaths[i].size () == 1)){

					//deal with discontinue traces
					matchingLength += traceLengths[j];

					//if (firstIndx1 == 932 && secondIndx1 == 1082)
					//	cout<<"matchingLength= "<<matchingLength<<"  loopLength= "<<loopLength<<endl;

					if (matchingLength <= loopLength){

						lastCluster2Indx	= tracePaths[j][tracesSize[j]];
						centroid2Indx		= clusters[lastCluster2Indx].size ()-1;
						trace2LastPoint		= clusters[lastCluster2Indx][centroid2Indx];

						if (getDistance(trace1LastPoint, trace2LastPoint) < THR){
							matchingLength		+= getDistance(trace1LastPoint, trace2LastPoint);

							//if (firstIndx1 == 932 && secondIndx1 == 1082)
							//	cout<<"matchingLength now = "<<matchingLength<<endl;

							if (matchingLength <= (loopLength + 5.0)){
								if (fabs(loopLength-matchingLength) < bestWeight){
									bestWeight = fabs(loopLength-matchingLength);
									bestTraceIndx1 = i;
									bestTraceIndx2 = j;
								}
							}
						}
					}
				}
				//if (firstEnd == 76) {getchar();}
			}
		}
	}

	//cout<<"    final bestWeight = "<<bestWeight<<endl;


	return bestWeight;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void findSticksEndsOnSkeleton(Map inSkeleton, vector<vector<Coordinate> >  ssEdges, pair<Coordinate, Coordinate> *ssEnds, string outPath){

	int i, j, irow, icol, islc, startIndx, endIndx, delta;
	Coordinate skeletonPnt;
	short ix, ix_, iy, iy_, iz, iz_;		//skeleton matrix indeces

	Protein trace;

	//convert 3d points to MAP indexing system
	for (i=0; i<ssEdges.size (); i++){
		for (j=0; j<ssEdges[i].size (); j++){
			ssEdges[i][j].x = (ssEdges[i][j].x - inSkeleton.hdr.xorigin)/inSkeleton.apixX;
			ssEdges[i][j].y = (ssEdges[i][j].y - inSkeleton.hdr.yorigin)/inSkeleton.apixY;
			ssEdges[i][j].z = (ssEdges[i][j].z - inSkeleton.hdr.zorigin)/inSkeleton.apixZ;
		}
	}

	//find the closest point for each stick End
	for (i=0; i<2*ssEdges.size (); i++){

		Atom tmpAtom;
		AminoAcid tmpAA;
		tmpAA.num = i+1;

		//Coordinate closestPnt;
		float	dist,// minDist = 9999.0,
				stickThirdLength = getDistance(ssEdges[i/2][0], ssEdges[i/2][ssEdges[i/2].size ()-1])/3;

		//convert stick End to MAP indexing
		if (i%2==0){
			startIndx = 0;
			delta = 1;
			//determine which point on stick near the endpoint to compare with
			for (j=0; j<ssEdges[i/2].size (); j++){
				dist = getDistance(ssEdges[i/2][startIndx], ssEdges[i/2][j]);

				if (dist > 2.5 || dist> stickThirdLength){
					endIndx = j;
					break;
				}
			}
		}
		else{
			startIndx = ssEdges[i/2].size ()-1;
			delta = -1;
			//determine to which point on stick near the endpoint to compare with
			for (j=startIndx; j>-1; j--){
				dist = getDistance(ssEdges[i/2][startIndx], ssEdges[i/2][j]);

				if (dist > 2.5 || dist > stickThirdLength){
					endIndx = j;
					break;
				}
			}
		}

		ix = ssEdges[i/2][startIndx].x - 8;
		ix_ = ssEdges[i/2][startIndx].x + 8;
		if (ix<0)
			ix = 0;
		if (ix_ > inSkeleton.numRows())
			ix_ = inSkeleton.numRows();

		iy = ssEdges[i/2][startIndx].y - 8;
		iy_ = ssEdges[i/2][startIndx].y + 8;
		if (iy < 0)
			iy = 0;
		if (iy_ > inSkeleton.numCols())
			iy_ = inSkeleton.numCols();

		iz = ssEdges[i/2][startIndx].z - 8;
		iz_ = ssEdges[i/2][startIndx].z + 8;
		if (iz < 0)
			iz = 0;
		if (iz_ > inSkeleton.numSlcs())
			iz_ = inSkeleton.numSlcs();

		bool stop = false;
		float proj;
		Coordinate nPnt;

		float *closestDist;
		Coordinate *closestPnt;
		closestDist = new float [abs(endIndx-startIndx)+1];
		closestPnt = new Coordinate [abs(endIndx-startIndx)+1];

		for (j=0; j<abs(endIndx-startIndx)+1; j++){
			closestDist[j] = 9999.0;
		}
		short cntr=0;
		for (j=startIndx; j!=endIndx+delta; j+=delta){
			for (irow=ix; irow<ix_; irow++){
				skeletonPnt.x = irow;
				for (icol=iy; icol<iy_; icol++){
					skeletonPnt.y = icol;
					for (islc=iz; islc<iz_; islc++){
						skeletonPnt.z = islc;
						if (inSkeleton.cube[irow][icol][islc] > 0.0){

							dist = getDistance(ssEdges[i/2][j], skeletonPnt);

							//cout<<"stick# : "<<i/2+1<<"  j= "<<j<<"  dist= "<<dist<<"  "<<skeletonPnt.x<<" "<<skeletonPnt.y<<" "<<skeletonPnt.z<<endl;
							if (closestDist[cntr] > dist){
								if (dist < 0.5){
									closestDist[cntr] = dist;
									closestPnt[cntr] = skeletonPnt;
									stop = true;
									break;
								}
								else{
									proj = linePointIntersection(ssEdges[i/2][0], ssEdges[i/2][ssEdges[i/2].size ()-1], skeletonPnt, nPnt);
									if (proj>=0 && proj<=1){
										closestDist[cntr] = dist;
										closestPnt[cntr] = skeletonPnt;
									}
								}
							}
						}
					}
					if (stop)
						break;
				}
				if (stop)
					break;
			}
			if (stop)
				break;
			cntr++;
		}

		//find the skeleton point within 1.5 Angstrom if any
		for (j=0; j<abs(endIndx-startIndx)+1; j++){
			if (closestDist[j] <= 1.5){
				closestPnt[0] = closestPnt[j];
				break;
			}

		}

		//if no skeleton point within 1.5 Angstrom.....find the closest point then
		if (j == abs(endIndx-startIndx)+1){
			for (j=1; j<abs(endIndx-startIndx)+1; j++){
				if (closestDist[j] < closestDist[0]){
					closestDist[0] = closestDist[j];
					closestPnt[0] = closestPnt[j];
				}
			}
		}

		//if no closest point around the end of the stick (within 8*8*8 cube) then return the original end of the stick as the closest point in the skeleton...
		if (closestDist[0] == 9999.0){
			//cout<<"could not find a close point.."<<endl;
			closestPnt[0] = ssEdges[i][startIndx];
		}

		//convert the point to PDB coordinate system
		closestPnt[0].x = closestPnt[0].x * inSkeleton.apixX + inSkeleton.hdr.xorigin;
		closestPnt[0].y = closestPnt[0].y * inSkeleton.apixY + inSkeleton.hdr.yorigin;
		closestPnt[0].z = closestPnt[0].z * inSkeleton.apixZ + inSkeleton.hdr.zorigin;


		if (startIndx == 0)
			ssEnds[i/2].first = closestPnt[0];
		else
			ssEnds[i/2].second = closestPnt[0];

		//copy the closest point to an atom coordinate to be written later to a veiwable file by chimera
		tmpAtom.coord = closestPnt[0];
		tmpAA.atoms.push_back(tmpAtom);
		trace.AAs.push_back(tmpAA);
	}

	//write local peaks into a pdb file
	//trace.writePDB(outPath + "_EndPoints.pdb", 1, trace.numOfAA());							//all points


	/*
	Protein densityPDB;

	//copy density voxels to a viewable PDB file
	Atom tmpAtom;
	AminoAcid tmpAA;
	tmpAtom.name = " O  ";
	tmpAA.chr3 = "HOH";
	for (irow=0; irow<inSkeleton.numRows(); irow++){
		for (icol=0; icol<inSkeleton.numCols(); icol++){
			for (islc=0; islc<inSkeleton.numSlcs(); islc++){
				if (inSkeleton.cube[irow][icol][islc] > 0.0){
					tmpAA.atoms.clear();
					tmpAtom.coord.x = irow*inSkeleton.apixX + inSkeleton.hdr.xorigin;
					tmpAtom.coord.y = icol*inSkeleton.apixY + inSkeleton.hdr.yorigin;
					tmpAtom.coord.z = islc*inSkeleton.apixZ + inSkeleton.hdr.zorigin;
					tmpAA.atoms.clear();
					tmpAA.atoms.push_back(tmpAtom);
					densityPDB.AAs.push_back(tmpAA);
				}
			}
		}
	}
	//write local peaks into a pdb file
	densityPDB.writePDB(outPath + "_skeletonPoints.pdb", 1, densityPDB.numOfAA());
	*/

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void deleteEdgeDensityNaive(vector<vector<Coordinate> > &clusters, vector<Coordinate> ssEdge, pair<Coordinate, Coordinate> stickEndsOnDensity, float radius){

	short j, i, k, startIndx = -1, endIndx=-1;

	//determine indeces of points to start deleting from
	for (j=1; j<ssEdge.size (); j++){
		if (getDistance(stickEndsOnDensity.first, ssEdge[j]) > 2){
			startIndx = j;
			break;
		}
	}
	for (j=ssEdge.size ()-2; j>-1; j--){
		if (getDistance(stickEndsOnDensity.second, ssEdge[j]) > 2){
			endIndx = j;
			break;
		}
	}

	for (i=startIndx; i<endIndx; i++){
		for (j=0; j<clusters.size (); j++){
			k=0;
			while (k<clusters[j].size ()){
				if (getDistance(clusters[j][k], ssEdge[i]) < radius){
					clusters[j].erase (clusters[j].begin () + k--);
				}
				k++;
			}
		}
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void deleteEdgesDensity(Map &inMRC, vector<vector<Coordinate> > &clusters, vector<vector<Coordinate> > ssEdges, pair<Coordinate, Coordinate>  *ssEnds, string outPath){

	int i,j, k, sourceCluster;
	findSticksEndsOnSkeleton(inMRC, ssEdges, ssEnds, outPath);

	/*
	 *		find shortest distances b/w each pair of clusters and then build adjMtrx
	 */
	float continuatyTHR = 1.1 * inMRC.apixX;
	short nClusters = clusters.size ();
	vector<vector<short> > adjMtrx;
	float **cDist;

	buildAdjMtrx(clusters, cDist, adjMtrx, inMRC.apixX, nClusters, continuatyTHR);


	float *pairDist;							//the maximum length of a path allowed b.w any two points
	float *tracesLengths;						//the cost to reach each cluster
	short *prevCluster;							//the prev cluster in the trace

	short *nodes;
	float **linksW;

	pairDist = new float [nClusters];
	prevCluster = new short [nClusters];
	tracesLengths = new float [nClusters];

	//initialization links weights
    nodes = new short[nClusters];
	linksW = AllocateDynamicArray<float> (nClusters, nClusters);			//create links weights
	for (i=0; i<nClusters;i++){
		for (j=0; j<nClusters; j++)
			linksW[i][j] = 99999999.0;
	}

	for (i=0; i<nClusters; i++){
		for (j=0; j<adjMtrx[i].size (); j++){
			linksW[i][adjMtrx[i][j]] = getDistance(clusters[i][clusters[i].size ()-1], clusters[adjMtrx[i][j]][clusters[adjMtrx[i][j]].size ()-1]);
		}
		pairDist[i] = 999.0;
	}

	pair<short, short> *endsClusters;
	endsClusters = new pair<short, short> [ssEdges.size ()];
	for (k=0; k<ssEdges.size (); k++){
		for (i=0; i<nClusters; i++){
			for (j=0; j<clusters[i].size(); j++){
				if (getDistance(clusters[i][j], ssEnds[k].first) <= 0.05){
					endsClusters[k].first = i;
				}

				if (getDistance(clusters[i][j], ssEnds[k].second) <= 0.05){
					endsClusters[k].second = i;
				}
			}
		}
	}

	vector<vector<short> > traces (ssEdges.size ());		//to save traces for each stick

	for (i=0; i<ssEdges.size (); i++){

		//cout<<ssEnds[i].first .x<<" "<<ssEnds[i].first .y<<" "<<ssEnds[i].first .z<< " | "<<ssEnds[i].second .x<<" "<<ssEnds[i].second .y<<" "<<ssEnds[i].second .z<<endl;
		//find the shortest distance from the beginning of the stick to the end

		//cout<<"Stick# "<<i+1<<" Length= "<<getDistance(ssEdges[i][0],ssEdges[i][ssEdges[i].size ()-1]);

		//find the cluster that first end belongs to
		sourceCluster = endsClusters[i].first;
		getShortestTraces(nodes, adjMtrx, linksW, pairDist, sourceCluster, prevCluster, tracesLengths);

		short u = endsClusters[i].second;
		traces[i].push_back (u);
		while (prevCluster[u] != -1){						//if prevCluster = -1 , means no path found in between....
			traces[i].push_back (prevCluster[u]);
			u = prevCluster[u];
		}

		//cout<<"  traceFound "<<prevCluster[endsClusters[i].second]<<"  trace Length = "<<tracesLengths[endsClusters[i].second]<<endl;

		//if the trace found too long.....then ignore it
		if (tracesLengths[endsClusters[i].second] > 2 * getDistance(ssEdges[i][0],ssEdges[i][ssEdges[i].size ()-1])){
			j=1;
			while (j<traces[i].size ()-1){
				traces[i].erase (traces[i].begin () + j--);				//delete the trace..keep only the two clusters in both ends
				j++;
			}
		}
		//in case no path was found in between the two ends that represent one stick  (prevCluster = -1) or too long path
		if (prevCluster[endsClusters[i].second] == -1 || tracesLengths[endsClusters[i].second] > 2 * getDistance(ssEdges[i][0],ssEdges[i][ssEdges[i].size ()-1])){
			vector<Coordinate> firstClusterEnd = clusters[endsClusters[i].first],
								secondClusterEnd = clusters[endsClusters[i].second];

			clusters[endsClusters[i].first].clear ();
			clusters[endsClusters[i].second].clear ();

			deleteEdgeDensityNaive(clusters, ssEdges[i], ssEnds[i], 2.5);

			//restore clusters have the two ends
			clusters[endsClusters[i].first] = firstClusterEnd;
			clusters[endsClusters[i].second] = secondClusterEnd;
		}
	}

	//print paths
	/*
	for (i=0; i<ssEdges.size (); i++){
		cout<<" sticks# "<<i+1<<"  : "<<" ";
		for (j=0; j<traces[i].size (); j++)
			cout<<traces[i][j]+1<<" ";
		if (traces[i].size () == 1)
			cout<<endsClusters[i].first+1;
		cout<<endl;
	}
	*/


	//delete clusters represent sticks on the skeleton
	for (i=0; i<ssEdges.size (); i++){
		for (j=1; j<traces[i].size ()-1; j++){			//don't delete clusters contain the two ends
			clusters[traces[i][j]].clear();
		}
	}


	//re-calculate centroids of clusters
	Coordinate centroid;
	for (j=0; j<clusters.size (); j++){
		if (clusters[j].size ()>0)
			clusters[j].pop_back ();

		if (clusters[j].size ()>1){
			//re-calculate the centroid
			centroid.x = centroid.y = centroid.z = 0.0;
			for (k=0; k<clusters[j].size (); k++){
				centroid.x += clusters[j][k].x;
				centroid.y += clusters[j][k].y;
				centroid.z += clusters[j][k].z;
			}

			centroid.x /= clusters[j].size ();
			centroid.y /= clusters[j].size ();
			centroid.z /= clusters[j].size ();

			clusters[j].push_back (centroid);
		}
	}
	//cout<<"Done..."<<endl;

	//delete points from the clusters where ends belong to that are outside the stick
	float proj;
	Coordinate nPnt;
	//Coordinate centroid;
	for (i=0; i<ssEdges.size (); i++){
		if (traces[i].size () < 4){
			//cout<<"Stick# : "<<i+1<<endl;
			//cout<<"First End belongs to cluster : "<<endsClusters[i].first+1<<" # of points= "<<clusters[endsClusters[i].first].size ()<<endl;

			//delete the centroid
			if (clusters[endsClusters[i].first].size())
				clusters[endsClusters[i].first].pop_back ();

			j=0;
			while (j<clusters[endsClusters[i].first].size()){				//last point is the centroid

				//proj = linePointIntersection(ssEdges[i][0], ssEdges[i][ssEdges[i].size ()-1], clusters[endsClusters[i].first][j], nPnt);
				proj = round(linePointIntersection(ssEnds[i].first, ssEnds[i].second, clusters[endsClusters[i].first][j], nPnt),3);

				//cout<<"   projection of point : "<<j+1<<" = "<<proj<<endl;

				if (proj >=0 && proj <= 1.0){			//inside the stick
					clusters[endsClusters[i].first].erase (clusters[endsClusters[i].first].begin () + j--);
				}
				j++;
			}

			//re-add end point
			//clusters[endsClusters[i].first].push_back (ssEnds[i].first);

			//re-calculate the centroid
			if (clusters[endsClusters[i].first].size () > 1){
				centroid.x = centroid.y = centroid.z = 0.0;

				for (j=0; j<clusters[endsClusters[i].first].size(); j++){
					centroid.x += clusters[endsClusters[i].first][j].x;
					centroid.y += clusters[endsClusters[i].first][j].y;
					centroid.z += clusters[endsClusters[i].first][j].z;
				}
				centroid.x /= clusters[endsClusters[i].first].size ();
				centroid.y /= clusters[endsClusters[i].first].size ();
				centroid.z /= clusters[endsClusters[i].first].size ();

				clusters[endsClusters[i].first].push_back (centroid);
			}

			//cout<<"Second End belongs to cluster : "<<endsClusters[i].second+1<<" # of points = "<<clusters[endsClusters[i].second].size ()<<endl;

			//delete the centroid
			if (clusters[endsClusters[i].second].size())
				clusters[endsClusters[i].second].pop_back ();

			j=0;
			while (j<clusters[endsClusters[i].second].size()){				//last point is the centroid

				//proj = linePointIntersection(ssEdges[i][0], ssEdges[i][ssEdges[i].size ()-1], clusters[endsClusters[i].second][j], nPnt);
				proj = round(linePointIntersection(ssEnds[i].first, ssEnds[i].second, clusters[endsClusters[i].second][j], nPnt),3);

				//cout<<"   projection of point : "<<j+1<<" = "<<proj<<endl;

				if (proj >=0 && proj <= 1.0){			//inside the stick
					clusters[endsClusters[i].second].erase (clusters[endsClusters[i].second].begin () + j--);
				}
				j++;
			}


			//re-add end point
			//clusters[endsClusters[i].second].push_back (ssEnds[i].second);

			//update centroid
			if (clusters[endsClusters[i].second].size () > 1){
				centroid.x = centroid.y = centroid.z = 0.0;

				for (j=0; j<clusters[endsClusters[i].second].size(); j++){
					centroid.x += clusters[endsClusters[i].second][j].x;
					centroid.y += clusters[endsClusters[i].second][j].y;
					centroid.z += clusters[endsClusters[i].second][j].z;
				}
				centroid.x /= clusters[endsClusters[i].second].size ();
				centroid.y /= clusters[endsClusters[i].second].size ();
				centroid.z /= clusters[endsClusters[i].second].size ();

				clusters[endsClusters[i].second].push_back (centroid);
			}
		}

	}

	//delete empty clusters
	i=0;
	while (i<clusters.size ()){
		if (clusters[i].empty ()){
			clusters.erase (clusters.begin() + i--);
		}
		i++;
	}

	Protein densityPDB;

    /*
	//copy density voxels to a viewable PDB file
	Atom tmpAtom;
	AminoAcid tmpAA;
	int seqNum = 1;
	tmpAtom.name = " O  ";
	tmpAA.chr3 = "HOH";

	for (i=0; i<clusters.size (); i++){
		for (j=0; j<clusters[i].size (); j++){

			tmpAA.atoms.clear();
			tmpAA.num = seqNum;
			tmpAtom.coord = clusters[i][j];
			tmpAA.atoms.push_back(tmpAtom);
			seqNum++;
			densityPDB.AAs.push_back(tmpAA);

		}
	}

	//write local peaks into a pdb file
	densityPDB.writePDB(outPath + "_skeletonPoints_withoutSticks.pdb", 1, densityPDB.numOfAA());
    */
	//getchar();getchar();
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 *		Find all traces b/w any two SSends...then select the best trace fit the number of AAs in the loop b/w the 2 SSends according to the link in the graph
 */
void setLoopsWeights(vector<vector<cell> > & graph,
                     Map inMRC, vector<vector<Coordinate> > ssEdges,
                     vector<SecondaryStruct>	seqSS,
                     vector<Coordinate> pnts,
                     float peakTHRg,
                     Coordinate **(&traceList),
                     float *(&traceLength),
                     short *(&traceNoOfPoints),
                     string outPath){
	int i, j, k, l;

	float	distTHR = 15.0,					//thresold distance b/w two ends of SSs
			continuatyTHR = 1.1 * inMRC.apixX;			//a threshold used as a cuttoff distance b/w any two points to consider them continued (non-disjoint)

	vector<vector<Coordinate> > clusters, newClusters;
	short nClusters;								//the number of initial clusters (without clusters for SS ends)

	float **cDist;									//centroid to centroid distances

	vector<vector<short> >	adjMtrx,				//adjacency matrix for initial clusters
							clx,					//data structure where cliques will be stored (the indices of clusters form the clique)
							clxAdjMtrx;				//the adjacency matrix for the cliques
	vector<vector< Coordinate> > clxClusters;		//the actual coordinates of points form the clique

	vector<vector<short> > tracePaths;				//all trace paths for a given start point...also could be all trace path b/w two end points
	vector<float>			traceLengths;			//the approximated length of each trace path
	float maxLength = 0.0;							//the maximum possible length of the loop

	Protein trace;
	/*
	 *		find local peaks locally....b/w each pair of sticks
	 */
	/*
	//delete density around Sticks

	inMRC.cleanVxls (ssEdges, 3.5);

	short stk1Indx, stk2Indx, stk1start, stk2start;
	Coordinate sIndx, eIndx;
	vector<Coordinate> tmpPnt;
	for(i=1; i<graph.size ()-1; i++){
		for (j=0; j<graph[i].size (); j++){
			//determine first stick
			stk1Indx = j/2;
			if (j%2)
				stk1start = 0;
			else
				stk1start = ssEdges[stk1Indx].size ()-1;

			for (int k=0; k<graph[i][j].outLinks .size (); k++){

				if (graph[i][j].outLinks[k].rowIndx != graph.size ()-1 &&				//ignore links to the end node
					graph[i][j].nAAloop [k] <= 15){										//and work only on short loops

					pnts.clear ();
					cout<<"finding local peaks b/w ["<<i<<","<<j<<"] and ["<<graph[i][j].outLinks [k].rowIndx<<","<<graph[i][j].outLinks [k].colIndx<<"]. nAA= "<<graph[i][j].nAAloop [k];
					//determine second sticks
					stk2Indx = graph[i][j].outLinks [k].colIndx/2;
					if (graph[i][j].outLinks [k].colIndx %2)
						stk2start = ssEdges[stk2Indx].size ()-1;
					else
						stk2start = 0;

					//detrmine the two points on density map
					sIndx = ssEdges[stk1Indx][stk1start];
					eIndx = ssEdges[stk2Indx][stk2start];

					cout<<"  dist1= "<<getDistance(graph[i][j].cTerminal , graph[graph[i][j].outLinks [k].rowIndx][graph[i][j].outLinks [k].colIndx].nTerminal);

					//convert points from XYZ to MAP indexing system
					sIndx.x = (int) (sIndx.x/inMRC.apixX - inMRC.hdr.xorigin/inMRC.apixX + 0.5);
					sIndx.y = (int) (sIndx.y/inMRC.apixY - inMRC.hdr.yorigin/inMRC.apixY + 0.5);
					sIndx.z = (int) (sIndx.z/inMRC.apixZ - inMRC.hdr.zorigin/inMRC.apixZ + 0.5);

					eIndx.x = (int) (eIndx.x/inMRC.apixX - inMRC.hdr.xorigin/inMRC.apixX + 0.5);
					eIndx.y = (int) (eIndx.y/inMRC.apixY - inMRC.hdr.yorigin/inMRC.apixY + 0.5);
					eIndx.z = (int) (eIndx.z/inMRC.apixZ - inMRC.hdr.zorigin/inMRC.apixZ + 0.5);


					cout<<"  dist2= "<<getDistance(sIndx, eIndx)<<endl;
					//find local peaks
					inMRC.localPeaks(sIndx, eIndx, graph[i][j].nAAloop [k], pnts, peakTHRg);

					//write local peaks into a pdb file
					trace = points2pdb(pnts, "HOH", " O  ");
					if (trace.numOfAA()){
						trace.writePDB(outPath + "_localPeaks.pdb", 1, trace.numOfAA());
						cout<<"done..."<<endl;
					//	getchar();
					//	getchar();
					}

					//cluster local peaks
					peakClustering(pnts, clusters, ssEdges[stk1Indx][stk1start], 3.5); //start from a random end...which is here the first SS and the first end

					//allocate and calculate min distances and build adjacency matrix
					buildAdjMtrx(clusters, cDist, adjMtrx, inMRC.apixX, continuatyTHR);

					//	find all cliques in the graph....using Bron-Kerbosch algorithm
					findallCliques(clusters, adjMtrx, clx, clxAdjMtrx, clxClusters);

					//add the two ends of SSs as two clx....save everything in newClusters
					clusters.clear ();
					nClusters = clxClusters.size ();

					tmpPnt.clear ();
					tmpPnt.push_back (ssEdges[stk1Indx][stk1start]);
					clxClusters.push_back (tmpPnt);
					tmpPnt.clear ();
					tmpPnt.push_back (ssEdges[stk2Indx][stk2start]);
					clxClusters.push_back (tmpPnt);

					//print cliques
					for (l=0; l<clxClusters.size(); l++){
						Protein trace;

						//write local peaks into a pdb file
						trace = points2pdb(clxClusters[l], "HOH", " O  ");
						trace.writePDB(outPath +"_Clique_" + toString(l+1) +".pdb", 1, trace.numOfAA());
					}
					cout<<"done printing..."<<endl;
					//recalculate adjacency matrix for cliques
					buildAdjMtrx(clxClusters, cDist, clxAdjMtrx, inMRC.apixX, continuatyTHR);

					//find trace paths for first end point
					maxLength = graph[i][j].nAAloop [k] * 3.8;
					tracePaths.clear ();
					pathsCosts.clear ();
					findTracePaths(clxAdjMtrx, clxClusters, cDist, nClusters + stk1Indx, nClusters, tracePaths, pathsCosts, maxLength);


					//print paths
					Protein pathTrace;
					cout<<"Number of paths = "<<tracePaths.size ()<<endl;
					for (l=0; l<tracePaths.size (); l++){

						pathTrace.initialize();
						cout<<"Path# "<<l+1<<": ";
						for (int m=0; m<tracePaths[l].size (); m++){
							cout<<"  "<<tracePaths[l][m]+1;
							trace = points2pdb(clxClusters[tracePaths[l][m]], "HOH", " O  ");
							pathTrace.append(trace, 0, trace.numOfAA()-1);
						}
						cout<<" total length = "<<pathsCosts[l]<<endl;
						pathTrace.writePDB(outPath+"_pathTrace_"+toString(l+1)+".pdb", 1, pathTrace.numOfAA());
						getchar();
					}

					getchar();
					getchar();
				}
			}
		}

	}
	*/

	/*
	 *		Cluster peak points...start from a random point (one SS end)...this could be not accurate as if you
	 *		cluster points starting from the end point you r working on.....
	 */

	//save density indeces (in XYZ coordinate) into pnts data structure
	Coordinate tmppnt;
	for (int ix=0; ix<inMRC.numRows (); ix++){
		tmppnt.x = round(ix*inMRC.apixX + inMRC.hdr.xorigin,3);
		for (int iy=0;iy<inMRC.numCols (); iy++){
			tmppnt.y = round(iy*inMRC.apixY + inMRC.hdr.yorigin,3);
			for (int iz=0; iz<inMRC.numSlcs (); iz++){
				if (inMRC.cube [ix][iy][iz] > 0.0){
					tmppnt.z = round(iz*inMRC.apixZ + inMRC.hdr.zorigin,3);

					pnts.push_back (tmppnt);
				}
			}
		}
	}

	//write local peaks into a pdb file
	//trace = points2pdb(pnts, "HOH", " O  ");
	//trace.writePDB(outPath + "_localPeaks.pdb", 1, trace.numOfAA());


	peakClustering(pnts, clusters, ssEdges[0][0], 2.0*inMRC.apixX); //start from a random end...which is here the first SS and the first end

	nClusters = clusters.size ();

	/*
	//print clusters (all points or centriod only)
	for (i=0; i<nClusters; i++){
		Protein trace;
		//write local peaks into a pdb file
		trace = points2pdb(clusters[i], "HOH", " O  ");
		trace.writePDB(outPath +"_Cluster_" + toString(i+1) +".pdb", 1, trace.numOfAA());							//all points
		//trace.writePDB(outPath +"_Cluster_" + toString(i+1) +".pdb", trace.numOfAA(), trace.numOfAA());				//centroid only
	}
	*/

	/*
	 *				Delete Density around edges sticks
	 */

	pair<Coordinate, Coordinate> *edgesEnds;
	edgesEnds = new pair<Coordinate, Coordinate> [ssEdges.size ()];

	deleteEdgesDensity(inMRC, clusters, ssEdges, edgesEnds, outPath);


	nClusters = clusters.size ();

	/*
	//print clusters (all points or centriod only)
	for (i=0; i<nClusters; i++){
		Protein trace;
		//write local peaks into a pdb file
		trace = points2pdb(clusters[i], "HOH", " O  ");
		trace.writePDB(outPath +"_Cluster_" + toString(i+1) +".pdb", 1, trace.numOfAA());							//all points
		//trace.writePDB(outPath +"_Cluster_" + toString(i+1) +".pdb", trace.numOfAA(), trace.numOfAA());				//centroid only
	}
	*/
	/*
	 *		find shortest distances b/w each pair of clusters and then build adjMtrx
	 */
	buildAdjMtrx(clusters, cDist, adjMtrx, inMRC.apixX, nClusters, continuatyTHR);


#ifdef _WIN				//if work under Windows
	clock_t start, finish;
	float tTime;
	start = clock();
/*
#else
	float tTime;
	struct timeval ti_start, ti_end;
	double time_di;
	cout<<"Time(0)= "<<time(0)<<endl;
	//start timer
	gettimeofday(&ti_start,0);
*/
#endif

	/*
	 *	find all cliques in the graph....using Bron-Kerbosch algorithm
	 */

	//cout<<"Finding Cliques....";
	findallCliques(clusters, adjMtrx, clx, clxAdjMtrx, clxClusters);
	//cout<<" Done"<<endl;

	//print final cliques and their associate local peaks clusters
	/*
	for (i=0; i<clx.size (); i++){
		cout<<"Clique# "<<i+1<<"  : ";
		for (j=0; j<clx[i].size(); j++){
			cout<<clx[i][j]+1<<" ";
		}
		cout<<endl;
	}
	*/

	//for SS ends....create a cluster of length = 1 contains one end point for each SS
	newClusters = clxClusters;
	clusters.clear ();
	vector<Coordinate> tmpCluster;
	for (i=0; i<ssEdges.size (); i++){

		tmpCluster.clear ();
		//tmpCluster.push_back (ssEdges[i][0]);
		tmpCluster.push_back (edgesEnds[i].first);
		newClusters.push_back (tmpCluster);


		tmpCluster.clear ();
		//tmpCluster.push_back (ssEdges[i][ssEdges[i].size ()-1]);
		tmpCluster.push_back (edgesEnds[i].second);
		newClusters.push_back (tmpCluster);
	}

	//write cliques (all Points all centroid only)
    /*
	for (i=0; i<newClusters.size(); i++){
		Protein trace;

		//write local peaks into a pdb file
		trace = points2pdb(newClusters[i], "HOH", " O  ");
		trace.writePDB(outPath +"_Clique_" + toString(i+1) +".pdb", 1, trace.numOfAA());						//all points
		//trace.writePDB(outPath +"_Clique_" + toString(i+1) +".pdb", trace.numOfAA(), trace.numOfAA());			//centroid only
	}
	*/


//	cout<<"Press any key.."<<endl;
//	getchar();
//	getchar();

	nClusters = clx.size ();			//old number of cliques before adding the end points


	//re-calculate Adjacency Mtrx and short distances matrix for cliques
	buildAdjMtrx(newClusters, cDist, clxAdjMtrx, inMRC.apixX, nClusters, continuatyTHR);

	//print Adjacency Matrix
	/*
	for (i=0; i<clxAdjMtrx.size (); i++){
		cout<<"Clique# "<<i+1<<"  : ";
		for (j=0; j< clxAdjMtrx[i].size (); j++){
			cout<<" "<<clxAdjMtrx[i][j]+1;
		}
		cout<<endl;
	}
	*/

	//getchar();
	//getchar();

	/*
	 *		find trace paths from density map "Skeleton" or "Local Peaks"
	 */
	Protein pathTrace;
	tracePaths.clear ();
	traceLengths.clear ();
	for (i=0; i<2*ssEdges.size (); i++){

		//determine the maximum length of a loop could get out from this end
		maxLength = 0.0;
		for (j=0; j<graph.size (); j++){
			for (k=0; k<graph[j][i].nAAloop .size (); k++){
				if (graph[j][i].nAAloop[k] > maxLength)
					maxLength = graph[j][i].nAAloop[k];
			}
		}

		//find trace paths for first end point
		findTracePaths(clxAdjMtrx, cDist, nClusters + i, nClusters, tracePaths, traceLengths, maxLength*3.8);
	}

	//cout<<"Number of Total Paths : "<<tracePaths.size()<<" and "<<traceLengths.size()<<endl;

	//print paths

	//print out cliques not considered in all paths
	/*
	vector<bool>  considered (newClusters.size (), false);
	for (i=0; i<tracePaths.size (); i++){
		for (j=0; j<tracePaths[i].size (); j++){
			considered[tracePaths[i][j]] = true;
		}
	}
	cout<<"Cliques were not considered in all traces are "<<endl;
	for (i=0; i<considered.size (); i++){
		if (!considered[i])
			cout<<"clique# "<<i+1<<endl;
	}
	*/

	//add traces discontinue (have gap) close to any one of SS end
	/*
	vector<short> newTrace;
	vector<vector<short> > newTraceSet;
	float newLength;
	vector<float> newTraceLengths;
	float ssGapTHR = 8.0;						//the maximum distance b/w the SS and the first cluster in the path to consider
	float gapDist = 9999.0;
	for (i=0; i<graph[0].size (); i++){

		//find how many traces out from a particular SS end
		int tracesCntr = 0;
		for (j=0; j<tracePaths.size (); j++){
			if (tracePaths[j][0] == i+nClusters && tracePaths[j].size() == 1)
				tracesCntr++;
		}

		if (tracesCntr == 1){				//only consider SS ends without any traces out from them...
			j=0;
			while (j<tracePaths.size ()){
				if (tracePaths[j][tracePaths[j].size ()-1] < nClusters){
					gapDist = getDistance(newClusters[tracePaths[j][tracePaths[j].size()-1]][newClusters[tracePaths[j][tracePaths[j].size()-1]].size()-1], newClusters[i+nClusters][newClusters[i+nClusters].size ()-1]);
					if (gapDist < ssGapTHR){
						//push the new path into tracePaths
						//first find where this SS starts in trace Paths
						newTrace = tracePaths[j];
						newLength = traceLengths[j];
						newTrace.push_back(i+nClusters);
						newLength += gapDist;

						newTraceSet.push_back(newTrace);
						newTraceLengths.push_back(newLength);

						//cout<<"I found a new path b/w graph col "<<i<<" and path# "<<j+1<< " old length = "<<traceLengths[j]<<" new Length= "<<newLength<<" gap= "<<gapDist<<endl;
					}
				}
				j++;
			}
		}
		//getchar();getchar();
	}
	cout<<"Number of traces found with gap next to SS ends : "<<newTraceSet.size ()<<endl;

	//save new traces found
	for (i=0; i<newTraceSet.size(); i++){
		for (k=0; k<tracePaths.size (); k++){
			if (tracePaths[k][0] == newTraceSet[i][0]){
				tracePaths.insert(tracePaths.begin() + k+1, newTraceSet[i]);
				traceLengths.insert(traceLengths.begin() + k + 1, newTraceLengths[i]);
				break;
			}
		}
		//push in the reverse direction
		newTrace.clear ();
		for (j=newTraceSet[i].size ()-1; j>-1; j--){
			newTrace.push_back(newTraceSet[i][j]);
		}

		for (k=0; k<tracePaths.size (); k++){
			if (tracePaths[k][0] == newTrace[0]){
				tracePaths.insert(tracePaths.begin() + k+1, newTrace);
				traceLengths.insert(traceLengths.begin() + k + 1, newTraceLengths[i]);
				break;
			}
		}
	}

	newTrace.clear();
	newTraceSet.clear();
	newTraceLengths.clear();
	*/

#ifdef _WIN
	finish = clock();
	cout<<"Time taken to find all traces "<<finish-start<<" ms."<<endl;
	tTime = finish-start;
/*
#else
	//stop timer
	gettimeofday(&ti_end,0);
	time_di = (ti_end.tv_sec-ti_start.tv_sec)*1000000 + ti_end.tv_usec - ti_start.tv_usec;
	cout<<"Time(0)= "<<time(0)<<endl;
	cout<<"Time taken to enumerate top paths is "<<(double) (time_di/1000)<<" ms "<<endl;

	tTime = (float) (time_di/1000);
*/
#endif

	//getchar();getchar();

	//delete upnormal paths that form sharp curves
	int nPadPaths = 0;
	float angle;
	bool pad = false;
	for (i=0; i<tracePaths.size (); i++){
		pad = false;


		if (traceLengths[i] != -1){

			//Check for sharp curves
			if (tracePaths[i].size () > 2){
				for (j=0; j<tracePaths[i].size ()-2; j++){
					angle = getAngleDegree(newClusters[tracePaths[i][j]][newClusters[tracePaths[i][j]].size ()-1], newClusters[tracePaths[i][j+1]][newClusters[tracePaths[i][j+1]].size()-1], newClusters[tracePaths[i][j+2]][newClusters[tracePaths[i][j+2]].size()-1]);
					if ( angle < 30.0){
						//cout<<"PAD#"<<i+1<<" ("<<angle<<"): ";
						//for (k=0; k<tracePaths[i].size (); k++){
						//	cout<<"  "<<tracePaths[i][k]+1;
						//}
						//tracePaths[i].clear ();

						traceLengths[i] = -1;
						nPadPaths++;
						pad = true;
						break;
					}
				}
			}


			/*
			//check for redundancy
			if (!pad){
				j=i+1;
				while (j < tracePaths.size () && tracePaths[j][0] == tracePaths[i][0]){
					if (tracePaths[i][tracePaths[i].size ()-1] == tracePaths[j][tracePaths[j].size ()-1]){			//same destination
						short numDiff = 0;		//num of different clusters in the trace
						for (int icur = 1; icur<tracePaths[i].size ()-1; icur++){
							for (int jcur=1; jcur<tracePaths[j].size ()-1; jcur++){
								if (tracePaths[i][icur] != tracePaths[j][jcur])
									numDiff++;
							}
						}
						if (tracePaths[i].size () > tracePaths[j].size ()){
							if ((float) numDiff/tracePaths[i].size () < 0.2){
								pad = true;
								//tracePaths[i].clear ();
								traceLengths[i] = -1;
								nPadPaths++;
								break;
							}
						}
						else{
							if ((float) numDiff/tracePaths[j].size () < 0.2){
								pad = true;
								//tracePaths[j].clear ();
								traceLengths[j] = -1;
								nPadPaths++;
								break;
							}
						}

					}
					j++;
				}
			}
			*/

		}

		//if (!pad){
		//	cout<<"GOOD#"<<i+1<<" : ";
		//	for (k=0; k<tracePaths[i].size (); k++){
		//		cout<<"  "<<tracePaths[i][k]+1;
		//	}
		//}
		//cout<<endl;
	}
	//cout<<"# of Pad Paths found = "<<nPadPaths<<endl;

	//copy tracePaths into Array structure to speed up the process
	short **tracePathsArr,
			*tracePathsSize;
	float *traceLengthsArr;
	int nTracePaths = tracePaths.size () - nPadPaths,
		nxtCntr=0;
	tracePathsArr = new short *[nTracePaths];
	traceLengthsArr = new float [nTracePaths];
	tracePathsSize = new short [nTracePaths];
	for (i=0; i<tracePaths.size (); i++){
		if (traceLengths[i] > -1){
			traceLengthsArr[nxtCntr] = traceLengths[i];
			tracePathsSize[nxtCntr] = tracePaths[i].size ()-1;
			tracePathsArr[nxtCntr] = new short [tracePaths[i].size ()];
			for (j=0; j<tracePaths[i].size (); j++){
				tracePathsArr[nxtCntr][j] = tracePaths[i][j];
			}
			nxtCntr++;
		}
	}
	tracePaths.clear ();
	traceLengths.clear ();

	//cout<<"Traces copied to new data structure..."<<endl;
	/*
	//write paths into viewable files by Chimera
	for (j=0; j<nTracePaths; j++){
		pathTrace.initialize();
		pathTrace.header.push_back("seq. of clusters in the path : ");
		for (k=0; k<=tracePathsSize[j]; k++){
			//cout<<"  "<<tracePaths[j][k]+1;
			trace = points2pdb(newClusters[tracePathsArr[j][k]], "GLY", " CA ");
			pathTrace.append(trace, trace.numOfAA()-1, trace.numOfAA()-1);
			pathTrace.header[0] += toString(tracePathsArr[j][k]+1);
			pathTrace.header[0] += " ";
		}
		pathTrace.header[0] += "  eLength= ";
		pathTrace.header[0] += toString(traceLengthsArr[j]);


		//cout<<" length= "<<traceLengths[j]<<endl;
		pathTrace.writePDB(outPath+"_pathTrace_"+toString(tracePathsArr[j][0]-nClusters)+"w"+toString(tracePathsArr[j][tracePathsSize[j]]-nClusters)+"_"+toString(j+1)+".pdb", 1, pathTrace.numOfAA(), 1);
	}
	*/

#ifdef _WIN
	start = clock();
/*
#else
	//start timer
	gettimeofday(&ti_start,0);
*/
#endif


	//set weights on the graph
	int bestTraceIndx1,bestTraceIndx2;
	float loopLength = 0;			//the length of the loop b/w the two SS
	for (i=1; i<graph.size ()-1; i++){
		for (j=0; j<graph[i].size (); j++){
			k=0;
			while (k<graph[i][j].outLinks .size ()){
			//for (k=0; k<graph[i][j].outLinks .size (); k++){
				if (graph[i][j].outLinks[k].rowIndx != graph.size ()-1){
					bestTraceIndx1 = -1;
					bestTraceIndx2 = -1;

					//cout<<"graph ["<<i<<" "<<j<<"] with ["<<graph[i][j].outLinks[k].rowIndx <<" "<<graph[i][j].outLinks[k].colIndx <<"]"<<endl;

					// calculate the length (in Angstrom) of the loop b/w two SS ends
					loopLength = 3.8 * graph[i][j].nAAloop [k];


					if (graph[i][j].outLinks [k].rowIndx != i+1){
						//calculate how many AA in the loop is a hlx AA
						for (l=i+1; l<graph[i][j].outLinks [k].rowIndx; l++){
							if (seqSS[l-1].type == 'H'){				//l-1 b/s i is shifted by 1 (first row is for start node)
								loopLength = loopLength - seqSS[l-1].nAA *3.8 + seqSS[l-1].nAA * ALPHA_RISE;
							}
						}
					}

					if (j%2==0)
						graph[i][j].outLinks [k].w = getBestPathTrace(newClusters, tracePathsArr, nTracePaths, tracePathsSize, traceLengthsArr, nClusters, loopLength, j+nClusters + 1, graph[i][j].outLinks[k].colIndx + nClusters, bestTraceIndx1, bestTraceIndx2);
					else
						graph[i][j].outLinks [k].w = getBestPathTrace(newClusters, tracePathsArr, nTracePaths, tracePathsSize, traceLengthsArr, nClusters, loopLength, j+nClusters - 1, graph[i][j].outLinks[k].colIndx + nClusters, bestTraceIndx1, bestTraceIndx2);

					//delete un-realistic links
					if (bestTraceIndx1 == -1 && graph[i][j].nAAloop[k] > 35){
						graph[i][j].outLinks.erase (graph[i][j].outLinks .begin () + k);
						graph[i][j].nAAloop .erase (graph[i][j].nAAloop .begin () + k);
						graph[i][j].traceNum .erase (graph[i][j].traceNum .begin () + k);
						k--;
					}
					else{

						graph[i][j].traceNum [k] = toString(bestTraceIndx1+1);
						graph[i][j].traceNum [k] += "-";
						graph[i][j].traceNum [k] += toString(bestTraceIndx2+1);

						//update inLink in the other node
						/*
						for (l=0; l<graph[graph[i][j].outLinks [k].rowIndx][graph[i][j].outLinks [k].colIndx].inLinks.size (); l++){
							if (graph[graph[i][j].outLinks [k].rowIndx][graph[i][j].outLinks [k].colIndx].inLinks[l].rowIndx == i){
								if (graph[graph[i][j].outLinks [k].rowIndx][graph[i][j].outLinks [k].colIndx].inLinks[l].colIndx == j){
									graph[graph[i][j].outLinks [k].rowIndx][graph[i][j].outLinks [k].colIndx].inLinks[l].w = graph[i][j].outLinks [k].w;
									break;
								}
							}
						}
						*/
					}

					//if (i==5 && j==14 && graph[i][j].outLinks[k].rowIndx == 7 && graph[i][j].outLinks[k].colIndx == 19) {getchar(); getchar();}
				}
				k++;
			}
			//getchar();getchar();

		}
	}

#ifdef _WIN
	finish = clock();
	cout<<"Time taken to re-weight links "<<finish-start<<" ms."<<endl;
	tTime += finish-start;
/*
#else
	//stop timer
	gettimeofday(&ti_end,0);
	time_di = (ti_end.tv_sec-ti_start.tv_sec)*1000000 + ti_end.tv_usec - ti_start.tv_usec;
	cout<<"Time(0)= "<<time(0)<<endl;
	cout<<"Time taken to re-weight links "<<(double) (time_di/1000)<<" ms "<<endl;

	tTime += (float) (time_di/1000);
*/
#endif

	//save traces found on the map
	traceList   = new Coordinate *[nTracePaths+1];    //traces  coordinates
	traceLength = new float [nTracePaths+1];          //lengths of traces
	traceNoOfPoints = new short [nTracePaths+1];          //number of points in each trace

	traceNoOfPoints[0] = nTracePaths+1;       //save number of traces

	for (i=1; i<nTracePaths+1; i++){
	    traceList[i] = new Coordinate [tracePathsSize[i-1]+1];
	    for (j=0; j<=tracePathsSize[i-1];j++){
	        traceList[i][j] = newClusters[tracePathsArr[i-1][j]][newClusters[tracePathsArr[i-1][j]].size()-1];        //get the centroid coordinate of that cluster
	    }
	    traceLength[i] = traceLengthsArr[i-1];
	    traceNoOfPoints[i]  = tracePathsSize[i-1]+1;
	}

    //delete old traces containers
    delete [] traceLengthsArr;
    delete [] *tracePathsArr;
    delete [] tracePathsArr;
    delete [] tracePathsSize;

	clusters.clear();
	newClusters.clear();
    delete [] *cDist;
    delete [] cDist;

	adjMtrx.clear();
    clx.clear();
    clxAdjMtrx.clear();
	clxClusters.clear();

	tracePaths.clear();
	traceLengths.clear();

	delete [] edgesEnds;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void fillSkeletonGaps(Map &inSkeleton, vector<vector<Coordinate> > ssEdges, string outDir){

	vector<Coordinate> skeletonPnts;
	vector<vector<Coordinate> > clusters;
	int nClusters, i, j, ix, iy, iz;
	float continuatyTHR = 1.1 ;//* inSkeleton.apixX;
	float **cDist;
	vector<vector<short> > adjMtrx;

	//save density indeces (in XYZ coordinate) into pnts data structure
	Coordinate tmppnt;
	for (ix=0; ix<inSkeleton.numRows (); ix++){
		tmppnt.x = ix;
		//tmppnt.x = ix*inMRC.apixX + inMRC.hdr.xorigin;
		for (iy=0;iy<inSkeleton.numCols (); iy++){
			tmppnt.y = iy;
			//tmppnt.y = iy*inMRC.apixY + inMRC.hdr.yorigin;
			for (iz=0; iz<inSkeleton.numSlcs (); iz++){
				if (inSkeleton.cube [ix][iy][iz] > 0.0){
					tmppnt.z = iz;
					//tmppnt.z = iz*inMRC.apixZ + inMRC.hdr.zorigin;

					skeletonPnts.push_back (tmppnt);
				}
			}
		}
	}

	//start point in indeces system
	Coordinate startPnt = ssEdges[0][0];
	startPnt.x = (startPnt.x - inSkeleton.hdr.xorigin)/inSkeleton.apixX;
	startPnt.y = (startPnt.y - inSkeleton.hdr.yorigin)/inSkeleton.apixY;
	startPnt.z = (startPnt.z - inSkeleton.hdr.zorigin)/inSkeleton.apixZ;

	//get clusters
	peakClustering(skeletonPnts, clusters, startPnt, inSkeleton.apixX); //start from a random end...which is here the first SS and the first end

	nClusters = clusters.size ();

	/*
	 *		find shortest distances b/w each pair of clusters and then build adjMtrx
	 */
	buildAdjMtrx(clusters, cDist, adjMtrx, inSkeleton.apixX, nClusters, continuatyTHR);


	float dist = 999.0;
	vector<int> iClusters;					//target clusters that we wanna re-connect with other clusters in the map....
	vector<pair<short, short> > nListRange;		//the starting (first) and Ending (second) indeces of iCluster[x] neighbors in clusterDistances

	//onePairInfo = new float [4];

	//find clusters with less number of neighbors
	for (i=0; i<nClusters; i++){
		if (adjMtrx[i].size ()<2){
			pair<short, short> tmpIndeces;

			iClusters.push_back (i);
			tmpIndeces.first = tmpIndeces.second = -1;
			nListRange.push_back (tmpIndeces);
		}
	}

	bool cont = true;
	while (cont){

		cout<<"================== one round ================="<<endl;
		cont = false;

		vector<vector<float> >  clusterDistances;		//[0] is isolated cluster  [1] the neighbor cluster  [2] trace length  [3] pairwise distance
		vector<float> onePairInfo (4, -1);

		for (i=0; i<iClusters.size (); i++){

			onePairInfo[0] = iClusters[i];
			nListRange[i].first = clusterDistances.size ();

			//cout<<endl<<onePairInfo[0]+1<<endl;

			//search for close clusters to connect them
			for (j=0; j<nClusters; j++){
				float minDist = 9999.9;
				if (iClusters[i]!=j){
					for (int k=0; k<clusters[iClusters[i]].size (); k++){
						for (int l=0; l<clusters[j].size (); l++){
							dist = getDistance(clusters[iClusters[i]][k], clusters[j][l]) * inSkeleton.apixX;
							if (dist <= minDist){
								minDist = dist;

								//k = clusters[iClusters[i]].size ();				//to break the outer loop
								//break;
							}
						}
					}
					if (minDist <= 15.0){
						onePairInfo[1] = j;
						onePairInfo[3] = dist;
						clusterDistances.push_back (onePairInfo);
					}
				}
			}
			nListRange[i].second = clusterDistances.size ()-1;
		}

		float *pairDist;							//the maximum length of a path allowed b.w any two points
		float *tracesLengths;						//the cost to reach each cluster
		short *prevCluster;							//the prev cluster in the trace

		short *nodes;
		float **linksW;

		pairDist = new float [nClusters];
		prevCluster = new short [nClusters];
		tracesLengths = new float [nClusters];

		//initialization links weights
		nodes = new short[nClusters];
		linksW = AllocateDynamicArray<float> (nClusters, nClusters);			//create links weights
		for (i=0; i<nClusters;i++)
			for (j=0; j<nClusters; j++)
				linksW[i][j] = 99999999.0;

		for (i=0; i<nClusters; i++){
			for (j=0; j<adjMtrx[i].size (); j++){
				linksW[i][adjMtrx[i][j]] = getDistance(clusters[i][clusters[i].size ()-1], clusters[adjMtrx[i][j]][clusters[adjMtrx[i][j]].size ()-1]);
			}
			pairDist[i] = 999.0;
		}

		for (i=0; i<iClusters.size (); i++){
			//find the cluster that first end belongs to
			int sourceCluster = iClusters[i];

			//cout<<"Cluster# "<<sourceCluster+1<<" should be connected to a close cluster [ "<<nListRange[i].first <<" "<<nListRange[i].second<<"  ]"<<endl;

			getShortestTraces(nodes, adjMtrx, linksW, pairDist, sourceCluster, prevCluster, tracesLengths);

			for (j=nListRange[i].first; j<=nListRange[i].second; j++){
				clusterDistances[j][2] = tracesLengths[(int) clusterDistances[j][1]];
			//	cout<<"   "<<clusterDistances[j][1]+1<<"  traceDist= "<<clusterDistances[j][2]<<"   straightDist= "<<clusterDistances[j][3]<<endl;
			}
		}

		//sort information according to straight distances
		for (i=0; i<clusterDistances.size (); i++){
			for (j=i+1; j<clusterDistances.size (); j++){
				if (clusterDistances[i][3] > clusterDistances[j][3]){
					//swap
					onePairInfo = clusterDistances[i];
					clusterDistances[i] = clusterDistances[j];
					clusterDistances[j] = onePairInfo;
				}
			}
		}

		for (i=0; i<clusterDistances.size (); i++){
			if (clusterDistances[i][2] >= 30.0){
				//add the the cluster to the neighbors of cluster clusterDistances[i][0]
				adjMtrx[clusterDistances[i][0]].push_back(clusterDistances[i][1]);
				adjMtrx[clusterDistances[i][1]].push_back(clusterDistances[i][0]);

				//fill density map b/w the two clusters
				//get the index of the two cluster on the density cube
				Coordinate indx1 = clusters[clusterDistances[i][0]][clusters[clusterDistances[i][0]].size ()-1],
							indx2 = clusters[clusterDistances[i][1]][clusters[clusterDistances[i][1]].size ()-1];

				//cout<<"indx1.x = "<<(int)indx1.x<<" y= "<<(int)indx1.y<<"  z= "<<(int)indx1.z<<"  indx2.x= "<<(int)indx2.x<<" y= "<<(int)indx2.y<<"  z= "<<(int)indx2.z<<endl;

				int xDelta = (int) indx2.x - (int) indx1.x;
				int yDelta = (int) indx2.y - (int) indx1.y;
				int zDelta = (int) indx2.z - (int) indx1.z;

				//cout<<"xDelta= "<<xDelta<<"  yDelta= "<<yDelta<<"  zDelta= "<<zDelta<<endl;

				int xSign = 1;
				int ySign = 1;
				int zSign = 1;

				if (xDelta <0)
					xSign = -1;
				if (yDelta < 0)
					ySign = -1;
				if (zDelta <0)
					zSign = -1;

				//cout<<"xSign= "<<xSign<<"  ySign= "<<ySign<<" zSign= "<<zSign<<endl; getchar();getchar();

				ix = (int) indx1.x;
				iy = (int) indx1.y;
				iz = (int) indx1.z;

				while (ix != (int) indx2.x){

					//cout<<"ix= "<<ix<<" iy= "<<iy<<" iz= "<<iz<<endl;

					inSkeleton.cube[ix][iy][iz] = 1.0;		//fill density
					ix += xSign;
				}
				//cout<<endl;

				//iy = (int) indx1.y + ySign;
				while (iy != (int) indx2.y){
					//cout<<"ix= "<<ix<<" iy= "<<iy<<" iz= "<<iz<<endl;
					inSkeleton.cube[ix][iy][iz] = 1.0;		//fill density
					iy += ySign;

				}

				//cout<<endl;
				//iz = (int) indx1.z + zSign;
				while (iz != (int) indx2.z){
					//cout<<"ix= "<<ix<<" iy= "<<iy<<" iz= "<<iz<<endl;
					inSkeleton.cube[ix][iy][iz] = 1.0;		//fill density
					iz += zSign;
				}

				//cout<<"    "<<clusterDistances[i][0]+1<<"  connected with "<<clusterDistances[i][1]+1<<endl;
				cont = true;
				break;

			}
		}
		//print information
		//for (i=0; i<clusterDistances.size (); i++){
		//	cout<<"Cluster1= "<<clusterDistances[i][0]+1<<"  cluster2= "<<clusterDistances[i][1]+1<<"  tDist= "<<clusterDistances[i][2]<<"  sDist= "<<clusterDistances[i][3]<<endl;
		//}
	}

	for (i=0; i<adjMtrx.size (); i++){
		cout<<"cluster# "<<i+1<<" : ";
		for (j=0; j<adjMtrx[i].size (); j++){
			cout<<adjMtrx[i][j]+1<<" ";
		}
		cout<<endl;
	}



	//write out the density after filling
	inSkeleton.write(outDir + "_outDensity.mrc");
	/*
	//print clusters (all points or centriod only)
	for (i=0; i<nClusters; i++){
		Protein trace;
		//convert cluster points to 3d XYZ coordinate system
		for (j=0; j<clusters[i].size (); j++){
			clusters[i][j].x = clusters[i][j].x*inSkeleton.apixX + inSkeleton.hdr.xorigin;
			clusters[i][j].y = clusters[i][j].y*inSkeleton.apixY + inSkeleton.hdr.yorigin;
			clusters[i][j].z = clusters[i][j].z*inSkeleton.apixZ + inSkeleton.hdr.zorigin;
		}
		//write local peaks into a pdb file
		trace = points2pdb(clusters[i], "HOH", " O  ");
		trace.writePDB(outDir +"_Cluster_" + toString(i+1) +".pdb", 1, trace.numOfAA());								//all points
		//trace.writePDB(outDir +"_Cluster_" + toString(i+1) +".pdb", trace.numOfAA(), trace.numOfAA());				//centroid only
	}
	*/
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Dijkstra algorithm to find the shortest trace b/w a SS end and all other SS ends
inline
void getShortestTraces(short *nodes, vector<vector<short> >	&adjMtrx, float **linksW, float *maxLength, int source, short *(&prev), float *(&D)){

	/*	Dijkstra Algorithm (graph, source)
		 2      for each vertex v in Graph:           // Initializations
		 3          dist[v] := infinity ;              // Unknown distance function from source to v
		 4          previous[v] := undefined ;         // Previous node in optimal path from source
		 5      end for ;
		 6      dist[source] := 0 ;                    // Distance from source to source
		 7      Q := the set of all nodes in Graph ;
				// All nodes in the graph are unoptimized - thus are in Q
		 8      while Q is not empty:                 // The main loop
		 9          u := vertex in Q with smallest dist[] ;
		10          if dist[u] = infinity:
		11              break ;                        // all remaining vertices are inaccessible from source
		12          fi ;
		13          remove u from Q ;
		14          for each neighbor v of u:         // where v has not yet been removed from Q.
		15              alt := dist[u] + dist_between(u, v) ;
		16              if alt < dist[v]:             // Relax (u,v,a)
		17                  dist[v] := alt ;
		18                  previous[v] := u ;
		19              fi  ;
		20          end for ;
		21      end while ;
		22      return dist[] ;

	*/


	//variables
	short i,j;
    int nNodes = adjMtrx.size ();


	for (i=0; i<nNodes; i++){
		D[i] = 999999999.0;					//intial distance b/w source any any other node is undefined (-1)
		prev[i] = -1;						//no previous node at the begining
		nodes[i] = i;						//indeces of nodes
	}
	//nodes[source] = -1;					//The prev of source in the trace is undefined

	D[source] = 0;							//dist from source to source is undefined
	int nNodesSolved = 0;
	while (nNodesSolved < nNodes){
		//cout<<"iteration "<<nNodesSolved+1<<endl;
		//find the closest node to source
		int minValue = 99999999.0;
		int minNode = -1;			//int minNode = start node;
		for (i = 0; i < nNodes; i++)
		{
			if (nodes[i] == -1)
				continue;
			if (D[i] < minValue)
			{
				minValue = D[i];
				minNode = i;
			}
		}

		if (minNode == -1 || D[minNode] == 999999999.0)			// all remaining vertices are inaccessible from source
			break;

		//cout<<"   Cluster w min dist= "<<minNode+1<<endl;
		nodes[minNode] = -1;			//remove the closest node (minNode) from the list of vertices
		for (i=0; i<adjMtrx[minNode].size (); i++){
			float nDist = D[minNode] + linksW[minNode][adjMtrx[minNode][i]];
			if (nDist < D[adjMtrx[minNode][i]] && nDist < maxLength[adjMtrx[minNode][i]]){
				D[adjMtrx[minNode][i]] = nDist;
				prev[adjMtrx[minNode][i]] = minNode;
			}
		}
		nNodesSolved++;
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 *		find the shortest trace b/w 2 SSends. the only condition is to fit the number of AA in the loop b/w these 2 SSends according to the graph link
 */
void setLoopsWeights_ShortestTrace(vector<vector<cell> > & graph, Map inMRC, vector<vector<Coordinate> > ssEdges, vector<SecondaryStruct>	seqSS, vector<Coordinate> pnts, float peakTHRg, string outPath){

	float **cDist;									//centroid Distance clusters (indeces and distance) ... each cluster will represent a node in a network later
	float continuatyTHR = 1.1 * inMRC.apixX;		//a threshold used as a cuttoff distance b/w any two points in any two clusters to consider them continued (non-disjoint)
	vector<vector<short> >	adjMtrx;				//adjacency matrix for clusters

	vector<vector< Coordinate> > clustersCoord;		//the actual coordinates of points form the clique
	short nClustersAll;								//number of all clusters (include SSEnds as clusters)

	Protein trace;									// a temp pdb file to write clusters into pdb files viewable by Chimera

	int i, j, k, l, m,
		nClusters;									//number of clusters found for local peaks

	short sourceCluster;							//the endex of source cluster cluster (node) for Dijkstra

	/*
	 *		Clustering peak points...start from a random point (one SS end)...this could be not accurate as if you
	 *		cluster points starting from the end point you r working on.....
	 */

	/*
	 *		if work on skeletonization.....comment this
	 */
	/*
		//find local peaks globally
	inMRC.localPeaksMap(pnts, ssEdges, 3.5, peakTHRg);

	*/

	/*
	 *		if do not work on skeletonization...comment this
	 */
	//save density indeces (in XYZ coordinate) into pnts data structure
	//*

	//remove density around sticks		... if you removed density around sticks manually...comment this statement
	//inMRC.cleanVxls (ssEdges, 3.5);

	Coordinate tmppnt;
	for (int ix=0; ix<inMRC.numRows (); ix++){
		for (int iy=0;iy<inMRC.numCols (); iy++){
			for (int iz=0; iz<inMRC.numSlcs (); iz++){
				if (inMRC.cube [ix][iy][iz] > 0.0){
					tmppnt.x = ix*inMRC.apixX + inMRC.hdr.xorigin;
					tmppnt.y = iy*inMRC.apixY + inMRC.hdr.yorigin;
					tmppnt.z = iz*inMRC.apixZ + inMRC.hdr.zorigin;

					pnts.push_back (tmppnt);
				}
			}
		}
	}

	//*/

	//write local peaks into a pdb file
	trace = points2pdb(pnts, "HOH", " O  ");
	trace.writePDB(outPath + "_localPeaks.pdb", 1, trace.numOfAA());


	//start from a random end...which is here the first SS and the first end..
	//if use very small number then you want each voxel in a seperate cluster
	peakClustering(pnts, clustersCoord, ssEdges[0][0], 2.0*inMRC.apixX);

	nClusters = clustersCoord.size ();

	/*
	 *				Delete Density around edges sticks
	 */

	pair<Coordinate, Coordinate> *edgesEnds;
	edgesEnds = new pair<Coordinate, Coordinate> [ssEdges.size ()];

	deleteEdgesDensity(inMRC, clustersCoord, ssEdges, edgesEnds, outPath);

	clustersCoord.clear ();

	//start from a random end...which is here the first SS and the first end..
	//if use very small number then you want each voxel in a seperate cluster
	peakClustering(pnts, clustersCoord, ssEdges[0][0], 0.1*inMRC.apixX);


	nClusters = clustersCoord.size ();


	//for SS ends....create a cluster of length = 1 contains one end point for each SS
	vector<Coordinate> tmpCluster;
	for (i=0; i<ssEdges.size (); i++){

		tmpCluster.clear ();
		tmpCluster.push_back (edgesEnds[i].first);
		clustersCoord.push_back (tmpCluster);


		tmpCluster.clear ();
		tmpCluster.push_back (edgesEnds[i].second);
		clustersCoord.push_back (tmpCluster);
	}

	nClustersAll = clustersCoord.size ();
	//print clusters (all points or centriod only)
	for (i=0; i<nClustersAll; i++){
		Protein trace;
		//write local peaks into a pdb file
		trace = points2pdb(clustersCoord[i], "HOH", " O  ");
		trace.writePDB(outPath +"_Cluster_" + toString(i+1) +".pdb", 1, trace.numOfAA());							//all points
		//trace.writePDB(outPath +"_Cluster_" + toString(i+1) +".pdb", trace.numOfAA(), trace.numOfAA());				//centroid only
	}


	/*
	 *		find shortest distances b/w each pair of clusters and then build adjMtrx
	 */
	buildAdjMtrx(clustersCoord, cDist, adjMtrx, inMRC.apixX, nClusters, continuatyTHR);


	float *pairDist;							//maximum distances (acoording to nAA loop) b/w any SSend and all other SSends ...used in local calculations
	float *loopLength;							//the length of the loop b.w two ssEnds...used to find traces locally
	float *tracesLengths;						//the cost to reach each cluster
	short *prevCluster;							//the prev cluster in the trace

	short *nodes;
	float **linksW;

	pairDist = new float [nClustersAll];
	prevCluster = new short [nClustersAll];
	tracesLengths = new float [nClustersAll];

	//initialization links weights
    nodes = new short[nClustersAll];
	linksW = AllocateDynamicArray<float> (nClustersAll, nClustersAll);			//create links weights
	for (i=0; i<nClustersAll;i++)
		for (j=0; j<nClustersAll; j++)
			linksW[i][j] = 99999999.0;

	for (i=0; i<nClustersAll; i++){
		for (j=0; j<adjMtrx[i].size (); j++){
			linksW[i][adjMtrx[i][j]] = getDistance(clustersCoord[i][clustersCoord[i].size ()-1], clustersCoord[adjMtrx[i][j]][clustersCoord[adjMtrx[i][j]].size ()-1]);
		}
	}

	/*
	 *		build the shortest distance b/w each pair of ssEnds (globally- independantly of links on the graph)
	 */
	/*
	short	**ssShortTraces;				//traces from each ssEnd to other clusters
	float	**ssShortLengths;				//lengths of traces from each ssEnd to each other clusters

	float maxLoopLength;

	//find the maximum loop length could ever found in the sequence
	//this can be ignored and choose a large number....actually this used to avoid traces longer than a particular length...this more useful when find traces locally
	short mxMissingSS=0;
	for (i=1; i<graph[0][0].outLinks .size (); i++){
		if (graph[0][0].outLinks[i].rowIndx > mxMissingSS)
			mxMissingSS = graph[0][0].outLinks[i].rowIndx ;
	}

	for (i=0; i+mxMissingSS<seqSS.size (); i+= mxMissingSS){
		if (seqSS[i+mxMissingSS].startIndx - seqSS[i].endIndx -1 > maxLoopLength)
			maxLoopLength = seqSS[i+mxMissingSS].startIndx - seqSS[i].endIndx -1;
	}
	maxLoopLength = (maxLoopLength + 2*MAX_SHIFT_ALLOWED) * 3.8;

	pairDist = new float [nClustersAll];
	prevCluster = new short [nClustersAll];		//no need to initialize ... will be initialized in Dijkstra algorithm
	tracesLengths = new float [nClustersAll];	//no need to initialize.....will be initialized in Dijkstra algorithm
	//set the maximum pair distance to other clusters (not ssEnds) to be the maximum distance found for this end
	for (k=0; k<nClustersAll; k++){
		pairDist[k] = maxLoopLength;
	}

	cout<<"MaxLoopLength = "<<maxLoopLength<<endl;
	//build traces
	ssShortTraces = AllocateDynamicArray<short> (ssEdges.size () * 2, nClustersAll);
	ssShortLengths = AllocateDynamicArray<float> (ssEdges.size () * 2, nClustersAll);

	for (i=0; i<ssEdges.size ()*2; i++){
		sourceCluster = i+nClusters;
		getShortestTraces(nodes, adjMtrx, linksW, pairDist, sourceCluster, prevCluster, tracesLengths);
		for (j=0; j<nClustersAll; j++){
			ssShortTraces[i][j] = prevCluster[j];
			ssShortLengths[i][j] = tracesLengths[j];
		}
	}

	//set links weights for the graph
	for (i=1;i<graph.size ()-1; i++){
		for (j=0; j<graph[i].size (); j++){
			//find the shortest trace paths b/w the pair (this ssEnd and all other clusters) but those shorter than the maximum length b/w this pair
			if (j%2==0)
				sourceCluster = j + 1;
			else
				sourceCluster = j - 1;
			//set weights of links according to shortest trace paths
			for (k=0; k<graph[i][j].outLinks .size (); k++){
				if (graph[i][j].outLinks[k].rowIndx < graph.size ()-1){
					//initial weight = nAA * 3.8 + 15.0
					graph[i][j].outLinks[k].w = graph[i][j].nAAloop [k] * 3.8 + 15.0;

					//for right now ... check complete paths
					if (ssShortLengths[sourceCluster][graph[i][j].outLinks[k].colIndx + nClusters] < graph[i][j].nAAloop [k]*3.8 + 5.0){			//add 5 AA to tolerate error
						if (fabs(graph[i][j].nAAloop [k]*3.8 - ssShortLengths[sourceCluster][graph[i][j].outLinks[k].colIndx + nClusters]) < graph[i][j].outLinks [k].w){
							graph[i][j].outLinks [k].w = fabs(graph[i][j].nAAloop [k]*3.8 - ssShortLengths[sourceCluster][graph[i][j].outLinks[k].colIndx + nClusters]);
							graph[i][j].traceNum[k] = toString(sourceCluster+1);
						}
					}
				}
			}
		}
	}

	FreeDynamicArray<short> (ssShortTraces);
	FreeDynamicArray<float> (ssShortLengths);
	delete [] pairDist;
	delete [] tracesLengths;
	delete [] prevCluster;
	delete [] nodes;
	FreeDynamicArray<float>(linksW);
	*/


	/*
	 *		build the shortest distance b/w each pair of ssEnds (locally- based on the links in the graph)
	 */

	for (i=1; i<graph.size ()-1; i++){
		for (j=0; j<graph[i].size (); j++){

			//cout<<"["<<i<<" "<<j<<"]  : "<<endl;
			//initialize max length of a trace could be found b/w two clusters (this ssEnd and all other clusters)
			for (k=0; k<nClustersAll; k++){
				pairDist[k] = 0;
			}

			loopLength = new float [graph[i][j].outLinks.size()];
			float maxDist=0;					//max length of a trace from this ssEnd to any other cluster
			for (l=0; l<graph[i][j].outLinks.size(); l++){
					// calculate the length (in Angstrom) of the loop b/w two SS ends
					loopLength[l] = 3.8 * (graph[i][j].nAAloop [l] + 1);				//add one amino acid to tolerate error

					if (graph[i][j].outLinks [l].rowIndx != i+1){
						//calculate how many AA in the loop is a hlx AA
						for (m=i+1; m<graph[i][j].outLinks [l].rowIndx; m++){
							if (seqSS[m-1].type == 'H'){				//m-1 b/s i is shifted by 1 (first row is for start node)
								loopLength[l] = loopLength[l] - seqSS[m-1].nAA *3.8 + seqSS[m-1].nAA * ALPHA_RISE;
							}
						}
					}

					if (loopLength[l] > pairDist[nClusters + graph[i][j].outLinks [l].colIndx]){
						pairDist[nClusters + graph[i][j].outLinks [l].colIndx] = loopLength[l];
						if (loopLength[l] > maxDist)
							maxDist = loopLength[l];					//set the maximum loop seen so far
					}
			}

			//set the maximum pair distance to other clusters (not ssEnds) to be the maximum distance found for this end
			for (k=0; k<nClusters; k++){
				if (pairDist[k] == 0)
					pairDist[k] = maxDist;
			}
			//print max possible lengths
			//cout<<"maxDist from cluster# "<<j+nClusters+1<<endl;
			//for (k=0; k<nClustersAll;k++)
			//	cout<<k+1<<" : "<<pairDist[k]<<endl;


			//get shortest Distances to other SS ends

			//find the shortest trace paths b/w the pair (this ssEnd and all other clusters) but those shorter than the maximum length b/w this pair
			if (j%2==0)
				sourceCluster = nClusters + j + 1;
			else
				sourceCluster = nClusters + j - 1;


			getShortestTraces(nodes, adjMtrx, linksW, pairDist, sourceCluster, prevCluster, tracesLengths);


			//set weights of links according to shortest trace paths
			for (k=0; k<graph[i][j].outLinks .size (); k++){
				if (graph[i][j].outLinks[k].rowIndx < graph.size ()-1){
					//cout<<"    ["<<graph[i][j].outLinks[k].rowIndx<<" "<<graph[i][j].outLinks[k].colIndx <<"]  "<<tracesLengths[graph[i][j].outLinks[k].colIndx + nClusters]<<endl;
					//initial weight = nAA * 3.8 + 15.0
					graph[i][j].outLinks[k].w = loopLength[k] + 15.0;

					//for right now ... check complete paths
					if (fabs(loopLength[k] - tracesLengths[graph[i][j].outLinks[k].colIndx + nClusters]) < graph[i][j].outLinks [k].w){
						graph[i][j].outLinks [k].w = fabs(loopLength[k] - tracesLengths[graph[i][j].outLinks[k].colIndx + nClusters]);
						graph[i][j].traceNum[k] = toString(sourceCluster+1);
					}
				}
			}
			delete [] loopLength;
		}
	}
	delete [] pairDist;
	delete [] tracesLengths;
	delete [] prevCluster;
	delete [] nodes;
	FreeDynamicArray<float>(linksW);


}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
iPath getShortestPath(cellSets **graphSets,
				vector<vector<cell> > & graph,
				unsigned int **cntr,								//set counters (contains number of sets in each node)
				short sRow,									//the node to find the shortest path from (given the indx of row and col)
				short sCol,
				short &nRslvd,								//number of nodes resolved in the path so far
				short *skipLinksRow,						//the list of nodes u r not allowd to visit
				short *skipLinksCol,
				short skipLinksSize,
				unsigned int visitedCols){							//the set of cols visited so far

	short nVisitedCols = nRslvd-1;
	short nSticks = graph[0].size ()/2;
	short fromlink, torow, tocol, tocolmod2, icol;
	unsigned int itoset;

	//find the shortest path
	//two cases
	//1. from node is the start node in original graph
	//2. from node is not the start node in original graph
	float shortestW = 9999.0;
	short rowIndx = -1;	//the index of the first row connected with the start node
	short colIndx = -1;
	unsigned int setIndx = -1;


	for (fromlink=0; fromlink<graph[sRow][sCol].outLinks .size (); fromlink++){

		bool found = false;
		for (int i=0; i<skipLinksSize; i++){
			if (graph[sRow][sCol].outLinks[fromlink].rowIndx == skipLinksRow[i] &&
				graph[sRow][sCol].outLinks[fromlink].colIndx == skipLinksCol[i]){
				found = true;
				break;
			}
		}

		if (!found){
			torow = graph[sRow][sCol].outLinks[fromlink].rowIndx;
			tocol = graph[sRow][sCol].outLinks[fromlink].colIndx/2;
			tocolmod2 = graph[sRow][sCol].outLinks[fromlink].colIndx % 2;

			//for all sets on the connected node
			for (itoset=0; itoset<cntr[torow][tocol];itoset++){

				if (graphSets[torow][tocol].length[itoset] + nVisitedCols == nSticks){	//the number of nodes should equal to number of sticks (cols)
					if ((graphSets[torow][tocol].sets [itoset] & visitedCols) == 0){

						if (graphSets[torow][tocol].weight[tocolmod2][itoset] + graph[sRow][sCol].outLinks [fromlink].w < shortestW){
							rowIndx = torow;
							colIndx = graph[sRow][sCol].outLinks[fromlink].colIndx;
							setIndx = itoset;
							shortestW = graphSets[torow][tocol].weight[tocolmod2][itoset]+ graph[sRow][sCol].outLinks [fromlink].w ;
						}
					}
				}
			}
		}
	}

	//find shortest path
	//short  *tmpShortestPathRow;
	//short  *tmpShortestPathCol;

	vector<pair<short,short> > tmpShortestPath (nSticks+2);

	//tmpShortestPath = new pair<short, short> [nSticks+2];			//+2 for start and end nodes
	//tmpShortestPathRow = new short [nSticks+2];						//+2 for start and end nodes
	//tmpShortestPathCol = new short [nSticks+2];						//+2 for start and end nodes

	//pair<short, short> tmpNode;

	if (rowIndx != -1){
		//add the first node in the path

		tmpShortestPath[nRslvd].first = rowIndx;
		tmpShortestPath[nRslvd++].second = colIndx;


		short colIndxDivide2;
		short colIndxMod2;
		short sColMod2;
		short sColDivide2;


		for (icol=0; icol<nSticks-nVisitedCols-1; icol++){
			colIndxDivide2 = colIndx/2;
			colIndxMod2 = colIndx%2;
			sColMod2 = tmpShortestPath[nRslvd-1].second%2;
			sColDivide2 = tmpShortestPath[nRslvd-1].second/2;

			rowIndx = graphSets[rowIndx][colIndxDivide2].nxtRow[colIndxMod2][setIndx];
			colIndx = graphSets[tmpShortestPath[nRslvd-1].first][colIndxDivide2].nxtCol[colIndxMod2][setIndx];
			setIndx = graphSets[tmpShortestPath[nRslvd-1].first][sColDivide2].nxtSetIndx[sColMod2][setIndx];

			tmpShortestPath[nRslvd].first = rowIndx;
			tmpShortestPath[nRslvd++].second = colIndx;
		}
	}


	iPath tmpPath;
	//tmpPath.pCandidateCol = new short [nSticks+2];						//+2 for start and end nodes
	//tmpPath.pCandidateRow = new short [nSticks+2];						//+2 for start and end nodes
    for (icol=0; icol<nSticks+2;icol++){
        //tmpPath.pCandidateRow[icol] = tmpShortestPathRow[icol];
        //tmpPath.pCandidateCol[icol] = tmpShortestPathCol[icol];
        tmpPath.pCandidate.push_back(tmpShortestPath[icol]);
    }

	tmpPath.w = shortestW;

    //clean up
    //delete [] tmpShortestPathRow;
    //delete [] tmpShortestPathCol;
    tmpShortestPath.clear();

	return tmpPath;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//find K-shortest Paths
void findKShortestPaths(cellSets **graphSets, vector<vector<cell> > &graph, short ***topoList, unsigned int **cntr, short nRows, short nSticks, short sRow, short sCol, int K){


#ifdef _WIN				//if work under Windows
	clock_t start, finish;
	start = clock();
/*
#else
	struct timeval ti_start, ti_end;
	double time_di;
	cout<<"Time(0)= "<<time(0)<<endl;
	//start timer
	gettimeofday(&ti_start,0);
*/
#endif

	int k;
	unsigned int icol;

	unsigned int preSet = 0;		//represents the set of current path....the cols included in the current path and we need to avoid

	//create data structure to save skip links for the deviation node...
	short *skipLinks_Row;
	short *skipLinks_Col;

	short maxSkipLinksSize = graph[0].size();				//maximum number of not allowed links is the number of nodes in eaxh row
	short skipLinkSize = 0;									//a pointer to the last link to skip

	skipLinks_Row = new short [maxSkipLinksSize];			//maximum number of not allowed links is the number of out links
	skipLinks_Col = new short [maxSkipLinksSize];


	/*
	 *		Shortest Path . . . top 1
	 */
	short nResolved=1;					//number of resolved nodes so far
	iPath bestPath = getShortestPath(graphSets, graph, cntr, sRow, sCol, nResolved, skipLinks_Row, skipLinks_Col, skipLinkSize, preSet);

	//add the start node
	//bestPath.pCandidateRow[0] = 0;
	//bestPath.pCandidateCol [0] = 0;


	//add end node
	//bestPath.pCandidateRow [nSticks+1] = graph.size ()-1;
	//bestPath.pCandidateCol [nSticks+1] = 0;

	bestPath.pCandidate[nSticks+1].first = graph.size ()-1;


	//bestPath.pCandidate.push_back(tmpRowCol);
    /*
	cout<<" shortest path "<<endl;
	for (int kk =0; kk<nSticks+2; kk++){
		//cout<<"["<<bestPath.pCandidateRow[kk]<<" , "<<bestPath.pCandidateCol[kk]<<"] ";
		cout<<"["<<bestPath.pCandidate[kk].first<<" , "<<bestPath.pCandidate[kk].second<<"] ";
	}
	cout<<" Wt = "<<bestPath.w<<endl;
	*/


	//find the top K-1 shortest Paths
	//initiate T and X sets
	//T is the list of best k solved shortest paths, X is the list of candidates of shortest paths
	iPath *T, *X;
	int Tcntr=0, Xcntr=0, i;
	T = new iPath[K];
	X = new iPath [K*nSticks+K];			//the maximum number of candidates in X would be K*nSticks becuase for each path from T the max possible paths could be generated is the number of edges in the path

	bestPath.parentPath = -1;
	bestPath.deviationNode = 0;

	//add shortest path to X list
	X[Xcntr++] = bestPath;

	k=1;

	while (k<=K){
		//find the shortest path among candidate paths in X
		float shortestW = 9999.0;
		int minX=-1;			//the index of shortest path in X
		for (int iX=0; iX<Xcntr; iX++){
			if (X[iX].w < shortestW){
				minX = iX;
				shortestW = X[iX].w;
			}
		}

		//no more valid paths
		if (minX == -1)
			break;

		//delete shortest path from X and add it to T
		bestPath = X[minX];

		//add best (shortest) path found in X to T...it is the next best (shortest) path
		T[Tcntr++] = bestPath;

		//eleminate this path by increasing the weight
		X[minX].w = 9999.0;

		//print best path in X
		//cout<<" Best Path "<<k<<" : ";
		//for (i=0; i<bestPath.pCandidate .size (); i++)
		//	cout<<"  "<<bestPath.pCandidate [i].second;
		//cout<<"  weight= "<<bestPath.w <<" deviation node= "<<bestPath.deviationNode <<" prnt= "<<bestPath.parentPath<<endl;

		//links that we should not consider (we should skip)
		skipLinkSize = 0;
		//skipLinks_Row [skipLinkSize] = bestPath.pCandidateRow [bestPath.deviationNode +1];
		//skipLinks_Col [skipLinkSize++] = bestPath.pCandidateCol [bestPath.deviationNode +1];
		skipLinks_Row [skipLinkSize] = bestPath.pCandidate [bestPath.deviationNode +1].first;
		skipLinks_Col [skipLinkSize++] = bestPath.pCandidate[bestPath.deviationNode +1].second;


		int parent = bestPath.parentPath;
		while (parent != -1){
			if (T[parent].pCandidate[bestPath.deviationNode ].first == bestPath.pCandidate[bestPath.deviationNode ].first &&
				T[parent].pCandidate[bestPath.deviationNode ].second  == bestPath.pCandidate[bestPath.deviationNode ].second){

				skipLinks_Row [skipLinkSize] = T[parent].pCandidate [bestPath.deviationNode+1].first;
				skipLinks_Col [skipLinkSize++] = T[parent].pCandidate [bestPath.deviationNode +1].second;

			}
			parent = T[parent].parentPath;
		}

		icol = bestPath.deviationNode;
		do{
			nResolved = icol+1;
			//for each node...find the shortest path to the end not passing through the original link in the parent
			iPath tmpPath;

			//save the preselected colms.....the colms visited before the deviation node
			for (i=1; i<=icol ; i++){
				preSet = preSet | (1 << bestPath.pCandidate [i].second/2);
			}

			//get the next best candidate deviates from bestPath
			tmpPath = getShortestPath(	graphSets,
										graph,
										cntr,
										bestPath.pCandidate [icol].first,
										bestPath.pCandidate[icol].second,
										nResolved,
										skipLinks_Row,
										skipLinks_Col,
										skipLinkSize,
										preSet);




			tmpPath.deviationNode = icol;
			tmpPath.parentPath = Tcntr-1;

			//insert first portion of the path
			//insert deviation node
			tmpPath.pCandidate [icol].first = bestPath.pCandidate [icol].first;
			tmpPath.pCandidate [icol].second = bestPath.pCandidate [icol].second;
			for (i=icol-1; i>=0; i--){
				tmpPath.pCandidate[i].first = bestPath.pCandidate[i].first;
				tmpPath.pCandidate[i].second = bestPath.pCandidate[i].second;

				//add weights
				for (int fromlink=0; fromlink<graph[tmpPath.pCandidate[i].first][tmpPath.pCandidate[i].second].outLinks.size(); fromlink++){
					if (graph[tmpPath.pCandidate[i].first][tmpPath.pCandidate[i].second].outLinks [fromlink].rowIndx == tmpPath.pCandidate[i+1].first &&
						graph[tmpPath.pCandidate[i].first][tmpPath.pCandidate[i].second].outLinks [fromlink].colIndx == tmpPath.pCandidate[i+1].second){
						tmpPath.w += graph[tmpPath.pCandidate[i].first][tmpPath.pCandidate[i].second].outLinks [fromlink].w;
						break;
					}
				}
			}

			//insert end node
			tmpPath.pCandidate[nSticks+1].first = bestPath.pCandidate[nSticks+1].first;
			tmpPath.pCandidate[nSticks+1].second = bestPath.pCandidate[nSticks+1].second;

			//print the path
            /*
			for (i=0; i<nSticks+2; i++){
				cout<<"  ["<<tmpPath.pCandidate[i].first<<" "<<tmpPath.pCandidate[i].second<<"] ";
			}
            cout<<"  weight= "<<tmpPath.w <<" dv= "<<tmpPath.deviationNode <<" prnt= "<<tmpPath.parentPath<<endl;
            */

			//save the candidate to X list
			X[Xcntr] = tmpPath;
			Xcntr++;
			//reset preSet
			preSet = 0;
			icol++;
			//reset skiplinks for the new deviation node
			skipLinkSize = 0;
			skipLinks_Row [skipLinkSize] = bestPath.pCandidate[icol+1].first;
			skipLinks_Col [skipLinkSize++] = bestPath.pCandidate[icol+1].second;

			//delete [] tmpPath.pCandidateCol;
			//delete [] tmpPath.pCandidateRow;

		}while(icol<nSticks);
		k+= 1;
		//getchar();
	}

#ifdef _WIN
	finish = clock();
	cout<<"Time taken to enumerate top "<<K<<" paths is "<<finish-start<<" ms."<<endl;
	start = clock();
/*
#else
	//stop timer
	gettimeofday(&ti_end,0);
	time_di = (ti_end.tv_sec-ti_start.tv_sec)*1000000 + ti_end.tv_usec - ti_start.tv_usec;
	cout<<"Time(0)= "<<time(0)<<endl;
	cout<<"Time taken to enumerate top "<<K<<" paths is "<<(double) (time_di/1000)<<" ms "<<endl;
*/
#endif
	for (int iT=0 ;iT<Tcntr; iT++){
		for (icol=1; icol<=nSticks; icol++){
			//saves the topology to topology list
			topoList[iT][icol-1][0] = T[iT].pCandidate[icol].first ;
			topoList[iT][icol-1][1] = T[iT].pCandidate[icol].second;
		}
	}
	//dispose memory 
	/*
	delete []X;
	delete []T;
	delete []skipLinks_Row;
	delete []skipLinks_Col;
	*/
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//finds the shortest path in the graph using a recurrence....saves all combinations from the cell down to the end node
//given is the index of start node (row and col) and the index of end node (row and col)
void gShortestRec(vector<vector<cell> > graph, short ***topoList, int K){

	short irow, icol;

	short nSticks = graph[0].size ()/2;
	short nRows = graph.size();

	//the indeces of START and END nodes
	short	STARTrow=0,
			STARTcol=0,
			ENDrow=nRows-1,
			ENDcol=0;
	//cntr saves the last indx in sets for each node
	unsigned int **cntr = AllocateDynamicArray<unsigned int> (nRows, nSticks);

	for (irow=0; irow<nRows; irow++)
		for (icol=0; icol<nSticks; icol++)
			cntr[irow][icol] = 0;

	//combination index used to hash the position of a combination
	unsigned int *combinationIndx ;

	//initialize data structures
	cellSets **graphSets = AllocateDynamicArray<cellSets> (nRows, nSticks); //number of cols is half number of cols in the original graph

	short fromrow, fromcol, fromlink, torow, tocol, gcol;			//counters
	unsigned int itoset, totalNumSets=0;							//total number of sets for the protein

	//estimate number of sets for each node
	vector<vector<vector<short> > > combSets (nRows, vector<vector< short> > (nSticks));
	combSets[nRows-1][0].push_back (0);
	for (irow=nRows-2; irow>0; irow--){
		for (icol=0; icol<graph[0].size (); icol++){
			for (int ilink=0; ilink<graph[irow][icol].outLinks .size (); ilink++){
				for (int icomb =0; icomb<combSets[graph[irow][icol].outLinks[ilink].rowIndx][graph[irow][icol].outLinks[ilink].colIndx/2].size (); icomb++){
					int k = combSets[graph[irow][icol].outLinks[ilink].rowIndx][graph[irow][icol].outLinks[ilink].colIndx/2][icomb] + 1;
					bool found = false;
					for (int iset=0; iset<combSets[irow][icol/2].size (); iset++){
						if (combSets[irow][icol/2][iset] == k){
							found = true;
							break;
						}
					}
					if (!found)
						combSets[irow][icol/2].push_back (k);
				}
			}
		}
	}
	vector<vector<vector<long> > > nSets (nRows, vector<vector< long> > (nSticks));
	nSets[0].resize (nSticks, vector<long> (1, 0));
	nSets[nRows-1].resize (nSticks, vector<long> (1, 0));
	for (irow=1; irow<nRows-1; irow++){
		for (icol=0; icol<nSticks; icol++){
			long nComb=0;		//number of total sets
			//cout<<"[ "<<irow<<" , "<<icol<<" ] : "<<endl;
			for (int iset=0; iset<combSets[irow][icol].size (); iset++){
				//cout<<"   "<<combSets[irow][icol][iset]<<endl;
				if (combSets[irow][icol][iset] + irow > nSticks)			//exclude short paths
					nComb += (nCombination(nSticks, combSets[irow][icol][iset]) - nCombination(nSticks-1, combSets[irow][icol][iset]));
			}
			nSets[irow][icol].push_back (nComb);
			totalNumSets += nComb;
			//cout<<"[ "<<irow<<" , "<<icol<<" ] : "<<nComb<<endl;
		}
	}

	combSets.clear ();
#ifndef _MPI_RUN
	cout<<"\tPress any key to start TopoDP..."<<endl;getchar();getchar();
	cout<<"\t\t==================== TopoDP ======================"<<endl;
	cout<<"\t\tReserve Memory and Initialize data structures...";
#endif

	unsigned int i;
	int k;

#ifdef _WIN				//if work under Windows
	clock_t start, finish;
	start = clock();
#else
	struct timeval ti_start, ti_end;
	double time_di;
	//start timer
	gettimeofday(&ti_start,0);
#endif

	//initialize indeces table
	unsigned int indecesTableSize = pow(double(2), (double) nSticks);
	combinationIndx = new unsigned int [indecesTableSize];

	for (i=0; i<indecesTableSize; i++)
		combinationIndx[i] = -1;


	//initialize some data structure
	for (irow=1; irow<ENDrow; irow++){
		for (icol=0; icol<nSticks; icol++){

			//cntr[irow][icol] = 0;

			graphSets[irow][icol].sets = new unsigned int [nSets[irow][icol][0]];
			graphSets[irow][icol].nxtCol = AllocateDynamicArray<char> (2, nSets[irow][icol][0]);
			graphSets[irow][icol].nxtRow = AllocateDynamicArray<char> (2, nSets[irow][icol][0]);
			graphSets[irow][icol].nxtSetIndx = AllocateDynamicArray<unsigned int> (2, nSets[irow][icol][0]);
			graphSets[irow][icol].length = new char [nSets[irow][icol][0]];
			graphSets[irow][icol].weight = AllocateDynamicArray<float> (2, nSets[irow][icol][0]);

			//initiate the value of other data structures
			for (i=0; i<nSets[irow][icol][0]; i++){

				graphSets[irow][icol].sets [i] = 0;

				graphSets[irow][icol].nxtCol[0][i] = -1;
				graphSets[irow][icol].nxtCol[1][i] = -1;

				graphSets[irow][icol].nxtRow [0][i] = -1;
				graphSets[irow][icol].nxtRow [1][i] = -1;

				graphSets[irow][icol].nxtSetIndx [0][i] = -1;
				graphSets[irow][icol].nxtSetIndx [1][i] = -1;

				graphSets[irow][icol].weight [0][i] = 9999.0;
				graphSets[irow][icol].weight [1][i] = 9999.0;

				graphSets[irow][icol].length [i] = 0;

			}
		}
	}

	/*
	for (irow=STARTrow+1; irow<ENDrow; irow++){
		for (icol=0; icol<nSticks; icol++){

			cout<<"["<<irow<<"]["<<icol<<"]     nSets= "<<nSets[irow][icol][0]<<endl;
			for (i=0; i<nSets[irow][icol][0]; i++){
				cout<<"  "<<i<<" weight[0]= "<<graphSets[irow][icol].weight [0][i]
				<<" weight[1]= "<<graphSets[irow][icol].weight [1][i]<<" length= "<<int(graphSets[irow][icol].length [i])
				<<" cntr= "<<cntr[irow][icol]<<endl;
			}

		}
	}
	*/

	nSets.clear ();
#ifndef _MPI_RUN
	cout<<"Done."<<endl;
	cout<<"\t\tTopoDP is building sets for nodes...";
#endif
	//the recurrence

	//First Deal with nodes connected with End node
	for (fromrow=ENDrow-1; fromrow>STARTrow; fromrow--){
		for (gcol=0; gcol<graph[fromrow].size (); gcol++){

			short gcolmod2 = gcol%2;
			fromcol = gcol/2;			//the index of the stick

			//get the combination ... is the set contains the current col
			int tmpSet = 1 << fromcol;
			unsigned int setIndx;

			for (fromlink=0; fromlink<graph[fromrow][gcol].outLinks.size (); fromlink++){
				if (graph[fromrow][gcol].outLinks[fromlink].rowIndx == ENDrow &&
					graph[fromrow][gcol].outLinks[fromlink].colIndx == ENDcol ){

					//no sets yet....initiate a set with only one col
					if (combinationIndx[tmpSet] != -1){
						//I have the combination in the list of sets
						setIndx = combinationIndx [tmpSet];
					}
					else{
						//we don't have it in the list yet...add it and save the index
						setIndx = cntr[fromrow][fromcol]++;
						combinationIndx[tmpSet] = setIndx;
						graphSets[fromrow][fromcol].sets[setIndx] = tmpSet;

					}

					// update (initiate) minimum cost and parent information (row, col, set indx)
					graphSets[fromrow][fromcol].nxtCol[gcolmod2][setIndx] = graph[fromrow][gcol].outLinks[fromlink].colIndx;
					graphSets[fromrow][fromcol].nxtRow[gcolmod2][setIndx] = torow;
					graphSets[fromrow][fromcol].nxtSetIndx[gcolmod2][setIndx] = -1;

					graphSets[fromrow][fromcol].weight[gcolmod2][setIndx] = graph[fromrow][gcol].outLinks[fromlink].w;
					graphSets[fromrow][fromcol].length[setIndx] = 1;			//the number of cols in the set is 1
				}
			}
			//re-initiate combination index
			if (gcolmod2 != 0)
				combinationIndx[tmpSet] = -1;


		}
	}

	for (fromrow=ENDrow-2; fromrow>STARTrow; fromrow--){

		for (gcol=0; gcol<graph[fromrow].size(); gcol++){		//graph col
			short gcolmod2 = gcol%2;
			fromcol = gcol/2;

			//test for all out links
			for (fromlink=0; fromlink<graph[fromrow][gcol].outLinks .size (); fromlink++){

				tocol = graph[fromrow][gcol].outLinks [fromlink].colIndx/2;		//the index of the col of the node on the tail of the link
				torow = graph[fromrow][gcol].outLinks [fromlink].rowIndx;		//the indx  of the row of the node on the tail of the link

				if (torow == ENDrow)	//if the node I am connected with is END node...then skip...the set is initiated before
					continue;

				short gtocol = graph[fromrow][gcol].outLinks[fromlink].colIndx%2;

				unsigned int combIndx=-1, setIndx;

				//for all sets on the connected node
				for (itoset=0; itoset<cntr[torow][tocol];itoset++){
					if ((((1 << fromcol) & graphSets[torow][tocol].sets [itoset]) == 0) &&				//exclude paths already have this col
						graphSets[torow][tocol].length[itoset] + fromrow>= nSticks){			//exclude unfeasable paths...paths that too short

						//exclude long paths
						if (graphSets[torow][tocol].length[itoset] + 1 > nSticks)
							continue;

						//check if the col I am connected with has valid paths
						if (graphSets[torow][tocol].weight [gtocol][itoset] == 9999.0)
							continue;

						//get the index of the set
						combIndx = graphSets[torow][tocol].sets[itoset] | (1 << fromcol);


						if (combinationIndx[combIndx] != -1){
							//set already in the list
							setIndx = combinationIndx[combIndx];
						}
						else{
							//set is not in the list of sets....add it and save its index
							setIndx = cntr[fromrow][fromcol]++;
							graphSets[fromrow][fromcol].sets[setIndx] = combIndx;
							combinationIndx [combIndx] = setIndx;
						}

						//update minimum cost and parent information (row, col, setIndx) accordingly
						if (graphSets[fromrow][fromcol].weight[gcolmod2][setIndx] > graphSets[torow][tocol].weight[gtocol][itoset] +  graph[fromrow][gcol].outLinks[fromlink].w){

							graphSets[fromrow][fromcol].length[setIndx] = graphSets[torow][tocol].length[itoset] + 1;
							graphSets[fromrow][fromcol].nxtCol[gcolmod2][setIndx] = graph[fromrow][gcol].outLinks [fromlink].colIndx;
							graphSets[fromrow][fromcol].nxtRow[gcolmod2][setIndx] = torow;
							graphSets[fromrow][fromcol].nxtSetIndx[gcolmod2][setIndx] = itoset;
							graphSets[fromrow][fromcol].weight[gcolmod2][setIndx] = graphSets[torow][tocol].weight[gtocol][itoset] +  graph[fromrow][gcol].outLinks[fromlink].w;
						}
					}
				}
			}

			//re-initialize indeces every two colms
			if (gcolmod2 != 0){
				for (k=0; k<cntr[fromrow][fromcol]; k++){
					combinationIndx [graphSets[fromrow][fromcol].sets [k]] = -1;				//re-initialize only used combinations
				}
			}
		}
	}

#ifdef _WIN
	finish = clock();
	cout<<"Done."<<endl;
	cout<<"\t\tTime taken to build all Sets is "<<finish-start<<" ms."<<endl;
/*
#else
	//stop timer
	gettimeofday(&ti_end,0);
	cout<<"\t\tDone."<<endl;
	time_di = (ti_end.tv_sec-ti_start.tv_sec)*1000000 + ti_end.tv_usec - ti_start.tv_usec;
	cout<<"\t\tTime taken to enumerate top "<<K<<" paths is "<<(double) (time_di/1000)<<" ms "<<endl;
*/
#endif


	/*
	//print information
	cout<<"++++++++++++++++++++++ print Information +++++++++++++"<<endl;
	int j;
	for (irow=0; irow<ENDrow; irow++){
		for (icol=0; icol<nSticks;icol++){
			cout<<"[ "<<irow<<" , "<<icol<<" ]  : "<<endl;
			for (i=0; i<cntr[irow][icol]; i++){
				cout<<"   set# "<<i+1<<" length= "<<int(graphSets[irow][icol].length [i])<<" : "<<graphSets[irow][icol].sets[i]<<endl;
				cout<<"     shortest path F (w= "<<graphSets[irow][icol].weight[0][i]<<") : ";
				cout<<"     shortest path B (w= "<<graphSets[irow][icol].weight[1][i]<<") : "<<endl;
			}
		}
		getchar();
	}
	*/
#ifndef _MPI_RUN
	cout<<"\t\tTopoDP will find and then print "<<K<<"-Shortest Paths Now..."<<endl;

	cout<<"\t\t=================================================="<<endl;
#endif

	findKShortestPaths(graphSets, graph, topoList, cntr, nRows, nSticks, STARTrow, STARTcol, K);
		
	
	//dispose memory
	FreeDynamicArray<cellSets> (graphSets);
	FreeDynamicArray<unsigned int> (cntr);
	delete [] combinationIndx;
	combSets.clear();


}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
//old implementation...gives same result...but less memory
//finds the shortest path in the graph using a recurrence....saves all combinations from the cell down to the end node
//given is the index of start node (row and col) and the index of end node (row and col)
void gShortestRec(vector<vector<cell> > &graph, short sRow, short sCol, short eRow, short eCol){


	int irow, icol;

	int level = 2;		//level we r working on since the start of process
	int nSticks = graph[0].size ()/2;

	//cntr saves the last indx in sets for each node
	vector<vector<long> > cntr(graph.size(), vector<long> (nSticks, 0));

	//initialize data structures
	vector<vector<cellSets> > graphSets (graph.size (), vector<cellSets> (nSticks)); //number of cols is half number of cols in the original graph

	//prepare the graph...
	//the row of start and end nodes should only has outcoming and incoming links from only start and end nodes repectively
	for (irow=0; irow<=sRow; irow++){
		for (icol=0; icol<graph[irow].size (); icol++){
			if (icol != sCol){
				//delete all out links ... no need to delete inlink from other nodes because the recurrence will not work on inlinks
				graph[irow][icol].outLinks .clear ();
			}
		}
	}
	for (irow=eRow; irow<graph.size (); irow++){
		for (icol=0; icol<graph[irow].size (); icol++){
			if (icol != sCol){
				//delete all out links ... no need to delete inlink from other nodes because the recurrence will not work on inlinks
				graph[irow][icol].outLinks .clear ();
			}
		}
	}


	//work on remaining rows
	short fromrow, fromcol, fromlink, torow, tocol, itoset, gcol;

	//initialize some data structure
	for (irow=sRow; irow<=eRow; irow++){
		for (icol=0; icol<nSticks; icol++){
			graphSets[irow][icol].combinationIndx.resize (nSticks, vector<long> (1,-1));
		}
	}
	for (fromrow=eRow-1; fromrow>sRow; fromrow--){
		for (gcol=0; gcol<graph[fromrow].size(); gcol++){		//graph col
			fromcol = gcol/2;

			//cout<<"============================="<<endl;
			//cout<<"row "<<fromrow<<" col "<<gcol<<" #links = "<<graph[fromrow][gcol].outLinks .size ()<<endl;
			for (fromlink=0; fromlink<graph[fromrow][gcol].outLinks .size (); fromlink++){

				tocol = graph[fromrow][gcol].outLinks [fromlink].colIndx/2;		//the index of the col of the node on the tail of the link
				torow = graph[fromrow][gcol].outLinks [fromlink].rowIndx;		//the indx  of the row of the node on the tail of the link


				bool emptySet = true;		//this flag to indicate those node connected with end node

				//cout<<"for link #"<<fromlink+1<<" tocol "<<graph[fromrow][gcol].outLinks[fromlink].colIndx<<" #sets= "<<graphSets[torow][tocol].sets .size ()<<endl;

				//for all sets on the connected node
				for (itoset=0; itoset<graphSets[torow][tocol].sets .size ();itoset++){
					//check if the combination is already saved on sets
					//prepare the combination so we find its indx
					vector<short> tmpComb;
					double combIndx=-1, setIndx=-1;
					bool found = false;


					//checking sets of to node
					//cout<<"checking { ";
					//for (int r=0; r<nSticks; r++)
					//	cout<<graphSets[torow][tocol].sets [itoset][r]<<" ";
					//cout<<" }"<<endl;

					if (!graphSets[torow][tocol].sets [itoset][fromcol]){// &&
						//graphSets[torow][tocol].length[itoset] + fromrow - sRow == nSticks  ){		//if the combination saved does not have from col

						//cout<<"col "<<fromcol<<" is not in the set..."<<endl;

						//check if the combination stored in to node is already stored in from node
						vector<bool> tmpSet = graphSets[torow][tocol].sets [itoset];
						tmpSet[fromcol] = true;

						//get the combination
						//cout<<"looking for { ";
						for (int i=0; i<nSticks; i++){
							if (tmpSet[i]){
								tmpComb.push_back (i+1);
							//	cout<<i+1<<" ";
							}
						}
						//cout<<" }";

						short combSize = tmpComb.size ();
						//find the index
						combIndx = getCombIndx(nSticks, combSize, tmpComb)-1;

						//cout<<"combSize= "<<combSize<<" combIndx= "<<combIndx<<endl;

						if (graphSets[fromrow][fromcol].combinationIndx.size () &&
							graphSets[fromrow][fromcol].combinationIndx[combSize-1].size() > combIndx &&
							graphSets[fromrow][fromcol].combinationIndx[combSize-1][combIndx] != -1){

							setIndx = graphSets[fromrow][fromcol].combinationIndx [combSize-1][combIndx];
							found = true;
						}
						if (!found){

							//cout<<" not found...combinationIndx["<<combSize-1<<"].size= "<<graphSets[fromrow][fromcol].combinationIndx [combSize-1].size ()<<endl;

							vector<short> emptyPath;

							setIndx = graphSets[fromrow][fromcol].sets.size ();

							graphSets[fromrow][fromcol].sets.push_back (tmpSet);
							graphSets[fromrow][fromcol].length .push_back (combSize);
							graphSets[fromrow][fromcol].fPath .push_back (emptyPath);
							graphSets[fromrow][fromcol].bPath .push_back (emptyPath);
							graphSets[fromrow][fromcol].fWeight .push_back (9999.0);
							graphSets[fromrow][fromcol].bWeight .push_back (9999.0);

							if (graphSets[fromrow][fromcol].combinationIndx [combSize-1].size () == 1)
								graphSets[fromrow][fromcol].combinationIndx [combSize-1].resize(nCombination(nSticks, combSize), -1);

							graphSets[fromrow][fromcol].combinationIndx [combSize-1][combIndx] = setIndx;

						}

						//cout<<"search result : "<<found<<" setIndx= "<<setIndx<<endl;

						//compare it with the one already stored
						if (gcol % 2 == 0){		//from node is forward
							if (graph[fromrow][gcol].outLinks [fromlink].colIndx % 2 == 0){		//to node is forward
								if (graphSets[fromrow][fromcol].fWeight[setIndx] > graphSets[torow][tocol].fWeight[itoset] +  graph[fromrow][gcol].outLinks[fromlink].w){
									graphSets[fromrow][fromcol].fPath [setIndx] = graphSets[torow][tocol].fPath[itoset];
									graphSets[fromrow][fromcol].fPath [setIndx].push_back (gcol);
									graphSets[fromrow][fromcol].fWeight[setIndx] = graphSets[torow][tocol].fWeight[itoset] +  graph[fromrow][gcol].outLinks[fromlink].w;
								}
							}
							else{	//to node is backward
								if (graphSets[fromrow][fromcol].fWeight[setIndx] > graphSets[torow][tocol].bWeight[itoset] +  graph[fromrow][gcol].outLinks[fromlink].w){
									graphSets[fromrow][fromcol].fPath [setIndx] = graphSets[torow][tocol].bPath[itoset];
									graphSets[fromrow][fromcol].fPath [setIndx].push_back (gcol);
									graphSets[fromrow][fromcol].fWeight[setIndx] = graphSets[torow][tocol].bWeight[itoset] +  graph[fromrow][gcol].outLinks[fromlink].w;
								}
							}
						}
						else{	//from node is backward
							if (graph[fromrow][gcol].outLinks [fromlink].colIndx % 2 == 0){		//to node is forward
								if (graphSets[fromrow][fromcol].bWeight[setIndx] > graphSets[torow][tocol].fWeight[itoset] +  graph[fromrow][gcol].outLinks[fromlink].w){
									graphSets[fromrow][fromcol].bPath [setIndx] = graphSets[torow][tocol].fPath[itoset];
									graphSets[fromrow][fromcol].bPath [setIndx].push_back (gcol);
									graphSets[fromrow][fromcol].bWeight[setIndx] = graphSets[torow][tocol].fWeight[itoset] +  graph[fromrow][gcol].outLinks[fromlink].w;
								}
							}
							else{	//to node is backward
								if (graphSets[fromrow][fromcol].bWeight[setIndx] > graphSets[torow][tocol].bWeight[itoset] +  graph[fromrow][gcol].outLinks[fromlink].w){
									graphSets[fromrow][fromcol].bPath [setIndx] = graphSets[torow][tocol].bPath[itoset];
									graphSets[fromrow][fromcol].bPath [setIndx].push_back (gcol);
									graphSets[fromrow][fromcol].bWeight[setIndx] = graphSets[torow][tocol].bWeight[itoset] +  graph[fromrow][gcol].outLinks[fromlink].w;
								}
							}
						}

					}
					emptySet = false;
				}

				if (emptySet){
					//cout<<"Empty set.."<<endl<<endl;
					//this case for nodes connected with end node (eRow, eCol) ... the only allowed dead end is the end node
					vector<short> tmpComb;
					if (torow < graph.size()-1){
						if (fromcol<tocol){
							tmpComb.push_back (fromcol+1);
							tmpComb.push_back (tocol+1);
						}
						else{
							tmpComb.push_back (tocol+1);
							tmpComb.push_back (fromcol+1);
						}
					}
					else{
						tmpComb.push_back(fromcol+1);
					}

					short combSize = tmpComb.size();

					//find the index of the combination
					double combIndx = getCombIndx(nSticks, combSize, tmpComb)-1;

					//cout<<"combSize= "<<combSize<<" {"<<tmpComb[0]<<"}  C("<<nSticks<<" , "<<combSize<<") combIndx= "<<combIndx<<endl;

					//initiate and update indeces
					//each row represents combination list...for example first row represent n choose 1, second row represents n choose 2 ...etc
					//total number of rows should be nSticks b/s the number of k could be choosen ranges from 1 to n which is equal to nSticks
					//later on...each row will be resized to be of size nCombination(n,k)
					//number of rows in combinationIndx is initialized to nSticks
					if (graphSets[fromrow][fromcol].combinationIndx[combSize-1].size()<nSticks)
						graphSets[fromrow][fromcol].combinationIndx[combSize-1].resize (nCombination(nSticks, combSize),-1);

					double setIndx = graphSets[fromrow][fromcol].sets.size();
					//save the index of the combination
					graphSets[fromrow][fromcol].combinationIndx[combSize-1][combIndx] = setIndx;

					//to save information for the combination
					//each row is a combination...first col is for combination member number 1, second col is for combination member number 2
					//if the number is included in the set then its value would be true
					//for example if we have a combination set {1,3,5} then col 0, 2, and 4 will be true..and other cols will be false
					//each time we will expand (resize) the rows by 1
					vector<bool> emptySet(nSticks, false);
					graphSets[fromrow][fromcol].sets.push_back (emptySet);//.resize (graphSets[fromrow][fromcol].sets.size() + 1, vector<bool> (nSticks, false));

					graphSets[fromrow][fromcol].sets[setIndx][fromcol] = true;
					if (combSize>1)
						graphSets[fromrow][fromcol].sets[setIndx][tocol] = true;

					//resize data structure to keep the length of the set
					graphSets[fromrow][fromcol].length.push_back(combSize);


					//resize data structure to save the shortest path for forward and backward nodes

					//graphSets[fromrow][fromcol].fPath.resize (graphSets[fromrow][fromcol].fPath.size()+1, vector<short> (combSize, -1));
					//graphSets[fromrow][fromcol].bPath.resize (graphSets[fromrow][fromcol].bPath.size()+1, vector<short> (combSize, -1));

					vector<short> emptyComb(combSize, -1);
					graphSets[fromrow][fromcol].fPath.push_back (emptyComb);
					graphSets[fromrow][fromcol].bPath.push_back (emptyComb);
					graphSets[fromrow][fromcol].fWeight.push_back (9999.0);
					graphSets[fromrow][fromcol].bWeight.push_back (9999.0);


					if ( gcol % 2 == 0){

						if (combSize>1){
							graphSets[fromrow][fromcol].fPath[setIndx][0] = graph[fromrow][gcol].outLinks[fromlink].colIndx;
							graphSets[fromrow][fromcol].fPath[setIndx][1] = gcol;
						}
						else
							graphSets[fromrow][fromcol].fPath[setIndx][0] = gcol;

						graphSets[fromrow][fromcol].fWeight[setIndx] = graph[fromrow][gcol].outLinks[fromlink].w;
					}
					else{

						if (combSize>1){
							graphSets[fromrow][fromcol].bPath[setIndx][0] = graph[fromrow][gcol].outLinks[fromlink].colIndx;
							graphSets[fromrow][fromcol].bPath[setIndx][1] = gcol;
						}
						else
							graphSets[fromrow][fromcol].bPath[setIndx][0] = gcol;

						graphSets[fromrow][fromcol].bWeight[setIndx] = graph[fromrow][gcol].outLinks[fromlink].w;
					}
				}
			}
		}
	}
	//find the shortest path
	//two cases
	//1. from node is the start node in original graph
	//2. from node is not the start node in original graph
	float shortest = 100000000;
	vector<short> shortestPath;
	if (sRow == 0){
		for (fromlink=0; fromlink<graph[sRow][sCol].outLinks .size (); fromlink++){

			torow = graph[sRow][sCol].outLinks[fromlink].rowIndx;
			tocol = graph[sRow][sCol].outLinks[fromlink].colIndx/2;
			//if from node is the start node...then we need just to save shortest paths with length equals to nSticks
			//for all sets on the connected node
			for (itoset=0; itoset<graphSets[torow][tocol].sets .size ();itoset++){
				if (graphSets[torow][tocol].length[itoset] == nSticks){		//check if it is equal to nSticks
					if (graphSets[torow][tocol].fWeight[itoset] < shortest){
						shortest = graphSets[torow][tocol].fWeight[itoset];
						shortestPath = graphSets[torow][tocol].fPath [itoset];
					}
					if (graphSets[torow][tocol].bWeight[itoset] < shortest){
						shortest = graphSets[torow][tocol].bWeight[itoset];
						shortestPath = graphSets[torow][tocol].bPath [itoset];
					}
				}
			}
		}
	}

	cout<<"shortest path is : "<<endl;
	for (int kk =shortestPath.size()-1;kk>=0; kk--)
		cout<<shortestPath[kk]<<" ";
	cout<<" weight = "<<shortest<<endl;
	getchar();
	getchar();

//	allSets.clear ();

}
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void estimateNumPaths(vector<vector<cell> > & g)
{
	int i, j, k, m;

	g[g.size()-1][0].nPaths = 1;		//set number of paths from end node which is equal to zero
	for (i=g.size ()-1;i>=0; i--)		//start from buttom
	{
		for (j=0; j<g[i].size(); j++)
		{
			for (k=0; k<g[i][j].outLinks .size(); k++)
				g[i][j].nPaths += g[g[i][j].outLinks [k].rowIndx][g[i][j].outLinks [k].colIndx].nPaths;

		}
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//given an axis...draw a line every dist and the axis number of points should be exactly as numPoints
//if the axis is taller trim it or if it is shorter extend it
void tuneAxis(vector<Coordinate> &axis, int numPoints, float dist){

	//axis is longer than number of points needed....trim it
	if (axis.size () > numPoints){
		axis.erase (axis.begin ()+numPoints, axis.end ());
	}
	//axis is shorter than the given number of points...need to draw many points equal to the difference
	if (axis.size ()<numPoints){
		int i = axis.size ();
		while (i<numPoints){
			axis.push_back (pointOnLine(axis[axis.size ()-2], axis[axis.size ()-1], dist));
			i++;
		}
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Protein mergeTwoHlces(vector<Coordinate>& axis, Protein &tmpProt, int nRise, int nAA){

	int i;

	// Now overlap the beginning of the generated structure with the curved line
	//set moving points which represent the axis of the generated structure
	vector<Coordinate> mPoints;
	setAxis(tmpProt, nRise, mPoints);

	//set points on target axis
    axis = extractCurves(axis, ALPHA_RISE, ALPHA_RISE * nRise, 1);

	if (axis.size () != mPoints.size()){
		tuneAxis(axis, mPoints.size (), ALPHA_RISE * nRise);
	}

    /*
    Protein curvd;
    curvd = points2pdb(axis);
    curvd.writePDB("axis.pdb", 1, curvd.numOfAA());
    */

	//set moving points
	Coordinate p;
	p.x	= axis[0].x - mPoints[0].x;
	p.y	= axis[0].y - mPoints[0].y;
	p.z	= axis[0].z - mPoints[0].z;

	tmpProt.translateBy(p);		//move start to start
	Vectors v1(axis[1],axis[0]), //v1(axis[axis.size()-1], axis[0]),
			v2(tmpProt.AAs[nRise-1].atoms[0].coord, axis[0]), //v2(tmpProt.AAs[tmpProt.numOfAA()-2].atoms[0].coord, axis[0]),
			normal;

	//get angle b/w two vectors...to make them on the same direction (roughly)
	double angle = v1.getAngleDegree(v2);

	//find the normal b/w two vectors
	normal = v1.cross(v2);

	normal += axis[0];

	tmpProt.rotate(0, tmpProt.numOfAA()-1, 0, axis[0], normal.getCoordinates(), angle);

	//tmpProt.writePDB("builtHlx.pdb",1,tmpProt.numOfAA());
	//reset moving points
	setAxis(tmpProt, nRise, mPoints);		//again after translation and rotation


	if (mPoints.size() == axis.size()){
		//overlap the moving points (mPoints) with target points nEdge by modifying tmpProt structure
		overlapFBCCD(tmpProt, nRise, mPoints, axis, 0.1, 100);
	}
	else
	{
		cout<<"=================================== in mergeTwoHlces(..........) ================="<<endl;
		cerr << "Number of moving points ( "<<mPoints.size()<<" ) and target points ( "<<axis.size()<<" ) are not equal."<< endl;
		cout<<"==================================================================================="<<endl;
		exit(1);
	}

	return tmpProt;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void traceSticksCollision(vector<vector<Coordinate> > &traceToAvoid, vector<Coordinate> &trace, int size){
	int i,j,k;

	//cout<<"remove duplicate points on trace...";
    //remove duplicate points on trace
    i=0;
    while (i<trace.size()-1){
        int j=i+1;
        while (j<trace.size()){
            if (getDistance(trace[j],trace[i])<0.5)
                trace.erase(trace.begin()+j--);
            j++;
        }
        i++;
    }
	//cout<<" Done.."<<endl;

/*
	bool collide = true;
	while (collide){
		collide = false;
		for (k=2; k<trace.size ()-2;k++){
			for (i=0; i<size; i++){
				for (j=1; j<traceToAvoid[i].size ()-1; j++){
					//get the distance b/w  the two line segments
					while (getDistLines(traceToAvoid[i][j], traceToAvoid[i][j+1], trace[k], trace[k+1]) < 5){
						collide = true;
						cout<<"Collision b/w the trace point# "<<k+1<<" and stick# "<<i+1<<endl;

						//move the point away from the stick segment one Angstrom at a time
						trace[k] = pointOnLine(traceToAvoid[i][j+1], trace[k], getDistance(traceToAvoid[i][j+1], trace[k])+1);
						trace[k+1] = pointOnLine(traceToAvoid[i][j+1], trace[k+1], getDistance(traceToAvoid[i][j+1], trace[k+1])+1);

						//getchar();
					}
				}
			}
		}
		cout<<endl;
	}
	//getchar();getchar();
*/
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double getRMSDSSE(vector<Protein> &portions, Protein &nativePDBFile, int *sIndeces)
{
	////////////get RMSD of SSW with native

	int i, length=0;
	double rmsd = 0;
	Coordinate p1,p2;
	//nX3 coordinates matrices
	double **native, **model;

	for (i=0;i<portions.size ();i++){
		length += portions[i].numOfAA();
	}
	//resize array of coordinates
	native = AllocateDynamicArray<double> (length,3);
	model =  AllocateDynamicArray<double> (length,3);


	//other data structure needed in calculating RMSD...see rmsd for details
	double mov_com[3], mov_to_ref[3], U[3][3];


	int gCntr=0;
	for (i=0;i<portions.size();i++)
	{
		int startIndx = sIndeces[i];

		for (int j=0; j<portions[i].numOfAA (); j++)
		{
			//cout<<portions[i].AAs[j].chr3<<" "<<portions[i].AAs[j].num<<" pdbfile startAA = "<<nativePDBFile.AAs [startIndx].chr3	\
				<<"  "<<nativePDBFile.AAs [startIndx].num<<" shift= "<<sticks[i].shift<<"  direction= "<<sticks[i].direction<<endl;
			p1 = portions[i].getAtomCoordinate(j, " CA ");
			p2 = nativePDBFile.getAtomCoordinate(startIndx, " CA ");

			native[gCntr][0] = p2.x;  native[gCntr][1] = p2.y;  native[gCntr][2] = p2.z;
			model[gCntr][0]  = p1.x;  model[gCntr][1]  = p1.y;  model[gCntr][2]  = p1.z;

			startIndx++;
			gCntr++;

			//rmsd += ((p1.x - p2.x) * (p1.x - p2.x)) + ((p1.y - p2.y) * (p1.y - p2.y)) + ((p1.z - p2.z) * (p1.z - p2.z));
		}
	}

	//one way of calculating RMSD....finding the rotation matrix first
	calculate_rotation_rmsd (native, model, length, mov_com, mov_to_ref, U, &rmsd);


	return rmsd;

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Protein buildFullModel(Protein & pdb, vector<vector<cell> > graph, vector<vector<Coordinate> > hEdges, short **topology, Coordinate **traceList, float *traceLength, short *traceNoOfPoints, short * shifts, short &sIndx, double &rmsd, Protein &allSSE){

	int i,j,k;

	/*
	//print top topologies info and traces info
	for (i=0; i<K; i++){
		cout<<"Topo# "<<i+1<<" : ";
		for (j=0; j<graph[0].size ()/2; j++){
			cout<<"["<<topTopos[i][j][0]<<" "<<topTopos[i][j][1]<<"] ";
		}
		cout<<endl;
	}

	for (i=1; i<traceNoOfPoints[0]; i++){
	    cout<<"trace# "<<i<<"  length= "<<traceLength[i]<<" noPoints= "<<traceNoOfPoints[i]<<endl;
	}
	*/

    size_t found;
    string key = "-";
    short path1, path2;
	short nRise = 4;
	const short nSticks= hEdges.size ();



 //   for (i=0; i<K; i++){

        Protein fullModel;									// a temp protein model to save the current topology model

		vector<short> stickIndx (nSticks);					//the indeces of the sticks for the model
		vector<short> hlxIndx (nSticks);					//the indeces of the hlces
		vector<short> stickDirection (nSticks);				//the direction of the stick
		vector<Protein> skeleton (nSticks);					//skeletons of sticks
		vector<short> startIndx (nSticks);					//starting indeces of sticks on the sequence
		vector<short> endIndx (nSticks);                    //the end index on the sequence
		vector<vector<short> > trace (nSticks-1, vector<short> (2));			//the traces for each index		traces[0][0] first sub-trace traces[0][1] second subtrace for first loop
		vector<short> traceNumPnts (nSticks-1);				//the number of points in each trace
		vector<float> traceNaa (nSticks-1);					//the expected length (in terms of num of AA) for each trace (trace 1 is b/w stick1 and 2)
		vector<short> loopNaa (nSticks-1);					//number of actual amino acids in the loop after assigning sticks to helices

		vector<vector<Coordinate> >  newEdges(nSticks);				//the new edges (follow the direction of the topology
		vector<vector<Coordinate> > tracePnts(nSticks-1);		    //the actual 3-d points of each trace



        for (j=0; j<nSticks-1; j++){	//for all sticks, find the traces b/w each two sticks

			Protein tmpLoop;
			//float nAALoopMap = 0;

			short rowIndx = topology[j][0],
					colIndx = topology[j][1];

			//set hlces and sticks indeces and other information
			stickIndx[j]	= (int) colIndx/2;
			hlxIndx[j]		= rowIndx;
			if (colIndx%2==0)
				stickDirection[j] = 1;
			else
				stickDirection[j] = -1;
			skeleton[j]		= graph[rowIndx][colIndx].skeleton;

			//set the index of the start AA of this stick on the sequence
			//cout<<"stick "<<stickIndx[j]+1<<" shift= "<<shifts[j]<<" vLeft= "<<graph[rowIndx][colIndx].validityLeft<<" vRight= "<<graph[rowIndx][colIndx].validityRight<<" start @ "<<graph[rowIndx][colIndx].startAAindx<<endl;
			if (shifts[j]<0){       //shift to the left
			    if (graph[rowIndx][colIndx].validityLeft >= abs(shifts[j])){
                    startIndx[j]	= graph[rowIndx][colIndx].startAAindx + shifts[j];
                    endIndx[j]      = startIndx[j] + skeleton[j].numOfAA()-1;

                    graph[topology[j+1][0]][topology[j+1][1]].validityLeft += abs(shifts[j]);       //add the valid (available) number of AAs to the left of the next stick
			    }
			    else{   //or shift the maximum allowed
                    startIndx[j]	= graph[rowIndx][colIndx].startAAindx - graph[rowIndx][colIndx].validityLeft;
                    endIndx[j]      = startIndx[j] + skeleton[j].numOfAA()-1;

                    graph[topology[j+1][0]][topology[j+1][1]].validityLeft += graph[rowIndx][colIndx].validityLeft;       //add the valid (available) number of AAs to the left of the next stick
			    }
			}
			else{
			    if (graph[rowIndx][colIndx].validityRight >= abs(shifts[j])){
                    startIndx[j]	= graph[rowIndx][colIndx].startAAindx + shifts[j];
                    endIndx[j]      = startIndx[j] + skeleton[j].numOfAA()-1;

                    graph[topology[j+1][0]][topology[j+1][1]].validityLeft -= shifts[j];
			    }
			    else{
                    startIndx[j]	= graph[rowIndx][colIndx].startAAindx + graph[rowIndx][colIndx].validityRight;
                    endIndx[j]      = startIndx[j] + skeleton[j].numOfAA()-1;

                    graph[topology[j+1][0]][topology[j+1][1]].validityLeft -= graph[rowIndx][colIndx].validityRight;
			    }
			}
			//startIndx[j]	= graph[rowIndx][colIndx].startAAindx;
			//endIndx[j]      = startIndx[j] + skeleton[j].numOfAA()-1;

            for (k=0; k<graph[rowIndx][colIndx].outLinks.size(); k++){
                if (graph[rowIndx][colIndx].outLinks[k].rowIndx==topology[j+1][0] &&         //find trace number to model the loop
                    graph[rowIndx][colIndx].outLinks[k].colIndx==topology[j+1][1]){

					//extract the traces numbers from the string that we have used to save traces for that link
					//and calculate the length of the trace in terms of number of AAs
                    found = graph[rowIndx][colIndx].traceNum[k].rfind(key);
                    if (found != graph[rowIndx][colIndx].traceNum[k].npos){

                        //get the index of first sub-trace...if the sub-trace is continous...then path1=path2
                        path1 = atoi(graph[rowIndx][colIndx].traceNum[k].substr(0, found).c_str());         //the index of the first sub-trace
                        path2 = atoi(graph[rowIndx][colIndx].traceNum[k].substr(found+1, graph[rowIndx][colIndx].traceNum[k].length()).c_str());    //the index of the second sub-trace

						trace[j][0] = path1;
						trace[j][1] = path2;

						//add points to the vector
						if (path1!=0 && path2!=0){
							int  mm;
							for (mm=0; mm<traceNoOfPoints[path1]; mm++){
								tracePnts[j].push_back (traceList[path1][mm]);
							}
							traceNaa[j] = traceLength[path1];

							if (path1 != path2){
								traceNaa[j] += traceLength[path2];

								for (int kk=traceNoOfPoints[path2]-1;kk>-1; kk--){
									tracePnts[j].push_back (traceList[path2][kk]);
								}

								//get the distance b/w the two sub-trace
								traceNaa[j] += getDistance(tracePnts[j][traceNoOfPoints[path1]-1], tracePnts[j][traceNoOfPoints[path1]]);
							}
							traceNumPnts[j] = tracePnts[j].size ();

							traceNaa[j] = (int) (traceNaa[j]/3.8 + 0.5);
						}
						break;
                    }
                }
            }
			//get the information of the actual loop...
			//number of amino acid to add from left
			//loopNaa[j] = (graph[rowIndx][colIndx].validityLeft  < MAX_SHIFT_ALLOWED) ? (graph[rowIndx][colIndx].validityLeft) : (MAX_SHIFT_ALLOWED);
			//number of amino acids to add from right
			//loopNaa[j] += (graph[graph[rowIndx][colIndx].outLinks[k].rowIndx][graph[rowIndx][colIndx].outLinks[k].colIndx].validityRight < MAX_SHIFT_ALLOWED) ? (graph[graph[rowIndx][colIndx].outLinks[k].rowIndx][graph[rowIndx][colIndx].outLinks[k].colIndx].validityRight ) : (MAX_SHIFT_ALLOWED);
			//number of amino acids in between
			//loopNaa[j] += graph[graph[rowIndx][colIndx].outLinks[k].rowIndx][graph[rowIndx][colIndx].outLinks[k].colIndx].startAAindx - graph[rowIndx][colIndx].startAAindx - graph[rowIndx][colIndx].skeleton .numOfAA();

		}
		//get the information of the last stick
		hlxIndx[j] = topology[j][0];
		stickIndx[j] = (int) topology[j][1]/2;
		if (stickIndx[j]%2==0)
			stickDirection[j] = 1;
		else
			stickDirection[j] = -1;
		skeleton[j]		= graph[hlxIndx[j]][topology[j][1]].skeleton;

        //cout<<"stick "<<stickIndx[j]+1<<" shift= "<<shifts[j]<<" vLeft= "<<graph[hlxIndx[j]][topology[j][1]].validityLeft<<" vRight= "<<graph[hlxIndx[j]][topology[j][1]].validityRight<<" start @ "<<graph[hlxIndx[j]][topology[j][1]].startAAindx<<endl;
		//set the index of the start AA of this stick on the sequence (for the last stick)
        if (shifts[j]<0){       //shift to the left
            if (graph[hlxIndx[j]][topology[j][1]].validityLeft >= abs(shifts[j])){
                startIndx[j]	= graph[hlxIndx[j]][topology[j][1]].startAAindx + shifts[j];
                endIndx[j]      = startIndx[j] + skeleton[j].numOfAA()-1;

                graph[topology[j-1][0]][topology[j-1][1]].validityRight--;       //add the valid (available) number of AAs to the left of the next stick
            }
            else{
                startIndx[j]	= graph[hlxIndx[j]][topology[j][1]].startAAindx - graph[hlxIndx[j]][topology[j][1]].validityLeft;
                endIndx[j]      = startIndx[j] + skeleton[j].numOfAA()-1;

                graph[topology[j-1][0]][topology[j-1][1]].validityRight -= graph[hlxIndx[j]][topology[j][1]].validityLeft;       //add the valid (available) number of AAs to the left of the next stick
            }
        }
        else{
            if (graph[hlxIndx[j]][topology[j][1]].validityRight >= abs(shifts[j])){
                startIndx[j]	= graph[hlxIndx[j]][topology[j][1]].startAAindx + shifts[j];
                endIndx[j]      = startIndx[j] + skeleton[j].numOfAA()-1;

                graph[topology[j-1][0]][topology[j-1][1]].validityRight++;
            }
            else{
                startIndx[j]	= graph[hlxIndx[j]][topology[j][1]].startAAindx + graph[hlxIndx[j]][topology[j][1]].validityRight;
                endIndx[j]      = startIndx[j] + skeleton[j].numOfAA()-1;

                graph[topology[j-1][0]][topology[j-1][1]].validityRight += graph[hlxIndx[j]][topology[j][1]].validityRight;
            }
        }
		//startIndx[j]	= graph[hlxIndx[j]][topology[j][1]].startAAindx;
		//endIndx[j]      = startIndx[j] + skeleton[j].numOfAA()-1;

        //resolve the overlap b/w consecutive sticks on the sequence
        //if the two consecutive sticks collide on the sequence....try to shift the one to the left and one to the right (if possible)..otherwise reduce the number of amino acid on each one
        short nAAOverlap;
        short ic;

        for (j=1; j<nSticks;j++){
            if (startIndx[j] <= endIndx[j-1]){
            //if (loopNaa[j-1]<0){
                //cout<<"overlap"<<endl;
                //cout<<"  stick# "<<stickIndx[j-1]+1<<" start @ "<<startIndx [j-1]<<" ands @ "<<endIndx[j-1]<<" length= "<<skeleton[j-1].numOfAA()<<endl;
                //cout<<"  stick# "<<stickIndx[j]+1<<" start @ "<<startIndx [j]<<" ands @ "<<endIndx[j]<<" length= "<<skeleton[j].numOfAA()<<endl;

                nAAOverlap = abs(endIndx[j-1]-startIndx[j])+1;//abs(pdb.hlces[hlxIndx[j-1]-1].endIndx - pdb.hlces[hlxIndx[j]-1].startIndx)+1;
                //try to shift the first one to the left
                //the index of the end of the stick in the left
                short lastIndx =-1;
                if (j>1)
                    lastIndx = startIndx[j-2]+skeleton[j-2].numOfAA()-1;
                else
                    lastIndx = 0;

                short nAABetween = startIndx[j-1] - lastIndx-1;
                if (nAABetween > nAAOverlap){
                    //cout<<"stick# "<<stickIndx[j-1]+1<<" shifted "<<nAAOverlap<<" to the left"<<endl;
                    //shift the first helix to the left
                    startIndx[j-1] -= nAAOverlap;
                    nAAOverlap = 0;

                }
                else{
                    if (nAABetween>0){
                        //cout<<"stick# "<<stickIndx[j-1]+1<<" shifted "<<nAABetween<<" to the left"<<endl;
                        startIndx[j-1] -= nAABetween;
                        nAAOverlap -= nAABetween;

                    }
                }
                //try to shift the second stick to the right
                //the index of the start of the next AA
                if (j==nSticks-1)
                    lastIndx = pdb.numOfAA()-1;
                else
                    lastIndx = startIndx[j+1];
                nAABetween = lastIndx - (startIndx[j]+skeleton[j].numOfAA()-1)-1;
                if (nAABetween > nAAOverlap){
                    //shift the first helix to the left
                    startIndx[j] += nAAOverlap;
                    //cout<<"stick# "<<stickIndx[j]+1<<" shifted "<<nAAOverlap<<" to the right"<<endl;
                }
                else{
                    if (nAABetween>0){
                        startIndx[j] += nAABetween;
                        //cout<<"stick# "<<stickIndx[j]+1<<" shifted "<<nAABetween<<" to the right"<<endl;
                    }
                }
                endIndx[j-1] = startIndx[j]-1;
            }
        }
        //update loop number of AAs
        for (j=0; j<nSticks-1; j++)
            loopNaa[j] = startIndx[j+1]-endIndx[j]-1;

        /*
		for (j=0; j<nSticks-1; j++){
			cout<<"hlx# "<<hlxIndx[j]<<" to Stck# "<<stickIndx[j]+1<<" and hlx# "<<hlxIndx[j+1]<<" to stck# "<<stickIndx[j+1]+1<<endl;
			cout<<"   stck# "<<stickIndx[j]+1<<" start at AA indx "<<startIndx[j]<<" and ends at "<<endIndx[j]<<endl;
			cout<<"   stck# "<<stickIndx[j+1]+1<<" starts at AA indx"<<startIndx[j+1]<<" and ends at "<<endIndx[j+1]<<endl;
			cout<<"   stick 1 nAA = "<<skeleton[j].numOfAA()<<" ("<<endIndx[j]-startIndx[j]+1<<") stick 2 nAA= "<<skeleton[j+1].numOfAA()<<" ("<<endIndx[j+1]-startIndx[j+1]+1<<")"<<endl;
			cout<<"   hlx 1 nAA = "<<pdb.hlces[hlxIndx[j]-1].nAA<<"  hlx2 nAA= "<<pdb.hlces[hlxIndx[j+1]-1].nAA<<endl;
			cout<<"   Trace found : "<<trace[j][0]<<"-"<<trace[j][1]<<" expected nAA= "<<traceNaa[j]<<endl;
			cout<<"   numOfAA in the loop on the Sequence is : "<<pdb.hlces[hlxIndx[j+1]-1].startIndx - pdb.hlces[hlxIndx[j]-1].endIndx-1<<endl;
			cout<<"   loop nAA : "<<loopNaa[j]<<endl;
		}
        */

		//update sticks edges according to the direction of the topology
		for (j=0; j<nSticks; j++){
			if (stickDirection[j]==-1){
				//reverse the edge
				for (int m=hEdges[stickIndx[j]].size ()-1;m>-1; m--)
					newEdges[j].push_back(hEdges[stickIndx[j]][m]);
				//cout<<"reverse stick found...stick# "<<stickIndx[j]+1<<" number of points in new Edge= "<<newEdges[j].size ()<<" in old edge= "<<hEdges[stickIndx[j]].size ()<<endl;
			}
			else
				newEdges[j] = hEdges[stickIndx[j]];
		}
		//search for consecutive short loops and build all hlces at once
		vector<vector<short> > mergedSticks;
		short mergedLastIndx = -1;
		vector<short> tmp;
		j=0;
		while (j<nSticks-1){
			//cout<<"loop# "<<j+1<<" with loop# "<<j+2<<" dist= "<<getDistance(newEdges[j][newEdges[j].size()-1],newEdges[j+1][0])<<endl;
			if (loopNaa[j] < 4 && (loopNaa[j] < 3 || getDistance(newEdges[j][newEdges[j].size()-1],newEdges[j+1][0])<10)){
				//check if this loop already pushed in the previous step
				if (mergedLastIndx > -1 && mergedSticks[mergedLastIndx][mergedSticks[mergedLastIndx].size()-1] == j){
					mergedSticks[mergedLastIndx].push_back(j+1);
					//cout<<"added "<<j+1<<" (stck# "<<stickIndx[j+1]+1<<")"<<endl;
				}
				else{
					//cout<<"pushed "<<j<<" ( stk# "<<stickIndx[j]+1<<") and "<<j+1<<" (stick# "<<stickIndx[j+1]+1<<" )"<<endl;
					tmp.clear();
					tmp.push_back(j);
					tmp.push_back(j+1);
					mergedSticks.push_back(tmp);
					mergedLastIndx += 1;
				}
			}
			else{
				//add the axis of the stick ... for collision detection later
				newEdges[j] = extractCurves(newEdges[j], ALPHA_RISE, ALPHA_RISE * nRise, 1);

			}
			j++;
		}
		//add the axis of the last helix
		if (mergedSticks.size()==0 || mergedSticks[mergedLastIndx][mergedSticks[mergedLastIndx].size()-1] != j)
			newEdges[j] = extractCurves(newEdges[j], ALPHA_RISE, ALPHA_RISE * nRise, 1);

		//cout<<"Merged Sticks..."<<endl;

		for (j=0; j<=mergedLastIndx; j++){
		    vector<Coordinate> axis;
		    short newHlxLength=0, kk;
			for (kk=0; kk<mergedSticks[j].size(); kk++){
				//cout<<stickIndx[mergedSticks[j][kk]]+1<<" ";
				for (ic=0;ic<newEdges[mergedSticks[j][kk]].size(); ic++){
				    //save all edges to one extended edge..save all points
				    axis.push_back(newEdges[mergedSticks[j][kk]][ic]);
				}
				newHlxLength += endIndx[mergedSticks[j][kk]]-startIndx[mergedSticks[j][kk]]+1;
				newHlxLength += loopNaa[mergedSticks[j][kk]];
			}
			//cout<<endl;
			//delete the last loop added
			newHlxLength -= loopNaa[mergedSticks[j][kk-1]];
			//cout<<" new long hlx length is : "<<newHlxLength<<endl;

			// build a random initial structure
			Protein tmpProt;
			AminoAcid tmpAA;

			tmpAA.chr3 = "GLY";
			tmpAA.num = 1;
			tmpAA.angles.phi = ALPHA_HLX_PHI_RIGHT;
			tmpAA.angles.psi = ALPHA_HLX_PSI_RIGHT;

			for(ic=0; ic<newHlxLength; ic++)
				tmpProt.AAs.push_back(tmpAA);
			//generate structure from a list of torsion angles
			tmpProt.torsion2coord();

			//set SS type flag
			short cntr=0;
			for (kk=0; kk<mergedSticks[j].size(); kk++){
				for (ic=startIndx[mergedSticks[j][kk]]; ic<=endIndx[mergedSticks[j][kk]]; ic++){
					tmpProt.AAs[cntr++].SStype = 'H';
				}
				if (loopNaa[mergedSticks[j][kk]]>0 && kk!= mergedSticks[j].size ()-1){
					for (ic=0; ic<loopNaa[mergedSticks[j][kk]]; ic++)
						tmpProt.AAs[cntr++].SStype = 'L';
				}
			}
			//or you can simply set all AA types to H.
			//for (ic=0; ic<newHlxLength; ic++)
			//	tmpProt.AAs[ic].SStype = 'H';

			Protein axisHlx;
			axisHlx = mergeTwoHlces(axis, tmpProt, nRise, newHlxLength);
            //axisHlx.writePDB("axisHelix"+toString(stickIndx[mergedSticks[j][0]]+1)+".pdb", 1, axisHlx.numOfAA());

			skeleton[mergedSticks[j][0]] = axisHlx;
			newEdges[mergedSticks[j][0]] = axis;
		}
		//cout<<"Done.."<<endl;
		//merge helices
		short nDeleted=0;
		for (j=0; j<mergedSticks.size(); j++){
		    //cout<<j<<" "<<mergedSticks[j][0]<<" "<<mergedSticks[j][mergedSticks[j].size()-1]-nDeleted<<endl;
		    //update information of the first hlx
		    endIndx[mergedSticks[j][0]-nDeleted] = endIndx[mergedSticks[j][mergedSticks[j].size()-1]-nDeleted];
		    if (mergedSticks[j][mergedSticks[j].size()-1]-nDeleted < nSticks-1){
                trace[mergedSticks[j][0]-nDeleted] = trace[mergedSticks[j][mergedSticks[j].size()-1]-nDeleted];
                traceNaa[mergedSticks[j][0]-nDeleted] = traceNaa[mergedSticks[j][mergedSticks[j].size()-1]-nDeleted];
                loopNaa[mergedSticks[j][0]-nDeleted] = loopNaa[mergedSticks[j][mergedSticks[j].size()-1]-nDeleted];
                tracePnts[mergedSticks[j][0]-nDeleted] = tracePnts[mergedSticks[j][mergedSticks[j].size()-1]-nDeleted];
                traceNumPnts[mergedSticks[j][0]-nDeleted] = traceNumPnts[mergedSticks[j][mergedSticks[j].size()-1]-nDeleted];
		    }
		    for (ic=1;ic<mergedSticks[j].size(); ic++){
		        //cout<<"  ic= "<<ic<<" indx to delete : "<<mergedSticks[j][ic]-nDeleted<<endl;
		        stickIndx.erase(stickIndx.begin()+mergedSticks[j][ic]-nDeleted);
		        hlxIndx.erase(hlxIndx.begin()+mergedSticks[j][ic]-nDeleted);
		        startIndx.erase(startIndx.begin()+mergedSticks[j][ic]-nDeleted);
		        endIndx.erase(endIndx.begin()+mergedSticks[j][ic]-nDeleted);
		        stickDirection.erase(stickDirection.begin()+mergedSticks[j][ic]-nDeleted);
		        skeleton.erase(skeleton.begin()+mergedSticks[j][ic]-nDeleted);
				newEdges.erase (newEdges.begin () +mergedSticks[j][ic]-nDeleted);
		        if (mergedSticks[j][ic]<nSticks-1){
                    loopNaa.erase(loopNaa.begin()+mergedSticks[j][ic]-nDeleted);
                    trace.erase(trace.begin()+mergedSticks[j][ic]-nDeleted);
                    traceNaa.erase(traceNaa.begin()+mergedSticks[j][ic]-nDeleted);
                    tracePnts.erase(tracePnts.begin()+mergedSticks[j][ic]-nDeleted);
                    traceNumPnts.erase(traceNumPnts.begin()+mergedSticks[j][ic]-nDeleted);
		        }
                nDeleted++;
		    }
		}

		//calculate the RMSD between SSE only
		int *sIndeces;
		sIndeces = new int [stickIndx.size()];
		for(i=0;i<stickIndx.size();i++)
			sIndeces[i] = startIndx[i];

		vector<Protein> SSE = skeleton;
		rmsd = getRMSDSSE(SSE,pdb,sIndeces);

		
		//write SSE before modeling the loop
		//Protein allSSE;
		//rename AAs
		for (i=0;i<stickIndx.size();i++){
			for (j=0;j<SSE[i].numOfAA();j++){
				SSE[i].renameAA(j, pdb.AAs[startIndx[i]+j].chr3);
				SSE[i].AAs[j].num = pdb.AAs[startIndx[i]+j].num;
			}
		}

		allSSE = SSE[0];
		for (i=1;i<stickIndx.size();i++){
			allSSE.append(SSE[i],0,SSE[i].numOfAA()-1, SSE[i].AAs[0].num - SSE[i-1].AAs[SSE[i-1].numOfAA()-1].num-1);
		}

        /*
		for (j=0; j<stickIndx.size()-1; j++){
			cout<<"hlx# "<<hlxIndx[j]<<" to Stck# "<<stickIndx[j]+1<<" and hlx# "<<hlxIndx[j+1]<<" to stck# "<<stickIndx[j+1]+1<<endl;
			cout<<"   stck# "<<stickIndx[j]+1<<" start at AA indx "<<startIndx[j]<<" and ends at "<<endIndx[j]<<endl;
			cout<<"   stck# "<<stickIndx[j+1]+1<<" starts at AA indx"<<startIndx[j+1]<<" and ends at "<<endIndx[j+1]<<endl;
			cout<<"   stick 1 nAA = "<<skeleton[j].numOfAA()<<" ("<<endIndx[j]-startIndx[j]+1<<") stick 2 nAA= "<<skeleton[j+1].numOfAA()<<" ("<<endIndx[j+1]-startIndx[j+1]+1<<")"<<endl;
			cout<<"   Trace found : "<<trace[j][0]<<"-"<<trace[j][1]<<" expected nAA= "<<traceNaa[j]<<endl;
			cout<<"   loop nAA : "<<loopNaa[j]<<endl;
		}
		*/

        //getchar();getchar();
		//model loops
		vector<vector<Coordinate> > traceToAvoid,orgTraceToAvoid;								//the list of points represents the part of the topology built so far o use for collision detection
		fullModel.append(skeleton[0], 0, skeleton[0].numOfAA()-1);

		//add the traces of hlces to avoid during the building of loop
		for (j=0; j<stickIndx.size (); j++){
			vector<Coordinate> tmp;
			//do not add first and last segments
			for (ic=1;ic<newEdges[j].size ()-2; ic++){
				tmp.push_back (newEdges[j][ic]);
			}
			//cout<<stickIndx[j]+1<<" ic= "<<ic<<" size= "<<newEdges[j].size ()<<endl;
			//for last segment ... the ength maynot be nRise*ALPHA_RISE, so, we may need to ignore the last two
			if (getDistance(newEdges[j][ic], newEdges[j][ic+1])> nRise*ALPHA_RISE-1 || newEdges[j].size ()<5)
				tmp.push_back (newEdges[j][ic]);

			traceToAvoid.push_back (tmp);

            /*
			//write the new edge
			Protein curvd;
			curvd = points2pdb(newEdges[j],"GLY", " CA ");
			curvd.writePDB("newStick_"+toString(stickIndx[j]+1)+".pdb",1 , curvd.numOfAA());

			curvd = points2pdb(tmp, "GLY", " CA ");
			curvd.writePDB("collTrace"+toString(stickIndx[j]+1)+".pdb",1 , curvd.numOfAA());
			*/
		}
		orgTraceToAvoid = traceToAvoid;

		//cout<<"fullModel number of AAs = "<<fullModel.numOfAA()<<endl;
		for (j=0; j<stickIndx.size()-1; j++){
           vector<char> aaTypes;
            for (ic=endIndx[j]+1; ic<startIndx[j+1]; ic++){
                aaTypes.push_back(pdb.AAs[ic].SStype);
            }

			//for the next stick..add a point close to the beginning so the loop does not collide at the beginning of the stick
			traceToAvoid[1].insert (traceToAvoid[1].begin (), pointOnLine(newEdges[j+1][1], newEdges[j+1][0], (nRise*ALPHA_RISE)/4));

		    if (trace[j][0]!=0){ //} && abs(traceNaa[j]-loopNaa[j])<traceNaa[j]*0.5){		//the length b/w trace and sequence is very different

				//resolve the collision b/w the trace and hlces sticks
				traceSticksCollision(orgTraceToAvoid, tracePnts[j], stickIndx.size ());

                //write the new edge 
                //Protein curvd;
                //curvd = points2pdb(tracePnts[j],"GLY", " CA ");
                //curvd.writePDB("loopTrace_"+toString(j+1)+".pdb",1 , curvd.numOfAA());
                

				fullModel = buildStructureOnTrace(fullModel, skeleton[j+1], 1, tracePnts[j], traceToAvoid, aaTypes);
		    }
		    else{

				//cout<<"modeled randomly..."<<endl;

				randomFBCCD(fullModel, skeleton[j+1], 1, aaTypes, traceToAvoid, 0.05);

		    }
			traceToAvoid.erase(traceToAvoid.begin ());

            //fullModel.writePDB("currentLoopTrace.pdb", 1, fullModel.numOfAA());

            //getchar();getchar();
        }

		//fullModel.writePDB("fullModel.pdb", 1 , fullModel.numOfAA());
        //getchar(); getchar();
    //}

    //cout<<"assign sequence to the structure... ";
     //assign sequence to the model
	int iiAA = sIndx = startIndx[0];
	for (i=0 ; i<fullModel.numOfAA(); i++)
	{
		fullModel.renameAA(i, pdb.AAs[iiAA].chr3);
		iiAA++;
	}
	//cout<<"Done. "<<iAA<<endl;

    return fullModel;

}

#endif
