#ifndef	TOPOLOGY_GRAPH_H
#define TOPOLOGY_GRAPH_H

#include <vector>
#include <bitset>
#include <queue>

#define _WIN

#ifdef _WIN
#include <time.h>					//used to find elapsed time
#else
#include <sys/time.h>
#endif

#include "skeleton_overall.h"
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
	//pair<short, short>  *pCandidate;
	short *pCandidateRow;
	short *pCandidateCol;
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
void gShortestRec(vector<vector<cell> > graph,  short sRow, short sCol, short eRow, short eCol, int K);	//find shortest path using a recurrence
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
	if ((nAASeq * parcent > nAAstick) ||		\
		(nAAstick * parcent > nAASeq))
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
		//cout<<setw(cntr-1)<<" "<<"Clique is found.. {";
		//for (int i=0; i<R.size (); i++)
		//	cout<<R[i]+1<<" ";
		//cout<<"}"<<endl;
		clx.push_back(R);
	}
	//choose a pivot node from P U X
	//I'll choose the node with maximum number of neighbors
	short u= -1;		//the index of the pivot node
	short maxN = 0;		//maximum number of neighbors
	short i,j,k;

	/*
	//print P
	cout<<setw(cntr)<<" "<<"P = {";
	for (i=0; i<P.size (); i++)
		cout<<P[i]+1<<" ";
	cout<<"}"<<endl;
	//Print R
	cout<<setw(cntr)<<" "<<"R = {";
	for (i=0;i<R.size ();i++)
		cout<<R[i]+1<<" ";
	cout<<"}"<<endl;
	//print X
	cout<<setw(cntr)<<" "<<"X = {";
	for (i=0; i<X.size (); i++)
		cout<<X[i]+1<<" ";
	cout<<"}"<<endl;
	*/

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
		//cout<<setw(cntr)<<" "<<i<<" : working on cluster : "<<P[i]+1<<endl;
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
			/*
			cout<<setw(cntr)<<" "<<"P intersect N("<<v+1<<")= {";
			for (j=0; j<tmpP.size (); j++)
				cout<<tmpP[j]+1<<" ";
			cout<<"}"<<endl;
			cout<<setw(cntr)<<" "<<"X intersect N("<<v+1<<")= {";
			for (j=0; j<tmpX.size ();j++)
				cout<<tmpX[j]+1<<" ";
			cout<<"}"<<endl;
			*/

			BronKerbosch(adjMtrx, clx, tmpR, tmpP, tmpX);

			//update P and X
			X.push_back (v);
			//--cntr;

			/*
			cout<<setw(cntr)<<" "<<"P is : {";
			for (int m=0; m<P.size (); m++)
				cout<<P[m]+1<<" ";
			cout<<"}"<<endl;
			cout<<"  deleting... "<<P[i]+1;
			*/
			P.erase (P.begin () + i);
			i--;
			
			//cout<<" Done.. i="<<i<<endl;
		}
		i++;
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

	//initiate sets
	vector<short>	P,			//the set contains vertices
					X,			//the NOT set
					R;			//the set will contain the maximum clique
	for (i=0; i<adjMtrx.size(); i++)
		P.push_back(i);

	//find cliques and save the indices of cluster form cliques in clx
	BronKerbosch(adjMtrx, clx, R,P,X);

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
	
	cout<<"Number of cliques found = "<<clx.size ()<<endl;
	for (i=0; i<clx.size (); i++){
		cout<<"Clique# "<<i+1<<" : ";
		for (int j=0; j<clx[i].size (); j++)
			cout<<clx[i][j]+1<<" ";
		cout<<endl;
	}
	
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

			cout<<"clqSze= "<<clqSze<<endl;

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
			cout<<"  Done.."<<endl;

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
		cout<<"# of cliques = "<<nClxArr<<"  Deleting duplicate cliques... ";

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
		cout<<"Done."<<endl;
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
	cout<<"# of cliques = "<<clx.size()<<" Deleting completely shared cliques... ";
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
	cout<<"Done."<<endl;

	//delete empty clusters
	cout<<"Deleting empty cliques... ";
	i=0;
	while (i<clx.size ()){
		if (clx[i].size () == 0){
			clx.erase (clx.begin () + i);
			i--;
		}
		i++;
	}
	cout<<"Done."<<endl;

	//print valid cliques
	
	cout<<"AFTER DELETING SHARED CLUSTERS.."<<endl;
	for (i=0; i<clx.size (); i++){
		cout<<"Clq# "<<i+1<<" : ";
		for (j=0; j<clx[i].size (); j++){
			cout<<clx[i][j]+1<<" ";
		}
		cout<<endl;
	}
	

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
					vector<vector<Coordinate> > &clusters,		//coordinate of points form clusters
					pair<short,float> **sDist,					//matrix has shortest distance b/w each pair of clusters and the indeces of points with shortest distances
					short startPnt,								//the starting point where you want to start ur tracing
					short nClusters,							//number of clusters without clusters for SS ends
					vector<vector<short> > &paths,				//trace paths found for density map
					vector<float> &pathsCosts,					//parallel vector with paths contains the cost for each trace path
					float maxLngth = 500.0){					//maximum length of path allowed

	vector<pair<short, short> > open;			//open set represetns clusters that are visited already and the number of children in the stack
	vector<float>				openLength;		//the length of the trace so far for each node in open
	vector<pair<short, short> > stack;			//clusters being processed and the parent for each cluster
	vector<bool> visitFlag (adjMtrx.size(), false);
	pair<short, short> tmp, top;
	vector<short> onePath;
	float pathLength = 0.0;
	short i, j;


	//push the starting point (the indx of starting point) to the stack
	tmp.first = startPnt;
	tmp.second = -1;				//no parent
	stack.push_back (tmp);

	while (!stack.empty ()){
		//pop the top of the stack and put its children (those clusters are satisfying threshold distance)
		top = stack[stack.size ()-1];

		cout<<top.first +1;			//cluster being processed
			
		//decrement number of neighbors in the stack for its parent
		if (top.second != -1){		//if it is not the start pnt
			open[top.second].second--;
			cout<<"  PathLength= "<<openLength[openLength.size ()-1];
			cout<<"     ( parent: "<<open[top.second].first+1<<"  p.Indx= "<<top.second<<"      #children: "<<open[top.second].second<<" ).";
		}

		cout<<endl;

		stack.pop_back ();

		//initialize number of children and push to open
		top.second = 0;
		open.push_back (top);
		visitFlag[top.first] = true;		//set visit flag

		//calculate the length of the path so far
		if (open.size () > 1){
			//get distance b/w last two clusters in open list
			short	openLastIndx = open.size () -1,
					addedClusterIndx = open[openLastIndx].first,
					centroid1Indx = clusters[addedClusterIndx].size()-1,

					clusterIndx = open[openLastIndx-1].first,
					centroid2Indx = clusters[clusterIndx].size()-1,

					point1Indx = sDist[addedClusterIndx][clusterIndx].first,						//the index of closest point to the other cluster
					point2Indx = sDist[clusterIndx][addedClusterIndx].first;

			//cout<<"addedCluster : "<<addedClusterIndx<<" centroid : "<<centroid1Indx<<" pnt1 : "<<point1Indx<<"  clustIndx : "<<clusterIndx<<" centroid2 : "<<centroid2Indx<<" pnt2 : "<<point2Indx<<endl;
		//	pathLength =	openLength[openLastIndx-1] + 
		//					getDistance(clusters[addedClusterIndx][centroid1Indx], clusters[addedClusterIndx][point1Indx]) +		//distance b/w cluster1 centroid and the closest point to cluster2	
		//					getDistance(clusters[addedClusterIndx][point1Indx], clusters[clusterIndx][point2Indx]) +				//distance b/w closest points in both clusters						
		//					getDistance(clusters[clusterIndx][point2Indx], clusters[clusterIndx][centroid2Indx]);					//distance b/w closest point in cluster#2 and its centroid

							//distance b/w centroids of the last two clusters in open list
			pathLength =	openLength[openLastIndx-1] + getDistance(clusters[addedClusterIndx][centroid1Indx], clusters[clusterIndx][centroid2Indx]);

			openLength.push_back (pathLength);
		}
		else
			openLength.push_back (0.0);				//path length for first node
		
		//don't continue with this node if you already exceeded the max length allowed
		if (pathLength > maxLngth && top.first < nClusters){
			cout<<" The Path exceeds the maximum length allowd which is : "<<maxLngth<<endl;
			visitFlag[open[open.size ()-1].first] = false;
			open.pop_back ();
			openLength.pop_back ();
			continue;
		}

		//if the cluster we r working on is one of SS ends. then do not generate its children b/s this would be a trace end
		if ( top.first < nClusters || top.first == startPnt){
			//generate children (neighbors) for the cluster
			for (i=0 ; i<adjMtrx[top.first].size(); i++){
				cout<<"  "<<adjMtrx[top.first][i]+1<<" ";

				if (!visitFlag[adjMtrx[top.first][i]]){
					tmp.first  = adjMtrx[top.first][i];
					tmp.second = open.size ()-1;			//save the index of the parent in open set
					stack.push_back (tmp);					//push neighbor to stack
					open[tmp.second].second++;				//increment number of neighbors in stack
					cout<<"pushed to stack.   parentIndx= "<<tmp.second<<" ( "<<open[tmp.second].first+1<<" )."<<endl;
				}
				else
					cout<<" found in open..."<<endl;
			}
		}

		
		//if we reach a leaf node...with no neighbor in the stack
		if (open[open.size ()-1].second<1){
			cout<<"  a leaf node has been reached "<<open[open.size ()-1].first+1<<" .... "<<endl;
			onePath.clear();
			for (i=0; i<open.size ()-1; i++)
				onePath.push_back(open[i].first);
			//push last cluster
			onePath.push_back (open[i].first );

			//save the path and its cost
			paths.push_back (onePath);
			pathsCosts.push_back (openLength[i]);

			//delete all nodes in open set with 0 neighbors
			i=open.size ()-1;
			while (i>=0 && open[i].second <=0){	
				cout<<"       "<<open[i].first +1<<" deleted from open"<<endl;
				visitFlag[open[i].first] = false;
				open.pop_back ();
				openLength.pop_back ();
				i--;
			}
			cout<<"    Open contents"<<endl;
			for (i=0; i<open.size ();i++)
				cout<<"("<<open[i].first+1 <<","<<open[i].second <<") ";
			cout<<endl;
		}
		else{
			//save a new path if we reach a stick end
			if (top.first>= nClusters && top.first != startPnt){
				onePath.clear ();
				cout<<"  find stick end "<<top.first+1<<" path saved"<<endl;
				for (i=0; i<open.size ()-1; i++){
					onePath.push_back(open[i].first);
				}

				//push last cluster
				onePath.push_back(open[i].first);

				//save the path and its cost
				paths.push_back (onePath);
				pathsCosts.push_back (openLength[i]);
			}
		}

		//getchar();
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void buildAdjMtrx(vector<vector<Coordinate> > &clusters, pair<short, float> **(&cDist), vector<vector<short> > &adjMtrx, float apix, float continuatyTHR){

	int i,j,k,l;
	int nClusters = clusters.size ();

	Coordinate cPnt, nPnt;			//the two points we are measuring distance b/w...in MAP indexing system

	cout<<"Allocate dynamic memory... of size "<<nClusters<<" X "<<nClusters<<".... "<<endl;

	//initiate cluster distance matrix....add also distance b/w SS ends and clusters
	cDist = AllocateDynamicArray<pair<short,float> > (nClusters, nClusters);		
	
	for (i=0; i<nClusters; i++){
		for (j=0; j<nClusters; j++){
			cDist[i][j].second = 9999.0;
			cDist[i][j].first = -1;
		}
			
	}

	cout<<"Done.."<<endl;

	/*
	 *		find minimum distance b/w each pair of clusters
	 */
	float dist;
	short p1Indx = -1, p2Indx;			//the indeces of the points have the minimum distance with the other cluster
	for (i=0; i<nClusters-1; i++){
		//cout<<i+1<<" : ";
		for (k= i+1; k<nClusters; k++){
			//cout<<"    :  "<<k+1<<endl;
			p1Indx = -1;
			p2Indx = -1;
			float minDist = 999999.0;
			for (j=0; j<clusters[i].size (); j++){
				cPnt = clusters[i][j];
				for (l=0; l<clusters[k].size (); l++){
					nPnt = clusters[k][l];
					dist  = getDistance(cPnt, nPnt)/apix;
					if (dist < minDist){
						minDist = dist;
						p1Indx = j;
						p2Indx = l;
					}
				}
			}
			cDist[i][k].second = cDist[k][i].second = minDist;
			cDist[i][k].first = p1Indx;
			cDist[k][i].first = p2Indx;


			//cout<<cDist[i][k].first<<" - "<<cDist[k][i].first<<" dist= "<<cDist[i][k].second<<" == "<<cDist[k][i].second<<endl;			
		}
		//cout<<endl;
	}

	cout<<"shortest distances is done.."<<endl;
	/*
	 *	create adjacency Matrix
	 */
	adjMtrx.clear ();
	adjMtrx.resize (nClusters);			//the size is : the initial clusters and the clusters for SSs' ends
	for (i=0 ; i<nClusters; i++){
		cout<<"Neigbors(cluster# "<<i+1<<")  : ";
		for (j=0; j<nClusters; j++){
			if (i!=j){
				if (cDist[i][j].second <= continuatyTHR){
					cout<<j+1<<"  ";
					adjMtrx[i].push_back (j);
				}
			}
		}
		cout<<endl;
	}
	cout<<"build adjacency matrix is done.."<<endl;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//select the best trace path among all paths b/w two SS ends
// this function assumes that all paths of a particular end are follow each other
float getBestPathTrace(vector<vector<Coordinate> > &clusters, vector<vector<short> > & tracePaths, vector<float> &traceLengths, short nClusters, float loopLength, int firstEnd, int secondEnd, short &bestTraceIndx1, short &bestTraceIndx2){
	int i,j,k;
	//get the starting and ending indeces of the paths belong to first End
	int firstIndx1 = -1, firstIndx2=-1, secondIndx1=-1, secondIndx2=-1;

	//always tracePaths[i][0] is the starting end
	i=0;
	while (i<tracePaths.size () && tracePaths[i][0] != firstEnd) i++;
	firstIndx1 = i;
	while (i<tracePaths.size () && tracePaths[i][0] == firstEnd) i++;
	firstIndx2 = i-1;

	i=0;
	while (i<tracePaths.size () && tracePaths[i][0] != secondEnd) i++;
	secondIndx1 = i;
	while (i<tracePaths.size () && tracePaths[i][0] == secondEnd) i++;
	secondIndx2 = i-1;

//	if (firstEnd == 78 && secondEnd == 82){
///		cout<<"["<<firstEnd+1<<" "<<secondEnd+1<<"] fIndx1= "<<firstIndx1<<" fIndx2= "<<firstIndx2<<" sIndx1= "<<secondIndx1<<" sIndx2= "<<secondIndx2<<" Lplngth= "<<loopLength<<endl;
//	}


	float bestWeight = fabs(getDistance(clusters[firstEnd][clusters[firstEnd].size ()-1], clusters[secondEnd][clusters[secondEnd].size ()-1]) - loopLength) + 15.0;		//add some weights to penalty of missing trace in between

	int lastCluster1Indx, lastCluster2Indx, centroid1Indx, centroid2Indx;
	Coordinate trace1LastPoint, trace2LastPoint;			//the coordinate of last point in both traces

	//cout<<"  initial bestWeight= "<<bestWeight<<endl;

/*
	//work only on full traces
	for (i=firstIndx1; i<=firstIndx2; i++){
		if (tracePaths[i][tracePaths[i].size ()-1] == secondEnd){

			if (fabs(loopLength-traceLengths[i]) < bestWeight){
				bestWeight = fabs(loopLength - traceLengths[i]);
				bestTraceIndx1 = i;
				bestTraceIndx2 = i;
			}
		}
	}
*/
	
	for (i=firstIndx1; i<=firstIndx2; i++){
		float matchingLength = 9999.0;
		if (tracePaths[i][tracePaths[i].size ()-1] == secondEnd){
			//cout<<"Found completer path  trace Length= "<<traceLengths[i]<<"   and loop length= "<<loopLength<<endl;
			if (traceLengths[i] <= loopLength){
				if (fabs(loopLength-traceLengths[i]) < bestWeight){
					bestWeight = fabs(loopLength - traceLengths[i]);
					bestTraceIndx1 = i;
					bestTraceIndx2 = i;
				}
				//bestWeight -= 10.0;					//reward a connected full path
			}
		}
		else{
			float THR = 15.0;			//discontinuty threshold distance (for real maps use 15. simulated maps use 8.0)
			lastCluster1Indx	= tracePaths[i][tracePaths[i].size ()-1];
			centroid1Indx		= clusters[lastCluster1Indx].size ()-1;
			trace1LastPoint		= clusters[lastCluster1Indx][centroid1Indx];

			for (j=secondIndx1; j<=secondIndx2; j++){
				
				//cout<<"checking path "<<i+1<<" with Path "<<j+1<<endl;

				matchingLength = traceLengths[i];
				if ((tracePaths[j][tracePaths[j].size ()-1] < nClusters) &&		//  || tracePaths[j].size () == 1) &&
					(tracePaths[i][tracePaths[i].size ()-1] < nClusters)){		//  || tracePaths[i].size () == 1)){

					//if (tracePaths[j].size () == 1 || tracePaths[i].size () == 1){
					//	THR = 8.0;
					//	cout<<"path # "<<i+1<<" and path # "<<j+1<<" gapDist= "<<getDistance(trace1LastPoint, clusters[tracePaths[j][tracePaths[j].size ()-1]][clusters[tracePaths[j][tracePaths[j].size ()-1]].size ()-1])<<endl;
					//	getchar();getchar();
					//}
					//else
					//	THR = 15.0;

					//deal with discontinue traces
					matchingLength += traceLengths[j];

					//cout<<"  matchingLength now = "<<matchingLength<<endl;

					if (matchingLength <= loopLength){

						lastCluster2Indx	= tracePaths[j][tracePaths[j].size ()-1];
						centroid2Indx		= clusters[lastCluster2Indx].size ()-1;
						trace2LastPoint		= clusters[lastCluster2Indx][centroid2Indx];

						//if (tracePaths[j].size () == 1 || tracePaths[i].size () == 1){
						//	cout<<" gapDist= "<<getDistance(trace1LastPoint, trace2LastPoint);
						//}

						if (getDistance(trace1LastPoint, trace2LastPoint) < THR){
							matchingLength		+= getDistance(trace1LastPoint, trace2LastPoint);
							
							//cout<<"  gap Distance= "<<getDistance(trace1LastPoint, trace2LastPoint)
							//	<<" P1: "<<trace1LastPoint.x<<" "<<trace1LastPoint.y<<" "<<trace1LastPoint.z<<" P2: "<<trace2LastPoint.x<<" "<<trace2LastPoint.y<<" "<<trace2LastPoint.z<<endl;

							//if (firstEnd == 78 && secondEnd == 82){
							//	cout<<"checking path "<<i+1<<" with Path "<<j+1<<endl;
							//	cout<<"  gap Distance= "<<getDistance(trace1LastPoint, trace2LastPoint)
							//		<<" P1: "<<trace1LastPoint.x<<" "<<trace1LastPoint.y<<" "<<trace1LastPoint.z<<" P2: "<<trace2LastPoint.x<<" "<<trace2LastPoint.y<<" "<<trace2LastPoint.z<<endl;
							//	getchar(); 
							//	getchar();
							//}
							//if (matchingLength > loopLength && (tracePaths[j].size () == 1 || tracePaths[i].size () == 1)){
							//	cout<<"  rejected. lngths ("<<matchingLength<<" "<<loopLength<<")";
							//}

							if (matchingLength <= loopLength){
								if (fabs(loopLength-matchingLength) < bestWeight){
									bestWeight = fabs(loopLength-matchingLength);
									bestTraceIndx1 = i;
									bestTraceIndx2 = j;
								}
								//if (tracePaths[j].size () == 1 || tracePaths[i].size () == 1)
								//	cout<<" lengths ("<<matchingLength<<" "<<loopLength<<")"<<" bst= "<<bestWeight;
							}
						}
					}
					//if (tracePaths[j].size () == 1 || tracePaths[i].size () == 1)
					//	cout<<endl;
				}
				//if (firstEnd == 78 && (secondEnd == 82)) { getchar(); getchar();}
			}
			//if (firstEnd == 78 && secondEnd == 82) getchar();
		}
	}

//	cout<<"    final bestWeight = "<<bestWeight<<endl;


	return bestWeight;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 *		Find all traces b/w any two SSends...then select the best trace fit the number of AAs in the loop b/w the 2 SSends according to the link in the graph
 */
void setLoopsWeights(vector<vector<cell> > & graph, Map inMRC, vector<vector<Coordinate> > ssEdges, vector<SecondaryStruct>	seqSS, vector<Coordinate> pnts, float peakTHRg, string outPath){
	int i, j, k, l;

	float	distTHR = 15.0,					//thresold distance b/w two ends of SSs
			continuatyTHR = 1.1 * inMRC.apixX;			//a threshold used as a cuttoff distance b/w any two points to consider them continued (non-disjoint)

	vector<vector<Coordinate> > clusters, newClusters;
	short nClusters;										//the number of initial clusters (without clusters for SS ends)	

	pair<short,float> **cDist;						//cluster shortest distances (indeces and distance)
	

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


	peakClustering(pnts, clusters, ssEdges[0][0], 2.0*inMRC.apixX); //start from a random end...which is here the first SS and the first end

	nClusters = clusters.size ();
	
	
	//print clusters (all points or centriod only)
	for (i=0; i<nClusters; i++){
		Protein trace;
		//write local peaks into a pdb file
		trace = points2pdb(clusters[i], "HOH", " O  ");
		trace.writePDB(outPath +"_Cluster_" + toString(i+1) +".pdb", 1, trace.numOfAA());							//all points
		//trace.writePDB(outPath +"_Cluster_" + toString(i+1) +".pdb", trace.numOfAA(), trace.numOfAA());				//centroid only
	}
	
	/*
	 *		find shortest distances b/w each pair of clusters and then build adjMtrx
	 */
	buildAdjMtrx(clusters, cDist, adjMtrx, inMRC.apixX, continuatyTHR);

	/*
	 *	find all cliques in the graph....using Bron-Kerbosch algorithm
	 */

	cout<<"Finding Cliques....";
	findallCliques(clusters, adjMtrx, clx, clxAdjMtrx, clxClusters);
	cout<<" Done"<<endl;

	//print final cliques and their associate local peaks clusters
	for (i=0; i<clx.size (); i++){
		cout<<"Clique# "<<i+1<<"  : ";
		for (j=0; j<clx[i].size(); j++){
			cout<<clx[i][j]+1<<" ";
		}
		cout<<endl;
	}

	//for SS ends....create a cluster of length = 1 contains one end point for each SS
	newClusters = clxClusters;
	clusters.clear ();
	vector<Coordinate> tmpCluster;
	for (i=0; i<ssEdges.size (); i++){

		tmpCluster.clear ();
		tmpCluster.push_back (ssEdges[i][0]);
		newClusters.push_back (tmpCluster);


		tmpCluster.clear ();
		tmpCluster.push_back (ssEdges[i][ssEdges[i].size ()-1]);
		newClusters.push_back (tmpCluster);
	}

	//print cliques (all Points all centroid only)
	for (i=0; i<newClusters.size(); i++){
		Protein trace;

		//write local peaks into a pdb file
		trace = points2pdb(newClusters[i], "HOH", " O  ");
		trace.writePDB(outPath +"_Clique_" + toString(i+1) +".pdb", 1, trace.numOfAA());						//all points
		//trace.writePDB(outPath +"_Clique_" + toString(i+1) +".pdb", trace.numOfAA(), trace.numOfAA());			//centroid only
	}

//	cout<<"Press any key.."<<endl;
//	getchar();
//	getchar();

	nClusters = clx.size ();


	//re-calculate Adjacency Mtrx and short distances matrix for cliques
	buildAdjMtrx(newClusters, cDist, clxAdjMtrx, inMRC.apixX, continuatyTHR);

	/*
	//update clxAdjMtrx ... by adding ssEdges ends
	clxAdjMtrx.resize (nClusters+2*ssEdges.size ());
	for (i=nClusters; i<newClusters.size (); i++){
	//	cout<<"Clq# "<<i+1<<endl;
		for (j=0; j<nClusters; j++){
		//	cout<<" with clq# "<<j+1;
			//find minimum distance between the two clusters
			for (k=0; k<newClusters[j].size (); k++){
				float dist = getDistance(newClusters[j][k], newClusters[i][0])/inMRC.apixX;
				if (dist <= continuatyTHR){
					clxAdjMtrx[i].push_back (j);
					clxAdjMtrx[j].push_back (i);
				//	cout<<" close.";
					break;
				}
			}
			//cout<<endl;
		}
	}
	cout<<"CLX_ADJ_MTRX"<<endl;
	for(i=0; i<clxAdjMtrx.size (); i++){
		cout<<"Clx# "<<i+1<<"  : ";
		for (j=0; j<clxAdjMtrx[i].size (); j++){
			cout<<clxAdjMtrx[i][j]+1<<" ";
		}
		cout<<endl;
	}
	*/
	
	getchar();
	getchar();

	/*
	 *		find trace paths from density map "Skeleton" or "Local Peaks"
	 */

	Protein pathTrace;

	tracePaths.clear ();
	traceLengths.clear ();
	for (i=0; i<2*ssEdges.size (); i+=2){
		
		//determine the maximum length of a loop could get out from this end
		maxLength = 0.0;
		for (j=0; j<graph.size (); j++){
			for (k=0; k<graph[j][i].nAAloop .size (); k++){
				if (graph[j][i].nAAloop[k] > maxLength)
					maxLength = graph[j][i].nAAloop[k];
			}
		}
		
		//find trace paths for first end point
		findTracePaths(clxAdjMtrx, newClusters, cDist, nClusters + i, nClusters, tracePaths, traceLengths, maxLength*3.8);

		//cout<<"Number of paths = "<<tracePaths.size ()<<endl;
		/*
		for (j=0; j<tracePaths.size (); j++){

			pathTrace.initialize();
			cout<<"Path# "<<j+1<<": ";
			for (k=0; k<tracePaths[j].size (); k++){
				cout<<"  "<<tracePaths[j][k]+1;
				trace = points2pdb(newClusters[tracePaths[j][k]], "GLY", " CA ");
				pathTrace.append(trace, trace.numOfAA()-1, trace.numOfAA()-1);						//centroid only							
			}
			cout<<"  length= "<<pathsCosts[j]<<endl;
			pathTrace.writePDB(outPath+"_pathTrace_"+toString(j+1)+".pdb", 1, pathTrace.numOfAA());
			//getchar();
		}
		*/
		

		//find tracePaths for the other end
		//tracePaths.clear ();
		//pathsCosts.clear ();
		pathTrace.initialize();

		//determine the maximum length of a loop could get out from this end
		maxLength = 0.0;
		for (j=0; j<graph.size (); j++){
			for (k=0; k<graph[j][i+1].nAAloop .size (); k++){
				if (graph[j][i].nAAloop[k] > maxLength)
					maxLength = graph[j][i].nAAloop[k];
			}
		}

		//find trace paths for second end point
		findTracePaths(clxAdjMtrx, newClusters, cDist, nClusters + i + 1, nClusters, tracePaths, traceLengths, maxLength*3.8);

		//cout<<"Number of paths = "<<tracePaths.size ()<<endl;
		/*
		for (j=0; j<tracePaths.size (); j++){
			pathTrace.initialize();

			cout<<"Path# "<<j+1<<": ";
			for (k=0; k<tracePaths[j].size (); k++){
				cout<<"  "<<tracePaths[j][k]+1;
				trace = points2pdb(newClusters[tracePaths[j][k]], "GLY", " CA ");
				pathTrace.append(trace, trace.numOfAA()-1, trace.numOfAA()-1);
			}
			cout<<" length= "<<pathsCosts[j]<<endl;
			pathTrace.writePDB(outPath+"_pathTrace_"+toString(j+1)+".pdb", 1, pathTrace.numOfAA());
			getchar();
		}
		*/

		//cout<<"Press any key.."<<endl;
		//getchar();
		//getchar();
	}

	//delete traces that have SS end in the middle
	/*
	i=0;
	while (i<tracePaths.size ()){
		for (j=1; j<=tracePaths[i].size ()-1; j++){
			if (tracePaths[i][j] >= nClusters){
				//delete trace
				tracePaths.erase (tracePaths.begin () + i);
				traceLengths.erase (traceLengths.begin () + i);
				i--;
				break;
			}
			j++;
		}
		i++;
	}
	*/

	cout<<"Number of Total Paths : "<<tracePaths.size()<<" and "<<traceLengths.size()<<endl;

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
	vector<short> newTrace;
	vector<vector<short> > newTraceSet;
	float newLength;
	vector<float> newTraceLengths;
	float ssGapTHR = 8.0;						//the maximum distance b/w the SS and the first cluster in the path to consider
	float gapDist = 9999.0;
	for (i=0; i<graph[0].size (); i++){

		//find how many traces out from a particular SS end
		short tracesCntr = 0;
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

						cout<<"I found a new path b/w graph col "<<i<<" and path# "<<j+1<< " old length = "<<traceLengths[j]<<" new Length= "<<newLength<<" gap= "<<gapDist<<endl;
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

	//getchar();getchar();
	
	
	

	//write paths into viewable files by Chimera
	for (j=0; j<tracePaths.size (); j++){
		cout<<"Path# "<<j+1<<": ";
		pathTrace.initialize();
		pathTrace.header.push_back("seq. of clusters in the path : ");
		for (k=0; k<tracePaths[j].size (); k++){
			cout<<"  "<<tracePaths[j][k]+1;
			trace = points2pdb(newClusters[tracePaths[j][k]], "GLY", " CA ");
			pathTrace.append(trace, trace.numOfAA()-1, trace.numOfAA()-1);
			pathTrace.header[0] += toString(tracePaths[j][k]+1);
			pathTrace.header[0] += " ";
		}
		pathTrace.header[0] += "  eLength= ";
		pathTrace.header[0] += toString(traceLengths[j]);

		
		cout<<" length= "<<traceLengths[j]<<endl;
		pathTrace.writePDB(outPath+"_pathTrace_"+toString(tracePaths[j][0]-nClusters)+"w"+toString(tracePaths[j][tracePaths[j].size ()-1]-nClusters)+"_"+toString(j+1)+".pdb", 1, pathTrace.numOfAA(), 1);
	}

	//set weights on the graph
	short bestTraceIndx1,bestTraceIndx2;
	float loopLength = 0;			//the length of the loop b/w the two SS
	for (i=1; i<graph.size ()-1; i++){
		for (j=0; j<graph[i].size (); j++){
			for (k=0; k<graph[i][j].outLinks .size (); k++){
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
						graph[i][j].outLinks [k].w = getBestPathTrace(newClusters, tracePaths, traceLengths, nClusters, loopLength, j+nClusters + 1, graph[i][j].outLinks[k].colIndx + nClusters, bestTraceIndx1, bestTraceIndx2);
					else
						graph[i][j].outLinks [k].w = getBestPathTrace(newClusters, tracePaths, traceLengths, nClusters, loopLength, j+nClusters - 1, graph[i][j].outLinks[k].colIndx + nClusters, bestTraceIndx1, bestTraceIndx2);


					graph[i][j].traceNum [k] = toString(bestTraceIndx1+1);
					graph[i][j].traceNum [k] += "-";
					graph[i][j].traceNum [k] += toString(bestTraceIndx2+1);

					//update inLink in the other node
					for (l=0; l<graph[graph[i][j].outLinks [k].rowIndx][graph[i][j].outLinks [k].colIndx].inLinks.size (); l++){
						if (graph[graph[i][j].outLinks [k].rowIndx][graph[i][j].outLinks [k].colIndx].inLinks[l].rowIndx == i){
							if (graph[graph[i][j].outLinks [k].rowIndx][graph[i][j].outLinks [k].colIndx].inLinks[l].colIndx == j){
								graph[graph[i][j].outLinks [k].rowIndx][graph[i][j].outLinks [k].colIndx].inLinks[l].w = graph[i][j].outLinks [k].w;
								break;
							}
						}
					}
				}
			}
			//getchar();getchar();
			
		}
	}


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

	pair<short,float> **cDist;						//shortest distances b/w clusters (indeces and distance) ... each cluster will represent a node in a network later
	float continuatyTHR = 4.1 * inMRC.apixX;		//a threshold used as a cuttoff distance b/w any two points in any two clusters to consider them continued (non-disjoint)
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
	peakClustering(pnts, clustersCoord, ssEdges[0][0], 0.1*inMRC.apixX); 

	nClusters = clustersCoord.size ();
	
	
	//for SS ends....create a cluster of length = 1 contains one end point for each SS
	vector<Coordinate> tmpCluster;
	for (i=0; i<ssEdges.size (); i++){

		tmpCluster.clear ();
		tmpCluster.push_back (ssEdges[i][0]);
		clustersCoord.push_back (tmpCluster);


		tmpCluster.clear ();
		tmpCluster.push_back (ssEdges[i][ssEdges[i].size ()-1]);
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
	buildAdjMtrx(clustersCoord, cDist, adjMtrx, inMRC.apixX, continuatyTHR);


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
	short  *tmpShortestPathRow;
	short  *tmpShortestPathCol;

	//tmpShortestPath = new pair<short, short> [nSticks+2];			//+2 for start and end nodes
	tmpShortestPathRow = new short [nSticks+2];						//+2 for start and end nodes
	tmpShortestPathCol = new short [nSticks+2];						//+2 for start and end nodes

	//pair<short, short> tmpNode;

	if (rowIndx != -1){
		//add the first node in the path

		tmpShortestPathRow[nRslvd] = rowIndx;
		tmpShortestPathCol[nRslvd++] = colIndx;
		

		short colIndxDivide2;
		short colIndxMod2;
		short sColMod2;
		short sColDivide2;


		for (icol=0; icol<nSticks-nVisitedCols-1; icol++){
			colIndxDivide2 = colIndx/2;
			colIndxMod2 = colIndx%2;
			sColMod2 = tmpShortestPathCol[nRslvd-1]%2;
			sColDivide2 = tmpShortestPathCol[nRslvd-1]/2;

			rowIndx = graphSets[rowIndx][colIndxDivide2].nxtRow[colIndxMod2][setIndx];
			colIndx = graphSets[tmpShortestPathRow[nRslvd-1]][colIndxDivide2].nxtCol[colIndxMod2][setIndx];
			setIndx = graphSets[tmpShortestPathRow[nRslvd-1]][sColDivide2].nxtSetIndx[sColMod2][setIndx];

			tmpShortestPathRow[nRslvd] = rowIndx;
			tmpShortestPathCol[nRslvd++] = colIndx;
		}
	}


	iPath tmpPath;

	tmpPath.pCandidateRow = tmpShortestPathRow;
	tmpPath.pCandidateCol = tmpShortestPathCol;

	tmpPath.w = shortestW;

	return tmpPath;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//find K-shortest Paths
void findKShortestPaths(cellSets **graphSets, vector<vector<cell> > &graph, unsigned int **cntr, short nRows, short nSticks, short sRow, short sCol, int K){

	
#ifdef _WIN				//if work under Windows
	clock_t start, finish;
	start = clock();
#else
	struct timeval ti_start, ti_end;
	double time_di;
	cout<<"Time(0)= "<<time(0)<<endl;
	//start timer
	gettimeofday(&ti_start,0);
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
	bestPath.pCandidateRow[0] = 0;
	bestPath.pCandidateCol [0] = 0;

	//add end node
	bestPath.pCandidateRow [nSticks+1] = graph.size ()-1;
	bestPath.pCandidateCol [nSticks+1] = 0;

	cout<<" shortest path "<<endl;
	for (int kk =0; kk<nSticks+2; kk++){
		cout<<"["<<bestPath.pCandidateRow[kk]<<" , "<<bestPath.pCandidateCol[kk]<<"] ";
	}
	cout<<" Wt = "<<bestPath.w<<endl;


	//find the top K-1 shortest Paths		
	//initiate T and X sets
	//T is the list of best k solved shortest paths, X is the list of candidates of shortest paths
	iPath *T, *X;
	int Tcntr=0, Xcntr=0, i;
	T = new iPath[K];
	X = new iPath [K*nSticks];			//the maximum number of candidates in X would be K*nSticks becuase for each path from T the max possible paths could be generated is the number of edges in the path

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
		skipLinks_Row [skipLinkSize] = bestPath.pCandidateRow [bestPath.deviationNode +1];
		skipLinks_Col [skipLinkSize++] = bestPath.pCandidateCol [bestPath.deviationNode +1];
		
		int parent = bestPath.parentPath;
		while (parent != -1){
			if (T[parent].pCandidateRow [bestPath.deviationNode ] == bestPath.pCandidateRow [bestPath.deviationNode ] &&
				T[parent].pCandidateCol [bestPath.deviationNode ]  == bestPath.pCandidateCol [bestPath.deviationNode ]){

				skipLinks_Row [skipLinkSize] = T[parent].pCandidateRow [bestPath.deviationNode+1];
				skipLinks_Col [skipLinkSize++] = T[parent].pCandidateCol [bestPath.deviationNode +1];

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
				preSet = preSet | (1 << bestPath.pCandidateCol [i]/2);
			}

			//get the next best candidate deviates from bestPath
			tmpPath = getShortestPath(	graphSets, 
										graph, 
										cntr, 
										bestPath.pCandidateRow [icol], 
										bestPath.pCandidateCol[icol], 
										nResolved,
										skipLinks_Row,
										skipLinks_Col,
										skipLinkSize,
										preSet);
			tmpPath.deviationNode = icol;
			tmpPath.parentPath = Tcntr-1;

			//insert first portion of the path
			//insert deviation node
			tmpPath.pCandidateRow [icol] = bestPath.pCandidateRow [icol];
			tmpPath.pCandidateCol [icol] = bestPath.pCandidateCol [icol];
			for (i=icol-1; i>=0; i--){
				tmpPath.pCandidateRow [i] = bestPath.pCandidateRow [i];
				tmpPath.pCandidateCol [i] = bestPath.pCandidateCol [i];

				//add weights
				for (int fromlink=0; fromlink<graph[tmpPath.pCandidateRow [i]][tmpPath.pCandidateCol [i]].outLinks.size(); fromlink++){
					if (graph[tmpPath.pCandidateRow [i]][tmpPath.pCandidateCol [i]].outLinks [fromlink].rowIndx == tmpPath.pCandidateRow [i+1] &&
						graph[tmpPath.pCandidateRow [i]][tmpPath.pCandidateCol [i]].outLinks [fromlink].colIndx == tmpPath.pCandidateCol [i+1]){
						tmpPath.w += graph[tmpPath.pCandidateRow [i]][tmpPath.pCandidateCol [i]].outLinks [fromlink].w;
						break;
					}
				}
			}

			//insert end node
			tmpPath.pCandidateRow [nSticks+1] = bestPath.pCandidateRow[nSticks+1];
			tmpPath.pCandidateCol [nSticks+1] = bestPath.pCandidateCol[nSticks+1];

			//print the path
			//for (i=0; i<nSticks+2; i++)
			//	cout<<"  ["<<tmpPath.pCandidateRow[i]<<" "<<tmpPath.pCandidateCol [i]<<"] ";
			//cout<<"  weight= "<<tmpPath.w <<" dv= "<<tmpPath.deviationNode <<" prnt= "<<tmpPath.parentPath <<endl;

			//save the candidate to X list
			X[Xcntr++] = tmpPath;
			//reset preSet
			preSet = 0;
			icol++;
			//reset skiplinks for the new deviation node
			skipLinkSize = 0;
			skipLinks_Row [skipLinkSize] = bestPath.pCandidateRow [icol+1];
			skipLinks_Col [skipLinkSize++] = bestPath.pCandidateCol [icol+1];

		}while(icol<nSticks);

		k+= 1;
		//getchar();
	}

#ifdef _WIN
	finish = clock();
	cout<<"Time taken to enumerate top "<<K<<" paths is "<<finish-start<<" ms."<<endl;
	start = clock();
#else
	//stop timer
	gettimeofday(&ti_end,0);
	time_di = (ti_end.tv_sec-ti_start.tv_sec)*1000000 + ti_end.tv_usec - ti_start.tv_usec;
	cout<<"Time(0)= "<<time(0)<<endl;
	cout<<"Time taken to enumerate top "<<K<<" paths is "<<(double) (time_di/1000)<<" ms "<<endl;
#endif

	//print top k paths
	for (int iT=0 ;iT<Tcntr; iT++){
		cout<<"path# "<<iT+1<<" : ";
		for (icol=1; icol<=nSticks; icol++)
			cout<<"["<<T[iT].pCandidateRow[icol]<<","<<T[iT].pCandidateCol [icol] <<"] ";
		cout<<"    Weight: "<<T[iT].w<<endl;
	}

	cout<<"Done.."<<endl;
	//dispose memory
	delete []X;
	delete []T;
//	delete []skipLinks_Row;
//	delete []skipLinks_Col;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//finds the shortest path in the graph using a recurrence....saves all combinations from the cell down to the end node
//given is the index of start node (row and col) and the index of end node (row and col)
void gShortestRec(vector<vector<cell> > graph, short sRow, short sCol, short eRow, short eCol, int K){


	short irow, icol;

	short nSticks = graph[0].size ()/2;
	short nRows = graph.size();				

	//cntr saves the last indx in sets for each node
	unsigned int **cntr = AllocateDynamicArray<unsigned int> (nRows, nSticks);

	//combination index used to hash the position of a combination
	unsigned int *combinationIndx ;				
	for (irow=0; irow<nRows;irow++)
		for (icol=0; icol<nSticks; icol++)
			cntr[irow][icol] = 0;

	//initialize data structures
	cellSets **graphSets = AllocateDynamicArray<cellSets> (nRows, nSticks); //number of cols is half number of cols in the original graph

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

	short fromrow, fromcol, fromlink, torow, tocol, gcol;			//counters
	unsigned int itoset, totalNumSets=0;							//total number of sets for the protein

	//estimate number of sets for each node
	vector<vector<vector<short> > > combSets (graph.size (), vector<vector< short> > (nSticks));
	combSets[graph.size ()-1][0].push_back (0);
	for (irow=graph.size ()-2; irow>0; irow--){
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
	vector<vector<vector<long> > > nSets (graph.size (), vector<vector< long> > (nSticks));
	nSets[0].resize (nSticks, vector<long> (1, 0));
	nSets[graph.size ()-1].resize (nSticks, vector<long> (1, 0));
	for (irow=1; irow<graph.size ()-1; irow++){
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
	cout<<"Press any key to start..."<<endl;getchar();getchar();
	cout<<"Initializing..."<<endl;
	unsigned int i;
	int k;

#ifdef _WIN				//if work under Windows
	clock_t start, finish;
	start = clock();
#else
	struct timeval ti_start, ti_end;
	double time_di;
	cout<<"Time(0)= "<<time(0)<<endl;
	//start timer
	gettimeofday(&ti_start,0);
#endif

	//initialize indeces table
	unsigned int indecesTableSize = pow(2, nSticks);
	combinationIndx = new unsigned int [indecesTableSize];

	for (i=0; i<indecesTableSize; i++)
		combinationIndx[i] = -1;


	//initialize some data structure
	for (irow=sRow+1; irow<eRow; irow++){		
		for (icol=0; icol<nSticks; icol++){
			//cout<<"     nSets= "<<nSets[irow][icol][0]<<endl;			
			
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

	nSets.clear ();

	
	cout<<"Start..."<<endl;	
	
	//the recurrence
	for (fromrow=eRow-1; fromrow>sRow; fromrow--){
		
		for (gcol=0; gcol<graph[fromrow].size(); gcol++){		//graph col
			short gcolmod2 = gcol%2;
			fromcol = gcol/2;
			
			//test for all out links
			for (fromlink=0; fromlink<graph[fromrow][gcol].outLinks .size (); fromlink++){

				tocol = graph[fromrow][gcol].outLinks [fromlink].colIndx/2;		//the index of the col of the node on the tail of the link
				torow = graph[fromrow][gcol].outLinks [fromlink].rowIndx;		//the indx  of the row of the node on the tail of the link
				short gtocol = graph[fromrow][gcol].outLinks[fromlink].colIndx%2;
				
				unsigned int combIndx=-1, setIndx;

				//for all sets on the connected node
				for (itoset=0; itoset<cntr[torow][tocol];itoset++){
									
					if ((((1 << fromcol) & graphSets[torow][tocol].sets [itoset]) == 0) &&				//exclude paths already have this col
						graphSets[torow][tocol].length[itoset] + fromrow - sRow >= nSticks){			//exclude unfeasable paths...paths that too short

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

				//if the node has no sets yet....then it is connected to END node or the eRow
				if (cntr[torow][tocol]<1){

					short combSize = 1;

					//get the combination ... is the set contains the current col
					int tmpSet = 1 << fromcol;

					if (torow < eRow-1){
						tmpSet = tmpSet | (1<<tocol);
						combSize = 2;
					}


					if (combinationIndx[tmpSet] != -1){
						//I have the combination in the list of sets
						setIndx = combinationIndx [tmpSet];
					}
					else{
						//we don't have it in the list...add it and save the index
						setIndx = cntr[fromrow][fromcol]++;
						combinationIndx[tmpSet] = setIndx;
						graphSets[fromrow][fromcol].sets[setIndx] = tmpSet;

					}
						
					// update minimum cost and parent information (row, col, set indx)
					//if (combSize>1){										
						graphSets[fromrow][fromcol].nxtCol[gcolmod2][setIndx] = graph[fromrow][gcol].outLinks[fromlink].colIndx;
						graphSets[fromrow][fromcol].nxtRow[gcolmod2][setIndx] = torow;
						graphSets[fromrow][fromcol].nxtSetIndx[gcolmod2][setIndx] = -1;
					//}
					
					graphSets[fromrow][fromcol].weight[gcolmod2][setIndx] = graph[fromrow][gcol].outLinks[fromlink].w;
					graphSets[fromrow][fromcol].length[setIndx] = combSize;
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
	cout<<"Time taken to build all Sets is "<<finish-start<<" ms."<<endl;
#else
	//stop timer
	gettimeofday(&ti_end,0);
	time_di = (ti_end.tv_sec-ti_start.tv_sec)*1000000 + ti_end.tv_usec - ti_start.tv_usec;
	cout<<"Time(0)= "<<time(0)<<endl;
	cout<<"Time taken to enumerate top "<<K<<" paths is "<<(double) (time_di/1000)<<" ms "<<endl;
#endif

	cout<<"End... "<<endl;

/*
	//print information
	cout<<"++++++++++++++++++++++ print Information +++++++++++++"<<endl;
	int j;
	for (irow=0; irow<nRows; irow++){
		for (icol=0; icol<nSticks;icol++){
			cout<<"[ "<<irow<<" , "<<icol<<" ]  : "<<endl;
			for (i=0; i<cntr[irow][icol]; i++){
				cout<<"    set# "<<i+1<<" length= "<<graphSets[irow][icol].length [i]<<" : ";
				for (j=0; j<graphSets[irow][icol].length[i]; j++){
					if (graphSets[irow][icol].sets[i].bitSet [j])
						cout<<j<<" ";
				}
				cout<<endl;
				cout<<"     shortest path F (w= "<<graphSets[irow][icol].weight[0][i]<<") : ";
			//	for (j=0; j<graphSets[irow][icol].fPath[i].size(); j++){
			//		if (graphSets[irow][icol].fPath[i][j] != -1)
			//			cout<<graphSets[irow][icol].fPath[i][j]<<" ";
			//	}
				cout<<endl;
				cout<<"     shortest path B (w= "<<graphSets[irow][icol].weight[1][i]<<") : ";
			//	for (j=0; j<graphSets[irow][icol].bPath[i].size(); j++){
			//		if (graphSets[irow][icol].bPath[i][j] != -1)
			//			cout<<graphSets[irow][icol].bPath[i][j]<<" ";
			//	}
				cout<<endl;
			}
		}
	}
*/
	findKShortestPaths(graphSets, graph, cntr, nRows, nSticks, sRow, sCol, K);

	//dispose memory
	FreeDynamicArray<cellSets> (graphSets);
	FreeDynamicArray<unsigned int> (cntr);
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
#endif