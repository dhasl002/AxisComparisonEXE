//fast-RNN v2.0
//    Copyright (C) 2012  Roberto J. López-Sastre (robertoj.lopez@uah.es)
//                        Daniel Oñoro-Rubio
//			  Víctor Carrasco-Valdelvira
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.


#ifndef AGGCLUSTERING_H
#define AGGCLUSTERING_H

//Includes

#include <fstream>
#include <vector>
#include <list>
#include <sys/time.h>
#include <stdlib.h>
#include <limits>
#include <string.h>
#include <math.h>
#include "Point.h"

using namespace std;



//MACRO
#define TIME_THIS(X,Y) \
{ \
  struct timeval t_ini, t_fin; \
  gettimeofday(&t_ini, NULL); \
  X; \
  gettimeofday(&t_fin, NULL); \
  Y=timeval_diff(&t_fin, &t_ini); \
}


typedef struct scluster{
  vector<double> centroid; //centroid
  vector<unsigned int> data_index; //index of vectors of this cluster
  double cvar; //cluster varianceabout:startpage
}cluster;


typedef list <vector <double> > matrix_data;


//element structure
typedef struct selement{
  list <cluster>::iterator it;
  bool mask;
}element;

typedef vector <element> fmap;

//Struct of candidates to nn
typedef struct scandidate{
  list <cluster>::iterator it;
  unsigned int index;
}candidate;

typedef list <candidate> slice;

//functions

vector<double> centroid_mul(double n, vector<double> centroid);
vector<double> centroid_plus(vector<double> A, vector<double> B);
vector<double> centroid_div(double n, vector<double> centroid);
vector<double> centroid_diff(vector<double> A, vector<double> B);
double squared_magnitude(vector<double>);


  char * data_file_path = NULL;
  char * out_file_path = NULL;
  unsigned int dim;
  unsigned int n;
  double thres;    //threshold for agglomerative clustering
  double thres_euclid;
  double ep; //threshold for slicing
  unsigned int nc; //num clusters
  int free_top;
  unsigned int COMP;
  unsigned long ndist;


vector<double> centroid_mul(double n, vector<double> centroid){
	int size = centroid.size();
	for (int i=0; i<size; i++)
		centroid[i] = centroid[i] * n;
	return centroid;
}

vector<double> centroid_plus(vector<double> A, vector<double> B){
	vector<double> centroid;
	int size = A.size();
	for (int i=0; i<size; i++)
		centroid.push_back(A[i] + B[i]);
	return centroid;
}

vector<double> centroid_div(double n, vector<double> centroid){
	int size = centroid.size();
	for (int i=0; i<size; i++)
		centroid[i] = centroid[i] / n;
	return centroid;
}

vector<double> centroid_diff(vector<double> A, vector<double> B){
	vector<double> centroid;
	int size = A.size();
	for (int i=0; i<size; i++)
		centroid.push_back(A[i] - B[i]);
	return centroid;
}

double squared_magnitude(vector<double> vec){
	int size = vec.size();
	double sum = 0;
	for (int i=0; i<size; i++)
		sum += vec[i]*vec[i];
	return sum;
}


double timeval_diff(struct timeval *a, struct timeval *b)
{
  return
    (double)(a->tv_sec + (double)a->tv_usec/1000000) -
    (double)(b->tv_sec + (double)b->tv_usec/1000000);
}


void agglomerate_clusters(cluster &C1, cluster &C2)
{
  //Agglomerate two clusters
  //C1=C1+C2

  unsigned int m = C2.data_index.size();
  unsigned int n = C1.data_index.size();
  double d_sum = double(n+m);

  //Copy index values
  for(unsigned int i=0; i<m; i++)
    C1.data_index.push_back(C2.data_index[i]);

  //update centroid
  C1.centroid=centroid_div(d_sum,centroid_plus(centroid_mul(double(n),C1.centroid),centroid_mul(double(m),C2.centroid)));

  //update variance
  vector <double> diff=centroid_diff(C1.centroid,C2.centroid);
  C1.cvar=((double(n)*C1.cvar)+(double(m)*C2.cvar)+(((double(n*m))/(d_sum))*(squared_magnitude(diff))))/(d_sum);

}

//Get Nearest Neighbor
void get_nn(cluster &C, list <cluster> &R, double &sim, list <cluster>::iterator &iNN)
{
  //Return the NN of cluster C in the list R.
  //iNN is the iterator of the NN in R, and sim is the similarity.

  unsigned int n = R.size();
  list <cluster>::iterator it, itBegin = R.begin(), itEnd = R.end();
  double d;
  vector <double> diff;

  if(n > 0)
  {
    //First iteration
    it = itBegin;
    diff = centroid_diff(C.centroid,(*it).centroid);
    sim=-(C.cvar + (*it).cvar + squared_magnitude(diff));
    iNN = it;
    it++;

    for(; it!=itEnd; it++)
    {
      diff = centroid_diff(C.centroid,(*it).centroid);

      d=-(C.cvar + (*it).cvar + squared_magnitude(diff));

      if(d>sim)
      {
	sim = d;
	iNN = it;
      }
    }
  }
  else
  {
    cout << "Warning: R is empty (function: get_nn)" << endl;
    sim = 0;
  }
}


void get_nn_in_slices(cluster &C, list <cluster> &X, slice &Si, slice &So, double &sim, list <cluster>::iterator &iNN, double limit)
{
  //Search within the interior and exterior slices, Si and So respectively
  double d;
  slice::iterator itS, endS;
  list <cluster>::iterator it;
  vector <double> diff;
  bool isfirst = true;

  //Search first in the interior slice
  if(Si.size() > 0)
  {
    endS=Si.end();
    //First iteration
    isfirst = false;
    itS= Si.begin();
    it = (*itS).it;
    diff = centroid_diff(C.centroid,(*it).centroid);
    sim=-(C.cvar + (*it).cvar + squared_magnitude(diff));
    iNN = it;

    itS++;
    for(; itS!=endS; itS++)
    {
      it = (*itS).it;
      diff = centroid_diff(C.centroid,(*it).centroid);
      d=-(C.cvar + (*it).cvar + squared_magnitude(diff));
      if(d>sim)
      {
	sim = d;
	iNN = it;
      }
    }

    //DEBUG
    ndist += Si.size();

    //Do we need to search in the exterior slice?
    if(sim >= limit)
    {
      return; //NO
    }
  }

  //Search in the exterior slice (if any)
  if(So.size()>0)
  {
    endS=So.end();
    if(isfirst)
    {
      //First iteration
      isfirst = false;
      itS= So.begin();
      it = (*itS).it;
      diff = centroid_diff(C.centroid,(*it).centroid);
      sim=-(C.cvar + (*it).cvar + squared_magnitude(diff));
      iNN = it;
    }

    for(itS=So.begin(); itS!=endS; itS++)
    {
      it = (*itS).it;
      diff = centroid_diff(C.centroid,(*it).centroid);
      d=-(C.cvar + (*it).cvar + squared_magnitude(diff));
      if(d>sim)
      {
	sim = d;
	iNN = it;
      }
    }

    //DEBUG
    ndist += So.size();
  }
}

int unsigned bsearchL(fmap &f_map, double d)
{
  //Search in f_map to the left
  int b = 0, c, t = f_map.size() - 1, aux_c;
  double q;
  bool end = false;

  // Move until we find a top value not erased
  while( !f_map[t].mask )
    t--;

  // Move until we find a base value not erased
  while( !f_map[b].mask )
    b++;

  while( ((t - b) > 1) && !end )
  {
    //the middle
    c = (b + t) >> 1;

    //erased
    if( !f_map[c].mask )
    {
       // Search a not erased element

       // Fitst iteration
       aux_c = c + 1;

       //Searching upward
       while( !f_map[aux_c].mask && (aux_c < t) )
	 aux_c++;


       //Do we need to search downward?
       if( !f_map[aux_c].mask || (aux_c >= t) )
       {
          aux_c = c - 1;

          // Searching downward
          while( !f_map[aux_c].mask && (aux_c > b) )
	    aux_c--;

       }


       if( aux_c == b )
	 end = true;
       else
	 c = aux_c;

    }//if erased

    if( !end && !f_map[c].mask )
    {
       cout << "bsearchL failed" << endl;
       exit(-1);
    }

    if(!end)
    {
      q = (*f_map[c].it).centroid[COMP];

      if(d < q)
	t = c;
      else if(d > q)
	b = c;
      else
	return c;

    }
  }

  if( !f_map[b].mask || !f_map[t].mask )
  {
    cout << "Error: bsearchL failed" << endl;
    exit(-1);
  }

  return ((d <= (*f_map[b].it).centroid[COMP]) ? b : t);

}



unsigned int bsearchR(fmap &f_map, double d)
{
  int b = 0, c, t = f_map.size() - 1, aux_c;
  double q;
  bool end = false;

  // Move until find a top not erased
  while( !f_map[t].mask )
    t--;

  // Move until find a base not erased
  while( !f_map[b].mask )
    b++;

  while( (t - b > 1) && !end )
  {

    c = (b + t) >> 1;

    //erased
    if( !f_map[c].mask )
    {
      aux_c = c + 1;

      // Search upward
      while( !f_map[aux_c].mask && (aux_c < t) )
	aux_c++;


      //Do we need to search downward
      if( !f_map[aux_c].mask || (aux_c >= t) )
      {
	aux_c = c - 1;

	// Search downward
	while( !f_map[aux_c].mask && (aux_c > b) )
	  aux_c--;
      }

       if( aux_c == b )
	 end = true;
       else
	 c = aux_c;
    }


    if( !end && !f_map[c].mask )
    {
       cout << "error: bsearchR failed" << endl;
       exit(-1);
    }



    if(!end)
    {
      q = (*f_map[c].it).centroid[COMP];
      if(d < q)
	t = c;
      else if(d > q)
	b = c;
      else
	return c;
    }
  }


  if( !f_map[b].mask || !f_map[t].mask )
  {
    cout << "error: bsearchR failed" << endl;
    exit(-1);
  }


  return ((d >= (*f_map[t].it).centroid[COMP]) ? t : b);
}

int bsearch(fmap &f_map, double d, unsigned int &b,unsigned int &t)
{
  unsigned int leng = f_map.size() - 1;
  b=0; t=leng;
  unsigned int c, aux_c;
  double q;
  bool end = false;

  //highest no deleted position
  while( (t > 0) && !f_map[t].mask )
    t--;
  //lowest no deleted position
  while( (b < leng) && !f_map[b].mask )
    b++;

  //Check conditions
  if( b > t )
    return -1;

  //Binary search
  while( ((t - b) > 1) && !end )
  {
    c = (b + t) >> 1;

    //is it erased?
    if( !f_map[c].mask )
    {
      aux_c = c + 1;

      // Search upward
      while( !f_map[aux_c].mask && (aux_c < t) )
      	aux_c++;


      // Do we have to search downward?
      if( !f_map[aux_c].mask || (aux_c >= t) )
      {
	aux_c = c-1;

	// Search downward
	while( !f_map[aux_c].mask && (aux_c > b) )
	  aux_c--;
      }

      if( aux_c == b )
      	end = true;
      else
	c = aux_c;
    }

    if( !end && !f_map[c].mask )
    {
       cout << "error: bsearch error" << endl;
       exit(-1);
    }


    if(!end)
    {
      q = (*f_map[c].it).centroid[COMP];
      if(d < q)
	t = c;
      else if(d > q)
	b = c;
      else
	return c;
    }
  }

  if( !f_map[b].mask || !f_map[t].mask )
  {
    cout << "error: bsearch error" << endl;
    exit(-1);
  }

  if(b == t)
    return ((*f_map[b].it).centroid[COMP] >= d) ? t : t+1;
  else
    return ((*f_map[b].it).centroid[COMP] >= d) ? b : t;

}


void init_slices(fmap &f_map, list <cluster> &X, slice &Si, slice &So, double V, double e)
{
  //Generate the slice in the space where the NN candidates are. The slice has 2e width.
  unsigned int min,max,bmax,bmin,tmax,tmin,i;
  candidate c;

  //Three slices? (recall: ep is the parameter for slicing)
  if(e > ep)
  {
    //Build interior slice
    min = bsearchL(f_map,V-ep);
    max = bsearchR(f_map,V+ep);

    for(i=min; i<=max; i++)
    {
      if( f_map[i].mask )
      {
        c.it = f_map[i].it;
        c.index = i;
        Si.push_back( c );
      }
    }

    if(min != 0) //generate bottom candidate list
    {
      bmax = min-1;
      //Build bottom slice
      bmin = bsearchL(f_map,V-e);

      for(i=bmin; i<=bmax; i++)
      {
	if( f_map[i].mask )
	{
	  c.it = f_map[i].it;
	  c.index = i;
	  So.push_back( c );
	}
      }
    }

    if(max != (f_map.size()-1))
    {
      tmin = max+1;

      //Build top slice
      tmax = bsearchR(f_map,V+e);

      for(i=tmin; i<=tmax; i++)
      {
	if( f_map[i].mask )
	{
	  c.it = f_map[i].it;
	  c.index = i;
	  So.push_back( c );
	}
      }
    }
  }
  else //only one slice
  {
    min = bsearchL(f_map,V-e);
    max = bsearchR(f_map,V+e);

    for(i=min; i<=max; i++)
    {
      if( f_map[i].mask )
      {
        c.it = f_map[i].it;
        c.index = i;
        Si.push_back( c );
      }
    }
  }
}

inline void erase_element(fmap &f_map, list <cluster> &X, list <cluster>::iterator it)
{
  int l = f_map.size(), i;

  for(i=0; i<l; i++)
  {
    if(f_map[i].mask)
      if(f_map[i].it == it)
	{
	  f_map[i].mask = false;
	  X.erase(it);

	  //Update free_top
	  free_top = (i>free_top) ? i : free_top;
	  return;
	}
  }

  cout << "error: erasing element" << endl;
  exit(-1);
}

void insert_element(fmap &f_map, list <cluster> &X, cluster &V)
{
  //Insert element V in X and update f_map
  element elem2insert, aux_elem;
  int pos=0;
  unsigned int b,t;
  bool update_free_top;

  //Push back the element in X
  X.push_back(V);

  // Initialize the element to insert
  elem2insert.mask = true;

  elem2insert.it = X.end(); elem2insert.it--;

  //Search for a position
  pos = bsearch(f_map,V.centroid[COMP],b,t);

  // f_map is empty
  if( -1 == pos )
  {
    pos = t;
    f_map[pos] = elem2insert;
    return;
  }

  //Is pos the last position?
  if( f_map.size() == (unsigned int) pos )
  {
    pos--;
    if(!f_map[pos].mask)
    {
      f_map[pos] = elem2insert;
      return;
    }

    //Insert downwards

    while(elem2insert.mask)
    {

      //Save current pos element
      aux_elem = f_map[pos];

      //Insert element in pos
      f_map[pos] = elem2insert;

      //Update elem2insert
      elem2insert = aux_elem;

      pos--;

    }
    free_top=pos;
    return;


  }

  if(pos >= free_top)
  {
    //Insert downwards
    while(elem2insert.mask)
    {
      pos--;
      //Save current pos element
      aux_elem = f_map[pos];

      //Insert element in pos
      f_map[pos] = elem2insert;

      //Update elem2insert
      elem2insert = aux_elem;
    }
    //update free_top?
    update_free_top = true;


  }
  else //upwards
  {
    if(f_map[pos].mask)
      if(V.centroid[COMP]>=(*f_map[pos].it).centroid[COMP])
	pos++;

    while(elem2insert.mask)
    {
      //Save current pos element
      aux_elem = f_map[pos];

      //Insert element in pos
      f_map[pos] = elem2insert;

      //Update elem2insert
      elem2insert = aux_elem;

      //Update index
      pos++;
    }

    //update free_top?
    update_free_top = (pos<free_top) ? false : true;


  }

  //Update free_top just in case it has been occupied
  if(update_free_top)
  {
    free_top = pos-1;
    while((free_top>0) &&(f_map[free_top].mask))
      free_top--;
  }

}

void init_map(fmap &f_map, list <cluster> &X)
{

  //Create forward map
  list <cluster>::iterator itX, endX;
  list <Point<list <cluster>::iterator> > aux_list;
  list <Point<list <cluster>::iterator> >::iterator it,itend;
  element eaux;

  //Get the list of iterators
  endX = X.end();
  for(itX=X.begin(); itX!=endX; itX++)
    aux_list.push_back(Point<list <cluster>::iterator>(itX, (*itX).centroid[COMP]));

  //Sorting
  aux_list.sort();

  eaux.mask = true; //Mask = true for all elements

  //Convert the sorted list to f_map
  itend = aux_list.end();
  for(it=aux_list.begin(); it!=itend; it++)
  {
    eaux.it = (*it).get_index(); //save the iterators
    f_map.push_back( eaux );
  }

  //Update free top
  free_top = -1;
}

unsigned int agg_clustering_fast_rnn(matrix_data &X, vector <unsigned int> &labels, double agg_thres,vector< vector<double> > &cluster_centre )
{
  //This function computes the Fast RNN (reciprocal nearest neighbors) clustering.

  //Chain for the clusters R
  list <cluster> R(n);
  list <cluster> ::iterator it, itEnd = R.end(), iNN, penult;
  matrix_data ::iterator Xit = X.begin();
  unsigned int Xindex = 0;
  double sim = 0;
  double l_agg_thres = (-1*agg_thres);
  double epsilon;
  slice Si,So; //slices with candidates
  bool RNNfound = false;


  //DEBUG
  ndist = 0;

  //Initialize list R - each vector in a separate cluster (start point of the algorithm)
  for(it=R.begin(); it!=itEnd; it++, Xit++, Xindex++)
  {
    //update index of the vector
    (*it).data_index.push_back(Xindex);
    //update centroid
    (*it).centroid = *Xit;
    //update variance
    (*it).cvar = 0;
  }

  //NN-Chain (the pair at the end of this chain is always a RNN)
  list <cluster> L;
  //Chain for similarities
  list <double>  Lsim;
  //chain for the (final) clusters C
  list <cluster> C;

  //Create forward map
  fmap f_map;
  init_map(f_map,R);


  //The algorithm starts with a random cluster
  srand( time(NULL) );
  unsigned int rp = (unsigned int)rand() % n; //random integer [0,n-1]



  //Add to L
  it = R.begin();
  advance(it,rp);
  L.push_back(*it);




  //R\rp -> delete the cluster in R and mark as erased in f_map
  erase_element(f_map, R, it);


  //First iteration
  if( R.size() > 0 )
  {
    //Get nearest neighbor
    get_nn(L.back(),R,sim,iNN);

    //DEBUG
    ndist += R.size();



    //Add to the NN chain
    L.push_back(*iNN); //add to L
    erase_element(f_map, R, iNN); //delete from R
    Lsim.push_back(sim);//add to Lsim


    //Only two clusters?
    if(R.size() == 0)
    {
      penult=L.end(); penult--; penult--;
      //check the similarity (last element)
      if(sim > l_agg_thres)
      {
	//Agglomerate clusters
	agglomerate_clusters(L.back(),*penult);
	C.push_back(L.back());
      }
      else
      {
	//Save in C separately
	C.push_back(*penult);
	C.push_back(L.back());
      }
      L.clear(); //free memory
    }
  }
  else //R is empty
  {
    if(L.size() == 1)
      C.push_back(L.back()); // Only one vector, only one cluster

    L.clear();  //free memory
  }
  //Main loop
  while(R.size() > 0)
  {
    RNNfound = false;

    //Clear slices
    Si.clear();
    So.clear();

    //Update epsilon with the last sim
    epsilon = sqrt(-1*Lsim.back());

    epsilon = (epsilon < ep) ? epsilon : ep;

    //Identify slices
    init_slices(f_map,R,Si,So,L.back().centroid[COMP],epsilon);


    if((Si.size()>0) || (So.size()>0)) //Search for a NN within the candidate list
    {

      get_nn_in_slices(L.back(),R,Si,So,sim,iNN,l_agg_thres);



      if(sim > Lsim.back()) //no RNN
      {

	//No RNNs, add s to the NN chain
	L.push_back(*iNN); //add to L
	erase_element(f_map, R, iNN); //delete from R
	Lsim.push_back(sim);//add to Lsim



	if(R.size()==0) //R has been emptied
	{
	  //check the last similarity
	  if(Lsim.back() > l_agg_thres)
	  {
	    //Agglomerate clusters
	    penult=L.end(); penult--; penult--;
	    agglomerate_clusters(L.back(),*penult);
	    insert_element(f_map,R,L.back());



	    //delete the last two elements in L
	    L.pop_back();
	    L.pop_back();

	    //delete similarities
	    Lsim.pop_back();
	    if(Lsim.size()>=1)
	      Lsim.pop_back();

	    //Initialize the chain with the following nearest neighbour
	    if(L.size() == 1)
	    {
	      //Get nearest neighbor
	      get_nn(L.back(),R,sim,iNN);


	      ndist += R.size();

	      //Add to the NN chain
	      L.push_back(*iNN); //add to L
	      erase_element(f_map, R, iNN); //delete from R
	      Lsim.push_back(sim);//add to Lsim


	      if(R.size()==0) //R has been emptied?
	      {
		penult=L.end(); penult--; penult--;

		//check the similarity
		if(Lsim.back() > l_agg_thres)
		{
		  //Agglomerate clusters
		  agglomerate_clusters(L.back(),*penult);
		  C.push_back(L.back()); //add the cluster to C
		}
		else
		{
		  //Save in C
		  C.push_back(*penult);
		  C.push_back(L.back());
		}
		break;//end main while
	      }
	    }
	  }
	  else
	  {
	    //Add the clusters to C (separate clusters)
	    itEnd = L.end();
	    for(it=L.begin(); it!=itEnd; it++)
	      C.push_back(*it);
	    break;//end main while
	  }
	}
      }
      else //A RNN
	RNNfound = true;
    }

    //RNN found
    if( RNNfound || ((Si.size() == 0)&&(So.size() == 0)))
    {



      if(Lsim.back() > l_agg_thres) //can they be agglomerated?
      {
	//Agglomerate clusters
	penult=L.end(); penult--; penult--;
	agglomerate_clusters(L.back(),*penult);


	insert_element(f_map,R,L.back());


	L.pop_back();
	L.pop_back();

	//delete similarities
	Lsim.pop_back();
	if(Lsim.size()>=1)
	  Lsim.pop_back();


	if(L.size() == 1)
	{
	  //Get nearest neighbor
	  get_nn(L.back(),R,sim,iNN);


	  //DEBUG
	  ndist += R.size();


	  //Add the NN chain
	  //Add to the NN chain
	  L.push_back(*iNN); //add to L
	  erase_element(f_map, R, iNN); //delete from R
	  Lsim.push_back(sim);//add to Lsim



	  if(R.size()==0) //R has been emptied?
	  {
	    penult = L.end(); penult--; penult--;
	    //check the similarity
	    if(Lsim.back() > l_agg_thres)
	    {
	      //Agglomerate clusters
	      agglomerate_clusters(L.back(),*penult);
	      C.push_back(L.back()); //add the cluster to C

	    }
	    else
	    {
	      //Save in C
	      C.push_back(*penult);
	      C.push_back(L.back());
	    }

	    break;
	  }
	}
      }
      else //discard this chain
      {
	//Add the clusters to C (separate clusters)
	itEnd = L.end();
	for(it=L.begin(); it!=itEnd; it++)
	  C.push_back(*it);

	L.clear();
      }
    }

    //Do we need to start a new chain?
    if( L.size() == 0 )
    {
      //Initialize a new chain
      Lsim.clear();

      //random point
      srand( time(NULL) );
      rp= rand() % R.size(); //random point


      //Add to L
      it = R.begin();
      advance(it,rp);
      L.push_back(*it);



      //R\rp -> delete the cluster in R and mark as erased in f_map
      erase_element(f_map, R, it);


      //First iteration
      if( R.size() > 0 )
      {
	//Get nearest neighbor
	get_nn(L.front(),R,sim,iNN);



	//DEBUG
	ndist += R.size();


	//Add to the NN chain
	L.push_back(*iNN); //add to L
	erase_element(f_map, R, iNN); //delete from R
	Lsim.push_back(sim);//add to Lsim



	//Only two clusters?
	if(R.size()==0)
	{
	  penult=L.end(); penult--; penult--;

	  //check the similarity (last element)
	  if(Lsim.back() > l_agg_thres)
	  {
	    //Agglomerate clusters
	    agglomerate_clusters(L.back(),*penult);
	    C.push_back(L.back());
	  }
	  else
	  {
	    //Save in C separately
	    C.push_back(*penult);
	    C.push_back(L.back());
	  }
	  break;//end main while
	}
      }
      else //R is empty
      {
	if(L.size() == 1)
	  C.push_back(L.front());
      }
    }
  }

  //Chain C contains all the clusters
  nc=C.size();//number of clusters

  if(nc > 0)
    dim = (C.front()).centroid.size(); // take the dimension from the fitst element in C

  cluster_centre.clear();//delete the content

  vector <double> aux_centroid(dim);

  itEnd = C.end();
  unsigned int c_label = 0, num_labels = 0, s;

  for(it = C.begin(); it!=itEnd; it++)
  {
    //convert from vcl to vnl
    for(s=0;s<dim;s++)
      aux_centroid[s]=(*it).centroid[s];

    //add centroid
    cluster_centre.push_back(aux_centroid);

    num_labels+=(*it).data_index.size();

    for(s=0; s<(*it).data_index.size(); s++)
      labels[(*it).data_index[s]]=c_label;

    c_label++;
  }

  //Were all the points asigned?.
  if(num_labels != n )
  {
    cout << "Warning: all the points were not assigned to a cluster!" << endl;
    cout << "Num. Labels = " << num_labels << " Num. Points = " << n << endl;
  }


  return nc;
}

//Get options from the command line.
int get_options(int argc,char *argv[])
{
	if(argc < 2)
	{
		cout << "Not enough arguments, please use -help for usage information.\n";
		exit(0);
	}
	if (strcmp(argv[1],"-help") == 0)
	{
		cout << "fast-RNN" << endl;

      		cout << "Usage: rnn [options]" << endl
	   	<< "-help: Prints this information" << endl
	   	<< "-data: path to the data file" <<endl
	   	<< "-dim: dimension of data" << endl
	   	<< "-comp: component for slicing" << endl
	   	<< "-n: number of vectors" << endl
	   	<< "-thres: threshold for the agglomertive clustering (default 0.4)" << endl
	  	<< "-e: epsilon (parameter for slicing)" << endl
	   	<< "-out: path to the out file" << endl;

		exit(0);
	}

	if(argc < 15)
	{
		cout << "Not enough arguments, please use -help for usage information.\n";

		exit(0);
	}
	else
	{
		for (int i = 1; i < argc; i=i+2)
		{
			if (i + 1 != argc)
			{
				if (strcmp(argv[i],"-data") == 0)
					data_file_path = argv[i + 1];
				else if (strcmp(argv[i],"-dim") == 0)
					dim=atoi(argv[i + 1]);
				else if (strcmp(argv[i],"-comp") == 0)
                    			COMP = atoi(argv[i + 1]);
				else if (strcmp(argv[i],"-n") == 0)
                    			n = atoi(argv[i + 1]);
				else if (strcmp(argv[i],"-thres") == 0)
				{
                    			thres_euclid = atof(argv[i + 1]);
					thres = thres_euclid * thres_euclid;
				}
				else if (strcmp(argv[i],"-e") == 0)
                    			ep = atof(argv[i + 1]);
				else if (strcmp(argv[i],"-out") == 0)
                    			out_file_path = argv[i + 1];
				else
				{
                    			std::cout << "Not enough or invalid arguments, please try again.\n";
                    			exit(0);
				}
			}
            	}

	}
  return 1;
}

int save_data(vector <unsigned int> &IDs,vector< vector<double> > &centres, float time, unsigned long num_distances)
{
  unsigned int i=0,k=0;

  //Out file
  ofstream ffile(out_file_path);

  if( ffile.is_open() )
  {
    //Write time
    ffile << time; ffile << endl;
    //Write num_distances
    ffile << num_distances; ffile << endl;

    //Write labels first
    for(k = 0; k < n; k++)
    {
      ffile << IDs[k];
      ffile << endl;
    }

    //Write number of clusters
    ffile << nc << endl;

    //Write centers
    for( i=0; i < nc; i++)
    {
      for(k=0; k<dim; k++)
	ffile << centres[i][k] << " ";

      ffile << endl;

    }
  }
  else
  {
    cout << "Error writing results file." << endl;
    exit(-1);
  }

  //close file
  ffile.close();
  return 1;
}


int read_data(matrix_data& data)
{
  unsigned int i=0,k=0;
  //Vector to read each descriptor
  vector <double> descriptor(dim);

  //Read data from file
  ifstream ffile(data_file_path);

  if( ffile.is_open() )
  {
    for(k = 0; k < n; k++)
    {
      for(i = 0; i < dim; i++)
	ffile >> descriptor[i];

      //Insert in data
      data.push_back(descriptor);

    }
  }
  else
  {
    cout << "Error reading data file. Please, check whether it exists or not." << endl;
    exit(-1);
  }

  //close file
  ffile.close();
  return 1;
}


#endif
