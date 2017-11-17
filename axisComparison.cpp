#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <windows.h>
#include <stdio.h>
#include <math.h>
#include <iomanip>   // format manipulation
#include <direct.h>
#include <string>     // std::string, std::stoi
#include <conio.h>
#include <cstdlib>
#include <stdlib.h>
#include "Eigen/Dense"

#include "include/protein.h"
#include "include/skeleton_overall.h"
#include "include/MRC.h"
#include "include/axis.h"

#define SSTR( x ) dynamic_cast< std::ostringstream & >(( std::ostringstream() << std::dec << x ) ).str()

using namespace std;
using Eigen::MatrixXd;
int regionQuery(int tempGroup[], int pointPlace, Map mrc, int groupNum, bool visited[], double threshold);
int regionQueryIt(int tempGroup[], int pointPlace, Map mrc, int groupNum, bool visited[], int groupCurrent[], int &num, bool notVisited[], double threshold);
int expandGroup(int tempGroup[], bool visited[], Map mrc, int groupCurrent[], int num, int groupNumber, bool notVisited[], double threshold);
void linearFit(Map mrc, int tempGroup[], int groupNumber, int total[], string path, double stepSize, bool one, int &currHel, Protein pdb, int helixOffset, int numSplit, vector<int>groupToSplitArr, bool acute);
void outputPoints(Map mrc, string path, double threshold);

int returnedNum = 0;
ofstream out;
vector<Axis> helTraceArray;

int main(int argc, char *argv[])  {
	Protein pdb;		//protein structure in xyz coordinate
	Protein pdb2;
	string oFile, path, pdbID, helixName, pdbFileName, mrcFileName, pdbOutName, mrcOutName, chainID, option, outputFilename, id_str, helBaseName, helEndName, shtBaseName, shtEndName, tempera, tempera2;
    string pathFolder = "";
	int i, j, just_a_number_holder;
    float rise = 1.5;					//the rise of the axis
	short radius = 6;
	double sradius = 2.5;
	int currHel = -1;
	int currSht = -1;
	int offset2 = 0;
    bool one = false;
    bool one2 = false;
    bool can = false;
    bool outExist = false;
    bool noHelix = false;
    bool noSheet = false;
    bool touched = false;
    bool mrcBool = false;
    int num = 0;
    ifstream inFile7;
    string pdbormrc;
    Map mrc;
    string file_helix_name = "";
    int groupNumber = 0;
    vector<Coordinate> axis;
    ofstream outputFile;
    int number_of_hel;
    double stepSize;
    stepSize = 0.1;
    vector<Axis> helTrueArray;
    vector<Axis> acuteHelix;
    vector<Axis> acuteSplitHelices;
    vector<int> splitNum;
    bool acute = false;
    vector<bool> shortHelix;
    vector<int> groupToSplitArr;
    string nextChain = "";
    string currentChain = "";
    //bool shortHelix[pdb.hlces.size()+pdb.sheets.size()];
    int helixOffset = 0;
    int strandOffset = 0;
    int numSplit = 0;
    bool moreChains = true;
    int chainNum = 0;
    bool uselessVar = false;
    string uselessString = "";
    int numtraceHels = 0;
    bool n = false;
ofstream outCoordinates100;

    ///open files and separate if necessary
    if(argc == 5)
    {
        pdbFileName = argv[1]; //"C:/Users/paulh/Desktop/DeskTop-5-18/Research/1FLP/1FLP";
        int longNameOffset = 1;
        while(pdbFileName.substr(pdbFileName.length()-longNameOffset,1) != "\\")
            longNameOffset++;
        longNameOffset--;
        //cout << longNameOffset;
        pdbFileName = pdbFileName + ".pdb";
        path = pdbFileName.substr(0, pdbFileName.length()-(longNameOffset+4));
        //cout << "\n" << path << endl;
        pdbID = pdbFileName.substr(pdbFileName.length()-(longNameOffset+4), longNameOffset);
        tempera = argv[2];//"C:/Users/paulh/Desktop/DeskTop-5-18/Research/1FLP/1FLP_altHLX_0";
        pdbormrc = tempera + ".mrc";
        inFile7.open(pdbormrc.c_str());
        if(inFile7)
        {
            one = true;
            currHel = 1;
            helBaseName = tempera.substr(0, tempera.length()-2);
            helEndName = tempera.substr(tempera.length()-2+1,tempera.length()-(tempera.length()-2+1));
            mrcBool = true;
            mrc.read( tempera + ".mrc");
            //mrc.printInfo();
            mrc.normalize();
            mrc.filterize(float(.0001/mrc.hdr.amax));
            //mrc.buildGradient(0);
            //mrc.buildTensor();
            //mrc.buildThickness(0);
            //mrc.EDT();
            //mrc.DR();
        }

        if(tempera != "Empty" && mrcBool != true)
        {
            for(int y = 0; y < tempera.length(); y++)
            {
                if(tempera.substr(tempera.length()-y, 1) == "0")
                {
                    currHel = 0;
                    helBaseName = tempera.substr(0, tempera.length()-y);
                    helEndName = tempera.substr(tempera.length()-y+1,tempera.length()-(tempera.length()-y+1));
                    y = tempera.length();

                }
                if(tempera.substr(tempera.length()-y, 1) == "1")
                {
                    currHel = 1;
                    helBaseName = tempera.substr(0, tempera.length()-y);
                    helEndName = tempera.substr(tempera.length()-y+1,tempera.length()-(tempera.length()-y+1));
                    y = tempera.length();
                    one = true;
                }
            }
        }
        else
        {
            noHelix = true;
        }
        one = true;///HARD CODED
        tempera2 = argv[3];//"Empty";
        if(tempera2 != "Empty")
        {
            string fileName = tempera2.c_str();
            fileName = fileName + ".pdb";
            ifstream cFile;
            cFile.open(fileName.c_str());
            pdb2.read(fileName.c_str(), uselessVar, uselessString);
            size_t lines_count2 =0;
            string templine3;
            while (std::getline(cFile , templine3))
            {
                if(templine3.substr(0,4) == "ATOM")
                    ++lines_count2;
            }

            cFile.close();
            int vs = 0;
            double total2 = 0;
            Coordinate last2;
            Coordinate x;
            int temp2 = 0;
            bool OneFile = true;
            while(vs < lines_count2-1)
            {
                if(vs != 0)
                {
                    last2 = x;
                    x = pdb2.AAs[0].atoms[vs].coord;
                    total2 = sqrt((last2.x-x.x)*(last2.x-x.x) + (last2.y-x.y)*(last2.y-x.y) + (last2.z-x.z)*(last2.z-x.z));
                    if(total2 > 1.5)
                    {
                        OneFile = false;
                        vs = lines_count2;
                    }
                }
                else
                    x = pdb2.AAs[0].atoms[vs].coord;
                vs++;
            }
            if(OneFile == true)
            {
                for(int y = tempera2.length()-4; y >0; y--)
                {
                    if(tempera2.substr(y, 1) == "/")
                    {
                        offset2 = y;
                        break;
                    }
                }
                for(int y = 0; y < tempera2.length()-4-offset2; y++)
                {
                    if(tempera2.substr(tempera2.length()-y, 1) == "0")
                    {
                        currSht = 0;
                        shtBaseName = tempera2.substr(0, tempera2.length()-y);
                        shtEndName = tempera2.substr(tempera2.length()-y+1,tempera2.length()-(tempera2.length()-y+1));
                        y = tempera2.length();
                        touched = true;
                    }
                    if(tempera2.substr(tempera2.length()-y, 1) == "1")
                    {
                        currSht = 1;
                        shtBaseName = tempera2.substr(0, tempera2.length()-y);
                        shtEndName = tempera2.substr(tempera2.length()-y+1,tempera2.length()-(tempera2.length()-y+1));
                        y = tempera2.length();
                        one2 = true;
                        touched = true;
                    }
                }
            }
            //wait for python code to delete old files
            Sleep(1000);
            //handles when input is only one file
            if(touched == false)
            {
                //declaring variables
                ifstream inFile;
                ofstream outFile;
                ofstream outFile2;
                int sheetit = 1;
                int it2 = 1;
                int temp1;
                double tempX = 0;
                double tempY = 0;
                double tempZ = 0;
                string temp2;
                double total = 0;
                int bt = 0;
                bool firstPass = true;
                Coordinate x;
                Coordinate last;
                //opening initial input and output files
                //inFile.open(tempera2.c_str(), ios::in);

                tempera2 = tempera2 + ".pdb";
                pdb2.read(tempera2.c_str(), uselessVar, uselessString);
                x = pdb2.AAs[0].atoms[0].coord;
                string outputFilename = path + "/output/sheet_" + toString(sheetit) + ".pdb";
                outFile.open(outputFilename.c_str(), ios::out);
                //read in values from input file
                //inFile >> temp2 >> temp1 >> temp2 >> temp2 >> temp2 >> temp1 >> valX >> valY >> valZ;
                tempX = x.x;
                tempY = x.y;
                tempZ = x.z;
                if(it2 <= 9)
                    outFile << fixed << setprecision(3) << "ATOM      " << it2 << "  H   HOH A   1      " << tempX << "  " << tempY << "  " << tempZ << endl;
                else
                    outFile << fixed << setprecision(3) << "ATOM     " << it2 << "  H   HOH A   1      " << tempX  << "  " << tempY << "  " << tempZ << endl;

                //check number of lines in file
                ifstream aFile (tempera2.c_str());
                size_t lines_count =0;
                string line;
                while (std::getline(aFile , line))
                {
                    if(line.substr(0,4) == "ATOM")
                        ++lines_count;
                }

                aFile.close();

                int lineCount[100];
                int strandCount = 0;
                for(int i = 0; i < 100; i++)
                {
                    lineCount[i] = 0;
                }
                ifstream bFile (tempera2.c_str());
                for(int i = 0; i < lines_count; i++)
                {
                    if(i!=0)
                    {
                        last = x;
                        x = pdb2.AAs[0].atoms[i].coord;
                        total = sqrt((last.x-x.x)*(last.x-x.x) + (last.y-x.y)*(last.y-x.y) + (last.z-x.z)*(last.z-x.z));
                    }
                    else
                    {
                        x = pdb2.AAs[0].atoms[i].coord;
                    }
                    if(total > 1.5)
                    {
                        strandCount++;
                        lineCount[strandCount]++;
                    }
                    else
                        lineCount[strandCount]++;
                }
                bFile.close();
                strandCount = 0;
                int currentCompare = 0;
                int div = 0;
                int remainder = 0;
                bool tooManyLines = false;
                //for(int j = 0; j < 10; j++)
                //    cout << lineCount[j] <<  endl;
                //cout << lines_count;
                while(bt <lines_count)
                {

                    if(bt == 0)
                    {
                        if(lineCount[strandCount] > 20)
                        {
                            tooManyLines = true;
                            div = lineCount[strandCount]/20;
                            remainder = lineCount[strandCount]%20;
                            if(remainder != 0)
                            {
                                currentCompare = div+1;
                                remainder--;
                            }
                            else
                                currentCompare = div;
                        }
                        else
                            currentCompare = it2+1;
                    }
                    //cout << it2 << " " << bt << " " << lineCount[strandCount] << "\n";
                    if(it2 == lineCount[strandCount])
                    {
                        //cout << "test" << lineCount[strandCount];
                        strandCount++;
                        if(lineCount[strandCount] > 20)
                        {
                            div = lineCount[strandCount]/20;
                            remainder = lineCount[strandCount]%20;
                            if(remainder != 0)
                            {
                                currentCompare = div+1;
                                remainder--;
                            }
                            else
                                currentCompare = div;
                        }
                        else
                            currentCompare = it2+1;
                    }

                    bt++;
                    it2++;
                    total = 0;

                    //set up total comparison
                    last = x;
                    x = pdb2.AAs[0].atoms[bt].coord;
                    total = sqrt((last.x-x.x)*(last.x-x.x) + (last.y-x.y)*(last.y-x.y) + (last.z-x.z)*(last.z-x.z));
                    //cout << total << endl;
                    if(total > 1.5 && it2 != 2)
                    {
                        sheetit++;
                        it2=1;
                        if(firstPass == true)
                        {
                            outFile.close();
                            firstPass = false;
                        }
                        else
                            outFile2.close();
                        if(bt != lines_count)
                        {
                            //cout << "bt:" << bt <<endl;
                            outputFilename = path + "/output/sheet_" + toString(sheetit) + ".pdb";
                            outFile2.open(outputFilename.c_str(), ios::out);
                            tempX = x.x;
                            tempY = x.y;
                            tempZ = x.z;
                            if(it2 < 10)
                                outFile2 << fixed << setprecision(3) << "ATOM      " << it2 << "  H   HOH A   1      " << tempX  << "  " << tempY << "  " << tempZ << endl;
                            else if(it2 > 9 && it2 < 100)
                                outFile2 <<fixed << setprecision(3) << "ATOM     " << it2 << "  H   HOH A   1      " << tempX  << "  " << tempY << "  " << tempZ << endl;
                            else
                                outFile2 <<fixed << setprecision(3) << "ATOM    " << it2 << "  H   HOH A   1      " << tempX  << "  " << tempY << "  " << tempZ << endl;
                        }
                    }
                    else
                    {
                        if(true)
                        {
                            if(firstPass == true)
                            {

                                tempX = x.x;
                                tempY = x.y;
                                tempZ = x.z;
                                if(it2 < 10)
                                    outFile <<fixed << setprecision(3) << "ATOM      " << it2 << "  H   HOH A   1      " << tempX  << "  " << tempY << "  " << tempZ << endl;
                                else if(it2 > 9 && it2 < 100)
                                    outFile <<fixed << setprecision(3) << "ATOM     " << it2 << "  H   HOH A   1      " << tempX  << "  " << tempY << "  " << tempZ << endl;
                                else
                                    outFile <<fixed << setprecision(3) << "ATOM    " << it2 << "  H   HOH A   1      " << tempX  << "  " << tempY << "  " << tempZ << endl;
                            }
                            else
                            {
                                tempX = x.x;
                                tempY = x.y;
                                tempZ = x.z;
                                if(it2 < 10)
                                    outFile2 <<fixed << setprecision(3) << "ATOM      " << it2 << "  H   HOH A   1      " << tempX  << "  " << tempY << "  " << tempZ << endl;
                                else if(it2 > 9 && it2 < 100)
                                    outFile2 <<fixed << setprecision(3) << "ATOM     " << it2 << "  H   HOH A   1      " << tempX  << "  " << tempY << "  " << tempZ << endl;
                                else
                                    outFile2 <<fixed << setprecision(3) << "ATOM    " << it2 << "  H   HOH A   1      " << tempX  << "  " << tempY << "  " << tempZ << endl;
                            }

                            if(tooManyLines == false)
                            {
                                cout << it2 << " " << currentCompare << " " << bt << "\n";
                                currentCompare++;
                            }
                            else
                            {
                                if(remainder != 0)
                                {
                                    currentCompare += div+1;
                                    remainder--;
                                }
                                else
                                    currentCompare += div;
                            }

                        }

                    }
                }

                tempera2 = path + "/output/sheet_1";
                for(int y = 0; y < tempera2.length()-4; y++)
                {
                    //cout << tempera2.substr(tempera2.length()-y, 1);
                    if(tempera2.substr(tempera2.length()-y, 1) == "1")
                    {
                        //cout << "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbb";
                        currSht = 1;
                        shtBaseName = tempera2.substr(0, tempera2.length()-y);
                        shtEndName = tempera2.substr(tempera2.length()-y+1,tempera2.length()-(tempera2.length()-y+1));
                        y = tempera2.length();
                        one2 = true;
                    }
                }
                currSht = 1;
                //shtBaseName = path + "SHT";
                //shtEndName = ".pdb";
                one2 = true;
                //cout << "aaaaaaaaa" << shtBaseName << "  " << shtEndName << endl;
                outFile2.close();
            }
            ///end of truing to split files
        }
        else
        {
            noSheet = true;
        }

        oFile = argv[4];//"Empty";
        outExist = true;
        if(oFile == "Empty")
        {
            outExist = false;
        }
    }
    else
    {
        cout<<"usage: "<< argv[0] <<" trueStructureFileLocation "<<" DetectedHelixFileLocation"<<" DetectedStickFileLocation "<< " LocationToCreateOutputFile "<<endl; //argv[0] is the program name
        cout<< "Only trueStructureFileLocation is required. Replace argument with 'Empty' if not using it." << endl;
        exit(1);
    }

    while(moreChains == true)
    {
        moreChains = false;
        if(chainNum == 0)
            pdb.read(pdbFileName.c_str(), moreChains, nextChain, "");
        else
            pdb.read(pdbFileName.c_str(), moreChains, nextChain, currentChain);
        currentChain = nextChain;

        for(i = 0; i < pdb.hlces.size(); i++)
        {
            shortHelix.push_back(false);
        }

        ///create and output True Helix axis from .pdb file
        for (i = 0 ; i < pdb.hlces.size(); i++)
        {
            Coordinate pdbEnd1;
            Coordinate pdbEnd2;
            just_a_number_holder = i + 1 + helixOffset;
            stringstream ss;
            ss << just_a_number_holder;
            id_str = ss.str();
            ///id_str = SSTR(i);
            /** id_str = static_cast<ostringstream*>( &(ostringstream() << i) )->str(); */

            outputFilename = path + "output/trueHelix" + id_str + ".pdb";
            ss.str( std::string() );
            ss.clear();
            axis.clear();
            pdb.setAxisSegments(pdb.hlces[i].startIndx, pdb.hlces[i].endIndx, rise, axis, shortHelix, i);
            //cerr << "Helix " << i+1 << endl;
            pdbEnd1 = pdb.AAs[pdb.hlces[i].startIndx].atoms[pdb.getAtomIndx(pdb.hlces[i].startIndx, "CA")].coord;
            pdbEnd2 = pdb.AAs[pdb.hlces[i].endIndx].atoms[pdb.getAtomIndx(pdb.hlces[i].endIndx, "CA")].coord;
            Axis myAxis;
            myAxis.axisPoints = axis;
            myAxis.addTrueEnds(pdbEnd1, pdbEnd2);
            myAxis.catmullRom(stepSize);
            myAxis.printAsPnts(outputFilename);
            //split helix if the angle is too bent
            //myAxis.angle() < 1.5708 && myAxis.angle() > .68
            if(myAxis.angle())
            {
                splitNum.push_back(i);

                acute = true;
                acuteHelix.push_back(myAxis);
                //splitNum.push_back(pdb.hlces.size()-1+numSplit);

                //first split
                Axis firstAxis;
                Axis secondAxis;
                firstAxis = myAxis;
                secondAxis = myAxis;

                firstAxis.splitFirst();
                secondAxis.splitSecond();

                acuteSplitHelices.push_back(firstAxis);
                acuteSplitHelices.push_back(secondAxis);

                helTrueArray.push_back(firstAxis);

                //outputFilename = path + "output/TEST1.pdb";
                //firstAxis.printAsPnts(outputFilename);
                //outputFilename = path + "output/TEST2.pdb";
                //secondAxis.printAsPnts(outputFilename);

                numSplit++;
                //cout << numSplit << endl;
            }
            else if(myAxis.axisPoints.size() > 500)
            {
                splitNum.push_back(i);

                acute = true;
                acuteHelix.push_back(myAxis);
                //splitNum.push_back(pdb.hlces.size()-1+numSplit);

                //first split
                Axis firstAxis;
                Axis secondAxis;
                firstAxis = myAxis;
                secondAxis = myAxis;

                firstAxis.splitFirstHalf();
                secondAxis.splitSecondHalf();

                acuteSplitHelices.push_back(firstAxis);
                acuteSplitHelices.push_back(secondAxis);

                helTrueArray.push_back(firstAxis);

                numSplit++;

            }
            else
               helTrueArray.push_back(myAxis);
        }

        helixOffset += pdb.hlces.size();
        chainNum++;
    }

    moreChains = true;
    chainNum = 0;

    while(moreChains == true)
    {
        moreChains = false;
        if(chainNum == 0)
            pdb.read(pdbFileName.c_str(), moreChains, nextChain, "");
        else
            pdb.read(pdbFileName.c_str(), moreChains, nextChain, currentChain);
        currentChain = nextChain;

        for(i = 0; i < pdb.sheets.size(); i++)
        {
            shortHelix.push_back(false);
        }

        ///create and output True Strand axes from .pdb file
        for (i = 0 ; i < pdb.sheets.size(); i++)
        {
            Coordinate pdbEnd1;
            Coordinate pdbEnd2;
            Coordinate last;
            just_a_number_holder = i + 1 + helixOffset + strandOffset;
            stringstream ss;
            ss << just_a_number_holder;
            id_str = ss.str();
            ///id_str = SSTR(i);
            /** id_str = static_cast<ostringstream*>( &(ostringstream() << i) )->str(); */

            outputFilename = path + "/output/trueSheet" + id_str + ".pdb";
            ss.str( std::string() );
            ss.clear();
            axis.clear();

            pdbEnd1 = pdb.AAs[pdb.sheets[i].endIndx].atoms[pdb.numOfAtoms(pdb.sheets[i].endIndx)-2].coord;
            pdbEnd2 = pdb.AAs[pdb.sheets[i].endIndx+1].atoms[pdb.getAtomIndx(pdb.sheets[i].endIndx, "N")].coord;
            last.x = (pdbEnd1.x+pdbEnd2.x)/2;
            last.y = (pdbEnd1.y+pdbEnd2.y)/2;
            last.z = (pdbEnd1.z+pdbEnd2.z)/2;
            pdb.setAxisSegments2(pdb.sheets[i].startIndx, pdb.sheets[i].endIndx, rise, axis, last, shortHelix, i+helixOffset);
            //cout << pdbEnd1.x << " " << pdbEnd1.y << " " << pdbEnd1.z << "\n";
            //cout << pdbEnd1.x << " " << pdbEnd1.y << " " << pdbEnd1.z << "\n\n";
            Axis myAxis;
            myAxis.axisPoints = axis;
            //myAxis.addTrueEnds(pdbEnd1, pdbEnd2);
            myAxis.catmullRom(stepSize);
            myAxis.printAsPnts(outputFilename);
            helTrueArray.push_back(myAxis);

        }

        strandOffset += pdb.sheets.size();
        chainNum++;
    }


    ifstream inFile_tracer;
    vector<Coordinate> helPoints;

        int offset = 0;
        stringstream sstring;
        sstring << currHel;
        string current_helix_str = sstring.str();
        if(mrcBool == false)
        {
            file_helix_name = helBaseName  + current_helix_str + helEndName + ".pdb";
            inFile_tracer.open(file_helix_name.c_str());
            ///output trace helices
            while(inFile_tracer) {
                Coordinate tmppnt;
                string line;
                sstring.str(std::string());
                sstring.clear();
                while(getline(inFile_tracer, line))
                {
                    line.erase(0,30);
                    sstring << line;
                    sstring >> tmppnt.x  >> tmppnt.y >> tmppnt.z;

                    //cerr << tmppnt.x << " " << tmppnt.y << " " << tmppnt.z << endl;
                    helPoints.push_back(tmppnt);
                    sstring.str(std::string());
                    sstring.clear();
                }
                Axis myAxis;
                myAxis.axisPoints = helPoints;
                myAxis.catmullRom(stepSize);
                if(one == false)
                {
                    if (currHel <= helixOffset)
                    {
                        outputFilename = path + "/output/traceHelix" + current_helix_str + ".pdb";
                    }
                    else
                        outputFilename = path + "/output/traceSheet" + current_helix_str + ".pdb";
                }
                if(one == true)
                {
                    if (currHel <= helixOffset+1)
                    {
                        outputFilename = path + "/output/traceHelix" + current_helix_str + ".pdb";
                    }
                    else
                        outputFilename = path + "/output/traceSheet" + current_helix_str + ".pdb";
                }

                //cerr << "Now writing to " << outputFilename <<endl;
                myAxis.printAsPnts2(outputFilename);
                helTraceArray.push_back(myAxis);
                helPoints.clear();
                inFile_tracer.close();

                numtraceHels++;
                currHel++;
                sstring << currHel;
                current_helix_str = sstring.str();
                file_helix_name = helBaseName  + current_helix_str + helEndName + ".pdb";
                //cout << "Helix #" << current_helix_str << endl << "-----------------" << endl;
                inFile_tracer.open(file_helix_name.c_str());
            }
        }
        else    ///if file is MRC, group and perform least squares
        {
            cout << "Grouping voxels ..." << endl;

            ///grouping voxels into groups that touch each other
            int *tempGroup;
            tempGroup = (int *)malloc(sizeof(int)*(mrc.numSlcs()*mrc.numCols()*mrc.numRows()));
            bool *visited;
            visited = (bool *)malloc(sizeof(bool)*(mrc.numSlcs()*mrc.numCols()*mrc.numRows()));
            bool *notVisited;
            notVisited = (bool *)malloc(sizeof(bool)*(mrc.numSlcs()*mrc.numCols()*mrc.numRows()));
            //bool notVisited[mrc.numSlcs()*mrc.numCols()*mrc.numRows()];
            int *groupCurrent;
            groupCurrent = (int *)malloc(sizeof(int)*(10000));
            int total[500];

            for(int i = 0; i < mrc.numSlcs()*mrc.numCols()*mrc.numRows(); i++)
            {
                if(i < 500)
                    total[i] = 0;
                if(i < 10000)
                    groupCurrent[i] = -1;
                tempGroup[i] = -1;
                visited[i] = false;
                notVisited[i] = false;
            }
            ///recursive method
            /*
            int numReturned = -1;
            for (long k=0; k<mrc.numSlcs(); k++)
                for (long j=0; j<mrc.numCols(); j++)
                    for (long i=0; i<mrc.numRows(); i++)
                    {
                        if (mrc.cube[i][j][k]>0 && visited[(k*mrc.numCols()*mrc.numRows()+j*mrc.numRows()+i)] == false)
                        {
                            numReturned = regionQuery(tempGroup,(k*mrc.numCols()*mrc.numRows()+j*mrc.numRows()+i), mrc, groupNumber, visited);
                            cout << numReturned << " ";
                            groupNumber++;
                            returnedNum = 0;
                        }
                    }
                    cout << "\n\n";
            */
            ///figure out the values of cube[][][]
            double low = 1;
            double high = 0;
            double temp = 0;
            double middle = 0;
            int iter = 0;
            for (long k=0; k<mrc.numSlcs(); k++)
                for (long j=0; j<mrc.numCols(); j++)
                    for (long i=0; i<mrc.numRows(); i++)
                    {
                        temp = mrc.cube[i][j][k];
                        if(temp != 0)
                        {
                            if(temp < low)
                                low = temp;
                            if(temp > high)
                                high = temp;
                        }
                    }
            middle = (high-low)/2;
            middle += low;

            int tempJ = 0;
            int tempK = 0;
            int tempI = 0;
            double threshold = .1/mrc.hdr.amax;//middle;

            ///iterative method
            int numReturned = -1;
            int groupSize = 0;
            for (long k=0; k<mrc.numSlcs(); k++)
                for (long j=0; j<mrc.numCols(); j++)
                    for (long i=0; i<mrc.numRows(); i++)
                    {
                        //cout << k << endl;
                        if (mrc.cube[i][j][k]> threshold && visited[(k*mrc.numCols()*mrc.numRows()+j*mrc.numRows()+i)] == false)
                        {
                            numReturned = regionQueryIt(tempGroup,(k*mrc.numCols()*mrc.numRows()+j*mrc.numRows()+i), mrc, groupNumber, visited, groupCurrent, num, notVisited, threshold);
                            ///if numReturned is too small, delete it TODO
                            groupSize = expandGroup(tempGroup, visited, mrc, groupCurrent, numReturned, groupNumber, notVisited, threshold);
                            for(int y = 0; y < 10000; y++)
                            {
                                if(groupCurrent[y] == -1)
                                {
                                    y = 10000;
                                }
                                else
                                {
                                    total[groupNumber]++;
                                    tempGroup[groupCurrent[y]] = groupNumber;
                                    groupCurrent[y] = -1;
                                }

                            }
                            if(groupSize != 0)
                                groupNumber++;
                            returnedNum = 0;
                            num = 0;
                            //cout << "YYYYYYYYESSSSSSSSSSSSS" << endl;
                        }
                    }
            //cout << groupNumber << endl;
            ///split group according to acute helices
            //find which from to split from true line
            if(acute == true)
            {
                for(int z = 0; z < numSplit; z++)
                {

                    double distance = 0;
                    double smallestDist = 99999999;
                    int correspondingGroup = 0;
                    double cubeX = 0;
                    double cubeY = 0;
                    double cubeZ = 0;
                    int groupToSplit = -2;
                    Coordinate point = acuteHelix[z].firstPoint();
                    for (long k=0; k<mrc.numSlcs(); k++)
                        for (long j=0; j<mrc.numCols(); j++)
                            for (long i=0; i<mrc.numRows(); i++)
                            {
                                if(mrc.cube[i][j][k] > threshold)
                                {
                                    cubeX = i*mrc.apixX+mrc.hdr.xorigin;
                                    cubeY = j*mrc.apixY+mrc.hdr.yorigin;
                                    cubeZ = k*mrc.apixZ+mrc.hdr.zorigin;
                                    distance = sqrt(pow((point.x-cubeX),2)+pow((point.y-cubeY),2)+pow((point.z-cubeZ),2));

                                    if(distance < smallestDist)
                                    {
                                        smallestDist = distance;
                                        correspondingGroup = k*mrc.numCols()*mrc.numRows()+j*mrc.numRows()+i;
                                    }
                                }
                            }
                    groupToSplit = tempGroup[correspondingGroup];
                    groupToSplitArr.push_back(tempGroup[correspondingGroup]);
                    //cout << groupToSplit << " " << total[groupNumber] << endl;

                    //split group
                    double firstDist = 0;
                    double secondDist = 0;
                    Coordinate firstPoint = acuteSplitHelices[0+z*2].lastPoint();
                    Coordinate secondPoint = acuteSplitHelices[1+z*2].lastPoint();
                    for (long k=0; k<mrc.numSlcs(); k++)
                            for (long j=0; j<mrc.numCols(); j++)
                                for (long i=0; i<mrc.numRows(); i++)
                                {
                                    //if this is the group to split compare distance of point to regroup
                                    if(tempGroup[(k*mrc.numCols()*mrc.numRows()+j*mrc.numRows()+i)] == groupToSplit)
                                    {
                                        cubeX = i*mrc.apixX+mrc.hdr.xorigin;
                                        cubeY = j*mrc.apixY+mrc.hdr.yorigin;
                                        cubeZ = k*mrc.apixZ+mrc.hdr.zorigin;
                                        firstDist = sqrt(pow((firstPoint.x-cubeX),2)+pow((firstPoint.y-cubeY),2)+pow((firstPoint.z-cubeZ),2));
                                        secondDist = sqrt(pow((secondPoint.x-cubeX),2)+pow((secondPoint.y-cubeY),2)+pow((secondPoint.z-cubeZ),2));

                                        if(secondDist < firstDist)
                                        {
                                            tempGroup[(k*mrc.numCols()*mrc.numRows()+j*mrc.numRows()+i)] = groupNumber;
                                            total[groupToSplit]--;
                                            total[groupNumber]++;
                                        }
                                    }
                                }
                                groupNumber++;
                }
            }
            //outputPoints(mrc, path, 46);
            linearFit(mrc, tempGroup, groupNumber, total, path, stepSize, one, currHel, pdb, helixOffset, numSplit, groupToSplitArr, acute);
        }
        //end of mrc stuff

        offset = currHel-1;
        stringstream s2string;
        s2string << currSht;
        string current_sheet_str = s2string.str();
        string file_sheet_name = shtBaseName  + current_sheet_str + shtEndName + ".pdb";
        //cout << file_sheet_name;
        inFile_tracer.open(file_sheet_name.c_str());

    ///output trace strands
    while(inFile_tracer) {
        Coordinate tmppnt;
        string line;
        s2string.str(std::string());
        s2string.clear();
        while(getline(inFile_tracer, line))
        {
            line.erase(0,30);
            s2string << line;
            s2string >> tmppnt.x  >> tmppnt.y >> tmppnt.z;

            //cerr << tmppnt.x << " " << tmppnt.y << " " << tmppnt.z << endl;
            helPoints.push_back(tmppnt);
            s2string.str(std::string());
            s2string.clear();
        }
        Axis myAxis;
        myAxis.axisPoints = helPoints;
        myAxis.catmullRom(stepSize);

        outputFilename = path + "/output/traceSheet" + current_sheet_str + ".pdb";

        //cerr << "Now writing to " << outputFilename <<endl;
        myAxis.printAsPnts2(outputFilename);
        helTraceArray.push_back(myAxis);
        helPoints.clear();
        inFile_tracer.close();

        currSht++;
        s2string << currSht;
        current_sheet_str = s2string.str();
        file_sheet_name = shtBaseName  + current_sheet_str + shtEndName + ".pdb";

        //cout << "Helix #" << current_helix_str << endl << "-----------------" << endl;
        inFile_tracer.open(file_sheet_name.c_str());
    }
    if(noSheet)
        number_of_hel = currHel-1;
    else
    {
        if(noHelix)
            number_of_hel = currSht-1;
        else
            number_of_hel = currHel-1+currSht-1;
    }

    // END OF INPUT


cerr << "Beginning matching." << endl;

int matchedAlternate[helixOffset + strandOffset];
double AlternateLateral[helixOffset + strandOffset];
double AlternateLongitudinal[helixOffset + strandOffset];

int matchedHels[helixOffset + strandOffset];
double shortHels[helixOffset + strandOffset];
int splitHels[splitNum.size()+1];
double splitLateral[splitNum.size()*2];
double splitLongitudinal[splitNum.size()*2];

//int distanceFromCenter[helixOffset + strandOffset];
for (int i = 0; i < helixOffset + strandOffset; i++)
{
    AlternateLateral[i] = 99999;
    AlternateLongitudinal[i] = 99999;
    matchedAlternate[i] = 99999;
    matchedHels[i] = 99999;
    shortHels[i] = -1;
    //distanceFromCenter[i] = 99999;
}

Displacement Dis;
 double distCent;
 Displacement Dis2;
 double distCent2;
 Displacement Dis3;
 Coordinate trueCenter;
 Coordinate traceCenter;
 double minDist;
 double SecondDist = 99999;
 int secondTrace = 0;
 int minLong = 0;
 //int holder = 0;
 int it = 0;
 int itOdd = 1;

 double comparedDistance = 0;

 if(mrcBool == false)
    comparedDistance = 9;
 else
    comparedDistance = 15;
 //if(tempera == "Empty")
 //{
 //    holder = helixOffset;
 //}
                    //= holder
 for (int currTrueHel = 0; currTrueHel < helixOffset + strandOffset; currTrueHel++) //for each true helix position
 {
    trueCenter.x = (helTrueArray[currTrueHel].axisPoints[0].x + helTrueArray[currTrueHel].axisPoints[helTrueArray[currTrueHel].axisPoints.size()-1].x)/2;
    trueCenter.y = (helTrueArray[currTrueHel].axisPoints[0].y + helTrueArray[currTrueHel].axisPoints[helTrueArray[currTrueHel].axisPoints.size()-1].y)/2;
    trueCenter.z = (helTrueArray[currTrueHel].axisPoints[0].z + helTrueArray[currTrueHel].axisPoints[helTrueArray[currTrueHel].axisPoints.size()-1].z)/2;
    cerr <<endl << "True center: " << trueCenter.x << " " << trueCenter.y << " " << trueCenter.z << endl;

    SecondDist = 99999;
     minDist = 3.0;
     for (int currTraceHel = 0; currTraceHel < number_of_hel; currTraceHel++) // check each trace for nearness
     {
         cerr << "Comparing trace " << currTraceHel << endl;
         traceCenter.x = (helTraceArray[currTraceHel].axisPoints[0].x + helTraceArray[currTraceHel].axisPoints[helTraceArray[currTraceHel].axisPoints.size()-1].x)/2;
         traceCenter.y = (helTraceArray[currTraceHel].axisPoints[0].y + helTraceArray[currTraceHel].axisPoints[helTraceArray[currTraceHel].axisPoints.size()-1].y)/2;
         traceCenter.z = (helTraceArray[currTraceHel].axisPoints[0].z + helTraceArray[currTraceHel].axisPoints[helTraceArray[currTraceHel].axisPoints.size()-1].z)/2;
        distCent = getDistance(trueCenter, traceCenter);

        if(shortHelix[currTrueHel] == true && distCent < minDist)
        {
            if(distCent < 10)
            {
                matchedHels[currTrueHel] = currTraceHel;
                minDist = distCent;
                shortHels[currTrueHel] = distCent;
            }
        }
        else if((distCent < comparedDistance))
        {

             helTrueArray[currTrueHel].alignDirectionHels(&helTraceArray[currTraceHel]);
             Dis = helTrueArray[currTrueHel].findDisplacement(&helTraceArray[currTraceHel], stepSize);
             distCent = Dis.crossDisplace;
             minLong = Dis.longDisplace;
             if (distCent < minDist)
             {
                if(matchedHels[currTrueHel] != 99999 && minDist < 1.5)
                {
                    matchedAlternate[currTrueHel] = matchedHels[currTrueHel];
                    AlternateLateral[currTrueHel] = minDist;
                    AlternateLongitudinal[currTrueHel] = minLong;
                }
                matchedHels[currTrueHel] = currTraceHel;
                minDist = distCent;
             }
             else if(distCent < 1.5)
             {
                 matchedAlternate[currTrueHel] = currTraceHel;
                 AlternateLateral[currTrueHel] = distCent;
                 AlternateLongitudinal[currTrueHel] = Dis.longDisplace;
             }
        }
        //for split helices
        if(it < splitNum.size() && currTrueHel == splitNum[it])
        {
            //outputFilename = path + "output/TEST1.pdb";
            //acuteSplitHelices[itOdd].printAsPnts(outputFilename);
            //outputFilename = path + "output/TEST2.pdb";
            //helTraceArray[currTraceHel].printAsPnts(outputFilename);
            acuteSplitHelices[itOdd].alignDirectionHels(&helTraceArray[currTraceHel]);
            Dis2 =  acuteSplitHelices[itOdd].findDisplacement(&helTraceArray[currTraceHel], stepSize);
            distCent2 = Dis2.crossDisplace;

            if(distCent2 < SecondDist)
            {
                SecondDist = distCent2;
                splitHels[it] = currTraceHel;
                splitLateral[it] = Dis2.crossDisplace;
                splitLongitudinal[it] = Dis2.longDisplace;
            }

        }
         cerr << "Done comparing trace " << currTraceHel << endl;
     }
     if(it < splitNum.size() && currTrueHel == splitNum[it])
     {
         it++;
         itOdd+=2;
     }

     if (matchedHels[currTrueHel] == 99999)
     {
         cout << "There is no trace detected for true Helix/Sheet " << currTrueHel+1 << endl;
     }
}


double score1 = 0;
double score2 = 0;
int newComparison = -1;
//make sure that all strands are only matched once
for(int i = 0; i < helixOffset + strandOffset; i++)
{
    //cout << matchedHels[i] << endl;
    for(int j = 0; j < helixOffset + strandOffset; j++)
    {
        //cout << i << " " << j << " " << matchedHels[i] << " " << matchedHels[j] << endl;
        if(i != j && matchedHels[i] != 99999 && matchedHels[j] != 99999 && matchedHels[i] == matchedHels[j])
        {
            helTrueArray[i].alignDirectionHels(&helTraceArray[matchedHels[i]]);
            Dis = helTrueArray[i].findDisplacement(&helTraceArray[matchedHels[i]], stepSize);
            score1 = Dis.crossDisplace;

            helTrueArray[j].alignDirectionHels(&helTraceArray[matchedHels[j]]);
            Dis = helTrueArray[j].findDisplacement(&helTraceArray[matchedHels[j]], stepSize);
            score2 = Dis.crossDisplace;

            if(score1 < score2)
            {
                matchedHels[j] = 99999;
                newComparison = j;
            }
            else if(score1 >= score2)
            {
                matchedHels[i] = 99999;
                newComparison = i;
            }

            SecondDist = 99999;
             minDist = 3.0;
             for (int currTraceHel = 0; currTraceHel < number_of_hel; currTraceHel++) // check each trace for nearness
             {
                 cerr << "Comparing trace " << currTraceHel << endl;
                 traceCenter.x = (helTraceArray[currTraceHel].axisPoints[0].x + helTraceArray[currTraceHel].axisPoints[helTraceArray[currTraceHel].axisPoints.size()-1].x)/2;
                 traceCenter.y = (helTraceArray[currTraceHel].axisPoints[0].y + helTraceArray[currTraceHel].axisPoints[helTraceArray[currTraceHel].axisPoints.size()-1].y)/2;
                 traceCenter.z = (helTraceArray[currTraceHel].axisPoints[0].z + helTraceArray[currTraceHel].axisPoints[helTraceArray[currTraceHel].axisPoints.size()-1].z)/2;
                distCent = getDistance(trueCenter, traceCenter);

                if((distCent < comparedDistance))
                {

                     helTrueArray[newComparison].alignDirectionHels(&helTraceArray[currTraceHel]);
                     Dis = helTrueArray[newComparison].findDisplacement(&helTraceArray[currTraceHel], stepSize);
                     distCent = Dis.crossDisplace;

                     if (distCent < minDist)
                     {
                        matchedHels[newComparison] = currTraceHel;
                        minDist = distCent;
                     }
                }

            }
            i = 0;
            j = 0;
        }
    }
}




///////////////////////////////
// END MATCHING              //
// BEGIN DISTANCE COMPARISON //
///////////////////////////////
Displacement displaced;
double axisLength = 0;
ofstream outFile;
if(outExist == true)
{
    outFile.open(oFile.c_str());
    outFile << "Using axis interpolation f " << stepSize << endl;
    outFile << "Helix#True, Helix#Trace, LengthTrueAxis, LengthDetectedAxis, TwoWayDistance, CrossDisplacement, LengthDisplacement, LengthErrorProportion, Specificity, Sensitivity, F1Score" << endl << endl;
}

/// //////////////////////////
/// Accuracy print out
int tot = 0;
double spec = 0;
double sen = 0;
int sens = 0;
vector<int> specArray;
for(int i = 0; i < number_of_hel; i++)
{
    specArray.push_back(0);
}
if(outExist == true)
{
    int nonHLXAA = PrintSpecificity(&helTraceArray, &specArray, pdb, outFile, numtraceHels);
    int nonSTRAA = PrintSpecificityStrands(&helTraceArray, &specArray, pdb, outFile, numtraceHels);
}

int it2 = 0;
bool matchedTwice = false;
for (int currentHel = 0; currentHel <  helixOffset + strandOffset; currentHel++)
{
    matchedTwice = false;
    for(int i = 0; i < helixOffset + strandOffset; i++)
    {
        if(matchedHels[i] == matchedHels[currentHel] && i != currentHel)
            matchedTwice = true;
    }
    //normal matching
    //cout << one << " " << currHel << endl;
    if(one == false)
    {
        if(matchedHels[currentHel] > 100 || (noSheet == true && currentHel >= helixOffset))
        {
            if(currentHel < helixOffset)
                cout << endl << "trueHelix #" << currentHel+1 << endl << "------------------------------------" << endl;
            else
                cout << endl << "trueSheet #" << currentHel+1 << endl << "------------------------------------" << endl;
        }
        else
        {
            if(currentHel < helixOffset)
            {
                cout << endl << "trueHelix #" << currentHel+1 << " --> traceHelix" << matchedHels[currentHel] << endl << "------------------------------------" << endl;
                if(matchedTwice == true)
                    cout << "This traceHelix was matched twice" << endl;
            }
            else
            {
                cout << endl << "trueSheet #" << currentHel+1 << " --> traceSheet" << matchedHels[currentHel]-(currHel+1)-1 << endl << "------------------------------------" << endl;
                if(matchedTwice == true)
                    cout << "This traceStrand was matched twice" << endl;
            }
        }
    }
    if(one == true)
    {
        if(matchedHels[currentHel] > 100  || (noSheet == true && currentHel >= helixOffset))
        {
            if(currentHel < helixOffset)
                cout << endl << "trueHelix #" << currentHel+1 << endl << "------------------------------------" << endl;
            else
                cout << endl << "trueSheet #" << currentHel+1 << endl << "------------------------------------" << endl;

        }
        else
        {
            if(currentHel < helixOffset)
            {
                cout << endl << "trueHelix #" << currentHel+1 << " --> traceHelix" << matchedHels[currentHel]+1 << endl << "------------------------------------" << endl;
                if(matchedTwice == true)
                    cout << "This traceHelix was matched twice" << endl;
            }
            else
            {
                cout << endl << "trueSheet #" << currentHel+1 << " --> traceSheet" << matchedHels[currentHel]-(currHel+1)+1 << endl << "------------------------------------" << endl;
                if(matchedTwice == true)
                    cout << "This traceStrand was matched twice" << endl;
            }

        }

    }

    if (matchedHels[currentHel] != 99999 && ((noSheet == false)||(noSheet == true && currentHel < helixOffset)))   // If there is a match
        {
            if(shortHelix[currentHel] == true)
            {
                cout << "Distance From Center:\t\t\t " << shortHels[currentHel] << endl;
                if(outExist == true)
                {
                    outFile << currentHel+1 << ", " << matchedHels[currentHel]+1 << ", ";
                    outFile << "Distance from Center = " << shortHels[currentHel] << endl;
                }
            }
            else
            {
               helTrueArray[currentHel].alignDirectionHels(&helTraceArray[matchedHels[currentHel]]);
                displaced = helTrueArray[currentHel].findDisplacement(&helTraceArray[matchedHels[currentHel]], stepSize);
                if(outExist == true)
                {
                    outFile << currentHel+1 << ", " << matchedHels[currentHel]+1 << ", ";
                    outFile << displaced.lengthAxis1 << ", " << displaced.lengthAxis2 << ", " << displaced.twoWayDistance << ", " << displaced.crossDisplace << ", " << displaced.longDisplace << ", " << displaced.lengthProportion << ", ";
                    if(matchedHels[currentHel] < numtraceHels)
                    {
                        tot = specificityBoundaries(pdb, &helTrueArray[currentHel], true, currentHel);
                        if(tot-(pdb.hlces.at(currentHel).endIndx-pdb.hlces.at(currentHel).startIndx) == 0)
                            spec = 1;
                        else
                        {
                            spec = 1.0-(double)((double)specArray[matchedHels[currentHel]]/(double)(tot-(pdb.hlces.at(currentHel).endIndx-pdb.hlces.at(currentHel).startIndx)));
                            outFile << "1-" << specArray[matchedHels[currentHel]] << "/" << tot << "-" << pdb.hlces.at(currentHel).endIndx << "-" << pdb.hlces.at(currentHel).startIndx << ", ";
                        }
                        outFile << spec << ", ";
                    }
                    else if(matchedHels[currentHel] < number_of_hel)
                    {
                        tot = specificityBoundaries(pdb, &helTrueArray[currentHel], false, currentHel);
                        if(tot-(pdb.sheets.at(currentHel-pdb.hlces.size()).endIndx-pdb.sheets.at(currentHel-pdb.hlces.size()).startIndx) == 0)
                            spec = 1;
                        else
                        {
                            spec = 1.0-(double)((double)specArray[matchedHels[currentHel]]/(double)(tot-(pdb.sheets.at(currentHel-pdb.hlces.size()).endIndx-pdb.sheets.at(currentHel-pdb.hlces.size()).startIndx)));
                            outFile << "1-" << specArray[matchedHels[currentHel]] << "/" << tot << "-" << pdb.sheets.at(currentHel-pdb.hlces.size()).endIndx << "-" << pdb.sheets.at(currentHel-pdb.hlces.size()).startIndx << ", ";
                        }
                        outFile << spec << ", ";
                    }
                    if(currentHel >= helixOffset)
                    {
                        sen = PrintSensitivityStrands(&helTraceArray[matchedHels[currentHel]], pdb.sheets[currentHel-helixOffset].startIndx, pdb.sheets[currentHel-helixOffset].endIndx, pdb, sradius, outFile);
                    }
                    else
                    {
                        sen = PrintSensitivity(&helTraceArray[matchedHels[currentHel]], pdb.hlces[currentHel].startIndx, pdb.hlces[currentHel].endIndx, pdb, sradius, outFile);
                    }
                    outFile << (2*((sen*spec)/(sen+spec)));
                    outFile << endl;
                }
                cout << "Lateral Discrepancy:\t\t\t " << displaced.crossDisplace << "\nLongitudinal Discrepancy:\t\t      " << displaced.longDisplace << endl;

            }
        }
    else {
        cout << "No trace corresponds enough to true helix/sheet " << currentHel+1 << endl << endl;
        double lengthAxis1 = helTrueArray[currentHel].getArcLength(helTrueArray[currentHel].axisPoints.begin(), helTrueArray[currentHel].axisPoints.end()-1);
        if(outExist == true)
            outFile << currentHel+1 << ", " << "N/A" << ", " << lengthAxis1 << ", " << "N/A" << ", " << "N/A" << ", " << "N/A" << ", " << "N/A" << endl;
    }

    //for purposely split helices
    if(it2 < splitNum.size() && currentHel == splitNum[it2])
    {
        if(one == false)
        {
            if(matchedHels[currentHel] > 100 || shortHelix[currentHel] == true)
            {
                cout << endl << "trueHelix #" << currentHel+1 << endl << "------------------------------------" << endl;
            }
            else
            {
                cout << endl << "Helix split because of sharp angle\n\n" << "trueHelix #" << currentHel+1 << " --> traceHelix" << splitHels[it2] << endl << "------------------------------------" << endl;
            }
        }
        if(one == true)
        {
            if(matchedHels[currentHel] > 100 || shortHelix[currentHel] == true)
            {
                cout << endl << "trueHelix #" << currentHel+1 << endl << "------------------------------------" << endl;

            }
            else
            {
                cout << endl << "Helix split because of sharp angle\n\n" << "trueHelix #" << currentHel+1 << " --> traceHelix" << splitHels[it2]+1 << endl << "------------------------------------" << endl;
            }

        }

        cout << "Lateral Discrepancy:\t\t\t" << splitLateral[it2] << "\nLongitudinal Discrepancy:\t\t      " << splitLongitudinal[it2] << endl;
        it2++;
    }

    //for cases where trace helix is split and two trace helices need to be matched to a true helix
    if(matchedAlternate[currentHel] != 99999)
    {
        helTrueArray[currentHel].alignDirectionHels(&helTraceArray[matchedAlternate[currentHel]]);
        displaced = helTrueArray[currentHel].findDisplacement(&helTraceArray[matchedAlternate[currentHel]], stepSize);

        helTrueArray[currentHel].alignDirectionHels(&helTraceArray[matchedHels[currentHel]]);
        Dis3 = helTrueArray[currentHel].findDisplacement(&helTraceArray[matchedHels[currentHel]], stepSize);
        if(one == false)
        {
            if(currentHel < helixOffset)
            {
                cout << endl << "traceHelix split due to flawed density map" << endl;
                cout << "If combined:\n" << "Lateral Discrepancy:\t\t\t " << (displaced.crossDisplace + Dis3.crossDisplace)/2 << "\nLongitudinal Discrepancy:\t\t      " << (displaced.longDisplace+Dis3.longDisplace-displaced.lengthAxis1)/2 << "\n\n";

                cout << "trueHelix #" << currentHel+1 << " --> traceHelix" << matchedAlternate[currentHel] << endl << "------------------------------------" << endl;
                cout << "Lateral Discrepancy:\t\t\t " << displaced.crossDisplace << "\nLongitudinal Discrepancy:\t\t      " << displaced.longDisplace << endl;
            }
            else
            {
                cout << endl << "trueSheet #" << currentHel+1 << " --> traceSheet" << matchedAlternate[currentHel]-offset-1 << endl << "------------------------------------" << endl;
                cout << "Lateral Discrepancy:\t\t\t " << displaced.crossDisplace << "\nLongitudinal Discrepancy:\t\t      " << displaced.longDisplace << endl;
            }
        }
        if(one == true)
        {
            if(currentHel < helixOffset)
            {
                cout << endl << "traceHelix split due to flawed density map" << endl;
                cout << "If combined:\n" << "Lateral Discrepancy:\t\t\t " << (displaced.crossDisplace + Dis3.crossDisplace)/2 << "\nLongitudinal Discrepancy:\t\t      " << (displaced.longDisplace+Dis3.longDisplace-displaced.lengthAxis1)/2 << "\n\n";

                cout << "trueHelix #" << currentHel+1 << " --> traceHelix" << matchedAlternate[currentHel]+1 << endl << "------------------------------------" << endl;
                cout << "Lateral Discrepancy:\t\t\t " << displaced.crossDisplace << "\nLongitudinal Discrepancy:\t\t      " << displaced.longDisplace << endl;
            }
            else
            {
                cout << endl << "trueSheet #" << currentHel+1 << " --> traceSheet" << matchedAlternate[currentHel]-offset << endl << "------------------------------------" << endl;
                cout << "Lateral Discrepancy:\t\t\t " << displaced.crossDisplace << "\nLongitudinal Discrepancy:\t\t      " << displaced.longDisplace << endl;
            }
        }
    }
}
//cout << "# of nonHLX AA: " << nonHLXAA << endl;
if(outExist == true)
    outFile.close();

return 0;
}

int regionQuery(int tempGroup[], int pointPlace, Map mrc, int groupNum, bool visited[], double threshold)
{
    int tempI = 0;
    int tempJ = 0;
    int tempK = 0;
    int offset = 0;
    tempK = pointPlace/(mrc.numCols()*mrc.numRows());
    tempJ = (pointPlace-(tempK*mrc.numCols()*mrc.numRows()))/mrc.numRows();
    tempI = (pointPlace-(tempK*mrc.numCols()*mrc.numRows())-(tempJ*mrc.numRows()));

    //center point
    tempGroup[pointPlace] = groupNum;
    ///i
    //first point
    if(mrc.cube[tempI][tempJ+1][tempK] > threshold && visited[((tempK)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI))] == false)
    {
        visited[((tempK)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI))] = true;
        tempGroup[((tempK)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI))] = groupNum;
        returnedNum++;
        regionQuery(tempGroup, ((tempK)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI)), mrc, groupNum, visited, threshold);
    }
    //second point
    if(mrc.cube[tempI][tempJ+1][tempK-1] > threshold && visited[((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI))] == false)
    {
        visited[((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI))] = true;
        tempGroup[((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI))] = groupNum;
        returnedNum++;
        regionQuery(tempGroup, ((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI)), mrc, groupNum, visited, threshold);
    }
    //third point
    if(mrc.cube[tempI][tempJ+1][tempK+1] > threshold && visited[((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI))] == false)
    {
        visited[((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI))] = true;
        tempGroup[((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI))] = groupNum;
        returnedNum++;
        regionQuery(tempGroup, ((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI)), mrc, groupNum, visited, threshold);
    }
    //fourth point
    if(mrc.cube[tempI][tempJ-1][tempK] > threshold && visited[((tempK)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI))] == false)
    {
        visited[((tempK)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI))] = true;
        tempGroup[((tempK)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI))] = groupNum;
        returnedNum++;
        regionQuery(tempGroup, ((tempK)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI)), mrc, groupNum, visited, threshold);
    }
    //fifth point
    if(mrc.cube[tempI][tempJ-1][tempK-1] > threshold && visited[((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI))] == false)
    {
        visited[((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI))] = true;
        tempGroup[((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI))] = groupNum;
        returnedNum++;
        regionQuery(tempGroup, ((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI)), mrc, groupNum, visited, threshold);
    }
    //sixth point
    if(mrc.cube[tempI][tempJ-1][tempK+1] > threshold && visited[((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI))] == false)
    {
        visited[((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI))] = true;
        tempGroup[((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI))] = groupNum;
        returnedNum++;
        regionQuery(tempGroup, ((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI)), mrc, groupNum, visited, threshold);
    }
    //seventh point
    if(mrc.cube[tempI][tempJ][tempK-1] > threshold && visited[((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI))] == false)
    {
        visited[((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI))] = true;
        tempGroup[((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI))] = groupNum;
        returnedNum++;
        regionQuery(tempGroup, ((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI)), mrc, groupNum, visited, threshold);
    }
    //eighth point
    if(mrc.cube[tempI][tempJ][tempK+1] > threshold && visited[((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI))] == false)
    {
        visited[((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI))] = true;
        tempGroup[((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI))] = groupNum;
        returnedNum++;
        regionQuery(tempGroup, ((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI)), mrc, groupNum, visited, threshold);
    }
    ///i+1
    //first point
    if(mrc.cube[tempI+1][tempJ+1][tempK] > threshold && visited[((tempK)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI+1))] == false)
    {
        visited[((tempK)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI+1))] = true;
        tempGroup[((tempK)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI+1))] = groupNum;
        returnedNum++;
        regionQuery(tempGroup, ((tempK)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI+1)), mrc, groupNum, visited, threshold);
    }
    //second point
    if(mrc.cube[tempI+1][tempJ+1][tempK-1] > threshold && visited[((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI+1))] == false)
    {
        visited[((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI+1))] = true;
        tempGroup[((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI+1))] = groupNum;
        returnedNum++;
        regionQuery(tempGroup, ((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI+1)), mrc, groupNum, visited, threshold);
    }
    //third point
    if(mrc.cube[tempI+1][tempJ+1][tempK+1] > threshold && visited[((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI+1))] == false)
    {
        visited[((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI+1))] = true;
        tempGroup[((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI+1))] = groupNum;
        returnedNum++;
        regionQuery(tempGroup, ((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI+1)), mrc, groupNum, visited, threshold);
    }
    //fourth point
    if(mrc.cube[tempI+1][tempJ-1][tempK] > threshold && visited[((tempK)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI+1))] == false)
    {
        visited[((tempK)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI+1))] = true;
        tempGroup[((tempK)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI+1))] = groupNum;
        returnedNum++;
        regionQuery(tempGroup, ((tempK)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI+1)), mrc, groupNum, visited, threshold);
    }
    //fifth point
    if(mrc.cube[tempI+1][tempJ-1][tempK-1] > threshold && visited[((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI+1))] == false)
    {
        visited[((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI+1))] = true;
        tempGroup[((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI+1))] = groupNum;
        returnedNum++;
        regionQuery(tempGroup, ((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI+1)), mrc, groupNum, visited, threshold);
    }
    //sixth point
    if(mrc.cube[tempI+1][tempJ-1][tempK+1] > threshold && visited[((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI+1))] == false)
    {
        visited[((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI+1))] = true;
        tempGroup[((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI+1))] = groupNum;
        returnedNum++;
        regionQuery(tempGroup, ((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI+1)), mrc, groupNum, visited, threshold);
    }
    //seventh point
    if(mrc.cube[tempI+1][tempJ][tempK-1] > threshold && visited[((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI+1))] == false)
    {
        visited[((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI+1))] = true;
        tempGroup[((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI+1))] = groupNum;
        returnedNum++;
        regionQuery(tempGroup, ((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI+1)), mrc, groupNum, visited, threshold);
    }
    //eighth point
    if(mrc.cube[tempI+1][tempJ][tempK+1] > threshold && visited[((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI+1))] == false)
    {
        visited[((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI+1))] = true;
        tempGroup[((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI+1))] = groupNum;
        returnedNum++;
        regionQuery(tempGroup, ((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI+1)), mrc, groupNum, visited, threshold);
    }
    //ninth point
    if(mrc.cube[tempI+1][tempJ][tempK] > threshold && visited[((tempK)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI+1))] == false)
    {
        visited[((tempK)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI+1))] = true;
        tempGroup[((tempK)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI+1))] = groupNum;
        returnedNum++;
        regionQuery(tempGroup, ((tempK)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI+1)), mrc, groupNum, visited, threshold);
    }
    ///i-1
    //first point
    if(mrc.cube[tempI-1][tempJ+1][tempK] > threshold && visited[((tempK)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI-1))] == false)
    {
        visited[((tempK)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI-1))] = true;
        tempGroup[((tempK)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI-1))] = groupNum;
        returnedNum++;
        regionQuery(tempGroup, ((tempK)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI-1)), mrc, groupNum, visited, threshold);
    }
    //second point
    if(mrc.cube[tempI-1][tempJ+1][tempK-1] > threshold && visited[((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI-1))] == false)
    {
        visited[((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI-1))] = true;
        tempGroup[((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI-1))] = groupNum;
        returnedNum++;
        regionQuery(tempGroup, ((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI-1)), mrc, groupNum, visited, threshold);
    }
    //third point
    if(mrc.cube[tempI-1][tempJ+1][tempK+1] > threshold && visited[((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI-1))] == false)
    {
        visited[((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI-1))] = true;
        tempGroup[((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI-1))] = groupNum;
        returnedNum++;
        regionQuery(tempGroup, ((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI-1)), mrc, groupNum, visited, threshold);
    }
    //fourth point
    if(mrc.cube[tempI-1][tempJ-1][tempK] > threshold && visited[((tempK)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI-1))] == false)
    {
        visited[((tempK)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI-1))] = true;
        tempGroup[((tempK)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI-1))] = groupNum;
        returnedNum++;
        regionQuery(tempGroup, ((tempK)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI-1)), mrc, groupNum, visited, threshold);
    }
    //fifth point
    if(mrc.cube[tempI-1][tempJ-1][tempK-1] > threshold && visited[((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI-1))] == false)
    {
        visited[((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI-1))] = true;
        tempGroup[((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI-1))] = groupNum;
        returnedNum++;
        regionQuery(tempGroup, ((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI-1)), mrc, groupNum, visited, threshold);
    }
    //sixth point
    if(mrc.cube[tempI-1][tempJ-1][tempK+1] > threshold && visited[((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI-1))] == false)
    {
        visited[((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI-1))] = true;
        tempGroup[((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI-1))] = groupNum;
        returnedNum++;
        regionQuery(tempGroup, ((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI-1)), mrc, groupNum, visited, threshold);
    }
    //seventh point
    if(mrc.cube[tempI-1][tempJ][tempK-1] > threshold && visited[((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI-1))] == false)
    {
        visited[((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI-1))] = true;
        tempGroup[((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI-1))] = groupNum;
        returnedNum++;
        regionQuery(tempGroup, ((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI-1)), mrc, groupNum, visited, threshold);
    }
    //eighth point
    if(mrc.cube[tempI-1][tempJ][tempK+1] > threshold && visited[((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI-1))] == false)
    {
        visited[((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI-1))] = true;
        tempGroup[((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI-1))] = groupNum;
        returnedNum++;
        regionQuery(tempGroup, ((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI-1)), mrc, groupNum, visited, threshold);
    }
    //ninth point
    if(mrc.cube[tempI-1][tempJ][tempK] > threshold && visited[((tempK)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI-1))] == false)
    {
        visited[((tempK)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI-1))] = true;
        tempGroup[((tempK)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI-1))] = groupNum;
        returnedNum++;
        regionQuery(tempGroup, ((tempK)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI-1)), mrc, groupNum, visited, threshold);
    }


    //for(int i = 0; i < 50; i++)
    //    cout << tempGroup[i] << " ";
    //cout << "\n\n";

    return returnedNum;
}

int regionQueryIt(int tempGroup[], int pointPlace, Map mrc, int groupNum, bool visited[], int groupCurrent[], int &num, bool notVisited[], double threshold)
{
    int returnedNum = 0;
    int tempI = 0;
    int tempJ = 0;
    int tempK = 0;
    int offset = 0;
    tempK = pointPlace/(mrc.numCols()*mrc.numRows());
    tempJ = (pointPlace-(tempK*mrc.numCols()*mrc.numRows()))/mrc.numRows();
    tempI = (pointPlace-(tempK*mrc.numCols()*mrc.numRows())-(tempJ*mrc.numRows()));
    int tempEq = 0;

    //center point
    //tempGroup[pointPlace] = groupNum;
    ///i
    //first point
    if(tempI > 0 && tempJ+1 > 0 && tempK > 0 && tempI < mrc.numRows() && tempJ+1 < mrc.numCols() && tempK < mrc.numSlcs())
    {
        tempEq = ((tempK)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI));
        if(notVisited[tempEq] == false && mrc.cube[tempI][tempJ+1][tempK] > threshold && visited[tempEq] == false )
        {
            notVisited[tempEq] = true;
            //tempGroup[((tempK)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI))] = groupNum;
            returnedNum++;
            groupCurrent[num] = tempEq;
            num++;
        }
    }
    //second point
    if(tempI > 0 && tempJ+1 > 0 && tempK-1 > 0 && tempI < mrc.numRows() && tempJ+1 < mrc.numCols() && tempK-1 < mrc.numSlcs())
    {
        tempEq = ((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI));
        if( notVisited[tempEq] == false && mrc.cube[tempI][tempJ+1][tempK-1] > threshold && visited[tempEq] == false)
        {
            notVisited[tempEq] = true;
            //tempGroup[((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI))] = groupNum;
            returnedNum++;
            groupCurrent[num] = tempEq;
            num++;
        }
    }
    //third point
    if(tempI > 0 && tempJ+1 > 0 && tempK+1 > 0 && tempI < mrc.numRows() && tempJ+1 < mrc.numCols() && tempK+1 < mrc.numSlcs())
    {
        tempEq = ((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI));
        if(notVisited[tempEq] == false && mrc.cube[tempI][tempJ+1][tempK+1] > threshold && visited[tempEq] == false)
        {
            notVisited[tempEq] = true;
            //tempGroup[((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI))] = groupNum;
            returnedNum++;
            groupCurrent[num] = tempEq;
            num++;
        }
    }
    //fourth point
    if(tempI > 0 && tempJ-1 > 0 && tempK > 0 && tempI < mrc.numRows() && tempJ-1 < mrc.numCols() && tempK < mrc.numSlcs())
    {
        tempEq = ((tempK)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI));
        if(notVisited[tempEq] == false && mrc.cube[tempI][tempJ-1][tempK] > threshold && visited[tempEq] == false)
        {
            notVisited[tempEq] = true;
            //tempGroup[((tempK)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI))] = groupNum;
            returnedNum++;
            groupCurrent[num] = tempEq;
            num++;
        }
    }
    //fifth point
    if(tempI > 0 && tempJ-1 > 0 && tempK-1 > 0 && tempI < mrc.numRows() && tempJ-1 < mrc.numCols() && tempK-1 < mrc.numSlcs())
    {
        tempEq = ((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI));
        if(notVisited[tempEq] == false && mrc.cube[tempI][tempJ-1][tempK-1] > threshold && visited[tempEq] == false)
        {
            notVisited[tempEq] = true;
            //tempGroup[((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI))] = groupNum;
            returnedNum++;
            groupCurrent[num] = tempEq;
            num++;
        }
    }
    //sixth point
    if(tempI > 0 && tempJ-1 > 0 && tempK+1 > 0 && tempI < mrc.numRows() && tempJ-1 < mrc.numCols() && tempK+1 < mrc.numSlcs())
    {
        tempEq = ((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI));
        if(notVisited[tempEq] == false && mrc.cube[tempI][tempJ-1][tempK+1] > threshold && visited[tempEq] == false)
        {
            notVisited[tempEq] = true;
            //tempGroup[((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI))] = groupNum;
            returnedNum++;
            groupCurrent[num] = tempEq;
            num++;
        }
    }
    //seventh point
    if(tempI > 0 && tempJ > 0 && tempK-1 > 0 && tempI < mrc.numRows() && tempJ < mrc.numCols() && tempK-1 < mrc.numSlcs())
    {
        tempEq = ((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI));
        if(notVisited[tempEq] == false && mrc.cube[tempI][tempJ][tempK-1] >threshold && visited[tempEq] == false)
        {
            notVisited[tempEq] = true;
            //tempGroup[((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI))] = groupNum;
            returnedNum++;
            groupCurrent[num] = tempEq;
            num++;
        }
    }
    //eighth point
    if(tempI > 0 && tempJ > 0 && tempK+1 > 0 && tempI < mrc.numRows() && tempJ < mrc.numCols() && tempK+1 < mrc.numSlcs())
    {
        tempEq = ((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI));
        if(notVisited[tempEq] == false && mrc.cube[tempI][tempJ][tempK+1] > threshold && visited[tempEq] == false)
        {
            notVisited[tempEq] = true;
            //tempGroup[((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI))] = groupNum;
            returnedNum++;
            groupCurrent[num] = tempEq;
            num++;
        }
    }
    ///i+1
    //first point
    if(tempI+1 > 0 && tempJ+1 > 0 && tempK > 0 && tempI+1 < mrc.numRows() && tempJ+1 < mrc.numCols() && tempK < mrc.numSlcs())
    {
        tempEq = ((tempK)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI+1));
        if(notVisited[tempEq] == false && mrc.cube[tempI+1][tempJ+1][tempK] > threshold && visited[tempEq] == false)
        {
            notVisited[tempEq] = true;
            //tempGroup[((tempK)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI+1))] = groupNum;
            returnedNum++;
            groupCurrent[num] = tempEq;
            num++;
        }
    }
    //second point
    if(tempI+1 > 0 && tempJ+1 > 0 && tempK-1 > 0 && tempI+1 < mrc.numRows() && tempJ+1 < mrc.numCols() && tempK-1 < mrc.numSlcs())
    {
        tempEq = ((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI+1));
        if(notVisited[tempEq] == false && mrc.cube[tempI+1][tempJ+1][tempK-1] > threshold && visited[tempEq] == false)
        {
            notVisited[tempEq] = true;
            //tempGroup[((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI+1))] = groupNum;
            returnedNum++;
            groupCurrent[num] = tempEq;
            num++;
        }
    }
    //third point
    if(tempI+1 > 0 && tempJ+1 > 0 && tempK+1 > 0 && tempI+1 < mrc.numRows() && tempJ+1 < mrc.numCols() && tempK+1 < mrc.numSlcs())
    {
        tempEq = ((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI+1));
        if(notVisited[tempEq] == false && mrc.cube[tempI+1][tempJ+1][tempK+1] > threshold && visited[tempEq] == false)
        {
            notVisited[tempEq] = true;
            //tempGroup[((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI+1))] = groupNum;
            returnedNum++;
            groupCurrent[num] = tempEq;
            num++;
        }
    }
    //fourth point
    if(tempI+1 > 0 && tempJ-1 > 0 && tempK > 0 && tempI+1 < mrc.numRows() && tempJ-1 < mrc.numCols() && tempK < mrc.numSlcs())
    {
        tempEq = ((tempK)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI+1));
        if(notVisited[tempEq] == false && mrc.cube[tempI+1][tempJ-1][tempK] > threshold && visited[tempEq] == false )
        {
            notVisited[tempEq] = true;
            //tempGroup[((tempK)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI+1))] = groupNum;
            returnedNum++;
            groupCurrent[num] = tempEq;
            num++;
        }
    }
    //fifth point
    if(tempI+1 > 0 && tempJ-1 > 0 && tempK-1 > 0 && tempI+1 < mrc.numRows() && tempJ-1 < mrc.numCols() && tempK-1 < mrc.numSlcs())
    {
        tempEq = ((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI+1));
        if(notVisited[tempEq] == false && mrc.cube[tempI+1][tempJ-1][tempK-1] > threshold && visited[tempEq] == false)
        {
            notVisited[tempEq] = true;
            //tempGroup[((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI+1))] = groupNum;
            returnedNum++;
            groupCurrent[num] = tempEq;
            num++;
        }
    }
    //sixth point
    if(tempI+1 > 0 && tempJ-1 > 0 && tempK+1 > 0 && tempI+1 < mrc.numRows() && tempJ-1 < mrc.numCols() && tempK+1 < mrc.numSlcs())
    {
        tempEq = ((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI+1));
        if(notVisited[tempEq] == false && mrc.cube[tempI+1][tempJ-1][tempK+1] > threshold && visited[tempEq] == false)
        {
            notVisited[tempEq] = true;
            //tempGroup[((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI+1))] = groupNum;
            returnedNum++;
            groupCurrent[num] = tempEq;
            num++;
        }
    }
    //seventh point
    if(tempI+1 > 0 && tempJ > 0 && tempK-1 > 0 && tempI+1 < mrc.numRows() && tempJ < mrc.numCols() && tempK-1 < mrc.numSlcs())
    {
        tempEq = ((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI+1));
        if(notVisited[tempEq] == false && mrc.cube[tempI+1][tempJ][tempK-1] > threshold && visited[tempEq] == false)
        {
            notVisited[tempEq] = true;
            //tempGroup[((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI+1))] = groupNum;
            returnedNum++;
            groupCurrent[num] = tempEq;
            num++;
        }
    }
    //eighth point
    if(tempI+1 > 0 && tempJ > 0 && tempK+1 > 0 && tempI+1 < mrc.numRows() && tempJ < mrc.numCols() && tempK+1 < mrc.numSlcs())
    {
        tempEq = ((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI+1));
        if(notVisited[tempEq] == false && mrc.cube[tempI+1][tempJ][tempK+1] > threshold && visited[tempEq] == false)
        {
            notVisited[tempEq] = true;
            //tempGroup[((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI+1))] = groupNum;
            returnedNum++;
            groupCurrent[num] = tempEq;
            num++;
        }
    }
    //ninth point
    if(tempI+1 > 0 && tempJ > 0 && tempK > 0 && tempI+1 < mrc.numRows() && tempJ < mrc.numCols() && tempK < mrc.numSlcs())
    {
        tempEq = ((tempK)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI+1));
        if(notVisited[tempEq] == false && mrc.cube[tempI+1][tempJ][tempK] > threshold && visited[tempEq] == false)
        {
            notVisited[tempEq] = true;
            //tempGroup[((tempK)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI+1))] = groupNum;
            returnedNum++;
            groupCurrent[num] = tempEq;
            num++;
        }
    }
    ///i-1
    //first point
    if(tempI-1 > 0 && tempJ+1 > 0 && tempK > 0 && tempI-1 < mrc.numRows() && tempJ+1 < mrc.numCols() && tempK < mrc.numSlcs())
    {
        tempEq = ((tempK)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI-1));
        if(notVisited[tempEq] == false && mrc.cube[tempI-1][tempJ+1][tempK] > threshold && visited[tempEq] == false)
        {
            notVisited[tempEq] = true;
            //tempGroup[((tempK)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI-1))] = groupNum;
            returnedNum++;
            groupCurrent[num] = tempEq;
            num++;
        }
    }
    //second point
    if(tempI-1 > 0 && tempJ+1 > 0 && tempK-1 > 0 && tempI-1 < mrc.numRows() && tempJ+1 < mrc.numCols() && tempK-1 < mrc.numSlcs())
    {
        tempEq = ((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI-1));
        if(notVisited[tempEq] == false && mrc.cube[tempI-1][tempJ+1][tempK-1] > threshold && visited[tempEq] == false)
        {
            notVisited[tempEq] = true;
            //tempGroup[((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI-1))] = groupNum;
            returnedNum++;
            groupCurrent[num] = tempEq;
            num++;
        }
    }
    //third point
    if(tempI-1 > 0 && tempJ+1 > 0 && tempK+1 > 0 && tempI-1 < mrc.numRows() && tempJ+1 < mrc.numCols() && tempK+1 < mrc.numSlcs())
    {
        tempEq = ((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI-1));
        if(notVisited[tempEq] == false && mrc.cube[tempI-1][tempJ+1][tempK+1] > threshold && visited[tempEq] == false)
        {
            notVisited[tempEq] = true;
            //tempGroup[((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ+1)*mrc.numRows()+(tempI-1))] = groupNum;
            returnedNum++;
            groupCurrent[num] = tempEq;
            num++;
        }
    }
    //fourth point
    if(tempI-1 > 0 && tempJ-1 > 0 && tempK > 0 && tempI-1 < mrc.numRows() && tempJ-1 < mrc.numCols() && tempK < mrc.numSlcs())
    {
        tempEq = ((tempK)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI-1));
        if(notVisited[tempEq] == false && mrc.cube[tempI-1][tempJ-1][tempK] > threshold && visited[tempEq] == false)
        {
            notVisited[tempEq] = true;
            //tempGroup[((tempK)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI-1))] = groupNum;
            returnedNum++;
            groupCurrent[num] = tempEq;
            num++;
        }
    }
    //fifth point
    if(tempI-1 > 0 && tempJ-1 > 0 && tempK-1 > 0 && tempI-1 < mrc.numRows() && tempJ-1 < mrc.numCols() && tempK-1 < mrc.numSlcs())
    {
        tempEq = ((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI-1));
        if(notVisited[tempEq] == false && mrc.cube[tempI-1][tempJ-1][tempK-1] > threshold && visited[tempEq] == false)
        {
            notVisited[tempEq] = true;
            //tempGroup[((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI-1))] = groupNum;
            returnedNum++;
            groupCurrent[num] = tempEq;
            num++;
        }
    }
    //sixth point
    if(tempI-1 > 0 && tempJ-1 > 0 && tempK+1 > 0 && tempI-1 < mrc.numRows() && tempJ-1 < mrc.numCols() && tempK+1 < mrc.numSlcs())
    {
        tempEq = ((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI-1));
        if(notVisited[tempEq] == false && mrc.cube[tempI-1][tempJ-1][tempK+1] > threshold && visited[tempEq] == false)
        {
            notVisited[tempEq] = true;
            //tempGroup[((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ-1)*mrc.numRows()+(tempI-1))] = groupNum;
            returnedNum++;
            groupCurrent[num] = tempEq;
            num++;
        }
    }
    //seventh point
    if(tempI-1 > 0 && tempJ > 0 && tempK-1 > 0 && tempI-1 < mrc.numRows() && tempJ < mrc.numCols() && tempK-1 < mrc.numSlcs())
    {
        tempEq = ((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI-1));
        if(notVisited[tempEq] == false && mrc.cube[tempI-1][tempJ][tempK-1] > threshold && visited[tempEq] == false)
        {
            notVisited[tempEq] = true;
            //tempGroup[((tempK-1)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI-1))] = groupNum;
            returnedNum++;
            groupCurrent[num] = tempEq;
            num++;
        }
    }
    //eighth point
    if(tempI-1 > 0 && tempJ > 0 && tempK+1 > 0 && tempI-1 < mrc.numRows() && tempJ < mrc.numCols() && tempK+1 < mrc.numSlcs())
    {
        tempEq = ((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI-1));
        if(notVisited[tempEq] == false && mrc.cube[tempI-1][tempJ][tempK+1] > threshold && visited[tempEq] == false)
        {
            notVisited[tempEq] = true;
            //tempGroup[((tempK+1)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI-1))] = groupNum;
            returnedNum++;
            groupCurrent[num] = tempEq;
            num++;
        }
    }
    //ninth point
    if(tempI-1 > 0 && tempJ > 0 && tempK > 0 && tempI-1 < mrc.numRows() && tempJ < mrc.numCols() && tempK < mrc.numSlcs())
    {
        tempEq = ((tempK)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI-1));
        if(notVisited[tempEq] == false && mrc.cube[tempI-1][tempJ][tempK] > threshold && visited[tempEq] == false)
        {
            notVisited[tempEq] = true;
            //tempGroup[((tempK)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI-1))] = groupNum;
            returnedNum++;
            groupCurrent[num] = tempEq;
            num++;
        }
    }
    //cout << num << " ";

    //for(int i = 0; i < 50; i++)
    //    cout << tempGroup[i] << " ";
    //cout << "\n\n";
    //cout << returnedNum << " ";
    return returnedNum;
}

int expandGroup(int tempGroup[], bool visited[], Map mrc, int groupCurrent[], int num, int groupNumber, bool notVisited[], double threshold)
{
    int tempI = 0;
    int tempJ = 0;
    int tempK = 0;
    int localCopy = num;
    for(int i = 0; i < localCopy; i++)
    {
        if(groupCurrent[i] != -1)
        {
            tempK = groupCurrent[i]/(mrc.numCols()*mrc.numRows());
            tempJ = (groupCurrent[i]-(tempK*mrc.numCols()*mrc.numRows()))/mrc.numRows();
            tempI = (groupCurrent[i]-(tempK*mrc.numCols()*mrc.numRows())-(tempJ*mrc.numRows()));
            if(visited[((tempK)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI))] == false)
            {
                visited[((tempK)*mrc.numCols()*mrc.numRows()+(tempJ)*mrc.numRows()+(tempI))] == true;
                localCopy += regionQueryIt(tempGroup,(tempK*mrc.numCols()*mrc.numRows()+tempJ*mrc.numRows()+tempI), mrc, groupNumber, visited, groupCurrent, num, notVisited, threshold);
                returnedNum = 0;
            }
        }
        else
        {
            cout << "WWWWWWWWWRRRRRRRRRRONG";
            i = localCopy;
        }


    }
    //cout << localCopy << ": ";
    return localCopy;
}

void linearFit(Map mrc, int tempGroup[], int groupNumber, int total[], string path, double stepSize, bool one, int &currHel, Protein pdb, int helixOffset, int numSplit, vector<int> groupToSplitArr, bool acute)
{
    int counter = 0;
    double meanX = 0;
    double meanY = 0;
    double meanZ = 0;
    int tempK = 0;
    int tempJ = 0;
    int tempI = 0;
    double minX = 9999999;
    double maxX = 0;
    double minY = 9999999;
    double maxY = 0;
    double minZ = 9999999;
    double maxZ = 0;
    int it = 0;
    int it2 = 0;
    vector<Coordinate> splitLastPoint;
    //currHel++;
    ///loop over every group
    for(int x = 0; x < groupNumber; x++)
    {
            using namespace Eigen;
            MatrixXf points(total[x],3);
            MatrixXf realPoints(total[x], 3);
            MatrixXf pointx(total[x],1);
            MatrixXf pointy(total[x],1);
            MatrixXf pointz(total[x],1);
            counter = 0;
            meanX = 0;
            meanY = 0;
            meanZ = 0;
            tempK = 0;
            tempJ = 0;
            tempI = 0;
            //put points in 3d matrix and find mean
            for(int i = 0; i < mrc.numSlcs()*mrc.numCols()*mrc.numRows(); i++)
            {
                if(tempGroup[i] == x)
                {

                    tempK = i/(mrc.numCols()*mrc.numRows());
                    tempJ = (i-(tempK*mrc.numCols()*mrc.numRows()))/mrc.numRows();
                    tempI = (i-(tempK*mrc.numCols()*mrc.numRows())-(tempJ*mrc.numRows()));
                    points(counter, 0) = tempI*mrc.apixX+mrc.hdr.xorigin;
                    points(counter, 1) = tempJ*mrc.apixY+mrc.hdr.yorigin;
                    points(counter, 2) = tempK*mrc.apixZ+mrc.hdr.zorigin;
                    realPoints(counter, 0) = tempI*mrc.apixX+mrc.hdr.xorigin;
                    realPoints(counter, 1) = tempJ*mrc.apixY+mrc.hdr.yorigin;
                    realPoints(counter, 2) = tempK*mrc.apixZ+mrc.hdr.zorigin;
                    pointx(counter, 0) = tempI*mrc.apixX+mrc.hdr.xorigin;
                    pointy(counter, 0) = tempJ*mrc.apixY+mrc.hdr.yorigin;
                    pointz(counter, 0) = tempK*mrc.apixZ+mrc.hdr.zorigin;

                    meanX += (tempI*mrc.apixX+mrc.hdr.xorigin);
                    if((tempI*mrc.apixX+mrc.hdr.xorigin) > maxX)
                        maxX = (tempI*mrc.apixX+mrc.hdr.xorigin);
                    if((tempI*mrc.apixX+mrc.hdr.xorigin) < minX)
                        minX = (tempI*mrc.apixX+mrc.hdr.xorigin);

                    meanY += (tempJ*mrc.apixY+mrc.hdr.yorigin);
                    if((tempJ*mrc.apixY+mrc.hdr.yorigin) > maxY)
                        maxY = (tempJ*mrc.apixY+mrc.hdr.yorigin);
                    if((tempJ*mrc.apixY+mrc.hdr.yorigin) < minY)
                        minY = (tempJ*mrc.apixY+mrc.hdr.yorigin);

                    meanZ += (tempK*mrc.apixZ+mrc.hdr.zorigin);
                    if((tempK*mrc.apixZ+mrc.hdr.zorigin) > maxZ)
                        maxZ = (tempK*mrc.apixZ+mrc.hdr.zorigin);
                    if((tempK*mrc.apixZ+mrc.hdr.zorigin) < minZ)
                        minZ = (tempK*mrc.apixZ+mrc.hdr.zorigin);

                    counter++;
                }
            }
            ///if group is large enough
            if(counter > 20)
            {
                meanX /= counter;
                meanY /= counter;
                meanZ /= counter;
                //subtract mean from points
                for(int i = 0; i < counter; i++)
                {
                    points(i, 0) -= meanX;
                    points(i, 1) -= meanY;
                    points(i, 2) -= meanZ;
                }
                //perform SVD
                MatrixXf V(3,3);
                JacobiSVD<MatrixXf> svd(points, ComputeThinU | ComputeThinV);
                V = svd.matrixV();

                //output files
                ofstream outCoordinates101;
                string fileName121 = "";
                stringstream ss1;
                ss1 << x;
                string str2 = ss1.str();
                fileName121 = path + "GroupRealHelix" + str2 +".pdb";
                outCoordinates101.open(fileName121.c_str());
                double result11 = 0;
                double result21 = 0;
                double result31 = 0;
                double minI = -30;
                double maxI = 0;
                double distance2 = 0;
                double minDistance = 999999999999999999999999;
                //cout << V(0,0) << "  " << V(0,1) << "  " << V(0,2) << endl;
                ///print strait line
                for(double i = -15; i < 15.01; i+=.1)
                {
                    minDistance = 999999999999999999999999;
                    result11 = meanX + V(0,0)*i;
                    result21 = meanY + V(1,0)*i;
                    result31 = meanZ + V(2,0)*i;
                    //find min distance to shorten axis
                    for(int z = 0; z < total[x]; z++)
                    {
                        distance2 = sqrt(pow((realPoints(z,0)-result11),2)+pow((realPoints(z,1)-result21),2)+pow((realPoints(z,2)-result31),2));
                        if(distance2 < minDistance)
                            minDistance = distance2;
                    }
                    if(result11 > minX && result11 < maxX && result21 > minY && result21 < maxY && result31 > minZ && result31 < maxZ)
                    {
                        if(minI == -30)
                            minI = i;
                        maxI = i;
                        outCoordinates101 << fixed << setprecision(3) << "ATOM      1 CA GLY A   1      " <<  result11 << "    " << result21 << "    " << result31 << endl;
                    }
                }
                outCoordinates101.close();

                ///compute curve
                MatrixXf t(total[x],3);
                double pointsDistance[counter];
                double distance = 0;
                for(int i = 0; i < counter; i++)
                {
                    pointsDistance[i] = 1.8 * pow(10,308);
                }
                for(int i = 0; i < counter; i++)
                {
                    for(double j = minI; j < (maxI+.01); j+=.1)
                    {
                        result11 = meanX + V(0,0)*j;
                        result21 = meanY + V(1,0)*j;
                        result31 = meanZ + V(2,0)*j;
                        distance = sqrt(pow((realPoints(i,0)-result11),2)+pow((realPoints(i,1)-result21),2)+pow((realPoints(i,2)-result31),2));
                        if(distance < pointsDistance[i])
                        {
                            pointsDistance[i] = distance;
                            t(i,0) = 1;
                            t(i,1) = j;
                            t(i,2) = j*j;
                        }
                    }
                }
                MatrixXf a(3,1);
                MatrixXf b(3,1);
                MatrixXf c(3,1);
                a = t.jacobiSvd(ComputeThinU | ComputeThinV).solve(pointx);
                b = t.jacobiSvd(ComputeThinU | ComputeThinV).solve(pointy);
                c = t.jacobiSvd(ComputeThinU | ComputeThinV).solve(pointz);

                ///outputCurve
                string fileName12 = "";
                stringstream ss;
                ss << x;
                string str = "";
                str = ss.str();
                double result1 = 0;
                double result2 = 0;
                double result3 = 0;
                Coordinate tmppnt;
                vector<Coordinate> helPoints;
                Axis myAxis;
                stringstream sstring;
                sstring << currHel;
                currHel++;
                string current_helix_str = sstring.str();

                for(double i = minI-5; i < maxI+5; i+=1)
                {
                    minDistance = 999999999999999999999999;
                    result1 = a(0,0)+a(1,0)*i+a(2,0)*i*i;
                    result2 = b(0,0)+b(1,0)*i+b(2,0)*i*i;
                    result3 = c(0,0)+c(1,0)*i+c(2,0)*i*i;
                    for(int z = 0; z < total[x]; z++)
                    {
                        distance2 = sqrt(pow((realPoints(z,0)-result1),2)+pow((realPoints(z,1)-result2),2)+pow((realPoints(z,2)-result3),2));
                        if(distance2 < minDistance)
                            minDistance = distance2;
                    }
                    tmppnt.x = result1;
                    tmppnt.y = result2;
                    tmppnt.z = result3;
                    if(minDistance < 1 )
                        helPoints.push_back(tmppnt);
                }
                if(acute == true && x >= groupNumber-numSplit)
                {
                    double dist1 = 0;
                    double dist2 = 0;
                    dist1 = sqrt(pow((splitLastPoint[it2].x-helPoints[0].x),2)+pow((splitLastPoint[it2].y-helPoints[0].y),2)+pow((splitLastPoint[it2].z-helPoints[0].z),2));
                    dist2 = sqrt(pow((splitLastPoint[it2].x-helPoints[helPoints.size()-1].x),2)+pow((splitLastPoint[it2].y-helPoints[helPoints.size()-1].y),2)+pow((splitLastPoint[it2].z-helPoints[helPoints.size()-1].z),2));

                    if(dist1 < dist2)
                    {
                        vector<Coordinate> temp;
                        temp = helPoints;

                        helPoints.clear();
                        for(int r = 5; r < temp.size(); r++)
                        {
                            helPoints.push_back(temp[r]);
                        }
                    }
                    else
                    {
                        helPoints.pop_back();
                        helPoints.pop_back();
                        helPoints.pop_back();
                        helPoints.pop_back();
                        helPoints.pop_back();
                    }

                    it2++;
                }

                if(acute == true && groupToSplitArr[it] == x)
                {
                    splitLastPoint.push_back(helPoints[helPoints.size()-1]);
                    it++;
                }
                myAxis.axisPoints = helPoints;
                myAxis.catmullRom(stepSize);

                fileName12 = path + "/output/traceHelix" + current_helix_str + ".pdb";

                myAxis.printAsPnts2(fileName12);
                helTraceArray.push_back(myAxis);
                helPoints.clear();
            }
    }
}

void outputPoints(Map mrc, string path, double threshold)
{
    ofstream outCoordinates100;
    string fileName120 = "";
    fileName120 = path + "GROUPALLLLLL.pdb";
    outCoordinates100.open(fileName120.c_str());
    double result10 = 0;
    double result20 = 0;
    double result30 = 0;
    for (long k=0; k<mrc.numSlcs(); k++)
        for (long j=0; j<mrc.numCols(); j++)
            for (long i=0; i<mrc.numRows(); i++)
            {
                if(mrc.cube[i][j][k] > threshold)
                {
                    result10 = i*mrc.apixX+mrc.hdr.xorigin;
                    result20 = j*mrc.apixY+mrc.hdr.yorigin;
                    result30 = k*mrc.apixZ+mrc.hdr.zorigin;
                    if(result10 < 10 && result10 > 0)
                    {
                        outCoordinates100 << fixed << setprecision(3) << "ATOM      1 CA GLY A   1        " <<  result10;
                        if(result20 < 10 && result20 > 0)
                        {
                            outCoordinates100 << fixed << setprecision(3) << "    " << result20;
                            if(result30 < 10 && result30 > 0)
                            {
                                outCoordinates100 << fixed << setprecision(3) << "    " << result30 << endl;
                            }
                            else if((result30 > 10 ^ result30 < 0) && (result30 > -10))
                            {
                                outCoordinates100 << fixed << setprecision(3) << "   " << result30 << endl;
                            }
                            else if((result30 > 10 && result30 < 0) || (result30 < -10))
                            {
                                outCoordinates100 << fixed << setprecision(3) << "  " << result30 << endl;
                            }
                        }
                        else if((result20 > 10 ^ result20 < 0) && (result20 > -10))
                        {
                            outCoordinates100 << fixed << setprecision(3) << "   " << result20;
                            if(result30 < 10 && result30 > 0)
                            {
                                outCoordinates100 << fixed << setprecision(3) << "    " << result30 << endl;
                            }
                            else if((result30 > 10 ^ result30 < 0) && (result30 > -10))
                            {
                                outCoordinates100 << fixed << setprecision(3) << "   " << result30 << endl;
                            }
                            else if((result30 > 10 && result30 < 0) || (result30 < -10))
                            {
                                outCoordinates100 << fixed << setprecision(3) << "  " << result30 << endl;
                            }
                        }
                        else if((result20 > 10 && result20 < 0) || (result20 < -10))
                        {
                            outCoordinates100 << fixed << setprecision(3) << "  " << result20;
                            if(result30 < 10 && result30 > 0)
                            {
                                outCoordinates100 << fixed << setprecision(3) << "    " << result30 << endl;
                            }
                            else if((result30 > 10 ^ result30 < 0) && (result30 > -10))
                            {
                                outCoordinates100 << fixed << setprecision(3) << "   " << result30 << endl;
                            }
                            else if((result30 > 10 && result30 < 0) || (result30 < -10))
                            {
                                outCoordinates100 << fixed << setprecision(3) << "  " << result30 << endl;
                            }
                        }
                    }
                    else if((result10 > 10 ^ result10 < 0) && (result10 > -10))
                    {
                        outCoordinates100 << fixed << setprecision(3) << "ATOM      1 CA GLY A   1       " <<  result10;
                        if(result20 < 10 && result20 > 0)
                        {
                            outCoordinates100 << fixed << setprecision(3) << "    " << result20;
                            if(result30 < 10 && result30 > 0)
                            {
                                outCoordinates100 << fixed << setprecision(3) << "    " << result30 << endl;
                            }
                            else if((result30 > 10 ^ result30 < 0) && (result30 > -10))
                            {
                                outCoordinates100 << fixed << setprecision(3) << "   " << result30 << endl;
                            }
                            else if((result30 > 10 && result30 < 0) || (result30 < -10))
                            {
                                outCoordinates100 << fixed << setprecision(3) << "  " << result30 << endl;
                            }
                        }
                        else if((result20 > 10 ^ result20 < 0) && (result20 > -10))
                        {
                            outCoordinates100 << fixed << setprecision(3) << "   " << result20;
                            if(result30 < 10 && result30 > 0)
                            {
                                outCoordinates100 << fixed << setprecision(3) << "    " << result30 << endl;
                            }
                            else if((result30 > 10 ^ result30 < 0) && (result30 > -10))
                            {
                                outCoordinates100 << fixed << setprecision(3) << "   " << result30 << endl;
                            }
                            else if((result30 > 10 && result30 < 0) || (result30 < -10))
                            {
                                outCoordinates100 << fixed << setprecision(3) << "  " << result30 << endl;
                            }
                        }
                        else if((result20 > 10 && result20 < 0) || (result20 < -10))
                        {
                            outCoordinates100 << fixed << setprecision(3) << "  " << result20;
                            if(result30 < 10 && result30 > 0)
                            {
                                outCoordinates100 << fixed << setprecision(3) << "    " << result30 << endl;
                            }
                            else if((result30 > 10 ^ result30 < 0) && (result30 > -10))
                            {
                                outCoordinates100 << fixed << setprecision(3) << "   " << result30 << endl;
                            }
                            else if((result30 > 10 && result30 < 0) || (result30 < -10))
                            {
                                outCoordinates100 << fixed << setprecision(3) << "  " << result30 << endl;
                            }
                        }
                    }
                    else if((result10 > 10 && result10 < 0) || (result10 < -10))
                    {
                        outCoordinates100 << fixed << setprecision(3) << "ATOM      1 CA GLY A   1      " <<  result10;
                        if(result20 < 10 && result20 > 0)
                        {
                            outCoordinates100 << fixed << setprecision(3) << "    " << result20;
                            if(result30 < 10 && result30 > 0)
                            {
                                outCoordinates100 << fixed << setprecision(3) << "    " << result30 << endl;
                            }
                            else if((result30 > 10 ^ result30 < 0) && (result30 > -10))
                            {
                                outCoordinates100 << fixed << setprecision(3) << "   " << result30 << endl;
                            }
                            else if((result30 > 10 && result30 < 0) || (result30 < -10))
                            {
                                outCoordinates100 << fixed << setprecision(3) << "  " << result30 << endl;
                            }
                        }
                        else if((result20 > 10 ^ result20 < 0) && (result20 > -10))
                        {
                            outCoordinates100 << fixed << setprecision(3) << "   " << result20;
                            if(result30 < 10 && result30 > 0)
                            {
                                outCoordinates100 << fixed << setprecision(3) << "    " << result30 << endl;
                            }
                            else if((result30 > 10 ^ result30 < 0) && (result30 > -10))
                            {
                                outCoordinates100 << fixed << setprecision(3) << "   " << result30 << endl;
                            }
                            else if((result30 > 10 && result30 < 0) || (result30 < -10))
                            {
                                outCoordinates100 << fixed << setprecision(3) << "  " << result30 << endl;
                            }
                        }
                        else if((result20 > 10 && result20 < 0) || (result20 < -10))
                        {
                            outCoordinates100 << fixed << setprecision(3) << "  " << result20;
                            if(result30 < 10 && result30 > 0)
                            {
                                outCoordinates100 << fixed << setprecision(3) << "    " << result30 << endl;
                            }
                            else if((result30 > 10 ^ result30 < 0) && (result30 > -10))
                            {
                                outCoordinates100 << fixed << setprecision(3) << "   " << result30 << endl;
                            }
                            else if((result30 > 10 && result30 < 0) || (result30 < -10))
                            {
                                outCoordinates100 << fixed << setprecision(3) << "  " << result30 << endl;
                            }
                        }
                    }

                }
            }
    outCoordinates100.close();
}
