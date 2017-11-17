#ifndef AXIS_H_INCLUDED
#define AXIS_H_INCLUDED

struct Displacement
{
    double lengthAxis1;
    double lengthAxis2;
    double lengthComparison;

    double twoWayDistance;

    double crossDisplace;
    double longDisplace;
    double lengthProportion;
};

class Axis{
public:
    inline void addPoint(Coordinate c){axisPoints.push_back(c);}
    void reverseOrder();
    void alignDirectionHels(Axis* qptr);
    Displacement findDisplacement(Axis* qptr, double stepSize);
    Displacement findCrossDisplacement(Axis* qptr, double stepSize);
    double getArcLength(vector<Coordinate>::iterator p1, vector<Coordinate>::iterator p2);
    void catmullRom(double stepSize);
    void addTrueEnds(Coordinate pdbEnd1, Coordinate pdbEnd2);
    double twoWayDistance(Axis* qptr);

    bool angle();
    void splitFirst();
    void splitSecond();
    void splitFirstHalf();
    void splitSecondHalf();
    Coordinate firstPoint();
    Coordinate lastPoint();

    void printAsPnts (string outputFileName);
    void printAsPnts2 (string outputFileName);
    //inline void changePoint(){}

// should be private but I'm not being that good.
    vector<Coordinate> axisPoints;



};

void Axis::addTrueEnds(Coordinate pdbEnd1, Coordinate pdbEnd2) // where End1 is the "far" end of the helix, at the end of axisPoints
{
    Vectors axisEnd1;
    Vectors axisEnd2;
    if(axisPoints.size() < 2)
        return;
    else if(axisPoints.size() == 2)
    {
        axisEnd1.set(axisPoints[1], axisPoints[0]);
        axisEnd2.set(axisPoints[0], axisPoints[1]);
    }
    else
    {
        axisEnd1.set(axisPoints[axisPoints.size()-1], axisPoints[axisPoints.size()-3]);
        axisEnd2.set(axisPoints[0], axisPoints[2]);
    }

    Coordinate tmppnt;
    if(getDistance(pdbEnd1, axisPoints[axisPoints.size()-1]) > getDistance(pdbEnd2, axisPoints[axisPoints.size()-1]))
    {
        tmppnt = pdbEnd1;
        pdbEnd1 = pdbEnd2;
        pdbEnd2 = tmppnt;
    }


    double multiplier = 1;
    double lengthNewEnd = 15.0;
    Coordinate newEnd;
    Coordinate tmpEnd;
    Coordinate bestEnd;
    Vectors tmpVect;
    double closestToPlane = 99999999;
    multiplier = lengthNewEnd/axisEnd1.length();
    axisEnd1 = axisEnd1.mul(multiplier);
    newEnd.x = axisPoints[axisPoints.size()-1].x + axisEnd1.getX();
    newEnd.y = axisPoints[axisPoints.size()-1].y + axisEnd1.getY();
    newEnd.z = axisPoints[axisPoints.size()-1].z + axisEnd1.getZ();

    for (double i = 0; i < 5.0; i+=0.1)
    {
        tmpEnd = pointOnLine(axisPoints[axisPoints.size()-1], newEnd, (float)i);
        tmpVect.set(tmpEnd, pdbEnd1);
        if(abs(axisEnd1.dot(tmpVect)) < closestToPlane)
        {
            closestToPlane = abs(axisEnd1.dot(tmpVect));
            bestEnd = tmpEnd;
        }
    }
    axisPoints.push_back(bestEnd);

/// Error is in this section
    closestToPlane = 99999999;
    multiplier = lengthNewEnd/axisEnd2.length();
    axisEnd2 = axisEnd2.mul(multiplier);
    newEnd.x = axisPoints[0].x + axisEnd2.getX();
    newEnd.y = axisPoints[0].y + axisEnd2.getY();
    newEnd.z = axisPoints[0].z + axisEnd2.getZ();

    for (double i = 0; i < 5.0; i+=0.1)
    {
        tmpEnd = pointOnLine(axisPoints[0], newEnd, (float)i);
        tmpVect.set(tmpEnd, pdbEnd2);
        //cerr << "Dot product is " << abs(axisEnd2.dot(tmpVect)) << endl;
        //cerr << "Closest to plane value was " << closestToPlane << endl;
        if(abs(axisEnd2.dot(tmpVect)) < closestToPlane)
        {
            closestToPlane = abs(axisEnd2.dot(tmpVect));
            bestEnd = tmpEnd;
            //cerr << "Best end is " << bestEnd.x << " " << bestEnd.y << " " << bestEnd.z << endl;
        }
        //else cerr << "End did not improve." << endl;
    }
    axisPoints.insert(axisPoints.begin(), bestEnd);
      // cerr << "Final choice is " << axisPoints[0].x << " " << axisPoints[0].y << " " << axisPoints[0].z << endl;
    return;
}

void Axis::reverseOrder()
{
    vector<Coordinate> tempAxisPoints;
    for (int i = axisPoints.size()-1; i >=0; i--)
        tempAxisPoints.push_back(axisPoints[i]);

    axisPoints.clear();
    for(int i = 0; i < tempAxisPoints.size(); i++)
        axisPoints.push_back(tempAxisPoints[i]);

    return;
}

bool Axis::angle()
{
    double V1X = 0;
    double V1Y = 0;
    double V1Z = 0;
    double V2X = 0;
    double V2Y = 0;
    double V2Z = 0;
    double magnitudeV1 = 0;
    double magnitudeV2 = 0;
    double dotProduct = 0;
    double angle = 0;


    V1X = axisPoints[0].x - axisPoints[1].x;
    V1Y = axisPoints[0].y - axisPoints[1].y;
    V1Z = axisPoints[0].z - axisPoints[1].z;
    V2X = axisPoints[axisPoints.size()-2].x - axisPoints[axisPoints.size()-1].x;
    V2Y = axisPoints[axisPoints.size()-2].y - axisPoints[axisPoints.size()-1].y;
    V2Z = axisPoints[axisPoints.size()-2].z - axisPoints[axisPoints.size()-1].z;
    magnitudeV1 = sqrt(pow(V1X,2)+pow(V1Y,2)+pow(V1Z,2));
    magnitudeV2 = sqrt(pow(V2X,2)+pow(V2Y,2)+pow(V2Z,2));
    dotProduct = V1X*V2X+V1Y*V2Y+V1Z*V2Z;
    angle = acos(dotProduct/(magnitudeV1*magnitudeV2));
    //cout << angle << endl;
    if(angle < 2.3 && angle > .7)
        return true;
    return false;

}

void Axis::splitFirst()
{
    double V1X = 0;
    double V1Y = 0;
    double V1Z = 0;
    double V2X = 0;
    double V2Y = 0;
    double V2Z = 0;
    double magnitudeV1 = 0;
    double magnitudeV2 = 0;
    double dotProduct = 0;
    double angle = 0;
    double orignalAngle = 0;
    V1X = axisPoints[0].x - axisPoints[1].x;
    V1Y = axisPoints[0].y - axisPoints[1].y;
    V1Z = axisPoints[0].z - axisPoints[1].z;
    magnitudeV1 = sqrt(pow(V1X,2)+pow(V1Y,2)+pow(V1Z,2));
    vector<Coordinate> tempPoints;
    tempPoints.push_back(axisPoints[0]);

    for(int i = 1; i < axisPoints.size(); i++)
    {
        V2X = axisPoints[i].x - axisPoints[i+1].x;
        V2Y = axisPoints[i].y - axisPoints[i+1].y;
        V2Z = axisPoints[i].z - axisPoints[i+1].z;

        magnitudeV2 = sqrt(pow(V2X,2)+pow(V2Y,2)+pow(V2Z,2));
        dotProduct = V1X*V2X+V1Y*V2Y+V1Z*V2Z;
        angle = acos(dotProduct/(magnitudeV1*magnitudeV2));
        tempPoints.push_back(axisPoints[i]);

        //cout << angle << endl;
        if(angle < 2.8 && angle > .34)
        {
            /*
            ofstream out;
            out.open("C:/Users/Malsa/Desktop/TESTTTT.pdb");
             for( int j = 0; j < i; j++ ) {
                                out << "ATOM  "
                                           << setw(5) << j
                                           << setw(4) << " CA "
                                           << "GLY "
                                           << "A"
                                           << setw(4) << j << "    "
                                           << setw(9) << fixed << setprecision(3) << axisPoints[j].x
                                           << setw(9) << fixed << setprecision(3) << axisPoints[j].y
                                           << setw(9) << fixed << setprecision(3) << axisPoints[j].z
                                           << endl;

                        }
            out.close();
            */
            i = axisPoints.size();
            axisPoints = tempPoints;
        }
    }
}

void Axis::splitSecond()
{
    double V1X = 0;
    double V1Y = 0;
    double V1Z = 0;
    double V2X = 0;
    double V2Y = 0;
    double V2Z = 0;
    double magnitudeV1 = 0;
    double magnitudeV2 = 0;
    double dotProduct = 0;
    double angle = 0;
    double orignalAngle = 0;
    V1X = axisPoints[axisPoints.size()-2].x - axisPoints[axisPoints.size()-1].x;
    V1Y = axisPoints[axisPoints.size()-2].y - axisPoints[axisPoints.size()-1].y;
    V1Z = axisPoints[axisPoints.size()-2].z - axisPoints[axisPoints.size()-1].z;
    magnitudeV1 = sqrt(pow(V1X,2)+pow(V1Y,2)+pow(V1Z,2));
    vector<Coordinate> tempPoints;
    tempPoints.push_back(axisPoints[axisPoints.size()-1]);

    for(int i = axisPoints.size()-2; i > 0; i--)
    {
        V2X = axisPoints[i].x - axisPoints[i-1].x;
        V2Y = axisPoints[i].y - axisPoints[i-1].y;
        V2Z = axisPoints[i].z - axisPoints[i-1].z;

        magnitudeV2 = sqrt(pow(V2X,2)+pow(V2Y,2)+pow(V2Z,2));
        dotProduct = V1X*V2X+V1Y*V2Y+V1Z*V2Z;
        angle = acos(dotProduct/(magnitudeV1*magnitudeV2));
        tempPoints.push_back(axisPoints[i]);
        //cout << angle << endl;
        if(angle < 2.1 && angle > .34)
        {
            /*
            ofstream out;
            out.open("C:/Users/Malsa/Desktop/TESTTTT2.pdb");
             for( int j = axisPoints.size()-1; j > i; j--) {
                                out << "ATOM  "
                                           << setw(5) << j
                                           << setw(4) << " CA "
                                           << "GLY "
                                           << "A"
                                           << setw(4) << j << "    "
                                           << setw(9) << fixed << setprecision(3) << axisPoints[j].x
                                           << setw(9) << fixed << setprecision(3) << axisPoints[j].y
                                           << setw(9) << fixed << setprecision(3) << axisPoints[j].z
                                           << endl;

                        }
            out.close();
            */
            i = 0;
            axisPoints = tempPoints;
        }
    }
}

void Axis::splitFirstHalf()
{
    vector<Coordinate> tempPoints;

    for(int i = 0; i < axisPoints.size()/2; i++)
    {
        tempPoints.push_back(axisPoints[i]);
    }
    axisPoints = tempPoints;
}

void Axis::splitSecondHalf()
{
    vector<Coordinate> tempPoints;

    for(int i = axisPoints.size(); i > axisPoints.size()/2; i--)
    {
        tempPoints.push_back(axisPoints[i]);
    }
    axisPoints = tempPoints;
}

Coordinate Axis::firstPoint()
{
    return axisPoints[0];
}

Coordinate Axis::lastPoint()
{
    return axisPoints[axisPoints.size()-1];
}

void Axis::printAsPnts (string outputFileName)
{
    ofstream out;
    out.open(outputFileName.c_str());
     for( int j = 0; j < axisPoints.size(); j++ ) {
                        out << "ATOM  "
                                   << setw(5) << j
                                   << setw(4) << " CA "
                                   << "GLY "
                                   << "A"
                                   << setw(4) << j << "    "
                                   << setw(9) << fixed << setprecision(3) << axisPoints[j].x
                                   << setw(9) << fixed << setprecision(3) << axisPoints[j].y
                                   << setw(9) << fixed << setprecision(3) << axisPoints[j].z
                                   << endl;

                }
    out.close();
}

void Axis::printAsPnts2 (string outputFileName)
{
    ofstream out;
    out.open(outputFileName.c_str());
     for( int j = 0; j < axisPoints.size(); j++ ) {
                        out << "ATOM  "
                                   << setw(5) << j
                                   << setw(4) << " H  "
                                   << "HOH "
                                   << "E"
                                   << setw(4) << j << "    "
                                   << setw(9) << fixed << setprecision(3) << axisPoints[j].x
                                   << setw(9) << fixed << setprecision(3) << axisPoints[j].y
                                   << setw(9) << fixed << setprecision(3) << axisPoints[j].z
                                   << endl;

                }
    out.close();
}

void Axis::catmullRom(double stepSize)
{
    struct Coordinate p0;
    struct Coordinate p1;
    struct Coordinate p2;
    struct Coordinate p3;
    double approxDistance = 0;
    int numInterpPoints = 0;
    vector<Coordinate> tempAxis;
    Coordinate v;
    //cout << axisPoints[0].x << " " << axisPoints[0].y << " " << axisPoints[0].z << endl;
    tempAxis.push_back(axisPoints[0]);
    axisPoints.insert(axisPoints.begin(), *(axisPoints.begin()));
    axisPoints.push_back(axisPoints[axisPoints.size()-1]);
    for (int j = 1; j < axisPoints.size()-2; j++)
        {
            approxDistance = getDistance(axisPoints[j], axisPoints[j+1]);
            //cerr << "Approximate distance: " << approxDistance << endl;
            numInterpPoints = approxDistance/stepSize;
            //cerr << "Number of interpolated points: " << numInterpPoints << endl;
            p0 = axisPoints[j-1];
            p1 = axisPoints[j];
            p2 = axisPoints[j+1];
            p3 = axisPoints[j+2];
            for (double t = 0.0; t < 1.00; t = t + 1.0/numInterpPoints)
            {
                double t2 = t*t;
                double t3 = t*t*t;
                v.x = 0.5 *((2 * p1.x) + (-p0.x + p2.x) * t + (2*p0.x - 5*p1.x + 4*p2.x - p3.x) * t2 + (-p0.x + 3*p1.x- 3*p2.x + p3.x) * t3);
                v.y = 0.5 *((2 * p1.y) + (-p0.y + p2.y) * t + (2*p0.y - 5*p1.y + 4*p2.y - p3.y) * t2 + (-p0.y + 3*p1.y- 3*p2.y + p3.y) * t3);
                v.z = 0.5 *((2 * p1.z) + (-p0.z + p2.z) * t + (2*p0.z - 5*p1.z + 4*p2.z - p3.z) * t2 + (-p0.z + 3*p1.z- 3*p2.z + p3.z) * t3);
                tempAxis.push_back(v);
            }
        }
    axisPoints.clear();

    for (int i = 1; i < tempAxis.size(); i++) ///
        axisPoints.push_back(tempAxis[i]);

    if (tempAxis.size() == 1 && axisPoints.size() == 0)
        axisPoints.push_back(tempAxis[0]);

    return;
}

void Axis::alignDirectionHels(Axis* qptr)
{
    vector<Coordinate> tempAxis;
    Vectors p (axisPoints[0], axisPoints[axisPoints.size()-1]);
    Vectors q1 (qptr->axisPoints[0], qptr->axisPoints[qptr->axisPoints.size()-1]);
    Vectors q2 (qptr->axisPoints[qptr->axisPoints.size()-1], qptr->axisPoints[0]);

    if (p.dot(q1) < p.dot(q2))
        {qptr->reverseOrder();
        cerr << "Reversed this helix." << endl;
        }

    return;
}

double Axis::getArcLength(vector<Coordinate>::iterator p1, vector<Coordinate>::iterator p2)
{
    double arcLength = 0;
    bool flip = false;
    vector<Coordinate>::iterator prevItr;
    vector<Coordinate>::iterator nextItr;

    int p1ind, p2ind, index;
    index = 0;
    p1ind = 0;
    p2ind = 0;

    for(prevItr = axisPoints.begin(); prevItr != axisPoints.end(); prevItr++)
    {
        index++;
        if (prevItr == p1)
            p1ind = index;
        if (prevItr == p2)
            p2ind = index;
    }


    for (nextItr = axisPoints.begin(); nextItr != p1; nextItr++) // if p1 comes after p2 in the list
        if(nextItr == p2)
        {
            prevItr = p2;
            p2 = p1;
            p1 = prevItr;
            flip = true;
            break;
        }

    nextItr = p1;
    nextItr++;
    index = 0;
    for (prevItr = p1; (prevItr != p2) && (nextItr != axisPoints.end()); prevItr++)
    {
        index++;
        //cerr << index << " steps currently calculated for arclength" << endl;
        arcLength += getDistance(*prevItr, *nextItr);
        nextItr++;
    }

    if (flip == true)
        arcLength *= -1;

    return arcLength;
}

Displacement Axis::findCrossDisplacement(Axis* qptr, double stepSize)
{
    Displacement displaced;
    double a, b;
    vector<Coordinate>::iterator pc, qc, pi, qi, qiTemp;
    Axis s;
    double distance = 0;
    double minDistance = 99999;
    double currentArc = 0;
    double closestArc = 99999;
    double sum = 0;
    bool startS = false;

    for (pi = axisPoints.begin(); pi != axisPoints.end(); pi++)
        for (qi = qptr->axisPoints.begin(); qi != qptr->axisPoints.end(); qi++)
        {
            distance = getDistance(*pi, *qi);
            if (distance < minDistance)
            {
                minDistance = distance;
                pc = pi;
                qc = qi;
            }
        }

    int closeIndexP = 0;
    int closeIndexQ = 0;
    for(pi = axisPoints.begin(); pi != axisPoints.end(); pi++)
        {
            if (pi == pc)
                break;
            closeIndexP++;
        }
    for(qi = qptr->axisPoints.begin(); qi != qptr->axisPoints.end(); qi++)
        {
            if (qi == qc)
                break;
            closeIndexQ++;
        }

    cerr << "Closest point on p is " << pc->x << " " << pc->y << " " << pc->z << " at index " << closeIndexP << endl;
    cerr << "Closest point on q is " << qc->x << " " << qc->y << " " << qc->z << " at index " << closeIndexQ << endl;
    cerr << endl;

    double a1 = -abs(getArcLength(axisPoints.begin(), pc));
    double a2 = -abs(qptr->getArcLength(qptr->axisPoints.begin(), qc));
    double b1 = abs(getArcLength(axisPoints.begin(), axisPoints.end()-1)) - abs(getArcLength(axisPoints.begin(), pc));
    double b2 = abs(qptr->getArcLength(qptr->axisPoints.begin(), qptr->axisPoints.end()-1)) - abs(qptr->getArcLength(qptr->axisPoints.begin(), qc));

    a = max(a1, a2);
    b = min(b1, b2);
    cerr << "b = " << b << endl;
    cerr << "a = " << a << endl;
    int testIndexP = 0;
    cerr << "Number of points in p is " << axisPoints.size() << endl;
    cerr << "Number of points in q is " << qptr->axisPoints.size() << endl << endl;
    for (pi = axisPoints.begin(); pi != axisPoints.end(); pi++)
    {
        closestArc = 999999;
        currentArc = getArcLength(pi, pc);
        //cerr << "Matching an arc on p of length " << currentArc << " at index " << testIndexP << endl;

        int testIndexQ = 0;
        for(qiTemp = qptr->axisPoints.begin(); qiTemp != qptr->axisPoints.end(); qiTemp++)
        {
            if (abs(qptr->getArcLength(qiTemp, qc) - currentArc) < closestArc)
            {
                closestArc = abs(qptr->getArcLength(qiTemp, qc) - currentArc);
                qi = qiTemp;
            }
        }
        //cerr << "Found the closest matching arc; the difference in arc lengths is " << closestArc << endl;
        if (closestArc < 2*stepSize) // iffy
            {
                for(qiTemp = qptr->axisPoints.begin(); qiTemp != qptr->axisPoints.end(); qiTemp++)
                {
                    if (qiTemp == qi)
                        break;
                        testIndexQ++;
                }
                //cerr << "Matched with an arc on q of length " << qptr->getArcLength(qi, qc) << " at index " << testIndexQ << endl << endl;
                startS = true;
                sum +=pow(getDistance(*pi,*qi), 2);
                s.axisPoints.push_back(*pi);
            }
        //else if (startS == false)
            //cerr << "p has started before q. Not starting line s yet." << endl << endl;
        else if (startS == true)
            {cerr << "q has ended before p. Ending line s." << endl << endl; break;}

        testIndexP++;
    }
    if (b-a == 0)
        cerr << "Insufficient length for cross-displacement comparison. Continuing..." << endl;
    displaced.crossDisplace = sqrt(sum*stepSize/(b-a));
    return displaced;
}

Displacement Axis::findDisplacement(Axis* qptr, double stepSize)
{
    cerr << "For true axis p and traced axis q: " << endl;
    Displacement displaced;
    double a, b;
    vector<Coordinate>::iterator pc, qc, pi, qi, qiTemp;
    Axis s;
    double distance = 0;
    double minDistance = 99999;
    double currentArc = 0;
    double closestArc = 99999;
    double sum = 0;
    bool startS = false;

    for (pi = axisPoints.begin(); pi != axisPoints.end(); pi++)
        for (qi = qptr->axisPoints.begin(); qi != qptr->axisPoints.end(); qi++)
        {
            distance = getDistance(*pi, *qi);
            if (distance < minDistance)
            {
                minDistance = distance;
                pc = pi;
                qc = qi;
            }
        }

    int closeIndexP = 0;
    int closeIndexQ = 0;
    for(pi = axisPoints.begin(); pi != axisPoints.end(); pi++)
        {
            if (pi == pc)
                break;
            closeIndexP++;
        }
    for(qi = qptr->axisPoints.begin(); qi != qptr->axisPoints.end(); qi++)
        {
            if (qi == qc)
                break;
            closeIndexQ++;
        }

    cerr << "Closest point on p is " << pc->x << " " << pc->y << " " << pc->z << " at index " << closeIndexP << endl;
    cerr << "Closest point on q is " << qc->x << " " << qc->y << " " << qc->z << " at index " << closeIndexQ << endl;
    cerr << endl;

    double a1 = -abs(getArcLength(axisPoints.begin(), pc));
    double a2 = -abs(qptr->getArcLength(qptr->axisPoints.begin(), qc));
    double b1 = abs(getArcLength(axisPoints.begin(), axisPoints.end()-1)) - abs(getArcLength(axisPoints.begin(), pc));
    double b2 = abs(qptr->getArcLength(qptr->axisPoints.begin(), qptr->axisPoints.end()-1)) - abs(qptr->getArcLength(qptr->axisPoints.begin(), qc));

    a = max(a1, a2);
    b = min(b1, b2);
    cerr << "b = " << b << endl;
    cerr << "a = " << a << endl;
    int testIndexP = 0;
    cerr << "Number of points in p is " << axisPoints.size() << endl;
    cerr << "Number of points in q is " << qptr->axisPoints.size() << endl << endl;
    for (pi = axisPoints.begin(); pi != axisPoints.end(); pi++)
    {
        closestArc = 999999;
        currentArc = getArcLength(pi, pc);
        //cerr << "Matching an arc on p of length " << currentArc << " at index " << testIndexP << endl;

        int testIndexQ = 0;
        for(qiTemp = qptr->axisPoints.begin(); qiTemp != qptr->axisPoints.end(); qiTemp++)
        {
            if (abs(qptr->getArcLength(qiTemp, qc) - currentArc) < closestArc)
            {
                closestArc = abs(qptr->getArcLength(qiTemp, qc) - currentArc);
                qi = qiTemp;
            }
        }
        //cerr << "Found the closest matching arc; the difference in arc lengths is " << closestArc << endl;
        if (closestArc < 2*stepSize) // iffy
            {
                for(qiTemp = qptr->axisPoints.begin(); qiTemp != qptr->axisPoints.end(); qiTemp++)
                {
                    if (qiTemp == qi)
                        break;
                        testIndexQ++;
                }
                //cerr << "Matched with an arc on q of length " << qptr->getArcLength(qi, qc) << " at index " << testIndexQ << endl << endl;
                startS = true;
                sum +=pow(getDistance(*pi,*qi), 2);
                s.axisPoints.push_back(*pi);
            }
        //else if (startS == false)
            //cerr << "p has started before q. Not starting line s yet." << endl << endl;
        else if (startS == true)
            {cerr << "q has ended before p. Ending line s." << endl << endl; break;}

        testIndexP++;
    }
    if (b-a == 0)
        cerr << "Insufficient length for cross-displacement comparison. Continuing..." << endl;
    displaced.crossDisplace = sqrt(sum*stepSize/(b-a));


    displaced.lengthAxis1 = getArcLength(axisPoints.begin(), axisPoints.end()-1);
    displaced.lengthAxis2 = qptr->getArcLength(qptr->axisPoints.begin(), qptr->axisPoints.end()-1);
    displaced.twoWayDistance = s.twoWayDistance(qptr);
    displaced.lengthComparison = s.getArcLength(s.axisPoints.begin(), s.axisPoints.end()-1);
    displaced.longDisplace = displaced.lengthAxis1 + displaced.lengthAxis2 - 2*displaced.lengthComparison;
        if(displaced.longDisplace < 0)
        {
            cerr << "Lines match through near the entirety of their lengths." << endl;
            displaced.lengthComparison = min(displaced.lengthAxis1, displaced.lengthAxis2);
            displaced.longDisplace = displaced.lengthAxis1 + displaced.lengthAxis2 - 2*displaced.lengthComparison;
        }
    cerr << "Length of p is " << displaced.lengthAxis1 << endl;
    cerr << "Length of q is " << displaced.lengthAxis2 << endl;
    cerr << "Length of s is " << displaced.lengthComparison << endl;

    cerr << "Cross displacement is " << displaced.crossDisplace << endl;
    cerr << "Longitudinal displacement is " << displaced.longDisplace << endl;

    displaced.lengthProportion = displaced.longDisplace/(displaced.longDisplace+b-a);
    cerr << "Proportion of incorrect length to combined length of axes is " << displaced.lengthProportion << endl;

    return displaced;
}

double Axis::twoWayDistance(Axis* qptr)
{
    double average1 = 0;
    double average2 = 0;
    double distance = 0;
    double closestDistance = 99999;

    for (int i = 0; i < axisPoints.size(); i++)
    {
        for (int j = 0; j < qptr->axisPoints.size(); j++)
        {
            if (getDistance(axisPoints[i], qptr->axisPoints[j]) < closestDistance)
                closestDistance = getDistance(axisPoints[i], qptr->axisPoints[j]);
        }
        average1 += closestDistance;
        closestDistance = 99999;
    }
    average1 /= axisPoints.size();

    for (int j = 0; j < qptr->axisPoints.size(); j++)
    {
        for(int i = 0; i < axisPoints.size(); i++)
            if(getDistance(axisPoints[i], qptr->axisPoints[j]) < closestDistance)
                closestDistance = getDistance(axisPoints[i], qptr->axisPoints[j]);
        average2 += closestDistance;
        closestDistance = 99999;
    }
    average2 /= qptr->axisPoints.size();

    distance = (average1 + average2)/2;

    return distance;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
// print out accuracy statistics
int PrintSpecificity(vector<Axis>* traces, vector<int>* specArray, Protein pdb, ofstream& out, int helixOffset)
{
    int hlxAA=0, fpHLXAA=0;  // # of total, fp-false positive
    int totalAA = pdb.numOfAA();
    bool matched = false;

    for (int n=0; n<pdb.numOfAA(); n++)
    {
        pdb.setBBCenter(n);  // BB center of this true sheet AA
        for (int m=0; m<helixOffset; m++)
        {
            matched = false;
            if(n < pdb.hlces[m].startIndx || n > pdb.hlces[m].endIndx)
            {
                for(int i = 0; i < traces->at(m).axisPoints.size(); i++)
                {
                    Coordinate p = traces->at(m).axisPoints[i];
                    if (getDistance(p, pdb.AAs[n].coord) <= 2.5) // was 3
                    {
                        if(matched == false)
                        {
                            fpHLXAA++;  // count for helix false positive
                            specArray->at(m)++;
                        }
                        matched = true;
                        break;
                    }
                }
            }
            else
                hlxAA++;
        }
    }
    out <<"   Specificity of Helix= "<<fixed<<setprecision(2)<<(double)(100.0-(double)((double)fpHLXAA*100/(double)(totalAA-hlxAA)))<<"%"<<endl;
    return (totalAA-hlxAA);
}

int PrintSpecificityStrands(vector<Axis>* traces, vector<int>* specArray, Protein pdb, ofstream& out, int helixOffset)
{
    int StrandAA=0, fpStrandAA=0;  // # of total, fp-false positive
    int totalAA = pdb.numOfAA();
    bool matched = false;

    for (int n=0; n<pdb.numOfAA(); n++)
    {
        pdb.setBBCenter(n);  // BB center of this true helix AA
        for (int m=helixOffset; m<traces->size(); m++)
        {
            matched = false;
            if(n < pdb.sheets[m].startIndx || n > pdb.sheets[m].endIndx)
            {
                for(int i = 0; i < traces->at(m).axisPoints.size(); i++)
                {
                    Coordinate p = traces->at(m).axisPoints[i];
                    if (getDistance(p, pdb.AAs[n].coord) <= 2.5) // was 3
                    {
                        if(matched == false)
                        {
                            fpStrandAA++;  // count for helix false positive
                            specArray->at(m)++;
                        }
                        matched = true;
                        break;
                    }
                }
            }
            else
                StrandAA++;
        }
    }
    out <<"   Specificity of Stands= "<<fixed<<setprecision(2)<<(double)(100.0-(double)((double)fpStrandAA*100/(double)(totalAA-StrandAA)))<<"%"<<endl;
    return (totalAA-StrandAA);
}

void ModifiedSpecificity(vector<Axis>* traces, vector<int>* specArray, Protein pdb, ofstream& out, int helixOffset)
{
    for(int i = 0; i < traces->size(); i++)
    {

    }
}

double PrintSensitivity(Axis* trace, int startIndx, int endIndx, Protein pdb, double radius, ofstream& out)
{
    int hlxAA=0, tpHLXAA=0; // # of total, tp-true positive
    int totalAA = pdb.numOfAA();
    bool matched = false;

    //cout << startIndx << " " << endIndx << endl;
    for (int n=startIndx; n<endIndx; n++) // Double check that we want a <= instead of a <
    {
        matched = false;

        //cout << pdb.AAs[n].SStype << endl;
        if (pdb.AAs[n].SStype == 'H')
        {
            hlxAA++;
            pdb.setBBCenter(n);  // BB center of this true helix AA

            for (double i = 0; i < trace->axisPoints.size(); i++)
            {
                    Coordinate p = trace->axisPoints[i];
                    if (getDistance(p, pdb.AAs[n].coord) <= radius)
                    {
                        tpHLXAA++; // count for helix true positive
                        matched = true;
                        break;
                    }
            }
        }
    }
    //cout << tpHLXAA << " " << hlxAA << endl;
    out <<fixed<<setprecision(2)<<(double)((double)tpHLXAA*100/(double)hlxAA)<<"%" << ", ";
    return ((double)tpHLXAA/(double)hlxAA);
}
double PrintSensitivityStrands(Axis* trace, int startIndx, int endIndx, Protein pdb, double radius, ofstream& out)
{
    int hlxAA=0, tpHLXAA=0; // # of total, tp-true positive
    int totalAA = pdb.numOfAA();
    bool matched = false;

    //cout << startIndx << " " << endIndx << endl;
    for (int n=startIndx; n<endIndx; n++) // Double check that we want a <= instead of a <
    {
        matched = false;

        //cout << pdb.AAs[n].SStype << endl;
        if (pdb.AAs[n].SStype == 'S')
        {
            hlxAA++;
            pdb.setBBCenter(n);  // BB center of this true helix AA

            for (double i = 0; i < trace->axisPoints.size(); i++)
            {
                    Coordinate p = trace->axisPoints[i];
                    if (getDistance(p, pdb.AAs[n].coord) <= radius)
                    {
                        tpHLXAA++; // count for helix true positive
                        matched = true;
                        break;
                    }
            }
        }
    }
    //cout << tpHLXAA << " " << hlxAA << endl;
    out <<fixed<<setprecision(2)<<(double)((double)tpHLXAA*100/(double)hlxAA)<<"%" << ", ";
    return ((double)tpHLXAA/(double)hlxAA);
}

int specificityBoundaries(Protein pdb, Axis* trueHels, bool helix, int current)
{
    int numInBox = 0;
    double changeFrontX = 0;
    double changeFrontY = 0;
    double changeFrontZ = 0;
    double changeLastX = 0;
    double changeLastY = 0;
    double changeLastZ = 0;

    changeFrontX = trueHels->axisPoints[5].x - trueHels->axisPoints[0].x;
    changeFrontY = trueHels->axisPoints[5].y - trueHels->axisPoints[0].y;
    changeFrontZ = trueHels->axisPoints[5].z - trueHels->axisPoints[0].z;
    changeLastX = trueHels->axisPoints[trueHels->axisPoints.size()-6].x - trueHels->axisPoints[trueHels->axisPoints.size()-1].x;
    changeLastY = trueHels->axisPoints[trueHels->axisPoints.size()-6].y - trueHels->axisPoints[trueHels->axisPoints.size()-1].y;
    changeLastZ = trueHels->axisPoints[trueHels->axisPoints.size()-6].z - trueHels->axisPoints[trueHels->axisPoints.size()-1].z;
    Axis temp;

    for(int i = 0; i < 15; i++)
    {
        Coordinate tempC;
        tempC.x = trueHels->axisPoints[0].x-changeFrontX*i;
        tempC.y = trueHels->axisPoints[0].y-changeFrontY*i;
        tempC.z = trueHels->axisPoints[0].z-changeFrontZ*i;
        temp.addPoint(tempC);
    }
    for(int i = 0; i < trueHels->axisPoints.size(); i++)
    {
        temp.addPoint(trueHels->axisPoints[i]);
    }
    for(int i = 0; i < 15; i++)
    {
        Coordinate tempC;
        tempC.x = trueHels->axisPoints[trueHels->axisPoints.size()-1].x-changeLastX*i;
        tempC.y = trueHels->axisPoints[trueHels->axisPoints.size()-1].y-changeLastY*i;
        tempC.z = trueHels->axisPoints[trueHels->axisPoints.size()-1].z-changeLastZ*i;
        temp.addPoint(tempC);
    }
    temp.catmullRom(.1);
    //string id_str;
    //stringstream ss;
    //ss << current;
    //id_str = ss.str();
    //string s = "C:\\Users\\dhaslam\\Desktop\\1CHD\\test" + id_str +".pdb";
    //temp.printAsPnts(s.c_str());

    for(int i = 0; i < pdb.numOfAA(); i++)
    {
        for(int j = 0; j < temp.axisPoints.size(); j++)
        {
            Coordinate p = temp.axisPoints[j];
            if (getDistance(p, pdb.AAs[i].coord) <= 10)
            {
                numInBox++;
                break;
            }
        }
    }

    return numInBox;



}

#endif // AXIS_H_INCLUDED
