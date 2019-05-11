#ifndef BUILDING_H
#define BUILDING_H

#include <vector>
#include <string>
#include <algorithm>
using namespace std;

class Building
{
public:
    enum StruType{RM1, RM2, URM, C1, C2, C3, W1, W2, S1, S2, S3, S4, S5, PC1, PC2, MH, UNKNOWN}; //Hazus structural type
    enum BldgOccupancy{office, education, healthcare, hospitality, residence, retail, warehouse, research, unknown};

    struct EDP{     //engineering demand parameters
        vector<double> IDR; // Inter-story drift ratio. size = nStory
        vector<double> PFA; // Peak floor acceleration. size = nStory+1
        double residual;    // residual drift
        //vector<double> PFV; // Peak floor velocity. size = nStory
        //vector<double> rotation; // rotation. size = nStory
    };

    //basic building info
    int id;
    string name;
    StruType strutype;
    int year;
    BldgOccupancy occupancy;
    int nStory;
    double storyheight;   //unit: m
    double area;		//story area. unit: m^2
    double im;  //intensity measure

    double clpsMedian;					//collapse median. unit: m/s^2
    double clpsDispersion;				//collapse dispersion


    double dampingRatio = 0.0;
    double damageCriteria[4] = {0,0,0,0};
    double T0;
    double T2;
	double **dispX;
	double **accX;
	double **dispY;
	double **accY;
    //output info
    vector <int> damagestates;  //ds for each story
    int maxDamage;  // max damage of this building
    EDP edps;
	int Ntime;
	double dt;
	double PGAX;
	double PGAY;
	vector <double> GMX;
	vector <double> GMY;

    Building();
	void getDrift();
	void getPFA();
	void getResidual();
	void CreateEDP(const char *filenameEDP);
	void WriteDisp(int BuildingID);
	void WriteAcc(int BuildingID);


};

#endif // BUILDING_H
