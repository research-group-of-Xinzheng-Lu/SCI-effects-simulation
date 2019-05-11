/********************************************************
* @file    : PostprocessSpeed_MDOF.cpp
* @brief   : 	
* @details : 	
* @author  : Chinler
* @date    : 2018-10-4
*********************************************************/
#include <iostream>
#include <windows.h>
#include <sstream>
#include <fstream>
#include<string> 
#include<vector>  
#include <numeric>
#include "building.h"

#define PI 3.1415927
using namespace std;
int countline(string filename)
{
	int line = 0;
	ifstream ipt(filename.c_str());
	char c;
	while (ipt.get(c))
	{
		if (c == '\n')
			line++;
	}
	return line + 1;

}
string getOutputName();
string& trim(string &s)
{
	if (s.empty())
	{
		return s;
	}
	s.erase(0, s.find_first_not_of(" "));
	s.erase(s.find_last_not_of(" ") + 1);
	return s;
}
int main()
{
	string fileMPI = ".//FILES_MPI//mpitest0.txt";
	string fileCtrl = "BRCtrl.txt";
	vector <int> ID; 
	vector <int> numBuilding; //每个并行文件中建筑数量
	int BRoutput[6] = { 1,1,1,1,1,1 };
	int line=countline(fileMPI)-1; //并行数量
	ifstream inMPI(fileMPI.c_str());
	int tempID, tempNum;
	double totalTime, step;
	double tempTime,temp;
	for (int i=0;i<line;i++)
	{
		inMPI >> tempID >> tempNum;
		ID.push_back(tempID);
		numBuilding.push_back(tempNum);
	}
	int numB = accumulate(numBuilding.begin(), numBuilding.end(), 0);//建筑总数
	Building *build;
	build = new Building[numB];

	/* -------read BuildingInfo.txt--------- */
	ifstream inBuildingInfo("BuildingInfo.txt");
	string tempString;
	inBuildingInfo >> temp;
	for (int i=0;i<numB;i++)
	{
		inBuildingInfo >> build[i].name >> temp>> build[i].nStory>> build[i].storyheight>> build[i].area>>temp>> temp;
		for (int j = 0; j < build[i].nStory+1; j++)
			getline(inBuildingInfo, tempString);
	}
	


	ifstream inCtrl(fileCtrl.c_str());
	if (inCtrl) //if file exits,then read the file
	{
		for (int i=0;i<6;i++)
		{
			inCtrl >> BRoutput[i];
		}		
	} 	
	/* -------read config.txt--------- */
	string outputFileName=getOutputName();
	string monitorFilename = ".//" + outputFileName + "//MONITOR.INFO";
	ifstream inMonitor(monitorFilename.c_str());
	int stepGap = 0;
	inMonitor >> totalTime >> step >> stepGap;
	//ifstream inConfig("Config.txt");
	//for (int i=0;i<7;i++)
	//{
	//	getline(inConfig, tempString);
	//}
	//inConfig>> totalTime>> step;
	for (int i = 0; i < numB; i++)
	{
		build[i].Ntime = int(totalTime / step/ stepGap) +1;
		build[i].dt = step*stepGap;
	}
		
	int TimeN= int(totalTime / step) - 1;
	

/*************************Read GM information 2018/10/5 15:32 Begin**************************/
	string AcIpX = ".//MONITOR//AcIpX-";
	string AcIpY = ".//MONITOR//AcIpY-";
	int indexx = 0;
	for (int i = 0; i < line; i++) //core
	{
		string fileAcIpX = AcIpX + to_string(i + 1) + ".txt";
		string fileAcIpY = AcIpY + to_string(i + 1) + ".txt";
		ifstream infileAcIpX(fileAcIpX.c_str());
		ifstream infileAcIpY(fileAcIpY.c_str());
		
		for (int j = 0; j <TimeN; j++)
		{
			infileAcIpX >> tempTime;
			infileAcIpY >> tempTime;
			for (int k = 0; k < numBuilding[i]; k++) //building
			{
				infileAcIpX >>temp;				
				build[indexx+k].GMX.push_back(temp);
				infileAcIpY >> temp;
				build[indexx + k].GMY.push_back(temp);
			}
		}
		indexx = indexx+ numBuilding[i];
	}
	for (int i=0;i<numB;i++)
	{
		build[i].PGAX=0.0;
		build[i].PGAY = 0.0;
		for (int j=0;j<build[i].Ntime;j++)
		{
			if (build[i].PGAX<abs(build[i].GMX[j]))
			{
				build[i].PGAX =abs(build[i].GMX[j]);
			}
			if (build[i].PGAY < abs(build[i].GMY[j]))
			{
				build[i].PGAY = abs(build[i].GMY[j]);
			}
		}
	}


/*************************Read GM information 2018/10/5 15:32 End****************************/

	/* -------read results--------- */
	string TDispX = ".//MONITOR//TDisX-";
	string TDispY = ".//MONITOR//TDisY-";
	string SDispX = ".//MONITOR//SDisX-";
	string SDispY = ".//MONITOR//SDisY-";
	int index = 0;
	if (BRoutput[0]==1)//读取每栋建筑顶点位移时程
	{
		for (int i=0;i<line;i++)
		{
			string fileDispX = TDispX + to_string(i + 1) + ".txt";
			string fileDispY = TDispY + to_string(i + 1) + ".txt";
			ifstream infileDispX(fileDispX.c_str());
			ifstream infileDispY(fileDispY.c_str());

		}
	} 
	else if (BRoutput[0] == 2)//读取每栋每层的位移时程
	{
		
		for (int i = 0; i < line; i++)
		{
			string fileDispX = SDispX + to_string(i + 1) + ".txt";
			string fileDispY = SDispY + to_string(i + 1) + ".txt";
			ifstream infileDispX(fileDispX.c_str());
			ifstream infileDispY(fileDispY.c_str());
			for (int j = 0; j < numBuilding[i]; j++)
			{
				build[index + j].dispX = new double *[build[index + j].nStory];
				build[index + j].dispY = new double *[build[index + j].nStory];
				for (int l = 0; l < build[index + j].nStory; l++)
				{
					build[index + j].dispX[l] = new double[build->Ntime];
					build[index + j].dispY[l] = new double[build->Ntime];
				}
			}
			for (int j = 0; j < TimeN; j++)
			{
				infileDispX >> tempTime;
				infileDispY >> tempTime;
				for (int k = 0; k < numBuilding[i]; k++)
				{
					for (int m = 0; m < build[index + k].nStory; m++)
					{
						infileDispX >> build[index + k].dispX[m][j];
						infileDispY >> build[index + k].dispY[m][j];
					}
				}
			}			
			index = index + numBuilding[i];
		}
		for (int i=0;i<numB;i++)
		{
			build[i].getDrift();
			build[i].getResidual();
		}

	}else
	{
		cout << "Error in Displacement!" << endl;
	}


/*************************Read accleration 2018/10/5 16:11 Begin**************************/
	string TAccX = ".//MONITOR//TAccX-";
	string TAccY = ".//MONITOR//TAccY-";
	string SAccX = ".//MONITOR//SAccX-";
	string SAccY = ".//MONITOR//SAccY-";
	int indexAcc = 0;
	if (BRoutput[1] == 1)//读取每栋建筑顶点位移时程
	{
		for (int i = 0; i < line; i++)
		{
			string fileAccX = TAccX + to_string(i + 1) + ".txt";
			string fileAccY = TAccY + to_string(i + 1) + ".txt";
			ifstream infileAccX(fileAccX.c_str());
			ifstream infileAccY(fileAccY.c_str());

		}
	}
	else if (BRoutput[1] == 2)//读取每栋每层的位移时程
	{

		for (int i = 0; i < line; i++)
		{
			string fileAccX = SAccX + to_string(i + 1) + ".txt";
			string fileAccY = SAccY + to_string(i + 1) + ".txt";
			ifstream infileAccX(fileAccX.c_str());
			ifstream infileAccY(fileAccY.c_str());
			for (int j = 0; j < numBuilding[i]; j++)
			{
				build[indexAcc + j].accX = new double *[build[indexAcc + j].nStory];
				build[indexAcc + j].accY = new double *[build[indexAcc + j].nStory];
				for (int l = 0; l < build[indexAcc + j].nStory; l++)
				{
					build[indexAcc + j].accX[l] = new double[build->Ntime];
					build[indexAcc + j].accY[l] = new double[build->Ntime];
				}
			}
			for (int j = 0; j < build->Ntime; j++)
			{
				infileAccX >> tempTime;
				infileAccY >> tempTime;
				for (int k = 0; k < numBuilding[i]; k++)
				{
					for (int m = 0; m < build[indexAcc + k].nStory; m++)
					{
						infileAccX >> build[indexAcc + k].accX[m][j];
						infileAccY >> build[indexAcc + k].accY[m][j];
					}
				}
			}
			indexAcc = indexAcc + numBuilding[i];
		}
		for (int i = 0; i < numB; i++)
		{
			build[i].getPFA();
		}

	}
	else
	{
		cout << "Error in accleration!" << endl;
	}
/*************************Read accleration 2018/10/5 16:11 End****************************/
	 

	ofstream opt("edps.txt");
	for (int i=0;i<numB;i++)
	{
		for (int j=0;j<build[i].nStory;j++)
		{
			opt << build[i].edps.IDR[j] << "\t";
		}
		for (int j = 0; j < build[i].nStory+1; j++)
		{
			opt << build[i].edps.PFA[j] << "\t";
		}
		opt <<  "\n";
		
	}
	for (int i = 0; i < numB; i++)
	{
		string fileEDP = to_string(i+1)+"-EDP.json";
		build[i].CreateEDP(fileEDP.c_str());
		build[i].WriteDisp(i + 1);
		build[i].WriteAcc(i + 1);
	}

	
	return 0;

}
string getOutputName()
{
	string outputName,tempString;
	ifstream fin("speed.input");
	int line = countline("speed.input");
	
	int index = 0;
	for (int i=0;i<line;i++)
	{
		char filename[8]="";
		fin.get(filename,8);
		if (strcmp(filename, "MONFILE")==0)
		{
			index = i;
			getline(fin, tempString);
			break;
		}
		getline(fin, tempString);
	}

	outputName = trim(tempString);
	return outputName;
}

