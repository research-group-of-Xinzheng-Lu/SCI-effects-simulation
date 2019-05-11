#include "building.h"
#include <fstream>
#include <jansson.h> // for Json
Building::Building()
{

}
void Building::getDrift()
{
	double **InterStory;
	double *MaxInter;
	InterStory = new double*[nStory];
	MaxInter = new double[nStory];
	for (int i = 0; i < nStory; i++)
		MaxInter[i] = 0.0;
	for (int i = 0; i < nStory; i++)
	{ 
		InterStory[i] = new double[Ntime];
		for (int j=0;j<Ntime;j++)
		{
			if (i==0)
			{
				InterStory[i][j] = abs(dispX[i][j])/storyheight;
			}
			else {
				InterStory[i][j] = abs(dispX[i][j]- dispX[i-1][j])/storyheight;
			}	
			if (MaxInter[i]<InterStory[i][j])
			{
				MaxInter[i] = InterStory[i][j];
			}

			if (i == 0)
			{
				InterStory[i][j] = abs(dispY[i][j]) / storyheight;
			}
			else {
				InterStory[i][j] = abs(dispY[i][j] - dispY[i - 1][j]) / storyheight;
			}
			if (MaxInter[i] < InterStory[i][j])
			{
				MaxInter[i] = InterStory[i][j];
			}
			
		}
		edps.IDR.push_back(MaxInter[i]);
	}
}
void Building::getPFA()
{
	double *MaxPFA;
	MaxPFA = new double[nStory];
	for (int i = 0; i < nStory; i++)
		MaxPFA[i] = 0.0;
	PGAX < PGAY ? edps.PFA.push_back(PGAY) : edps.PFA.push_back(PGAX);
	for (int i = 0; i < nStory; i++)
	{
		for (int j = 0; j < Ntime; j++)
		{
			if (MaxPFA[i] < abs(accX[i][j]))
			{
				MaxPFA[i] = abs(accX[i][j]);
			}
			if (MaxPFA[i] < abs(accY[i][j]))
			{
				MaxPFA[i] = abs(accY[i][j]);
			}
			
		}
		edps.PFA.push_back(MaxPFA[i]);
	}
}
void Building::getResidual()
{
	double res = 0.0;
	for (int i = 0; i < nStory; i++)
	{
		if (i==0)
		{
			if (res < abs(dispX[i][Ntime - 1]) / storyheight)
			{
				res = abs(dispX[i][Ntime - 1]) / storyheight;
			}
		}else{
			if (res < abs(dispX[i][Ntime - 1] - dispX[i - 1][Ntime - 1]) / storyheight)
		{
			res = abs(dispX[i][Ntime - 1] - dispX[i - 1][Ntime - 1]) / storyheight;
		}
		}

		if (i == 0)
		{
			if (res < abs(dispY[i][Ntime - 1]) / storyheight)
			{
				res = abs(dispY[i][Ntime - 1]) / storyheight;
			}
		}
		else {
			if (res < abs(dispY[i][Ntime - 1] - dispY[i - 1][Ntime - 1]) / storyheight)
			{
				res = abs(dispY[i][Ntime - 1] - dispY[i - 1][Ntime - 1]) / storyheight;
			}
		}

	}
	edps.residual = res;
}

void Building::CreateEDP(const char *filenameEDP)
{
	// create output JSON object
	json_t *rootEDP = json_object();
	json_t *eventArray = json_array(); // for each analysis event
	// add the EDP for the event
	json_t *eventObj = json_object();
	json_object_set(eventObj, "name", json_string("name"));
	json_t *responsesArray = json_array(); // for each analysis event
	for (int i = 0; i < nStory; i++) {
		json_t *response = json_object();
		json_object_set(response, "type", json_string("max_drift"));
		json_object_set(response, "cline", json_integer(1));
		json_object_set(response, "floor1", json_integer(i + 1));
		json_object_set(response, "floor2", json_integer(i + 2));
		json_object_set(response, "scalar_data", json_real(edps.IDR[i]));
		json_array_append(responsesArray, response);
	}
	for (int i = 0; i < nStory+1; i++) {
		json_t *response = json_object();
		json_object_set(response, "type", json_string("max_abs_acceleration"));
		json_object_set(response, "cline", json_integer(1));
		json_object_set(response, "floor", json_integer(i + 1));
		json_object_set(response, "scalar_data", json_real(edps.PFA[i]));
		json_array_append(responsesArray, response);
	}
	json_t *response = json_object();
	json_object_set(response, "type", json_string("residual_disp"));
	json_object_set(response, "cline", json_integer(1));
	json_object_set(response, "floor", json_integer(nStory + 1));
	json_object_set(response, "scalar_data", json_real(edps.residual));
	json_object_set(eventObj, "responses", responsesArray);
	json_array_append(responsesArray, response);
	json_array_append(eventArray, eventObj);
	json_object_set(rootEDP, "EngineeringDemandParameters", eventArray);
	json_dump_file(rootEDP, filenameEDP, 0);
}



void Building::WriteDisp(int BuildingID)
{
	string filenameDispX = to_string(BuildingID) + "-DispX.txt";
	string filenameDispY = to_string(BuildingID) + "-DispY.txt";
	ofstream foutDispX(filenameDispX);
	ofstream foutDispY(filenameDispY);
	double dtX = -dt;
	double dtY = -dt;
	for (int i=0;i<Ntime;i++)
	{
		dtX = dtX + dt;
		dtY = dtY + dt;
		foutDispX << dtX  << "\t";
		foutDispY << dtY  << "\t";
		for (int j = 0; j < nStory; j++)
		{
			foutDispX << dispX[j][i] << "\t";
			foutDispY << dispY[j][i] << "\t";
		}
		foutDispX << "\n";
		foutDispY << "\n";
	}

}

void Building::WriteAcc(int BuildingID)
{
	string filenameAccX = to_string(BuildingID) + "-AccX.txt";
	string filenameAccY = to_string(BuildingID) + "-AccY.txt";
	ofstream foutAccX(filenameAccX);
	ofstream foutAccY(filenameAccY);
	double dtX = -dt;
	double dtY = -dt;
	for (int i = 0; i < Ntime; i++)
	{
		dtX = dtX + dt;
		dtY = dtY + dt;
		foutAccX << dtX << "\t";
		foutAccY << dtY << "\t";
		for (int j = 0; j < nStory; j++)
		{
			foutAccX << accX[j][i] << "\t";
			foutAccY << accY[j][i] << "\t";
		}
		foutAccX << "\n";
		foutAccY << "\n";
	}

}

