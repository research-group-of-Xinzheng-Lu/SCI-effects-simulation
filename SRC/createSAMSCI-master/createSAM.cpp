
#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <cctype>
#include "HazusSAM_Generator.h"
#include "Building.h"
#include <string>

int main()
{

//  if (argc != 5) {
//    printf("ERROR: correct usage: createSAM fileNameBIM fileNameEVENT fileNameSAM fileNameSAMSCI\n");
//    exit(0);
//  }
  int num=0;
  std::cout<<"Please input the number of buildings:"<<"\n";
  std::cin>>num;


//  char *filenameBIM = argv[1];
//  char *filenameEVENT = argv[2];
//  char *filenameSAM = argv[3];
  string fileNameSAMSCI="BuildingInfo.txt";
  ofstream outfile(fileNameSAMSCI.c_str());//ios::app表示在原文件末尾追加
  outfile<<num<<endl;
  HazusSAM_Generator* aim = new HazusSAM_Generator();

  for (int i=1;i<num+1;i++)
  {
      string filenameBIM=to_string(i)+"-BIM.json";
      Building *theBuilding = new Building();
      theBuilding->readBIM(filenameBIM.c_str());
      aim->CalcBldgPara(theBuilding);
      //theBuilding->writeSAM(filenameSAM);
      theBuilding->writeSAMSCI(fileNameSAMSCI.c_str());
  }
  for (int i=1;i<num+1;i++)
  {
      string filenameBIM=to_string(i)+"-BIM.json";
      Building *theBuilding = new Building();
      theBuilding->readBIM(filenameBIM.c_str());
      aim->CalcBldgPara(theBuilding);
      //theBuilding->writeSAM(filenameSAM);
      theBuilding->writeSAMdamage(fileNameSAMSCI.c_str());
  }
  delete aim;
  return 0;
}

