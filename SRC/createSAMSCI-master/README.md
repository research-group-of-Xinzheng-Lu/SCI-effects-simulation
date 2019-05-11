# createSAMSCI
Create input file for regional SCI computation based on the SimCenter Workflow


##input file

1-BIM.json 2-BIM.json 3-BIM.json 4-BIM.json 5-BIM.json ...


Input the number of buildings to run the program.



##Output file

BuildingInfo.txt


NOTE: the applications require the jansson lib be installed. jansson is 
a free BSD lib recomended online for dealing with reading/writing json files.
It is assumed installed in /usr/local/jansson as seen in Makefiles.


