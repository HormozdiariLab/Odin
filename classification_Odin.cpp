#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct samples{
	char sampleName[50];
	int caseControl; //1: case, 0: control
	int trainOrTest; //1: test, 0: train
	double values[50000];
	double distanceToCenter;
	int coveredInThreeSampleLoF;
	int filterSample;
}samples;

samples listSamples[2000];
int totalSamples;
int numberGene;
double distanceMatrix[2000][2000];
int casesInSphereCenter[2000]; // the center is the id 
double radiusOfSphere[2000];

#define maxCloRow 25000
double weightsCalcultaed[maxCloRow];


int centerId;
int casesInSphereId[1000];
int numberCasesInSphere=0;
int controlInSphereId[100];
int numberControlInSphere=0;
int sphereCount;

typedef struct spheres{
	double geneWeights[50000];
	char centerSampleName[50];
	int centerId; // the center id from the samples (listSamples)
	double radius; // the radius
	int numberCasesInside;
}spheres;

spheres listSpheres[20];



int compare (const void * a, const void * b)
{

double xx = *(double*)a; 
double yy = *(double*)b;
if (xx<yy) return -1;
if (xx>yy) return 1;
return 0;

//  return ( *(double*)a - *(double*)b );
}

int classify_NN(int sampleId)
{
double distance=0;
int returnValue=0;
double totalCount=0;
double distanceArray[totalSamples];
double tempScore[totalSamples];
double minDistanceControl=0;
int caseInside, maxCaseInside;


for (int count=0; count<totalSamples; count++)
{
	distanceArray[count]=0;
}

maxCaseInside=0;

for (int count=0; count<1; count++)
{
distance=0;
minDistanceControl=10000000; 
caseInside=0;
	for (int count2=0; count2<totalSamples; count2++)
	{
	distance=0;
	distanceArray[count2]=0;
	tempScore[count2]=0;
		//if (listSamples[count2].trainOrTest==0)
		for (int count3=0; count3<numberGene; count3++)
		{
		
			//distance=distance+listSpheres[count].geneWeights[count3]*fabs(listSamples[sampleId].values[count3]-listSamples[count2].values[count3]);
			distance=distance+1*fabs(listSamples[sampleId].values[count3]-listSamples[count2].values[count3]);
		
		}
		
		distanceArray[count2]=distance;	
		tempScore[count2]=distance;

		if((listSamples[count2].trainOrTest==0)&&(listSamples[count2].caseControl==0))
		{
			if (minDistanceControl>distance)
				minDistanceControl=distance;
		} 

	}

	

	for (int count2=0; count2<totalSamples; count2++)
	{
		if ((listSamples[count2].trainOrTest==0)&&(listSamples[count2].caseControl==1))
		{
			if (distanceArray[count2]<minDistanceControl)
			{
				//printf("%s %lf ", listSamples[count2].sampleName, distanceArray[count2]);
				caseInside++;
			}
		}
	}

	//printf("\n");

	if (maxCaseInside<caseInside)
		maxCaseInside=caseInside;

	distance=0;

	for (int count3=0; count3<numberGene; count3++)
	{
		distance=distance+listSpheres[count].geneWeights[count3]*fabs(listSamples[sampleId].values[count3]-listSamples[listSpheres[count].centerId].values[count3]);
	}


	//if (distanceArray[listSpheres[count].centerId]<listSpheres[count].radius)
	if (distance<listSpheres[count].radius)
	{
		if (caseInside>0)
			printf("SamplesInsideCenterBest %i %s %i %lf %lf\n", count, listSamples[sampleId].sampleName, caseInside, distance, listSpheres[count].radius );
	}//else printf("SamplesInsideCenterNotBest %i %s %i %lf %lf\n", count, listSamples[sampleId].sampleName, caseInside, distance, listSpheres[count].radius);

	
}


//printf("SamplesInsideCenter %s %i %lf\n", listSamples[sampleId].sampleName, maxCaseInside, minDistanceControl);
return 0;

}


int classify(int sampleId)
{
double  distance=0;
int returnValue=0;
double totalCount=0;


	for (int count=0; count<3; count++)
	{
	distance=0;
		for (int count2=0; count2<numberGene; count2++)
		{
			distance=distance+listSpheres[count].geneWeights[count2]*fabs(  listSamples[listSpheres[count].centerId].values[count2]-listSamples[sampleId].values[count2]);
		}

	//printf("distance %lf %lf\n", distance, listSpheres[count].radius);


	//if ((count == 0 && distance<1.5*listSpheres[count].radius)||(count==1 && distance<0.9*listSpheres[count].radius))
	if ((distance<listSpheres[count].radius))
	{
			
		//printf("distance %lf %lf %s %s %i\n", distance, listSpheres[count].radius, listSamples[sampleId].sampleName, listSamples[listSpheres[count].centerId].sampleName, listSpheres[count].numberCasesInside);
		totalCount=totalCount+(1-((double)distance/(double)listSpheres[count].radius))*listSpheres[count].numberCasesInside;
		returnValue=1;	
	}
	}

//if (returnValue==1)
//{
//	printf("CASE_PR %s %lf\n", listSamples[sampleId].sampleName, totalCount);
//}

return returnValue;
	
}



int calculateNumCasesInsideSphere()
{
double distance=0;
//printf("L199 %i\n", sphereCount);
	for (int count=0; count<totalSamples; count++)
	{
	if ((listSamples[count].caseControl==1)&&(listSamples[count].trainOrTest==0)) 
		{	
			for (int count2=0; count2<sphereCount; count2++)
			{
				distance=0;
				for (int count3=0; count3<numberGene; count3++)
				{
					distance=distance+listSpheres[count2].geneWeights[count3]*fabs(listSamples[listSpheres[count2].centerId].values[count3] - listSamples[count].values[count3]);

				}
				if (distance<listSpheres[count2].radius)
				{
					listSpheres[count2].numberCasesInside++;
				}
			}
		}
	}

	for (int count=0; count<sphereCount; count++)
	{
		//printf("Center %i\t%s\t%i\n", count, listSpheres[count].centerSampleName, listSpheres[count].numberCasesInside);
	}



}




int classification()
{
int isCase;
int countCases=0;
int countControl=0;
int casePredicedtCase=0;
int controlPredictedCase=0;
int coveredByThreeGeneFind=0;
int coveredByThreeGeneNotFind=0;
	for (int count=0; count<totalSamples; count++)
	{
		if (listSamples[count].trainOrTest==1)
		{
			//isCase=classify(count);
			isCase=classify_NN(count);
			if (listSamples[count].caseControl==1)
			{
		
				countCases++;
				if (isCase==1)
				{
					casePredicedtCase++;
				
					if (listSamples[count].coveredInThreeSampleLoF==1)
						coveredByThreeGeneFind++;			
				}else{
					
					if (listSamples[count].coveredInThreeSampleLoF==1)
						coveredByThreeGeneNotFind++;			
				}
			}else
			if (listSamples[count].caseControl==0)
			{
				countControl++;
				if (isCase==1)
					controlPredictedCase++;
			}
			//if (isCase==1)
			//printf("Predicted:%i\tReal:%i\n", isCase, listSamples[count].caseControl);
		}
	}

//printf("FINAL Result %i %i %f %i %i OnlyMyModel:%f OnlyMyModelTPR:%f UnionTwoModels:%f %i %i %f %f\n", countCases, casePredicedtCase,(float)casePredicedtCase/(float)countCases,  coveredByThreeGeneFind, coveredByThreeGeneNotFind, (float)(casePredicedtCase-coveredByThreeGeneFind)/(float)countCases, (float)(casePredicedtCase-coveredByThreeGeneFind)/(float)(countCases-coveredByThreeGeneFind-coveredByThreeGeneNotFind), (float)(casePredicedtCase+coveredByThreeGeneNotFind)/(float)countCases, countControl, controlPredictedCase,(float)controlPredictedCase/(float)countControl, (((float)casePredicedtCase/(float)countCases)/((float)controlPredictedCase/(float)countControl)) );

}



int main(int argv, char *argc[])
{
	FILE *fp=fopen(argc[1],"r"); // All Samples calculated
	numberGene=atoi(argc[2]);
	totalSamples=atoi(argc[3]);
	FILE *fp2=fopen(argc[4],"r");// the classification spheres (weights, centers, radius)
	FILE *fp3=fopen(argc[5],"r"); //Every samples name grouped as training (0) or testing sample (1)
	FILE *fp4=fopen(argc[6],"r");// the 3 gene LoF


	int sampleId, trainOrTest;
	char sampleName[50];
	for (int count=0; count<totalSamples; count++)
	{
		fscanf(fp, "%s\t%i", listSamples[count].sampleName, &(listSamples[count].caseControl));
		//printf("%s\n", listSamples[count].sampleName);
		listSamples[count].filterSample=0;
		for (int count2=0; count2<numberGene; count2++)
		{
			fscanf(fp,"\t%lf", &(listSamples[count].values[count2]));
		}
		fscanf(fp,"\n");

	}



	fscanf(fp2,"%i\n", &sphereCount);	
	for (int count=0; count<sphereCount; count++)
	{
		for (int count2=0; count2<numberGene; count2++)
		{
			fscanf(fp2, "%lf ", &(listSpheres[count].geneWeights[count2]));
		}
		fscanf(fp2, "%s %lf\n", listSpheres[count].centerSampleName, &(listSpheres[count].radius));
		//printf("L233 %s %lf\n", listSpheres[count].centerSampleName, listSpheres[count].radius);
		for (int count2=0; count2<totalSamples; count2++)
		{
			if (strcmp(listSamples[count2].sampleName, listSpheres[count].centerSampleName)==0)
			{
				listSpheres[count].centerId=count2;	
			}
		}
	}




	while(fscanf(fp3,"%s\t%i\n", sampleName, &trainOrTest)!=EOF)
	{
		for (int count=0; count<totalSamples; count++)
		{
			if (strcmp(listSamples[count].sampleName, sampleName)==0)
			{
				listSamples[count].trainOrTest=trainOrTest;
			}
		}
	}


	while(fscanf(fp4,"%s\n", sampleName)!=EOF)
	{
		for (int count=0; count<totalSamples; count++)
		{
			if (strcmp(sampleName, listSamples[count].sampleName)==0)
				listSamples[count].coveredInThreeSampleLoF=1;
		}
	}

	calculateNumCasesInsideSphere();
	classification();

}
