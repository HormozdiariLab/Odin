#include <ilcplex/cplex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
/* Total number of dimensions allowed after dimension reduction is <10. In the LP the summation of weights of dimensions should be less than 5 (this can be changed as user wants)*/

FILE *fpOut;

typedef struct samples{
	char sampleName[50];
	int caseControl; // 1: case, 0: control
	int trainOrTest; //1: test, 0: train 
	double values[50000];
	double distanceToCenter;
	int filterSample;
}samples;

samples listSamples[3000];
int totalSamples;
int numberGene;
double distanceMatrix[3000][3000];
int casesInSphereCenter[3000]; // the center is the id 
double radiusOfSphere[3000];

#define maxCloRow 25000
double weightsCalcultaed[maxCloRow];


int centerId;
int casesInSphereId[1000];
int numberCasesInSphere=0;
int controlInSphereId[100];
int numberControlInSphere=0;



typedef struct spheres{
	double weights[50000];
	int centerId;
	double radius;
	int numCases;
}spheres;


spheres listSpheres[20];






int calculateWeights(int centerId, int *casesToConsider, int *controlToConsider, int numCase, int numControl)
{

	
int solstat;
double objval;
double *x = NULL;
double *pi = NULL;
double *slack = NULL;
double *dj = NULL;
int cur_numrows, cur_numcols;
int status=0;
int NumRows, NumCols, NumNZ;


double obj[maxCloRow];
double lb[maxCloRow];
double ub[maxCloRow];
char *colname[maxCloRow];
int rmatbeg[maxCloRow];
int *rmatind;
double *rmatval;
double rhs[maxCloRow];
char sense[maxCloRow];
char *rowname[maxCloRow];
double cofficient1, cofficient2, cofficient3;
int totalControlCount=0;


rmatind=(int *)malloc(maxCloRow*3000*sizeof(int));
rmatval=(double *)malloc(maxCloRow*3000*sizeof(double));


char *strings;

CPXENVptr env = NULL;
CPXLPptr lp = NULL;
env = CPXopenCPLEX (&status);
status = CPXsetintparam (env, CPXPARAM_ScreenOutput, CPX_ON);
status = CPXsetintparam (env, CPXPARAM_Read_DataCheck, CPX_ON);
lp = CPXcreateprob (env, &status, "lpex1");

//printf("L62 %i\n", numberGene);

status = CPXchgobjsen (env, lp, CPX_MIN);

//printf("L66\n");

for (int w=0; w<numberGene; w++)
{
strings=(char *)malloc(100*sizeof(char));
//printf("W %i\n", w);
cofficient1=0;
cofficient2=0;
cofficient3=0;
totalControlCount=0;

	for (int cases=0; cases<numCase; cases++)
	{
		cofficient1=cofficient1+fabs(listSamples[casesToConsider[cases]].values[w]-listSamples[centerId].values[w]);	
	}

	for (int samples=0; samples<totalSamples; samples++)
	{
		if ((listSamples[samples].caseControl==0) && (listSamples[samples].trainOrTest==0) && (listSamples[samples].filterSample==0))
		{
			for (int cases=0; cases<numCase; cases++)
			{
				cofficient3=cofficient3+fabs(listSamples[casesToConsider[cases]].values[w]-listSamples[samples].values[w]);
				totalControlCount++;
			}
				cofficient3=cofficient3+fabs(listSamples[samples].values[w]-listSamples[centerId].values[w]);
				totalControlCount++;
		}
	}




obj[w]=((double)cofficient1/(double)numCase)-((double)cofficient3/(double)totalControlCount);
//obj[w]=((double)cofficient1)-((double)cofficient3);


lb[w]=0;
ub[w]=1;

//printf("L86\n");
sprintf(strings, "w%i\0", w);
colname[w]=strings;

}


//printf("L87\n");

obj[numberGene]=0;
lb[numberGene]=0;
ub[numberGene]=CPX_INFBOUND;
//ub[numberGene]=1;
colname[numberGene]="r\0";
status = CPXnewcols (env, lp, numberGene+1, obj, lb, ub, NULL, colname);


int rmatbegCount=0;
int constraintCount=0;


for (int cases=0; cases<numCase; cases++)
{
	rmatbeg[constraintCount]=rmatbegCount;
	strings=(char *)malloc(100*sizeof(char));
	sprintf(strings, "c%i\0", constraintCount);
	rowname[constraintCount]=strings;	
	for (int w=0; w<numberGene; w++)
	{
		rmatind[rmatbegCount]=w;
		rmatval[rmatbegCount]=fabs(listSamples[casesToConsider[cases]].values[w]-listSamples[centerId].values[w]);
		rmatbegCount++;
	} 
	
		rmatind[rmatbegCount]=numberGene;
		rmatval[rmatbegCount]=-1;
		rmatbegCount++;

	sense[constraintCount]='L';
	rhs[constraintCount]=0;


	constraintCount++;
}

bool filter;
for (int allSamples=0; allSamples<totalSamples; allSamples++)
{
filter=false;
	if ((listSamples[allSamples].caseControl==0) && (listSamples[allSamples].trainOrTest==0))// it is a controls
	{
	if (filter==false)
	{
		//printf("ConstCount %i %i\n", constraintCount, rmatbegCount); 
		rmatbeg[constraintCount]=rmatbegCount;
		strings=(char *)malloc(100*sizeof(char));
		sprintf(strings, "d%i\0", constraintCount);
		
		rowname[constraintCount]=strings;
		for (int w=0; w<numberGene; w++)
		{
			rmatind[rmatbegCount]=w;
			rmatval[rmatbegCount]=-1*fabs(listSamples[allSamples].values[w]-listSamples[centerId].values[w]);
			rmatbegCount++;
		}

		rmatind[rmatbegCount]=numberGene;
		rmatval[rmatbegCount]=1;
		rmatbegCount++;

		sense[constraintCount]='L';
		rhs[constraintCount]=0;  
		constraintCount++;

		

	}


	}
} 


rmatbeg[constraintCount]=rmatbegCount;
rowname[constraintCount]="t";
for (int w=0; w<numberGene; w++)
{
	rmatind[rmatbegCount]=w;
	rmatval[rmatbegCount]=1;
	rmatbegCount++;
} 
sense[constraintCount]='L';
rhs[constraintCount]=5;
constraintCount++; 


status = CPXaddrows (env, lp, 0, constraintCount, rmatbegCount, rhs, sense, rmatbeg, rmatind, rmatval, NULL, rowname);
//status = populatebyrow (env, lp); 






status = CPXlpopt (env, lp); 

cur_numrows = CPXgetnumrows (env, lp);
cur_numcols = CPXgetnumcols (env, lp);

x = (double *) malloc (cur_numcols * sizeof(double));
slack = (double *) malloc (cur_numrows * sizeof(double));
dj = (double *) malloc (cur_numcols * sizeof(double));
pi = (double *) malloc (cur_numrows * sizeof(double));
status = CPXsolution (env, lp, &solstat, &objval, x, pi, slack, dj);


for (int count=0; count<cur_numcols; count++)
{
if (x[count]>0)
	printf("Col %i %lf\n", count, x[count]);

weightsCalcultaed[count]=x[count];
}

//status = CPXwriteprob (env, lp, "lpex1.lp", NULL);

free(rmatind);
free(rmatval);
free(x);
free(slack);
free(dj);
free(pi);
CPXfreeprob(env, &lp);
 CPXcloseCPLEX (&env);
}








int compare (const void * a, const void * b)
{

double xx = *(double*)a; 
double yy = *(double*)b;
if (xx<yy) return -1;
if (xx>yy) return 1;
return 0;

}


double calDist(double *geneWeights)
{

	double radisusControls[3000];
	double returnValue=0;
	double distance=0;
	int numControls=0;
	double radius;
	int countDistanceBigZero=0;

	numberCasesInSphere=0;
	numberControlInSphere=0;



	for (int count=0; count<3000; count++)
	{
		radisusControls[count]=0;
	}



	for (int count=0; count<totalSamples; count++)
	{
	if (listSamples[count].caseControl==0)
		numControls++;
		for (int count2=count+1; count2<totalSamples; count2++)
		{
		
		distance=0;
			for (int numGene=0; numGene<numberGene; numGene++)
			{
				distance=distance+geneWeights[numGene]*fabs(listSamples[count].values[numGene]-listSamples[count2].values[numGene]);
			}
		
		if (distance>0)
			countDistanceBigZero++;


		distanceMatrix[count][count2]=distance;
		distanceMatrix[count2][count]=distance;
		}
	casesInSphereCenter[count]=0;
	}

	int maxCaseInCenter=0;


	int casePickedAsCenter=-1;

	for (int count=0; count<totalSamples; count++)
	{
		numControls=0;
		if ((listSamples[count].caseControl==1) && (listSamples[count].filterSample==0) && (listSamples[count].trainOrTest==0))
		{
			for (int tempC=0; tempC<3000; tempC++)
			{
				radisusControls[tempC]=0;
			}

			for (int count2=0; count2<totalSamples; count2++)
			{
				if ((listSamples[count2].caseControl==0) && (listSamples[count2].filterSample==0) && (listSamples[count2].trainOrTest==0))
				{
					radisusControls[numControls]=distanceMatrix[count][count2];
					numControls++;
				}	
			}
			qsort(radisusControls, numControls, sizeof(double), compare);
			radiusOfSphere[count]=radisusControls[0];
			for (int count2=0; count2<totalSamples; count2++)
			{
				if ((listSamples[count2].caseControl==1) && (distanceMatrix[count][count2]<radiusOfSphere[count]) && (listSamples[count2].filterSample==0) && (listSamples[count2].trainOrTest==0))
				{
					casesInSphereCenter[count]++;
				}
			}
			if (maxCaseInCenter<casesInSphereCenter[count])
			{
				casePickedAsCenter=count;
				maxCaseInCenter=casesInSphereCenter[count];
				printf("Max cases %i\n", casesInSphereCenter[count]);
			}		
		}
	}
	

	bool notCovered=true;

	//for (int count=0; count<totalSamples; count++)
	//{
		if ((casePickedAsCenter>-1) && (casesInSphereCenter[casePickedAsCenter]==maxCaseInCenter) && (notCovered==true) && (listSamples[casePickedAsCenter].caseControl==1) && (listSamples[casePickedAsCenter].trainOrTest==0))
		{
			notCovered=false;
			centerId=casePickedAsCenter;
			printf("%s\n", listSamples[casePickedAsCenter].sampleName);
			radius=radiusOfSphere[casePickedAsCenter];
			for (int count2=0; count2<totalSamples; count2++)
			{
				if ((distanceMatrix[casePickedAsCenter][count2]<radiusOfSphere[casePickedAsCenter]) && (listSamples[count2].filterSample==0) && (listSamples[count2].trainOrTest==0))
				{
					if (listSamples[count2].caseControl==0)
					{
						controlInSphereId[numberControlInSphere]=count2;
						numberControlInSphere++;
					}else if (listSamples[count2].caseControl==1)
					{
						casesInSphereId[numberCasesInSphere]=count2;
						numberCasesInSphere++;
					}
					printf("%s\t%lf\n", listSamples[count2].sampleName, distanceMatrix[casePickedAsCenter][count2]);
				}
			}
		}
	//}

printf("%i %i %lf %i\n", numberCasesInSphere, numberControlInSphere, radius, countDistanceBigZero);
//calculateWeights(centerId, casesInSphereId, controlInSphereId, numberCasesInSphere, numberControlInSphere);
return radius;


}



int main(int argv, char *argc[])
{
	FILE *fp=fopen(argc[1],"r");
	numberGene=atoi(argc[2]);
	totalSamples=atoi(argc[3]);
	FILE *fp2=fopen(argc[4],"r");// Every samples name grouped as training (0) or testing sample (1)
	FILE *fp4=fopen(argc[5],"r");//Samples to remove 
	fpOut = fopen(argc[6],"w");


	int sampleId;
	int trainOrTest;
	char sampleName[50];
	double radius;
	for (int count=0; count<totalSamples; count++)
	{
		fscanf(fp, "%s\t%i", listSamples[count].sampleName, &(listSamples[count].caseControl));
		listSamples[count].filterSample=0;
		for (int count2=0; count2<numberGene; count2++)
		{
			fscanf(fp,"\t%lf", &(listSamples[count].values[count2]));
			//printf("%lf\n", listSamples[count].values[count2]);
		}
		fscanf(fp,"\n");

	}


	while(fscanf(fp2,"%s\t%i\n", sampleName, &trainOrTest)!=EOF)
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
			if (strcmp(listSamples[count].sampleName, sampleName)==0)
			{
				listSamples[count].filterSample=1;
			}
		}
	}



	for (int i=0; i<numberGene; i++)
	{
		weightsCalcultaed[i]=1;
	}


double *tempWeights;
tempWeights=(double *)malloc(sizeof(double)*numberGene);


fprintf(fpOut, "1\n"); 
int countX=0;
for (int numCenters=0; numCenters<1; numCenters++)
{

	while(numberCasesInSphere>=listSpheres[numCenters].numCases && countX<10) 
	{
	printf("CENTER NUM %i\n", numCenters);
	printf("CENTER NUM %i\n", numCenters);
	printf("CENTER NUM %i\n", numCenters);
	printf("CENTER NUM %i\n", numCenters);
	radius=calDist(weightsCalcultaed);

		

		      for (int count2=0; count2<numberCasesInSphere; count2++)
        		{
               			 printf("Samples Being Covered %i %s\n", countX, listSamples[casesInSphereId[count2]].sampleName);
        		}



		if (numberCasesInSphere>0)
		{

	
			if ((numberCasesInSphere>=listSpheres[numCenters].numCases) && (countX>0))
			{
				listSpheres[numCenters].radius=radius;
				listSpheres[numCenters].centerId=centerId;
				listSpheres[numCenters].numCases=numberCasesInSphere;
				for (int count2=0; count2<numberGene; count2++)
				{
					listSpheres[numCenters].weights[count2]=weightsCalcultaed[count2];
				}
			}
			

			/*for (int count=0; count<numberGene; count++)
			{
				fprintf(fpOut,"%lf ", weightsCalcultaed[count]);
			}
				fprintf(fpOut,"%s %lf\n", listSamples[centerId].sampleName, radius);

			*/

			calculateWeights(centerId, casesInSphereId, controlInSphereId, numberCasesInSphere, numberControlInSphere);
			



			for (int count2=0; count2<numberGene; count2++)
			{
				tempWeights[count2]=weightsCalcultaed[count2];
				if (tempWeights[count2]>0)
					printf("Temp Weights %lf\n", tempWeights[count2]);
			}	

			qsort(tempWeights, numberGene, sizeof(double), compare);

			for (int count2=0; count2<numberGene; count2++)
			{
				if (weightsCalcultaed[count2]<=tempWeights[numberGene-10])
				{
					weightsCalcultaed[count2]=0;
				}else 
					printf("Non Zero Weights %lf %lf\n", weightsCalcultaed[count2], tempWeights[10]);
			}
		};
		countX++;
	}


	for (int count=0; count<numberGene; count++)
	{
		fprintf(fpOut,"%lf ", listSpheres[numCenters].weights[count]);
	}
		fprintf(fpOut,"%s %lf\n", listSamples[listSpheres[numCenters].centerId].sampleName, listSpheres[numCenters].radius);



	for (int count2=0; count2<numberGene; count2++)
	{
		weightsCalcultaed[count2]=listSpheres[numCenters].weights[count2];
	}
	radius=calDist(weightsCalcultaed);



	listSamples[listSpheres[numCenters].centerId].filterSample=1;
	
	for (int count2=0; count2<numberCasesInSphere; count2++)
	{	
		printf("Samples Being Covered %s\n", listSamples[casesInSphereId[count2]].sampleName);
		listSamples[casesInSphereId[count2]].filterSample=1;
	}
	
	for (int count2=0; count2<numberGene; count2++)
	{
		weightsCalcultaed[count2]=1;
	}

}
fclose(fpOut);
}
