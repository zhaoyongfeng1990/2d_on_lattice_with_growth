#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "lattice2d.h"

//save a binary snapshot, containing all the information of the system at current time
void snapshot(const latticeStruct *__restrict__ lattice, const ecoli *__restrict__ ecoliList, const status *__restrict__ pstatus, const int fileIndex)
{
	char filename[10];  //file name will be like s0000.bin
	char buffer[5];
	sprintf(buffer, "%04d", fileIndex);
	strcpy(filename, "s");
	strcat(filename, buffer);
	strcat(filename, ".bin");
	FILE* pOutputFile = fopen(filename, "wb");

	fwrite(pstatus, sizeof(status), 1, pOutputFile);
	fwrite(&(lattice->LatticeDim), sizeof(short), 1, pOutputFile);
	fwrite(lattice->lattice, sizeof(short), lattice->NumSite, pOutputFile);
	fwrite(ecoliList, sizeof(ecoli), pstatus->NumOfCells, pOutputFile);
	fclose(pOutputFile);
}

void snapshotLattices(const latticeStruct *__restrict__ lattice, const ecoli *__restrict__ ecoliList, const status *__restrict__ pstatus, const int fileIndex)
{
	char filename[10];  //file name will be like s0000.bin
	char buffer[5];
	sprintf(buffer, "%04d", fileIndex);
	strcpy(filename, "l");
	strcat(filename, buffer);
	strcat(filename, ".bin");
	FILE* pOutputFile = fopen(filename, "wb");

	fwrite(pstatus, sizeof(status), 1, pOutputFile);
	fwrite(&(lattice->LatticeDim), sizeof(short), 1, pOutputFile);
	fwrite(lattice->lattice, sizeof(short), lattice->NumSite, pOutputFile);
	fclose(pOutputFile);
}

void snapshotParticles(const latticeStruct *__restrict__ lattice, const ecoli *__restrict__ ecoliList, const status *__restrict__ pstatus, const int fileIndex)
{
	char filename[10];  //file name will be like s0000.bin
	char buffer[5];
	sprintf(buffer, "%04d", fileIndex);
	strcpy(filename, "p");
	strcat(filename, buffer);
	strcat(filename, ".bin");
	FILE* pOutputFile = fopen(filename, "wb");

	fwrite(pstatus, sizeof(status), 1, pOutputFile);
	fwrite(&(lattice->LatticeDim), sizeof(short), 1, pOutputFile);
	fwrite(ecoliList, sizeof(ecoli), pstatus->NumOfCells, pOutputFile);
	fclose(pOutputFile);
}

//read the time and parameters of the system from a binary snapshot
void readSnapshotHead(latticeStruct *__restrict__ lattice, status *__restrict__ pstatus, const int fileIndex)
{
	char filename[10];
	char buffer[5];
	sprintf(buffer, "%04d", fileIndex);
	strcpy(filename, "s");
	strcat(filename, buffer);
	strcat(filename, ".bin");
	FILE* pInputFile = fopen(filename, "rb");

	fread(pstatus, sizeof(status), 1, pInputFile);
	fread(&(lattice->LatticeDim), sizeof(short), 1, pInputFile);

	fclose(pInputFile);
}

//read the data of lattices and cells from a binary snapshot
void readSnapshotData(latticeStruct *__restrict__ lattice, ecoli *__restrict__ ecoliList, const status *__restrict__ pstatus, const int fileIndex)
{
	char filename[10];
	char buffer[5];
	sprintf(buffer, "%04d", fileIndex);
	strcpy(filename, "s");
	strcat(filename, buffer);
	strcat(filename, ".bin");
	FILE* pInputFile = fopen(filename, "rb");

	fseek(pInputFile, sizeof(status), SEEK_SET);
	fread(&(lattice->LatticeDim), sizeof(short), 1, pInputFile);
	fread(lattice->lattice, sizeof(short), lattice ->NumSite, pInputFile);
	fread(ecoliList, sizeof(ecoli), pstatus->NumOfCells, pInputFile);

	fclose(pInputFile);
}

void recoverStructures(latticeStruct *__restrict__ lattice, ecoli *__restrict__ ecoliList, status *__restrict__ pstatus, movingList *__restrict__ pmlist, int *__restrict__ DensityProfile,/* circleRegion *pcirc,*/ histogram *__restrict__ phist)
{
	pmlist->NumOfMovingCellsEdge = 0;

	//for circular interaction region
	//pcirc->circleIdx = (int*)malloc((pstatus->interactionRange + 1) * sizeof(int));
	//int setNumSitesInCircle = 0;
	//pcirc->circleIdx[0] = setInteractionRange;
	//for (int i = 1; i <= setInteractionRange; ++i)
	//{
	//	for (int j = 1; j <= setInteractionRange; ++j)
	//	{
	//		if (i*i+j*j>setInteractionRange*setInteractionRange)
	//		{
	//			pcirc->circleIdx[i] = j-1;
	//			break;
	//		}
	//		++setNumSitesInCircle;
	//	}
	//}
	//setNumSitesInCircle *= 4;
	//setNumSitesInCircle += 4 * setInteractionRange + 1;
	//pcirc->UnitDensity = 1.0 / (double)setNumSitesInCircle;;
	lattice->NumSite = lattice->LatticeDim*lattice->LatticeDim;
	lattice->obstacles = (bool*)malloc(lattice->NumSite*sizeof(int));
	lattice->lattice = (unsigned short*)malloc(lattice->NumSite * sizeof(unsigned short));
	lattice->pLattice = (int**)malloc(lattice->NumSite * sizeof(int*));
	lattice->whereInHist = (int*)malloc(lattice->NumSite * sizeof(int));
	phist->maxNumOfbins = 10 * (int)pstatus->capacity;
	phist->binCapacity = (int*)malloc(phist->maxNumOfbins * sizeof(int));
	phist->bins = (int**)malloc(phist->maxNumOfbins * sizeof(int*));
	phist->cellNum = (int *)malloc(phist->maxNumOfbins * sizeof(int));
	for (int i = 0; i < phist->maxNumOfbins; ++i)
	{
		phist->bins[i] = NULL;
		phist->cellNum[i] = 0;
		phist->binCapacity = 0;
	}
	for (int i = 0; i < lattice->NumSite; ++i)
	{
		lattice->pLattice[i] = NULL;
		lattice->whereInHist[i] = -1;
	}
	for (int i = 0; i < pstatus->NumOfCells; ++i)
	{
		const int idx = ecoliList[i].posY*lattice->LatticeDim + ecoliList[i].posX;   //index in lattice array
		if (lattice->pLattice[idx] == NULL)
		{
			lattice->pLattice[idx] = (int*)malloc((size_t)pstatus->capacity * 20 * sizeof(int));
		}
		lattice->pLattice[idx][lattice->lattice[idx]] = i;
		ecoliList[i].whereInpLattice = lattice->lattice[idx];
		++lattice->lattice[idx];
	}
	generateDensityProfile(lattice, ecoliList, pstatus/*, pcirc*/, phist);
	//checkIfBlocked(pstatus, lattice, ecoliList, pmlist);
}