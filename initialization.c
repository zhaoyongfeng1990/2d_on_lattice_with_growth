#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "lattice2d.h"
#include "mt64.h"


//uniformly sample a new direction of bacteria.
void sampleDirection(ecoli *__restrict__ cell);

//generate a lattices with obstruction. Value is the number of particles, and -1 denotes blocked site.
void putObstacles(latticeStruct *__restrict__ lattice, const status *__restrict__ pstatus)
{
	const double GelConcentration = pstatus->GelConcentration;
	const int NumSite = lattice->NumSite;

	for (int i = 0; i < NumSite; ++i)
	{
		lattice->pLattice[i] = NULL;
		lattice->whereInHist[i] = -1;
		lattice->lattice[i] = 0;
		lattice->obstacles[i] = (genrand64_real1() < GelConcentration);
	}
}

//uniformly sample a new direction of bacteria.
void sampleDirection(ecoli *__restrict__ cell)
{
	const int randomDirection = (int)(genrand64_int64() & 3); //a random number from 0 to 7
	switch (randomDirection)
	{
	case 0:
		cell->directionX = 1;
		cell->directionY = 0;
		break;
	case 1:
		cell->directionX = -1;
		cell->directionY = 0;
		break;
	case 2:
		cell->directionX = 0;
		cell->directionY = 1;
		break;
	case 3:
		cell->directionX = 0;
		cell->directionY = -1;
		break;
	}
}

//put particles on lattices.
void putParticles(latticeStruct *__restrict__ lattice, movingList *__restrict__ pmlist, hoppingList *__restrict__ phlist, ecoli *__restrict__ ecoliList, status *__restrict__ pstatus)
{
	const int NumOfCells = pstatus->NumOfCells;
	const int LatticeDim = lattice->LatticeDim;
	int NumOfMovingCellsEdge = 0;
	int NumOfHoppingCellsEdge = 0;

	for (int i = 0; i < NumOfCells; ++i)
	{
		int idx; //random coordinates and index in lattice array
		int randomX, randomY;

		randomX = (int)floor(LatticeDim * genrand64_real2());
		randomY = (int)floor(LatticeDim * genrand64_real2());
		idx = randomY * LatticeDim + randomX; //index in lattice array

		//Non-interacting particles: comment this line
		//lattice[idx]=i+1;   //the site is occupied
		//////
		if (lattice->pLattice[idx] == NULL)
		{
			lattice->pLattice[idx] = (int *)malloc((size_t)pstatus->capacity * 20 * sizeof(int));
		}
		lattice->pLattice[idx][lattice->lattice[idx]] = i;
		ecoliList[i].whereInpLattice = lattice->lattice[idx];
		++lattice->lattice[idx];
		ecoliList[i].deltaX = 0;
		ecoliList[i].deltaY = 0;

		pmlist->movingEcoliEdge[NumOfMovingCellsEdge] = i;
		ecoliList[i].whereInMove = NumOfMovingCellsEdge;
		++NumOfMovingCellsEdge;
		phlist->movingEcoliEdge[NumOfHoppingCellsEdge] = i;
		ecoliList[i].whereInHop = NumOfHoppingCellsEdge;
		++NumOfHoppingCellsEdge;

		ecoliList[i].posX = randomX;
		ecoliList[i].posY = randomY;
		sampleDirection(ecoliList + i);
		ecoliList[i].ifMoving = 1;
		ecoliList[i].ifTumbling = 0;
	}
	pmlist->NumOfMovingCellsEdge = NumOfMovingCellsEdge;
	phlist->NumOfMovingCellsEdge = NumOfHoppingCellsEdge;
}

//put all the particles at the origin, with random directions
void putParticlesGaussian(latticeStruct *__restrict__ lattice, movingList *__restrict__ pmlist, hoppingList *__restrict__ phlist, ecoli *__restrict__ ecoliList, status *__restrict__ pstatus)
{
	const int NumOfCells = pstatus->NumOfCells;
	const int LatticeDim = lattice->LatticeDim;
	int NumOfMovingCellsEdge = 0;
	int NumOfHoppingCellsEdge = 0;
	//const int idx = (LatticeDim + 1)*LatticeDim / 2;
	//lattice[idx] = NumOfCells;
	const double diffusivity = sqrt(0.5*pstatus->HoppingRate*pstatus->HoppingRate*(1+pstatus->RotDiff/pstatus->RecoverRate) / pstatus->TumblingRate/(1+pstatus->TumblingRate/pstatus->RecoverRate)/pstatus->RotDiff*pstatus->RecoverRate/pstatus->GrowthRate );
	const int Module=LatticeDim-1;
	for (int i = 0; i < NumOfCells; ++i)
	{
		int idx; //random coordinates and index in lattice array
		int randomX, randomY;

		const double rand1 = sqrt(-2 * log(genrand64_real3()));
		const double rand2 = genrand64_real3();
		const double nrand1 = rand1 * cos(2 * pi * rand2) * diffusivity;
		const double nrand2 = rand1 * sin(2 * pi * rand2) * diffusivity;
		randomX = ((int)floor(nrand1) + LatticeDim / 2)&Module;
		randomY = ((int)floor(nrand2) + LatticeDim / 2)&Module;
		idx = randomY * LatticeDim + randomX; //index in lattice array

		if (lattice->pLattice[idx] == NULL)
		{
			lattice->pLattice[idx] = (int *)malloc((size_t)pstatus->capacity * 20 * sizeof(int));
		}
		lattice->pLattice[idx][lattice->lattice[idx]] = i;
		ecoliList[i].whereInpLattice = lattice->lattice[idx];
		++lattice->lattice[idx];
		ecoliList[i].deltaX = 0;
		ecoliList[i].deltaY = 0;

		pmlist->movingEcoliEdge[NumOfMovingCellsEdge] = i;
		ecoliList[i].whereInMove = NumOfMovingCellsEdge;
		++NumOfMovingCellsEdge;
		phlist->movingEcoliEdge[NumOfHoppingCellsEdge] = i;
		ecoliList[i].whereInHop = NumOfHoppingCellsEdge;
		++NumOfHoppingCellsEdge;

		ecoliList[i].posX = randomX;
		ecoliList[i].posY = randomY;
		sampleDirection(ecoliList + i);
		ecoliList[i].ifMoving = 1;
		ecoliList[i].ifTumbling = 0;
	}
	pmlist->NumOfMovingCellsEdge = NumOfMovingCellsEdge;
	phlist->NumOfMovingCellsEdge = NumOfHoppingCellsEdge;
}

//initialize pstatus with the pre-setting values
void setDefaultStatus(status *__restrict__ pstatus, movingList *__restrict__ pmlist, tumblingList *__restrict__ ptlist, hoppingList *__restrict__ phlist, /* circleRegion *pcirc,*/ histogram *__restrict__ phist, latticeStruct *__restrict__ lattice)
{
	FILE *parafile = fopen("parameters.txt", "r");
	double setTotalTime;
	short setLatticeDim;
	uint64_t setMaximumNumOfCells;
	uint64_t setNumOfCells;
	//short setInteractionRange;
	double setGelConcentration;
	double setTumblingRate;
	double setHoppingRate;
	double setTumbDuration;
	double setRotDiff;
	double setGrowthRate;
	short setCapacity;
	
	fscanf(parafile, "%lf\n", &setTotalTime);
	fscanf(parafile, "%hd\n", &setLatticeDim);
	//fscanf(parafile, "%d\n", &setInteractionRange);
	fscanf(parafile, "%ld\n", &setMaximumNumOfCells);
	fscanf(parafile, "%ld\n", &setNumOfCells);
	fscanf(parafile, "%lf\n", &setGelConcentration);
	fscanf(parafile, "%lf\n", &setTumblingRate);
	fscanf(parafile, "%lf\n", &setHoppingRate);
	fscanf(parafile, "%lf\n", &setTumbDuration);
	fscanf(parafile, "%lf\n", &setRotDiff);
	fscanf(parafile, "%lf\n", &setGrowthRate);
	fscanf(parafile, "%hd\n", &setCapacity);
	fscanf(parafile, "%ld\n", &(pstatus->seed));

	fclose(parafile);
	setNumOfCells *= (int)(0.5*setHoppingRate*setHoppingRate*(1+setRotDiff*setTumbDuration) / setTumblingRate/(1+setTumblingRate*setTumbDuration)/setTumbDuration/setRotDiff / setGrowthRate *2*pi);
	int setNumSite = setLatticeDim * setLatticeDim;
	pstatus->time = 0;
	pstatus->totalTime = setTotalTime;
	pstatus->GelConcentration = setGelConcentration;
	pstatus->TumblingRate = setTumblingRate;
	pstatus->HoppingRate = setHoppingRate;
	//pstatus->interactionRange = setInteractionRange;
	pstatus->NumOfCells = setNumOfCells;
	pstatus->MaximumNumOfCells = setMaximumNumOfCells;
	pstatus->GrowthRate = setGrowthRate;
	pstatus->capacity = setCapacity;
	pstatus->totalDeathRate = 0;
	pstatus->RecoverRate = 1.0 / setTumbDuration;
	pstatus->RotDiff = setRotDiff;
	//pstatus->deathRate = (double*)malloc(setMaximumNumOfCells * sizeof(double));
	pmlist->movingEcoliEdge = (int *)malloc(setMaximumNumOfCells * sizeof(int));
	pmlist->NumOfMovingCellsEdge = 0;
	ptlist->movingEcoliEdge = (int *)malloc(setMaximumNumOfCells * sizeof(int));
	ptlist->NumOfMovingCellsEdge = 0;
	phlist->movingEcoliEdge = (int *)malloc(setMaximumNumOfCells * sizeof(int));
	phlist->NumOfMovingCellsEdge = 0;

	//for circular interaction region
	//pcirc->circleIdx = (int*)malloc((setInteractionRange + 1)*sizeof(int));
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
	//pcirc->UnitDensity = 1.0 / (double)setNumSitesInCircle;
	lattice->LatticeDim = setLatticeDim;
	lattice->NumSite = setNumSite;
	lattice->lattice = (unsigned short *)malloc(setNumSite * sizeof(unsigned short));
	lattice->obstacles = (bool *)malloc(setNumSite * sizeof(bool));
	//lattice->DensityProfile = (int *)malloc(setNumSite * sizeof(int));
	lattice->pLattice = (int **)malloc(setNumSite * sizeof(int *));
	lattice->whereInHist = (int *)malloc(setNumSite * sizeof(int));
	phist->maxNumOfbins = 20 * (int)setCapacity;
	phist->binCapacity = (int *)malloc(phist->maxNumOfbins * sizeof(int));
	phist->bins = (int **)malloc(phist->maxNumOfbins * sizeof(int *));
	phist->cellNum = (int *)malloc(phist->maxNumOfbins * sizeof(int));
	for (int i = 0; i < phist->maxNumOfbins; ++i)
	{
		phist->bins[i] = NULL;
		phist->cellNum[i] = 0;
		phist->binCapacity[i] = 0;
	}
}

void destructeStatus(movingList *__restrict__ pmlist, tumblingList *__restrict__ ptlist, hoppingList *__restrict__ phlist, /* circleRegion *pcirc,*/ latticeStruct *__restrict__ lattice, histogram *__restrict__ phist)
{
	free(pmlist->movingEcoliEdge);
	free(ptlist->movingEcoliEdge);
	free(phlist->movingEcoliEdge);
	free(phist->cellNum);
	for (int i = 0; i < phist->maxNumOfbins; ++i)
	{
		if (phist->bins[i] != NULL)
		{
			free(phist->bins[i]);
		}
	}
	free(phist->bins);
	free(lattice->lattice);
	//free(lattice->DensityProfile);
	for (int i = 0; i < lattice->NumSite; ++i)
	{
		if (lattice->pLattice[i] != NULL)
		{
			free(lattice->pLattice[i]);
		}
	}
	free(lattice->pLattice);
	free(lattice->whereInHist);
	//free(pcirc->circleIdx);
}

void checkHistBinCapacity(histogram *__restrict__ phist, const int binIdx)
{
	if (phist->bins[binIdx] == NULL)
	{
		phist->bins[binIdx] = (int *)malloc(100 * sizeof(int));
		phist->binCapacity[binIdx] = 100;
	}
	else if (phist->cellNum[binIdx] == phist->binCapacity[binIdx])
	{
		phist->binCapacity[binIdx] *= 2;
		int *temp = (int *)malloc(phist->binCapacity[binIdx] * sizeof(int));
		for (int j = 0; j < phist->cellNum[binIdx]; ++j)
		{
			temp[j] = phist->bins[binIdx][j];
		}
		free(phist->bins[binIdx]);
		phist->bins[binIdx] = temp;
	}
	if (binIdx > phist->largestNonZeroBin)
	{
		phist->largestNonZeroBin = binIdx;
	}
}