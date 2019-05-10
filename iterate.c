#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "lattice2d.h"
#include "mt64.h"
//#include <omp.h>

//simulate one gillespie step, return denote the type of event: 0:tumble, 1:hop, 2:growth, 3:death, -1:something is wrong!
int iterate(latticeStruct *__restrict__ lattice, ecoli **__restrict__ pecoliList, status *__restrict__ pstatus, movingList *__restrict__ pmlist, tumblingList *__restrict__ ptlist, hoppingList *__restrict__ phlist,/* const circleRegion *pcirc,*/ histogram *__restrict__ phist)
{
	ecoli *ecoliList = *pecoliList;
	//local constants
	const double GrowthRate = pstatus->GrowthRate;
	const double capacity = pstatus->capacity;
	const int NumOfCells = pstatus->NumOfCells;
	const int LatticeDim = lattice->LatticeDim;
	const int Module = lattice->LatticeDim - 1;
	const int NumSite = lattice->NumSite;
	//Calculate the death rate of each particles, and do summation in the same loop
	//const double UnitDensity = 1.0;// pcirc->UnitDensity;
	double totaldeathRate = pstatus->totalDeathRate * GrowthRate / capacity;

	//local constants again
	const double TumblingRate = pstatus->TumblingRate;
	const double HoppingRate = pstatus->HoppingRate;
	const double RecoverRate = pstatus->RecoverRate;
	const double RotDiff = pstatus->RotDiff;
	int NumOfMovingCellsEdge = pmlist->NumOfMovingCellsEdge;
	int NumOfTumblingCellsEdge = ptlist->NumOfMovingCellsEdge;
	int NumOfHoppingCellsEdge = phlist->NumOfMovingCellsEdge;
	const double hoppingThres = NumOfMovingCellsEdge * HoppingRate;					  //Threshold for particles to hop
	const double tumblingThres = NumOfHoppingCellsEdge * TumblingRate + hoppingThres; //Threshold for particles to tumble
	const double recoverThres = NumOfTumblingCellsEdge * RecoverRate + tumblingThres;
	const double rotdiffThres = NumOfTumblingCellsEdge * RotDiff + recoverThres;
	const double growthThres = NumOfCells * GrowthRate + rotdiffThres; //Threshold for particles to divide
	const double totalProb = totaldeathRate + growthThres;

	const double r1 = genrand64_real3();
	const double tau = -log(r1) / totalProb; //time increment
	pstatus->time += tau;
	//finding the event which is going to happen
	double randomNum = genrand64_real3() * totalProb;

	int *movingEcoliEdge = pmlist->movingEcoliEdge;
	int *tumblingEcoliEdge = ptlist->movingEcoliEdge;
	int *hoppingEcoliEdge = phlist->movingEcoliEdge;
	//The conditions are ordered such that the most frequent types comes first
	if (randomNum < hoppingThres) //hopping event
	{
		const int idx = floor(randomNum / HoppingRate);
		const int i = movingEcoliEdge[idx];

		const int posX = ecoliList[i].posX;
		const int posY = ecoliList[i].posY;
		const char directionX = ecoliList[i].directionX;
		const char directionY = ecoliList[i].directionY;
		const int index = posY * LatticeDim + posX;
		updateProfile_Hop(lattice, pstatus/*, pcirc*/, posX, posY, directionX, directionY); //Update density profile

		const int oldNumber = lattice->lattice[index] - 1;
		--phist->cellNum[oldNumber];
		lattice->whereInHist[phist->bins[oldNumber][phist->cellNum[oldNumber]]] = lattice->whereInHist[index];
		phist->bins[oldNumber][lattice->whereInHist[index]] = phist->bins[oldNumber][phist->cellNum[oldNumber]];

		if (oldNumber > 0)
		{
			const int newNumber = oldNumber - 1;
			checkHistBinCapacity(phist, newNumber);
			phist->bins[newNumber][phist->cellNum[newNumber]] = index;
			lattice->whereInHist[index] = phist->cellNum[newNumber];
			++phist->cellNum[newNumber];
		}
		--lattice->lattice[index]; //the particle number in the old site reduced by 1

		ecoliList[lattice->pLattice[index][lattice->lattice[index]]].whereInpLattice = ecoliList[i].whereInpLattice;
		lattice->pLattice[index][ecoliList[i].whereInpLattice] = lattice->pLattice[index][lattice->lattice[index]];

		while (phist->cellNum[phist->largestNonZeroBin] == 0)
		{
			--phist->largestNonZeroBin;
		}

		//the mod operation is to make periodic boundary condition
		//plus with LatticeDim helps eliminate the possibility of getting a negative value
		const int newX = (posX + directionX) & Module;
		const int newY = (posY + directionY) & Module;
		const int newIndex = newY * LatticeDim + newX;

		if (lattice->pLattice[newIndex] == NULL)
		{
			lattice->pLattice[newIndex] = (int *)malloc((size_t)pstatus->capacity * 20 * sizeof(int));
		}
		lattice->pLattice[newIndex][lattice->lattice[newIndex]] = i;
		ecoliList[i].whereInpLattice = lattice->lattice[newIndex];

		const int newNumber2 = lattice->lattice[newIndex];
		if (newNumber2 > 0)
		{
			const int oldNumber2 = newNumber2 - 1;
			--phist->cellNum[oldNumber2];
			lattice->whereInHist[phist->bins[oldNumber2][phist->cellNum[oldNumber2]]] = lattice->whereInHist[newIndex];
			phist->bins[oldNumber2][lattice->whereInHist[newIndex]] = phist->bins[oldNumber2][phist->cellNum[oldNumber2]];
		}
;
		checkHistBinCapacity(phist, newNumber2);

		phist->bins[newNumber2][phist->cellNum[newNumber2]] = newIndex;
		lattice->whereInHist[newIndex] = phist->cellNum[newNumber2];
		++phist->cellNum[newNumber2];

		++lattice->lattice[newIndex]; //and the new site has a new particle

		ecoliList[i].posX = newX;
		ecoliList[i].posY = newY;
		ecoliList[i].deltaX += directionX;
		ecoliList[i].deltaY += directionY;

		if (lattice->obstacles[newIndex] == 1) //now cell cannot move, deleting from corresponding movingEcoli array
		{
			ecoliList[i].ifMoving = 0;
			const char dX = ecoliList[i].directionX;
			const char dY = ecoliList[i].directionY;
			ecoliList[i].lastDirectCode = 3 * dX + dY;
			--NumOfMovingCellsEdge;
			ecoliList[movingEcoliEdge[NumOfMovingCellsEdge]].whereInMove = ecoliList[i].whereInMove;
			movingEcoliEdge[ecoliList[i].whereInMove] = movingEcoliEdge[NumOfMovingCellsEdge];

			pmlist->NumOfMovingCellsEdge = NumOfMovingCellsEdge;
		}
		return 1;
	}
	else if (randomNum < tumblingThres) //tumbling event
	{
		const int idx = hoppingEcoliEdge[(int)floor((randomNum - hoppingThres) / TumblingRate)];
		ecoliList[idx].ifTumbling = 1;
		ecoliList[idx].whereInTumb = NumOfTumblingCellsEdge;
		tumblingEcoliEdge[NumOfTumblingCellsEdge] = idx;
		++NumOfTumblingCellsEdge;
		ptlist->NumOfMovingCellsEdge = NumOfTumblingCellsEdge;

		--NumOfHoppingCellsEdge;
		ecoliList[hoppingEcoliEdge[NumOfHoppingCellsEdge]].whereInHop = ecoliList[idx].whereInHop;
		hoppingEcoliEdge[ecoliList[idx].whereInHop] = hoppingEcoliEdge[NumOfHoppingCellsEdge];
		phlist->NumOfMovingCellsEdge = NumOfHoppingCellsEdge;

		if (ecoliList[idx].ifMoving)
		{
			ecoliList[idx].ifMoving=0;
			--NumOfMovingCellsEdge;
			ecoliList[movingEcoliEdge[NumOfMovingCellsEdge]].whereInMove = ecoliList[idx].whereInMove;
			movingEcoliEdge[ecoliList[idx].whereInMove] = movingEcoliEdge[NumOfMovingCellsEdge];
			pmlist->NumOfMovingCellsEdge = NumOfMovingCellsEdge;
		}
		return 0;
	}
	else if (randomNum < recoverThres)  //recovery event
	{
		const int idx = tumblingEcoliEdge[(int)floor((randomNum - tumblingThres) / RecoverRate)];

		ecoliList[idx].ifTumbling = 0;
		--NumOfTumblingCellsEdge;
		ecoliList[tumblingEcoliEdge[NumOfTumblingCellsEdge]].whereInTumb = ecoliList[idx].whereInTumb;
		tumblingEcoliEdge[ecoliList[idx].whereInTumb] = tumblingEcoliEdge[NumOfTumblingCellsEdge];
		ptlist->NumOfMovingCellsEdge = NumOfTumblingCellsEdge;

		ecoliList[idx].whereInHop = NumOfHoppingCellsEdge;
		hoppingEcoliEdge[NumOfHoppingCellsEdge] = idx;
		++NumOfHoppingCellsEdge;
		phlist->NumOfMovingCellsEdge = NumOfHoppingCellsEdge;
		
		const int posX = ecoliList[idx].posX;
		const int posY = ecoliList[idx].posY;
		ecoliList[idx].ifMoving=1;
		if (lattice->obstacles[posY * LatticeDim + posX] == 1)
		{
		  const int dX = ecoliList[idx].directionX;
		  const int dY = ecoliList[idx].directionY;
		  const int directionCode = 3 * dX + dY;
		  ecoliList[idx].ifMoving = (directionCode != ecoliList[idx].lastDirectCode);
		}
		if(ecoliList[idx].ifMoving)
		{
			ecoliList[idx].whereInMove = NumOfMovingCellsEdge;
			movingEcoliEdge[NumOfMovingCellsEdge] = idx;
			++NumOfMovingCellsEdge;
			pmlist->NumOfMovingCellsEdge = NumOfMovingCellsEdge;
		}
		return 4;
	}
	else if (randomNum < rotdiffThres)  //rotation event
	{
		const int idx = tumblingEcoliEdge[(int)floor((randomNum - recoverThres) / RotDiff)];
		char directChange = (char)((genrand64_int64() & 1) << 1) - 1;

		ecoliList[idx].directionX = (!ecoliList[idx].directionX)*directChange;
		ecoliList[idx].directionY = (!ecoliList[idx].directionY)*directChange;
		return 5;
	}
	else if (randomNum < growthThres)
	{
		const int idx = (int)floor((randomNum - rotdiffThres) / GrowthRate);
		cellDivide(lattice, pecoliList, pstatus, idx, pmlist, ptlist, phlist/*, pcirc*/, phist);
		return 2;
	}
	else
	{
		randomNum -= growthThres;
		randomNum *= (capacity/GrowthRate);

		int i = phist->largestNonZeroBin;
		while (phist->cellNum[i] * (i + 1) * (i + 1) < randomNum)
		{
			randomNum -= phist->cellNum[i] * (i + 1) * (i + 1);
			--i;
		}
		const int latticeIdx = (int)floor(randomNum / (i + 1) / (i + 1));
		const int latticePos = phist->bins[i][latticeIdx];
		randomNum -= latticeIdx * (i + 1) * (i + 1);
		const int cellIdx = (int)floor(randomNum / (i + 1));
		const int idx = lattice->pLattice[latticePos][cellIdx];

		//int idx=NumOfCells;
		//uint64_t partialProb = 0;
		//do
		//{
		//	--idx;
		//	//printf("test\n");
		//	const int posX = ecoliList[idx].posX;
		//	const int posY = ecoliList[idx].posY;
		//	const int idxLattice = posX + posY*LatticeDim;
		//	partialProb += lattice->lattice[idxLattice];
		//	//partialProb += DensityProfile[idxLattice];
		//	//partialProb += pstatus->deathRate[idx];
		//} while (partialProb < randomNum);
		//printf("test \n");
		cellDeath(lattice, ecoliList, pstatus, idx, pmlist, ptlist, phlist/*, pcirc*/, phist);
		//printf("Dead: %d\n",idx);
		return 3;
	}
	return -1;
}
