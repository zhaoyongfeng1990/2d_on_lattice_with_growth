#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "lattice2d.h"
#include "mt64.h"

//update the density profile after a new cell appears at (posX, posY)
void updateProfile_add(const latticeStruct *__restrict__ lattice, status *__restrict__ pstatus, /* const circleRegion *pcirc,*/ const int posX, const int posY);

//update the density profile after a new cell disappears from (posX, posY)
void updateProfile_minus(const latticeStruct *__restrict__ lattice, status *__restrict__ pstatus, /* const circleRegion *pcirc,*/ const int posX, const int posY);

//let the No. cellIdx cell to divide, without blocking interaction, the new cell occupies the same site with its mother
void cellDivide(latticeStruct *__restrict__ lattice, ecoli **__restrict__ ecoliList, status *__restrict__ pstatus, const int cellIdx, movingList *__restrict__ pmlist, tumblingList *__restrict__ ptlist, hoppingList *__restrict__ phlist, /* const circleRegion *pcirc,*/ histogram *__restrict__ phist)
{
	ecoli *clist = *ecoliList;
	int cMaxNum = pstatus->MaximumNumOfCells;
	const int cCellNum = pstatus->NumOfCells;

	//If the array is not long enough, make a larger one.
	if (cCellNum == cMaxNum)
	{
		ecoli *newList = (ecoli *)malloc(2 * cMaxNum * sizeof(ecoli));
		for (int iter = 0; iter < cMaxNum; ++iter)
		{
			newList[iter] = clist[iter];
		}
		free(clist);
		cMaxNum *= 2;
		pstatus->MaximumNumOfCells = cMaxNum;
		*ecoliList = newList;
		clist = newList;
		int *newMovingCells = (int *)malloc(cMaxNum * sizeof(int));
		for (int i = 0; i < pmlist->NumOfMovingCellsEdge; ++i)
		{
			newMovingCells[i] = pmlist->movingEcoliEdge[i];
		}
		free(pmlist->movingEcoliEdge);
		pmlist->movingEcoliEdge = newMovingCells;

		int *newTumblingCells = (int *)malloc(cMaxNum * sizeof(int));
		for (int i = 0; i < ptlist->NumOfMovingCellsEdge; ++i)
		{
			newTumblingCells[i] = ptlist->movingEcoliEdge[i];
		}
		free(ptlist->movingEcoliEdge);
		ptlist->movingEcoliEdge = newTumblingCells;

		int *newHoppingCells = (int *)malloc(cMaxNum * sizeof(int));
		for (int i = 0; i < phlist->NumOfMovingCellsEdge; ++i)
		{
			newHoppingCells[i] = phlist->movingEcoliEdge[i];
		}
		free(phlist->movingEcoliEdge);
		phlist->movingEcoliEdge = newHoppingCells;
	}

	//Some codes in comments are for simulations with blocking interactions. In that case, the daughter
	//should appear at an empty site beside its mother. If there's no empty site, it will not be born.

	const int nextPosX = clist[cellIdx].posX;
	const int nextPosY = clist[cellIdx].posY;
	updateProfile_add(lattice, pstatus /*, pcirc*/, nextPosX, nextPosY); //update density profile

	const int LatticeDim = lattice->LatticeDim;
	const int idx = nextPosY * LatticeDim + nextPosX;

	if (lattice->pLattice[idx] == NULL)
	{
		lattice->pLattice[idx] = (int *)malloc((size_t)pstatus->capacity * 20 * sizeof(int));
	}
	lattice->pLattice[idx][lattice->lattice[idx]] = cCellNum;

	const int newNumber = lattice->lattice[idx];
	if (newNumber > 0)
	{
		const int oldNumber = newNumber - 1;
		--phist->cellNum[oldNumber];
		lattice->whereInHist[phist->bins[oldNumber][phist->cellNum[oldNumber]]] = lattice->whereInHist[idx];
		phist->bins[oldNumber][lattice->whereInHist[idx]] = phist->bins[oldNumber][phist->cellNum[oldNumber]];
	}
	checkHistBinCapacity(phist, newNumber);

	phist->bins[newNumber][phist->cellNum[newNumber]] = idx;
	lattice->whereInHist[idx] = phist->cellNum[newNumber];
	++phist->cellNum[newNumber];

	clist[cCellNum].deltaX = 0;
	clist[cCellNum].deltaY = 0;
	if (clist[cellIdx].ifMoving)
	{
		clist[cCellNum].whereInMove = pmlist->NumOfMovingCellsEdge;
		pmlist->movingEcoliEdge[pmlist->NumOfMovingCellsEdge] = cCellNum;
		++pmlist->NumOfMovingCellsEdge;
	}
	if (clist[cellIdx].ifTumbling)
	{
		clist[cCellNum].whereInTumb = ptlist->NumOfMovingCellsEdge;
		ptlist->movingEcoliEdge[ptlist->NumOfMovingCellsEdge] = cCellNum;
		++ptlist->NumOfMovingCellsEdge;
	}
	else
	{
		clist[cCellNum].whereInHop = phlist->NumOfMovingCellsEdge;
		phlist->movingEcoliEdge[phlist->NumOfMovingCellsEdge] = cCellNum;
		++phlist->NumOfMovingCellsEdge;
	}
	clist[cCellNum].whereInpLattice = lattice->lattice[idx];
	clist[cCellNum].posX = nextPosX;
	clist[cCellNum].posY = nextPosY;
	clist[cCellNum].directionX = clist[cellIdx].directionX;
	clist[cCellNum].directionY = clist[cellIdx].directionY;
	clist[cCellNum].lastDirectCode = clist[cellIdx].lastDirectCode;
	clist[cCellNum].ifMoving = clist[cellIdx].ifMoving;
	clist[cCellNum].ifTumbling = clist[cellIdx].ifTumbling;

	++lattice->lattice[idx];
	++pstatus->NumOfCells;
	return;
}

//let the No. cellIdx cell to die
void cellDeath(latticeStruct *__restrict__ lattice, ecoli *__restrict__ ecoliList, status *__restrict__ pstatus, const int cellIdx, movingList *__restrict__ pmlist, tumblingList *__restrict__ ptlist, hoppingList *__restrict__ phlist, /* const circleRegion *pcirc,*/ histogram *__restrict__ phist)
{
	const int posX = ecoliList[cellIdx].posX;
	const int posY = ecoliList[cellIdx].posY;
	updateProfile_minus(lattice, pstatus /*, pcirc*/, posX, posY); //update density profile

	const int LatticeDim = lattice->LatticeDim;
	const int idx = posY * LatticeDim + posX;

	const int oldNumber = lattice->lattice[idx] - 1;
	--phist->cellNum[oldNumber];
	lattice->whereInHist[phist->bins[oldNumber][phist->cellNum[oldNumber]]] = lattice->whereInHist[idx];
	phist->bins[oldNumber][lattice->whereInHist[idx]] = phist->bins[oldNumber][phist->cellNum[oldNumber]];

	if (oldNumber > 0)
	{
		const int newNumber = oldNumber - 1;
		checkHistBinCapacity(phist, newNumber);
		phist->bins[newNumber][phist->cellNum[newNumber]] = idx;
		lattice->whereInHist[idx] = phist->cellNum[newNumber];
		++phist->cellNum[newNumber];
	}

	--lattice->lattice[idx];
	int idxToExchange = lattice->pLattice[idx][lattice->lattice[idx]];
	lattice->pLattice[idx][ecoliList[cellIdx].whereInpLattice] = idxToExchange;
	ecoliList[idxToExchange].whereInpLattice = ecoliList[cellIdx].whereInpLattice;

	//modifying movingEcoliEdge array, we have three particles to care about:
	//the one at the end of ecoliList array, the one at the end of movingEcoliEdge array,
	//and the dead one.

	while (phist->cellNum[phist->largestNonZeroBin] == 0)
	{
		--phist->largestNonZeroBin;
	}

	if (ecoliList[cellIdx].ifMoving)
	{
		int idxToExchange = pmlist->movingEcoliEdge[pmlist->NumOfMovingCellsEdge - 1];
		pmlist->movingEcoliEdge[ecoliList[cellIdx].whereInMove] = idxToExchange;
		ecoliList[idxToExchange].whereInMove = ecoliList[cellIdx].whereInMove;
		--pmlist->NumOfMovingCellsEdge;
	}
	if (ecoliList[cellIdx].ifTumbling)
	{
		--ptlist->NumOfMovingCellsEdge;
		int idxToExchange = ptlist->movingEcoliEdge[ptlist->NumOfMovingCellsEdge];
		ptlist->movingEcoliEdge[ecoliList[cellIdx].whereInTumb] = idxToExchange;
		ecoliList[idxToExchange].whereInTumb = ecoliList[cellIdx].whereInTumb;
	}
	else
	{
		--phlist->NumOfMovingCellsEdge;
		int idxToExchange = phlist->movingEcoliEdge[phlist->NumOfMovingCellsEdge];
		phlist->movingEcoliEdge[ecoliList[cellIdx].whereInHop] = idxToExchange;
		ecoliList[idxToExchange].whereInHop = ecoliList[cellIdx].whereInHop;
	}
	//If cellIdx is the last cell in ecoliList, the exchange will destroy the data structure
	if (cellIdx != pstatus->NumOfCells - 1)
	{
		ecoliList[cellIdx] = ecoliList[pstatus->NumOfCells - 1];

		const int idx = ecoliList[cellIdx].posX + ecoliList[cellIdx].posY * LatticeDim;
		lattice->pLattice[idx][ecoliList[cellIdx].whereInpLattice] = cellIdx;

		if (ecoliList[cellIdx].ifMoving)
		{
			pmlist->movingEcoliEdge[ecoliList[cellIdx].whereInMove] = cellIdx;
		}
		if (ecoliList[cellIdx].ifTumbling)
		{
			ptlist->movingEcoliEdge[ecoliList[cellIdx].whereInTumb] = cellIdx;
		}
		else
		{
			phlist->movingEcoliEdge[ecoliList[cellIdx].whereInHop] = cellIdx;
		}
	}
	--pstatus->NumOfCells;
}

//calculate the coarse grained density profile, only need once at the beginning of the simulation
void generateDensityProfile(latticeStruct *__restrict__ lattice, const ecoli *__restrict__ ecoliList, status *__restrict__ pstatus, /* const circleRegion *pcirc,*/ histogram *__restrict__ phist)
{
	//const int NumSite = pstatus->NumSite;
	//for (int i = 0; i < NumSite; ++i)
	//{
	//	DensityProfile[i] = 0;
	//}

	//const int interactionRange = pstatus->interactionRange;
	const int LatticeDim = lattice->LatticeDim;
	//const int* circleIdx = pcirc->circleIdx;
	const int NumOfCells = pstatus->NumOfCells;
	//const double UnitDensity = pcirc->UnitDensity;
	////usually the initial number of particles is less than the number of sites
	////otherwise we can loop over the lattice sites
	//for (int i = 0; i < NumOfCells; i++)
	//{
	//	const int posX = ecoliList[i].posX;
	//	const int posY = ecoliList[i].posY;
	//	for (int x = posX - interactionRange; x <= posX + interactionRange; ++x)
	//	{
	//		int boundaryOffset = circleIdx[abs(x - posX)];
	//		for (int y = posY - boundaryOffset; y <= posY + boundaryOffset; ++y)
	//		{
	//			int coX = (x + LatticeDim) % LatticeDim;
	//			int coY = (y + LatticeDim) % LatticeDim;
	//			int idx = coY*LatticeDim + coX;
	//			++DensityProfile[idx];
	//		}
	//	}
	//}
	phist->largestNonZeroBin = -1;
	uint64_t totaldeathRate = 0;
	for (int i = 0; i < lattice->NumSite; ++i)
	{
		//totaldeathRate += DensityProfile[idx];
		if (lattice->lattice[i] > 0)
		{
			totaldeathRate += lattice->lattice[i] * lattice->lattice[i];
			int latticeidx = lattice->lattice[i] - 1;

			checkHistBinCapacity(phist, latticeidx);
			phist->bins[latticeidx][phist->cellNum[latticeidx]] = i;
			lattice->whereInHist[i] = phist->cellNum[latticeidx];
			++phist->cellNum[latticeidx];
		}
	}
	pstatus->totalDeathRate = totaldeathRate;
}

//update the density profile after a new cell appears at (posX, posY)
void updateProfile_add(const latticeStruct *__restrict__ lattice, status *__restrict__ pstatus, /* const circleRegion *pcirc,*/ const int posX, const int posY)
{
	//const int interactionRange = pstatus->interactionRange;
	const int LatticeDim = lattice->LatticeDim;
	//const int* circleIdx = pcirc->circleIdx;
	//const double UnitDensity = pcirc->UnitDensity;

	////only the sites in a small region need to be updated
	//int count = 0;
	//for (int y = posY - interactionRange; y <= posY + interactionRange; ++y)
	//{
	//	const int boundaryOffset = circleIdx[abs(y - posY)];
	//	const int coY = ((y + LatticeDim) % LatticeDim)*LatticeDim;
	//	for (int x = posX - boundaryOffset; x <= posX + boundaryOffset; ++x)
	//	{
	//		const int coX = (x + LatticeDim) % LatticeDim;
	//		const int idx = coY + coX;
	//		++DensityProfile[idx];
	//		if (lattice[idx]>0)
	//		{
	//			count+=lattice[idx];
	//		}
	//	}
	//}
	//const int idx = posY*LatticeDim + posX;
	//pstatus->totalDeathRate += (count + DensityProfile[idx]);

	//For local density
	const int idx = posY * LatticeDim + posX;
	pstatus->totalDeathRate += (2 * lattice->lattice[idx] + 1);
}

//update the density profile after a new cell disappears from (posX, posY)
void updateProfile_minus(const latticeStruct *__restrict__ lattice, status *__restrict__ pstatus, /* const circleRegion *pcirc,*/ const int posX, const int posY)
{
	//const int interactionRange = pstatus->interactionRange;
	const int LatticeDim = lattice->LatticeDim;
	//const int* circleIdx = pcirc->circleIdx;
	//const double UnitDensity = pcirc->UnitDensity;

	////only the sites in a small region need to be updated
	//int count = 0;
	//for (int y = posY - interactionRange; y <= posY + interactionRange; ++y)
	//{
	//	const int boundaryOffset = circleIdx[abs(y - posY)];
	//	const int coY = ((y + LatticeDim) % LatticeDim)*LatticeDim;
	//	for (int x = posX - boundaryOffset; x <= posX + boundaryOffset; ++x)
	//	{
	//		const int coX = (x + LatticeDim) % LatticeDim;
	//		const int idx = coY + coX;
	//		--DensityProfile[idx];
	//		if (lattice[idx]>0)
	//		{
	//			count+=lattice[idx];
	//		}
	//	}
	//}
	//const int idx = posY*LatticeDim + posX;
	//pstatus->totalDeathRate -= (count + DensityProfile[idx]);

	//For local density
	const int idx = posY * LatticeDim + posX;
	pstatus->totalDeathRate -= (2 * lattice->lattice[idx] - 1);
}

//update the density profile after a cell at (posX, posY) hopping with direction (directionX, directionY), magic inside
void updateProfile_Hop(const latticeStruct *__restrict__ lattice, status *__restrict__ pstatus, /* const circleRegion *pcirc,*/ const int posX, const int posY, const int directionX, const int directionY)
{
	//const int interactionRange = pstatus->interactionRange;
	const int LatticeDim = lattice->LatticeDim;
	const int Module = lattice->LatticeDim - 1;
	//const int* circleIdx = pcirc->circleIdx;
	//const double UnitDensity = pcirc->UnitDensity;

	//only the sites near the boundary need to be updated
	//take care of the direction
	//	int count = 0;
	//	if (directionY == 1)
	//	{
	//		for (int x = posX - interactionRange; x <= posX + interactionRange; ++x)
	//		{
	//			const int boundaryOffset = circleIdx[abs(x - posX)];
	//			const int y = posY - boundaryOffset;
	//
	//			const int coX = (x + LatticeDim) % LatticeDim;
	//			const int coY = (y + LatticeDim) % LatticeDim;
	//
	//			const int idx = coY*LatticeDim + coX;
	//
	//			--DensityProfile[idx];
	//
	//			if (lattice[idx] > 0)
	//			{
	//				count -= lattice[idx];
	//			}
	//
	//			const int y2 = posY + boundaryOffset + 1;
	//			const int coY2 = (y2 + LatticeDim) % LatticeDim;
	//			const int idx2 = coY2*LatticeDim + coX;
	//
	//			++DensityProfile[idx2];
	//
	//			if (lattice[idx2] > 0)
	//			{
	//				count += lattice[idx2];
	//			}
	//		}
	//	}
	//	else if (directionY == -1)
	//	{
	//		for (int x = posX - interactionRange; x <= posX + interactionRange; ++x)
	//		{
	//			const int boundaryOffset = circleIdx[abs(x - posX)];
	//			const int y = posY - boundaryOffset - 1;
	//
	//			const int coX = (x + LatticeDim) % LatticeDim;
	//			const int coY = (y + LatticeDim) % LatticeDim;
	//
	//			const int idx = coY*LatticeDim + coX;
	//
	//			++DensityProfile[idx];
	//
	//			if (lattice[idx] > 0)
	//			{
	//				count += lattice[idx];
	//			}
	//
	//			const int y2 = posY + boundaryOffset;
	//			const int coY2 = (y2 + LatticeDim) % LatticeDim;
	//			const int idx2 = coY2*LatticeDim + coX;
	//
	//			--DensityProfile[idx2];
	//
	//			if (lattice[idx2] > 0)
	//			{
	//				count -= lattice[idx2];
	//			}
	//		}
	//}
	//	else if (directionX == 1)
	//	{
	//		for (int y = posY - interactionRange; y <= posY + interactionRange; ++y)
	//		{
	//			const int boundaryOffset = circleIdx[abs(y - posY)];
	//			const int x = posX - boundaryOffset;
	//
	//			const int coX = (x + LatticeDim) % LatticeDim;
	//			const int coY = ((y + LatticeDim) % LatticeDim)*LatticeDim;
	//
	//			const int idx = coY + coX;
	//
	//			--DensityProfile[idx];
	//
	//			if (lattice[idx]>0)
	//			{
	//				count -= lattice[idx];
	//			}
	//
	//			const int x2 = posX + boundaryOffset + 1;
	//			const int coX2 = (x2 + LatticeDim) % LatticeDim;
	//			const int idx2 = coY + coX2;
	//
	//			++DensityProfile[idx2];
	//
	//			if (lattice[idx2]>0)
	//			{
	//				count += lattice[idx2];
	//			}
	//		}
	//	}
	//	else
	//	{
	//		for (int y = posY - interactionRange; y <= posY + interactionRange; ++y)
	//		{
	//			const int boundaryOffset = circleIdx[abs(y - posY)];
	//			const int x = posX - boundaryOffset-1;
	//
	//			const int coX = (x + LatticeDim) % LatticeDim;
	//			const int coY = ((y + LatticeDim) % LatticeDim)*LatticeDim;
	//
	//			const int idx = coY + coX;
	//
	//			++DensityProfile[idx];
	//
	//			if (lattice[idx]>0)
	//			{
	//				count += lattice[idx];
	//			}
	//
	//			const int x2 = posX + boundaryOffset;
	//			const int coX2 = (x2 + LatticeDim) % LatticeDim;
	//			const int idx2 = coY + coX2;
	//
	//			--DensityProfile[idx2];
	//
	//			if (lattice[idx2]>0)
	//			{
	//				count -= lattice[idx2];
	//			}
	//		}
	//	}
	//	const int index = posY*LatticeDim + posX;
	//	const int newX = (posX + directionX + LatticeDim) % LatticeDim;
	//	const int newY = (posY + directionY + LatticeDim) % LatticeDim;
	//	const int newIndex = newY*LatticeDim + newX;
	//	pstatus->totalDeathRate += (count + DensityProfile[newIndex] - DensityProfile[index]);

	//For local density
	const int index = posY * LatticeDim + posX;
	const int newX = (posX + directionX) & Module;
	const int newY = (posY + directionY) & Module;
	const int newIndex = newY * LatticeDim + newX;
	pstatus->totalDeathRate += 2 * (lattice->lattice[newIndex] - lattice->lattice[index] + 1);
}
