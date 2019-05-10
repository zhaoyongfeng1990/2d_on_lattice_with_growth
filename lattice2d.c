#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
//#include <omp.h>

#include "lattice2d.h"
#include "mt64.h"

int main(int argc, const char *argv[])
{
	clock_t it, ft;
	it = clock();

	status cstatus; //the structure for status of current system
	//circleRegion circR;
	movingList movingEcoli;
	tumblingList tumblingEcoli;
	hoppingList hoppingEcoli;
	latticeStruct lattice;
	histogram hist;
	setDefaultStatus(&cstatus, &movingEcoli, &tumblingEcoli, &hoppingEcoli /*, &circR*/, &hist, &lattice);
	uint64_t seed=cstatus.seed;
	init_genrand64(seed);
	
	//initializeOldSnapshotParticles(&cstatus, &movingEcoli, &circR, 224);

	//int* DensityProfile = NULL;// (int*)malloc(lattice.NumSite * sizeof(int));
	ecoli *ecoliList = (ecoli *)malloc(cstatus.MaximumNumOfCells * sizeof(ecoli));
	double *radicalDensity = (double *)malloc(lattice.LatticeDim / 2 * sizeof(double));
	double *count = (double *)malloc(lattice.LatticeDim / 2 * sizeof(double));

	putObstacles(&lattice, &cstatus);
	//putParticles(lattice, ecoliList, &cstatus);
	putParticlesGaussian(&lattice, &movingEcoli, &hoppingEcoli, ecoliList, &cstatus);
	//putParticlesLine(lattice, ecoliList, &cstatus);
	//readOldSnapshotParticles(lattice, ecoliList, &cstatus, 224);
	generateDensityProfile(&lattice, ecoliList, &cstatus /*, &circR*/, &hist);
	//checkIfBlocked(&cstatus, &lattice, ecoliList, &movingEcoli);

	//snapshot(&lattice, ecoliList, &cstatus, 0);
	double totalTime = cstatus.totalTime;
	//int fileIdx = 225;
	int fileIdx = 1;
	FILE *outfile = fopen("result.txt", "w");
	int step = 0;
	int timeIntv = 0;
	int HopEventTime = 0;
	int TumbEventTime = 0;
	int GrowthEventTime = 0;
	int DeathEventTime = 0;
	int RecoverEventTime = 0;
	int RotDiffEventTime = 0;
	int failure = 0;
	while (cstatus.time < totalTime)
	{
		//++step;
		int ifHopping = iterate(&lattice, &ecoliList, &cstatus, &movingEcoli, &tumblingEcoli, &hoppingEcoli /*, &circR*/, &hist);
		switch (ifHopping)
		{
		case 0:
			++TumbEventTime;
			break;
		case 1:
			++HopEventTime;
			break;
		case 2:
			++GrowthEventTime;
			break;
		case 3:
			++DeathEventTime;
			break;
		case 4:
			++RecoverEventTime;
			break;
		case 5:
			++RotDiffEventTime;
			break;
		default:
			++failure;
			break;
		}

		// //if(1==ifHopping)
		// //{
		//   //fprintf(outfile, "%f ", cstatus.time);
		//   //double sumx=0;
		//   //double sumy=0;
		//   // #pragma omp parallel for reduction(+: sumx, sumy)
		//   // for (int i = 0; i < cstatus.NumOfCells; ++i)
		//   // {
		//   //   int dx=ecoliList[i].deltaX;
		//   //   int dy=ecoliList[i].deltaY;
		//   //   sumx+=dx*dx;
		//   //   sumy+=dy*dy;
		//   //   //fprintf(outfile, "%d ", dx*dx+dy*dy);
		//   // }
		//   // sumx=sumx/cstatus.NumOfCells;
		//   // sumy=sumy/cstatus.NumOfCells;
		//   // fprintf(outfile, "%f %f \n", sumx, sumy);
		//   //++timeIntv;
		//   //if (1000==timeIntv)
		//   //{
		//   //  timeIntv=0;
		//   //}
		// //}
		if (cstatus.time > (double)fileIdx*2)
		{
			//aveQBilinear(&cstatus, DensityProfile, radicalDensity, count);

			aveQBilinear(&cstatus, &lattice, radicalDensity, count);

			/*for (int j = 0; j < cstatus.LatticeDim / 2; ++j)
			{
				radicalDensity[j] = 0;
			}
			for (int i = 0; i < cstatus.LatticeDim; ++i)
			{
				radicalDensity[0] += DensityProfile[i*cstatus.LatticeDim + cstatus.LatticeDim / 2];
				for (int j = 1; j < cstatus.LatticeDim/2; ++j)
				{
					radicalDensity[j] += DensityProfile[i*cstatus.LatticeDim + cstatus.LatticeDim / 2 + j];
					radicalDensity[j] += DensityProfile[i*cstatus.LatticeDim + cstatus.LatticeDim / 2 - j];
				}
			}*/
			for (int i = 0; i < lattice.LatticeDim / 2; ++i)
			{
				//fprintf(outfile, "%e ", radicalDensity[i]*circR.UnitDensity);
				fprintf(outfile, "%e ", radicalDensity[i]);
			}
			fprintf(outfile, "\n");

			//	snapshotParticles(&lattice, ecoliList, &cstatus, fileIdx);
			//snapshotLattices(&lattice, ecoliList, &cstatus, fileIdx);
			++fileIdx;
			printf("%d: %d %d %d %d %d %d %d\n",fileIdx, HopEventTime, TumbEventTime, GrowthEventTime, DeathEventTime, RecoverEventTime, RotDiffEventTime, failure);
		}
		//for (int i = 0; i < cstatus.NumOfCells; ++i)
		//printf("%f, %d, %d \n", cstatus.time, ecoliList[i].deltaX, ecoliList[i].deltaY);
	}

	//for (int i = 0; i < cstatus.NumOfCells; ++i)
	//{
	//	fprintf(outfile, "%d, %d \n", ecoliList[i].deltaX, ecoliList[i].deltaY);
	//}

	fclose(outfile);
	/////////////////////

	/////////////////////l
	destructeStatus(&movingEcoli, &tumblingEcoli, &hoppingEcoli, /* &circR,*/ &lattice, &hist);
	free(ecoliList);
	//free(DensityProfile);
	free(radicalDensity);
	free(count);

	ft = clock();
	printf("%f\n", (double)(ft - it) / CLOCKS_PER_SEC);
	return 0;
}
