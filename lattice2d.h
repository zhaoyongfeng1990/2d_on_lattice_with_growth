//the default settings
// #define setLatticeDim 2048
// #define setDelta 1  //unit: um, useless for now.
// #define setNumOfCells 1000
// #define setGelConcentration 0.3
// #define setTumblingRate 0
// #define setHoppingRate 10

// static const int setLatticeDim=512*4;
// static const double setDelta=1;  //unit: um, useless for now.
// static const int setNumOfCells=1000;
// static const double setGelConcentration=0.3;
// static const double setTumblingRate=10;
// static const double setHoppingRate=10;
//
// //two constant for speed
// static const int setNumSite=setLatticeDim*setLatticeDim;
// static const double setdiagHoppingRate=setHoppingRate/1.414213562373095;
#include <stdint.h>
#include <stdbool.h>

static const double pi = 3.14159265358979323846264338327950288419716939937510;

typedef struct
{
	int deltaX; //displacement
	int deltaY;
	int whereInMove; //the index of this cell in movingEcoli array
	int whereInTumb;
	int whereInHop;
	int whereInpLattice;
	short posX; //position
	short posY;
	signed char directionX; //direction
	signed char directionY;
	signed char lastDirectCode;
	bool ifMoving; //bool value, true if cell is moving
	bool ifTumbling;
} ecoli;

typedef struct
{
	double time;
	double totalTime;
	double GelConcentration; //volume concentration of the gel
	double TumblingRate;
	double HoppingRate;
	double RecoverRate;
	double RotDiff;
	double GrowthRate;
	uint64_t totalDeathRate;
	uint64_t seed;
	uint64_t NumOfCells;			//number of cells
	uint64_t MaximumNumOfCells;  //maximum number of cells
	//short interactionRange; //defines the range of local density
	short capacity;			//environment capacity
							//double* deathRate;		  //list of death rates
} status;					//status and parameters of the system

/*
typedef struct
{
	short* circleIdx;			  //a list of the boundary of a circle, to avoid multiple calculations
	double UnitDensity;		  //the contribution to density of a single particle
} circleRegion;
*/

typedef struct
{
	int *movingEcoliEdge;	 //list of index of cells moving along edge
	int NumOfMovingCellsEdge; //number of cells moving along edge
} movingList;

typedef movingList tumblingList;

typedef movingList hoppingList;

typedef struct
{
	bool *obstacles;
	unsigned short *lattice;
	//int* DensityProfile;
	int **pLattice;
	int *whereInHist;
	int NumSite;	  //total number of sites
	short LatticeDim; //number of sites in a row
} latticeStruct;

typedef struct
{
	int **bins;
	int *cellNum;
	int *binCapacity;
	int maxNumOfbins;
	int largestNonZeroBin;
} histogram;

//generate a lattices with obstruction. Value is the number of particles, and -1 denotes blocked site.
void putObstacles(latticeStruct *__restrict__ lattice, const status *__restrict__ pstatus);

//put particles on lattices. lattice[idx] denotes the number of particles.
void putParticles(latticeStruct *__restrict__ lattice, movingList *__restrict__ pmlist, hoppingList *__restrict__ phlist, ecoli *__restrict__ ecoliList, status *__restrict__ pstatus);

//put all the particles at the origin, with random directions
void putParticlesGaussian(latticeStruct *__restrict__ lattice, movingList *__restrict__ pmlist, hoppingList *__restrict__ phlist, ecoli *__restrict__ ecoliList, status *__restrict__ pstatus);

//generate movingEcoliEdge
//void checkIfBlocked(const status *__restrict__ pstatus, const latticeStruct *__restrict__ lattice, ecoli *__restrict__ ecoliList, movingList *__restrict__ pmlist);

//modifying movingEcoliEdge after one event
//void checkIfBlockedSingleCell(const status *__restrict__ pstatus, const latticeStruct *__restrict__ lattice, ecoli *__restrict__ ecoliList, const int idx, movingList *__restrict__ pmlist);

//initialize pstatus with the pre-setting values
void setDefaultStatus(status *__restrict__ pstatus, movingList *__restrict__ pmlist, tumblingList *__restrict__ ptlist, hoppingList *__restrict__ phlist,/* circleRegion *pcirc,*/ histogram *__restrict__ phist, latticeStruct *__restrict__ lattice);

//destruction function
void destructeStatus(movingList *__restrict__ pmlist, tumblingList *__restrict__ ptlist, hoppingList *__restrict__ phlist, /* circleRegion *pcirc,*/ latticeStruct *__restrict__ lattice, histogram *__restrict__ phist);

//save a binary snapshot, containing all the information of the system at current time
void snapshot(const latticeStruct *__restrict__ lattice, const ecoli *__restrict__ ecoliList, const status *__restrict__ pstatus, const int fileIndex);

//save a binary snapshot for lattice only
void snapshotLattices(const latticeStruct *__restrict__ lattice, const ecoli *__restrict__ ecoliList, const status *__restrict__ pstatus, const int fileIndex);

//save a binary snapshot for particles only
void snapshotParticles(const latticeStruct *__restrict__ lattice, const ecoli *__restrict__ ecoliList, const status *__restrict__ pstatus, const int fileIndex);

//read the time and parameters of the system from a binary snapshot
void readSnapshotHead(latticeStruct *__restrict__ lattice, status *__restrict__ pstatus, const int fileIndex);

//read the data of lattices and cells from a binary snapshot
void readSnapshotData(latticeStruct *__restrict__ lattice, ecoli *__restrict__ ecoliList, const status *__restrict__ pstatus, const int fileIndex);

//simulate one gillespie step, return denote the type of event: 0:tumble, 1:hop, 2:blocked
int iterate(latticeStruct *__restrict__ lattice, ecoli **__restrict__ pecoliList, status *__restrict__ pstatus, movingList *__restrict__ pmlist, tumblingList *__restrict__ ptlist, hoppingList *__restrict__ phlist,/* const circleRegion *pcirc,*/ histogram *__restrict__ phist);

//let the No. cellIdx cell to divide, without blocking interaction, the new cell occupies the same site with its mother
void cellDivide(latticeStruct *__restrict__ lattice, ecoli **__restrict__ ecoliList, status *__restrict__ pstatus, const int cellIdx, movingList *__restrict__ pmlist, tumblingList *__restrict__ ptlist, hoppingList *__restrict__ phlist,/* const circleRegion *pcirc,*/ histogram *__restrict__ phist);

//let the No. cellIdx cell to die
void cellDeath(latticeStruct *__restrict__ lattice, ecoli *__restrict__ ecoliList, status *__restrict__ pstatus, const int cellIdx, movingList *__restrict__ pmlist, tumblingList *__restrict__ ptlist, hoppingList *__restrict__ phlist,/* const circleRegion *pcirc,*/ histogram *__restrict__ phist);

//calculate the coarse grained density profile, only need once at the beginning of the simulation
void generateDensityProfile(latticeStruct *__restrict__ lattice, const ecoli *__restrict__ ecoliList, status *__restrict__ pstatus,/* const circleRegion *pcirc,*/ histogram *__restrict__ phist);

//update the density profile after a cell at (posX, posY) hopping with direction (directionX, directionY), magic inside
void updateProfile_Hop(const latticeStruct *__restrict__ lattice, status *__restrict__ pstatus,/* const circleRegion *pcirc,*/ const int posX, const int posY, const int directionX, const int directionY);

void aveQBilinear(const status *__restrict__ pstatus, const latticeStruct *__restrict__ lattice, double *__restrict__ radicalDensity, double *__restrict__ count);
void checkHistBinCapacity(histogram *__restrict__ phist, const int binIdx);
void recoverStructures(latticeStruct *__restrict__ lattice, ecoli *__restrict__ ecoliList, status *__restrict__ pstatus, movingList *__restrict__ pmlist, int *__restrict__ DensityProfile,/* circleRegion *pcirc,*/ histogram *__restrict__ phist);
