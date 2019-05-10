#include "lattice2d.h"
#include <math.h>

void aveQBilinear(const status *__restrict__ pstatus, const latticeStruct *__restrict__ lattice, double *__restrict__ radicalDensity, double *__restrict__ count)
{
	//Local variables
	const unsigned short* DensityProfile = lattice->lattice;
	const int LatticeDim = lattice->LatticeDim;
	const int cqsize = LatticeDim / 2; //qsize is the number of different q value samples.
	const double cqmax = cqsize - 1.0;

	radicalDensity[0] = DensityProfile[(LatticeDim + 1)*LatticeDim / 2];
	count[0] = 0;
	for (int i = 1; i < cqsize; ++i)
	{
		count[i] = 0;
		radicalDensity[i] = 0;
	}

	//We go though every cell (and its reflection about y axis simultaneously) in one quadrant of q space, find all the interpolation q value, calculate the integration of interpolation function over the arc and the lenght of arc.

	for (int itercol = 1; itercol < cqsize; ++itercol)
	{
		//The y index of the 4 vertices of the cell
		const int ky1 = itercol;
		const int ky2 = itercol + 1;
		const int refky1 = LatticeDim - ky1;
		const int refky2 = LatticeDim - ky2;
		const int dy1 = cqsize - ky2;
		const int dy2 = cqsize - ky1;
		const int dy12 = dy1*dy1;
		const int dy22 = dy2*dy2;
		for (int iterrow = 1; iterrow < cqsize; ++iterrow)
		{
			//The x index of the 4 vertices of the cell and its reflection about y axis.
			const int kx1 = iterrow;
			const int kx2 = iterrow + 1;
			const int refkx1 = LatticeDim - kx1;
			const int refkx2 = LatticeDim - kx2;
			const int dx1 = cqsize - kx2;
			const int dx2 = cqsize - kx1;
			const int dx12 = dx1*dx1;
			const int dx22 = dx2*dx2;

			//The distance between the upper-right corner and center
			const float dist1 = sqrtf(dx12 + dy12);
			//The distance between the lower-left corner and center
			const float dist2 = sqrtf(dx22 + dy22);

			//Calculate the range of interpolation points in this cell
			const int maxqidx = (int)ceilf(dist2);
			const int minqidx = (int)ceilf(dist1);

			//The distance between center and the other two vertices, will be used later
			const int dist3 = (dx22 + dy12);
			const int dist4 = (dx12 + dy22);

			//For each interpolation point.
			for (int iterq = minqidx; iterq < maxqidx; ++iterq)
			{
				if (iterq >= cqsize)
				{
					continue;
				}
				if (iterq == 0)
				{
					continue;
				}
				//To store the coordinates of the two ends of the arc in the cell
				double px[2];
				double py[2];

				//Current interpolation point, in the unit of qstep
				//bool ifdist3 = (dist3 >= iterq);
				//bool ifdist4 = (dist4 >= iterq);
				//px[0] = ifdist3*sqrt(iterq*iterq - dy1 * dy1) + dx2*(!ifdist3);
				//px[1] = dx1*ifdist4 + (!ifdist4)*sqrt(iterq*iterq - dy2 * dy2);
				//py[0] = (!ifdist3)*sqrt(iterq*iterq - dx2 * dx2) + dy1*ifdist3;
				//py[1] = ifdist4*sqrt(iterq*iterq - dx1 * dx1)+ dy2*(!ifdist4);


				//Calculate the coordinates of the ends of the arc
				const int q2 = iterq*iterq;
				if (dist3 >= q2)
				{
					py[0] = dy1;
					px[0] = sqrt(q2 - dy1 * dy1);
				}
				else
				{
					px[0] = dx2;
					py[0] = sqrt(q2 - dx2 * dx2);
				}
				
				if (dist4 >= q2)
				{
					px[1] = dx1;
					py[1] = sqrt(q2 - dx1 * dx1);
				}
				else
				{
					py[1] = dy2;
					px[1] = sqrt(q2 - dy2 * dy2);
				}

				//Calculate the angle that correspond to the arc, by solving the triangle
				const double inviterq = 1.0 / iterq;
				const double dist = (px[0] - px[1])*(px[0] - px[1]) + (py[0] - py[1])*(py[0] - py[1]);
				const double cosdt = 1.0 - 0.5 * dist * inviterq * inviterq;
				const double dtheta = acos(cosdt);
				const double dy1dtheta = dy1*dtheta;
				const double dx1dtheta = dx1*dtheta;

				const double cost1 = px[0] * inviterq;
				const double cost2 = px[1] * inviterq;
				const double sint1 = py[0] * inviterq;
				const double sint2 = py[1] * inviterq;
				const double iterqdcost = iterq*(cost2 - cost1);
				const double iterqdsint = iterq*(sint2 - sint1);
				const double dcos2t = 2 * (cost2*cost2 - cost1*cost1);
				const double factoru12u11 = -(iterqdcost + dy1dtheta);
				const double factoru21u11 = iterqdsint - dx1dtheta;
				const double factor4t = dx1*(iterqdcost + dy1dtheta) - iterq*iterq*dcos2t *0.25 - dy1*iterqdsint;
				const double factoru11 = dtheta - factoru12u11 - factoru21u11 + factor4t;
				const double factoru12 = factoru12u11 - factor4t;
				const double factoru21 = factoru21u11 - factor4t;

				double u22 = DensityProfile[kx1 + LatticeDim * ky1];
				double u12 = DensityProfile[kx2 + LatticeDim * ky1];
				double u21 = DensityProfile[kx1 + LatticeDim * ky2];
				double u11 = DensityProfile[kx2 + LatticeDim * ky2];

				//The integration has been done analytically, this is just substitution.
				radicalDensity[iterq] += u11*factoru11 + factoru12*u12 + factoru21*u21 + factor4t*u22;

				u12 = DensityProfile[refkx2 + LatticeDim * ky1];
				u22 = DensityProfile[refkx1 + LatticeDim * ky1];
				u11 = DensityProfile[refkx2 + LatticeDim * ky2];
				u21 = DensityProfile[refkx1 + LatticeDim * ky2];

				radicalDensity[iterq] += u11*factoru11 + factoru12*u12 + factoru21*u21 + factor4t*u22;

				u21 = DensityProfile[kx1 + LatticeDim * refky2];
				u11 = DensityProfile[kx2 + LatticeDim * refky2];
				u22 = DensityProfile[kx1 + LatticeDim * refky1];
				u12 = DensityProfile[kx2 + LatticeDim * refky1];

				radicalDensity[iterq] += u11*factoru11 + factoru12*u12 + factoru21*u21 + factor4t*u22;

				u11 = DensityProfile[refkx2 + LatticeDim * refky2];
				u21 = DensityProfile[refkx1 + LatticeDim * refky2];
				u12 = DensityProfile[refkx2 + LatticeDim * refky1];
				u22 = DensityProfile[refkx1 + LatticeDim * refky1];

				radicalDensity[iterq] += u11*factoru11 + factoru12*u12 + factoru21*u21 + factor4t*u22;

				count[iterq] += dtheta * 4.0;
			}
		}
	}
	for (int i = 1; i < cqsize; ++i)
	{
		radicalDensity[i] /= count[i];
	}
}
