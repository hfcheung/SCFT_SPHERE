#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "struct.h"
#include "solver.h"
#include "calc.h"

#define DEL(k, i, n) del[(i) + scftInfo->size * (n) + scftInfo->nHist * scftInfo->size * (k)]
#define OUTS(k, i, n) outs[(i) + scftInfo->size * (n) + scftInfo->nHist * scftInfo->size * (k)]

#define U(n, m) up[(m - 1) + (nRec - 1) * (n - 1)]
#define V(n) vp[n - 1]
#define A(n) ap[n - 1]

#define PI 3.1415926535897932385

void removeAvgFld(SCFT_INFO *scftInfo, double *wA, double *wB)
{
    // Setting the spatial average of the fields to zero
    double sA = 0.0, sB = 0.0;

#pragma omp parallel for reduction(+ \
                                   : sA, sB)
    for (int idx = 0; idx < scftInfo->size; ++idx)
    {
        sA += scftInfo->wts[idx] * wA[idx];
        sB += scftInfo->wts[idx] * wB[idx];
    }

    sA /= 4 * PI;
    sB /= 4 * PI;

#pragma omp parallel for
    for (int idx = 0; idx < scftInfo->size; ++idx)
    {
        wA[idx] -= sA;
        wB[idx] -= sB;
    }
}

void calcFreeE(SCFT_INFO *scftInfo)
{
    scftInfo->incompMax = 0.0;

    scftInfo->freeS = -log(scftInfo->Q);

    double freeAB1 = 0.0, freeW1 = 0.0;

#pragma omp parallel for reduction(+ \
                                   : freeAB1, freeW1)
    for (int idx = 0; idx < scftInfo->size; ++idx)
    {
        freeAB1 += scftInfo->wts[idx] * scftInfo->hAB * scftInfo->phA[idx] * scftInfo->phB[idx];
        freeW1 += scftInfo->wts[idx] * (scftInfo->phA[idx] * scftInfo->wA[idx] +
                                        scftInfo->phB[idx] * scftInfo->wB[idx] +
                                        scftInfo->ksi[idx] * (1.0 - scftInfo->phA[idx] - scftInfo->phB[idx]));
    }

    for (int idx = 0; idx < scftInfo->size; ++idx)
    {
        double t = fabs(1 - scftInfo->phA[idx] - scftInfo->phB[idx]);
        if (t > scftInfo->incompMax)
        {
            scftInfo->incompMax = t;
        }
    }

    freeAB1 /= 4 * PI;
    freeW1 /= 4 * PI;

    scftInfo->freeAB = freeAB1;
    scftInfo->freeW = (-1.0) * freeW1;

    scftInfo->freeE = scftInfo->freeAB + scftInfo->freeW + scftInfo->freeS;
}

void updateFldsHist(SCFT_INFO *scftInfo, double *del, double *outs, double *wADiff, double *wBDiff,
                    double *wANew, double *wBNew)
{
    // Only update the flds stored at "posUpdate"

    // for (int idx = 0; idx < scftInfo->size; ++idx)
    // {
    //     DEL(0, idx, scftInfo->posUpdate) = wADiff[idx];
    //     DEL(1, idx, scftInfo->posUpdate) = wBDiff[idx];
    //     OUTS(0, idx, scftInfo->posUpdate) = wANew[idx];
    //     OUTS(1, idx, scftInfo->posUpdate) = wBNew[idx];
    // }

    // scftInfo->posUpdate = (scftInfo->posUpdate + 1) % scftInfo->nHist;

    for (int idxHist = scftInfo->nHist - 1; idxHist > 0; --idxHist)
    {
        for (int idx = 0; idx < scftInfo->size; ++idx)
        {
            DEL(0, idx, idxHist) = DEL(0, idx, idxHist - 1);
            DEL(1, idx, idxHist) = DEL(1, idx, idxHist - 1);
            OUTS(0, idx, idxHist) = OUTS(0, idx, idxHist - 1);
            OUTS(1, idx, idxHist) = OUTS(1, idx, idxHist - 1);
        }
    }

    for (int idx = 0; idx < scftInfo->size; ++idx)
    {
        DEL(0, idx, 0) = wADiff[idx];
        DEL(1, idx, 0) = wBDiff[idx];
        OUTS(0, idx, 0) = wANew[idx];
        OUTS(1, idx, 0) = wBNew[idx];
    }
}

void AndersonMixing(SCFT_INFO *scftInfo, int nRec, double *del, double *outs)
{
    gsl_matrix_view uGnu;
    gsl_vector_view vGnu, aGnu;
    gsl_permutation *p;

    double *up, *vp, *ap;
    int s;

    up = (double *)malloc(sizeof(double) * (nRec - 1) * (nRec - 1));
    vp = (double *)malloc(sizeof(double) * (nRec - 1));
    ap = (double *)malloc(sizeof(double) * (nRec - 1));

    // Calculate the U-matrix and the V-vector
    for (int idxHist1 = 1; idxHist1 < nRec; ++idxHist1)
    {
        // Skip "posUpdate"

        // if (idxHist1 == scftInfo->posUpdate)
        //     continue;
        // if (idxHist1 > scftInfo->posUpdate)
        //     idxHist1 -= 1;

        V(idxHist1) = 0.0;
        double v = 0.0;

#pragma omp parallel for reduction(+ \
                                   : v)
        for (int idx = 0; idx < scftInfo->size; ++idx)
        {
            // v += scftInfo->wts[idx] * (DEL(0, idx, scftInfo->posUpdate) - DEL(0, idx, idxHist1)) * DEL(0, idx, scftInfo->posUpdate);
            // v += scftInfo->wts[idx] * (DEL(1, idx, scftInfo->posUpdate) - DEL(1, idx, idxHist1)) * DEL(1, idx, scftInfo->posUpdate);
            v += scftInfo->wts[idx] * (DEL(0, idx, 0) - DEL(0, idx, idxHist1)) * DEL(0, idx, 0);
            v += scftInfo->wts[idx] * (DEL(1, idx, 0) - DEL(1, idx, idxHist1)) * DEL(1, idx, 0);
        }

        V(idxHist1) = v;

        for (int idxHist2 = idxHist1; idxHist2 < nRec; ++idxHist2)
        {
            // if (idxHist2 == scftInfo->posUpdate)
            //     continue;
            // if (idxHist2 > scftInfo->posUpdate)
            //     idxHist2 -= 1;

            U(idxHist1, idxHist2) = 0.0;
            double u = 0.0;

#pragma omp parallel for reduction(+ \
                                   : u)
            for (int idx = 0; idx < scftInfo->size; ++idx)
            {
                // u += scftInfo->wts[idx] * (DEL(0, idx, scftInfo->posUpdate) - DEL(0, idx, idxHist1)) * (DEL(0, idx, scftInfo->posUpdate) - DEL(0, idx, idxHist2));
                // u += scftInfo->wts[idx] * (DEL(1, idx, scftInfo->posUpdate) - DEL(1, idx, idxHist1)) * (DEL(1, idx, scftInfo->posUpdate) - DEL(1, idx, idxHist2));
                u += scftInfo->wts[idx] * (DEL(0, idx, 0) - DEL(0, idx, idxHist1)) * (DEL(0, idx, 0) - DEL(0, idx, idxHist2));
                u += scftInfo->wts[idx] * (DEL(1, idx, 0) - DEL(1, idx, idxHist1)) * (DEL(1, idx, 0) - DEL(1, idx, idxHist2));
            }

            U(idxHist1, idxHist2) = u;
            U(idxHist2, idxHist1) = u;
        }
    }

    // Calculate A - uses GNU LU decomposition for U A = V
    uGnu = gsl_matrix_view_array(up, nRec - 1, nRec - 1);
    vGnu = gsl_vector_view_array(vp, nRec - 1);
    aGnu = gsl_vector_view_array(ap, nRec - 1);

    p = gsl_permutation_alloc(nRec - 1);

    gsl_linalg_LU_decomp(&uGnu.matrix, p, &s);
    gsl_linalg_LU_solve(&uGnu.matrix, p, &vGnu.vector, &aGnu.vector);

    // Update fields
    for (int idx = 0; idx < scftInfo->size; ++idx)
    {
        // scftInfo->wA[idx] = OUTS(0, idx, scftInfo->posUpdate);
        // scftInfo->wB[idx] = OUTS(1, idx, scftInfo->posUpdate);
        scftInfo->wA[idx] = OUTS(0, idx, 0);
        scftInfo->wB[idx] = OUTS(1, idx, 0);

        // for (int idxHist = 0; idxHist < nRec; ++idxHist)
        for (int idxHist = 1; idxHist < nRec; ++idxHist)
        {
            // Skip "posUpdate"

            // if (idxHist == scftInfo->posUpdate)
            //     continue;
            // if (idxHist > scftInfo->posUpdate)
            //     idxHist -= 1;

            // scftInfo->wA[idx] += A(idxHist) * (OUTS(0, idx, idxHist) - OUTS(0, idx, scftInfo->posUpdate));
            // scftInfo->wB[idx] += A(idxHist) * (OUTS(1, idx, idxHist) - OUTS(1, idx, scftInfo->posUpdate));
            scftInfo->wA[idx] += A(idxHist) * (OUTS(0, idx, idxHist) - OUTS(0, idx, 0));
            scftInfo->wB[idx] += A(idxHist) * (OUTS(1, idx, idxHist) - OUTS(1, idx, 0));
        }
    }

    gsl_permutation_free(p);

    free(ap);
    free(vp);
    free(up);
}

// Simple mixing or Anderson mixing
void calcNewFlds(SCFT_INFO *scftInfo, int itr)
{
    double *wANew, *wBNew, *wADiff, *wBDiff;
    double errDiff = 0.0, errW = 0.0;
    // double err1 = 0.0;

    wANew = (double *)malloc(sizeof(double) * scftInfo->size);
    wBNew = (double *)malloc(sizeof(double) * scftInfo->size);
    wADiff = (double *)malloc(sizeof(double) * scftInfo->size);
    wBDiff = (double *)malloc(sizeof(double) * scftInfo->size);

#pragma omp parallel for
    for (int idx = 0; idx < scftInfo->size; ++idx)
    {
        wANew[idx] = scftInfo->hAB * scftInfo->phB[idx] + scftInfo->ksi[idx];
        wBNew[idx] = scftInfo->hAB * scftInfo->phA[idx] + scftInfo->ksi[idx];
    }
    if (scftInfo->isRemoveAvgFld == 1)
    {
        removeAvgFld(scftInfo, wANew, wBNew);
    }

#pragma omp parallel for reduction(+ \
                                   : errDiff, errW)
    for (int idx = 0; idx < scftInfo->size; ++idx)
    {
        wADiff[idx] = wANew[idx] - scftInfo->wA[idx];
        wBDiff[idx] = wBNew[idx] - scftInfo->wB[idx];

        errDiff += scftInfo->wts[idx] * (wADiff[idx] * wADiff[idx] + wBDiff[idx] * wBDiff[idx]);
        errW += scftInfo->wts[idx] * (scftInfo->wA[idx] * scftInfo->wA[idx] +
                                      scftInfo->wB[idx] * scftInfo->wB[idx]);
        // err1 += scftInfo->wts[idxLat] * (wADiff[idx] * wADiff[idx] + wBDiff[idx] * wBDiff[idx]);
    }

    scftInfo->err = sqrt(errDiff / errW);

    // err1 *= 2 * PI / scftInfo->nphi;
    // err1 /= 4 * PI;
    // err1 /= scftInfo->hAB;
    // scftInfo->err = err1;

    updateFldsHist(scftInfo, scftInfo->del, scftInfo->outs, wADiff, wBDiff, wANew, wBNew);

    if (scftInfo->err > scftInfo->errThreshold || itr <= scftInfo->stepThreshold)
    {
        // Simple mixing
#pragma omp parallel for
        for (int idx = 0; idx < scftInfo->size; ++idx)
        {
            scftInfo->wA[idx] += scftInfo->acc * wADiff[idx];
            scftInfo->wB[idx] += scftInfo->acc * wBDiff[idx];
        }
    }
    else
    {
        // Anderson mixing
        // printf("*** Anderson mixing ***\n");
        // FILE *fp = fopen("printout.txt", "a");
        // fprintf(fp, "*** Anderson mixing ***\n");
        // fclose(fp);

        int nRec = (itr - 1) < scftInfo->nHist ? (itr - 1) : scftInfo->nHist;
        AndersonMixing(scftInfo, nRec, scftInfo->del, scftInfo->outs);
    }

    if (scftInfo->isRemoveAvgFld == 1)
    {
        removeAvgFld(scftInfo, scftInfo->wA, scftInfo->wB);
    }

    free(wANew);
    free(wBNew);
    free(wADiff);
    free(wBDiff);
}

void calcJointDistr(SCFT_INFO *scftInfo)
{
    int offset = 0;
    double *jointCurr, *qFwdJoint, *qBwdJoint;

    for (int count = 0; count <= scftInfo->blockNum; ++count)
    {
        offset += count != scftInfo->blockNum ? scftInfo->blockInfo[count].offset
                                              : scftInfo->blockInfo[count - 1].ns * scftInfo->size;
        jointCurr = scftInfo->joint + count * scftInfo->size;
        qFwdJoint = scftInfo->qFwd + offset;
        qBwdJoint = scftInfo->qBwd + offset;

        for (int idx = 0; idx < scftInfo->size; ++idx)
        {
            jointCurr[idx] = qFwdJoint[idx] * qBwdJoint[idx] / scftInfo->Q;
        }
    }
}