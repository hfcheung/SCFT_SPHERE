#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include <shtns.h>
#include <gsl/gsl_integration.h>

#include "struct.h"
#include "solver.h"

#define PI 3.1415926535897932385
#define LM_IDX(lmax, l, m) (((((unsigned short)(m)) * (2 * lmax + 2 - ((m) + 1))) >> 1) + (l))

// Intergral for Gaussian grid
// nphi: equally spaced nodes in longitude, spanning the range of angle
//       between 0 (included) and 2PI (excluded)
// nlat: gauss nodes in latitude, spanning angles from 0 to PI
//       renormalized to [-1, 1] Gauss-Legendre quadrature

// Integral for regular grid
// nphi: equally spaced nodes in longitude, spanning the range of angle
//       between 0 (included) and 2PI (excluded)
// nlat: equally spaced nodes in longitude, spanning the range of angle
//       between 0 (included) and PI (included)
//       renormalized to [-1, 1] Fejer quadrature

void randomInitFlds(SCFT_INFO *scftInfo)
{
#pragma omp parallel for
    for (int idx = 0; idx < scftInfo->size; ++idx)
    {
        scftInfo->wA[idx] = scftInfo->hAB * scftInfo->fB * (1 + 0.4 * (rand() / (double)RAND_MAX - 0.5));
        scftInfo->wB[idx] = scftInfo->hAB * scftInfo->fA * (1 + 0.4 * (rand() / (double)RAND_MAX - 0.5));
    }
}

void readInitFlds(SCFT_INFO *scftInfo, char *initFile, int flag)
{
    FILE *fp = fopen(initFile, "r");
    if (fp == NULL)
    {
        printf("Can't open the file \"%s\".\n", initFile);
        exit(EXIT_FAILURE);
    }
    if (flag == 0)
    {
        for (int idx = 0; idx < scftInfo->size; ++idx)
        {
            fscanf(fp, "%lf %lf %*[^\n]\n", &scftInfo->phA[idx], &scftInfo->phB[idx]);
            scftInfo->wA[idx] = scftInfo->hAB * scftInfo->phB[idx] + 0.4 * (rand() / (double)RAND_MAX - 0.5);
            scftInfo->wB[idx] = scftInfo->hAB * scftInfo->phA[idx] + 0.4 * (rand() / (double)RAND_MAX - 0.5);
        }
    }
    else
    {
        for (int idx = 0; idx < scftInfo->size; ++idx)
        {
            fscanf(fp, "%lf %lf %lf %lf\n", &scftInfo->phA[idx], &scftInfo->phB[idx],
                   &scftInfo->wA[idx], &scftInfo->wB[idx]);
        }
    }
    fclose(fp);
}

void calcWts(SCFT_INFO *scftInfo)
{
    if (scftInfo->isGauss == 1)
    {
        gsl_integration_fixed_workspace *wksp; // A workspace for computing integrals
        const gsl_integration_fixed_type *T = gsl_integration_fixed_legendre;
        wksp = gsl_integration_fixed_alloc(T, scftInfo->nlat, -1.0, 1.0, 0.0, 0.0);

        for (int idxPhi = 0; idxPhi < scftInfo->nphi; ++idxPhi)
        {
            for (int idxLat = 0; idxLat < scftInfo->nlat; ++idxLat)
            {
                int idx = idxPhi * scftInfo->nlat + idxLat;
                scftInfo->wts[idx] = wksp->weights[idxLat];
                scftInfo->wts[idx] *= 2 * PI / scftInfo->nphi;
            }
        }

        gsl_integration_fixed_free(wksp);
    }
    else
    {
        for (int idxPhi = 0; idxPhi < scftInfo->nphi; ++idxPhi)
        {
            for (int idxLat = 0; idxLat < scftInfo->nlat; ++idxLat)
            {
                int idx = idxPhi * scftInfo->nlat + idxLat;
                scftInfo->wts[idx] = idxLat == 0 ? 0.5 : idxLat == scftInfo->nlat - 1 ? 0.5
                                                                                      : 1.0;
                scftInfo->wts[idx] *= sin(idxLat * PI / (scftInfo->nlat - 1));
                scftInfo->wts[idx] *= PI / (scftInfo->nlat - 1) * 2 * PI / scftInfo->nphi;
            }
        }
    }

    // Or: shtns_gauss_wts(shtns, scftInfo->wts);
}

void calcKsi(SCFT_INFO *scftInfo)
{
#pragma omp parallel for
    for (int idx = 0; idx < scftInfo->size; ++idx)
    {
        scftInfo->ksi[idx] = 0.5 * (scftInfo->wA[idx] + scftInfo->wB[idx] - scftInfo->hAB);
    }
}

void calcWds(SCFT_INFO *scftInfo, double *w, double ds)
{
#pragma omp parallel for
    for (int idx = 0; idx < scftInfo->size; ++idx)
    {
        scftInfo->wds2[idx] = exp(-0.5 * w[idx] * ds);
    }
}

void calcLds(SCFT_INFO *scftInfo, double ds)
{
    // m dominant

    for (int idxL = 0; idxL <= scftInfo->lmax; ++idxL)
    {
        double t = -1.0 * ds * idxL * (idxL + 1) / (scftInfo->r * scftInfo->r);
        for (int idxM = 0; idxM <= idxL && idxM <= scftInfo->mmax; ++idxM)
        {
            scftInfo->lds[LM_IDX(scftInfo->lmax, idxL, idxM)] = exp(t);
        }
    }
}

void solveDiffEqn(SCFT_INFO *scftInfo, int solveCount)
{
    double *q, *w;
    int offset, ns;
    double ds;
    SOLVE_DIR dir;

    q = scftInfo->solveScheme[solveCount].q;
    offset = scftInfo->solveScheme[solveCount].offset;
    ns = scftInfo->solveScheme[solveCount].ns;
    ds = scftInfo->solveScheme[solveCount].ds;
    dir = scftInfo->solveScheme[solveCount].dir;

    q += offset;

    w = scftInfo->solveScheme[solveCount].type == ATYPE ? scftInfo->wA : scftInfo->wB;

    calcWds(scftInfo, w, ds);
    calcLds(scftInfo, ds);

    complex double *Slm; // Spherical harmonics coefficients
    double *in;          // Real space: theta, phi
    double *qCurr = (dir == FWD) ? q : q + ns * scftInfo->size;

    shtns_cfg shtns; // Handle to sht transform configuration
    // shtns_verbose(1);

    if (scftInfo->isGauss)
    {
        shtns = shtns_init(sht_gauss, scftInfo->lmax, scftInfo->mmax, 1, scftInfo->nlat, scftInfo->nphi);
    }
    else
    {
        shtns = shtns_init(sht_reg_dct, scftInfo->lmax, scftInfo->mmax, 1, scftInfo->nlat, scftInfo->nphi);
    }

    shtns_use_threads(scftInfo->threadNum);

    // Allocate spatial fields
    in = (double *)shtns_malloc(NSPAT_ALLOC(shtns) * sizeof(double));
    // Allocate SH representations
    Slm = (complex double *)shtns_malloc(scftInfo->nlm * sizeof(complex double));

    for (int idx = 0; idx < scftInfo->size; ++idx)
    {
        in[idx] = scftInfo->qInit[idx];
        qCurr[idx] = scftInfo->qInit[idx];
    }

    for (int idxS = 0; idxS < ns; ++idxS)
    {
#pragma omp parallel for
        for (int idx = 0; idx < scftInfo->size; ++idx)
        {
            in[idx] *= scftInfo->wds2[idx];
        }

        // SH analysis
        spat_to_SH(shtns, in, Slm); // Warning: this destroys the spatial data

        for (int idxL = 0; idxL <= scftInfo->lmax; ++idxL)
        {
            for (int idxM = 0; idxM <= idxL && idxM <= scftInfo->mmax; ++idxM)
            {
                Slm[LM_IDX(scftInfo->lmax, idxL, idxM)] *= scftInfo->lds[LM_IDX(scftInfo->lmax, idxL, idxM)];
            }
        }

        // SH synthesis
        SH_to_spat(shtns, Slm, in);

#pragma omp parallel for
        for (int idx = 0; idx < scftInfo->size; ++idx)
        {
            in[idx] *= scftInfo->wds2[idx];
        }

        qCurr = (dir == FWD) ? qCurr + scftInfo->size : qCurr - scftInfo->size;

        for (int idx = 0; idx < scftInfo->size; ++idx)
        {
            qCurr[idx] = in[idx];
        }
    }

    shtns_free(in);
    shtns_free(Slm);

    shtns_destroy(shtns);
}

void calcQ(SCFT_INFO *scftInfo)
{
    double sum = 0.0;

#pragma omp parallel for reduction(+ \
                                   : sum)
    for (int idx = 0; idx < scftInfo->size; ++idx)
    {
        sum += scftInfo->wts[idx] * scftInfo->qBwd[idx]; 
    }
    sum /= 4 * PI;
    scftInfo->Q = sum;
}

void calcPh(SCFT_INFO *scftInfo)
{
    for (int idx = 0; idx < scftInfo->size; ++idx)
    {
        scftInfo->phA[idx] = 0.0;
        scftInfo->phB[idx] = 0.0;
    }

    for (int count = 0; count < scftInfo->blockNum; ++count)
    {
        int offset = scftInfo->blockInfo[count].offset;
        int offsetS = offset / scftInfo->size;
        int ns = scftInfo->blockInfo[count].ns;
        double ds = scftInfo->blockInfo[count].ds;
        int coef = scftInfo->blockInfo[count].coef;
        BLOCK_TYPE type = scftInfo->blockInfo[count].type;

        for (int idxS = offsetS; idxS <= offsetS + ns; ++idxS)
        {
            double factor = idxS == offsetS ? 0.5
                                            : (idxS == offsetS + ns ? 0.5 : 1.0);
            factor /= scftInfo->Q;
            factor *= coef * ds;

            if (type == ATYPE)
            {
#pragma omp parallel for
                for (int idx = 0; idx < scftInfo->size; ++idx)
                {
                    scftInfo->phA[idx] += factor * scftInfo->qFwd[idxS * scftInfo->size + idx] *
                                          scftInfo->qBwd[idxS * scftInfo->size + idx];
                }
            }
            else
            {
#pragma omp parallel for
                for (int idx = 0; idx < scftInfo->size; ++idx)
                {
                    scftInfo->phB[idx] += factor * scftInfo->qFwd[idxS * scftInfo->size + idx] *
                                          scftInfo->qBwd[idxS * scftInfo->size + idx];
                }
            }
        }
    }
}
