#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include <omp.h>

#include "struct.h"
#include "solver.h"
#include "calc.h"
#include "topo.h"

#define PI 3.1415926535897932385

int main(int argc, char *argv[])
{
    clock_t start, end;
    double tRun = 0.0;

    time_t ts;
    int iseed = time(&ts);
    srand(iseed);

    SCFT_INFO scftInfo;

    FILE *fp;

    char type[20], branch[20];
    char gridType[20];
    char initMode[20];
    char initFile[20];
    char isRemoveAvgFldStr[20];

    char ATypeStr[2] = "A";
    char BTypeStr[2] = "B";
    char leftTypeStr[2] = "L";
    char rightTypeStr[2] = "R";
    char yesStr[2] = "Y";
    char noStr[2] = "N";

    char randomMode[20] = "random";
    char phiMode[20] = "phi";
    char omegaMode[20] = "omega";

    fp = fopen("params", "r");

    fscanf(fp, "%*[^\n]\n");
    fscanf(fp, "fA %lf\n", &scftInfo.fA);
    fscanf(fp, "hAB %lf\n", &scftInfo.hAB);
    fscanf(fp, "r %lf\n", &scftInfo.r);

    fscanf(fp, "%*[^\n]\n");
    fscanf(fp, "blockNum %d\n", &scftInfo.blockNum);
    fscanf(fp, "%*[^\n]\n");
    fscanf(fp, "%*[^\n]\n");
    for (int count = 0; count < scftInfo.blockNum; ++count)
    {
        fscanf(fp, "%d %d %s %d %d %s %lf %lf\n", &scftInfo.blockInfo[count].start, &scftInfo.blockInfo[count].end,
               &type, &scftInfo.blockInfo[count].multi, &scftInfo.blockInfo[count].coef,
               &branch, &scftInfo.blockInfo[count].len, &scftInfo.blockInfo[count].ds);

        if (scftInfo.blockInfo[count].start > scftInfo.blockNum ||
            scftInfo.blockInfo[count].end > scftInfo.blockNum ||
            scftInfo.blockInfo[count].start < 0 ||
            scftInfo.blockInfo[count].end < 0 ||
            scftInfo.blockInfo[count].start >= scftInfo.blockInfo[count].end)
        {
            printf("Invalid \"start\" or \"end\".\n");
            exit(EXIT_FAILURE);
        }

        if (!(scftInfo.blockInfo[count].multi > 0 && scftInfo.blockInfo[count].multi < 10))
        {
            printf("Invalid blockMulti.\n");
            exit(EXIT_FAILURE);
        }

        if (!(scftInfo.blockInfo[count].coef > 0 && scftInfo.blockInfo[count].coef < 200))
        {
            printf("Invalid blockMulti.\n");
            exit(EXIT_FAILURE);
        }

        if (scftInfo.blockInfo[count].len < 0 || scftInfo.blockInfo[count].len > 1)
        {
            printf("Invalid blockLen.\n");
            exit(EXIT_FAILURE);
        }

        if (strcmp(type, ATypeStr) == 0)
        {
            scftInfo.blockInfo[count].type = ATYPE;
        }
        else if (strcmp(type, BTypeStr) == 0)
        {
            scftInfo.blockInfo[count].type = BTYPE;
        }
        else
        {
            printf("Invalid blockType.\n");
            exit(EXIT_FAILURE);
        }

        if (strcmp(branch, leftTypeStr) == 0)
        {

            scftInfo.blockInfo[count].branch = LEFT;
        }
        else if (strcmp(branch, rightTypeStr) == 0)
        {
            scftInfo.blockInfo[count].branch = RIGHT;
        }
        else
        {
            printf("Invalid branchType.\n");
            exit(EXIT_FAILURE);
        }
    }

    fscanf(fp, "%*[^\n]\n");
    fscanf(fp, "lmax %d\n", &scftInfo.lmax);
    fscanf(fp, "mmax %d\n", &scftInfo.mmax);
    fscanf(fp, "nlat %d\n", &scftInfo.nlat);
    fscanf(fp, "nphi %d\n", &scftInfo.nphi);
    fscanf(fp, "gridType %s\n", &gridType);

    fscanf(fp, "%*[^\n]\n");
    fscanf(fp, "initMode %s\n", &initMode);
    fscanf(fp, "initFile %s\n", &initFile);

    fscanf(fp, "%*[^\n]\n");
    fscanf(fp, "itrMax %d\n", &scftInfo.itrMax);
    fscanf(fp, "threadNum %d\n", &scftInfo.threadNum);
    fscanf(fp, "removeAvgFld %s\n", &isRemoveAvgFldStr);

    fscanf(fp, "%*[^\n]\n");
    fscanf(fp, "%*[^\n]\n");
    fscanf(fp, "acc %lf\n", &scftInfo.acc);

    fscanf(fp, "%*[^\n]\n");
    fscanf(fp, "%*[^\n]\n");
    fscanf(fp, "errThreshold %lf\n", &scftInfo.errThreshold);
    fscanf(fp, "stepThreshold %lf\n", &scftInfo.stepThreshold);
    fscanf(fp, "nHist %d\n", &scftInfo.nHist);

    fscanf(fp, "%*[^\n]\n");
    fscanf(fp, "phoutStep %d\n", &scftInfo.phoutStep);
    fscanf(fp, "printoutStep %d\n", &scftInfo.printoutStep);

    fclose(fp);

    omp_set_num_threads(scftInfo.threadNum);

    scftInfo.fB = 1.0 - scftInfo.fA;
    scftInfo.size = scftInfo.nlat * scftInfo.nphi;

    scftInfo.nlm = (scftInfo.lmax + 1) * (scftInfo.mmax + 1);
    scftInfo.nlm -= (scftInfo.mmax * (scftInfo.mmax + 1)) / 2;

    checkLenSum(&scftInfo);
    calcTotalNs(&scftInfo);

    scftInfo.phA = (double *)malloc(sizeof(double) * scftInfo.size);
    scftInfo.phB = (double *)malloc(sizeof(double) * scftInfo.size);
    scftInfo.wA = (double *)malloc(sizeof(double) * scftInfo.size);
    scftInfo.wB = (double *)malloc(sizeof(double) * scftInfo.size);
    scftInfo.ksi = (double *)malloc(sizeof(double) * scftInfo.size);

    scftInfo.qInit = (double *)malloc(sizeof(double) * scftInfo.size);
    scftInfo.qFwd = (double *)malloc(sizeof(double) * scftInfo.size * scftInfo.totalNs);
    scftInfo.qBwd = (double *)malloc(sizeof(double) * scftInfo.size * scftInfo.totalNs);

    scftInfo.wts = (double *)malloc(sizeof(double) * scftInfo.size);

    scftInfo.wds2 = (double *)malloc(sizeof(double) * scftInfo.size);
    scftInfo.lds = (double *)malloc(sizeof(double) * scftInfo.nlm);

    scftInfo.del = (double *)malloc(sizeof(double) * scftInfo.size * scftInfo.nHist * 2);
    scftInfo.outs = (double *)malloc(sizeof(double) * scftInfo.size * scftInfo.nHist * 2);

    scftInfo.joint = (double *)malloc(sizeof(double) * scftInfo.size * (scftInfo.blockNum + 1));

    char gaussType[20] = "gauss";
    char regularType[20] = "regular";

    if (strcmp(gridType, gaussType) == 0)
    {
        scftInfo.isGauss = 1;
    }
    else if (strcmp(gridType, regularType) == 0)
    {
        scftInfo.isGauss = 0;
    }
    else
    {
        printf("Invalid gridType.\n");
        exit(EXIT_FAILURE);
    }

    if (strcmp(initMode, randomMode) == 0)
    {
        randomInitFlds(&scftInfo);
    }
    else if (strcmp(initMode, phiMode) == 0)
    {
        readInitFlds(&scftInfo, initFile, 0);
    }
    else if (strcmp(initMode, omegaMode) == 0)
    {
        readInitFlds(&scftInfo, initFile, 1);
    }
    else
    {
        printf("Invalid initMode.\n");
        exit(EXIT_FAILURE);
    }

    if (strcmp(isRemoveAvgFldStr, yesStr) == 0)
    {
        scftInfo.isRemoveAvgFld = 1;
    }
    else if (strcmp(isRemoveAvgFldStr, noStr) == 0)
    {
        scftInfo.isRemoveAvgFld = 0;
    }
    else
    {
        printf("Invalid isRemoveAvgFld answer.\n");
        exit(EXIT_FAILURE);
    }

    if (scftInfo.isRemoveAvgFld == 1)
    {
        removeAvgFld(&scftInfo, scftInfo.wA, scftInfo.wB);
    }

    calcWts(&scftInfo);
    findSolveScheme(&scftInfo);

    fp = fopen("printout.txt", "w");
    fprintf(fp, "fA = %.3f, fB = %.3f, hAB = %.3f, r = %.3f\n", scftInfo.fA,
            scftInfo.fB, scftInfo.hAB, scftInfo.r);
    fprintf(fp, "lmax = %d, mmax = %d, nlat = %d, nphi = %d, size = %d, gridType = %s\n", scftInfo.lmax, scftInfo.mmax,
            scftInfo.nlat, scftInfo.nphi, scftInfo.size, gridType);
    fprintf(fp, "initMode = %s, initFile = %s\n", initMode, initFile);
    fprintf(fp, "\n");

    fprintf(fp, "blockSeq\n");
    fprintf(fp, "start end type multi coef dir len ds\n");

    for (int count = 0; count < scftInfo.blockNum; ++count)
    {
        printf("%d %d %d %d %d %d %.6f %.6f\n", scftInfo.blockInfo[count].start, scftInfo.blockInfo[count].end,
               scftInfo.blockInfo[count].type, scftInfo.blockInfo[count].multi, scftInfo.blockInfo[count].coef,
               scftInfo.blockInfo[count].branch, scftInfo.blockInfo[count].len, scftInfo.blockInfo[count].ds);
        fprintf(fp, "%d %d %d %d %d %d %.6f %.6f\n", scftInfo.blockInfo[count].start, scftInfo.blockInfo[count].end,
                scftInfo.blockInfo[count].type, scftInfo.blockInfo[count].multi, scftInfo.blockInfo[count].coef,
                scftInfo.blockInfo[count].branch, scftInfo.blockInfo[count].len, scftInfo.blockInfo[count].ds);
    }
    fprintf(fp, "\n");

    fprintf(fp, "solveScheme\n");
    fprintf(fp, "type offsetS ns dir\n");

    for (int solveCount = 0; solveCount < scftInfo.blockNum * 2; ++solveCount)
    {
        fprintf(fp, "step %d\n", solveCount);
        fprintf(fp, "%d %d %d %d\n", scftInfo.solveScheme[solveCount].type,
                scftInfo.solveScheme[solveCount].offset / scftInfo.size, scftInfo.solveScheme[solveCount].ns,
                scftInfo.solveScheme[solveCount].dir);
        for (int j = 0; j <= 2 * scftInfo.blockNum; ++j)
        {
            fprintf(fp, "%d ", scftInfo.solveScheme[solveCount].initPow[j]);
        }
        fprintf(fp, "\n");
    }

    fprintf(fp, "\n\n");

    fprintf(fp, "               freeE           freeAB            freeW           freeS             err            incompMax\n");
    fclose(fp);

    int itr = 0;
    do
    {
        start = clock();

        runSolveScheme(&scftInfo);

        calcQ(&scftInfo);

        calcPh(&scftInfo);

        calcKsi(&scftInfo);

        calcFreeE(&scftInfo);

        calcNewFlds(&scftInfo, itr);

        printf("itr %-4d %16.9e %16.9e %16.9e\n", itr, scftInfo.freeE, scftInfo.err, scftInfo.incompMax);

        if (itr % scftInfo.phoutStep == 0)
        {
            fp = fopen("phout.txt", "w");
            for (int idx = 0; idx < scftInfo.size; ++idx)
            {
                fprintf(fp, "%.6f %.6f %.6f %.6f\n", scftInfo.phA[idx], scftInfo.phB[idx],
                        scftInfo.wA[idx], scftInfo.wB[idx]);
            }
            fclose(fp);
        }

        if (itr % scftInfo.printoutStep == 0)
        {
            fp = fopen("printout.txt", "a");
            fprintf(fp, "itr %-4d %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n", itr, scftInfo.freeE, scftInfo.freeAB,
                    scftInfo.freeW, scftInfo.freeS, scftInfo.err, scftInfo.incompMax);
            fclose(fp);
        }
        itr++;

        end = clock();
        tRun += (double)(end - start) / CLOCKS_PER_SEC;
        printf("%f s passed.\n", tRun);

        if (itr >= scftInfo.itrMax || (scftInfo.err < 1e-6 && scftInfo.incompMax < 1e-6))
        {
            fp = fopen("phout.txt", "w");
            for (int idx = 0; idx < scftInfo.size; ++idx)
            {
                fprintf(fp, "%.6f %.6f %.6f %.6f\n", scftInfo.phA[idx], scftInfo.phB[idx],
                        scftInfo.wA[idx], scftInfo.wB[idx]);
            }
            fclose(fp);

            fp = fopen("printout.txt", "a");
            fprintf(fp, "\n");
            fprintf(fp, "itr %-4d %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n", itr, scftInfo.freeE, scftInfo.freeAB,
                    scftInfo.freeW, scftInfo.freeS, scftInfo.err, scftInfo.incompMax);
            fclose(fp);

            fp = fopen("joint.txt", "w");
            calcJointDistr(&scftInfo);
            for (int idx = 0; idx < scftInfo.size; ++idx)
            {
                for (int count = 0; count <= scftInfo.blockNum; ++count)
                {
                    int offset = count * scftInfo.size;
                    fprintf(fp, "%.6f ", scftInfo.joint[idx + offset]);
                }
                fprintf(fp, "\n");
            }
            fclose(fp);

            break;
        }

    } while (1);

    free(scftInfo.phA);
    free(scftInfo.phB);
    free(scftInfo.wA);
    free(scftInfo.wB);
    free(scftInfo.ksi);
    free(scftInfo.qInit);
    free(scftInfo.qFwd);
    free(scftInfo.qBwd);
    free(scftInfo.wts);
    free(scftInfo.wds2);
    free(scftInfo.lds);
    free(scftInfo.del);
    free(scftInfo.outs);
    free(scftInfo.joint);

    return 0;
}