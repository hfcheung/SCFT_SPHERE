#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "solver.h"
#include "struct.h"
#include "topo.h"

void checkLenSum(SCFT_INFO *scftInfo)
{
    double sA = 0.0, sB = 0.0;

    for (int count = 0; count < scftInfo->blockNum; ++count)
    {
        double s = scftInfo->blockInfo[count].len * scftInfo->blockInfo[count].coef;
        if (scftInfo->blockInfo[count].type == ATYPE)
        {
            sA += s;
        }
        else
        {
            sB += s;
        }
    }

    if (fabs(sA + sB - 1.0) > 1e-3)
    {
        printf("The total length of the chain is not 1.\n");
        exit(EXIT_FAILURE);
    }

    if (fabs(sA - scftInfo->fA) > 1e-3)
    {
        printf("The total length of the A blocks is not %f\n.", scftInfo->fA);
        exit(EXIT_FAILURE);
    }

    if (fabs(sB - scftInfo->fB) > 1e-3)
    {
        printf("The total length of the B blocks is not %f\n.", scftInfo->fB);
        exit(EXIT_FAILURE);
    }
}

void calcTotalNs(SCFT_INFO *scftInfo)
{
    int totalNs = 0;

    // Consider the case when f / ds is not an integer

    for (int count = 0; count < scftInfo->blockNum; ++count)
    {
        scftInfo->blockInfo[count].ns = (int)ceil((scftInfo->blockInfo[count].len / scftInfo->blockInfo[count].ds));
        scftInfo->blockInfo[count].ds = scftInfo->blockInfo[count].len / scftInfo->blockInfo[count].ns;
        scftInfo->blockInfo[count].offset = count == 0 ? 0
                                                       : scftInfo->blockInfo[count - 1].offset +
                                                             (scftInfo->blockInfo[count - 1].ns + 1) * scftInfo->size;
        totalNs += scftInfo->blockInfo[count].ns;
    }

    totalNs += scftInfo->blockNum;
    scftInfo->totalNs = totalNs;
}

void findSolveScheme(SCFT_INFO *scftInfo)
{
    // Find out the solving scheme

    int okFlagFwd[scftInfo->blockNum];
    int okFlagBwd[scftInfo->blockNum];
    int okCount = 0;

    for (int count = 0; count < scftInfo->blockNum; ++count)
    {
        okFlagFwd[count] = 0;
        okFlagBwd[count] = 0;
    }

    for (int count1 = 0; count1 < scftInfo->blockNum * 2; ++count1)
    {
        for (int count2 = 0; count2 <= scftInfo->blockNum * 2; ++count2)
        {
            scftInfo->solveScheme[count1].initPow[count2] = 0;
        }
    }

    int okFlag = 0, initFlag = 0;

    while (okCount < scftInfo->blockNum * 2)
    {
        for (int count = 0; count < scftInfo->blockNum; ++count)
        {
            int start = scftInfo->blockInfo[count].start;
            int end = scftInfo->blockInfo[count].end;

            if (okFlagFwd[count] != 1)
            {

                // Forward

                okFlag = 1;
                initFlag = 1;

                for (int itr = 0; itr < scftInfo->blockNum; ++itr)
                {
                    if (scftInfo->blockInfo[itr].end == start)
                    {
                        initFlag = 0;

                        if (okFlagFwd[itr] != 1)
                        {
                            okFlag = 0;
                            break;
                        }
                    }

                    if (scftInfo->blockInfo[itr].start == start && scftInfo->blockInfo[itr].branch == LEFT)
                    {
                        initFlag = 0;

                        if (okFlagBwd[itr] != 1)
                        {
                            okFlag = 0;
                            break;
                        }
                    }
                }

                // okFlag == 1 : can be solved

                if (okFlag == 1)
                {
                    scftInfo->solveScheme[okCount].q = scftInfo->qFwd;
                    scftInfo->solveScheme[okCount].type = scftInfo->blockInfo[count].type;
                    scftInfo->solveScheme[okCount].offset = scftInfo->blockInfo[count].offset;
                    scftInfo->solveScheme[okCount].ns = scftInfo->blockInfo[count].ns;
                    scftInfo->solveScheme[okCount].ds = scftInfo->blockInfo[count].ds;
                    scftInfo->solveScheme[okCount].dir = FWD;

                    if (initFlag == 1)
                    {
                        scftInfo->solveScheme[okCount].initPow[0] = -1;
                    }
                    else
                    {
                        for (int itr = 0; itr < scftInfo->blockNum; ++itr)
                        {
                            if (scftInfo->blockInfo[itr].end == start)
                            {
                                if (scftInfo->blockInfo[itr].branch == RIGHT)
                                {
                                    scftInfo->solveScheme[okCount].initPow[itr + 1] =
                                        scftInfo->blockInfo[itr].multi;
                                }
                                else
                                {
                                    scftInfo->solveScheme[okCount].initPow[itr + 1] = 1;
                                }
                            }

                            if (scftInfo->blockInfo[itr].start == start && scftInfo->blockInfo[itr].branch == LEFT)
                            {
                                if (itr == count) // Same block
                                {
                                    scftInfo->solveScheme[okCount].initPow[itr + scftInfo->blockNum + 1] =
                                        scftInfo->blockInfo[itr].multi - 1;
                                }
                                else
                                {
                                    scftInfo->solveScheme[okCount].initPow[itr + scftInfo->blockNum + 1] =
                                        scftInfo->blockInfo[itr].multi;
                                }
                            }
                        }
                    }

                    okFlagFwd[count] = 1;
                    okCount++;
                }
            }

            if (okFlagBwd[count] != 1)
            {

                // Backward

                okFlag = 1;
                initFlag = 1;

                for (int itr = 0; itr < scftInfo->blockNum; ++itr)
                {
                    if (scftInfo->blockInfo[itr].start == end)
                    {
                        initFlag = 0;

                        if (okFlagBwd[itr] != 1)
                        {
                            okFlag = 0;
                            break;
                        }
                    }

                    if (scftInfo->blockInfo[itr].end == end && scftInfo->blockInfo[itr].branch == RIGHT)
                    {
                        initFlag = 0;

                        if (okFlagFwd[itr] != 1)
                        {
                            okFlag = 0;
                            break;
                        }
                    }
                }

                // okFlag == 1 : can be solved

                if (okFlag == 1)
                {
                    scftInfo->solveScheme[okCount].q = scftInfo->qBwd;
                    scftInfo->solveScheme[okCount].type = scftInfo->blockInfo[count].type;
                    scftInfo->solveScheme[okCount].offset = scftInfo->blockInfo[count].offset;
                    scftInfo->solveScheme[okCount].ns = scftInfo->blockInfo[count].ns;
                    scftInfo->solveScheme[okCount].ds = scftInfo->blockInfo[count].ds;
                    scftInfo->solveScheme[okCount].dir = BWD;

                    if (initFlag == 1)
                    {
                        scftInfo->solveScheme[okCount].initPow[0] = -1;
                    }
                    else
                    {
                        for (int itr = 0; itr < scftInfo->blockNum; ++itr)
                        {
                            if (scftInfo->blockInfo[itr].start == end)
                            {
                                if (scftInfo->blockInfo[itr].branch == LEFT)
                                {
                                    scftInfo->solveScheme[okCount].initPow[itr + scftInfo->blockNum + 1] =
                                        scftInfo->blockInfo[itr].multi;
                                }
                                else
                                {
                                    scftInfo->solveScheme[okCount].initPow[itr + scftInfo->blockNum + 1] = 1;
                                }
                            }

                            if (scftInfo->blockInfo[itr].end == end && scftInfo->blockInfo[itr].branch == RIGHT)
                            {
                                if (itr == count) // Same block
                                {
                                    scftInfo->solveScheme[okCount].initPow[itr + 1] =
                                        scftInfo->blockInfo[itr].multi - 1;
                                }
                                else
                                {
                                    scftInfo->solveScheme[okCount].initPow[itr + 1] =
                                        scftInfo->blockInfo[itr].multi;
                                }
                            }
                        }
                    }

                    okFlagBwd[count] = 1;
                    okCount++;
                }
            }
        }
    }
}

void runSolveScheme(SCFT_INFO *scftInfo)
{
    for (int solveCount = 0; solveCount < 2 * scftInfo->blockNum; ++solveCount)
    {
        for (int idx = 0; idx < scftInfo->size; ++idx)
        {
            scftInfo->qInit[idx] = 1.0;
        }

        if (scftInfo->solveScheme[solveCount].initPow[0] != -1)
        {
            for (int count = 1; count <= scftInfo->blockNum * 2; ++count)
            {
                int power = scftInfo->solveScheme[solveCount].initPow[count];

                if (power == 0)
                {
                    continue;
                }

                double *q;

                // Consider AB_2 and A_2B

                q = count <= scftInfo->blockNum ? scftInfo->qFwd + scftInfo->blockInfo[count - 1].offset + scftInfo->blockInfo[count - 1].ns * scftInfo->size
                                                : scftInfo->qBwd + scftInfo->blockInfo[count - scftInfo->blockNum - 1].offset;

                if (power == 1)
#pragma omp parallel for
                    for (int idx = 0; idx < scftInfo->size; ++idx)
                    {
                        scftInfo->qInit[idx] *= q[idx];
                    }
                else
                {
#pragma omp parallel for
                    for (int idx = 0; idx < scftInfo->size; ++idx)
                    {
                        scftInfo->qInit[idx] *= pow(q[idx], power);
                    }
                }
            }
        }

        solveDiffEqn(scftInfo, solveCount);
    }
}