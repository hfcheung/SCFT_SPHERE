#ifndef SOLVER_H
#define SOLVER_H

#include "struct.h"

void randomInitFlds(SCFT_INFO *scftInfo);
void readInitFlds(SCFT_INFO *scftInfo, char *initFile, int flag);
void calcWts(SCFT_INFO* scftInfo);
void calcKsi(SCFT_INFO *scftInfo);
void calcWds(SCFT_INFO *scftInfo, double *w, double ds);
void calcLds(SCFT_INFO *scftInfo, double ds);
void solveDiffEqn(SCFT_INFO *scftInfo, int solveCount);
void calcQ(SCFT_INFO *scftInfo);
void calcPh(SCFT_INFO *scftInfo);

#endif // SOLVER_H