#ifndef STRUCT_H
#define STRUCT_H

typedef enum BLOCK_TYPE {
    ATYPE = 1, BTYPE = 2
} BLOCK_TYPE;

typedef enum BRANCH_TYPE {
    LEFT = 1, RIGHT = 2
} BRANCH_TYPE;

typedef enum SOLVE_DIR {
    FWD = 1, BWD = -1
} SOLVE_DIR;

typedef struct BLOCK_INFO {
    int start;
    int end;
    BLOCK_TYPE type;
    int multi;
    int coef;
    BRANCH_TYPE branch;
    double len;
    int ns;
    double ds;
    int offset;
} BLOCK_INFO;

typedef struct SOLVE_SCHEME
{
    double *q;
    BLOCK_TYPE type;
    int initPow[40];
    int offset;
    int ns;
    double ds;
    SOLVE_DIR dir;
} SOLVE_SCHEME;

typedef struct SCFT_INFO
{
    // AB-type multiblock copolymer params
    double fA, fB;
    double hAB;

    int blockNum;
    BLOCK_INFO blockInfo[20];
    SOLVE_SCHEME solveScheme[20];

    // SCFT params and vars
    int lmax;                               // Maximum degree of spherical harmonics
    int mmax;                               // Maximum order of spherical harmonics
    int nlat;                               // Number of points in the latitude direction
    int nphi;                               // Number of points in the longitude direction
    int size;                               // Number of the grid points (size = nlat * nphi)
    int nlm;                                // Number of spherical harmonic components
    int totalNs;                            // Number of segments 
    int isGauss;                            // Gaussian grid (isGauss = 1) or regular grid (isGauss = 0)

    double r;                               // Radius of the sphere

    int itrMax;                             // Maximum iteration step
    int threadNum;                          // Number of threads used

    double Q;                               // Single-chain partition function
    double freeE, freeAB, freeW, freeS;     // Free energy and the related contributions
    double err, incompMax;                  // Error, incompressiblity coefficient
    int isRemoveAvgFld;                     // Flag for removing the average of the fields

    double acc;                             // Acceptance of the calculated new fields in simple mixing

    double errThreshold;                    // Threshold for the err to enter Anderson mixing
    double stepThreshold;                   // Threshold for the iteration step to enter Anderson mixing
    int nHist;                              // Number of "history" in Anderson mixing

    int phoutStep;
    int printoutStep;

    // AB diblock copolymer arrs
    
    double *phA, *phB;
    double *wA, *wB;
    double *ksi;

    double *qInit;
    double *qFwd, *qBwd;

    // SCFT arrs
    double *wts;

    double *wds2;
    double *lds;

    // Anderson mixing
    int posUpdate;
    double *del, *outs;

    double *joint;

} SCFT_INFO;

#endif // STRUCT_H