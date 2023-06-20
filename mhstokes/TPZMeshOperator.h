#ifndef STOKESEXAMPLE_H
#define STOKESEXAMPLE_H

#include "ProblemData.h"
#include <TPZAnalyticSolution.h>

class TPZMeshOperator{
public:
    enum SimSpace {EStandard, EConstant};
    
    static TPZGeoMesh* CreateGMesh(ProblemData* simData);
    static void InsertLambdaGEl(ProblemData* simData, TPZGeoMesh* gmesh);
    
    static TPZCompMesh* CreateCMeshV(ProblemData* simData, TPZGeoMesh* gmesh);
    
    static TPZCompMesh* CreateCmeshP(ProblemData* simData, TPZGeoMesh* gmesh);
    
    static TPZCompMesh* CreateCmeshMv(ProblemData* simData, TPZGeoMesh* gmesh);
    
    static TPZCompMesh* CreateCmeshMp(ProblemData* simData, TPZGeoMesh* gmesh);
    
    static TPZMultiphysicsCompMesh* CreateMultiPhysicsMesh(ProblemData* simData, TPZGeoMesh* gmesh);
    
    static void CondenseElements(TPZMultiphysicsCompMesh* cmesh_m);
    
    static void InsertBCInterfaces(TPZMultiphysicsCompMesh* cmesh_m, ProblemData* simData, TPZGeoMesh* gmesh);
    
    static void InsertInterfaces(TPZMultiphysicsCompMesh* cmesh_m, ProblemData* simData, TPZGeoMesh* gmesh);
    
    static void PrintGeoMesh(TPZGeoMesh* gmesh);

    static void PrintCompMesh(TPZVec<TPZCompMesh*> CMeshVec);
<<<<<<< HEAD
    
    static void PrintCompMesh(TPZCompMesh* cmesh);
    
    static void GenerateMshFile(ProblemData* geoFile);
||||||| parent of 9e51454 (Axisymmetric working)
=======

    static void PrintCompMesh(TPZCompMesh* cmesh);
>>>>>>> 9e51454 (Axisymmetric working)
};
#endif
