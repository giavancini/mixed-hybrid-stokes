#include <iostream>
#include <fstream>
#include <string>
#include <TPZMultiphysicsCompMesh.h>
#include <TPZLinearAnalysis.h>
#include <TPZSSpStructMatrix.h>
#include <pzstepsolver.h>
#include <pzfstrmatrix.h>
#include <pzlog.h>
#include <TPZSSpStructMatrix.h>

#include"ProblemData.h"
#include "SimpleExample.h"
#include "TPZMeshOperator.h"


int main(){
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG("Stokes.cfg");
#endif
    
    bool printdata = false;
    
    std::string filepath = "../DataInput/";
    std::string filenamejson =  "LidDrivenFlow.json";

    ProblemData simData;
    simData.ReadJson(filepath+filenamejson);

    TPZGeoMesh* gmesh = TPZMeshOperator::CreateGMesh(&simData);
    TPZCompMesh* cmesh_v = TPZMeshOperator::CreateCMeshV(&simData, gmesh);
    TPZCompMesh* cmesh_p = TPZMeshOperator::CreateCmeshP(&simData, gmesh);

    TPZMultiphysicsCompMesh* cmesh_m = TPZMeshOperator::CreateMultiPhysicsMesh(&simData, gmesh);

    if(printdata){
        simData.Print();
        TPZMeshOperator::PrintMesh(gmesh, cmesh_v, cmesh_p, cmesh_m);
    }

    TPZLinearAnalysis an(cmesh_m, true);
    TPZSSpStructMatrix<> strmat(cmesh_m);

    strmat.SetNumThreads(6);
    an.SetStructuralMatrix(strmat);

    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);

    an.Assemble();
    
    an.Solve();
    
    //vtk export
    TPZVec<std::string> scalarVars(1), vectorVars(1);
    scalarVars[0] = "Pressure";
    vectorVars[0] = "Velocity";
    
    an.DefineGraphMesh(simData.Dim(),scalarVars,vectorVars,"StokesSolution.vtk");
    constexpr int resolution{0};
    an.PostProcess(resolution);

    std::cout << "\n\nSimulation finished without errors :) \n\n";
	return 0;
}
