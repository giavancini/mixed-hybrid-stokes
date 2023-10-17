#include <iostream>
#include <fstream>
#include <string>

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "TPZMultiphysicsCompMesh.h"
#include "pzmatred.h"
#include "pzvec.h"
#include "TPZGenGrid2D.h"
#include "TPZGenGrid3D.h"
#include "TPZGeoMeshTools.h"
#include "TPZVTKGeoMesh.h"
#include "pzlog.h"
#include "TPZDarcy2Spaces.h"
#include "TPZNullMaterial.h"
#include "TPZNullMaterialCS.h"
#include "pzintel.h"
#include "pzcompel.h"
#include "TPZLinearAnalysis.h"
#include "TPZMatRedStructMatrix.h"
#include "pzfstrmatrix.h"

enum EMatid {EDomain,EBCbot,EBCright,EBCtop,EBCleft,EBCfront,EBCback};

TPZCompMesh* CreateH1CMesh(TPZGeoMesh* gmesh, const int pord);

TPZCompMesh* CreateGCMesh(TPZGeoMesh* gmesh);

TPZCompMesh* CreateMultiPhysicsMesh(TPZGeoMesh *gmesh, const TPZVec<TPZCompMesh*>& meshes);

int main()
{
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG("Stokes.cfg");
#endif

    TPZGeoMesh* gmesh = nullptr;

    const int dim = 3;

    const int division = 1;

    // lower corner
    const TPZManVector<REAL, 3> minX = {0., 0., 0.};

    // upper corner
    const TPZManVector<REAL, 3> maxX = {1., 1., 1.};

    //Creating geometric element
    if (dim == 2)
    {
        const TPZManVector<int, 2> ndiv = {division, division};
        TPZManVector<int,5> matids = {EDomain,EBCbot,EBCright,EBCtop,EBCleft};
        MMeshType meshType = MMeshType::EQuadrilateral;
        bool createBoundEls = true;
        gmesh = TPZGeoMeshTools::CreateGeoMeshOnGrid(dim, minX, maxX, matids, ndiv, meshType, createBoundEls);
        std::ofstream filegvtk("GMesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filegvtk, true);
    }
    else if (dim == 3)
    {
        const TPZManVector<int, 3> ndiv = {division, division, division};
        TPZManVector<int,7> matids = {EDomain,EBCbot,EBCright,EBCtop,EBCleft,EBCfront,EBCback};
        MMeshType meshType = MMeshType::EHexahedral;
        bool createBoundEls = true;
        gmesh = TPZGeoMeshTools::CreateGeoMeshOnGrid(dim, minX, maxX, matids, ndiv, meshType, createBoundEls);
        std::ofstream filegvtk("GMesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filegvtk, true);
    }

    //creating h1 mesh
    TPZCompMesh* cmesh_h1 = CreateH1CMesh(gmesh, 2);

    //creating constant mesh
    TPZCompMesh* cmesh_g = CreateGCMesh(gmesh);

    //creating multiphysics mesh
    TPZManVector<TPZCompMesh*,2> cmeshes = {cmesh_h1, cmesh_g};
    TPZCompMesh* cmesh_m = CreateMultiPhysicsMesh(gmesh, cmeshes);

    std::ofstream VTKCompMeshFile(cmesh_m->Name() + ".vtk");

    TPZVTKGeoMesh::PrintCMeshVTK(cmesh_m, VTKCompMeshFile);

    TPZLinearAnalysis analisys(cmesh_m, RenumType::ENone);

    TPZMatRedStructMatrix<> strmat(cmesh_m);
    strmat.Set

    analisys.SetStructuralMatrix(strmat);



    return 0;
}

TPZCompMesh* CreateH1CMesh(TPZGeoMesh* gmesh, const int pord)
{
    TPZCompMesh* cmesh = new TPZCompMesh(gmesh);
    cmesh->SetName("CMesh_H1");
    const int dim = gmesh->Dimension();
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(pord);
    cmesh->SetAllCreateFunctionsContinuous();
    
    // Domain elas mat
    TPZNullMaterial<>* mat = new TPZNullMaterial<>(EDomain);
    // mat->SetExactSol(elas->ExactSolution(), 2);
    // mat->SetForcingFunction(elas->ForceFunc(), 2);
    cmesh->InsertMaterialObject(mat);
    
    // BC null mat
    TPZFMatrix<STATE> val1(3,3,0.);
    TPZManVector<STATE> val2(3,0.);
    
    const int diri = 0, neu = 1;
    auto* BCCondBot = mat->CreateBC(mat, EBCbot, neu, val1, val2);
//    BCCond0->SetForcingFunctionBC(elas->ExactSolution(), 2);
    cmesh->InsertMaterialObject(BCCondBot);
    
    auto* BCCondRight = mat->CreateBC(mat, EBCright, neu, val1, val2);
//    BCCond0->SetForcingFunctionBC(elas->ExactSolution(), 2);
    cmesh->InsertMaterialObject(BCCondRight);

    auto* BCCondTop = mat->CreateBC(mat, EBCtop, neu, val1, val2);
//    BCCond0->SetForcingFunctionBC(elas->ExactSolution(), 2);
    cmesh->InsertMaterialObject(BCCondTop);

    auto* BCCondLeft = mat->CreateBC(mat, EBCleft, neu, val1, val2);
//    BCCond0->SetForcingFunctionBC(elas->ExactSolution(), 2);
    cmesh->InsertMaterialObject(BCCondLeft);

    if (dim == 3)
    {
        auto* BCCondFront = mat->CreateBC(mat, EBCfront, neu, val1, val2);
//    BCCond0->SetForcingFunctionBC(elas->ExactSolution(), 2);
        cmesh->InsertMaterialObject(BCCondFront);

        auto* BCCondBack = mat->CreateBC(mat, EBCback, neu, val1, val2);
//    BCCond0->SetForcingFunctionBC(elas->ExactSolution(), 2);
        cmesh->InsertMaterialObject(BCCondBack);
    }
    
    // Constructs mesh
    cmesh->AutoBuild();
    
    return cmesh;
}

TPZCompMesh* CreateGCMesh(TPZGeoMesh* gmesh)
{
    TPZCompMesh* cmesh = new TPZCompMesh(gmesh);
    cmesh->SetName("CMesh_G");
    const int dim = gmesh->Dimension();
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(0);
    cmesh->SetAllCreateFunctionsDiscontinuous();
    
    // Domain elas mat
    TPZNullMaterial<>* mat = new TPZNullMaterial<>(EDomain);
    cmesh->InsertMaterialObject(mat);
    
    // Constructs mesh
    cmesh->AutoBuild();

    for (TPZCompEl* el : cmesh->ElementVec())
    {
        TPZCompElDisc* compeldisc = dynamic_cast<TPZCompElDisc*>(el);
        if (!compeldisc) DebugStop();
        compeldisc->SetConnectIndex(0,0);
    }

    cmesh->ComputeNodElCon();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
}

TPZCompMesh* CreateMultiPhysicsMesh(TPZGeoMesh *gmesh, const TPZVec<TPZCompMesh*>& meshes)
{
    TPZMultiphysicsCompMesh *cmesh_m = new TPZMultiphysicsCompMesh(gmesh);

    cmesh_m->SetName("CMesh_M");
    const int dim = gmesh->Dimension();
    cmesh_m->SetDefaultOrder(meshes[0]->GetDefaultOrder());
    cmesh_m->SetAllCreateFunctionsMultiphysicElem();

    // Creating Materials

    // 1. For domain
    TPZDarcy2Spaces* mat = new TPZDarcy2Spaces(EDomain, dim);
    cmesh_m->InsertMaterialObject(mat);

    // 2. Boundary Conditions
    TPZFMatrix<STATE> val1(3, 3, 0.);
    TPZManVector<STATE> val2(3, 0.);

    const int diri = 0, neu = 1;
    auto* BCCondBot = mat->CreateBC(mat, EBCbot, neu, val1, val2);
//    BCCond0->SetForcingFunctionBC(elas->ExactSolution(), 2);
    cmesh_m->InsertMaterialObject(BCCondBot);
    
    auto* BCCondRight = mat->CreateBC(mat, EBCright, neu, val1, val2);
//    BCCond0->SetForcingFunctionBC(elas->ExactSolution(), 2);
    cmesh_m->InsertMaterialObject(BCCondRight);

    auto* BCCondTop = mat->CreateBC(mat, EBCtop, neu, val1, val2);
//    BCCond0->SetForcingFunctionBC(elas->ExactSolution(), 2);
    cmesh_m->InsertMaterialObject(BCCondTop);

    auto* BCCondLeft = mat->CreateBC(mat, EBCleft, neu, val1, val2);
//    BCCond0->SetForcingFunctionBC(elas->ExactSolution(), 2);
    cmesh_m->InsertMaterialObject(BCCondLeft);

    if (dim == 3)
    {
        auto* BCCondFront = mat->CreateBC(mat, EBCfront, neu, val1, val2);
//    BCCond0->SetForcingFunctionBC(elas->ExactSolution(), 2);
        cmesh_m->InsertMaterialObject(BCCondFront);

        auto* BCCondBack = mat->CreateBC(mat, EBCback, neu, val1, val2);
//    BCCond0->SetForcingFunctionBC(elas->ExactSolution(), 2);
        cmesh_m->InsertMaterialObject(BCCondBack);
    }

    TPZManVector<int, 2> active_approx_spaces(2, 1);

    cmesh_m->BuildMultiphysicsSpace(active_approx_spaces, meshes);
    cmesh_m->AdjustBoundaryElements();
    cmesh_m->CleanUpUnconnectedNodes();
    cmesh_m->LoadReferences();

    return cmesh_m;
}
