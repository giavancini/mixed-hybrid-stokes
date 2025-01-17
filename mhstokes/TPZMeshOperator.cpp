#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

#include <pzgmesh.h>
#include <TPZGenGrid2D.h>
#include <TPZVTKGeoMesh.h>
#include <TPZGmshReader.h>
#include <TPZLinearAnalysis.h>
#include <pzfmatrix.h>
#include <TPZGeoElement.h>
#include <TPZInterfaceEl.h>
#include <TPZMultiphysicsInterfaceEl.h>
#include <TPZMultiphysicsCompMesh.h>
#include <TPZGeoLinear.h>
#include <TPZNullMaterial.h>
#include <TPZNullMaterialCS.h>
#include <pzintel.h>
#include <pzelementgroup.h>
#include <pzcondensedcompel.h>
#include <filesystem>

#include "TPZInterfaceAxisymStokesMaterial.h"
#include "TPZInterface1dStokesMaterial.h"
#include "TPZInterface1dFlux.h"
#include "TPZStokesMaterial.h"
#include "TPZAxisymStokesMaterial.h"
#include "TPZ1dStokesMaterial.h"
#include "TPZMixedLinearElasticMaterial.h"
#include "TPZMeshOperator.h"
#include "ProblemData.h"

void TPZMeshOperator::GenerateMshFile(ProblemData *simData)
{
    std::string command;

    //    if(std::filesystem::exists(file)){
    //        std::cout << "\n ERROR: THERE IS NO SUCH FILE: " << simData->MeshName() << std::endl;
    //        DebugStop();
    //    }

    command = "gmsh " + simData->MeshName() + ".geo -o " + simData->MeshName() + ".msh -algo del2d -2 -order 1";
    system(command.c_str());
}

TPZGeoMesh *TPZMeshOperator::CreateGMesh(ProblemData *simData)
{

    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->SetName("GeoMesh");

    if (simData->CreateMshFile())
    {
        GenerateMshFile(simData);
    }

    TPZGmshReader reader;

    reader.GeometricGmshMesh(simData->MeshName() + ".msh", gmesh);

    TPZMeshOperator::InsertLambdaGEl(simData, gmesh);

    gmesh->BuildConnectivity();

    return gmesh;
}

void TPZMeshOperator::InsertLambdaGEl(ProblemData *simData, TPZGeoMesh *gmesh)
{
    int64_t nEl = gmesh->NElements();

    //We insert a lambda element either between two neightbour elements or at r=0 for axisymmetric problems

    for (int64_t el = 0; el < nEl; el++)
    {
        TPZGeoEl *geoEl = gmesh->Element(el);

        if (!geoEl)
            continue;
        if (geoEl->HasSubElement())
            continue;
        if (geoEl->Dimension() != gmesh->Dimension())
            continue;

        int nside = geoEl->NSides();

        if (simData->DomainVec().size() != 0)
        {
            for (int side = 0; side < nside; side++)
            {
                if (geoEl->SideDimension(side) != gmesh->Dimension() - 1)
                    continue;

                TPZGeoElSide geoElSide(geoEl, side);
                TPZGeoElSide neighbour = geoElSide.Neighbour();

                if (neighbour == geoElSide)
                    continue;
                if (neighbour.Element()->HasSubElement())
                    continue;

                while (neighbour != geoElSide)
                {

                    if (neighbour.Element()->Dimension() == gmesh->Dimension() - 1)
                    {
                        int neighbourMatId = neighbour.Element()->MaterialId();

                        break;
                    }

                    neighbour = neighbour.Neighbour();
                }

                if (neighbour == geoElSide)
                    TPZGeoElBC(geoElSide, simData->LambdaID());
            }
        }
        
        
    }

    for (int64_t el = 0; el < nEl; el++)
    {
        TPZGeoEl *geoEl = gmesh->Element(el);

        if (!geoEl)
            continue;
        if (geoEl->HasSubElement())
            continue;

        int matID = geoEl->MaterialId();

        if (simData->Axisymmetric() && simData->AxisymmetryDomainVec().size() != 0) //in this case, r=0 belongs to the domain
            for (auto axisymmetryDomain : simData->AxisymmetryDomainVec())
            {
                if (matID == axisymmetryDomain.matID) //Axisymmetric tube so we create a lambda material
                    TPZGeoElBC gbc(geoEl, 2, simData->AxiLambdaID());
            }
    }

}

void TPZMeshOperator::InsertBCInterfaces(TPZMultiphysicsCompMesh *cmesh_m, ProblemData *simData, TPZGeoMesh *gmesh)
{
    TPZManVector<int, 2> Interfaces(2, 0);
    Interfaces[0] = simData->InterfaceID();
    Interfaces[1] = -simData->InterfaceID();

    if (!gmesh)
        DebugStop();

    int64_t nel = gmesh->NElements();

    TPZVec<int> IDVec(simData->TangentialBCs().size(), 0);
    for (int i = 0; i < simData->TangentialBCs().size(); i++)
    {
        IDVec[i] = simData->TangentialBCs()[i].matID;
    }

    //For tangential boundary conditions
    for (auto const &BcMatID : IDVec)
    {
        for (int64_t el = 0; el < nel; el++)
        {
            TPZGeoEl *gel = gmesh->Element(el);
            int meshDim = gmesh->Dimension();
            int matID = gel->MaterialId();

            if (matID != BcMatID)
                continue;

            int nsides = gel->NSides();
            TPZGeoElSide gelSide(gel, nsides - 1);
            TPZCompElSide celSide = gelSide.Reference();

            TPZStack<TPZGeoElSide> neighbourSet;
            gelSide.AllNeighbours(neighbourSet);

            int64_t nneighs = neighbourSet.size();

            TPZManVector<int64_t, 1> LeftElIndex(1, 0), RightElIndex(1, 0);
            LeftElIndex[0] = 0;
            RightElIndex[0] = 1;

            for (int stack_i = 0; stack_i < nneighs; stack_i++)
            {
                TPZGeoElSide neigh = neighbourSet[stack_i];
                int neighMatID = neigh.Element()->MaterialId();
                TPZCompElSide celNeigh = neigh.Reference();

                int64_t neighIndex = neigh.Element()->Index();

                if (neigh.Element()->Dimension() != meshDim)
                    continue;

                if (neigh.Element()->HasSubElement())
                {
                    // Check if it is working in the case with refined meshes
                    DebugStop();
                }
                else
                {

                    TPZGeoElBC gbc(gelSide, Interfaces[0]);

                    TPZMultiphysicsInterfaceElement *interElem = new TPZMultiphysicsInterfaceElement(*cmesh_m, gbc.CreatedElement(), celNeigh, celSide);
                    interElem->SetLeftRightElementIndices(LeftElIndex, RightElIndex);
                }
            }
        }
    }

    //For the axisymmetric tube
    for (int64_t el = 0; el < nel; el++)
    {
        TPZGeoEl *gel = gmesh->Element(el);
        int meshDim = gmesh->Dimension();
        int matID = gel->MaterialId();

        if (matID == simData->AxiLambdaID()) 
        {
            int nsides = gel->NSides();
            TPZGeoElSide gelSide(gel, nsides - 1);
            TPZCompElSide celSide = gelSide.Reference();

            TPZStack<TPZGeoElSide> neighbourSet;
            gelSide.AllNeighbours(neighbourSet);

            int64_t nneighs = neighbourSet.size();

            TPZManVector<int64_t, 1> LeftElIndex(1, 0), RightElIndex(1, 0);
            LeftElIndex[0] = 0;
            RightElIndex[0] = 1;

            for (int stack_i = 0; stack_i < nneighs; stack_i++)
            {
                TPZGeoElSide neigh = neighbourSet[stack_i];
                int neighMatID = neigh.Element()->MaterialId();
                TPZCompElSide celNeigh = neigh.Reference();

                int64_t neighIndex = neigh.Element()->Index();

                if (neighMatID == simData->AxisymmetryDomainVec()[0].matID) //Creating an interface element between the 1d lambda and 1d stokes element
                {
                    if (neigh.Element()->HasSubElement())
                    {
                        // Check if it is working in the case with refined meshes
                        DebugStop();
                    }
                    else
                    {
                        TPZGeoElBC gbc(gelSide, simData->AxiInterfaceID());

                        TPZMultiphysicsInterfaceElement *interElem = new TPZMultiphysicsInterfaceElement(*cmesh_m, gbc.CreatedElement(), celNeigh, celSide);
                        interElem->SetLeftRightElementIndices(LeftElIndex, RightElIndex);
                    }
                }
                else if (neighMatID == simData->DomainVec()[0].matID) //Creating an interface element between the 1d lambda and 2d stokes element
                {
                    if (neigh.Element()->HasSubElement())
                    {
                        // Check if it is working in the case with refined meshes
                        DebugStop();
                    }
                    else
                    {
                        TPZGeoElBC gbc(gelSide, simData->InterfaceID());

                        TPZMultiphysicsInterfaceElement *interElem = new TPZMultiphysicsInterfaceElement(*cmesh_m, gbc.CreatedElement(), celNeigh, celSide);
                        interElem->SetLeftRightElementIndices(LeftElIndex, RightElIndex);
                    }
                }
            }
        }

        else if (matID == simData->AxisymmetryDomainVec()[0].matID) //Creating an interface element between 1d stokes element and 2d stokes element for normal flux computation
        {
            int nsides = gel->NSides();
            TPZGeoElSide gelSide(gel, nsides - 1);
            TPZCompElSide celSide = gelSide.Reference();

            TPZStack<TPZGeoElSide> neighbourSet;
            gelSide.AllNeighbours(neighbourSet);

            int64_t nneighs = neighbourSet.size();

            TPZManVector<int64_t, 1> LeftElIndex(1, 0), RightElIndex(1, 0);
            LeftElIndex[0] = 0;
            RightElIndex[0] = 1;

            for (int stack_i = 0; stack_i < nneighs; stack_i++)
            {
                TPZGeoElSide neigh = neighbourSet[stack_i];
                int neighMatID = neigh.Element()->MaterialId();
                TPZCompElSide celNeigh = neigh.Reference();

                int64_t neighIndex = neigh.Element()->Index();

                if (neighMatID != simData->DomainVec()[0].matID)
                    continue;

                if (neigh.Element()->HasSubElement())
                {
                    // Check if it is working in the case with refined meshes
                    DebugStop();
                }
                else
                {
                    TPZGeoElBC gbc(gelSide, simData->FluxInterfaceID());

                    TPZMultiphysicsInterfaceElement *interElem = new TPZMultiphysicsInterfaceElement(*cmesh_m, gbc.CreatedElement(), celNeigh, celSide);
                    interElem->SetLeftRightElementIndices(LeftElIndex, RightElIndex);
                }
            }
        }
    }
}

void TPZMeshOperator::InsertInterfaces(TPZMultiphysicsCompMesh *cmesh_m, ProblemData *simData, TPZGeoMesh *gmesh)
{
    TPZManVector<int, 2> Interfaces(2, 0);
    Interfaces[0] = simData->InterfaceID();
    Interfaces[1] = -simData->InterfaceID();

    int dim = cmesh_m->Dimension();

    if (!gmesh)
        DebugStop();

    int nInterfaceCreated = 0;

    int matfrom = simData->LambdaID();

    int64_t nel = gmesh->NElements();

    for (int64_t el = 0; el < nel; el++)
    {
        TPZGeoEl *gel = gmesh->Element(el);
        int meshdim = gmesh->Dimension();
        int matid = gel->MaterialId();

        if (matid != matfrom)
            continue;
        if (gel->HasSubElement())
            continue;

        int nsides = gel->NSides();
        TPZGeoElSide gelside(gel, nsides - 1);
        TPZCompElSide celside = gelside.Reference();

        TPZStack<TPZGeoElSide> neighbourSet;
        gelside.AllNeighbours(neighbourSet);

        gelside.LowerLevelCompElementList2(1);

        int64_t nneighs = neighbourSet.size();

        TPZManVector<int64_t, 3> LeftElIndices(1, 0.), RightElIndices(1, 0);
        LeftElIndices[0] = 0;
        RightElIndices[0] = 1;

        for (int stack_i = 0; stack_i < nneighs; stack_i++)
        {
            TPZGeoElSide neigh = neighbourSet[stack_i];
            int neighMatID = neigh.Element()->MaterialId();
            TPZCompElSide celneigh = neigh.Reference();

            if (!celside)
                DebugStop();

            if (neigh.Element()->HasSubElement())
            {
                DebugStop();
            }
            else
            {

                int64_t neighIndex = neigh.Element()->Index();

                if (neigh.Element()->Dimension() != meshdim)
                    continue;

                TPZGeoElBC gbc(gelside, Interfaces[stack_i]);

                TPZMultiphysicsInterfaceElement *interElem = new TPZMultiphysicsInterfaceElement(*cmesh_m, gbc.CreatedElement(), celneigh, celside);
                interElem->SetLeftRightElementIndices(LeftElIndices, RightElIndices);
                nInterfaceCreated++;
            }
        }
    }

    std::cout << __PRETTY_FUNCTION__ << "Number of Interfaces Created " << nInterfaceCreated << std::endl;
}

TPZCompMesh *TPZMeshOperator::CreateCMeshV(ProblemData *simData, TPZGeoMesh *gmesh)
{
    TPZCompMesh *cmesh_v = new TPZCompMesh(gmesh);
    cmesh_v->SetName("CMesh_V");
    
    std::set<int> materialIDs;

    // domain's material (2D or 3D)
    if (simData->DomainVec().size() != 0)
    {
        cmesh_v->SetDefaultOrder(simData->VelpOrder());

        cmesh_v->SetDimModel(simData->Dim());

        if (simData->HdivType() == EConstant)
        {
            cmesh_v->ApproxSpace().SetHDivFamily(HDivFamily::EHDivConstant);
        }
        else if (simData->HdivType() == EStandard)
        {
            cmesh_v->ApproxSpace().SetHDivFamily(HDivFamily::EHDivStandard);
        }

        cmesh_v->SetAllCreateFunctionsHDiv();

        auto *mat = new TPZNullMaterial<>(simData->DomainVec()[0].matID);
        cmesh_v->InsertMaterialObject(mat);

        materialIDs.insert(simData->DomainVec()[0].matID);

        // boundary conditions' material
        TPZFMatrix<STATE> val1(1, 1, 0.);
        TPZManVector<STATE> val2(1, 0.);

        for (const auto &bc : simData->NormalBCs())
        {
            val2 = bc.value;

            auto BCmat = mat->CreateBC(mat, bc.matID, bc.type, val1, val2);
            cmesh_v->InsertMaterialObject(BCmat);
            materialIDs.insert(bc.matID);
        }

        cmesh_v->AutoBuild(materialIDs);

        // setting the approximation order for the volume elements
        int64_t ncEl = cmesh_v->NElements();
        for (int64_t cEl = 0; cEl < ncEl; cEl++)
        {
            TPZCompEl *compEl = cmesh_v->Element(cEl);

            // only in those elements whose dimension equals to the simulation dim
            if (compEl->Dimension() == simData->Dim())
            {
                // dynamic casting the compEl object to use the ForceSideOrder function
                TPZInterpolatedElement *intercEl = dynamic_cast<TPZInterpolatedElement *>(compEl);

                // checking if the dynamic cast exists
                if (!intercEl)
                    continue;

                // finally using the desired function
                intercEl->ForceSideOrder(compEl->Reference()->NSides() - 1, simData->VelpOrder() + 1);
            }
        }

        int numEq2d = cmesh_v->NEquations();
        gmesh->ResetReference();
    }

    if (simData->Axisymmetric() && simData->AxisymmetryDomainVec().size() != 0) //in this case, we have a 1d mesh at r=0
    {
        cmesh_v->SetDefaultOrder(simData->VelpOrder());

        cmesh_v->SetDimModel(1);

        cmesh_v->ApproxSpace().SetHDivFamily(HDivFamily::EHDivStandard);

        cmesh_v->SetAllCreateFunctionsHDiv();

        // axisymmetry domain's material
        auto *Axisymmetricmat = new TPZNullMaterial<>(simData->AxisymmetryDomainVec()[0].matID);
        cmesh_v->InsertMaterialObject(Axisymmetricmat);

        materialIDs.clear();
        materialIDs.insert(simData->AxisymmetryDomainVec()[0].matID);

        // axisymmetry boundary conditions' material
        TPZFMatrix<STATE> val1(1, 1, 0.);
        TPZManVector<STATE> val2(1, 0.);
        for (const auto &bc : simData->AxisymmetryBCs())
        {
            val2 = bc.value;

            auto BCmat = Axisymmetricmat->CreateBC(Axisymmetricmat, bc.matID, bc.type, val1, val2);
            cmesh_v->InsertMaterialObject(BCmat);
            materialIDs.insert(bc.matID);
        }

        cmesh_v->AutoBuild(materialIDs);
        std::cout << cmesh_v->NElements() << std::endl;
        
        for (int64_t el = 0; el < cmesh_v->NElements(); el++)
        {
            auto cel = cmesh_v->Element(el);
            auto gel = cel->Reference();
            if (gel->Dimension() == 1 && cel->NConnects() == 3)
            {
                
                int64_t ncon = cel->NConnects();
                for (int64_t i = 0; i < ncon; i++)
                {
                    TPZConnect &newnod = cel->Connect(i);
                    newnod.SetLagrangeMultiplier(3);
                }
            }
        }
    }

    // expanding the solution vector
    cmesh_v->ExpandSolution();
    
    if ((simData->CondensedElements() && simData->MeshVector().size() != 4) || simData->MeshVector().size() < 2)
        DebugStop();

    simData->MeshVector()[0] = cmesh_v;

    CheckSideOrientOfCompEl(simData, gmesh);

    return cmesh_v;
}

TPZCompMesh *TPZMeshOperator::CreateCmeshP(ProblemData *simData, TPZGeoMesh *gmesh)
{
    TPZCompMesh *cmesh_p = new TPZCompMesh(gmesh);
    cmesh_p->SetName("CMesh_P");

    std::set<int> materialIDs;

    if (simData->DomainVec().size() != 0) //if true, there is a 2D/3D domain
    {
        cmesh_p->SetDimModel(simData->Dim());

        if (simData->HdivType() == EConstant)
        {
            cmesh_p->SetAllCreateFunctionsDiscontinuous();
            cmesh_p->SetDefaultOrder(0);
        }
        else if (simData->HdivType() == EStandard)
        {
            cmesh_p->SetDefaultOrder(simData->VelpOrder() + 1);
            cmesh_p->SetAllCreateFunctionsContinuous();
        }

        cmesh_p->ApproxSpace().CreateDisconnectedElements(true);

        // domain's material
        auto *mat = new TPZNullMaterial<>(simData->DomainVec()[0].matID);
        cmesh_p->InsertMaterialObject(mat);

        materialIDs.insert(simData->DomainVec()[0].matID);

        cmesh_p->AutoBuild(materialIDs);
        gmesh->ResetReference();

        materialIDs.clear();

        // matlambda traction material
        auto matLambda = new TPZNullMaterial<>(simData->LambdaID());
        matLambda->SetNStateVariables(simData->Dim() - 1); //In 3D, lambda has 2 state variables (one at each tangential direction)
        cmesh_p->InsertMaterialObject(matLambda);

        materialIDs.insert(simData->LambdaID());

        // traction on boundary material
        for (const auto &bc : simData->TangentialBCs())
        {
            auto matLambdaBC = new TPZNullMaterial<>(bc.matID);
            matLambdaBC->SetNStateVariables(simData->Dim() - 1);
            cmesh_p->InsertMaterialObject(matLambdaBC);

            materialIDs.insert(bc.matID);
        }

        if (simData->TracpOrder() > 0)
        {
            cmesh_p->SetAllCreateFunctionsContinuous();
            cmesh_p->ApproxSpace().CreateDisconnectedElements(true);
        }
        else
        {
            cmesh_p->SetAllCreateFunctionsDiscontinuous();
            cmesh_p->ApproxSpace().CreateDisconnectedElements(true);
        }

        cmesh_p->SetDefaultOrder(simData->TracpOrder());
        cmesh_p->SetDimModel(simData->Dim()-1);
        cmesh_p->AutoBuild(materialIDs);

        int64_t ncon = cmesh_p->NConnects();
        for (int64_t i = 0; i < ncon; i++)
        {
            TPZConnect &newnod = cmesh_p->ConnectVec()[i];
            newnod.SetLagrangeMultiplier(1);
        }

        //gmesh->ResetReference();
    }
    
    if (simData->Axisymmetric() && simData->AxisymmetryDomainVec().size() != 0) //in this case, we have a 1d mesh at r=0
    {
        cmesh_p->SetDimModel(1);

        cmesh_p->SetDefaultOrder(simData->VelpOrder());
        cmesh_p->SetAllCreateFunctionsContinuous();
        cmesh_p->ApproxSpace().CreateDisconnectedElements(true);

        // axisymmetry domain's material
        auto *Axisymmetricmat = new TPZNullMaterial<>(simData->AxisymmetryDomainVec()[0].matID);
        cmesh_p->InsertMaterialObject(Axisymmetricmat);

        materialIDs.clear();
        materialIDs.insert(simData->AxisymmetryDomainVec()[0].matID);

        cmesh_p->AutoBuild(materialIDs);

        int64_t ncon = cmesh_p->NConnects();
        for (int64_t i = 0; i < ncon; i++)
        {
            TPZConnect &newnod = cmesh_p->ConnectVec()[i];
            unsigned char lagrange = newnod.LagrangeMultiplier();
            unsigned char zero = '\000';
            if (lagrange == zero)
            {
                newnod.SetLagrangeMultiplier(2);
            }
        }

        gmesh->ResetReference();

        materialIDs.clear();

        // matlambda traction material
        auto matLambda = new TPZNullMaterial<>(simData->AxiLambdaID());
        cmesh_p->InsertMaterialObject(matLambda);
        materialIDs.insert(simData->AxiLambdaID());
        if (simData->TracpOrder() > 0)
        {
            cmesh_p->SetAllCreateFunctionsContinuous();
            cmesh_p->ApproxSpace().CreateDisconnectedElements(true);
        }
        else
        {
            cmesh_p->SetAllCreateFunctionsDiscontinuous();
            cmesh_p->ApproxSpace().CreateDisconnectedElements(true);
        }

        cmesh_p->SetDefaultOrder(simData->TracpOrder());
        cmesh_p->SetDimModel(1);
        cmesh_p->AutoBuild(materialIDs);

        ncon = cmesh_p->NConnects();
        for (int64_t i = 0; i < ncon; i++)
        {
            TPZConnect &newnod = cmesh_p->ConnectVec()[i];
            unsigned char lagrange = newnod.LagrangeMultiplier();
            unsigned char zero = '\000';
            if (lagrange == zero)
            {
                newnod.SetLagrangeMultiplier(1);
            }
        }
    }

    cmesh_p->ExpandSolution();

    simData->MeshVector()[1] = cmesh_p;

    return cmesh_p;
}

TPZCompMesh *TPZMeshOperator::CreateCmeshMv(ProblemData *simData, TPZGeoMesh *gmesh)
{
    TPZCompMesh *cmesh_Mv = new TPZCompMesh(gmesh);

    cmesh_Mv->SetName("CMesh_MV");
    cmesh_Mv->SetDefaultOrder(0);
    cmesh_Mv->SetDimModel(simData->Dim());

    cmesh_Mv->SetAllCreateFunctionsDiscontinuous();

    auto material_vM = new TPZNullMaterial<>(simData->DomainVec()[0].matID);
    cmesh_Mv->InsertMaterialObject(material_vM);

    int64_t ncel = cmesh_Mv->NElements();

    for (int64_t i = 0; i < ncel; i++)
    {
        TPZCompEl *compEl = cmesh_Mv->ElementVec()[i];

        if (!compEl)
            continue;

        TPZInterfaceElement *faceEl = dynamic_cast<TPZInterfaceElement *>(compEl);

        if (faceEl)
            DebugStop();
    }

    cmesh_Mv->AutoBuild();
    cmesh_Mv->AdjustBoundaryElements();
    cmesh_Mv->CleanUpUnconnectedNodes();

    simData->MeshVector()[3] = cmesh_Mv;

    return cmesh_Mv;
}

TPZCompMesh *TPZMeshOperator::CreateCmeshMp(ProblemData *simData, TPZGeoMesh *gmesh)
{
    TPZCompMesh *cmesh_Mp = new TPZCompMesh(gmesh);

    cmesh_Mp->SetName("CMesh_MP");
    cmesh_Mp->SetDefaultOrder(0);
    cmesh_Mp->SetDimModel(simData->Dim());

    cmesh_Mp->SetAllCreateFunctionsDiscontinuous();

    auto material_pM = new TPZNullMaterial<>(simData->DomainVec()[0].matID);
    cmesh_Mp->InsertMaterialObject(material_pM);

    int64_t ncel = cmesh_Mp->NElements();

    for (int64_t i = 0; i < ncel; i++)
    {
        TPZCompEl *compEl = cmesh_Mp->ElementVec()[i];

        if (!compEl)
            continue;

        TPZInterfaceElement *faceEl = dynamic_cast<TPZInterfaceElement *>(compEl);

        if (faceEl)
            DebugStop();
    }

    cmesh_Mp->AutoBuild();
    cmesh_Mp->AdjustBoundaryElements();
    cmesh_Mp->CleanUpUnconnectedNodes();

    simData->MeshVector()[2] = cmesh_Mp;

    return cmesh_Mp;
}

TPZMultiphysicsCompMesh *TPZMeshOperator::CreateMultiPhysicsMesh(ProblemData *simData, TPZGeoMesh *gmesh)
{
    TPZMultiphysicsCompMesh *cmesh_m = new TPZMultiphysicsCompMesh(gmesh);

    cmesh_m->SetName("CMesh_M");

    cmesh_m->SetDefaultOrder(simData->VelpOrder());
    cmesh_m->SetAllCreateFunctionsMultiphysicElem();

    // Creating Materials

    // 1. For domain
    if (simData->DomainVec().size() != 0)
    {
        TPZStokesMaterial *material = simData->Axisymmetric() ? new TPZAxisymStokesMaterial(simData->DomainVec()[0].matID, simData->Dim(), simData->DomainVec()[0].viscosity) : new TPZStokesMaterial(simData->DomainVec()[0].matID, simData->Dim(), simData->DomainVec()[0].viscosity);

        //TPZMixedLinearElasticMaterial* material = new TPZMixedLinearElasticMaterial(simData->DomainVec()[0].matID, simData->Dim(), 1.0, 0.0, AnalysisType::EPlaneStress);

        bool hasAnalyticSolution = true;
        if (hasAnalyticSolution && dynamic_cast<TPZStokesMaterial*>(material))
        {
            // TAxisymmetricStokesAnalytic *analyticSol = new TAxisymmetricStokesAnalytic();
            // analyticSol->fExactSol = TAxisymmetricStokesAnalytic::ESlidingCouetteFlow;
            // analyticSol->fL = 2.0;
            // analyticSol->fRe = 2.0;
            // analyticSol->fRi = 1.0;
            // analyticSol->fp = -1.0; 
            // analyticSol->fvel = 1.0;
            // analyticSol->fviscosity = simData->DomainVec()[0].viscosity;
            // material->SetExactSol(analyticSol->ExactSolution(), 0);
            // material->SetForcingFunction(analyticSol->ForceFunc(), 0);
        }
        else if (hasAnalyticSolution && dynamic_cast<TPZMixedLinearElasticMaterial*>(material))
        {
            if (simData->Dim() == 2)
            {
                TElasticity2DAnalytic *analyticSol = new TElasticity2DAnalytic();
                analyticSol->fProblemType = TElasticity2DAnalytic::EDefState::EBend;
                analyticSol->gE = 200000.0;
                analyticSol->gPoisson = 0.0;
                analyticSol->fPlaneStress = 1;
                material->SetExactSol(analyticSol->ExactSolution(), 3);
                material->SetForcingFunction(analyticSol->ForceFunc(), 3);
            }
            else if (simData->Dim() == 3)
            {
                TElasticity3DAnalytic *analyticSol = new TElasticity3DAnalytic();
                analyticSol->fProblemType = TElasticity3DAnalytic::EDefState::EBend;
                analyticSol->fE = 1.0;
                analyticSol->fPoisson = 0.0;
                material->SetExactSol(analyticSol->ExactSolution(), 3);
                material->SetForcingFunction(analyticSol->ForceFunc(), 3);
            }
            else
            {
                std::cout << "Dimension wrong!" << std::endl;
                DebugStop();
            }
        }

        cmesh_m->InsertMaterialObject(material);

        // 2. Boundary Conditions
        TPZFMatrix<STATE> val1(3, 3, 0.);
        TPZManVector<STATE> val2(3, 0.);

        for (const auto &bc : simData->NormalBCs())
        {
            val2 = bc.value;

            TPZBndCond *matBC = material->CreateBC(material, bc.matID, bc.type, val1, val2);
            auto matBC2 = dynamic_cast<TPZBndCondT<STATE>*>(matBC);
            if (material->HasExactSol())
                matBC2->SetForcingFunctionBC(material->ExactSol(),5);
            cmesh_m->InsertMaterialObject(matBC);
        }

        for (const auto &bc : simData->TangentialBCs())
        {
            val2 = bc.value;

            TPZBndCond *matBC = material->CreateBC(material, bc.matID, bc.type, val1, val2);
            auto matBC2 = dynamic_cast<TPZBndCondT<STATE>*>(matBC);
            if (material->HasExactSol())
                matBC2->SetForcingFunctionBC(material->ExactSol(),5);
            cmesh_m->InsertMaterialObject(matBC);
        }

        // 3 - Material for tangential traction
        TPZNullMaterialCS<> *matLambda = new TPZNullMaterialCS<>(simData->LambdaID());
        double dim = simData->Axisymmetric()? 1 : simData->Dim()-1;
        matLambda->SetDimension(dim);
        matLambda->SetNStateVariables(dim);
        cmesh_m->InsertMaterialObject(matLambda);

        // 4 - Material for interfaces (Inner)
        TPZInterfaceAxisymStokesMaterial *matInterfaceLeft = new TPZInterfaceAxisymStokesMaterial(simData->InterfaceID(), simData->Dim()-1);
        matInterfaceLeft->SetMultiplier(1.);
        cmesh_m->InsertMaterialObject(matInterfaceLeft);

        TPZInterfaceAxisymStokesMaterial *matInterfaceRight = new TPZInterfaceAxisymStokesMaterial(-simData->InterfaceID(), simData->Dim()-1);
        matInterfaceRight->SetMultiplier(-1.);
        cmesh_m->InsertMaterialObject(matInterfaceRight);
    }

    if (simData->Axisymmetric() && simData->AxisymmetryDomainVec().size() != 0) //in this case, we have 1d element at r=0
    {
        // 5. 1d Axisymmetric tube at r=0
        TPZ1dStokesMaterial *material = new TPZ1dStokesMaterial(simData->AxisymmetryDomainVec()[0].matID, simData->AxisymmetryDomainVec()[0].viscosity, simData->AxisymmetryDomainVec()[0].radius);
        cmesh_m->InsertMaterialObject(material);

        // 6. Boundary Conditions
        TPZFMatrix<STATE> val1(3, 3, 0.);
        TPZManVector<STATE> val2(3, 0.);

        for (const auto &bc : simData->AxisymmetryBCs())
        {
            val2 = bc.value;

            TPZBndCond *matBC = material->CreateBC(material, bc.matID, bc.type, val1, val2);
            cmesh_m->InsertMaterialObject(matBC);
        }

        // 7 - Material for tangential traction
        TPZNullMaterialCS<> *matLambda = new TPZNullMaterialCS<>(simData->AxiLambdaID());
        matLambda->SetDimension(1);
        matLambda->SetNStateVariables(1);
        cmesh_m->InsertMaterialObject(matLambda);

        // 8. - Material for 1d stokes interface
        TPZInterface1dStokesMaterial* mat1dInterface = new TPZInterface1dStokesMaterial(simData->AxiInterfaceID(), simData->AxisymmetryDomainVec()[0].viscosity, simData->AxisymmetryDomainVec()[0].radius);
        mat1dInterface->SetMultiplier(1.0);
        cmesh_m->InsertMaterialObject(mat1dInterface);

        // 9. - Material for 1d flux interface
        TPZInterface1dFlux* mat1dFluxInterface = new TPZInterface1dFlux(simData->FluxInterfaceID(), simData->AxisymmetryDomainVec()[0].viscosity, simData->AxisymmetryDomainVec()[0].radius);
        mat1dFluxInterface->SetMultiplier(1.0);
        cmesh_m->InsertMaterialObject(mat1dFluxInterface);
    }

    // Creating computational elements that will manage the mesh approximation space: Aparentemente não faz nada????????
    int64_t ncel = cmesh_m->NElements();
    for (int i = 0; i < ncel; i++)
    {
        TPZCompEl *compEl = cmesh_m->ElementVec()[i];
        if (!compEl)
            continue;
        TPZInterfaceElement *facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if (facel)
            DebugStop();
    }

    TPZManVector<int, 2> active_approx_spaces(simData->MeshVector().size(), 1);

    cmesh_m->BuildMultiphysicsSpace(active_approx_spaces, simData->MeshVector());
    cmesh_m->AdjustBoundaryElements();
    cmesh_m->CleanUpUnconnectedNodes();
    cmesh_m->LoadReferences();

    TPZMeshOperator::InsertBCInterfaces(cmesh_m, simData, gmesh);
    TPZMeshOperator::InsertInterfaces(cmesh_m, simData, gmesh);

    if (simData->CondensedElements())
    {
        cmesh_m->SetName("CMesh_M_BeforeCond");
        cmesh_m->ComputeNodElCon();
        PrintCompMesh(cmesh_m);
    }

    return cmesh_m;
}

void TPZMeshOperator::CondenseElements(ProblemData *simData, TPZMultiphysicsCompMesh *cmesh_m)
{
    int64_t ncompEl = cmesh_m->ElementVec().NElements();
    int dim = cmesh_m->Reference()->Dimension();

    std::set<int64_t> externalNode;
    std::vector<int64_t> groupIndex;
    TPZStack<TPZElementGroup *> elGroups;
    int count = 0;

    int domain1d_matid = simData->AxisymmetryDomainVec()[0].matID;
    int lambda1d_matid = simData->AxiLambdaID();
    int bc1d_matid[2] = {simData->AxisymmetryBCs()[0].matID, simData->AxisymmetryBCs()[1].matID};

    auto cmesh_v = cmesh_m->MeshVector()[0];

    // Creating the element groups for the domain
    for (int64_t el = 0; el < ncompEl; el++)
    {
        TPZCompEl *compEl = cmesh_m->Element(el);

        if (compEl->Dimension() != dim)
            continue;

        TPZMultiphysicsElement *multEl = dynamic_cast<TPZMultiphysicsElement *>(compEl);
        int64_t numSpaces = multEl->NMeshes();

        if (numSpaces < 4)
            DebugStop();

        int64_t numConnectExt = numSpaces - 3;
        int nConnect = multEl->NConnects();

        for (int ic = nConnect - numConnectExt; ic < nConnect; ic++)
        {
            int64_t conIndex = compEl->ConnectIndex(ic);
            externalNode.insert(conIndex);
        }

        count++;
        groupIndex.resize(count);
        groupIndex[count - 1] = compEl->Index();

        TPZElementGroup *groupEl = new TPZElementGroup(*cmesh_m);
        elGroups.Push(groupEl);
        elGroups[count - 1]->AddElement(compEl);
    }

    //Creating element group for axisymmetric tube
    // if (simData->Axisymmetric() && simData->AxisymmetryDomainVec().size() != 0)
    // {
    //     for (int64_t el = 0; el < ncompEl; el++)
    //     {
    //         TPZCompEl *compEl = cmesh_m->Element(el);

    //         if (!compEl) continue;

    //         TPZGeoEl* gel = compEl->Reference();
    //         int matid = gel->MaterialId();
            

    //         if (matid == domain1d_matid || matid == lambda1d_matid || matid == bc1d_matid[0] || matid == bc1d_matid[1])
    //         {
    //             TPZMultiphysicsElement *multEl = dynamic_cast<TPZMultiphysicsElement *>(compEl);
    //             int64_t numSpaces = multEl->NMeshes();

    //             int nConnect = multEl->NConnects();

    //             if (matid == bc1d_matid[1]) //Not condensing 1 pressure dof
    //             {
    //                 int64_t conIndex = compEl->ConnectIndex(0);
    //                 externalNode.insert(conIndex);
    //             }

    //             count++;
    //             groupIndex.resize(count);
    //             groupIndex[count - 1] = compEl->Index();

    //             TPZElementGroup *groupEl = new TPZElementGroup(*cmesh_m);
    //             elGroups.Push(groupEl);
    //             elGroups[count - 1]->AddElement(compEl);
    //         }
    //     }
    // }

    // Inserting interfaces and boundary conditions

    for (int64_t el = 0; el < ncompEl; el++)
    {
        TPZCompEl *compEl = cmesh_m->Element(el);

        TPZMultiphysicsInterfaceElement *interEl = dynamic_cast<TPZMultiphysicsInterfaceElement *>(compEl);

        if (interEl)
        {
            TPZCompEl *leftEl = interEl->LeftElement();
            TPZGeoEl *leftGel = leftEl->Reference();
            int left_matid = leftGel->MaterialId();

            if (leftEl->Dimension() != dim)
                continue;

            int64_t leftIndex = leftEl->Index();
            for (int64_t iEl = 0; iEl < groupIndex.size(); iEl++)
            {
                if (leftIndex == groupIndex[iEl])
                {
                    elGroups[iEl]->AddElement(compEl);
                }
            }
        }

        if (!compEl)
            continue;
        //
        //        TPZGeoEl* geoEl = compEl->Reference();
        //
        //        if(geoEl->Dimension()==dim-1){
        //            TPZBndCond* elBC = dynamic_cast<TPZBndCond*>(compEl->Material());
        //            if (!elBC) continue;
        //
        //            TPZStack<TPZCompElSide> compElStack;
        //            TPZGeoElSide geoElSide(geoEl, geoEl->NSides()-1);
        //
        //            geoElSide.EqualLevelCompElementList(compElStack, 0, 0);
        //
        //            for(auto& compElStackIndex: compElStack){
        //                if(compElStackIndex.Reference().Element()->Dimension()==dim){
        //                    int64_t IndexBC = compElStackIndex.Element()->Index();
        //
        //                    for(int64_t iEl=0; iEl<groupIndex.size(); iEl++){
        //                        if(IndexBC==groupIndex[iEl]){
        //                            elGroups[iEl]->AddElement(compEl);
        //                        }
        //                    }
        //                }
        //            }
        //        }
    }

    cmesh_m->ComputeNodElCon();

    for (auto it = externalNode.begin(); it != externalNode.end(); it++)
    {
        int64_t coIndex = *it;
        cmesh_m->ConnectVec()[coIndex].IncrementElConnected();
    }

    //     Creating  condensed elements
    int64_t nenvel = elGroups.NElements();
    for (int64_t iEnv = 0; iEnv < nenvel; iEnv++)
    {
        TPZElementGroup *elGroup = elGroups[iEnv];
        new TPZCondensedCompElT<STATE>(elGroup);
    }

    cmesh_m->SetName("CMesh_M_Condensed");

    cmesh_m->CleanUpUnconnectedNodes();
    cmesh_m->ExpandSolution();
}

void TPZMeshOperator::PrintGeoMesh(TPZGeoMesh *gmesh)
{
    std::cout << "\nPrinting geometric mesh in .txt and .vtk formats...\n";

    std::ofstream VTKGeoMeshFile(gmesh->Name() + ".vtk");
    std::ofstream TextGeoMeshFile(gmesh->Name() + ".txt");

    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, VTKGeoMeshFile);
    gmesh->Print(TextGeoMeshFile);
}

void TPZMeshOperator::PrintCompMesh(TPZCompMesh *cmesh)
{
    std::cout << "\nPrinting multiphysics mesh in .txt and .vtk formats...\n";

    std::ofstream VTKCompMeshFile(cmesh->Name() + ".vtk");
    std::ofstream TextCompMeshFile(cmesh->Name() + ".txt");

    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, VTKCompMeshFile);
    cmesh->Print(TextCompMeshFile);
}

void TPZMeshOperator::PrintCompMesh(TPZVec<TPZCompMesh *> CMeshVec)
{
    std::cout << "\nPrinting computational meshes in .txt and .vtk formats...\n";

    for (auto &cmesh : CMeshVec)
    {

        std::ofstream VTKCompMeshFile(cmesh->Name() + ".vtk");
        std::ofstream TextCompMeshFile(cmesh->Name() + ".txt");

        TPZVTKGeoMesh::PrintCMeshVTK(cmesh, VTKCompMeshFile);
        cmesh->ComputeNodElCon();
        cmesh->Print(TextCompMeshFile);
    }
}

void TPZMeshOperator::CheckSideOrientOfCompEl(ProblemData* simData, TPZGeoMesh* gmesh)
{
    if (simData->AxisymmetryDomainVec().size() != 0)
    {
        std::set<int> bcIds;
        for(auto& bc : simData->AxisymmetryBCs())
            bcIds.insert(bc.matID);

        int matid1d = simData->AxisymmetryDomainVec()[0].matID;

        for(TPZGeoEl* gel : gmesh->ElementVec())
        {
            const int gelmatid = gel->MaterialId();

            if (gelmatid != matid1d) //Check only 1d elements
                continue;

            int nside = gel->NSides(0); //Check only the vertices
            
            TPZInterpolatedElement* intel = dynamic_cast<TPZInterpolatedElement*>(gel->Reference());
            if (!intel) DebugStop();

            for (int is = 0; is < nside; is++)
            {
                const int sideorientgel = intel->GetSideOrient(is);
                TPZGeoElSide gelside(gel,is);
                //TPZGeoElSide neig = gelside.HasNeighbour(matid1d);
                TPZGeoElSide neig = gelside.Neighbour();

                for (; neig != gelside; neig++)
                {
                    int neigmatid = neig.Element()->MaterialId();

                    if (neigmatid == matid1d) //vertex is has an 1d element as neighbour
                    {
                        TPZInterpolatedElement* intelneig = dynamic_cast<TPZInterpolatedElement*>(neig.Element()->Reference());
                        const int sideorientneig = intelneig->GetSideOrient(neig.Side());

                        if (sideorientgel * sideorientneig != -1)
                        {
                            intelneig->SetSideOrient(neig.Side(), -sideorientgel);
                        }
                    }
                    else if (bcIds.find(neigmatid) != bcIds.end()) //vertex is has a bc as neighbour
                    {
                        intel->SetSideOrient(is, 1);
                    }
                }
                
                // if (neig && (gelside != neig))
                // {
                //     TPZInterpolatedElement* intelneig = dynamic_cast<TPZInterpolatedElement*>(neig.Element()->Reference());
                //     const int sideorientneig = intelneig->GetSideOrient(neig.Side());

                //     if (sideorientgel * sideorientneig != -1)
                //     {
                //         intelneig->SetSideOrient(neig.Side(), -sideorientgel);
                //     }
                // }

                // neig = gelside.HasNeighbour(bcIds);
                // if (neig)
                // {
                //     intel->SetSideOrient(is, 1);
                // }
            }
        }
    }
}
