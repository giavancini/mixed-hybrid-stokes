#include <pzfmatrix.h>
#include <TPZMaterial.h>
#include <TPZMatBase.h>
#include <TPZMaterialData.h>
#include <TPZMaterialDataT.h>
#include <TPZMatCombinedSpaces.h>
#include <TPZMatInterfaceCombinedSpaces.h>
#include <TPZLagrangeMultiplier.h>
#include <math.h>

#ifndef TPZINTERFACEMATERIAL_H
#define TPZINTERFACEMATERIAL_H

class TPZInterfaceMaterial : public TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>, TPZMatInterfaceCombinedSpaces<STATE>>,
    
    public TPZLagrangeMultiplierBase
    
{
    
    using TBase = TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>,TPZMatInterfaceCombinedSpaces<STATE>>;
    
    bool IsLagrangeMult() override{return true;};
    
protected:
    int fdimension;
    int fVindex = 1;
    int fPindex = 0;
    int fNStateVariables = -1;
    
    REAL fBigNumber = pow(10,std::numeric_limits<STATE>::max_digits10*2/3);

public:
    /// Creates a material object
    TPZInterfaceMaterial(int matID, int dimension);
    
    /// Destructor
    ~TPZInterfaceMaterial();
    
    int Dimension() const override {return fdimension;}
    
    void ContributeInterface(const TPZMaterialDataT<STATE>& data,
                             const std::map<int, TPZMaterialDataT<STATE>>& dataleft,
                             const std::map<int, TPZMaterialDataT<STATE>>& dataright, REAL
                             weight,
                             TPZFMatrix<STATE>& ek,
                             TPZFMatrix<STATE>& ef) override;
    
    void ContributeBCInterface(const TPZMaterialDataT<STATE> &data,
                               const std::map<int, TPZMaterialDataT<STATE>> &dataleft,
                               REAL weight,
                               TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef,
                               TPZBndCondT<STATE> &bc) override;
    
    void FillDataRequirementsInterface(TPZMaterialDataT<STATE> &data,
                                       std::map<int, TPZMaterialDataT<STATE>> &datavec_left,
                                       std::map<int, TPZMaterialDataT<STATE>> &datavec_right) override;
    
    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                    REAL weight,
                    TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;

    void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                      REAL weight, TPZFMatrix<STATE> &ek,
                      TPZFMatrix<STATE> &ef,
                      TPZBndCondT<STATE> &bc) override;
    
    void SolutionInterface(const TPZMaterialDataT<STATE> &data,
                           const std::map<int, TPZMaterialDataT<STATE>> &dataleftvec,
                           const std::map<int, TPZMaterialDataT<STATE>> &datarightvec,
                           int var, TPZVec<STATE> &Solout) override;


    void SolutionInterface(const TPZMaterialDataT<STATE> &data,
                           const std::map<int, TPZMaterialDataT<STATE>> &dataleftvec,
                           const std::map<int, TPZMaterialDataT<STATE>> &datarightvec,
                           int var, TPZVec<STATE> &Solout,
                           TPZCompEl *left,TPZCompEl *right) override;
    
    void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                  int var, TPZVec<STATE> &sol) override;

    int NStateVariables() const override
    {return fNStateVariables;}
    
};
#endif
