/**
 * @file
 * @brief Contains the TPZDarcy2Spaces material class that has 1 H1 space and 1 constant space
 * @author Giovane Avancini
 * @date 22/09/2023
 */

#ifndef TPZDARCY2SPACES_H
#define TPZDARCY2SPACES_H

#include "DarcyFlow/TPZHybridDarcyFlow.h"

class TPZDarcy2Spaces : public TPZHybridDarcyFlow
{

public:
    /**
     * @brief Default constructor
     */
    TPZDarcy2Spaces();

    /**
	 * @brief Class constructor
	 * @param [in] id material id
	 * @param [in] dim problem dimension
	 */
    TPZDarcy2Spaces(int id, int dim);

    TPZDarcy2Spaces(const TPZDarcy2Spaces &copy) : TPZHybridDarcyFlow(copy) {}

    TPZDarcy2Spaces &operator=(const TPZDarcy2Spaces &copy)
    {
        TPZHybridDarcyFlow::operator=(copy);
        return *this;
    }
    /**
     * @brief Returns the problem dimension
     */
    [[nodiscard]] int Dimension() const override {
        return this->fDim;
    }

    /**
     * @brief Returns the number of state variables
     */
    [[nodiscard]] int NStateVariables() const override { return 1; }


    /**
	 * @brief Returns a 'std::string' with the name of the material
	 */
    [[nodiscard]] std::string Name() const override { return "TPZDarcy2Spaces"; }

    virtual int NEvalErrors()  const override {return 4;}

    /** @name Contribute */
    /** @{ */
    /**
     * @brief It computes a contribution to the stiffness matrix
     * and load vector at one integration point.
     * @param[in] datavec stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ek is the element matrix
     * @param[out] ef is the rhs vector
     */
    virtual void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                            REAL weight,TPZFMatrix<STATE> &ek,
                            TPZFMatrix<STATE> &ef) override;
    /**@}*/

    /** @name ContributeBC
        @ingroup Contribute*/
    /**@{*/
    /**
     * @brief It computes a contribution to the stiffness matrix
     * and load vector at one BC integration point.
     * @param[in] datavec stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ek is the element matrix
     * @param[out] ef is the rhs vector
     * @param[in] bc is the boundary condition material
     */
    virtual void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                              REAL weight, TPZFMatrix<STATE> &ek,
                              TPZFMatrix<STATE> &ef,
                              TPZBndCondT<STATE> &bc) override;

    /**@}*/
    /** @brief Returns the solution associated with a given index
        based on the finite element approximation.
        @param[in] datavec Stores all the input data.
        @param[in] var Index of the queried solution
        @param[out] sol FEM Solution at the integration point
    */
    virtual void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                          int var, TPZVec<STATE> &sol) override {}
    /**
     * @brief Returns an integer associated with a post-processing variable name
     * @param [in] name string containing the name of the post-processing variable. Ex: "Pressure".
     */
    [[nodiscard]] int VariableIndex(const std::string &name) const override;

    /**
     * @brief Returns an integer with the dimension of a post-processing variable
     * @param [in] var index of the post-processing variable, according to TPZDarcyFlow::VariableIndex method.
     */
    [[nodiscard]] int NSolutionVariables(int var) const override;

    //! @name Error
    /** @{*/
    /*!
      \brief Calculates the error at a given point x.
      \param[in] datavec input data
      \param[out] errors The calculated errors.
     */
    virtual void Errors(const TPZVec<TPZMaterialDataT<STATE>> &data,
                        TPZVec<REAL> &errors) override;

    /**
     * @brief Returns an unique class identifier
     */
    [[nodiscard]] int ClassId() const override;
    
    /** @brief Writes this object to the TPZStream buffer. Include the classid if `withclassid = true` */
    virtual void Write(TPZStream &buf, int withclassid) const override
    {
        
    }
    
    /** @brief Reads an objects from the TPZStream buffer. */
    virtual void Read(TPZStream &buf, void *context) override
    {
        
    }


    /**
     * @brief Creates another material of the same type
     */
    [[nodiscard]] TPZMaterial *NewMaterial() const override;

    /**
     * @brief Prints data associated with the material.
     */
    void Print(std::ostream & out) const override;

    /** @brief Creates an associated boundary condition.
     @param[in] reference The volumetric material associated with the BC.
     @param[in] id Boundary condition identifier.
     @param[in] type Type of the boundary condition.
     @param[in] val1 Value to be set at the element matrix.
     @param[in] val2 Value to be set at the rhs vector.
    */
    // virtual TPZBndCondT<STATE>* CreateBC(TPZMaterial *reference,
    //                                     int id, int type,
    //                                     const TPZFMatrix<STATE> &val1,
    //                                     const TPZVec<STATE> &val2) override
    // {
    //     return TPZDarcyFlow::CreateBC(reference,id,type,val1,val2);
    // }

};


#endif //TPZDARCYFLOW_H