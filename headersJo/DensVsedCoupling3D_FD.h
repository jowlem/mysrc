/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2017 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/* Main author: Orestis Malaspinas
 */

#ifndef DENS_VSED_PROCESSOR_3D_H
#define DENS_VSED_PROCESSOR_3D_H

#include "core/globalDefs.h"
#include "core/block3D.h"
#include "atomicBlock/dataProcessor3D.h"
#include "atomicBlock/blockLattice3D.h"

namespace plb {


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Jonathan
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template< typename T1, typename T2  >
class DensVsedProcessor3D :
    public BoxProcessingFunctional3D_SS<T1, T2>
{
public:

    DensVsedProcessor3D(T1 rhoP_, T1 convers_);

    virtual void process( Box3D domain,
                          ScalarField3D<T1>& v_sed,
                          ScalarField3D<T2>& Density );
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
        modified[1] = modif::nothing;
    }
    virtual DensVsedProcessor3D<T1, T2>* clone() const;

private:
  T1 rhoP; T1 convers;

};


/*
template< typename T,
          template<typename U1> class FluidDescriptor,
          template<typename U2> class DensityDescriptor
        >
class SmagorinskyBoussinesqThermalProcessor3D :
    public BoxProcessingFunctional3D_LL<T,FluidDescriptor,T,DensityDescriptor>
{
public:
    SmagorinskyBoussinesqThermalProcessor3D (
            T gravity_, T T0_, T deltaTemp_, Array<T,FluidDescriptor<T>::d> dir_,
            T cSmagoFluid_, T cSmagoTemp_ )

    virtual void process( Box3D domain,
                          BlockLattice3D<T,FluidDescriptor>& fluid,
                          BlockLattice3D<T,DensityDescriptor>& temperature );
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
        modified[1] = modif::staticVariables;
    }
    virtual SmagorinskyBoussinesqThermalProcessor3D<T,FluidDescriptor,DensityDescriptor>* clone() const;
private:
    T gravity, T0, deltaTemp;
    Array<T,FluidDescriptor<T>::d> dir;
    T cSmagoFluid, T cSmagoTemp;
};
*/

}  // namespace plb

#endif  // BOUSSINESQ_THERMAL_PROCESSOR_3D_H
