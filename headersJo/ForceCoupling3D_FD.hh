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

#ifndef BUOYANT_THERMAL_PROCESSOR_3D_FD_HH
#define BUOYANT_THERMAL_PROCESSOR_3D_FD_HH

#include "ForceCoupling3D_FD.h"
#include "atomicBlock/blockLattice3D.h"
#include "core/util.h"
#include "finiteDifference/finiteDifference3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

namespace plb {


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Jonathan
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template< typename T, template<typename U> class FluidDescriptor>
ScalarBuoyanTermProcessor3D<T,FluidDescriptor>::
        ScalarBuoyanTermProcessor3D(T gravity_, T rho0_, T rhoP_, Array<T,FluidDescriptor<T>::d> dir_)
    :  gravity(gravity_), rho0(rho0_), rhoP(rhoP_),
       dir(dir_)
{
    // We normalize the direction of the force vector.
    T normDir = std::sqrt(VectorTemplate<T,FluidDescriptor>::normSqr(dir));
    for (pluint iD = 0; iD < FluidDescriptor<T>::d; ++iD) {
        dir[iD] /= normDir;
    }
}


template< typename T, template<typename U> class FluidDescriptor >
void ScalarBuoyanTermProcessor3D<T,FluidDescriptor>::process (
        Box3D domain,
        BlockLattice3D<T,FluidDescriptor>& fluid,
        ScalarField3D<T>& volfracfield )
{
    typedef FluidDescriptor<T> D;
    enum {
        forceOffset = FluidDescriptor<T>::ExternalField::forceBeginsAt
    };
    Dot3D offset = computeRelativeDisplacement(fluid, volfracfield);

    Array<T,D::d> gravOverrho0 (
            gravity*dir[0]/rho0,
            gravity*dir[1]/rho0,
            gravity*dir[2]/rho0 );


    for (plint iX=domain.x0; iX<=domain.x1; ++iX)
    {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY)
        {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ)
            {

                T localVolfrac = volfracfield.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                // Computation of the Boussinesq force
                T *force = fluid.get(iX,iY,iZ).getExternal(forceOffset);
                // volfracfield is the order-0 moment of the advection-diffusion lattice.
                const T diffT = rhoP-rho0;
                for (pluint iD = 0; iD < D::d; ++iD)
                {
			                 force[iD] = -localVolfrac * gravOverrho0[iD] * diffT * 0.002;
                }
            }
        }
    }
}

template< typename T, template<typename U> class FluidDescriptor>
ScalarBuoyanTermProcessor3D<T,FluidDescriptor>*
    ScalarBuoyanTermProcessor3D<T,FluidDescriptor>::clone() const
{
    return new ScalarBuoyanTermProcessor3D<T,FluidDescriptor>(*this);
}


/*
template< typename T,
          template<typename U1> class FluidDescriptor,
          template<typename U2> class volfracfieldDescriptor
        >
SmagorinskyBoussinesqThermalProcessor3D<T,FluidDescriptor,volfracfieldDescriptor>::
        SmagorinskyBoussinesqThermalProcessor3D(T gravity_, T T0_, T deltaTemp_, Array<T,FluidDescriptor<T>::d> dir_,
        T cSmagoFluid_, T cSmagoTemp_ )
    :  gravity(gravity_), T0(T0_), deltaTemp(deltaTemp_),
       dir(dir_),
       cSmagoFluid(cSmagoFluid_),
       cSmagoTemp(cSmagoTemp_),
       nu0(nu0_), d0(d0_)
{
    // We normalize the direction of the force vector.
    T normDir = std::sqrt(VectorTemplate<T,FluidDescriptor>::normSqr(dir));
    for (pl2uint iD = 0; iD < FluidDescriptor<T>::d; ++iD) {
        dir[iD] /= normDir;
    }
}

template< typename T,
          template<typename U1> class FluidDescriptor,
          template<typename U2> class volfracfieldDescriptor
        >
void SmagorinskyBoussinesqThermalProcessor3D<T,FluidDescriptor,volfracfieldDescriptor>::process (
        Box3D domain,
        BlockLattice3D<T,FluidDescriptor>& fluid,
        BlockLattice3D<T,volfracfieldDescriptor>& volfracfield )
{
    typedef FluidDescriptor<T> D;
    enum {
        velOffset   = volfracfieldDescriptor<T>::ExternalField::velocityBeginsAt,
        forceOffset = FluidDescriptor<T>::ExternalField::forceBeginsAt
    };
    Dot3D offset = computeRelativeDisplacement(fluid, volfracfield);

    Array<T,D::d> gravOverDeltaTemp (
            gravity*dir[0]/deltaTemp,
            gravity*dir[1]/deltaTemp,
            gravity*dir[2]/deltaTemp );

    for (plint iX=domain.x0; iX<=domain.x1; ++iX)
    {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY)
        {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ)
            {
                // Velocity coupling
                T *u = volfracfield.get(iX+offset.x,iY+offset.y,iZ+offset.z).getExternal(velOffset);
                Array<T,FluidDescriptor<T>::d> vel;
                fluid.get(iX,iY,iZ).computeVelocity(vel);
                vel.to_cArray(u);

                // Computation of the Boussinesq force
                T *force = fluid.get(iX,iY,iZ).getExternal(forceOffset);
                // volfracfield is the order-0 moment of the advection-diffusion lattice.
                //   You can compute it with the method computeDensity().
                T localvolfracfield = volfracfield.get(iX+offset.x,iY+offset.y,iZ+offset.z).computeDensity();
                const T diffT = localvolfracfield - T0;
                for (pluint iD = 0; iD < D::d; ++iD)
                {
                    force[iD] = gravOverDeltaTemp[iD] * diffT;
                }
            }
        }
    }
}

template< typename T,
          template<typename U1> class FluidDescriptor,
          template<typename U2> class volfracfieldDescriptor
        >
SmagorinskyBoussinesqThermalProcessor3D<T,FluidDescriptor,volfracfieldDescriptor>*
    SmagorinskyBoussinesqThermalProcessor3D<T,FluidDescriptor,volfracfieldDescriptor>::clone() const
{
    return new SmagorinskyBoussinesqThermalProcessor3D<T,FluidDescriptor,volfracfieldDescriptor>(*this);
}

*/

}  // namespace plb

#endif  // BOUSSINESQ_THERMAL_PROCESSOR_3D_HH
