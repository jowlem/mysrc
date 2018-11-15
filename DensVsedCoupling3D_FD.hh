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

#ifndef DENS_VSED_PROCESSOR_3D_HH
#define DENS_VSED_PROCESSOR_3D_HH

#include "DensVsedCoupling3D_FD.h"
#include "atomicBlock/blockLattice3D.h"
#include "core/util.h"
#include "finiteDifference/finiteDifference3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

namespace plb {


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Jonathan
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template< typename T1, typename T2 >
DensVsedProcessor3D<T1, T2>::DensVsedProcessor3D(T1 rhoP_, T1 convers_)
        : rhoP(rhoP_), convers(convers_)
{

}


template< typename T1, typename T2 >
void DensVsedProcessor3D<T1, T2>::process (
        Box3D domain,
        ScalarField3D<T1>& v_sed,
        ScalarField3D<T2>& Density )
{

    Dot3D offset = computeRelativeDisplacement(v_sed, Density);
    T1 g=9.81;
    T1 d=63e-6;
    T1 mu=1e-3;


    for (plint iX=domain.x0; iX<=domain.x1; ++iX)
    {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY)
        {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ)
            {

                T1 dens = Density.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                T1 vel_sed = v_sed.get(iX+offset.x,iY+offset.y,iZ+offset.z);

                if(vel_sed>0) vel_sed=convers*(0.5*d*d*g*(rhoP-dens))/(9*mu);
                else vel_sed=0;

                v_sed.get(iX+offset.x,iY+offset.y,iZ+offset.z)=vel_sed;


            }
        }
    }
}

template< typename T1, typename T2 >
DensVsedProcessor3D<T1,T2>*
    DensVsedProcessor3D<T1,T2>::clone() const
{
    return new DensVsedProcessor3D<T1,T2>(*this);
}


/*
template< typename T,
          template<typename U1> class FluidDescriptor,
          template<typename U2> class DensityDescriptor
        >
SmagorinskyBoussinesqThermalProcessor3D<T,FluidDescriptor,DensityDescriptor>::
        SmagorinskyBoussinesqThermalProceArray<T,FluidDescriptor<T>::d> dir;ssor3D(T gravity_, T T0_, T deltaTemp_, Array<T,FluidDescriptor<T>::d> dir_,
        T cSmagoFluid_, T cSmagoTemp_ )
    :  gravity(gravity_), T0(T0_), deltaTemp(deltaTemp_),
       dir(dir_),
       cSmagoFluid(cSmagoFluid_),
       cSmagoTemp(cSmagoTemp_),
       nu0(nu0_), d0(d0_)
{
    // We normalize the direction of the force vector.
    T normDir = std::sqrt(VectorTemplate<T,FluidDescriptor>::normSqr(dir));
    for (pluint iD = 0; iD < FluidDescriptor<T>::d; ++iD) {
        dir[iD] /= normDir;
    }
}

template< typename T,
          template<typename U1> class FluidDescriptor,
          template<typename U2> class DensityDescriptor
        >
void SmagorinskyBoussinesqThermalProcessor3D<T,FluidDescriptor,DensityDescriptor>::process (
        Box3D domain,
        BlockLattice3D<T,FluidDescriptor>& fluid,
        BlockLattice3D<T,DensityDescriptor>& Density )
{
    typedef FluidDescriptor<T> D;
    enum {
        velOffset   = DensityDescriptor<T>::ExternalField::velocityBeginsAt,
        forceOffset = FluidDescriptor<T>::ExternalField::forceBeginsAt
    };
    Dot3D offset = computeRelativeDisplacement(fluid, Density);

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
                T *u = Density.get(iX+offset.x,iY+offset.y,iZ+offset.z).getExternal(velOffset);
                Array<T,FluidDescriptor<T>::d> vel;
                fluid.get(iX,iY,iZ).computeVelocity(vel);
                vel.to_cArray(u);

                // Computation of the Boussinesq force
                T *force = fluid.get(iX,iY,iZ).getExternal(forceOffset);
                // Density is the order-0 moment of the advection-diffusion lattice.
                //   You can compute it with the method computeDensity().
                T localDensity = Density.get(iX+offset.x,iY+offset.y,iZ+offset.z).computeDensity();
                const T diffT = localDensity - T0;
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
          template<typename U2> class DensityDescriptor
        >
SmagorinskyBoussinesqThermalProcessor3D<T,FluidDescriptor,DensityDescriptor>*
    SmagorinskyBoussinesqThermalProcessor3D<T,FluidDescriptor,DensityDescriptor>::clone() const
{
    return new SmagorinskyBoussinesqThermalProcessor3D<T,FluidDescriptor,DensityDescriptor>(*this);
}

*/

}  // namespace plb

#endif  // BOUSSINESQ_THERMAL_PROCESSOR_3D_HH
