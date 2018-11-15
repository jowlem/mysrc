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

#ifndef VEL_COUPLING_PROCESSOR_3D_HH
#define VEL_COUPLING_PROCESSOR_3D_HH

#include "Velocity_Coupling_FD.h"
#include "atomicBlock/blockLattice3D.h"
#include "core/util.h"
#include "finiteDifference/finiteDifference3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

namespace plb {


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Jonathan
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template< typename T1, typename T2, int nDim >
Velocity_CouplingProcessor3D<T1, T2, nDim>::Velocity_CouplingProcessor3D()

{
}


template< typename T1, typename T2, int nDim >
void Velocity_CouplingProcessor3D<T1, T2, nDim>::process (
        Box3D domain,
        ScalarField3D<T1>& v_sed,
        TensorField3D<T2, nDim>& velocity )
{
    //typedef DensityDescriptor<T> D;

    Dot3D offset = computeRelativeDisplacement(v_sed, velocity);


    for (plint iX=domain.x0; iX<=domain.x1; ++iX)
    {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY)
        {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ)
            {



                T1 vel_sed = v_sed.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                Array<T1, 3> vel=velocity.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                vel[2] -= vel_sed;
                velocity.get(iX+offset.x,iY+offset.y,iZ+offset.z)=vel;



            }
        }
    }
}

template< typename T1, typename T2, int nDim >
Velocity_CouplingProcessor3D<T1, T2, nDim>*
    Velocity_CouplingProcessor3D<T1, T2, nDim>::clone() const
{
    return new Velocity_CouplingProcessor3D<T1, T2, nDim>(*this);
}


/*
template< typename T,
          template<typename U1> class DensityDescriptor,
          template<typename U2> class VolFracDescriptor
        >
SmagorinskyBoussinesqThermalProcessor3D<T,DensityDescriptor,VolFracDescriptor>::
        SmagorinskyBoussinesqThermalProcessor3D(T gravity_, T T0_, T deltaTemp_, Array<T,DensityDescriptor<T>::d> dir_,
        T cSmagoFluid_, T cSmagoTemp_ )
    :  gravity(gravity_), T0(T0_), deltaTemp(deltaTemp_),
       dir(dir_),
       cSmagoFluid(cSmagoFluid_),
       cSmagoTemp(cSmagoTemp_),
       nu0(nu0_), d0(d0_)
{
    // We normalize the direction of the force vector.
    T normDir = std::sqrt(VectorTemplate<T,DensityDescriptor>::normSqr(dir));
    for (pluint iD = 0; iD < DensityDescriptor<T>::d; ++iD) {
        dir[iD] /= normDir;
    }
}

template< typename T,
          template<typename U1> class DensityDescriptor,
          template<typename U2> class VolFracDescriptor
        >
void SmagorinskyBoussinesqThermalProcessor3D<T,DensityDescriptor,VolFracDescriptor>::process (
        Box3D domain,
        BlockLattice3D<T,DensityDescriptor>& fluid,
        BlockLattice3D<T,VolFracDescriptor>& temperature )
{
    typedef DensityDescriptor<T> D;
    enum {
        velOffset   = VolFracDescriptor<T>::ExternalField::velocityBeginsAt,
        forceOffset = DensityDescriptor<T>::ExternalField::forceBeginsAt
    };
    Dot3D offset = computeRelativeDisplacement(fluid, temperature);

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
                T *u = temperature.get(iX+offset.x,iY+offset.y,iZ+offset.z).getExternal(velOffset);
                Array<T,DensityDescriptor<T>::d> vel;
                fluid.get(iX,iY,iZ).computeVelocity(vel);
                vel.to_cArray(u);

                // Computation of the Boussinesq force
                T *force = fluid.get(iX,iY,iZ).getExternal(forceOffset);
                // Temperature is the order-0 moment of the advection-diffusion lattice.
                //   You can compute it with the method computeDensity().
                T localTemperature = temperature.get(iX+offset.x,iY+offset.y,iZ+offset.z).computeDensity();
                const T diffT = localTemperature - T0;
                for (pluint iD = 0; iD < D::d; ++iD)
                {
                    force[iD] = gravOverDeltaTemp[iD] * diffT;
                }
            }
        }
    }
}

template< typename T,
          template<typename U1> class DensityDescriptor,
          template<typename U2> class VolFracDescriptor
        >
SmagorinskyBoussinesqThermalProcessor3D<T,DensityDescriptor,VolFracDescriptor>*
    SmagorinskyBoussinesqThermalProcessor3D<T,DensityDescriptor,VolFracDescriptor>::clone() const
{
    return new SmagorinskyBoussinesqThermalProcessor3D<T,DensityDescriptor,VolFracDescriptor>(*this);
}

*/

}  // namespace plb

#endif  // BOUSSINESQ_THERMAL_PROCESSOR_3D_HH
