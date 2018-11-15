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

#ifndef DVF_PROCESSOR_3D_HH
#define DVF_PROCESSOR_3D_HH

#include "VolFracDensityCoupling.h"
#include "atomicBlock/blockLattice3D.h"
#include "core/util.h"
#include "finiteDifference/finiteDifference3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

namespace plb {


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Jonathan
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template< typename T,
          template<typename U1> class DensityDescriptor,
          template<typename U2> class VolFracDescriptor
        >
DVFCouplingProcessor3D<T,DensityDescriptor,VolFracDescriptor>::
        DVFCouplingProcessor3D(T rhoP_, T convers_)
    :  rhoP(rhoP_), convers(convers_)
{


}


template< typename T,
          template<typename U1> class DensityDescriptor,
          template<typename U2> class VolFracDescriptor
        >
void DVFCouplingProcessor3D<T,DensityDescriptor,VolFracDescriptor>::process (
        Box3D domain,
        BlockLattice3D<T,DensityDescriptor>& density,
        BlockLattice3D<T,VolFracDescriptor>& volumefrac )
{
    //typedef DensityDescriptor<T> D;
    enum {
        velOffset   = VolFracDescriptor<T>::ExternalField::velocityBeginsAt,
    };
    Dot3D offset = computeRelativeDisplacement(density, volumefrac);
    Array<T,3> v_sed(0., 0., 0.);

    for (plint iX=domain.x0; iX<=domain.x1; ++iX)
    {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY)
        {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ)
            {
                // Velocity coupling
                T *u = volumefrac.get(iX+offset.x,iY+offset.y,iZ+offset.z).getExternal(velOffset);


                T dens = density.get(iX,iY,iZ).computeDensity();

                const T g=9.81;
                const T deq=63e-6;
                const T Cd=124.73;

                T prevVolfrac=volumefrac.get(iX+offset.x,iY+offset.y,iZ+offset.z).computeDensity();
                if (prevVolfrac<0){
                    volumefrac.get(iX+offset.x,iY+offset.y,iZ+offset.z).defineDensity(0);
                }

                T localVolfrac = volumefrac.get(iX+offset.x,iY+offset.y,iZ+offset.z).computeDensity();
                if(localVolfrac > 0)
			          v_sed[2] = convers*sqrt((4*g*deq*(rhoP-dens))/(3*dens*Cd));
                else v_sed[2] = 0;

                v_sed.to_cArray(u);
            }
        }
    }
}

template< typename T,
          template<typename U1> class DensityDescriptor,
          template<typename U2> class VolFracDescriptor
        >
DVFCouplingProcessor3D<T,DensityDescriptor,VolFracDescriptor>*
    DVFCouplingProcessor3D<T,DensityDescriptor,VolFracDescriptor>::clone() const
{
    return new DVFCouplingProcessor3D<T,DensityDescriptor,VolFracDescriptor>(*this);
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
