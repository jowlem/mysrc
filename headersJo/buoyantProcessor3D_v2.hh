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

#ifndef BUOYANT_THERMAL_PROCESSOR_3Dv2_HH
#define BUOYANT_THERMAL_PROCESSOR_3Dv2_HH

#include "buoyantProcessor3D_v2.h"
#include "atomicBlock/blockLattice3D.h"
#include "core/util.h"
#include "finiteDifference/finiteDifference3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

namespace plb {


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Jonathan
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template< typename T,
          template<typename U1> class FluidDescriptor,
          template<typename U2> class TemperatureDescriptor
        >
BuoyanTermProcessor3Dv2<T,FluidDescriptor,TemperatureDescriptor>::
        BuoyanTermProcessor3Dv2(T gravity_, T rho0_, T rhoP_, Array<T,FluidDescriptor<T>::d> dir_)
    :  gravity(gravity_), rho0(rho0_), rhoP(rhoP_),
       dir(dir_)
{
    // We normalize the direction of the force vector.
    T normDir = std::sqrt(VectorTemplate<T,FluidDescriptor>::normSqr(dir));
    for (pluint iD = 0; iD < FluidDescriptor<T>::d; ++iD) {
        dir[iD] /= normDir;
    }
}


template< typename T,
          template<typename U1> class FluidDescriptor,
          template<typename U2> class TemperatureDescriptor
        >
void BuoyanTermProcessor3Dv2<T,FluidDescriptor,TemperatureDescriptor>::process (
        Box3D domain,
        BlockLattice3D<T,FluidDescriptor>& fluid,
        BlockLattice3D<T,TemperatureDescriptor>& temperature )
{
    typedef FluidDescriptor<T> D;
    enum {
        velOffset   = TemperatureDescriptor<T>::ExternalField::velocityBeginsAt,
        forceOffset = FluidDescriptor<T>::ExternalField::forceBeginsAt
    };
    Dot3D offset = computeRelativeDisplacement(fluid, temperature);

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
                // Velocity coupling
                T *u = temperature.get(iX+offset.x,iY+offset.y,iZ+offset.z).getExternal(velOffset);
                Array<T,TemperatureDescriptor<T>::d> v_sedi;
                temperature.get(iX+offset.x,iY+offset.y,iZ+offset.z).computeVelocity(v_sedi);
                Array<T,FluidDescriptor<T>::d> vel;
                fluid.get(iX,iY,iZ).computeVelocity(vel);
                T prevVolfrac=temperature.get(iX+offset.x,iY+offset.y,iZ+offset.z).computeDensity();
                if (prevVolfrac<0){
                    temperature.get(iX+offset.x,iY+offset.y,iZ+offset.z).defineDensity(0);
                }
		            T localVolfrac = temperature.get(iX+offset.x,iY+offset.y,iZ+offset.z).computeDensity();
		            if(localVolfrac > 0) vel[2] -= v_sedi[2];
                else vel[2] -= 0;
                vel.to_cArray(u);

                // Computation of the Boussinesq force
                T *force = fluid.get(iX,iY,iZ).getExternal(forceOffset);
                // Temperature is the order-0 moment of the advection-diffusion lattice.
                //   You can compute it with the method computeDensity().
                //T localVolfrac = temperature.get(iX+offset.x,iY+offset.y,iZ+offset.z).computeDensity();
                const T diffT = rhoP-rho0;
                for (pluint iD = 0; iD < D::d; ++iD)
                {
			if(localVolfrac > 0) force[iD] = -localVolfrac * gravOverrho0[iD] * diffT * 0.002;
      else force[iD] = 0.;
		//	force[iD] = 0.;

                }
            }
        }
    }
}

template< typename T,
          template<typename U1> class FluidDescriptor,
          template<typename U2> class TemperatureDescriptor
        >
BuoyanTermProcessor3Dv2<T,FluidDescriptor,TemperatureDescriptor>*
    BuoyanTermProcessor3Dv2<T,FluidDescriptor,TemperatureDescriptor>::clone() const
{
    return new BuoyanTermProcessor3Dv2<T,FluidDescriptor,TemperatureDescriptor>(*this);
}


/*
template< typename T,
          template<typename U1> class FluidDescriptor,
          template<typename U2> class TemperatureDescriptor
        >
SmagorinskyBoussinesqThermalProcessor3D<T,FluidDescriptor,TemperatureDescriptor>::
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
          template<typename U2> class TemperatureDescriptor
        >
void SmagorinskyBoussinesqThermalProcessor3D<T,FluidDescriptor,TemperatureDescriptor>::process (
        Box3D domain,
        BlockLattice3D<T,FluidDescriptor>& fluid,
        BlockLattice3D<T,TemperatureDescriptor>& temperature )
{
    typedef FluidDescriptor<T> D;
    enum {
        velOffset   = TemperatureDescriptor<T>::ExternalField::velocityBeginsAt,
        forceOffset = FluidDescriptor<T>::ExternalField::forceBeginsAt
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
                Array<T,FluidDescriptor<T>::d> vel;
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
          template<typename U1> class FluidDescriptor,
          template<typename U2> class TemperatureDescriptor
        >
SmagorinskyBoussinesqThermalProcessor3D<T,FluidDescriptor,TemperatureDescriptor>*
    SmagorinskyBoussinesqThermalProcessor3D<T,FluidDescriptor,TemperatureDescriptor>::clone() const
{
    return new SmagorinskyBoussinesqThermalProcessor3D<T,FluidDescriptor,TemperatureDescriptor>(*this);
}

*/

}  // namespace plb

#endif  // BOUSSINESQ_THERMAL_PROCESSOR_3D_HH
