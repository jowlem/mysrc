#ifndef TAG_PROCESSOR_3D_HH
#define TAG_PROCESSOR_3D_HH

#include "TagProcessor3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "core/util.h"
#include "finiteDifference/finiteDifference3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "particles/particleProcessingFunctional3D.h"
#include "particles/particleField3D.h"
#include "dataProcessors/metaStuffFunctional3D.h"
#include "core/geometry3D.h"
#include "core/plbDebug.h"
#include "core/blockStatistics.h"
#include "atomicBlock/atomicBlock3D.h"
#include <algorithm>

namespace plb {


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Jonathan
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template< typename T, template<typename U> class ParticleDescriptor>
TagProcessor3D<T,ParticleDescriptor>::TagProcessor3D()
{   }



template< typename T, template<typename U> class ParticleDescriptor>
void TagProcessor3D<T,ParticleDescriptor>::processGenericBlocks (
        Box3D domain,
        std::vector<AtomicBlock3D*> blocks )
{


    PLB_PRECONDITION( blocks.size()==2 );
    ParticleField3D<T,ParticleDescriptor>& particleField
      = *dynamic_cast<ParticleField3D<T,ParticleDescriptor>*>(blocks[0]);
    BlockLattice3D<T,ParticleDescriptor>& fluid =
          *dynamic_cast<BlockLattice3D<T,ParticleDescriptor>*>(blocks[1]);
    std::vector<Particle3D<T,ParticleDescriptor>*> particles;

    Dot3D offset = computeRelativeDisplacement(fluid, particleField);
    #ifdef PLB_DEBUG
        Box3D bbox(fluid.getBoundingBox());
    #endif
    particleField.findParticles(domain, particles);
    T tag=1;
    Array<T,3> pos;
    Dot3D loc(fluid.getLocation());

    for (pluint iParticle=0; iParticle<particles.size(); ++iParticle) {

                pos=particles[iParticle]->getPosition() - Array<T,3>(loc.x,loc.y,loc.z);

                T localVolfrac = fluid.get( (plint) pos[0]+offset.x, (plint) pos[1]+offset.y, (plint) pos[2]+offset.z).computeDensity();

                if(localVolfrac > 0.2) tag = 1;
                else tag = 0;

                particles[iParticle]->setTag(tag);

                }

}


template< typename T, template<typename U> class ParticleDescriptor>

        void TagProcessor3D<T,ParticleDescriptor>::getTypeOfModification (
                std::vector<modif::ModifT>& modified ) const
        {
            modified[0] = modif::dynamicVariables;  // Particle field.
            modified[1] = modif::nothing;  // Fluid.
        }

template< typename T, template<typename U> class ParticleDescriptor>
TagProcessor3D<T,ParticleDescriptor>*
    TagProcessor3D<T,ParticleDescriptor>::clone() const
{
    return new TagProcessor3D<T,ParticleDescriptor>(*this);
}



}  // namespace plb

#endif  // BOUSSINESQ_THERMAL_PROCESSOR_3D_HH
