
#ifndef TAG_PROCESSOR_3D_H
#define TAG_PROCESSOR_3D_H

#include "core/globalDefs.h"
#include "core/block3D.h"
#include "atomicBlock/dataProcessor3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "core/functions.h"
#include "atomicBlock/dataProcessingFunctional3D.h"
#include "atomicBlock/reductiveDataProcessingFunctional3D.h"
#include "atomicBlock/atomicContainerBlock3D.h"

namespace plb {


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Jonathan
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template< typename T, template<typename U> class ParticleDescriptor>
class TagProcessor3D : public BoxProcessingFunctional3D
{
public:

    TagProcessor3D();

    virtual void processGenericBlocks( Box3D domain,
                          std::vector<AtomicBlock3D*> fields );
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual TagProcessor3D<T,ParticleDescriptor>* clone() const;

};

}

#endif
