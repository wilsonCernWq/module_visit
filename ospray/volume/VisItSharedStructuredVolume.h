// ======================================================================== //
// Copyright 2009-2017 Intel Corporation                                    //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#pragma once

#include "volume/structured/StructuredVolume.h"

namespace ospray {
  namespace visit {

    //! \brief A concrete implementation of the StructuredVolume class
    //!  in which the voxel data is laid out in memory in XYZ order and
    //!  provided via a shared data buffer.
    //!
    struct OSPRAY_SDK_INTERFACE VisItSharedStructuredVolume 
      : public StructuredVolume
    {
    public:
      VisItSharedStructuredVolume();
      virtual ~VisItSharedStructuredVolume();

      //! A string description of this class.
      virtual std::string toString() const override;

      //! Allocate storage and populate the volume,
      //! called through the OSPRay API.
      virtual void commit() override;

      //! Copy voxels into the volume at the given index; not allowed on
      //!  VisItSharedStructuredVolume.
      virtual int setRegion(const void *source,
			    const vec3i &index,
			    const vec3i &count) override
      { return false; };

      //! Complete volume initialization (only on first commit).
      virtual void finish();

    protected:
      //! Create the equivalent ISPC volume container.
      void createEquivalentISPC() override;

      //! Called when a dependency of this object changes.
      void dependencyGotChanged(ManagedObject *object) override;
      
      //! The voxelData object upon commit().
      Data *voxelData {nullptr};

      //! The number of voxels
      size_t voxelCount {0};
      
      //! The volume dimensions
      vec3i dims {0,0,0};

      //! The global volume bounding box
      box3f globalBoundingBox {vec3f{0.f,0.f,0.f},vec3f{0.f,0.f,0.f}};

      //! The global volume clipping box
      box3f globalClippingBox {vec3f{0.f,0.f,0.f},vec3f{0.f,0.f,0.f}};
	    
      //! flags
      bool useGridAccelerator {false}; //! flag to use grid accelerator
    };

  };
} // ::ospray
