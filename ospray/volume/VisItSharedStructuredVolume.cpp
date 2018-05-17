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

//ospray common
#include "common/Data.h"
//ospray
#include "VisItModuleCommon.h"
//header
#include "VisItSharedStructuredVolume.h"
//ISPC
#include "VisItSharedStructuredVolume_ispc.h"
#include "StructuredVolume_ispc.h"

namespace ospray {
  namespace visit {
    //----------------------------------------------------------------------//
    VisItSharedStructuredVolume::VisItSharedStructuredVolume()
    {
      if (::ospray::visit::CheckVerbose())
      {
        std::cout << "[ospray] New VisItSharedStructuredVolume" << std::endl;
      }
    }
    //----------------------------------------------------------------------//
    VisItSharedStructuredVolume::~VisItSharedStructuredVolume()
    {
      if (voxelData) voxelData->unregisterListener(this);
    }
    //----------------------------------------------------------------------//
    std::string VisItSharedStructuredVolume::toString() const
    {
      return ("ospray::VisItSharedStructuredVolume<" + voxelType + ">");
    }
    //----------------------------------------------------------------------//
    void VisItSharedStructuredVolume::finish() 
    {
      // Make the voxel value range visible to the application.
      if (!hasParam("voxelRange"))
      {
	setParam("voxelRange", StructuredVolume::voxelRange);
      }
      else {
	StructuredVolume::voxelRange =
	  getParam2f("voxelRange", StructuredVolume::voxelRange);
      }
      // Build empty space skipping grid
      if (this->useGridAccelerator)
      {
	buildAccelerator();
      }
      // Volume finish actions.
      Volume::finish();
    }
    //----------------------------------------------------------------------//
    void VisItSharedStructuredVolume::commit()
    {
      // Create the equivalent ISPC volume container.
      if (ispcEquivalent == nullptr) createEquivalentISPC();
      // StructuredVolume commit actions.
      // In order to completely disable grid accelerator, we need an ugly hack
      bool initialized = finished;
      finished = true;
      StructuredVolume::commit();
      finished = initialized;
      // Complete volume initialization (only on first commit).
      if (!finished)
      {
	finish();
	finished = true;
      }
      
      // Check if use grid accelerator before building the volume
      useGridAccelerator = (bool)getParam1i("useGridAccelerator", 0);
      if (::ospray::visit::CheckVerbose())
      {
        std::cout << "[ospray] using grid accelerator = "
                  << useGridAccelerator << std::endl;
      }
      ispc::VisItSSV_updateGridAccelerator(ispcEquivalent, 
					   (int)useGridAccelerator);

      // Get global volume bounding box
      globalBoundingBox = box3f(getParam3f("volumeGlobalBoundingBoxLower",
					   vec3f(0.f)),
				getParam3f("volumeGlobalBoundingBoxUpper",
					   vec3f(0.f)));
      ispc::VisItSSV_updateBoundingBox(ispcEquivalent, (const ispc::box3f &)globalBoundingBox);

      // Get global volume bounding box
      globalClippingBox = box3f(getParam3f("volumeGlobalClippingBoxLower",
					   vec3f(0.f)),
				getParam3f("volumeGlobalClippingBoxUpper",
					   vec3f(0.f)));
      ispc::VisItSSV_updateClippingBox(ispcEquivalent, (const ispc::box3f &)globalClippingBox);
    }
    //----------------------------------------------------------------------//
    void VisItSharedStructuredVolume::createEquivalentISPC()
    {
      // Get the voxel type.
      voxelType = getParamString("voxelType", "unspecified");
      const OSPDataType ospVoxelType = getVoxelType();
      if (ospVoxelType == OSP_UNKNOWN)
	throw std::runtime_error("unrecognized voxel type");
      // Get the volume dimensions.
      dims = getParam3i("dimensions", vec3i(0));
      if (reduce_min(dims) <= 0)
	throw std::runtime_error("invalid volume dimensions");
      // Get the voxel data.
      voxelData = (Data *) getParamObject("voxelData", nullptr);
      if (voxelData && !(voxelData->flags & OSP_DATA_SHARED_BUFFER)) {
	postStatusMsg(1)
	  << "WARNING: The voxel data buffer was not created with the "
	  << "OSP_DATA_SHARED_BUFFER flag";
      }
      // The voxel count.
      voxelCount = (size_t) dims.x * (size_t) dims.y * (size_t) dims.z;
      
      // Create an ISPC VisItSharedStructuredVolume object and assign
      // type-specific function pointers.
      ispcEquivalent
          = ispc::VisItSSV_createInstance(this,
					  ospVoxelType,
					  (const ispc::vec3i &) dims,
					  voxelData->data);
      // Listen for changes to voxelData.
      if (voxelData) { voxelData->registerListener(this); }
    }
    //----------------------------------------------------------------------//
    void VisItSharedStructuredVolume::dependencyGotChanged
        (ManagedObject *object) {
      // Rebuild volume accelerator when voxelData is committed.
      if (object == voxelData && ispcEquivalent && useGridAccelerator) {
        StructuredVolume::buildAccelerator();
      }
    }
    //----------------------------------------------------------------------//
    // A volume type with XYZ storage order. The voxel data is provided by the
    // application via a shared data buffer.
    OSP_REGISTER_VOLUME(VisItSharedStructuredVolume,
                        visit_shared_structured_volume);

  };
} // ::ospray
