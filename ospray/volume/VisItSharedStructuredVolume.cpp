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
    void VisItSharedStructuredVolume::finish() {}
    //----------------------------------------------------------------------//
    void VisItSharedStructuredVolume::commit()
    {           
      // Create the equivalent ISPC volume container.
      if (ispcEquivalent == nullptr) createEquivalentISPC();
      // StructuredVolume commit actions.
      // In order to completely disable grid accelerator, we need an ugly hack
      bool trueFinished = StructuredVolume::finished;
      StructuredVolume::finished = true;
      StructuredVolume::commit();
      StructuredVolume::finished = trueFinished;
      // Complete volume initialization (only on first commit).
      if (!StructuredVolume::finished)
      {
        // Make the voxel value range visible to the application.
        if (findParam("voxelRange") == NULL)
        {
          set("voxelRange", StructuredVolume::voxelRange);
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
        StructuredVolume::finished = true;
      }
    }
    //----------------------------------------------------------------------//
    void VisItSharedStructuredVolume::createEquivalentISPC()
    {
      // Check if use grid accelerator before building the volume
      useGridAccelerator = (bool)getParam1i("useGridAccelerator", 0);
      if (::ospray::visit::CheckVerbose())
      {
        std::cout << "[ospray] using grid accelerator = "
                  << useGridAccelerator << std::endl;
      }
      // Get the voxel type.
      voxelType = getParamString("voxelType", "unspecified");
      const OSPDataType ospVoxelType = getVoxelType();
      exitOnCondition(ospVoxelType == OSP_UNKNOWN, "Unrecognized voxel type");
      // Get the volume dimensions.
      dims = getParam3i("dimensions", vec3i(0));
      exitOnCondition(reduce_min(dims) <= 0, "Invalid volume dimensions");
      // Get the voxel data.
      voxelData = (Data *) getParamObject("voxelData", nullptr);
      if (voxelData)
      {
        exitOnCondition(!(voxelData->flags & OSP_DATA_SHARED_BUFFER),
             "The voxel data buffer was not created with the "
             "OSP_DATA_SHARED_BUFFER flag");
      }
      // The voxel count.
      voxelCount = (size_t) dims.x * (size_t) dims.y * (size_t) dims.z;

      // Get global volume bounding box
      globalBoundingBox = box3f(getParam3f("volumeGlobalBoundingBoxLower",
					   vec3f(0.f)),
				getParam3f("volumeGlobalBoundingBoxUpper",
					   vec3f(0.f)));

      // Get global volume bounding box
      globalClippingBox = box3f(getParam3f("volumeGlobalClippingBoxLower",
					   vec3f(0.f)),
				getParam3f("volumeGlobalClippingBoxUpper",
					   vec3f(0.f)));
      std::cout << globalBoundingBox.lower.x << " "
		<< globalBoundingBox.lower.y << " "
		<< globalBoundingBox.lower.z << " | "
		<< globalBoundingBox.upper.x << " "
		<< globalBoundingBox.upper.y << " "
		<< globalBoundingBox.upper.z << std::endl;
      
      // Create an ISPC VisItSharedStructuredVolume object and assign
      // type-specific function pointers.
      ispcEquivalent
          = ispc::VisItSSV_createInstance(this,
					  (int)useGridAccelerator,
					  ospVoxelType,
					  (const ispc::box3f &)
					  globalClippingBox,
					  (const ispc::box3f &)
					  globalBoundingBox,
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
