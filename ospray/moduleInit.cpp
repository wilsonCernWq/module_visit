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

/*! \file ospray/moduleInit \brief Defines the module initialization callback */

#include "VisItModuleCommon.h"
#include "moduleInit_ispc.h"
#include "ospcommon/utility/getEnvVar.h"
#include "volume/VisItSharedStructuredVolume.h"
#include <cmath>

/*! _everything_ in the ospray core universe should _always_ be in the
  'ospray' namespace. */
namespace ospray {
  
  /*! though not required, it is good practice to put any module into
    its own namespace (isnide of ospray:: ). Unlike for the naming of
    library and init function, the naming for this namespace doesn't
    particularlly matter. E.g., 'bilinearPatch', 'module_blp',
    'bilinar_patch' etc would all work equally well. */
  namespace visit {
    
    /*! the actual module initialization function. This function gets
      called exactly once, when the module gets first loaded through
      'ospLoadModule'. Notes:

      a) this function does _not_ get called if the application directly
      links to libospray_module_<modulename> (which it
      shouldn't!). Modules should _always_ be loaded through
      ospLoadModule. 

      b) it is _not_ valid for the module to do ospray _api_ calls
      inside such an intiailzatoin function. Ie, you can _not_ do a
      ospLoadModule("anotherModule") from within this function (but
      you could, of course, have this module dynamically link to the
      other one, and call its init function)

      c) in order for ospray to properly resolve that symbol, it
      _has_ to have extern C linkage, and it _has_ to correspond to
      name of the module and shared library containing this module
      (see comments regarding library name in CMakeLists.txt)
    */

    bool verbose = false;
    extern "C" void ospray_init_module_visit()
    {
      /* initialize global variables */
      using ospcommon::utility::getEnvVar;
      auto OSPRAY_VERBOSE = 
	getEnvVar<int>("OSPRAY_VERBOSE").value_or(0);
      ::ospray::visit::verbose = OSPRAY_VERBOSE > 0; 
      if (::ospray::visit::CheckVerbose()) 
      {
	std::cout << "[ospray] initializing the 'visit' module" << std::endl;
      }
    }

    /* nothing to do, actually - this is only an example */
    extern "C" void Experiment()
    {
      std::cout << "[ospray] experiment_visit from CXX" << std::endl;
      ispc::ISPC_Experiment();
    }

    extern "C" void ComposeBackground(int *screen,
				      int *compositedImageExtents,
				      int  compositedImageWidth,
				      int  compositedImageHeight,
				      float *compositedImageBuffer,
				      unsigned char *opaqueImageColor,
				      float         *opaqueImageDepth,
				      unsigned char *&imgFinal)
    {
      const int batch = TILE_SIZE * TILE_SIZE;
      const int psize = screen[0] * screen[1];
      const int tasks = std::ceil(psize / batch);
      tasking::parallel_for(tasks, [=](int taskIndex) {
	  ispc::ISPC_ComposeBackground(taskIndex * batch,
				       std::min(taskIndex * batch + batch,
						psize),
				       compositedImageExtents,
				       compositedImageWidth,
				       compositedImageHeight,
				       compositedImageBuffer,
				       screen[0], screen[1],
				       (ospray::uint8*)opaqueImageColor,
				       opaqueImageDepth,
				       (ospray::uint8*)imgFinal);

	});
    }
    
    extern "C" void BlendFrontToBack(const int   *blendExtents,
				     const int   *srcExtents,
				     const float *srcImage,
				     const int   *dstExtents,
				     float      *&dstImage)
    {  
      // image sizes
      const int srcX = srcExtents[1] - srcExtents[0];
      const int srcY = srcExtents[3] - srcExtents[2];
      const int dstX = dstExtents[1] - dstExtents[0];
      const int dstY = dstExtents[3] - dstExtents[2];
      // determin the region to blend
      const int startX = 
	std::max(std::max(blendExtents[0], srcExtents[0]), dstExtents[0]);
      const int startY = 
	std::max(std::max(blendExtents[2], srcExtents[2]), dstExtents[2]);
      const int endX = 
	std::min(std::min(blendExtents[1], srcExtents[1]), dstExtents[1]);
      const int endY = 
	std::min(std::min(blendExtents[3], srcExtents[3]), dstExtents[3]);
      // render
      const int W = endX - startX;
      const int H = endY - startY;
      const int psize = W * H;
      const int batch = TILE_SIZE * TILE_SIZE;
      const int tasks = std::ceil(psize / batch);      
      tasking::parallel_for(tasks, [=](int taskIndex) {
	  ispc::ISPC_BlendFrontToBack(taskIndex * batch,
				      std::min(taskIndex * batch + batch,
					       psize),
				      startX, startY, W, H,
				      srcX, srcY, dstX, dstY,
				      blendExtents,
				      srcExtents, srcImage,
				      dstExtents, dstImage);
	});
    }

    extern "C" void BlendBackToFront(const int   *blendExtents,
				     const int   *srcExtents,
				     const float *srcImage,
				     const int   *dstExtents,
				     float      *&dstImage)
    {  
      // image sizes
      const int srcX = srcExtents[1] - srcExtents[0];
      const int srcY = srcExtents[3] - srcExtents[2];
      const int dstX = dstExtents[1] - dstExtents[0];
      const int dstY = dstExtents[3] - dstExtents[2];
      // determin the region to blend
      const int startX = 
	std::max(std::max(blendExtents[0], srcExtents[0]), dstExtents[0]);
      const int startY = 
	std::max(std::max(blendExtents[2], srcExtents[2]), dstExtents[2]);
      const int endX = 
	std::min(std::min(blendExtents[1], srcExtents[1]), dstExtents[1]);
      const int endY = 
	std::min(std::min(blendExtents[3], srcExtents[3]), dstExtents[3]);
      // render
      const int W = endX - startX;
      const int H = endY - startY;
      const int psize = W * H;
      const int batch = TILE_SIZE * TILE_SIZE;
      const int tasks = std::ceil(psize / batch);
      tasking::parallel_for(tasks, [=](int taskIndex) {
	  ispc::ISPC_BlendBackToFront(taskIndex * batch,
				      std::min(taskIndex * batch + batch,
					       psize),
				      startX, startY, W, H,
				      srcX, srcY, dstX, dstY,
				      blendExtents,
				      srcExtents, srcImage,
				      dstExtents, dstImage);
	});
    }

  }; // ::ospray::visit

}; // ::ospray
  
