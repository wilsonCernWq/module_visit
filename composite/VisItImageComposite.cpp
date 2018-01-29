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

#include "VisItImageComposite.h"
#include "VisItImageComposite_ispc.h"

#include "common/OSPCommon.h"
#include "ospcommon/tasking/parallel_for.h"

#include <vector>
#include <limits>
#include <cmath>
#include <map>
#include <mutex>

using namespace ospray;

namespace visit {

  /* nothing to do, actually - this is only an example */
  extern "C" void Experiment()
  {
    std::cout << "[ospray] experiment_visit from CXX" << std::endl;
    ispc::ISPC_Experiment();
  }

  extern "C" void CompositeBackground(int *screen,
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
    const int tasks = std::ceil(static_cast<float>(psize) / 
				static_cast<float>(batch));
    tasking::parallel_for(tasks, [=](int taskIndex) {
	ispc::ISPC_CompositeBackground(taskIndex * batch,
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
    const int tasks = std::ceil(static_cast<float>(psize) / 
				static_cast<float>(batch));
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
    const int tasks = std::ceil(static_cast<float>(psize) / 
				static_cast<float>(batch));
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
  
