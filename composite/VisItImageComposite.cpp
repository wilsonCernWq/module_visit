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

// #include <IceT.h>
// #include <IceTMPI.h>

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
    const int tasks = std::ceil(static_cast<float>(psize) / 
				static_cast<float>(batch));
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

  // // ====================================================================== //
  // //                                                                        
  // // ====================================================================== //
  // bool icetInitialized = false;
  // int icetScreen[2];
  // IceTContext icetContext;
  // IceTCommunicator icetComm;        
  // IceTInt icetMPIRank, icetMPISize;
  // IceTDouble  icetMatProj[16] = {
  //   1.0, 0.0, 0.0, 0.0,
  //   0.0, 1.0, 0.0, 0.0,
  //   0.0, 0.0, 1.0, 0.0,
  //   0.0, 0.0, 0.0, 1.0
  // };
  // IceTDouble  icetMatMV[16] = {
  //   1.0, 0.0, 0.0, 0.0,
  //   0.0, 1.0, 0.0, 0.0,
  //   0.0, 0.0, 1.0, 0.0,
  //   0.0, 0.0, 0.0, 1.0
  // };
  // IceTFloat   icetBgColor[4] = {0.0, 0.0, 0.0, 0.0};

  // struct Image {
  //   const float* data;
  //   int extents[4];
  //   float depth;
  //   void SetTile(const float* d, const int e[4], const float& z) {
  //     data = d;
  //     extents[0] = e[0];
  //     extents[1] = e[1];
  //     extents[2] = e[2];
  //     extents[3] = e[3];
  //     depth = z;
  //   }
  // } icetImageLocal;

  // extern "C" void CompositeInit(int W, int H)
  // {
  //   if (!icetInitialized) 
  //   {
  //     icetScreen[0] = W;
  //     icetScreen[1] = H;
  //     // Create
  //     icetComm = icetCreateMPICommunicator(MPI_COMM_WORLD);
  //     icetContext = icetCreateContext(icetComm);
  //     icetGetIntegerv(ICET_RANK, &icetMPIRank);
  //     icetGetIntegerv(ICET_NUM_PROCESSES, &icetMPISize);
  //     icetDiagnostics(ICET_DIAG_FULL);

  //     // Setup IceT for alpha-blending compositing
  //     icetCompositeMode(ICET_COMPOSITE_MODE_BLEND);
  //     icetSetColorFormat(ICET_IMAGE_COLOR_RGBA_FLOAT);
  //     icetSetDepthFormat(ICET_IMAGE_DEPTH_NONE);
  //     icetDiagnostics(ICET_DIAG_ERRORS | ICET_DIAG_WARNINGS);
  //     icetEnable(ICET_ORDERED_COMPOSITE);
  //     icetDiagnostics(ICET_DIAG_ERRORS | ICET_DIAG_WARNINGS);
  //     icetDisable(ICET_INTERLACE_IMAGES);
  //     // Safety
  //     MPI_Barrier(MPI_COMM_WORLD);
  //     //
  //     icetResetTiles();
  //     icetAddTile(0, 0, W, H, 0);
  //     icetPhysicalRenderSize(W, H);
  //     //
  //     icetStrategy(ICET_STRATEGY_SEQUENTIAL);
  //     //icetSingleImageStrategy(ICET_SINGLE_IMAGE_STRATEGY_TREE);
  //     //icetSingleImageStrategy(ICET_SINGLE_IMAGE_STRATEGY_RADIXK);
  //     icetSingleImageStrategy(ICET_SINGLE_IMAGE_STRATEGY_BSWAP);
  //     //
  //     icetInitialized = true;
  //   }
  // }
  
  // extern "C" void CompositeSetTile(const float* data, const int*e, 
  // 				   const float& z, float*& output)
  // {
  //   // Setup tiles
  //   icetImageLocal.SetTile(data, e, z);

  //   // Gather depths
  //   std::vector<float> all_depths(icetMPISize);
  //   std::vector<IceTInt> all_orders(icetMPISize);
  //   MPI_Allgather(&(icetImageLocal.depth), 1, MPI_FLOAT, 
  // 		  all_depths.data(), 1, MPI_FLOAT, MPI_COMM_WORLD);

  //   // Sort the rank in compositingOrder                      
  //   std::multimap<float,int> ordered_depths;
  //   for (int i = 0; i < icetMPISize; i++)
  //   {
  //     ordered_depths.insert(std::pair<float, int>(all_depths[i], i)); 
  //   }

  //   // Setup IceT
  //   int i = 0;
  //   for (auto it = ordered_depths.begin(); it != ordered_depths.end(); ++it){
  //     all_orders[i] = (*it).second;
  //     i++;
  //   }
  //   icetCompositeOrder(all_orders.data());  // front to back
  //   icetBoundingBoxf((static_cast<float>(e[0])/static_cast<float>(icetScreen[0]) - 0.5f) * 2,
  // 		     (static_cast<float>(e[1])/static_cast<float>(icetScreen[0]) - 0.5f) * 2,
  // 		     (static_cast<float>(e[2])/static_cast<float>(icetScreen[1]) - 0.5f) * 2,
  // 		     (static_cast<float>(e[3])/static_cast<float>(icetScreen[1]) - 0.5f) * 2,
  // 		     0.0, 0.0);
    
  //   // Compose
  //   IceTImage result = icetCompositeImage
  //     (icetImageLocal.data, NULL, e, icetMatProj, icetMatMV, icetBgColor);
  //   MPI_Barrier(MPI_COMM_WORLD);
  //   output = new float[icetScreen[0] * icetScreen[1] * 4];
  //   icetImageCopyColorf(result, output, ICET_IMAGE_COLOR_RGBA_FLOAT);
  // }

}; // ::ospray::visit
  
