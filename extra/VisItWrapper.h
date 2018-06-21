// ======================================================================== //
// Copyright Qi WU                                                          //
//   Scientific Computing and Image Institution                             //
//   University of Utah                                                     //
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
// this file will be installed so can expose new functions to the users

#pragma once

#include <ospray/ospray.h>

namespace ospray {
  namespace visit {

    /** 
     * Helper Functions
     */
    template<typename T> void ospray_rm(T& obj) {
      if (!obj) { ospRelease(obj); obj = NULL; }
    }

    template<typename T> struct Object {
    protected:
      bool init;
      T    self;
    public:
      Object() { init = false; self = NULL; }
      virtual ~Object() { Delete(); }
      T operator*() { return self; }
      void Delete() { ospray_rm(self); init = false; }
    };

    /**
     * Transfer Function Wrapper
     */
    struct TransferFunction : public Object<OSPTransferFunction> {
    public:
      TransferFunction() : Object<OSPTransferFunction>() {}
      virtual ~TransferFunction() {};
      void Set(const void *table, const unsigned int size,
               const double datamin, const double datamax);      
    };





/*  
    class Volume
    {
    public:
      Volume(const int i, OSPVisItPatcheses* p) {
        initialized = false;
        parent      = p;
        index       = i;
        pixelnum    = 0;
        model  = NULL;
        volume = NULL;
        voxels = NULL;
        framebuffer = NULL;
        background  = NULL;
      }
      ~Volume() { Delete(); }
      void Delete() {
        OSPRAY_SAFE_RM(model);
        OSPRAY_SAFE_RM(volume);
        OSPRAY_SAFE_RM(voxels);
        OSPRAY_SAFE_RM(framebuffer);
        OSPRAY_SAFE_RM(background);
      }
      // create volume and model
      void Set(const int   voxelEnumVTK,
               const void *voxelDataPtr,
               const double *X, const double *Y, const double *Z,
               const int    nX, const int    nY, const int    nZ,
               const double pbox[6],
               const double bbox[6]);

      // create framebuffer
      void Render(OSPRenderer,
                  const unsigned int W,
                  const unsigned int H,
                  const int Xs, const int Ys,
                  const float* maxDepthBuffer,
                  const int    maxDepthStride,
                  float*);

    private:
      bool initialized;
      OSPVisItPatcheses *parent;
      int                index;    // patch index
      unsigned int       pixelnum; // number of pixels
      OSPModel           model;
      OSPVolume          volume;
      OSPData            voxels;
      OSPFrameBuffer     framebuffer;
      OSPTexture2D       background;
    };
*/
  
/*     OSPVisItPatches() { */
/*       enableShading = false; */
/*       enableDVR     = false; */
/*       enableGridAccelerator = false; */
/*       enablePreIntegration = true; */
/*       Ks    = 1.; */
/*       Ns    = 15.; */
/*       samplingRate  = 3.; */
/*       scaling[0] = scaling[1] = scaling[2] = 1.f; */
/*       scaledGlobalBBoxLower[0] = 0.; */
/*       scaledGlobalBBoxLower[1] = 0.; */
/*       scaledGlobalBBoxLower[2] = 0.; */
/*       scaledGlobalBBoxupper[0] = 0.; */
/*       scaledGlobalBBoxupper[1] = 0.; */
/*       scaledGlobalBBoxupper[2] = 0.; */
/*     }; */
/*     ~OSPVisItPatches() { Delete(); }; */
/*     void Delete() { */
/*       /\* for (int i = 0; i < patches.size(); ++i) { *\/ */
/*       /\*     patches[i]->Delete(); *\/ */
/*       /\* } *\/ */
/*     }; */

/*     OSPVisItPatches::Volume* GetPatch(int id) { return &patches[id]; }; */
/*     void RegisterPatch(int id) { */
/*       if (patches.find(id) == patches.end()) */
/*         patches[id] = Volume(id, this));     */
/*   };    */

/*   // parameters */
/*   void SetFlagDVR(const bool f) { enableDVR = f; /\* not used yet *\/}; */
/*   void SetFlagShading(const bool f) { enableShading = f; }; */
/*   void SetFlagGridAccelerator(const bool f) { enableGridAccelerator = f; }; */
/*   void SetFlagPreIntegration(const bool f) { enablePreIntegration = f; }; */
/*   void SetSamplingRate(const double v) { samplingRate = v; }; */
/*   void SetKs(const double v) { Ks = v; }; */
/*   void SetNs(const double v) { Ns = v; }; */
/*   void SetDataBBox(const double v[6]) { */
/*     dataBBox[0] = v[0]; /\* x_min *\/ */
/*     dataBBox[1] = v[1]; /\* x_max *\/ */
/*     dataBBox[2] = v[2]; /\* y_min *\/ */
/*     dataBBox[3] = v[3]; /\* y_max *\/ */
/*     dataBBox[4] = v[4]; /\* z_min *\/ */
/*     dataBBox[5] = v[5]; /\* z_max *\/ */
/*   }    */
/*   void SetScaling(const double v[3]) { /\* called after SetDataBBox *\/ */
/*     scaling[0] = s[0]; */
/*     scaling[1] = s[1]; */
/*     scaling[2] = s[2];	     */
/*     scaledGlobalBBoxLower[0] = dataBBox[0] * scaling[0]; */
/*     scaledGlobalBBoxLower[1] = dataBBox[2] * scaling[1]; */
/*     scaledGlobalBBoxLower[2] = dataBBox[4] * scaling[2]; */
/*     scaledGlobalBBoxUpper[0] = dataBBox[1] * scaling[0]; */
/*     scaledGlobalBBoxUpper[1] = dataBBox[3] * scaling[1]; */
/*     scaledGlobalBBoxUpper[2] = dataBBox[5] * scaling[2]; */
/*   }; */
    
/* private: */
    
/*   bool enableDVR;     // Distributed Volume Renderer */
/*   bool enableShading; */
/*   bool enableGridAccelerator; */
/*   bool enablePreIntegration; */
/*   double samplingRate; */
/*   double Ks; */
/*   double Ns; */
/*   double scaling[3]; */
/*   double dataBBox[6]; */
/*   double scaledGlobalBBoxLower[3]; */
/*   double scaledGlobalBBoxUpper[3]; */

/* public: */

/*   OSPVisItPatch::TransferFunction tfn;     */
/*   std::map<int, OSPVisItPatch> patches; */

/* }; */

};
};