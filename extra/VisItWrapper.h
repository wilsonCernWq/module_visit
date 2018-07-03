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

#include "VisItWrapperCore.h"
#include <string>

namespace ospray {
namespace visit {

  template<typename _CoreType, typename _OSPType> struct Manipulator {
  protected:
    typedef _CoreType CoreType;
    typedef _OSPType  OSPType;
    _CoreType *core;
  public:
    Manipulator(_CoreType& other) : core{&other} {}
    _OSPType   operator* () { return *(*core); }
    _CoreType* operator->() { return &(*core); }
  };
  
  /**
   * Transfer Function Wrapper
   */
  struct TransferFunction
    : public Manipulator<TransferFunctionCore, OSPTransferFunction>
  {
  public:
    TransferFunction(CoreType& other);
    void Set(const void *table, const unsigned int size,
             const double datamin, const double datamax);      
  };

  /**
   * Camera Wrapper
   */
  struct Camera
    : public Manipulator<CameraCore, OSPCamera>
  {
  public:
    Camera(CameraCore& other);
    double GetWindowExts(const int i) const { return core->windowExts[i]; }
    void Set(const bool ortho,
             const double camera_p[3], 
             const double camera_f[3], 
             const double camera_u[3], 
             const double fovy, 
             const double pan_ratio[2],
             const double zoom_ratio,
             const double near_clip,
             const double canvas_size[2],
             const int screen_size[2],
             const int tile_extents[4]);
    void SetScreen(const double xMin, const double xMax,
                   const double yMin, const double yMax);
  };

  /**
   * Light Wrapper
   */
  struct Light
    : public Manipulator<LightCore, OSPLight>
  {
  public:
    Light(LightCore& other);
    void Set(const bool ambient, const double i, 
             const double c, const double* d = NULL);
    void Set(const bool ambient, const double i, 
             const double cr, const double cg, const double cb,
             const double* d = NULL);
    void Set(const bool ambient, const double i, 
             const double c[3], const double* d = NULL);
  };

  /**
   * Renderer Wrapper
   */
  struct Renderer
    : public Manipulator<RendererCore, OSPRenderer>
  {
  public:
    Renderer(RendererCore& other);
    void Init();
    void ResetLights();
    Light AddLight();
    void FinalizeLights();
    void SetBackgroundBuffer(const unsigned char * color,
                             const float* depth, 
                             const int size[2]);
    void Set(const int aoSamples, const int spp, 
             const bool oneSidedLighting,
             const bool shadowsEnabled,
             const bool aoTransparencyEnabled);
    void Set(OSPCamera osp_camera);
    void Set(OSPModel osp_world);
  };

  /**
   * Model Wrapper
   */
  struct Model
    : public Manipulator<ModelCore, OSPModel>
  {
  public:
    Model(ModelCore& other);
    void Reset();
    void Init();
    void Set(OSPVolume);
  };

  /**
   * Volume Wrapper
   */
  struct Volume 
    : public Manipulator<VolumeCore, OSPVolume>
  {
  public:
    Volume(VolumeCore& other);
    bool Init(const std::string volume_type, 
              const OSPDataType data_type, const std::string data_char,
              const size_t data_size, const void* data_ptr);
    void Set(const bool useGridAccelerator, const bool adaptiveSampling,
             const bool preIntegration, const bool singleShade, 
             const bool gradientShadingEnabled, const double samplingRate, 
             const double Ks, const double Ns,
             const double *X, const double *Y, const double *Z, 
             const int nX, const int nY, const int nZ,
             const double dbox[6], const double cbox[6], 
             const osp::vec3f& global_upper,
             const osp::vec3f& global_lower,
             const osp::vec3f& scale,
             OSPTransferFunction tfn);
    static void ComputeGhostBounds(bool bound[6], 
                                   const unsigned char *ghosts, 
                                   const int gnX, 
                                   const int gnY, 
                                   const int gnZ);
  };

  /**
   * FrameBuffer Wrapper
   */
  struct FrameBuffer
    : public Manipulator<FrameBufferCore, OSPFrameBuffer>
  {
  public:
    FrameBuffer(FrameBufferCore& other);
    void Render(const int tile_w, const int tile_h,
                const int tile_x, const int tile_y,
                const int global_stride, 
                const float* global_depth,
                OSPRenderer renderer,
                float*& dest);
  };

  /**
   * Context
   */
  struct Context {
  public:
    ContextCore* core;
    Context(ContextCore& other);
    void InitPatch(const int patchID);
    void SetupPatch(const int patchID,
                    const OSPDataType data_type, 
                    const std::string data_char,
                    const size_t data_size, 
                    const void* data_ptr,
                    const double *X, const double *Y, const double *Z, 
                    const int nX, const int nY, const int nZ,
                    const double dbox[6], const double cbox[6]);
    void RenderPatch(const int patchID,
                     const float xMin, const float xMax, 
                     const float yMin, const float yMax,
                     const int tile_w, const int tile_h,
                     float*& dest); 
  };

};
};
