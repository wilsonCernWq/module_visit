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

// ************************************************************************
// We expose this in header because iy will be called in other components
// where we dont have direct library linkage
// ************************************************************************

#pragma once

#include <ospray/ospray.h>
#include <string>
#include <vector>
#include <map>

using std::string;

namespace ospray {
namespace visit {
  
  /** 
   * Helper Functions
   */
  template<typename T> void ospray_rm(T& obj) {
    if (!obj) { ospRelease(obj); obj = NULL; }
  }
  
  /**
   * Abstraction of an object
   */
  template<typename T> struct Object {
    bool init;
    T    self;
    Object() { init = false; self = NULL; }
    virtual ~Object() { ospray_rm(self); init = false; }
    T operator*() { return self; }
  };
  
  /**
   * Transfer Function Wrapper
   */
  struct TransferFunctionCore : public Object<OSPTransferFunction> {
    TransferFunctionCore() : Object<OSPTransferFunction>() {}
  };

  /**
   * Camera Wrapper
   */
  struct CameraCore : public Object<OSPCamera> {
    bool   orthographic;
    int    windowExts[4];
    int    screenSize[2];
    double pan[2]; // pan ratio [0, 1]
    double zoom;   // zoom factor
    CameraCore() : Object<OSPCamera>() {
      orthographic = false;
      windowExts[0] = windowExts[1] = 0;
      windowExts[2] = windowExts[3] = 0;
      screenSize[0] = screenSize[1] = 0;
      pan[0] = pan[1] = 0.0;
      zoom = 1.0;
    }
  };

  /**
   * Light Wrapper
   */
  struct LightCore : public Object<OSPLight> {
    bool isAmbient;
    LightCore() : Object<OSPLight>() {
      isAmbient = false;
    }
  };  

  /**
   * Renderer Wrapper
   */
  struct RendererCore : public Object<OSPRenderer> {
    OSPData                lightData;
    std::vector<LightCore> lightList;
    RendererCore() : Object<OSPRenderer>() {
      lightData = NULL;
    }
    ~RendererCore() { ospray_rm(lightData); }
  };

  /**
   * Model Wrapper
   */
  struct ModelCore : public Object<OSPModel> {
    ModelCore() : Object<OSPModel>() {}
  };

  /**
   * Volume Wrapper
   */
  struct VolumeCore : public Object<OSPVolume> {
    string      volumeType;
    OSPDataType dataType;
    size_t      dataSize;
    const void* dataPtr;
    VolumeCore() : Object<OSPVolume>() {
      volumeType = "";
      dataType = OSP_UCHAR; /* just give it a value */
      dataSize = 0;
      dataPtr  = NULL;
    }
  };
  
  /**
   * Framebuffer Wrapper
   */
  struct FrameBufferCore : public Object<OSPFrameBuffer> {
    FrameBufferCore() : Object<OSPFrameBuffer>() {}
  };

  /**
   * Now we define a PatchCore
   */
  struct Patch {
    VolumeCore      volume;
    ModelCore       model;
    FrameBufferCore fb;
  };

  /**
   * And a ContextCore
   */
  struct ContextCore {
    // data
    string varname;
    std::map<int, Patch> patches;
    CameraCore           camera;
    RendererCore         renderer;
    TransferFunctionCore tfn;
    // flags
    bool oneSidedLighting;       /* renderer */
    bool shadowsEnabled;         /* renderer */
    bool aoTransparencyEnabled;  /* renderer */
    bool useGridAccelerator;     /*  volume  */
    bool adaptiveSampling;       /*  volume  */
    bool preIntegration;         /*  volume  */
    bool singleShade;            /*  volume  */
    bool gradientShadingEnabled; /*  volume  */
    // other parameters
    double Ks;
    double Ns;
    double samplingRate;
    int aoSamples;
    int spp;
    double scale[3];
    double gbbox[6];
    // (shared, dont delete here)
    const unsigned char *bgColorBuffer;  // color channel for the backplatte
    const float         *bgDepthBuffer;  // depth channel for the backplatte 
    int                  bgSize[2];      // channel buffer size
    ContextCore() {
      varname = "";
      oneSidedLighting       = false;
      shadowsEnabled         = false;
      aoTransparencyEnabled  = false;
      useGridAccelerator     = false;
      adaptiveSampling       = false;
      preIntegration         = false;
      singleShade            = false;
      gradientShadingEnabled = false;
      Ks = 1.0; Ns = 20;
      samplingRate = 3.0;
      aoSamples = 0;
      spp = 1;
      scale[0] = scale[1] = scale[2] = 1.f;
      gbbox[0] = gbbox[1] = gbbox[2] = 0.f;
      gbbox[3] = gbbox[4] = gbbox[5] = 0.f;
      bgSize[0] = bgSize[1] = 0;
    }
  };

};
};
