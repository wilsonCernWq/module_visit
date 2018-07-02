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
#include <vector>

namespace ospray {
namespace visit {

  /**
   * Shared Variables
   */
  struct SharedVariables {
    
  };
  
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
    double intensity;
    double color[3];
    double direction[3];
    LightCore() : Object<OSPLight>() {
      isAmbient = false;
      intensity = 1.0;
      color[0] = color[1] = color[2] = 1.0;
      direction[0] = direction[1] = 0.0;
      direction[2] = 1.0;
    }
  };  

  /**
   * Renderer Wrapper
   */
  struct RendererCore : public Object<OSPRenderer> {
    OSPData                lightData;
    std::vector<LightCore> lightList;
    int  aoSamples;
    int  spp;
    bool enableOneSidedLighting;
    bool enableShadowsEnabled;
    bool enableAoTransparencyEnabled;
    // (shared, never delete)
    const unsigned char *bgColorBuffer;  // color channel for the backplatte
    const float         *bgDepthBuffer;  // depth channel for the backplatte 
    int bgSize[2]; // buffer extents (minX, maxX, minY, max)  
    RendererCore() : Object<OSPRenderer>() {
      lightData = NULL;
      aoSamples = 0;
      spp       = 1;
      enableOneSidedLighting      = false;
      enableShadowsEnabled        = false;
      enableAoTransparencyEnabled = false;
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
  VolumeCore() : Object<OSPVolume>() {}
  };
  
};
};
