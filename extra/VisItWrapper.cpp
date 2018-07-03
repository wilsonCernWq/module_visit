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

#include "VisItWrapper.h"
#include "VisItModuleCommon.h"

#include <ospcommon/vec.h>
#include <cmath>
#include <vector>

using namespace ospcommon;

namespace ospray {
namespace visit {

  template<typename T> 
  void ospray_check(const T& obj, const std::string s) {
    if (!obj) { 
      if (CheckVerbose()) {
        std::cerr << s << " is invalid" << std::endl; 
      }
      throw std::runtime_error(s + " is invalid");
    }
  }

  // =======================================================================//
  //
  // =======================================================================//
  struct Color { float R,G,B,A; };

  // =======================================================================//
  //
  // =======================================================================//
  TransferFunction::TransferFunction(TransferFunctionCore& other)
    : Manipulator<CoreType, OSPType>(other) {}
  void TransferFunction::Set(const void *_table,
                             const unsigned int size, 
                             const double datamin, 
                             const double datamax) 
  {
    // initialize it once
    if (!core->init) 
    {
      ospray_rm(core->self);
      core->self = ospNewTransferFunction("piecewise_linear");
      ospray_check(core->self, "transfer function");
      core->init = true;
    }
    // create OSP data
    const Color* table = reinterpret_cast<const Color*>(_table);
    std::vector<vec3f> cdata;
    std::vector<float> odata;
    for (unsigned int i = 0; i < size; ++i)
    {
      cdata.emplace_back(table[i].R, table[i].G, table[i].B);
      odata.emplace_back(table[i].A);
    }
    OSPData osp_cdata = ospNewData(cdata.size(), OSP_FLOAT3, cdata.data());
    OSPData osp_odata = ospNewData(odata.size(), OSP_FLOAT,  odata.data());
    ospray_check(osp_cdata,"TFN color data");
    ospray_check(osp_odata,"TFN opacity data");
    // commit
    ospSetData(core->self, "colors",    osp_cdata);
    ospSetData(core->self, "opacities", osp_odata);
    ospSet2f(core->self, "valueRange", datamin, datamax);
    ospCommit(core->self);
    // cleanup
    ospray_rm(osp_cdata);
    ospray_rm(osp_odata);
  }

  // =======================================================================//
  //
  // =======================================================================//
  Camera::Camera(CameraCore& other) 
    : Manipulator<CoreType, OSPType>(other) {}
  void Camera::Set(const bool ortho,
                   const double camera_p[3], 
                   const double camera_f[3], 
                   const double camera_u[3], 
                   const double fovy,
                   const double pan_ratio[2],
                   const double zoom_ratio,
                   const double near_clip,
                   const double canvas_size[2],
                   const int screen_size[2], /* overall screen size */
                   const int tile_extents[4] /*      tile size      */)
  {
    // create camera
    if (!core->init || (ortho != core->orthographic))
    {
      core->orthographic = ortho;
      ospray_rm(core->self);
      core->self = core->orthographic ? 
        ospNewCamera("orthographic") : 
        ospNewCamera("perspective");
      ospray_check(core->self, "camera");
      core->init = true;
    }
    // compute camera
    core->zoom = zoom_ratio;
    core->pan[0] = pan_ratio[0];
    core->pan[1] = pan_ratio[1];
    core->screenSize[0] = screen_size[0];
    core->screenSize[1] = screen_size[1];
    // commit
    ospSet1f(core->self, "nearClip", near_clip);
    ospSet3f(core->self, "pos", 
             camera_p[0], 
             camera_p[1], 
             camera_p[2]);
    ospSet3f(core->self, "dir", 
             camera_f[0] - camera_p[0], 
             camera_f[1] - camera_p[1], 
             camera_f[2] - camera_p[2]);
    ospSet3f(core->self, "up",
             camera_u[0],
             camera_u[1], 
             camera_u[2]);
    ospSet1f(core->self, "aspect", 
             static_cast<double>(screen_size[0]) /
             static_cast<double>(screen_size[1]));
    if (!ortho) ospSet1f(core->self, "fovy", fovy);
    else ospSet1f(core->self, "height", canvas_size[1]);
    SetScreen(tile_extents[0], tile_extents[1],
              tile_extents[2], tile_extents[3]);
  }
  void Camera::SetScreen(const double xMin, const double xMax,
                         const double yMin, const double yMax) 
  {
    core->windowExts[0] = std::max(static_cast<int>(std::round(xMin)), 0);
    core->windowExts[1] = std::min(static_cast<int>(std::round(xMax)), 
                                   core->screenSize[0]);
    core->windowExts[2] = std::max(static_cast<int>(std::round(yMin)), 0);
    core->windowExts[3] = std::min(static_cast<int>(std::round(yMax)), 
                                   core->screenSize[1]);
    const double r_xl = xMin/core->screenSize[0] - core->pan[0];
    const double r_yl = yMin/core->screenSize[1] - core->pan[1];
    const double r_xu = xMax/core->screenSize[0] - core->pan[0];
    const double r_yu = yMax/core->screenSize[1] - core->pan[1];
    ospSet2f(core->self, "imageStart",
             (r_xl - 0.5) / core->zoom + 0.5,
             (r_yl - 0.5) / core->zoom + 0.5);
    ospSet2f(core->self, "imageEnd",
             (r_xu - 0.5) / core->zoom + 0.5,
             (r_yu - 0.5) / core->zoom + 0.5);
    ospCommit(core->self);
  }

  // =======================================================================//
  //
  // =======================================================================//
  Light::Light(LightCore& other) 
    : Manipulator<CoreType, OSPType>(other) {}
  void Light::Set(const bool ambient, const double i, 
                  const double c, const double* d)
  {
    Set(ambient, i, c, c, c, d);
  }
  void Light::Set(const bool ambient, const double i, 
                  const double cr, const double cg, const double cb,
                  const double* d)
  {
    const double c[3] = {cr, cg, cb};
    Set(ambient, i, c, d);
  }
  void Light::Set(const bool ambient, const double i, 
                  const double c[3], const double* d)
  {
    // create light
    if (!core->init || ambient != core->isAmbient) {
      ospray_rm(core->self);
      core->self = ospNewLight2("scivis", "distant");
      ospray_check(core->self, "light");
      ospSet1i(core->self, "isVisible", 0);
      ospSet1f(core->self, "angularDiameter", 0.53f);
      core->init = true;
    }
    // commit light
    ospSet1f(core->self, "intensity", i);
    ospSet3f(core->self, "color", c[0], c[1], c[2]);
    if (!ambient) {
      ospSet3f(core->self, "direction", d[0], d[1], d[2]);
    }
    ospCommit(core->self);
  }

  // =======================================================================//
  //
  // =======================================================================//
  Renderer::Renderer(RendererCore& other) 
    : Manipulator<CoreType, OSPType>(other) {}
  void Renderer::Init()
  {
    if (!core->init) {
      ospray_rm(core->self);
      core->self = ospNewRenderer("scivis");
      ospray_check(core->self, "SCIVIS Renderer");
      core->init = true;
    }
  }
  void Renderer::ResetLights()
  {
    if (!core->init) { Init(); }
    ospray_rm(core->lightData);
    core->lightList.clear();
  }
  Light Renderer::AddLight()
  {
    core->lightList.emplace_back();
    return Light(core->lightList.back());
  }
  void Renderer::FinalizeLights()
  {
    std::vector<OSPLight> osp_light_list;
    for (auto& l : core->lightList) { osp_light_list.emplace_back(*l); }
    ospray_rm(core->lightData);
    core->lightData = ospNewData(osp_light_list.size(), OSP_OBJECT, 
                                 osp_light_list.data());
    ospray_check(core->lightData, "light list");
  }
  void Renderer::SetBackgroundBuffer(const unsigned char* color, 
                                     const float* depth, 
                                     const int size[2])
  {
    core->bgColorBuffer = color;
    core->bgDepthBuffer = depth;
    core->bgSize[0] = size[0];
    core->bgSize[1] = size[1];
  }
  void Renderer::Set(const int aoSamples, const int spp, 
                     const bool oneSidedLighting,
                     const bool shadowsEnabled,
                     const bool aoTransparencyEnabled)
  {
    if (!core->init) { Init(); }   
    ospSet1i(core->self, "aoSamples", aoSamples);
    ospSet1i(core->self, "spp", spp);
    ospSet1i(core->self, "oneSidedLighting", oneSidedLighting);
    ospSet1i(core->self, "shadowsEnabled", shadowsEnabled);
    ospSet1i(core->self, "aoTransparencyEnabled", aoTransparencyEnabled);
    ospSetData(core->self, "lights", core->lightData);
    ospCommit(core->self);
  }
  void Renderer::Set(OSPCamera osp_camera)
  {
    if (!core->init) { Init(); }
    ospSetObject(core->self, "camera", osp_camera);
    ospCommit(core->self);
  }
  void Renderer::Set(OSPModel osp_world)
  {
    if (!core->init) { Init(); }
    ospSetObject(core->self, "model", osp_world);
    ospCommit(core->self);
  }
    
  // =======================================================================//
  //
  // =======================================================================//
  Model::Model(ModelCore& other)
    : Manipulator<CoreType, OSPType>(other) {}
  void Model::Reset()
  {
    core->init = false;
  }
  void Model::Init()
  {
    if (!core->init) {
      ospray_rm(core->self);
      core->self = ospNewModel();
      ospray_check(core->self, "Model");
      core->init = true;
    }
  }
  void Model::Set(OSPVolume volume)
  {
    if (!core->init) { Init(); }
    ospAddVolume(core->self, volume);
    ospCommit(core->self);
  }

  // =======================================================================//
  //
  // =======================================================================//
  Volume::Volume(VolumeCore& other)
    : Manipulator<CoreType, OSPType>(other) {}
  bool Volume::Init(const std::string volume_type, 
                    const OSPDataType data_type, const std::string data_char,
                    const size_t data_size, const void* data_ptr)
  {
    if (!core->init || 
        volume_type != core->volumeType ||
        data_type   != core->dataType   ||
        data_size   != core->dataSize   ||
        data_ptr    != core->dataPtr)
    {
      core->volumeType = volume_type;
      core->dataType = data_type;
      core->dataSize = data_size;
      core->dataPtr  = data_ptr;
      ospray_rm(core->self);
      core->self = ospNewVolume(volume_type.c_str());
      ospray_check(core->self, volume_type);
      if (volume_type == "visit_shared_structured_volume" ||
          volume_type == "shared_structured_volume") 
      {
        OSPData osp_data = ospNewData(data_size, data_type,
                                      data_ptr, OSP_DATA_SHARED_BUFFER);
        ospSetString(core->self, "voxelType", data_char.c_str());
        ospSetData(core->self, "voxelData", osp_data);
        ospray_rm(osp_data);
      }
      core->init = true;
      return true;
    }
    return false;
  }
  void Volume::Set(const bool useGridAccelerator, const bool adaptiveSampling,
                   const bool preIntegration, const bool singleShade, 
                   const bool gradientShadingEnabled, const double samplingRate, 
                   const double Ks, const double Ns,
                   const double *X, const double *Y, const double *Z, 
                   const int nX, const int nY, const int nZ,
                   const double dbox[6], const double cbox[6], 
                   const osp::vec3f& global_upper,
                   const osp::vec3f& global_lower,
                   const osp::vec3f& scale,
                   OSPTransferFunction tfn)
  {
    const vec3i dims(nX, nY, nZ);
    const vec3f data_lower(vec3f(dbox[0], dbox[1], dbox[2]) * (const vec3f&)scale);
    const vec3f data_upper(vec3f(dbox[3], dbox[4], dbox[5]) * (const vec3f&)scale);
    const vec3f clip_lower(vec3f(cbox[0], cbox[1], cbox[2]) * (const vec3f&)scale);
    const vec3f clip_upper(vec3f(cbox[3], cbox[4], cbox[5]) * (const vec3f&)scale);
    const vec3f spacing = (data_upper - data_lower)/((const vec3f)dims - 1.0f);
    ospSetVec3f(core->self, "volumeGlobalBoundingBoxLower", 
                (const osp::vec3f&) global_upper);
    ospSetVec3f(core->self, "volumeGlobalBoundingBoxUpper",
                (const osp::vec3f&) global_lower);
    ospSetVec3f(core->self, "volumeClippingBoxLower", (const osp::vec3f&)clip_lower);
    ospSetVec3f(core->self, "volumeClippingBoxUpper", (const osp::vec3f&)clip_upper);
    ospSetVec3f(core->self, "gridSpacing", (const osp::vec3f&)spacing);
    ospSetVec3f(core->self, "gridOrigin",  (const osp::vec3f&)data_lower);
    ospSetVec3i(core->self, "dimensions",  (const osp::vec3i&)dims);
    ospSet1f(core->self, "samplingRate", samplingRate);
    ospSet3f(core->self, "Ks", Ks, Ks, Ks);
    ospSet1f(core->self, "Ns", Ns);
    ospSet1i(core->self, "gradientShadingEnabled", (int)gradientShadingEnabled);
    ospSet1i(core->self, "useGridAccelerator", (int)useGridAccelerator);
    ospSet1i(core->self, "adaptiveSampling", (int)adaptiveSampling);
    ospSet1i(core->self, "preIntegration", (int)preIntegration);
    ospSet1i(core->self, "singleShade", (int)singleShade);
    ospSetObject(core->self, "transferFunction", tfn);
    ospCommit(core->self);
  }  
  void Volume::ComputeGhostBounds(bool r[6], 
                                  const unsigned char *g, 
                                  const int nX, 
                                  const int nY, 
                                  const int nZ)
  {
    for (int y = 1; y < (nY-1); ++y) {
      for (int z = 1; z < (nZ-1); ++z) {
        if (!r[0]) if (g[z*nY*nX+y*nX       ] != 0) r[0] = true;
        if (!r[3]) if (g[z*nY*nX+y*nX+(nX-1)] != 0) r[3] = true;
        if (r[0] && r[3]) { break; }
      }
    }
    for (int x = 1; x < (nX-1); ++x) {
      for (int z = 1; z < (nZ-1); ++z) {
        if (!r[1]) if (g[z*nY*nX          +x] != 0) r[1] = true;
        if (!r[4]) if (g[z*nY*nX+(nY-1)*nX+x] != 0) r[4] = true;
        if (r[1] && r[4]) { break; }
      }
    }
    for (int x = 1; x < (nX-1); ++x) {
      for (int y = 1; y < (nY-1); ++y) {
        if (!r[2]) if (g[             y*nX+x] != 0) r[2] = true;
        if (!r[5]) if (g[(nZ-1)*nY*nX+y*nX+x] != 0) r[5] = true; 
        if (r[2] && r[5]) { break; }
      }
    }
  }

  // =======================================================================//
  //
  // =======================================================================//
  FrameBuffer::FrameBuffer(FrameBufferCore& other)
    : Manipulator<CoreType, OSPType>(other) {}
  void FrameBuffer::Render(const int tile_w, const int tile_h,
                           const int tile_x, const int tile_y,
                           const int    global_stride, 
                           const float* global_depth,
                           OSPRenderer renderer,
                           float*& dest)
  {
    const vec2i fb_size(tile_w, tile_h);
    // prepare the maxDepthDexture
    {
      // The reason I use round(r * (N-1)) instead of floor(r * N) is that
      // during the composition phase, there will be a wired offset between
      // rendered image and the background, which is about one pixel in size.
      // Using round(r * (N - 1)) can remove the problem
      std::vector<float> local_depth(tile_w * tile_h);
      for (int i = 0; i < tile_w; ++i) {
        for (int j = 0; j < tile_h; ++j) {
          local_depth[i + j * tile_w] = 
            global_depth[tile_x + i + (tile_y + j) * global_stride];
        }
      }
      OSPTexture2D maxDepthTexture
        = ospNewTexture2D((const osp::vec2i&)fb_size, OSP_TEXTURE_R32F,
                          local_depth.data(), OSP_TEXTURE_FILTER_NEAREST);    
      ospCommit(maxDepthTexture);
      ospSetObject(renderer, "maxDepthTexture", maxDepthTexture);
      ospCommit(renderer);
      ospray_rm(maxDepthTexture);
    }    
    // do the rendering
    // ALWAYS create a new framebuffer
    {
      ospray_rm(core->self);
      core->self = ospNewFrameBuffer((const osp::vec2i&)fb_size, OSP_FB_RGBA32F, OSP_FB_COLOR);
      ospray_check(core->self, "framebuffer");
      ospRenderFrame(core->self, renderer, OSP_FB_COLOR);
      const float* image = (float*)ospMapFrameBuffer(core->self, OSP_FB_COLOR);
      std::copy(image, image + (tile_w * tile_h) * 4, dest);
      ospUnmapFrameBuffer(image, core->self);
      ospray_rm(core->self);      
    }
  }

  // =======================================================================//
  //
  // =======================================================================//
  Context::Context(ContextCore& other)
    : core{&other} {}
  void Context::InitPatch(const int patchID)
  {
    if (core->patches.find(patchID) == core->patches.end()) {
        core->patches[patchID] = PatchCore();
    }
  }
  void Context::SetupPatch(const int patchID,
                           const OSPDataType data_type, 
                           const std::string data_char,
                           const size_t data_size, 
                           const void* data_ptr,
                           const double *X, const double *Y, const double *Z, 
                           const int nX, const int nY, const int nZ,
                           const double dbox[6], const double cbox[6])

  {
    Volume volume(core->patches[patchID].volume);
    Model  model(core->patches[patchID].model);
    bool reset = volume.Init("visit_shared_structured_volume",
                             data_type, data_char, data_size, data_ptr);
    volume.Set(core->useGridAccelerator, core->adaptiveSampling,
               core->preIntegration, core->singleShade,
               core->gradientShadingEnabled,
               core->samplingRate, 
               core->specularKs, core->specularNs,
               X, Y, Z, nX, nY, nZ, dbox, cbox,
               core->bbox.upper,
               core->bbox.lower,
               core->scale,
               *(core->tfn));
    if (reset) {
      model.Reset();
      model.Init();
      model.Set(*(core->patches[patchID].volume));
    }
  }
  void Context::RenderPatch(const int patchID,
                            const float xMin, const float xMax, 
                            const float yMin, const float yMax,
                            const int tile_w, const int tile_h,
                            float*& dest)
  {
    Camera camera(core->camera);
    Renderer renderer(core->renderer);
    FrameBuffer fb(core->patches[patchID].fb);
    camera.SetScreen(xMin, xMax, yMin, yMax);
    renderer.Set(*(core->patches[patchID].model));
    renderer.Set(*(core->camera));
    fb.Render(tile_w, tile_h, camera.GetWindowExts(0), camera.GetWindowExts(2),
              core->renderer.bgSize[0], core->renderer.bgDepthBuffer,
              *(core->renderer), dest);
  }
  
};
};
