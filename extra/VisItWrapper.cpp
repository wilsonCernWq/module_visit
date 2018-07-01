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
    : core{&other} {}
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
    : core{&other} {}
  void Camera::Set(const bool ortho,
                   const double camera_p[3], 
                   const double camera_f[3], 
                   const double camera_u[3], 
                   const double fovy,
                   const double pan_ratio[2],
                   const double zoom_ratio, 
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
    core->pan[0] = pan_ratio[0] * zoom_ratio;
    core->pan[1] = pan_ratio[1] * zoom_ratio;
    core->screenSize[0] = screen_size[0];
    core->screenSize[1] = screen_size[1];
    // commit
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
    : core{&other} {}
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
    core->intensity = i;
    ospSet1f(core->self, "intensity", i);
    core->color[0] = c[0];
    core->color[1] = c[1];
    core->color[2] = c[2];
    ospSet3f(core->self, "color", c[0], c[1], c[2]);
    if (!ambient) {
      core->direction[0] = d[0];
      core->direction[1] = d[1];
      core->direction[2] = d[2];
      ospSet3f(core->self, "direction", d[0], d[1], d[2]);
    }
    ospCommit(core->self);
  }

  // =======================================================================//
  //
  // =======================================================================//
  Renderer::Renderer(RendererCore& other) 
    : core{&other} {}
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
  void Renderer::Set(const int ao_samples, const int spp, 
                     const bool flag_oneSidedLighting,
                     const bool flag_shadowsEnabled,
                     const bool flag_aoTransparencyEnabled)
  {
    if (!core->init) { Init(); }
    core->aoSamples = ao_samples;
    core->spp       = spp;
    core->enableOneSidedLighting      = flag_oneSidedLighting;
    core->enableShadowsEnabled        = flag_shadowsEnabled;
    core->enableAoTransparencyEnabled = flag_aoTransparencyEnabled;
    ospSet1i(core->self, "aoSamples", ao_samples);
    ospSet1i(core->self, "spp", spp);
    ospSet1i(core->self, "oneSidedLighting", flag_oneSidedLighting);
    ospSet1i(core->self, "shadowsEnabled", flag_shadowsEnabled);
    ospSet1i(core->self, "aoTransparencyEnabled", flag_aoTransparencyEnabled);
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
};
};



/*

void 
OSPVisItPatches::Volume::SetupParameters(const int   voxelEnumVTK,
					 const void *voxelDataPtr,
					 const double *X,
					 const double *Y,
					 const double *Z, 
					 const int    nX,
					 const int    nY,
					 const int    nZ,
					 const double pbox[6],
					 const double bbox[6])
{
    // determine volume data type
    std::string voxelTypeOSP;
    std::string voxelEnumOSP;
    ospout << "[avtOSPRayCommon] current voxel type: ";
    switch (voxelEnumVTK) {
    case (VTK_UNSIGNED_CHAR):
	voxelTypeOSP = "uchar";
	voxelEnumOSP = OSP_UCHAR;
	break;
    case (VTK_SHORT):
	voxelTypeOSP = "short";
	voxelEnumOSP = OSP_SHORT;
	break;
    case (VTK_UNSIGNED_SHORT):
	voxelTypeOSP = "ushort";
	voxelEnumOSP = OSP_USHORT;
	break;
    case (VTK_FLOAT):
	voxelTypeOSP = "float";
	voxelEnumOSP = OSP_FLOAT;
	break;
    case (VTK_DOUBLE):
	voxelTypeOSP = "double";
	voxelEnumOSP = OSP_DOUBLE;
	break;
    default:
	ospray::Exception("[avtOSPRayCommon] ERROR: "
			  "Unsupported voxel type "
			  + std::to_string(voxelEnumVTK) +
			  "found.");
    }
    ospout << voxelTypeOSP << std::endl;

    // create volume and model
    if (true)
    {
	OSPRAY_SAFE_RM(model);
	model = ospNewModel();
	OSPRAY_CHECK(model, "");   
	OSPRAY_SAFE_RM(volume);
	volume = ospNewVolume("visit_shared_structured_volume");
	OSPRAY_CHECK(volume, "visit_shared_structured_volume");
	OSPRAY_SAFE_RM(voxels);
	voxels = ospNewData(nX * nY * nZ, voxelEnumOSP,
			    voxelDataPtr, OSP_DATA_SHARED_BUFFER);
	OSPRAY_CHECK(voxels,"");
	initialized = true;
    }

    // commit volume
    ospSetObject(volume, "transferFunction", *(parent->tfn));
    ospSetString(volume, "voxelType", voxelTypeOSP.c_str());
    ospSetData(volume, "voxelData", voxels);
    ospSet1i(volume, "gradientShadingEnabled", parent->enableShading);
    ospSet1i(volume, "useGridAccelerator", parent->enableGridAccelerator);
    ospSet1i(volume, "preIntegration", parent->enablePreIntegration);
    ospSet1i(volume, "adaptiveSampling", 0);
    ospSet1i(volume, "singleShade", 0);
    ospSet1f(volume, "samplingRate", parent->samplingRate);
    ospSet3f(volume, "Ks", parent->Ks, parent->Ks, parent->Ks);
    ospSet1f(volume, "Ns", parent->Ns);
    ospSet3i(volume, "dimensions", nX, nY, nZ);
    ospSet3f(volume, "gridOrigin",
	     pbox[0] * parent->scaling[0],
	     pbox[1] * parent->scaling[1],
	     pbox[2] * parent->scaling[2]);
    ospSet3f(volume, "gridSpacing",
	     parent->scaling[0]*(pbox[3]-pbox[0])/((double)nX-1.),
	     parent->scaling[1]*(pbox[4]-pbox[1])/((double)nY-1.),
	     parent->scaling[2]*(pbox[5]-pbox[2])/((double)nZ-1.));
    ospSet3f(volume, "volumeClippingBoxLower",
	     bbox[0] * parent->scaling[0],
	     bbox[1] * parent->scaling[1],
	     bbox[2] * parent->scaling[2]);
    ospSet3f(volume, "volumeClippingBoxUpper",
	     bbox[3] * parent->scaling[0],
	     bbox[4] * parent->scaling[1],
	     bbox[5] * parent->scaling[2]);
    ospSet3f(volume, "volumeGlobalBoundingBoxLower",
	     scaledGlobalBBoxLower[0],
	     scaledGlobalBBoxLower[1],
	     scaledGlobalBBoxLower[2]);
    ospSet3f(volume, "volumeGlobalBoundingBoxUpper",
	     scaledGlobalBBoxUpper[0],
	     scaledGlobalBBoxUpper[1],
	     scaledGlobalBBoxUpper[2]);
    ospCommit(volume);
    ospAddVolume(model, volume);
    ospCommit(model);
}

void 
OSPVisItPatches::Volume::PrepareToRender(OSPRenderer renderer,
					 const unsigned int W,
					 const unsigned int H,
					 const int Xs, const int Ys,
					 const float* maxDepthBuffer,
					 const int    maxDepthStride)
{
    // create maximum depth texture
    const osp::vec2i fbSize;
    fbSize.x = W; fbSize.y = H; pixelnum = W * H;
    std::vector<float> bgData(pixelnum);
    for (int i = 0; i < W; ++i) {
    	for (int j = 0; j < H; ++j) {
    	    maxDepth[i + j * W] =
		maxDepthBuffer[Xs + i + (Ys + j) * maxDepthStride];
    	}
    }

    // create background texture
    OSPRAY_SAFE_RM(background);
    background = ospNewTexture2D(fbSize, OSP_TEXTURE_R32F, 
				 maxDepth.data(),
				 OSP_TEXTURE_FILTER_NEAREST);
    OSPRAY_CHECK(background,"");
    ospCommit(background);

    // create framebuffer
    ospSetObject(renderer, "maxDepthTexture", background);
    ospCommit(renderer);
    OSPRAY_SAFE_RM(background);

    // create framebuffer
    OSPRAY_SAFE_RM(framebuffer);
    framebuffer = ospNewFrameBuffer(fbSize, OSP_FB_RGBA32F, OSP_FB_COLOR);
    OSPRAY_CHECK(framebuffer,"");
}

void
OSPVisItPatches::Volume::Render(OSPRenderer renderer, float* output)
{
    StackTimer t0("OSPVisItPatches Renders One Patch");
    {
	std::cout << "OSPVisItPatches::ospRenderFrame start" << std::endl;
	ospRenderFrame(framebuffer, renderer, OSP_FB_COLOR);
	std::cout << "OSPVisItPatches::ospRenderFrame done" << std::endl;
    }
    float* pixels = (float*) ospMapFrameBuffer(framebuffer, OSP_FB_COLOR);
    std::copy(pixels, pixels + pixelnum * 4, output);
    ospUnmapFrameBuffer(pixels, framebuffer);	    
}

// ****************************************************************************
//
// OSPCamera
//
// ****************************************************************************
// ****************************************************************************
//
// OSPVisItRenderer
//
// ****************************************************************************

void OSPVisItRenderer::Light::SetIntensity(const double v)
{
    intensity = v;
    ospSet1f(self, "intensity", intensity);
    ospCommit(self);
}

void OSPVisItRenderer::Light::SetColor(const double v[3])
{
    color[0] = v[0];
    color[1] = v[1];
    color[2] = v[2];
    ospSet3f(self, "color", color[0], color[1], color[2]);
    ospCommit(self);
}

void OSPVisItRenderer::AmbientLight::Init()
{
    if (!initialized) {
	OSPRAY_SAFE_RM(self);
	self = ospNewLight2("scivis", "ambient");
	OSPRAY_CHECK(self,"");
	ospSet1i(aLight, "isVisible", 0);
	ospCommit(self);
	initialized = true;
    }
}

void OSPVisItRenderer::DistantLight::Init()
{
    if (!initialized) {
	OSPRAY_SAFE_RM(self);
	self = ospNewLight2("scivis", "distant");
	OSPRAY_CHECK(self,"");
	ospSet1i(aLight, "isVisible", 0);
	ospSet1f(sLight, "angularDiameter", 0.53f);	    
	ospCommit(self);
	initialized = true;
    }
}

void OSPVisItRenderer::DistantLight::SetDirection(const double v[3])
{
    direction[0] = v[0];
    direction[1] = v[1];
    direction[2] = v[2];
    ospSet3f(self, "direction", direction[0], direction[1], direction[2]);
    ospCommit(self);
}
void OSPVisItRenderer::Init() 
{
    if (!initialized) {
	OSPRAY_SAFE_RM(self);
	self = ospNewRenderer("scivis");
	OSPRAY_CHECK(self,"scivis rencerer");
	initialized = true;
    }
}
void OSPVisItRenderer::SetMaxDepthBuffer(float* b, const int x, const int y)
{
    maxDepthBuffer = b;
    maxDepthSize[0] = x;
    maxDepthSize[1] = y;
}
void OSPVisItRenderer::SetCamera(OSPCamera camera)
{
    ospSetObject(self, "camera", camera);
}
void OSPVisItRenderer::SetModel(OSPModel world)
{
    ospSetObject(self, "model",  world);
}
void OSPVisItRenderer::SetLight(double v[3], const double a, const double d)
{
    aLight.SetIntensity(a);
    dLights.resize(2);
    dLights[0].SetDirection(v);
    dLights[0].SetIntensity(d);   /// diffuse light 
    dLights[1].SetDirection(v);
    dLights[1].SetIntensity(1.5); // sun light 
    std::vector<OSPLight> list;
    list.push_back(*aLight);
    for (int i = 0; i < dLights.size(); ++i) {
	list.pusb_back(*(dLights[i]));
    }
    OSPRAY_SAFE_RM(lightData);
    lightData = ospNewData(list.size(), OSP_OBJECT, list.data());
    OSPRAY_CHECK(lightData);
    ospSetData(self, "lights", lightData);
}
void OSPVisItRenderer::Commit()
{    
    ospSet1f(self, "bgColor",   0.f);
    ospSet1i(self, "aoSamples", aoSamples);
    ospSet1i(self, "spp",       spp);
    ospSet1i(self, "oneSidedLighting",      enableOneSidedLighting);
    ospSet1i(self, "shadowsEnabled",        enableShadowsEnabled);
    ospSet1i(self, "aoTransparencyEnabled", enableAoTransparencyEnabled);
    ospCommit(self);
}

// ****************************************************************************
//  Struct:  OSPContext
//
//  Purpose:
//
//
//  Programmer:  
//  Creation:   
//
// ****************************************************************************

void OSPContext_ErrorFunc(OSPError, const char* msg) { 
    osperr << "#osp: (rank " << PAR_Rank() << ")" 
           << msg; 
}
void OSPContext_StatusFunc(const char* msg) { 
    osperr << "#osp: (rank " << PAR_Rank() << ")" 
           << msg; 
}
bool OSPVisItContext::initialized = false;
void OSPVisItContext::Init(int numThreads) 
{ 
    if (!OSPVisItContext::initialized) 
    {
	// check hostname
#ifdef __unix__
	char hname[200];
	gethostname(hname, 200);
        ospout << "[ospray] on host >> " << hname << "<<" << std::endl;;
#endif
	// initialize ospray
        ospout << "[ospray] Initialize OSPRay";
	OSPDevice device = ospNewDevice();
	OSPRAY_CHECK(device,"");
	if (!device)
	{
	    osperr << "[avtOSPRayCommon] ERROR: "
		   << "can't create ospray device"
		   << std::endl;
	}
	// setup debug && number of threads (this can only be hard-coded)
	if (DebugStream::Level5())
	{
	    ospout << " debug mode";
	    ospDeviceSet1i(device, "debug", 0);
	}
	if (numThreads > 0)
	{
	    ospout << " numThreads: " << numThreads;
	    ospDeviceSet1i(device, "numThreads", numThreads);
	}
	ospout << std::endl;	
	ospDeviceSetErrorFunc(device, OSPContext_ErrorFunc);
	ospDeviceSetStatusFunc(device, OSPContext_StatusFunc);
	ospDeviceCommit(device);
	ospSetCurrentDevice(device);
	// load ospray module
	OSPError err = ospLoadModule("visit");
	if (err != OSP_NO_ERROR)
	{
	    osperr << "[avtOSPRayCommon] ERROR: "
		   << "can't load visit module"
		   << std::endl;
	}
    }
    OSPVisItContext::initialized = true;
}

// We use this function to minimize interface
void OSPVisItContext::Render(float xMin, float xMax, float yMin, float yMax,
			     int imgWidth, int imgHeight,
			     float*& dest, OSPVisItPatches* volume) 
{
    camera.SetScreen(xMin, xMax, yMin, yMax);
    renderer.SetModel(volume->GetWorld());
    renderer.SetCamera(camera.camera);
    volume->InitFB(imgWidth, imgHeight);
    volume->RenderFB();
}

*/
