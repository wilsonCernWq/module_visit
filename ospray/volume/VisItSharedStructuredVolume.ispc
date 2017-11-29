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

#include "volume/VisItSharedStructuredVolume.ih"

/*! get a voxel from given volume type */

#define template_getVoxel(type)                                              \
/* --------------------------------------------------------------------------\
// versions for pure 32-bit addressing. volume *MUST* be smaller than 2G     \
// ------------------------------------------------------------------------*/\
inline void SSV_getVoxel_##type##_32(void *uniform _self,                    \
                                     const varying vec3i &index,             \
                                     varying float &value)                   \
{                                                                            \
  /* Cast to the actual Volume subtype. */                                   \
  VisItSharedStructuredVolume *uniform self                                  \
      = (VisItSharedStructuredVolume *uniform)_self;                         \
                                                                             \
  /* Cast to the actual voxel type. */                                       \
  const type *uniform voxelData = (const type *uniform)self->voxelData;      \
  const uint32 addr = index.x +                                              \
      self->super.dimensions.x*(index.y + self->super.dimensions.y*index.z); \
                                                                             \
  /* The voxel value at the given index. */                                  \
  value = voxelData[addr];                                                   \
}                                                                            \
/* --------------------------------------------------------------------------\
// versions for 64/32-bit addressing. volume itself can be larger than       \
// 2G, but each slice must be within the 2G limit.                           \
// ------------------------------------------------------------------------*/\
inline void SSV_getVoxel_##type##_64_32(void *uniform _self,                 \
                                        const varying vec3i &index,          \
                                        varying float &value)                \
{                                                                            \
  /* Cast to the actual Volume subtype. */                                   \
  VisItSharedStructuredVolume *uniform self                                  \
      = (VisItSharedStructuredVolume *uniform)_self;                         \
                                                                             \
  const uniform uint8 *uniform basePtr =                                     \
      (const uniform uint8 *uniform)self->voxelData;                         \
                                                                             \
  /* iterate over slices, then do 32-bit gather in slice */                  \
  const uint32 ofs = index.x + self->super.dimensions.x * index.y;           \
  foreach_unique (z in index.z) {                                            \
    const uniform uint64 byteOffset = z * self->bytesPerSlice;               \
    const uniform type *uniform sliceData                                    \
      = (const uniform type *uniform )(basePtr + byteOffset);                \
    value = sliceData[ofs];                                                  \
  }                                                                          \
}                                                                            \
/* --------------------------------------------------------------------------\
// versions for full 64-bit addressing, for all dimensions or slice size     \
// ------------------------------------------------------------------------*/\
inline void SSV_getVoxel_##type##_64(void *uniform _self,                    \
                                     const varying vec3i &index,             \
                                     varying float &value)                   \
{                                                                            \
  /* Cast to the actual Volume subtype. */                                   \
  VisItSharedStructuredVolume *uniform self                                  \
      = (VisItSharedStructuredVolume *uniform)_self;                         \
                                                                             \
  const uint64 index64 = (uint64)index.x + self->super.dimensions.x *        \
      ((int64)index.y + self->super.dimensions.y * ((uint64)index.z));       \
  const uint32 hi28 = index64 >> 28;                                         \
  const uint32 lo28 = index64 & ((1<<28)-1);                                 \
                                                                             \
  foreach_unique (hi in hi28) {                                              \
    const uniform uint64 hi64 = hi;                                          \
    const type *uniform base = ((const type *)self->voxelData) + (hi64<<28); \
    value = base[lo28];                                                      \
  }                                                                          \
}

template_getVoxel(uint8);
template_getVoxel(int16);
template_getVoxel(uint16);
template_getVoxel(float);
template_getVoxel(double);
#undef template_getVoxel


/*! read a _typed_ value from an address that's given by an
    *BYTE*-offset relative to a base array. note that even though we
    assume that the offset is already in bytes (ie, WITHOUT scaling by
    the width of the array data type), the the base pointer must
    *STILL* be of the proper type for this macro to figure out the
    type of data it's reading from this address 
*/
#define template_accessArray(type)                                           \
inline float accessArrayWithOffset(const type *uniform basePtr,              \
                                   const varying uint32 offset)              \
{                                                                            \
  uniform uint8 *uniform base = (uniform uint8 *uniform)basePtr;             \
  return *((uniform type *)(base+offset));                                   \
}                                                                            \
inline float accessArrayWithOffset(const type *uniform basePtr,              \
                                   const uniform uint64 baseOfs,             \
                                   const varying uint32 offset)              \
{                                                                            \
  uniform uint8 *uniform base = (uniform uint8 *uniform)(basePtr);           \
  return *((uniform type *)((base+baseOfs)+offset));                         \
}                                                                            \
  
template_accessArray(uint8);
template_accessArray(int16);
template_accessArray(uint16);
template_accessArray(float);
template_accessArray(double);
#undef template_accessArray


/*! perform trilinear interpolation for given sample. unlike old way
  of doing this (a single computesample on the StructuredVolume level
  that calls the virtual 'getSample()' of the volume layout) this
  function will directly do all the addressing for the getSample
  (inlined), and thus be about 50% faster (wall-time, meaning even
  much faster in pure sample speed) 
*/
#define template_sample(type)                                                \
inline float SSV_sample_##type##_32(void *uniform _self,                     \
				    const vec3f &worldCoordinates)	     \
{                                                                            \
  /* Cast to the actual Volume subtype. */                                   \
  VisItSharedStructuredVolume *uniform self                                  \
      = (VisItSharedStructuredVolume *uniform)_self;                         \
                                                                             \
  /* Transform the sample location into the local coordinate system. */      \
  vec3f localCoordinates;                                                    \
  self->super.transformWorldToLocal(&self->super,                            \
                                    worldCoordinates,                        \
                                    localCoordinates);                       \
                                                                             \
  /* Coordinates outside the volume are clamped to the volume bounds. */     \
  const vec3f clampedLocalCoordinates                                        \
      = clamp(localCoordinates, make_vec3f(0.0f),                            \
              self->super.localCoordinatesUpperBound);                       \
                                                                             \
  /* Lower and upper corners of the box straddling the voxels to be          \
   * interpolated.                                                           \
   */									     \
  const vec3i voxelIndex_0 = to_int(clampedLocalCoordinates);                \
                                                                             \
  /* Fractional coordinates within the lower corner voxel used during        \
   * interpolation.                                                          \
   */									     \
  const vec3f frac = clampedLocalCoordinates - to_float(voxelIndex_0);       \
                                                                             \
  const uint32 voxelOfs                                                      \
    = voxelIndex_0.x * self->voxelOfs_dx                                     \
    + voxelIndex_0.y * self->voxelOfs_dy                                     \
    + voxelIndex_0.z * self->voxelOfs_dz;                                    \
  const type *uniform voxelData = (const type *uniform)self->voxelData;      \
  const uniform uint64 ofs000 = 0;                                           \
  const uniform uint64 ofs001 = self->bytesPerVoxel;                         \
  const float val000  = accessArrayWithOffset(voxelData,ofs000,voxelOfs);    \
  const float val001  = accessArrayWithOffset(voxelData,ofs001,voxelOfs);    \
  const float val00   = val000 + frac.x * (val001 - val000);                 \
                                                                             \
  const uniform uint64 ofs010 = self->bytesPerLine;                          \
  const uniform uint64 ofs011 = self->bytesPerLine+self->bytesPerVoxel;      \
  const float val010  = accessArrayWithOffset(voxelData,ofs010,voxelOfs);    \
  const float val011  = accessArrayWithOffset(voxelData,ofs011,voxelOfs);    \
  const float val01   = val010 + frac.x * (val011 - val010);                 \
                                                                             \
  const uniform uint64 ofs100 = self->bytesPerSlice;                         \
  const uniform uint64 ofs101 = ofs100 + ofs001;                             \
  const float val100  = accessArrayWithOffset(voxelData,ofs100,voxelOfs);    \
  const float val101  = accessArrayWithOffset(voxelData,ofs101,voxelOfs);    \
  const float val10   = val100 + frac.x * (val101 - val100);                 \
                                                                             \
  const uniform uint64 ofs110 = ofs100 + ofs010;                             \
  const uniform uint64 ofs111 = ofs100 + ofs011;                             \
  const float val110  = accessArrayWithOffset(voxelData,ofs110,voxelOfs);    \
  const float val111  = accessArrayWithOffset(voxelData,ofs111,voxelOfs);    \
  const float val11   = val110 + frac.x * (val111 - val110);                 \
                                                                             \
  /* Interpolate the voxel values. */                                        \
  const float val0    = val00  + frac.y * (val01  - val00);                  \
  const float val1    = val10  + frac.y * (val11  - val10);                  \
  const float val     = val0   + frac.z * (val1   - val0 );                  \
                                                                             \
  return val;                                                                \
}                                                                            \
                                                                             \
inline float SSV_sample_##type##_64_32(void *uniform _self,                  \
				       const vec3f &worldCoordinates)        \
{                                                                            \
  /* Cast to the actual Volume subtype. */                                   \
  VisItSharedStructuredVolume *uniform self                                  \
      = (VisItSharedStructuredVolume *uniform)_self;                         \
                                                                             \
  /* Transform the sample location into the local coordinate system. */      \
  vec3f localCoordinates;                                                    \
  self->super.transformWorldToLocal(&self->super,                            \
                                    worldCoordinates,                        \
                                    localCoordinates);                       \
                                                                             \
  /* Coordinates outside the volume are clamped to the volume bounds. */     \
  const vec3f clampedLocalCoordinates                                        \
    = clamp(localCoordinates, make_vec3f(0.0f),                              \
            self->super.localCoordinatesUpperBound);                         \
                                                                             \
  /* Lower and upper corners of the box straddling the voxels to be          \
   * interpolated.							     \
   */									     \
  const vec3i voxelIndex_0 = to_int(clampedLocalCoordinates);                \
                                                                             \
  /* Fractional coordinates within the lower corner voxel used during        \
   * interpolation.							     \
   */  									     \
  const vec3f frac = clampedLocalCoordinates - to_float(voxelIndex_0);       \
                                                                             \
  varying float ret = 0.f;                                                   \
  foreach_unique (sliceID in voxelIndex_0.z)  {                              \
    const uint32 voxelOfs                                                    \
      = voxelIndex_0.x * self->voxelOfs_dx                                   \
      + voxelIndex_0.y * self->voxelOfs_dy;                                  \
    const type *uniform voxelData                                            \
      = (const type *uniform)((uniform uint8*uniform)self->voxelData         \
                              +sliceID*self->bytesPerSlice);                 \
    const uniform uint64 ofs000 = 0;                                         \
    const uniform uint64 ofs001 = self->bytesPerVoxel;                       \
    const float val000  = accessArrayWithOffset(voxelData,ofs000,voxelOfs);  \
    const float val001  = accessArrayWithOffset(voxelData,ofs001,voxelOfs);  \
    const float val00   = val000 + frac.x * (val001 - val000);               \
                                                                             \
    const uniform uint64 ofs010 = self->bytesPerLine;                        \
    const uniform uint64 ofs011 = self->bytesPerLine+self->bytesPerVoxel;    \
    const float val010  = accessArrayWithOffset(voxelData,ofs010,voxelOfs);  \
    const float val011  = accessArrayWithOffset(voxelData,ofs011,voxelOfs);  \
    const float val01   = val010 + frac.x * (val011 - val010);               \
                                                                             \
    const uniform uint64 ofs100 = self->bytesPerSlice;                       \
    const uniform uint64 ofs101 = ofs100 + ofs001;                           \
    const float val100  = accessArrayWithOffset(voxelData,ofs100,voxelOfs);  \
    const float val101  = accessArrayWithOffset(voxelData,ofs101,voxelOfs);  \
    const float val10   = val100 + frac.x * (val101 - val100);               \
                                                                             \
    const uniform uint64 ofs110 = ofs100 + ofs010;                           \
    const uniform uint64 ofs111 = ofs100 + ofs011;                           \
    const float val110  = accessArrayWithOffset(voxelData,ofs110,voxelOfs);  \
    const float val111  = accessArrayWithOffset(voxelData,ofs111,voxelOfs);  \
    const float val11   = val110 + frac.x * (val111 - val110);               \
                                                                             \
    /* Interpolate the voxel values. */                                      \
    const float val0    = val00  + frac.y * (val01  - val00);                \
    const float val1    = val10  + frac.y * (val11  - val10);                \
    const float val     = val0   + frac.z * (val1   - val0 );                \
    ret = val;                                                               \
  }                                                                          \
  return ret;                                                                \
}

template_sample(uint8)
template_sample(int16)
template_sample(uint16)
template_sample(float)
template_sample(double)
#undef template_sample

// ray.time is set to interval length of intersected sample
inline void 
VisItSSV_stepRay(void *uniform _volume, 
		 varying Ray &ray, 
		 const varying float samplingRate)
{
  // Cast to the actual Volume subtype.
  uniform VisItSharedStructuredVolume *uniform volume =
    (uniform VisItSharedStructuredVolume *uniform) _volume;

  // The recommended step size for ray casting based volume renderers.
  const varying float step =
    volume->super.super.samplingStep * rcpf(samplingRate);
  const varying float rcp_step = rcpf(step);

  // Tentatively advance the ray.
  ray.t0 += step;
  ray.time = step;
  if (ray.t0 >= ray.t) return;

  // Intersect volume bounding box.
  varying float t0, t1;
  intersectBox(ray, volume->localBoundingBox, t0, t1);

  // Clip against volume clipping box (if specified).
  if (ne(volume->localClippingBox.lower,
  	 volume->localClippingBox.upper)) 
  {
    varying float tClip0, tClip1;
    intersectBox(ray, volume->localClippingBox, tClip0, tClip1);
    t0 = max(t0, tClip0);
    t1 = min(t1, tClip1);
  }

  // Advance ray if necessary
  if (ray.t0 <= t0) ray.t0 += max(0.f, ceil((t0    - ray.t0) * rcp_step)) * step;
  if (ray.t0 >= t1) ray.t0 += max(0.f, ceil((ray.t - ray.t0) * rcp_step)) * step;
}

inline void 
VisItSSV_intersectIsosurface(void *uniform _volume, 
			     uniform float *uniform isovalues, 
			     uniform int numIsovalues, 
			     varying Ray &ray)
{
  // Cast to the actual Volume subtype.
  uniform VisItSharedStructuredVolume *uniform volume =
    (uniform VisItSharedStructuredVolume *uniform) _volume;
  
  // The nominal step size for ray casting based volume renderers, not
  // considering the sampling rate.
  const uniform float step = volume->super.super.samplingStep;
  const uniform float rcp_step = rcpf(step);

  // Tentatively advance the ray.
  ray.t0 += step;
  if (ray.t0 >= ray.t) return;

  // Intersect volume bounding box.
  varying float t0, t1;
  intersectBox(ray, volume->localBoundingBox, t0, t1);

  // Clip against volume clipping box (if specified).
  if (ne(volume->localClippingBox.lower,
  	 volume->localClippingBox.upper)) 
  {
    varying float tClip0, tClip1;
    intersectBox(ray, volume->localClippingBox, tClip0, tClip1);
    t0 = max(t0, tClip0);
    t1 = min(t1, tClip1);
  }

  // Advance ray if necessary
  if (ray.t0 <= t0) ray.t0 += max(0.f, ceil((t0    - ray.t0) * rcp_step)) * step;
  if (ray.t0 >= t1) ray.t0 += max(0.f, ceil((ray.t - ray.t0) * rcp_step)) * step;
}

inline varying vec3f 
VisItSSV_computeGradient(void *uniform _volume, 
			 const varying vec3f &worldCoordinates)
{
  // Cast to the actual Volume subtype.
  StructuredVolume *uniform volume = (StructuredVolume *uniform) _volume;

  // Gradient step in each dimension (world coordinates).
  const uniform vec3f gradientStep = volume->gridSpacing * 0.25f;

  // The gradient will be computed using central differences.
  varying vec3f gradient;

  // Forward differences

  // Sample at gradient location.
  float sample = volume->super.sample(volume, worldCoordinates);

  // Gradient magnitude in the X direction.
  varying vec3f gradient_dx =
    worldCoordinates + make_vec3f(gradientStep.x, 0.0f, 0.0f);
  gradient.x = volume->super.sample(volume, gradient_dx) - sample;

  // Gradient magnitude in the Y direction.
  varying vec3f gradient_dy =
    worldCoordinates + make_vec3f(0.0f, gradientStep.y, 0.0f);
  gradient.y = volume->super.sample(volume, gradient_dy) - sample;

  // Gradient magnitude in the Z direction.
  varying vec3f gradient_dz =
    worldCoordinates + make_vec3f(0.0f, 0.0f, gradientStep.z);
  gradient.z = volume->super.sample(volume, gradient_dz) - sample;

  // This approximation may yield image artifacts.
  return(gradient / gradientStep);
}

inline void 
VisItSSV_setGlobalClippingBox
(uniform VisItSharedStructuredVolume *uniform self,
 const uniform box3f &globalClippingBox)
{
  // Set actual clipping box
  self->localClippingBox = self->super.super.volumeClippingBox;
  self->globalClippingBox = globalClippingBox;
  self->super.super.volumeClippingBox = globalClippingBox;
}

inline void 
VisItSSV_setGlobalBoundingBox
(uniform VisItSharedStructuredVolume *uniform self,
 const uniform box3f &globalBoundingBox)
{
  // Set actual bounding box
  if (ne(globalBoundingBox.lower,
	 globalBoundingBox.upper)) 
  {
    self->localBoundingBox = self->super.super.boundingBox;
    self->globalBoundingBox = globalBoundingBox;
    self->super.super.boundingBox = globalBoundingBox;
  }
}

inline void 
VisItSSV_selectGridAccelerator
(uniform VisItSharedStructuredVolume *uniform self,
 const uniform int useGridAccelerator)
{
  // Register
  self->useGridAccelerator = useGridAccelerator;
  if (useGridAccelerator == 0) {
    // Don't use grid acceleration
    self->super.super.stepRay = self->na_stepRay;
    self->super.super.intersectIsosurface = self->na_intersectIsosurface;
  }
  else {
    // Use grid acceleration
    self->super.super.stepRay = self->ga_stepRay;
    self->super.super.intersectIsosurface = self->ga_intersectIsosurface;
  }
}

void VisItSharedStructuredVolume_Constructor
(uniform VisItSharedStructuredVolume *uniform self,
 void *uniform cppEquivalent,
 const uniform int voxelType,
 const uniform vec3i &dimensions,
 const void *uniform voxelData)
{
  // Call parent constructor
  StructuredVolume_Constructor(&self->super, cppEquivalent, dimensions);
  
  // Backup grid accelerator functions
  self->ga_stepRay                  = self->super.super.stepRay;
  self->ga_intersectIsosurface      = self->super.super.intersectIsosurface;
  self->na_stepRay                  = VisItSSV_stepRay;
  self->na_intersectIsosurface      = VisItSSV_intersectIsosurface;
  self->super.super.computeGradient = VisItSSV_computeGradient;

  // Calculate volume size
  uniform uint64 bytesPerVoxel;

  if (voxelType == OSP_UCHAR)
    bytesPerVoxel = sizeof(uniform uint8);
  else if (voxelType == OSP_SHORT)
    bytesPerVoxel = sizeof(uniform int16);
  else if (voxelType == OSP_USHORT)
    bytesPerVoxel = sizeof(uniform uint16);
  else if (voxelType == OSP_FLOAT)
    bytesPerVoxel = sizeof(uniform float);
  else if (voxelType == OSP_DOUBLE)
    bytesPerVoxel = sizeof(uniform double);
  else {
    print("#osp:shared_structured_volume: unknown voxel type\n");
    return;
  }

  const uniform uint64 bytesPerLine = bytesPerVoxel * dimensions.x;
  const uniform uint64 bytesPerSlice = bytesPerLine * dimensions.y;
  const uniform uint64 bytesPerVolume = bytesPerSlice * dimensions.z;

  self->voxelType     = (OSPDataType) voxelType;
  self->voxelData     = voxelData;
  self->bytesPerVoxel = bytesPerVoxel;
  self->bytesPerLine  = bytesPerLine;
  self->bytesPerSlice = bytesPerSlice;
  self->voxelOfs_dx   = bytesPerVoxel;
  self->voxelOfs_dy   = bytesPerLine;
  self->voxelOfs_dz   = bytesPerSlice;

  if (bytesPerVolume <= (1ULL<<30)) {
    //print("#osp:shared_structured_volume: using 32-bit mode\n");
    // in this case, we know ALL addressing can be 32-bit.

    if (voxelType == OSP_UCHAR) {
      self->super.getVoxel = SSV_getVoxel_uint8_32;
      self->super.super.sample = SSV_sample_uint8_32;
    } else if (voxelType == OSP_SHORT) {
      self->super.getVoxel = SSV_getVoxel_int16_32;
      self->super.super.sample = SSV_sample_float_32;
    } else if (voxelType == OSP_USHORT) {
      self->super.getVoxel = SSV_getVoxel_uint16_32;
      self->super.super.sample = SSV_sample_float_32;
    } else if (voxelType == OSP_FLOAT) {
      self->super.getVoxel = SSV_getVoxel_float_32;
      self->super.super.sample = SSV_sample_float_32;
    } else if (voxelType == OSP_DOUBLE) {
      self->super.getVoxel = SSV_getVoxel_double_32;
      self->super.super.sample = SSV_sample_double_32;
    }

  } else if (bytesPerSlice <= (1ULL << 30)) {
    //print("#osp:shared_structured_volume: using 64/32-bit mode\n");
    // in this case, we know we can do 32-bit addressing within a
    // slice, but need 64-bit arithmetic to get slice begins

    if (voxelType == OSP_UCHAR) {
      self->super.getVoxel = SSV_getVoxel_uint8_64_32;
      self->super.super.sample = SSV_sample_uint8_64_32;
    } else if (voxelType == OSP_SHORT) {
      self->super.getVoxel = SSV_getVoxel_int16_64_32;
      self->super.super.sample = SSV_sample_int16_64_32;
    } else if (voxelType == OSP_USHORT) {
      self->super.getVoxel = SSV_getVoxel_uint16_64_32;
      self->super.super.sample = SSV_sample_uint16_64_32;
    } else if (voxelType == OSP_FLOAT) {
      self->super.getVoxel = SSV_getVoxel_float_64_32;
      self->super.super.sample = SSV_sample_float_64_32;
    } else if (voxelType == OSP_DOUBLE) {
      self->super.getVoxel = SSV_getVoxel_double_64_32;
      self->super.super.sample = SSV_sample_double_64_32;
    }
  } else {
    //print("#osp:shared_structured_volume: using 64-bit mode\n");
    // in this case, even a single slice is too big to do 32-bit
    // addressing, and we have to do 64-bit throughout

    if (voxelType == OSP_UCHAR)
      self->super.getVoxel = SSV_getVoxel_uint8_64;
    else if (voxelType == OSP_SHORT)
      self->super.getVoxel = SSV_getVoxel_int16_64;
    else if (voxelType == OSP_USHORT)
      self->super.getVoxel = SSV_getVoxel_uint16_64;
    else if (voxelType == OSP_FLOAT)
      self->super.getVoxel = SSV_getVoxel_float_64;
    else if (voxelType == OSP_DOUBLE)
      self->super.getVoxel = SSV_getVoxel_double_64;
  }
}

export void
VisItSSV_updateBoundingBox(void *uniform _self,			   
			   const uniform box3f &globalBoundingBox)
{
  uniform VisItSharedStructuredVolume *uniform self =
    (uniform VisItSharedStructuredVolume *uniform) _self;
  VisItSSV_setGlobalBoundingBox(self, globalBoundingBox);
}

export void
VisItSSV_updateClippingBox(void *uniform _self,			   
			   const uniform box3f &globalClippingBox)
{
  uniform VisItSharedStructuredVolume *uniform self =
    (uniform VisItSharedStructuredVolume *uniform) _self;
  VisItSSV_setGlobalClippingBox(self, globalClippingBox);
}

export void
VisItSSV_updateGridAccelerator(void *uniform _self, 
			       const uniform int useGridAccelerator)
{
  uniform VisItSharedStructuredVolume *uniform self =
    (uniform VisItSharedStructuredVolume *uniform) _self;
  VisItSSV_selectGridAccelerator(self, useGridAccelerator);
}

export void *uniform 
VisItSSV_createInstance(void *uniform cppEquivalent, 
			const uniform int voxelType,
			const uniform vec3i &dimensions, 
			const void *uniform voxelData)
{
  // Allocate the volume container
  VisItSharedStructuredVolume *uniform volume = 
    uniform new uniform VisItSharedStructuredVolume;
  
  // Create Volume
  VisItSharedStructuredVolume_Constructor(volume,
					  cppEquivalent,
					  voxelType,
					  dimensions,
					  voxelData);
  return volume;
}
