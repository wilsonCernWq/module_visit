// ======================================================================== //
// Copyright 2016-2018 Qi WU                                                //
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

#include "common/OSPCommon.ih"

export void ISPC_Experiment()
{
  print("[ospray] experiment_visit from ISPC\n");
}

export void ISPC_ComposeBackground(const uniform int start,
				   const uniform int end,
				   const uniform int *uniform fgRegion,
				   const uniform int fgW,
				   const uniform int fgH,
				   const uniform float *uniform fgBuffer,
				   const uniform int bgW,
				   const uniform int bgH,
				   const uniform uint8 *uniform bgColor,
				   const uniform float *uniform         bgDepth,
				   uniform uint8 *uniform output)
{
  foreach(k = start ... end)
  {
    int x = k % bgW;
    int y = k / bgW;
    int bgIdx = y * bgW + x;
    int fgIdx = (y - fgRegion[2]) * fgW + (x - fgRegion[0]);
    if ((x >= fgRegion[0] && x < fgRegion[1]) &&
	(y >= fgRegion[2] && y < fgRegion[3]))
    {
      if (fgBuffer[fgIdx * 4 + 3] == 0)
      {
	output[bgIdx * 3 + 0] = bgColor[bgIdx * 3 + 0];
	output[bgIdx * 3 + 1] = bgColor[bgIdx * 3 + 1];
	output[bgIdx * 3 + 2] = bgColor[bgIdx * 3 + 2];
      }
      else
      {
	float alpha = (1.0 - fgBuffer[fgIdx * 4 + 3]);
	output[bgIdx * 3 + 0] =
	  clamp(bgColor[bgIdx * 3 + 0] * alpha + fgBuffer[fgIdx * 4 + 0] * 255.f, 0.f, 255.f);
	output[bgIdx * 3 + 1] =
	  clamp(bgColor[bgIdx * 3 + 1] * alpha + fgBuffer[fgIdx * 4 + 1] * 255.f, 0.f, 255.f);
	output[bgIdx * 3 + 2] =
	  clamp(bgColor[bgIdx * 3 + 2] * alpha + fgBuffer[fgIdx * 4 + 2] * 255.f, 0.f, 255.f);
      }
    }
    else
    {
      output[bgIdx * 3 + 0] = bgColor[bgIdx * 3 + 0];
      output[bgIdx * 3 + 1] = bgColor[bgIdx * 3 + 1];
      output[bgIdx * 3 + 2] = bgColor[bgIdx * 3 + 2];
    }
  }
}

export void ISPC_BlendFrontToBack(const uniform int start,
				  const uniform int end,
				  const uniform int startX,
				  const uniform int startY,
				  const uniform int W, const uniform int H,
				  const uniform int srcX,
				  const uniform int srcY,
				  const uniform int dstX,
				  const uniform int dstY,
				  const uniform int   * uniform blendExtents,
				  const uniform int   * uniform srcExtents,
				  const uniform float * uniform srcImage,
				  const uniform int   * uniform dstExtents,
				  uniform float       * uniform dstImage)
{  
  foreach(k = start ... end)
  {
    // get indices
    int x = k % W + startX;
    int y = k / W + startY;
    int srcIndex = (srcX * (y-srcExtents[2]) + x-srcExtents[0]) * 4;
    int dstIndex = (dstX * (y-dstExtents[2]) + x-dstExtents[0]) * 4;    
    // front to back compositing
    if (dstImage[dstIndex + 3] < 1.0f) {
      float trans = 1.0f - dstImage[dstIndex + 3];
      dstImage[dstIndex+0] = 
        clamp(srcImage[srcIndex+0] * trans + dstImage[dstIndex+0],
	      0.0f, 1.0f);
      dstImage[dstIndex+1] = 
	clamp(srcImage[srcIndex+1] * trans + dstImage[dstIndex+1],
	      0.0f, 1.0f);
      dstImage[dstIndex+2] = 
	clamp(srcImage[srcIndex+2] * trans + dstImage[dstIndex+2],
	      0.0f, 1.0f);
      dstImage[dstIndex+3] = 
	clamp(srcImage[srcIndex+3] * trans + dstImage[dstIndex+3],
	      0.0f, 1.0f);
    }
  }
}

export void ISPC_BlendBackToFront(const uniform int start,
				  const uniform int end,
				  const uniform int startX,
				  const uniform int startY,
				  const uniform int W, const uniform int H,
				  const uniform int srcX,
				  const uniform int srcY,
				  const uniform int dstX,
				  const uniform int dstY,
				  const uniform int   * uniform blendExtents,
				  const uniform int   * uniform srcExtents,
				  const uniform float * uniform srcImage,
				  const uniform int   * uniform dstExtents,
				  uniform float       * uniform dstImage)
{  
  foreach(k = start ... end)
  {
    // get indices
    int x = k % W + startX;
    int y = k / W + startY;
    int srcIndex = (srcX * (y-srcExtents[2]) + x-srcExtents[0]) * 4;
    int dstIndex = (dstX * (y-dstExtents[2]) + x-dstExtents[0]) * 4;    
    // back to front compositing	    
    float trans = 1.0f - srcImage[srcIndex + 3];
    dstImage[dstIndex+0] = 
      clamp(dstImage[dstIndex+0] * trans + srcImage[srcIndex+0],
	    0.0f, 1.0f);
    dstImage[dstIndex+1] = 
      clamp(dstImage[dstIndex+1] * trans + srcImage[srcIndex+1],
	    0.0f, 1.0f);
    dstImage[dstIndex+2] = 
      clamp(dstImage[dstIndex+2] * trans + srcImage[srcIndex+2],
	    0.0f, 1.0f);
    dstImage[dstIndex+3] = 
      clamp(dstImage[dstIndex+3] * trans + srcImage[srcIndex+3],
	    0.0f, 1.0f);
  }
}

