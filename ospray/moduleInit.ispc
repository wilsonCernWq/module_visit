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
				   const uniform int count,
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
  foreach(k = start ... (start + count))
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
