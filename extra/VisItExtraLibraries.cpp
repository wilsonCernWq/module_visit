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

#include "VisItModuleCommon.h"
#include "VisItExtraLibraries.h"
#include <ospcommon/vec.h>
#include <cmath>
#include <vector>

using namespace ospcommon;

namespace ospray {
namespace visit {

  void ComputeGhostBounds(bool r[6], 
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
  
};
};
