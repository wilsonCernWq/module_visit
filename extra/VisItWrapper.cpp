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
        throw std::runtime_error(s + "is invalid");
      }
    }

    struct Color { float R,G,B,A; };

    void TransferFunction::Set(const void *_table,
                             const unsigned int size, 
                             const double datamin, 
                             const double datamax) 
    {
      // initialize it once
      if (!init) 
      {
        ospray_rm(self);
        self = ospNewTransferFunction("piecewise_linear");
        ospray_check(self,"");
        init = true;
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
      ospSetData(self, "colors",    osp_cdata);
      ospSetData(self, "opacities", osp_odata);
      ospSet2f(self, "valueRange", datamin, datamax);
      ospCommit(self);
      // cleanup
      ospray_rm(osp_cdata);
      ospray_rm(osp_odata);
    }

  };
};
