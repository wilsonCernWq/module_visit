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

#include <cstdlib>

namespace ospray {
  namespace visit {
    
    // *************************************************//
    //                                                  //
    // OSPRay defines following environmental variables //
    //                                                  //
    // OSPRAY_DEBUG                                     //
    // OSPRAY_THREADS                                   //
    // OSPRAY_SET_AFFINITY                              //
    // OSPRAY_LOG_LEVEL                                 //
    // OSPRAY_LOG_OUTPUT                                //
    //                                                  //
    // We define one more environmental variable here   //
    //                                                  //
    // OSPRAY_VERBOSE                                   //
    //                                                  //
    // *************************************************//
    inline bool InitVerbose(const unsigned int level) 
    { 
      const char* env_verbose = std::getenv("OSPRAY_VERBOSE");
      if (env_verbose) {
        if (atoi(env_verbose) > 0) { return true; }
      }
      const char* env_debug = std::getenv("OSPRAY_DEBUG");
      if (env_debug) {
        if (atoi(env_debug) > 0) { return true; }
      }
      const char* env_log_level = std::getenv("OSPRAY_LOG_LEVEL");
      if (env_log_level) {
        if (atoi(env_log_level) > 0) { return true; }
      }
      return false;
    }
    /* There are 6 levels of debug output
     * 
     * level 0 : output in normal mode
     * level 5 : output everything related to ospray
     */
    inline bool CheckVerbose(const unsigned int level = 5)
    { 
      static bool verbose[6] = { 
        true,
        InitVerbose(1),
        InitVerbose(2),
        InitVerbose(3),
        InitVerbose(4),
        InitVerbose(5)
      };
      return verbose[level < 6 ? level : 5];
    }

  };
};

