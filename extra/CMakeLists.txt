## ======================================================================== ##
## Copyright 2009-2017 Intel Corporation                                    ##
##                                                                          ##
## Licensed under the Apache License, Version 2.0 (the "License");          ##
## you may not use this file except in compliance with the License.         ##
## You may obtain a copy of the License at                                  ##
##                                                                          ##
##     http://www.apache.org/licenses/LICENSE-2.0                           ##
##                                                                          ##
## Unless required by applicable law or agreed to in writing, software      ##
## distributed under the License is distributed on an "AS IS" BASIS,        ##
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. ##
## See the License for the specific language governing permissions and      ##
## limitations under the License.                                           ##
## ======================================================================== ##

INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}/ospray
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

OSPRAY_CREATE_LIBRARY(ospray_module_visit_common
  VisItImageComposite.ispc
  VisItImageComposite.cpp
  VisItImageComposite.h
  VisItWrapper.cpp
  VisItWrapper.h
  LINK
  ospray
  ospray_common
  ospray_module_ispc
  )

INSTALL(
  FILES
  VisItWrapperCore.h
  VisItWrapper.h
  VisItImageComposite.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ospray/visit
  )
