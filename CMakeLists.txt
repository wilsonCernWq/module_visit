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

# Define an option to enable/disable this module.
#
# Though not actually required we guard each module 'moduleName' with
# a OSPRAY_MODULE_MODULENAME flag that enables/disables its building.
#
OPTION(OSPRAY_MODULE_VISIT "Build visit-ospray module" OFF)
IF (OSPRAY_MODULE_VISIT)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
  ADD_SUBDIRECTORY(ospray)
  ADD_SUBDIRECTORY(extra)
  INSTALL(
    FILES
    VisItModuleCommon.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ospray/visit)
ENDIF ()
