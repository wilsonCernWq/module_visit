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

# this subdirectory builds the actual ospray plugin into a library
# "libospray_module_<modulename>.so".
#
# This shared library can depend on libospray, but should _not_ assume
# that it is living in the same address space as the application: When
# running in MPI mode ospray will, upon encountering the
# ospLoadModule("<moduleName>"), actually load and use this library on
# every worker node (which will typically not run the application). As
# such, the library should _not_ "commuincate" with the application
# through global variables or globally defind functions, and instead
# communicate solely through the ospNewGeometry(...), ospSet...(),
# ospCommit(...) etc functionality provided by ospray.


# build the actual shared library that defines the geometry.
#
# Note the name is important: In order for ospray to properly find and
# initialize a module referenced by a call to
# "ospLoadModule(<moduleName>) this module _has_ to
#
# a) be called libospray_module_<modulename>.so, and
# b) contain a (extern C linkage) initializatoin routine named
#    void ospray_init_module_<moduleName>()
#

OSPRAY_CREATE_LIBRARY(ospray_module_visit

  # volume without marcocell
  volume/VisItSharedStructuredVolume.cpp
  volume/VisItSharedStructuredVolume.ispc

  # and finally, the module init code (not doing much, but must be there)
  moduleInit.cpp
  
  # this depends on ospray core:
  LINK
  ospray
  ospray_common
  )

INSTALL(
  FILES 
  VisItModuleCommon.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ospray/visit
)

# INSTALL(
#   FILES 
#   include/cpp.h
#   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ospray/visit
# )

# INSTALL(
#   FILES
#   include/cpp/Camera.h
#   include/cpp/Data.h
#   include/cpp/Device.h
#   include/cpp/FrameBuffer.h
#   include/cpp/Geometry.h
#   include/cpp/Light.h
#   include/cpp/ManagedObject.h
#   include/cpp/Material.h
#   include/cpp/Model.h
#   include/cpp/PixelOp.h
#   include/cpp/Renderer.h
#   include/cpp/Texture2D.h
#   include/cpp/TransferFunction.h
#   include/cpp/Volume.h
#   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ospray/visit/cpp
# )
