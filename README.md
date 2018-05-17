module_visit: Module for Data-Distributed Volume Rendering in VisIt
==========================================================================================================================

This module provides a new volume type for data	distributed rendering. The new volume provided is called `visit_shared_structured_volume`, which is very similar to the built-in `shared_structured_volume`. In particular it only takes three more additional parameters:

| Type | Name                         | Default  | Description                                         |
| ---- |:----------------------------:| --------:| ---------------------------------------------------:|
| bool | useGridAccelerator           | True     | Whether the empty space skipping grid will be built |
| vec3 | volumeGlobalBoundingBoxLower | disabled | Lower coordinate of the global volume               |
| vec3 | volumeGlobalBoundingBoxUpper | disabled | Upper coordinate of the global volume               |

The `useGridAccelerator` describes whether the empty-space-skipping grid will be build. The `volumeGlobalBoundingBoxLower` and `volumeGlobalBoundingBoxUpper` defines the global bounding box of the volume. 




