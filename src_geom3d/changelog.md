Ref: https://ww2.mathworks.cn/matlabcentral/fileexchange/24484-geom3d?s_tid=srchtitle
change log for geom3d

geom3d release 2011.10.13
=========================

New Features
- added function vectorCross3d (thanks to Sven Holcombe)
- added function triangleArea3d
- added function trimeshSurfaceArea
- added drawPolygon3d (mostly the same as drawPolyline3d)

Several bug fixes and improvements in code speed, thanks to Sven Holcombe


geom3d release 2011.06.30
=========================

Important changes:
- Package has been splitted up into 'geom3d' and 'meshes3d'.
- Representation of geometrical shapes (3D circles, ellipsoids...) now uses
    degrees for angles

New features
- added function distanceLine3d
- added functions inertiaEllipsoid and drawEllipsoid
- added functions intersectLinePolygon3d, intersectRayPolygon3d, and
    intersectLineTriangle3d
- added function eulerAngleToRotation3d
- added drawing functions: drawTorus, drawCube, drawCuboid
- added function randomPointInBox3d
- added functions sph2cart2d and cart2sph2d, that work with degrees
    
Bug fixes
- fixed bugs in intersectEdgePlanes
- fixed bugs in isPerpendicular3d
- enhanced createPlane
 
 
geom3d release 2011.01.10
=========================

New features
- new functions for meshes and polyhedra:
    + meshEdgeLength
    + meshDihedralAngles
    + meshEdgeFaces
    + meshSurfaceArea
    + polyhedronMeanBreadth
    + checkMeshAdjacentFaces
    + drawFaceNormals
    + meshFace
    
- new function createBasisTransform3d to compute coordinate changes between
    basis
    
- added vectorAngle3d, to compute angle between 2 3D vectors
    
- added midPoint3d, to compute middle points of either 2 points or a 3D edge

- added functions for 3D rotations:
    + createRotation3dLineAngle
    + rotation3dAxisAndAngle
    + rotation3dToEulerAngles

Regressions
- function 'local2global3d' moved to private directory

Bug fixes
- most meshes now have their faces oriented with normals pointing outwards of
    the meshes
- fix bug in drawPlane for planes touching current axis at one corner
- fix bug in faceNormal for cell array of faces
- fix bug in transformPoint3d for large arrays

Code
- creation of private directory, for utilitary functions


geom3d release 2010.07.28
=========================

New features
- added createEulerAnglesRotation
- added functions for management of 3D boxes (isothetic cuboids)
- added drawPolyline3d
- added drawAxisCube function
- added recenterTransform3d, to change invariant point of 3D transforms.

Enhancements
- added support for several lines in clipLine3d
- added support for multiple inputs in distancePointPlane
- updated doc for polyhedra, for function localToGlobal3d
- updated several drawing functions

Bug fixes
- fixed several bugs in 3D rotations
- fixed bugs in drawCircle3d, drawCircleArc3d, drawEllipse3d.
- fixed bugs in drawCylinder


geom3d release 2009-06-17
=========================

* new features
- added transform functions:
    + transformLine3d, transformVector3d
- added intersectLineCylinder
- added 'localToGlobal3d': transform from local coordinates to global
    coordinates, using a translation vector and 3 rotation angles
    
* Changes
- renamed some functions to avoid name conflicts and ambiguities
    translation3d -> createTranslation3d
    rotationOx -> createRotationOx
    rotationOy -> createRotationOy
    rotationOz -> createRotationOz
    scale3d -> createScaling3d
    vecnorm3d -> vectorNorm3d
    normalize3d -> normalizeVector3d
    
* Compatibility considerations
- changed convention for angle in rotationOx, rotationOy and rotationOz.
- changed convention for dihedral angle, to be consistent with other references 
