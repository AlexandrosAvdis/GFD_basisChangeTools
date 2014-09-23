from math import sqrt, pi, sin, asin, cos, tan, atan2, acos, degrees, radians
import numpy as np

def cartesian_2_sphericalPolar(positionVectorCartesian):
    '''Convert Cartesian coordinates to spherical-polar coordinates, in radians.

The origin of the Cartesian frame of reference is located at the centre of the sphere, the positive half of x-axis goes through 0 deg E, 0 deg N, the positive half of y-axis goes through 90 deg E, 0 deg N and the positive half of the z-axis goes through the North Pole equivalent.'''
    #Extract the Cartesian vector components.
    x = positionVectorCartesian[0]
    y = positionVectorCartesian[1]
    z = positionVectorCartesian[2]
    #Calculate the radius.
    radius = sqrt(x**2 + y**2 + z**2)
    #Calculate azimuthal (theta) and zenith (phi) angles
    theta = acos(z/radius)
    phi = atan2(y,x)
    return [radius, theta, phi]

def sphericalPolar_2_cartesian(positionVectorSphericalPolar):
    '''Convert spherical-polar coordinates (in radians), into Cartesian coordinates.

The origin of the Cartesian frame of reference is located at the centre of the sphere, the positive half of x-axis goes through 0 deg E, 0 deg N, the positive half of y-axis goes through 90 deg E, 0 deg N and the positive half of the z-axis goes through the North Pole equivalent.'''
    import numpy as np
    [radius, theta, phi] = positionVectorSphericalPolar
    x = radius*np.sin(theta)*np.cos(phi)
    y = radius*np.sin(theta)*np.sin(phi)
    z = radius*np.cos(theta)
    return [x, y, z]

def cartesian_2_lonlatradius(positionVectorCartesian):
    '''Convert Cartesian coordinates on a sphere to longitude-latitude-radius. Longitude and latitude are returned in degrees.

The origin of the Cartesian frame of reference is located at the centre of the sphere, the positive half of x-axis goes through 0 deg E, 0 deg N, the positive half of y-axis goes through 90 deg E, 0 deg N and the positive half of the z-axis goes through the North Pole equivalent.'''
    #Calculate azimuthal (phi), polar (theta) angles and distance from origin.
    [radius, theta, phi] = cartesian_2_sphericalPolar(positionVectorCartesian)
    #Calculate longitude and latitude
    lon = degrees(phi)
    lat = degrees(pi/2 - theta)
    
    positionVectorLonlat = [lon, lat, radius]
    return positionVectorLonlat

def lonlatradius_2_cartesian(positionVectorLonLatRad):
    '''Convert longitude-latitude-radial coordinates on the surface of the Earth (in degrees) to Cartesian coordinates. Longitude and latitude must be in degrees.

The origin of the Cartesian frame of reference is located at the centre of the sphere, the positive half of x-axis goes through 0 deg E, 0 deg N, the positive half of y-axis goes through 90 deg E, 0 deg N and the positive half of the z-axis goes through the North Pole equivalent.'''
    #Calculate spherical-polar coordinates form longitude-latitude-radius.
    [radius, theta, phi] = lonlatradius_2_sphericalPolar(positionVectorLonLatRad)
    #Calculate Cartesian coordinates from spherical-polar coordinates.
    [x, y, z] = sphericalPolar_2_cartesian([radius, theta, phi])
    return [x, y, z]

def lonlatradius_2_sphericalPolar(positionVectorLonLatRad):
    '''Convert longitude-latitude-radial coordinates on the surface of the Earth (in degrees) to Cartesian coordinates. Longitude and latitude must be in degrees, the azimuthal and polar angles are returned in radians.'''
    [longitude, latitude, radius] = positionVectorLonLatRad
    #Calculate azimuthal (phi), polar (theta) angles.
    phi = np.radians(longitude)
    theta = pi/2. - np.radians(latitude)
    return [radius, theta, phi]

def cartesian_2_polarStereographic(positionVectorCartesian, surfaceRadius=6.37101e+06, southPoleCoordinates=None):
    '''Convert Cartesian coordinates on a sphere to polar stereographic, where the stereographic plane is tangent to the north pole. The output coordinates are also normalised with the sphere diameter. If the southPoleCoordinates argument is specified the input coordinates basis is rotated.'''
    x = np.double(positionVectorCartesian[0])
    y = np.double(positionVectorCartesian[1])
    z = np.double(positionVectorCartesian[2])
    #Calculate output vertical coordinate
    Z = (sqrt(x**2 + y**2 + z**2) - np.double(surfaceRadius))
    #If the South pole coordinates are specified, rotate the inputi
    # coordinate basis.
    if southPoleCoordinates != None:
        [x,y,z] = rotateCartesianBasis([x,y,z], southPoleCoordinates)
    #Convert Cartesian coordinates into unit-sphere Cartesian Coordinates
    x = x/(2.*np.double(surfaceRadius))
    y = y/(2.*np.double(surfaceRadius))
    z = z/(2.*np.double(surfaceRadius))
    #Convert unit-diameter-sphere cordinates into polar stereographic, carefull with
    # the South pole.
    if z <= -np.double(0.5):
        X = np.nan
        Y = np.nan
    else:
        X = (x/(z + 0.5))
        Y = (y/(z + 0.5))
    return np.double([X,Y,Z])

def rotateCartesianBasis(positionVectorCartesian, southPoleCoordinates):
    '''Function rotating the planet-centered Cartesian basis of the input vector. The rotation angles are specified vie entering the coordinates of the point that will end on the south pole after the rotation.'''
    southPoleLon = southPoleCoordinates[0]
    southPoleLat = southPoleCoordinates[1]
    angle1 = np.radians(-(90.+southPoleLon))
    angle2 = np.radians((90.+southPoleLat))
    #Rotate about the z-axis
    rotationMatrix1 = np.array([[ np.cos(angle1), np.sin(angle1),  0.0],
                                [-np.sin(angle1), np.cos(angle1),  0.0],
                                [0.0,                  0.0,                  1.0]
                               ])
    #Rotation about the x-axis
    rotationMatrix2 = np.array([[1.0,          0.0,                  0.0],
                                [0.0,          np.cos(angle2),  np.sin(angle2)],
                                [0.0,         -np.sin(angle2),  np.cos(angle2)]
                               ])
    rotationMatrix = np.dot(rotationMatrix2,rotationMatrix1)
    rotatedCoords = np.dot(rotationMatrix, np.array(positionVectorCartesian))
    return rotatedCoords.tolist()

def lonlatradius_2_polarStereographic(positionVectorLonLatRad, surfaceRadius=6.37101e+06, southPoleCoordinates=None):
    '''Convert longitude-latitude-radial coordinates on the surface of the Earth (in degrees) to polar stereographic, where the stereographic plane is tangent to the north pole. The output coordinates are also normalised with the sphere diameter.'''
    [x,y,z] = lonlatradius_2_cartesian(positionVectorLonLatRad)
    [X,Y,Z] = cartesian_2_polarStereographic([x,y,z], surfaceRadius, southPoleCoordinates)
    return [X,Y,Z]

def cartesian_2_unitDisk(positionVectorCartesian, surfaceRadius=6.37101e+06, southPoleCoordinates=None):
    '''Convert Cartesian coordinates on a sphere to a plane unit disk, where the plane is tangent to the north pole.'''
    #If the South pole coordinates are specified, rotate the input
    # coordinate basis.
    if southPoleCoordinates != None:
        [x,y,z] = rotateCartesianBasis(positionVectorCartesian, southPoleCoordinates)
    else:
        x = positionVectorCartesian[0]
        y = positionVectorCartesian[1]
        z = positionVectorCartesian[2]
    [lon,lat,rad] = cartesian_2_lonlatradius([x,y,z])
    [X,Y,Z] = lonlatradius_2_unitDisk([lon,lat,rad], surfaceRadius)
    return [X,Y,Z]

def lonlatradius_2_unitDisk(positionVectorLonLatRad, surfaceRadius=6.37101e+06, southPoleCoordinates=None):
    '''Convert longitude-latitude-radial coordinates on the surface of the Earth (in degrees) to a flat, unit disk. Similar to the North-polar Stereographic projection, but the "South" pole maps to points on the periphery of the unit circle.'''
    #If the South pole coordinates are specified, rotate the input
    # coordinate basis.
    if southPoleCoordinates != None:
        [x,y,z] = lonlatradius_2_cartesian(positionVectorLonLatRad)
        [x,y,z] =  rotateCartesianBasis([x,y,z], southPoleCoordinates)
        [longitude, latitude, radius] = cartesian_2_lonlatradius([x,y,z])
    else:
        longitude = positionVectorLonLatRad[0]
        latitude = positionVectorLonLatRad[1]
        radius = positionVectorLonLatRad[2]
    alpha = longitude
    r = (-1./180.)*latitude + 0.5
    X = r*np.cos(np.radians(alpha))
    Y = r*np.sin(np.radians(alpha))
    Z = surfaceRadius - radius
    return [X,Y,Z]

def transform_tensor_sphericalPolar_2_cartesian(positionVectorSpherical, tensor):
    '''Function changing the basis of a tensor from zonal-meridional-radial basis to a Cartesian basis.

The origin of the Cartesian frame of reference is located at the centre of the sphere, the positive half of x-axis goes through 0 deg E, 0 deg N, the positive half of y-axis goes through 90 deg E, 0 deg N and the positive half of the z-axis goes through the North Pole equivalent.'''
    import numpy as np
    from numpy import linalg
    #extract distance from origin, polar (theta) angles and azimuthal (phi) angles.
    [radius, theta, phi] = positionVectorSpherical
    #Evaluate components of rotation matrices.
    transformationMatrix =\
       np.array([\
                 [np.sin(theta)*np.cos(phi),   np.cos(theta)*np.cos(phi),  -np.sin(phi)],\
                 [np.sin(theta)*np.sin(phi),   np.cos(theta)*np.sin(phi),   np.cos(phi)],\
                 [np.cos(theta),              -np.sin(theta),              0]])
    transposedTransformationMatrix = transformationMatrix.transpose()
    #Calculate the components of the tensor in the reference system.
    transformed_Tensor = np.dot(transformationMatrix, np.array(tensor))
    transformed_Tensor = np.dot(transformed_Tensor,transposedTransformationMatrix)
    return transformed_Tensor

def transform_tensor_cartesian_2_sphericalPolar(positionVectorCartesian, tensor):
    '''Function transforming the components of a tensor from a Cartesian basis to a zonal-meridional-radial basis.

The origin of the Cartesian frame of reference is located at the centre of the sphere, the positive half of x-axis goes through 0 deg E, 0 deg N, the positive half of y-axis goes through 90 deg E, 0 deg N and the positive half of the z-axis goes through the North Pole equivalent.'''
    import numpy as np
    from numpy import linalg
    #Calculate azimuthal (theta) and zenith (phi) angles and distance from origin
    [radius, theta, phi] = cartesian_2_sphericalPolar(positionVectorCartesian)
    #Evaluate components of rotation matrices.
    transformationMatrix =\
       np.array([\
                 [np.sin(theta)*np.cos(phi),   np.sin(theta)*np.sin(phi),   np.cos(theta)],\
                 [np.cos(theta)*np.cos(phi),   np.cos(theta)*np.sin(phi),  -np.sin(theta)],\
                 [-np.sin(phi),                np.cos(phi),                 0]])
    transposedTransformationMatrix = transformationMatrix.transpose()
    #Calculate the components of the tensor in the reference system.
    transformed_Tensor = np.dot(transformationMatrix, np.array(tensor))
    transformed_Tensor = np.dot(transformed_Tensor,transposedTransformationMatrix)
    return transformed_Tensor

def transform_tensor_sphericalPolar_2_lon_lat_rad(tensor):
    '''Function transforming the components of a tensor from a spherical-polar basis to a zonal-meridional-radial basis.'''
    import numpy as np
    #The requested transformation is just a reflection followed by a change in the order
    # of the components in order to get a right-handed system.
    transformationMatrix = np.array([[ 0.0, 0.0, 1.0 ],\
                                     [ 0.0,-1.0, 0.0 ],\
                                     [ 1.0, 0.0, 0.0 ]])
    transformed_Tensor = np.dot(np.dot(transformationMatrix, np.array(tensor)), transformationMatrix)
    return transformed_Tensor

def transform_tensor_lon_lat_rad_2_sphericalPolar(tensor):
    '''Function transforming the components of a tensor from a zonal-meridional-radial basis to a spherical-polar basis.'''
    import numpy as np
    #The requested transformation is just a reflection followed by a change in the order
    # of the components in order to get a right-handed system.
    transformationMatrix = np.array([[ 0.0, 0.0, 1.0 ],\
                                     [ 0.0,-1.0, 0.0 ],\
                                     [ 1.0, 0.0, 0.0 ]])
    transformed_Tensor = np.dot(np.dot(transformationMatrix, np.array(tensor)), transformationMatrix)
    return transformed_Tensor

def transform_tensor_cartesian_2_lon_lat_rad(positionVectorCartesian, tensor):
    '''Function transforming the components of a tensor from a Cartesian basis to a zonal-meridional-radial basis.

The origin of the Cartesian frame of reference is located at the centre of the sphere, the positive half of x-axis goes through 0 deg E, 0 deg N, the positive half of y-axis goes through 90 deg E, 0 deg N and the positive half of the z-axis goes through the North Pole equivalent.'''
    import numpy as np
    #Transform from Cartesian into spherical-polar
    transformed_Tensor = transform_tensor_cartesian_2_sphericalPolar(positionVectorCartesian, tensor)
    #Transform from spherical-polar into longitude-latitude-radius.
    transformed_Tensor = transform_tensor_sphericalPolar_2_lon_lat_rad(transformed_Tensor)
    return transformed_Tensor

def transform_tensor_lon_lat_rad_2_cartesian(positionVectorLonLatRad, tensor):
    '''Function transforming the components of a tensor from a Cartesian basis to a zonal-meridional-radial basis.

The origin of the Cartesian frame of reference is located at the centre of the sphere, the positive half of x-axis goes through 0 deg E, 0 deg N, the positive half of y-axis goes through 90 deg E, 0 deg N and the positive half of the z-axis goes through the North Pole equivalent.'''
    import numpy as np
    #Transform coordinates from lon-lat-rad into spherical-polar.
    positionVectorSpericalPolar = lonlatradius_2_sphericalPolar(positionVectorLonLatRad)
    #Transform tensor from lon-lat-rad into spherical-polar.
    transformed_Tensor = transform_tensor_lon_lat_rad_2_sphericalPolar(tensor)
    #Transform spherical-polar into Cartesian.
    transformed_Tensor = transform_tensor_sphericalPolar_2_cartesian(positionVectorSpericalPolar, transformed_Tensor)
    return transformed_Tensor

def transform_vector_sphericalPolar_2_cartesian(positionVectorSpherical, vector):
    '''Function transforming the components of a vector from a spherical-polar basis to a Cartesian basis. The input position vector must be given as [radius, polar angle, azimuthal angle], all angles specified in radians.

The origin of the Cartesian frame of reference is located at the centre of the sphere, the positive half of x-axis goes through 0 deg E, 0 deg N, the positive half of y-axis goes through 90 deg E, 0 deg N and the positive half of the z-axis goes through the North Pole equivalent.'''
    import numpy as np
    #extract distance from origin, polar (theta) angles and azimuthal (phi) angles.
    [radius, theta, phi] = positionVectorSpherical
    #Evaluate components of rotation matrices.
    transformationMatrix =\
       np.array([\
                 [np.sin(theta)*np.cos(phi),   np.cos(theta)*np.cos(phi),  -np.sin(phi)],\
                 [np.sin(theta)*np.sin(phi),   np.cos(theta)*np.sin(phi),   np.cos(phi)],\
                 [np.cos(theta),              -np.sin(theta),              0]])
    #Calculate the components of the tensor in the reference system.
    transformed_Vector = np.dot(transformationMatrix, np.array(vector))
    return transformed_Vector

def transform_vector_cartesian_2_sphericalPolar(positionVectorCartesian, vector):
    '''Function transforming the components of a vector from a Cartesian basis to a spherical-polar basis.

The origin of the Cartesian frame of reference is located at the centre of the sphere, the positive half of x-axis goes through 0 deg E, 0 deg N, the positive half of y-axis goes through 90 deg E, 0 deg N and the positive half of the z-axis goes through the North Pole equivalent.'''
    import numpy as np
    #Calculate distance from origin, polar (theta) angles and azimuthal (phi) angles.
    [radius, theta, phi] = cartesian_2_sphericalPolar(positionVectorCartesian)
    #Evaluate components of rotation matrices.
    transformationMatrix =\
       np.array([\
                 [np.sin(theta)*np.cos(phi),   np.sin(theta)*np.sin(phi),   np.cos(theta)],\
                 [np.cos(theta)*np.cos(phi),   np.cos(theta)*np.sin(phi),  -np.sin(theta)],\
                 [-np.sin(phi),                np.cos(phi),                 0]])
    #Calculate the components of the tensor in the reference system.
    transformed_Vector = np.dot(transformationMatrix, np.array(vector))
    return transformed_Vector

def transform_vector_sphericalPolar_2_lon_lat_rad(vector):
    '''Function transforming the components of a vector from a spherical-polar basis to a zonal-meridional-radial basis.'''
    import numpy as np
    #The requested transformation is just a reflection followed by a change in the order
    # of the components in order to get a right-handed system.
    transformationMatrix = np.array([[ 0.0, 0.0, 1.0 ],\
                                     [ 0.0,-1.0, 0.0 ],\
                                     [ 1.0, 0.0, 0.0 ]])
    transformed_Vector = np.dot(transformationMatrix, np.array(vector))
    return transformed_Vector

def transform_vector_lon_lat_rad_2_sphericalPolar(vector):
    '''Function transforming the components of a vector from a zonal-meridional-radial basis to a spherical-polar basis.'''
    import numpy as np
    #The requested transformation is just a reflection followed by a change in the order
    # of the components in order to get a right-handed system.
    transformationMatrix = np.array([[ 0.0, 0.0, 1.0 ],\
                                     [ 0.0,-1.0, 0.0 ],\
                                     [ 1.0, 0.0, 0.0 ]])
    transformed_Vector = np.dot(transformationMatrix, np.array(vector))
    return transformed_Vector

def transform_vector_cartesian_2_lon_lat_rad(positionVectorCartesian, vector):
    '''Function transforming the components of a vector from a Cartesian basis to a zonal-meridional-radial basis.

The origin of the Cartesian frame of reference is located at the centre of the sphere, the positive half of x-axis goes through 0 deg E, 0 deg N, the positive half of y-axis goes through 90 deg E, 0 deg N and the positive half of the z-axis goes through the North Pole equivalent.'''
    #Calculate spherical-polar components of the vector.
    transformed_Vector = transform_vector_cartesian_2_sphericalPolar(positionVectorCartesian, vector)
    #Calculate zonal, meridional and radial components of the vector.
    transformed_Vector = transform_vector_sphericalPolar_2_lon_lat_rad(transformed_Vector)
    return transformed_Vector

def transform_vector_lon_lat_rad_2_cartesian(positionVectorRadLonLat, vector):
    '''Function transforming the components of a vector from a zonal-meridional-radial basis into a Cartesian basis.

The origin of the Cartesian frame of reference is located at the centre of the sphere, the positive half of x-axis goes through 0 deg E, 0 deg N, the positive half of y-axis goes through 90 deg E, 0 deg N and the positive half of the z-axis goes through the North Pole equivalent.'''
    #Transform coordinates from longitude-latitude-radius into spherical-polar
    positionVectorSpherical = lonlatradius_2_sphericalPolar(positionVectorRadLonLat)
    #Transform vector from longitude-latitude-radius into spherical-polar
    transformed_Vector = transform_vector_lon_lat_rad_2_sphericalPolar(vector)
    #Transform from spherical-polar into Cartesian
    transformed_Vector = transform_vector_sphericalPolar_2_cartesian(positionVectorSpherical, vector)
    return transformed_Vector
