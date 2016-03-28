import MyLib_geometry as lib_geometry

def pbc_vector(vector,box):

    """ vector en numpy format y box en mylab format """

    if not box.__dna__:
        raise Exception('Box not in my format')

    return lib_geometry.pbc(vector,box.box,box.inv,box.ortho)


def distance_two_points(point1,point2,box=None):

    if box:

        if not box.__dna__:
            raise Exception('Box not in my format')

        return lib_geometry.distance_two_points(point1,point2,1,box.box,box.inv,box.ortho)

def angle_three_points(point1,point2,point3,box=None):

    if box:

        if not box.__dna__:
            raise Exception('Box not in my format')

        return lib_geometry.angle_three_points(point1,point2,point3,1,box.box,box.inv,box.ortho)


def angle_three_points_projxy(point1,point2,point3,box=None):

    if box:

        if not box.__dna__:
            raise Exception('Box not in my format')

        return lib_geometry.angle_three_points_projxy(point1,point2,point3,1,box.box,box.inv,box.ortho)


