from metrics import Metrics
import matplotlib.pyplot as plt
from numpy import meshgrid

class ObjectError(Exception): pass

class BaseType(type):

    def __repr__(cls):
        return cls.__name__

class Point(metaclass=BaseType):


    """
        Basic implementation of a Point using 3 variables x, y, z

        Point is the baseclass for all the objects defined next
        i.e. Vector, Line and Plane

        initialize using Point(x, y, z) -> returns a Point object with
                                           spatial position x, y, z

    """


    def __init__(self, x, y, z):

        self.x= x
        self.y= y
        self.z= z
        self.xyz= x, y, z
        self._xyz= round(x, 3), round(y, 3), round(z, 3)

    @staticmethod
    def origin():

        return Point(0, 0, 0)

    def P2Vec(self):

        return Vector(*self.xyz)

    def get_distance(self, other):

        __x= (self.x - other.x)**2
        __y= (self.y - other.y)**2
        __z= (self.z - other.z)**2

        return (__x + __y + __z)**0.5

    def __eq__(self, other):

        return (self.x==other.x and self.y==other.y and self.z==other.z)

    def __neg__(self):

        return Point(-self.x, -self.y, -self.z)

    def __sub__(self, other):

        X, Y, Z= self.x-other.x, self.y-other.y, self.z-other.z
        return Vector(X, Y, Z)

    def __repr__(self):

        return "Point({}, {}, {})".format(round(self.x, 3), round(self.y, 3), round(self.z, 3))


class Vector(Point, Metrics, metaclass=BaseType):


    """
        Basic implementation of a vector using 3 variables x, y, z
        interactions available with another vector.

        Vector is abstracted from another class Point and adds directional geometry
        to the Point class

        initialize using Vector(x, y, z) -> returns a Vector object

    """

    def __init__(self, x_coord, y_coord, z_coord):

        Point.__init__(self, x_coord, y_coord, z_coord)
        self.length= self.get_distance(Point.origin())
        self.unit_vector= self.get_unit_vector()
        self.str_d= "Vector({}, {}, {})".format(*self.xyz)

    def get_unit_vector(self):

        try:
            if round(self.length) == 1:
                return self
            else:
                unit_x= round(self.x/self.length, 3)
                unit_y= round(self.y/self.length, 3)
                unit_z= round(self.z/self.length, 3)
                return Vector(unit_x, unit_y, unit_z)

        except ZeroDivisionError:
            return None

    def equal_dir(self, other):

        svec, ovec= self.unit_vector, other.unit_vector
        return (svec.x==ovec.x and svec.y==ovec.y and svec.z==ovec.z)

    def __neg__(self):

        return Vector(-self.x, -self.y, -self.z)

    def __add__(self, other):

        X, Y, Z= self.x+other.x, self.y+other.y, self.z+other.z
        return Vector(X, Y, Z)

    def __sub__(self, other):

        X, Y, Z= self.x-other.x, self.y-other.y, self.z-other.z
        return Vector(X, Y, Z)

    def __mul__(self, other):

        try:
            return Vector(self.x*other, self.y*other, self.z*other)
        except:
            return Vector(self.x*other.x, self.y*other.y, self.z*other.z)

    def __rmul__(self, other):

        try:
            return Vector(self.x*other, self.y*other, self.z*other)
        except:
            return Vector(self.x*other.x, self.y*other.y, self.z*other.z)

    def __repr__(self):

        return "Vector({}, {}, {}, d={})".format(self.x, self.y, self.z, round(self.length, 3))

    def Vec2P(self):

        return Point(*self.xyz)

    @staticmethod
    def dot(*vectors):

        x=y=z=1
        for v in vectors:
            x*= v.x
            y*= v.y
            z*= v.z

        return x+y+z

    @staticmethod
    def cross(vec1, vec2):

        X= (vec1.y*vec2.z)-(vec1.z*vec2.y)
        Y= (vec1.z*vec2.x)-(vec1.x*vec2.z)
        Z= (vec1.x*vec2.y)-(vec1.y*vec2.x)

        try:
            return Vector(X, Y, Z)
        except ZeroDivisionError:
            return Point(X, Y, Z)

class Line(metaclass=BaseType):


    """
        Basic implementation of a line using 2 points or a point-vector pair.
        interactions available with a line and a point.

        initialize using Line(Point1, Vector) -> returns a line object

    """


    lamb= "\u03BB"

    def __init__(self, point, vector):

        self.ref= point
        self.dir= vector

    @staticmethod
    def from_points(p1, p2):

        v1= Vector(p1.x, p1.y, p1.z)
        dir= Vector(p2.x, p2.y, p2.z) - v1
        return Line(p1, dir)

    def get_point(self, k=0):

        r= self.ref
        p= Vector(r.x, r.y, r.z) + (k*self.dir)
        return Point(p.x, p.y, p.z)

    def interact_point(self, point):

        res= {}

        dot_dir_p= Vector.dot(self.dir, point)
        dot_dir_l= Vector.dot(self.dir, self.ref)
        dot_dir_dir= Vector.dot(self.dir, self.dir)
        l= (dot_dir_p - dot_dir_l)/dot_dir_dir
        FOP= self.get_point(l)
        distance= round(FOP.get_distance(point), 3)
        res["DISTANCE"]= distance
        res["FOOT"]= FOP
        if distance == 0:
            image= FOP

        else:
            image= ((2*FOP.P2Vec())-point.P2Vec()).Vec2P()
            del_p= point.P2Vec()-self.ref.P2Vec()
            N= Vector.cross(del_p, self.dir)
            d= Vector.dot(-self.ref.P2Vec(), N)
            plane= Plane(*N.xyz, d)
            res["PLANE"]= plane._xyzd

        res["IMAGE"]= image

        return res

    def interact_line(self, other):

        res= {}

        cos_theta= Vector.dot(self.dir, other.dir)/(self.dir.length*other.dir.length)
        res["cos_theta"]= round(cos_theta, 3)

        x1, y1, z1= self.ref.xyz
        a1, b1, c1= self.dir.xyz
        x2, y2, z2= other.ref.xyz
        a2, b2, c2= other.dir.xyz

        set1= (a2, -a1, x2-x1)
        set2= (b2, -b1, y2-y1)
        set3= (c2, -c1, z2-z1)

        try:
            mu, la= Metrics.solve2D(set1, set2)
            POI1= self.get_point(la)
            POI2= other.get_point(mu)
            if POI1==POI2:
                res["POI"]= POI1
            else:
                cross_dir= Vector.cross(self.dir, other.dir)
                nume= Vector.dot(other.ref.P2Vec()-self.ref.P2Vec(), cross_dir)
                deno= cross_dir.length
                res["SKEW"]= True
                res["MIN DISTANCE"]= round(abs(nume/deno), 3)

        except ZeroDivisionError:
            if self.ref==other.ref:
                raise ObjectError("Similar Objects found !!!")
            else:
                res["PARALLEL"]= True
                res["MIN DISTANCE"]= self.interact_point(other.ref)["DISTANCE"]

        return res

    def __repr__(self):

        return "Line({} + {}{})".format(self.ref._xyz, Line.lamb, self.dir._xyz)

class Plane(Vector, Point, metaclass=BaseType):

    """
        Basic implementation of a plane using 4 variables x, y, z, d
        interactions available with a plane, line and a point.

        initialize using Plane(x, y, z, d) -> returns a plane object

    """

    def __init__(self, normal_x, normal_y, normal_z, d):

        Vector.__init__(self, normal_x, normal_y, normal_z)
        self.d= d
        self.xyzd= *self.xyz, self.d
        self._xyzd= round(self.x, 3), round(self.y, 3), round(self.z, 3), round(self.d, 3)

    @staticmethod
    def from_xyzd(x, y, z, d):

        return Plane(x, y, z, d)

    @staticmethod
    def from_points(p1, p2, p3):

        direction= Vector.cross(p3-p1, p2-p1)
        d= Vector.dot(direction, p1) or direction.dot(p2)

        return Plane.from_xyzd(*direction.xyz, -d)

    def interact_point(self, point):

        res= {}

        D= Vector.dot(self.unit_vector, point)+(self.d/self.length)

        if D != 0:
            __x= self.x *(D/self.length)
            __y= self.y *(D/self.length)
            __z= self.z *(D/self.length)

            Foot_x= point.x- __x
            Foot_y= point.y- __y
            Foot_z= point.z- __z
            Foot= Point(Foot_x, Foot_y, Foot_z)

            Image_x= point.x- (2*__x)
            Image_y= point.y- (2*__y)
            Image_z= point.z- (2*__z)
            Image= Point(Image_x, Image_y, Image_z)

        else:
            Foot= Image= point

        res["DISTANCE"]= abs(round(D, 3))
        res["FOOT"]= Foot
        res["IMAGE"]= Image

        return res

    def interact_plane(self, other):

        res= {}

        if (self.unit_vector==other.unit_vector):
            if self.xyz == other.xyz:
                d= other.d
            else:
                l= self.x/other.x
                d= other.d*l

            D= abs(round((self.d-d)/self.length, 3))

        else:
            D= 0
            line_dir= Vector.cross(self, other)
            try:
                c1= self.x, self.y, self.d
                c2= other.x, other.y, other.d
                x, y= Metrics.solve2D(c1, c2)
                line_ref= Point(x, y, 0)

            except ZeroDivisionError:
                try:
                    c1= self.y, self.z, self.d
                    c2= other.y, other.z, other.d
                    y, z= Metrics.solve2D(c1, c2)
                    line_ref= Point(0, y, z)

                except ZeroDivisionError:
                    c1= self.x, self.z, self.d
                    c2= other.x, other.z, other.d
                    x, z= Metrics.solve2D(c1, c2)
                    line_ref= Point(x, 0, z)

            line= Line(line_ref, line_dir)
            res["LOI"]= line

        res["DISTANCE"]= D
        return res

    def interact_line(self, line):

        res= {}

        dot_ref= Vector.dot(self, line.ref)
        dot_dir= Vector.dot(self, line.dir)
        if dot_dir != 0:
            l= (self.d-dot_ref)/dot_dir
            res["INTERSECTION"]= line.get_point(l)

        else:
            deno= self.length
            nume= Vector.dot(self, line.ref) -self.d
            d= abs(nume/deno)
            res["DISTANCE"]= d

        return res

    def __repr__(self):

        return "PLANE({}, {}, {}, {})".format(self.x, self.y, self.z, self.d)

    def render(self, mesh_range=5):

        fig= plt.figure()
        ax= fig.add_subplot(111, projection='3d')

        xx, yy= meshgrid(range(mesh_range), range(mesh_range))
        z= (self.d-(self.x*xx)-(self.y*yy))/self.z

        ax.plot_surface(xx, yy, z, alpha=1)
        plt.scatter(0, 0, 0, label="origin")
        plt.show()



