import numpy as np
import matplotlib.pyplot as pyplot
from enum import Enum

#epsilon which defines if a point is different from other
eps= 1e-10

#Types of boundary conditions
class bc(Enum):
    NONE = 0
    DIRICHLET = 1
    NEUMANN = 2



class Point:
    def __init__(self,x ,y):
        self.x = x
        self.y = y
        
    def __str__(self):
        return f"({self.x}, {self.y})"
    
        


class Node:
    value_bc = 0.0
    def __init__(self, x, y, id, bc):
        self.point = Point(x,y)
        self.id = id
        self.bc = bc

    def __str__(self):
        if self.bc == bc.NONE:
            return f"Coordinate: {self.point}\nid: {self.id}\nBoundary condition: {self.bc.name}"
        else:
            return f"Coordinate: {self.point}\nid: {self.id}\nBoundary condition: {self.bc.name}\nValue: {self.value_bc}"
    
    def get_point(self):
        return self.point
    
    def get_id(self):
        return self.id
    
    def get_bc(self):
        return self.bc, self.value_bc
    
    def set_bc(self, bc_type, value):
        self.bc = bc_type
        self.value_bc = value

    



#Square domain with lower_bound low and high_bound high.
class Domain:
    node_list = []
    num_nodes_x = 0
    num_nodes_y = 0
    num_nodes = 0

    def __init__(self, low, high, h):
        self.low = low
        self.high = high
        self.h = h

        self.num_nodes_x = int((self.high.x - self.low.x)/h) + 1
        self.num_nodes_y = int((self.high.y - self.low.y)/h) + 1
        self.num_nodes = self.num_nodes_x * self.num_nodes_y


    def __str__(self):
        return f"Low point: {self.low}\nHigh point: {self.high}\n h: {self.h}\n Number of points: {self.num_nodes}"
    
    def get_num_nodes_x(self):
        return self.num_nodes_x
    
    def get_num_nodes_y(self):
        return self.num_nodes_y
    
    def get_num_nodes(self):
        return self.num_nodes
    
    def get_node_by_id(self, id):
        return self.node_list[id]

    #creates a node list with natural ordering
    def create_node_list(self):
        k = 0
        for i in range(self.num_nodes_y):
            y = self.low.y + i*h

            for j in range(self.num_nodes_x):
                x = self.low.x + j*h
                self.node_list.append(Node(x,y, k, bc.NONE))
                k += 1

    #Creates a box with low point low and high point high. every point inside the box will have the boundary condition inputed. 
    # The f parameter is a function (of a point) that corresponds to that bc
    def add_bc(self, low, high, bc_type, f):
        low_eps = Point(low.x - eps, low.y - eps)
        high_eps = Point(high.x + eps, high.y + eps)
        for i, it in enumerate(self.node_list):
            p = it.get_point()
            if (low_eps.x < p.x < high_eps.x) and (low_eps.y < p.y < high_eps.y):
                    it.set_bc(bc_type, f(p))
                    


class Stencil:
    c_id = -1 #center
    t_id = -1 #top
    b_id = -1 #bottom
    l_id = -1 #left
    r_id = -1 #right
    def __init__(self, domain):
        self.domain = domain
        
    def center_node(self, center_id):
        self.c_id = center_id


    #should only be called when there is already a center point
    def assemble_stencil(self):
        c_node = self.domain.get_node_by_id(self.c_id)
        bc_type, _ = c_node.get_bc()
        if bc_type == bc.DIRICHLET:
            return
        else:
            nodes_x = self.domain.get_num_nodes_x()
            if self.c_id - nodes_x > 0:
                self.b_id = self.c_id - nodes_x()
            if self.c_id + nodes_x < self.domain.get_num_nodes():
                self.t_id = self.c_id + nodes_x
            if (self.c_id + 1)%nodes_x != 0:
                self.r_id = self.c_id + 1
            if (self.c_id + 1)%nodes_x != 1:
                self.l_id = self.c_id - 1

    def reset_stencil(self):
        self.c_id = -1 #center
        self.t_id = -1 #top
        self.b_id = -1 #bottom
        self.l_id = -1 #left
        self.r_id = -1 #right


h = 0.1
#bc of right and left
def f(point):
    return 0.0

#upper bc
def g(p):
    return 0.0

#bottom bc
def s(p):
    return np.e * np.sin(2 * np.pi * p.x)


d = Domain(Point(0.0, -1.0), Point(1.0, 1.0), h)
d.create_node_list()
d.add_bc(Point(0.0, -1.0), Point(1.0, -1.0), bc.NEUMANN, s)
d.add_bc(Point(0.0, 1.0), Point(1.0, 1.0), bc.NEUMANN, g)
d.add_bc(Point(0.0, -1.0), Point(0.0, 1.0), bc.DIRICHLET, f)
d.add_bc(Point(1.0, -1.0), Point(1.0, 1.0), bc.DIRICHLET, f)
for i in range(d.num_nodes):
    print(d.node_list[i])
