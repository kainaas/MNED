import numpy as np
import matplotlib.pyplot as pyplot
from enum import Enum

###################################################
#NOTA 1 Essa implementação talvez tenha ficado grande demais pra esse proble e a implementação do ex_6 na pasta lista_3 esteja mais clean
#Isso aqui é algo mais geral, mas também mais ilegível. Deixo a vc a escolha do código

#NOTA 2 eu nao tinha visto que tinha uma condição de Robin e eu acho q vai ser um inferno implementar

#NOTA 3 acho melhor apenas modificar valores do ex_6.py

#NOTA 4 eu nao testei essa versão
##################################################

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
        return self.bc
    
    def get_bc_value(self):
        return self.value_bc
    
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

    def get_h(self):
        return self.h

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
                    


class Stencil: #Specific for this problem
    ids = -1*np.ones(5) #by order: center, top, bottom, left, right
    weights_general = np.zeros(5) #weights of interior points (should not change)
    weights = np.zeros(5) #changes if the node is bc or interior
    rhs = 0.0
    def __init__(self, domain):
        self.domain = domain
        self.h = self.domain.get_h()
        
    def center_node(self, center_id):
        self.ids[0] = center_id


    def set_weight_center(self, w):
        self.weights_general[0] = w

    def set_weight_top(self, w):
        self.weights_general[1] = w

    def set_weight_bottom(self, w):
        self.weights_general[2] = w

    def set_weight_left(self, w):
        self.weights_general[3] = w

    def set_weight_right(self, w):
        self.weights_general[4] = w


    #should only be called when there is already a center point
    def assemble_stencil(self):
        c_node = self.domain.get_node_by_id(self.ids[0])
        bc_type = c_node.get_bc()
        bc_value = c_node.get_bc_value()
        
        if bc_type == bc.DIRICHLET:
            self.rhs = bc_value
            self.weights[0] = 1

        elif bc_type == bc.NONE:
            nodes_x = self.domain.get_num_nodes_x()

            if self.ids[0] + nodes_x < nodes_x: #defines the top id
                self.ids[1] = self.ids[0] + nodes_x
                if self.ids[1] == bc.DIRICHLET:
                    self.rhs -= self.weights_general[1] * self.ids[1].get_bc_value()
                else:
                    self.weights[1]  = self.weights_general[1]
            
            if self.ids[0] - nodes_x > 0: #defines the bottom id
                self.ids[2] = self.ids[0] - nodes_x()
                if self.ids[2] == bc.DIRICHLET:
                    self.rhs -= self.weights_general[2] * self.ids[2].get_bc_value()
                else:
                    self.weights[2]  = self.weights_general[2]
            
            if (self.c_id + 1)%nodes_x != 1: #defines the left id
                self.ids[3] = self.ids[0] - 1
                if self.ids[3] == bc.DIRICHLET:
                    self.rhs -= self.weights_general[3] * self.ids[3].get_bc_value()
                else:
                    self.weights[3]  = self.weights_general[3]

            if (self.c_id + 1)%nodes_x != 0: #defines the right id
                self.ids[4] = self.ids[0] + 1
                if self.ids[4] == bc.DIRICHLET:
                    self.rhs -= self.weights_general[4] * self.ids[1].get_bc_value()
                else:
                    self.weights[4]  = self.weights_general[4]

        else: #TODO Neumann conditions
            a = 1#complete here

        return self.ids, self.weights, self.rhs 

    def reset_stencil(self):
        self.ids = -1*np.ones(5)
        self.weights = np.zeros(5)
        self.rhs = 0.0




class System:
    def __init__(self, domain, stencil, rhs_func):
        self.domain = domain
        self.stn = stencil
        num_nodes = self.domain.get_num_nodes()
        self.A = np.zeros((num_nodes, num_nodes))
        self.b = np.array([rhs_func(domain.list_nodes[i].point) for i in range(domain.get_num_nodes())]) #TODO encapsular isso?

    def assemble(self): #TODO: talvez implementar isso com produto de kronecker e remodelar toda a logica seja uma opcao
        self.stn.reset()
        for i, it in enumerate(self.domain.node_list):
            self.stn.center_node(i)
            ids, weights, rhs = self.stn.assemble_stencil()
            for j in range(5):
                if ids[j] != -1:
                    self.A[i, j] = weights[j]
            self.b[i] = rhs

        
                





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
