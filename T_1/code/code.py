import numpy as np
import matplotlib.pyplot as pyplot
from enum import Enum


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
    def __init__(self, x, y, id, bc):
        self.point = Point(x,y)
        self.id = id
        self.bc = bc

    def __str__(self):
        return f"Coordinate: {self.point}\nid: {self.id}\nBoundary condition: {self.bc.name}"
    
    def up(self, num_nodes_x):
        return self.id + num_nodes_x
    
    def down(self, num_nodes_y):
        return self.id - num_nodes_y
    
    def left(self):
        return self.id - 1
    
    def right(self):
        return self.id + 1
    



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
    
    #creates a node list with natural ordering
    def create_node_list(self):
        k = 0
        for i in range(self.num_nodes_y):
            y = self.low.y + i*h

            for j in range(self.num_nodes_x):
                x = self.low.x + j*h
                self.node_list.append(Node(x,y, k, bc.NONE))
                k += 1

h = 0.1

d = Domain(Point(0.0, -1.0), Point(1.0, 1.0), h)
d.create_node_list()
for i in range(d.num_nodes):
    print(d.node_list[i])