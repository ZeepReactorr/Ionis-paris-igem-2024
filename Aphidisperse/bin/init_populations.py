from dataclasses import dataclass
from shapely import Point, Polygon
import math as m
import itertools

@dataclass 
class Object:
    position : Point
    
    def __init__(self, position):
        self.position = position

    def update(self):
        pass

distance = 0.2

@dataclass
class Population(Object):
    local_density : float
    treshold_date : int

    def update(self):
        if self.local_density >= 1:
            self.local_density = 0.9
            
            poly_proba = [Population(Point(self.position.x + distance*m.sin(m.radians(45*i)), 
                                       self.position.y + distance*m.cos(m.radians(45*i))), 0.1, self.treshold_date) for i in range(0, 8)]
            out = list(itertools.chain.from_iterable([[self], poly_proba]))
            return out

        else:
            self.local_density *= 1.16576#value to change by the result of the core equation
            return [self]
        
class Culture(Object):
    _type : str
    _status : str

    def update(self):
        pass