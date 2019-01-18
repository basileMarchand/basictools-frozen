
# -*- coding: utf-8 -*-
#python octree implementation
# Code © Spencer Krum June 2011
# Released underl GPLv3 See LICENSE file in this repository

class node():
    """
    Class to be a node in my octree
    """

    def __init__(self,parent, Xupperlimit, Yupperlimit, Zupperlimit, Xlowerlimit, Ylowerlimit, Zlowerlimit):
        self.parent = parent
        self.Xupperlimit = Xupperlimit
        self.Yupperlimit = Yupperlimit
        self.Zupperlimit = Zupperlimit
        self.Xlowerlimit = Xlowerlimit
        self.Ylowerlimit = Ylowerlimit
        self.Zlowerlimit = Zlowerlimit
        self.Xcenter = (self.Xupperlimit + self.Xlowerlimit)/2.
        self.Ycenter = (self.Yupperlimit + self.Ylowerlimit)/2.
        self.Zcenter = (self.Zupperlimit + self.Xlowerlimit)/2.

        self.posXposYposZ = None
        self.posXposYnegZ = None
        self.posXnegYposZ = None
        self.posXnegYnegZ = None
        self.negXposYposZ = None
        self.negXposYnegZ = None
        self.negXnegYposZ = None
        self.negXnegYnegZ = None

        self.value = None

    def __str__(self):
        res = ""
        if self.value is not None:
            res += "("+str(len(self.value))+")\n"
            return res

        res += "(A)" if self.posXposYposZ is None else str(self.posXposYposZ)
        res += "(B)" if self.posXposYnegZ is None else str(self.posXposYnegZ)
        res += "(C)" if self.posXnegYposZ is None else str(self.posXnegYposZ)
        res += "(D)" if self.posXnegYnegZ is None else str(self.posXnegYnegZ)
        res += "(E)" if self.negXposYposZ is None else str(self.negXposYposZ)
        res += "(F)" if self.negXposYnegZ is None else str(self.negXposYnegZ)
        res += "(G)" if self.negXnegYposZ is None else str(self.negXnegYposZ)
        res += "(H)" if self.negXnegYnegZ is None else str(self.negXnegYnegZ)


        return res


    def add(self, payload, coord, level):

        """
        Create a subnode
        """


        if level != 0 :
            level -= 1
            #Determine quadrant
            if coord[0] <= self.Xcenter:
                Xupperlimit = self.Xcenter
                Xlowerlimit = self.Xlowerlimit

                #negX
                if coord[1] <= self.Ycenter:
                    #negY
                    Yupperlimit = self.Ycenter
                    Ylowerlimit = self.Ylowerlimit

                    if coord[2] <= self.Zcenter:
                        #negZ
                        Zupperlimit = self.Zcenter
                        Zlowerlimit = self.Zlowerlimit
                        if self.negXnegYnegZ is None:
                            self.negXnegYnegZ = node(self, Xupperlimit, Yupperlimit, Zupperlimit, Xlowerlimit, Ylowerlimit, Zlowerlimit)
                        self.negXnegYnegZ.add(payload, coord, level)
                    else:
                        #posZ
                        Zupperlimit = self.Zupperlimit
                        Zlowerlimit = self.Zcenter
                        if self.negXnegYposZ is None:
                            self.negXnegYposZ = node(self, Xupperlimit, Yupperlimit, Zupperlimit, Xlowerlimit, Ylowerlimit, Zlowerlimit)
                        self.negXnegYposZ.add(payload, coord, level)
                else:
                    #posY
                    if coord[2] <= self.Zcenter:
                        #negZ
                        Yupperlimit = self.Yupperlimit
                        Zupperlimit = self.Zcenter
                        Ylowerlimit = self.Ycenter
                        Zlowerlimit = self.Zlowerlimit
                        if self.negXposYnegZ is None:
                            self.negXposYnegZ = node(self, Xupperlimit, Yupperlimit, Zupperlimit, Xlowerlimit, Ylowerlimit, Zlowerlimit)
                        self.negXposYnegZ.add(payload, coord, level)

                    else:
                        #posZ
                        Yupperlimit = self.Yupperlimit
                        Zupperlimit = self.Zupperlimit
                        Ylowerlimit = self.Ycenter
                        Zlowerlimit = self.Zcenter
                        if self.negXposYposZ is None:
                            self.negXposYposZ = node(self, Xupperlimit, Yupperlimit, Zupperlimit, Xlowerlimit, Ylowerlimit, Zlowerlimit)
                        self.negXposYposZ.add(payload, coord, level)


            else:
                #posX
                if coord[1] <= self.Ycenter:
                    #negY
                    if coord[2] <= self.Zcenter:
                        #negZ
                        Xupperlimit = self.Xupperlimit
                        Yupperlimit = self.Ycenter
                        Zupperlimit = self.Zcenter
                        Xlowerlimit = self.Xcenter
                        Ylowerlimit = self.Ylowerlimit
                        Zlowerlimit = self.Zlowerlimit
                        if self.posXnegYnegZ is None:
                            self.posXnegYnegZ = node(self, Xupperlimit, Yupperlimit, Zupperlimit, Xlowerlimit, Ylowerlimit, Zlowerlimit)
                        self.posXnegYnegZ.add(payload, coord, level)

                    else:
                        #posZ
                        Xupperlimit = self.Xupperlimit
                        Yupperlimit = self.Ycenter
                        Zupperlimit = self.Zupperlimit
                        Xlowerlimit = self.Xcenter
                        Ylowerlimit = self.Ylowerlimit
                        Zlowerlimit = self.Zcenter
                        if self.posXnegYposZ is None:
                            self.posXnegYposZ = node(self, Xupperlimit, Yupperlimit, Zupperlimit, Xlowerlimit, Ylowerlimit, Zlowerlimit)
                        self.posXnegYposZ.add(payload, coord, level)

                else:
                    #posY
                    if coord[2] <= self.Zcenter:
                        #negZ
                        Xupperlimit = self.Xupperlimit
                        Yupperlimit = self.Yupperlimit
                        Zupperlimit = self.Zcenter
                        Xlowerlimit = self.Zcenter
                        Ylowerlimit = self.Ycenter
                        Zlowerlimit = self.Zlowerlimit
                        if self.posXposYnegZ is None:
                            self.posXposYnegZ = node(self, Xupperlimit, Yupperlimit, Zupperlimit, Xlowerlimit, Ylowerlimit, Zlowerlimit)
                        self.posXposYnegZ.add(payload, coord, level)

                    else:
                        #posZ
                        Xupperlimit = self.Xupperlimit
                        Yupperlimit = self.Yupperlimit
                        Zupperlimit = self.Zupperlimit
                        Xlowerlimit = self.Xcenter
                        Ylowerlimit = self.Ycenter
                        Zlowerlimit = self.Zcenter
                        if self.posXposYposZ is None:
                            self.posXposYposZ = node(self, Xupperlimit, Yupperlimit, Zupperlimit, Xlowerlimit, Ylowerlimit, Zlowerlimit)
                        self.posXposYposZ.add(payload, coord, level)

        else :
            #if level == 0:
            try:
                self.value.append((coord,payload))
            except AttributeError:
                self.value = []
                self.value.append((coord,payload))







class Octree():
    """
    class to hold the whole tree
    """

    def __init__(self, Xmax, Ymax, Zmax, Xmin, Ymin, Zmin, root_coords=(0,0,0), maxiter=7):
        self.Xmax = Xmax
        self.Ymax = Ymax
        self.Zmax = Xmax
        self.Xmin = Xmin
        self.Ymin = Ymin
        self.Zmin = Zmin
        self.root_coords = root_coords
        self.maxiter = maxiter

        self.root = node('root', Xmax, Ymax, Zmax, Xmin, Ymin, Zmin)

    def add_item(self, payload, coord):
        """
        Create recursively create subnodes until maxiter is reached
        then deposit payload in that node
        """

        self.root.add(payload, coord, self.maxiter)

    def find_within_range(self, center, size, shape):
        """
        Return payloads and coordinates of every payload within
        a specified area
        """
        """
        When shape is cube:
        Search space is defined as the cubic region where each face is 'size'
        distance directly away from the center.
        """
        """
        Should support "cube", "sphere", "doughnut"
        """
        if shape == "cube":
            return self.find_within_range_cube( center, size)

    def find_within_range_cube(self, center, size):
            """
            This deals with things around the center of a node in a box shape
            with a radius of 'size'
            It would be totally good to make a spere search space
            """
            """
            It works by making the (correct I think) assumption that the box
            shape has 8 vertices. To determine the overlap between the search
            space and the node area is hard. We can start by identifying which
            of the 8 children of the current node overlap with the search space.
            The assumption is that if the search box overlaps in any part with
            a child, then the most extreme vertex will be within the child.

            Let me put it 2 dimensions with ascii art:

            Fig1. A Quadtree Node Space with children space labeled

            1,-1--------1,0--------1,1
            |           |           |
            |  child_4  |  child_1  |
            |           |           |
           -1,0---------0,0--------1,0
            |           |           |
            |  child_3  |  child_2  |
            |           |           |
          (-1,-1)-----(-1,0)-----(-1,1)


            The four children of this node are:

            child_1 has the following attributes

            . node name is posXposY
            . initial value is Null
            . Once recursively filled out, its value is another node
            . it represents the rectangular area between ((0,0),(1,0)) and ((1,0), (1,1))
            . its center is at the center of the rectangular area it represents


            Children _2, _3, _4 are much the same

            Fig2. Quadtree Node Space with superimposed search space


            /-----------------------\
            |           |           |
            |           |           |
            |     |-------|         |
            |-----|-------|---------|
            |     |     | |         |
            |     |-------|         |
            |           |           |
            \-----------------------/


            Fig3. Fig2 Zoomed

            /---------------------------------------------------------------\
            |                                                ||   <child_1> |
            |    <child_4>                                   ||  <point A>  |
            |                                                ||        |    |
            |                <selection space>               ||        V    |
            |    0+++++++++++++++++++++++++++++++++++++++++++++++++++++0    |
            |    +                                           ||<exact  +    |
            |    +                             <x-axis>      || center>+    |
            |====+===========================================%%========+====|
            |    +                                           ||        +    |
            |    +                                           ||        +    |
            |    +                                           ||        +    |
            |    +                                          ^||        +    |
            |    +                                          y||        +    |
            |    +                                          -||        +    |
            |    +                                          a||        +    |
            |    +                                          x||        +    |
            |    +                                          i||        +    |
            |    +                                          s||        +    |
            |    +                                          v||        +    |
            |    +                                           ||        +    |
            |    +                                           ||        +    |
            |    0+++++++++++++++++++++++++++++++++++++++++++++++++++++0    |
            |             <child_3>                          ||   <child_2> |
            \---------------------------------------------------------------/


            Fig4. Quadtree Node Space with different superimposed search space


            /-----------------------\
            |           |           |
            |           |           |
            |     |---| |           |
            |-----|---|-------------|
            |     |   | |           |
            |     |---| |           |
            |           |           |
            \-----------------------/


            Fig5. Fig2 Zoomed

            /---------------------------------------------------------------\
            |                                                ||             |
            |    <child_4>                                   ||  <child_1>  |
            |                                    <point B>   ||             |
            |                <selection space>   |           ||             |
            |    0+++++++++++++++++++++++++++0 <--           ||             |
            |    +                           +               ||<exact       |
            |    +                           + <x-axis>      || center>     |
            |====+===========================+===============%%=============|
            |    +                           +               ||             |
            |    +                           +               ||             |
            |    +                           +               ||             |
            |    +                           +              ^||             |
            |    +                           +              y||             |
            |    +                           +              -||             |
            |    +                           +              a||             |
            |    +                           +              x||             |
            |    +                           +              i||             |
            |    +                           +              s||             |
            |    +                           +              v||             |
            |    +                           +               ||             |
            |    +                           +               ||             |
            |    0+++++++++++++++++++++++++++0               ||             |
            |             <child_3>                          ||   <child_2> |
            \---------------------------------------------------------------/



            What we see in Fig2/3 vs Figi4/5 is illustrated by point A, point B
            Points A, B are the points that are the most positive in both X and
            Y. We see that if and only if this 'upper rightmost' point is within
            child 1, then there is space in the selection overlapping part of
            child 1. This also follows for children_2, _3, _4 and generalizes to
            the 3D space of the octree.

            So our litmus test for 'does the selection overlap with a
            childspace?' is whether or not the most extreme(bad wording I know)
            part of that selection is within the childspace.

            """

            list_list = []
            list_list.append([self.root])
            for level in range(self.maxiter):
                list_list.append([])

            if hasattr(size,"__iter__"):
                print(size)
                print(type(size))
                Xedge_max = center[0] + size[0]
                Xedge_min = center[0] - size[0]
                Yedge_max = center[1] + size[1]
                Yedge_min = center[1] - size[1]
                Zedge_max = center[2] + size[2]
                Zedge_min = center[2] - size[2]
            else:
                Xedge_max = center[0] + size
                Xedge_min = center[0] - size
                Yedge_max = center[1] + size
                Yedge_min = center[1] - size
                Zedge_max = center[2] + size
                Zedge_min = center[2] - size

            corner0 = (Xedge_max, Yedge_max, Zedge_max)
            corner1 = (Xedge_max, Yedge_max, Zedge_min)

            corner2 = (Xedge_max, Yedge_min, Zedge_max)
            corner3 = (Xedge_max, Yedge_min, Zedge_min)

            corner4 = (Xedge_min, Yedge_max, Zedge_max)
            corner5 = (Xedge_min, Yedge_max, Zedge_min)
            corner6 = (Xedge_min, Yedge_min, Zedge_max)
            corner7 = (Xedge_min, Yedge_min, Zedge_min)
            #print list_list
            for level in range(self.maxiter):
                for node in list_list[level]:
                    if Xedge_max > node.Xcenter :
                        if corner0[1] > node.Ycenter :
                            if node.posXposYposZ is not None:
                                if  corner0[2] > node.Zcenter:
                                #if corner0[0] > node.Xcenter and corner0[1] > node.Ycenter  and corner0[2] > node.Zcenter:
                                #table = ((corner0[0] > node.Xcenter),(corner0[1] > node.Ycenter) ,(corner0[2] > node.Zcenter))
                                #if not False in table :
                                    list_list[level+1].append(node.posXposYposZ)
                            if node.posXposYnegZ is not None:
                                if corner1[2] <= node.Zcenter:
                                #if corner1[0] > node.Xcenter and corner1[1] > node.Ycenter and corner1[2] <= node.Zcenter:
                                #table = ((corner1[0] > node.Xcenter),(corner1[1] > node.Ycenter) ,(corner1[2] <= node.Zcenter))
                                #if not False in table:
                                    list_list[level+1].append(node.posXposYnegZ)
                        if corner2[1] <= node.Ycenter :
                            if node.posXnegYposZ is not None:
                                if corner2[2] > node.Zcenter:
                                #if corner2[0] > node.Xcenter and corner2[1] <= node.Ycenter and corner2[2] > node.Zcenter:
                                #table = ((corner2[0] > node.Xcenter),(corner2[1] <= node.Ycenter) ,(corner2[2] > node.Zcenter))
                                #if not False in table:
                                    list_list[level+1].append(node.posXnegYposZ)
                            if node.posXnegYnegZ is not None:
                                if corner3[2] <= node.Zcenter:
                                #if corner3[0] > node.Xcenter and corner3[1] <= node.Ycenter and corner3[2] <= node.Zcenter:
                                #table = ((corner3[0] > node.Xcenter),(corner3[1] <= node.Ycenter) ,(corner3[2] <= node.Zcenter))
                                #if not False in table:
                                    list_list[level+1].append(node.posXnegYnegZ)


                    if corner4[0] <= node.Xcenter:
                        if corner4[1] > node.Ycenter:
                            if node.negXposYposZ is not None:
                                if  corner4[2] > node.Zcenter:
                                #if corner4[0] <= node.Xcenter and corner4[1] > node.Ycenter and corner4[2] > node.Zcenter:
                                #table = ((corner4[0] <= node.Xcenter),(corner4[1] > node.Ycenter) ,(corner4[2] > node.Zcenter))
                                #if not False in table:
                                    list_list[level+1].append(node.negXposYposZ)
                            if node.negXposYnegZ is not None:
                                if corner5[2] <= node.Zcenter:
                                #if corner5[0] <= node.Xcenter and corner5[1] > node.Ycenter and corner5[2] <= node.Zcenter:
                                #table = ((corner5[0] <= node.Xcenter),(corner5[1] > node.Ycenter) ,(corner5[2] <= node.Zcenter))
                                #if not False in table:
                                    list_list[level+1].append(node.negXposYnegZ)

                        if corner6[1] <= node.Ycenter:
                            if node.negXnegYposZ is not None:
                                if corner6[2] > node.Zcenter:
                                #if corner6[0] <= node.Xcenter and corner6[1] <= node.Ycenter and corner6[2] > node.Zcenter:
                                #table = ((corner6[0] <= node.Xcenter),(corner6[1] <= node.Ycenter) ,(corner6[2] > node.Zcenter))
                                #if not False in table:
                                    list_list[level+1].append(node.negXnegYposZ)
                            if node.negXnegYnegZ is not None:
                                if corner7[2] <= node.Zcenter:
                                #if corner7[0] <= node.Xcenter and corner7[1] <= node.Ycenter and corner7[2] <= node.Zcenter:
                                #table = ((corner7[0] <= node.Xcenter),(corner7[1] <= node.Ycenter) ,(corner7[2] <= node.Zcenter))
                                #if not False in table:
                                    list_list[level+1].append(node.negXnegYnegZ)


            #flaten the output
            return [ item for sublist in list_list[-1] for item in  sublist.value]


def CheckIntegrity():

    print( "Creating octree")
    tree = Octree(100,100,100, -100, -100, -100)
    print( "inserting node")
    tree.add_item("derp1", (90.34251,10.1234,10.9876))
    print( "Great success")
    print( "inserting node")
    tree.add_item("derp2", (10.34251,10.1234,10.9876))
    print( "Great success")
    print( "inserting node")
    tree.add_item("derp3", (-10.34251,10.1234,10.9876))
    print( "Great success")
    print( "inserting node")
    tree.add_item("derp4", (10.34251,-10.1234,10.9876))
    print( "Great success")
    print( "inserting node")
    tree.add_item("derp5", (10.34251,10.1234,-10.9876))
    print( "Great success")

    print(tree.root)
    print("-------")
    #get some data
    entries = tree.find_within_range((0,0,0), 40, "cube")
    entries = tree.find_within_range((0,0,0), [40]*3, "cube")
    cpt =0
    for i in entries:
      print(cpt),
      print(i)
      cpt += 1

    print(len(i))
    if len(entries) != 4:# pragma: no cover
        raise(Exception("Error") )

    return "ok"

if __name__ == '__main__':
    print((CheckIntegrity()))# pragma: no cover

