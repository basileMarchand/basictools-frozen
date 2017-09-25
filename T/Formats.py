# -*- coding: utf-8 -*-

import numpy as np

__author__ = "Felipe Bordeu"
from BasicTools.Helpers.TextFormatHelper import TFormat
# name comp is reserved for the u_comp = [u_1 u_2 u_3] = [u v w]

class BaseTensor(object):
    def __init__(self):
        self.cores = [];

    def GetAxisNames(self):
        res = [];
        for c in self.cores:
            res.extend(c.GetAxisNames())
        return list(set(res))

    def GetNumberOfDims(self):
        return len(self.GetAxisNames())

    def AddSubTensor(self,Gs):
        import copy
        self.cores.append(copy.deepcopy(Gs))



    def Restriction(self,part={}):
        """ Particularisation of the canonic Tensor for some dimension , part is a dictionary with the names and the index to particularize"""
        " part = {'x':[1,3,4], 'comp':[0]}"
        import math
        import copy as cp
        res = cp.deepcopy(self)
        for c in res.cores:
            index = []
            for n,cpt in zip(c.names,list(range(c.ndims))):

                if list(part.keys()).count(n) :
                    #c.array = c.array.slice()
                    index.append(slice(part[n],int(part[n]+math.copysign(1,part[n])), int(math.copysign(1,part[n])) ) )
                else:
                    index.append(slice(None))
            c.array = c.array[index];
        return res



    def __str__(self):
        res =   str(type(self)).split(".")[-1].split("'")[0] + "\n"
        TFormat().II();
        CI = getattr(self, "CheckIntegrity", None)
        if CI != None:
            if self.CheckIntegrity(verbose=False):
                res +=TFormat().GetIndent() + 'Integrity OK \n';
            else:
                res += TFormat.InRed(TFormat().GetIndent() +'Integrity NOT OK !!!')+'\n';
        res += TFormat().GetIndent() + "GetAxisNames : " + str(", ".join(self.GetAxisNames()))
        cpt = 0
        for c in self.cores:
            res += "\n"+ TFormat().GetIndent() +str(cpt) +") "
            TFormat.II();
            res += c.__str__()
            TFormat.DI();
            cpt +=1

        TFormat.DI();
        return res



class FullTensor(BaseTensor):
    """ a numpy darray with names """
    def __init__(self, nRanks=2 , data = None):
        super(type(self),self).__init__()
        self.ndims = nRanks;
        if data is  None:
            self.array = np.ndarray(shape=(1,)*nRanks,dtype=float);
        else:
            self.array = data;
        self.names = ['']*nRanks;

    def SetRanksNames(self, names):
        assert(len(names) == self.ndims)
        self.names = names

    def GetAxisNames(self):
        return self.names;

    def __str__(self):
        res =  'FullTensor : '
        if  len(self.names):
            res += "" + ", ".join(self.names)
            res += " ("+ (",".join([str(x) for x in self.array.shape]) )+")"
        return res

class SparseBaseTensor(BaseTensor):
    def __init__(self):
        super(SparseBaseTensor,self).__init__()


    def CheckIntegrity(self,verbose=True):
        # first we get the list of all the Spaces
        sizes = {};
        # now for each SubTensor we check the size of the ranks
        res = True;
        for c in self.cores:
            for n,s in zip(c.names,c.array.shape):
                if n in sizes:
                    if (sizes[n] != s) :
                        res = False;
                        if  verbose:
                            print("Tensor Integrity check fail! : dimension '" + n + "' is inconsistent")

                else:
                   sizes[n] = s
        if res and verbose : print ("Tensor Integrity ok")
        return res
    def Reduction(self,redRanks=[]):
        import numpy as np
        # first we check if the reduction will generate a full tensor
        # we never do the reduction in the comp
        # for the moment we dont know how to treat the comp

        dic = 'abcdefghijklmnopqrstuvwxyz'

        # generation of the dictionary
        spaces = self.GetAxisNames()
        #print(spaces)
        dictionary = dict(list(zip(spaces, dic[0:len(spaces)])))

        # to count the number of
        singleDouble = {};
        script = ''
        for G in self.cores:
             script += "".join([dictionary[x] for x in  G.names  ])
             script += ","
             for y in G.names:
                 if y in singleDouble:
                     singleDouble[y] += 1
                 else:
                     singleDouble[y] = 1

        #print([x for x in singleDouble if singleDouble[x] == 1])
        script = script[0:-1] + "->"

        # Generation of the operation
        # the independent index are concatenated in the same order
        #print('-----')
        #print(singleDouble)
        fullDim = 0;
        names = []
        for G in self.cores:
            for y in G.names:
                if singleDouble[y] == 1:
                    script += dictionary[y]
                    fullDim +=1;
                    singleDouble[y] = 0
                    names.append(y);

        print(script)
        data = np.einsum(script, *[ G.array for G in self.cores]);

        res = FullTensor(nRanks = fullDim)
        res.array = data
        res.SetRanksNames(names);
        return res


class CanonicTensor(SparseBaseTensor):
    def __init__(self):
        super(CanonicTensor,self).__init__()

    def AddSubTensor(self,data, name=None):
        if name is None:
            super(CanonicTensor, self).AddSubTensor(data)
            return

        G1 = FullTensor(nRanks=2)
        G1.SetRanksNames([name,'alpha1'])
        G1.array = data;
        super(CanonicTensor, self).AddSubTensor(G1)

    def __add__(self, other):
        if type(self) != type(other):
            return TensorSum(self,other)
        else:
            import copy
            res = copy.copy(self)
            if self is other :
                res *= 2
                return res
            else:
                res += other
                return res

    def __iadd__(self,other):
        if type(self) != type(other):
            raise Exception("can't do an inplace operation with different types of objects")# pragma: no cover
        else:
            if self is other :
                return self*2
            else:
                for G1,G2 in zip(self.cores,other.cores):
                    x = np.concatenate((G1.array, G2.array), axis=1)
                    G1.array = x
                return self

    def __mul__(self,other):

        if type(other) == float or type(other) == int:
            import copy
            res = copy.copy(self)
            res *= other
            return res
        else :
            raise Exception("Don't know how to multiply a tensor and a " + str(type(other)))


    def __imul__(self,other):
        if type(other) == float or type(other) == int:
            self.cores[0].array *=  other
            return self
        else:
            raise Exception("Don't know how to multiply a tensor and a " + str(type(other)) ) # pragma: no cover

class TensorTrain(SparseBaseTensor):
    def __init__(self):
        super(type(self),self).__init__()

    def AddSubTensor(self,data, name=None):
        if name is None:
            super(TensorTrain, self).AddSubTensor(data)
            return

        #if empty
        if len(self.cores) == 0 and len(data.shape) != 2: # pragma: no cover
            raise Exception("The first core must be of dimensionality 2")

        if len(self.cores) > 1  and self.cores[-1].GetNumberOfDims() == 2 :# pragma: no cover
            raise Exception("Can't add a new core (last core is only of rank 2)")

        if len(self.cores) == 0 :
            G1 = FullTensor(nRanks=2)
            G1.SetRanksNames([name,'alpha1'])
            G1.array = data;
            super(TensorTrain, self).AddSubTensor(G1)
            return

        if len(self.cores) > 0 :
            G1 = FullTensor(nRanks=len(data.shape))
            names = [name]
            if len(data.shape) == 3:
                names.append('alpha'+str(len(self.cores)+1))

            names.append(self.cores[-1].names[1])
            print( names)
            G1.SetRanksNames(names)
            G1.array = data;
            super(TensorTrain, self).AddSubTensor(G1)
            return

 # for the moment we use the reduction from the base class, slow but robust
#    def Reduction(self,redRanks=[]):
#        #for now only works fo
#        mat = [ G.array for G in self.cores]
#        mat.reverse()
#        print([ x.shape for x in mat])
#        data = reduce(np.dot,mat)
#        res = FullTensor(nRanks = len(self.cores))
#        res.array = data
#        res.SetRanksNames([X.names[0] for X in self.cores]);
#        return res

class TensorSum(BaseTensor):
    def __init__(self,first,second):
        super(type(self),self).__init__()
        if issubclass(type(first), SparseBaseTensor):
            self.AddSubTensor(first)
        else:
            raise Exception("TensorSum: First Object not a SparseBaseTensor class") # pragma: no cover

        if issubclass(type(second), SparseBaseTensor):
            self.AddSubTensor(second)
        else:
            raise Exception("TensorSum: Second Object not a SparseBaseTensor class") # pragma: no cover

    def GetAxisNames(self):
        res = [];
        for c in self.cores:
            res.extend(c.GetAxisNames())
        return list(set(res))
    def __add__(self,other):
        import copy
        res = copy.deepcopy(self)
        res += other
        return res
    def __iadd__(self,other):
        if type(other) == TensorSum:
            self.cores.extend(other.cores)
            return self
        elif issubclass(type(other), SparseBaseTensor):
            self.AddSubTensor(other)
            return self
        else:
            raise Exception("TensorSum.__iadd__: Don't know how add type  " + str(type(other)))


def CheckIntegrity():
    import numpy as np

    # Test Full tensors
    ft = FullTensor(data=np.ndarray((2,3)));
    name = ft.GetAxisNames();
    disp = ft.__str__()

    # Test SparseBaseTensor CheckIntegrity
    spt = SparseBaseTensor();
    for n in [0,1]:
        G = FullTensor(nRanks=2);
        G.SetRanksNames(["X"+str(n) ,'alpha1'])
        G.array = np.ndarray((2,3+n))
        spt.AddSubTensor(G)

    spt.CheckIntegrity()
    print(spt)

    spt.cores[1].array = np.ndarray((2,3))

    spt.Restriction({'X1':int(0)})
    spt.Reduction()


    # Test SparseBaseTensor CheckIntegrity
    ct = CanonicTensor()
    ct.AddSubTensor(spt.cores[0])
    ct.AddSubTensor(spt.cores[1].array,'X1')

    ct2 = CanonicTensor()
    ct2.AddSubTensor(spt.cores[0])
    ct2.AddSubTensor(spt.cores[1].array,'X1')

    # sum
    ct = ct + ct + ct2
    ct += ct

    # product
    ct = ct *2
    ct *= 2
    #TODO
    try :
        ct = ct * ct
        return 'Not ok' # pragma: no cover
    except Exception as e:
        pass

    #TODO
    try :
        ct = ct * ct2
        return 'Not ok' # pragma: no cover
    except Exception as e:
        pass
    # heterogenius
    tt = TensorTrain()
    print("1")
    print(tt)
    tt.AddSubTensor(spt.cores[0].array,'x0')
    print("2")
    print(tt)
    tt.AddSubTensor(np.empty([1, 1, 3]) ,'X1')
    print("3")
    print(tt)
    tt.AddSubTensor(np.empty([1, 1]) ,'X2')
    print("4")
    print(tt)



    su =  ct + tt
    print(su)
    su += ct2
    print(su)
    su = su +su
    print(su)
    try:
        su += 3
        return'Not OK' # pragma: no cover
    except Exception as e:
        pass

    return 'OK'

if __name__ == '__main__':
    CheckIntegrity() # pragma: no cover