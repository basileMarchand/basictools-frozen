# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

from BasicTools.FE.KR.KRBase import KRBase

class KRDofs(KRBase):
    def __init__(self):
        super(KRDofs,self).__init__()

    def GenerateEquations(self,meshI=None,fields=None,CH=None,meshII=None,fieldsII=None):
        if CH is None:
            from BasicTools.Linalg.ConstraintsHolder import ConstraintsHolder
            CH = ConstraintsHolder()

        ei,ej,ev = self.GetIsJsVs()
        CH.AddEquationsFromIJV(ei,ej,ev)
        return CH

    def GetIsJsVs(self):
        raise(Exception("Function not implemeneted, please subclass"))
        #example
        #ei= [0, 0, 1, 1, 2, 2, 3, 3]
        #ej= [1, 4, 2, 7,10,15, 9,12]
        #ev= [1,-1, 1,-1, 1,-1, 1,-1]
        return

def CheckIntegrity(GUI=False):

    class KRDofsTest(KRDofs):
        def __init__(self):
            super(KRDofsTest,self).__init__()

        def GetIsJsVs(self):
            ei= [0, 0, 1, 1, 2, 2, 3, 3]
            ej= [1, 4, 2, 7,10,15, 9,12]
            ev= [1,-1, 1,-1, 1,-1, 1,-1]
            return ei, ej, ev

    obj = KRDofsTest()
    ch = obj.GenerateEquations()
    ch.SetNumberOfDofs(20)
    ch.Compact()
    mat, dofs = ch.ToSparse()
    print(mat.toarray())
    print(dofs)
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(GUI=True))