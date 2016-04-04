import numpy as np

class GramSchmidt():
    """Modified Gram-Schmidt orthonormalization routine, with user-provided scalar product
    
    """
    
    def __init__(self):
        self.sizeVec = None
        self.nbVec = None
        self.u = None
        self.e = None

    def Init(self, v):
        shape = v.shape
        self.nbVec   = shape[0]
        self.sizeVec = shape[1]
        self.u = np.zeros((self.nbVec,self.sizeVec), dtype = float)
        self.e = np.zeros((self.nbVec,self.sizeVec), dtype = float)
        return

    def ScalProd(self, A, v, u):
        return np.dot(v, A.dot(u))

    def Proj(self, A, v, u):
        return (self.ScalProd(A, v, u)/self.ScalProd(A, u, u))*u

    def Compute(self, A, v):
        self.Init(v)
        for i in xrange(self.nbVec):
            self.u[i] = v[i]
            for j in xrange(i):
                self.u[i] = self.u[i] - self.Proj(A, self.u[i], self.u[j])
            self.e[i] = self.u[i]/np.sqrt(self.ScalProd(A, self.u[i], self.u[i]))
        return

    def Check(self, A):
        import OTTools.Helpers.TextFormatHelper as TFH
        matCheck = np.zeros((self.nbVec,self.nbVec), dtype = float)
        for i in xrange(self.nbVec):
            for j in xrange(self.nbVec):
                matCheck[i,j] = self.ScalProd(A, self.e[i], self.e[j])
        matCheck -= np.eye(self.nbVec)
        
        print("Gram-Schmidt relative error = "+ TFH.TFormat.GoodBad(np.linalg.norm(matCheck)/np.linalg.norm(np.eye(self.nbVec)),1.e-10))
        return
        

def CheckIntegrity():
    vec = np.array([np.random.rand(10), np.random.rand(10), np.random.rand(10)])
    A = np.eye(10)
    GS = GramSchmidt()
    GS.Compute(A,vec)
    GS.Check(A)
    return 'ok'
        
        
if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover   

    
    
    
    
