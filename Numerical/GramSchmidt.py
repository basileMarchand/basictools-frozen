import numpy as np

class GramSchmidt(object):
    """Modified Gram-Schmidt orthonormalization routine, with user-provided scalar product
    
    """
    
    def __init__(self, A, v):
        self.A = A
        self.v = v
        shape = v.shape
        self.nbVec   = shape[0]
        self.sizeVec = shape[1]
        self.u = np.zeros((self.nbVec,self.sizeVec), dtype = float)
        self.e = np.zeros((self.nbVec,self.sizeVec), dtype = float)

    def ScalProd(self, v, u):
        return np.dot(v, self.A.dot(u))

    def Proj(self, v, u):
        return (self.ScalProd(v, u)/self.ScalProd(u, u))*u

    def Compute(self):
        for i in xrange(self.nbVec):
            self.u[i] = self.v[i]
            for j in xrange(i):
                self.u[i] = self.u[i] - self.Proj(self.u[i], self.u[j])
            self.e[i] = self.u[i]/np.sqrt(self.ScalProd(self.u[i], self.u[i]))
        return

    def Check(self):
        import OTTools.Helpers.TextFormatHelper as TFH
        matCheck = np.zeros((self.nbVec,self.nbVec), dtype = float)
        for i in xrange(self.nbVec):
            for j in xrange(self.nbVec):
                matCheck[i,j] = self.ScalProd(self.e[i], self.e[j])
        matCheck -= np.eye(self.nbVec)
        
        error = np.linalg.norm(matCheck)/np.linalg.norm(np.eye(self.nbVec))
        print("Gram-Schmidt relative error = "+ TFH.TFormat.GoodBad(error,error<1.e-10))
        return
        

def CheckIntegrity():
    v  = np.random.rand(3,10)
    A  = np.eye(10)
    GS = GramSchmidt(A, v)
    GS.Compute()
    GS.Check()
    return 'ok'
        
        
if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover   

    
    
    
    
