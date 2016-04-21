import numpy as np

        
def TruncatedSymetricSVD(a, epsilon):        
    
    eigenValues, eigenVectors = np.linalg.eigh(a,UPLO='L')

    idx = eigenValues.argsort()[::-1]   
    eigenValues = eigenValues[idx]
    eigenVectors = eigenVectors[:,idx]

    id_max = 0
    bound = (epsilon**2)*eigenValues[0]
    for e in eigenValues:
      if e > bound:
        id_max += 1
    id_max2 = 0
    bound = (1-epsilon**2)*np.sum(eigenValues)
    temp = 0
    for e in eigenValues:
      temp += e
      if temp < bound:
        id_max2 += 1
    id_max = max(id_max, id_max2)

    return eigenValues[0:id_max], eigenVectors[:,0:id_max]

    
def CheckIntegrity():
    a = np.random.rand(10,10)
    a = np.dot(a.T,a)
    TruncatedSymetricSVD(a, 1.e-6)
    return 'ok'
        
        
if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover   

    
    
    
    
