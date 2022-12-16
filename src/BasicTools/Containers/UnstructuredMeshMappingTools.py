import math
import numpy as np

from scipy.optimize import fmin_l_bfgs_b

from BasicTools.Containers.Filters import ElementFilter
from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh
from BasicTools.NumpyDefs import PBasicFloatType

class FoldOverFreeMaps():

    def __init__(self,mesh: UnstructuredMesh, boundaryEF:ElementFilter):
        self.eps = 1
        self.pdim = mesh.GetPointsDimensionality()
        self.edim = mesh.GetElementsDimensionality()
        self.nbNodes = mesh.GetNumberOfNodes()
        self.mesh = mesh
        self.nodeMask = np.zeros((self.nbNodes,self.pdim),dtype=bool)
        for name, data , ids in boundaryEF:
            nbIds = np.unique(data.connectivity[ids,:])
            self.nodeMask[nbIds,:] = True
        from BasicTools.FE.IntegrationsRules import NodalEvaluationP1
        from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceP1
        self.space = LagrangeSpaceP1
        self.ip = NodalEvaluationP1
        self.spaceAtIP = { name:self.space[name].SetIntegrationRule(ipv,ipw) for name,(ipv,ipw) in  NodalEvaluationP1.items()}
        self.theta = 0.5

        self.callback = None

    def CallCallback(self):
        if self.callback is not None:
            self.callback(self.mesh)

    def Start(self):
        n =0
        while True: # outer L-BFGS loop, Alg. 1, line 2
            [Fprev, grad] = self.ComputePotential(self.mesh.nodes.flatten()) # compute $F(U^k, \varepsilon^k)$
            #print(f"Fprev {Fprev}")#; exit()
            [newNodes, F, _] = fmin_l_bfgs_b(self.ComputePotential, self.mesh.nodes.flatten(),factr=1e12)  # inner L-BFGS loop, compute $F(U^{k+1}, \varepsilon^k)$
            #newNodes = mesh.nodes.flatten() - grad/(np.linalg.norm(grad))*0.01/energy.eps

            #print(newNodes)
            newNodes.shape = self.mesh.nodes.shape
            self.mesh.nodes = newNodes
            mindet = self.MinDetJacobian()
            mu = min(F/Fprev, 9e-1)*(mindet + math.sqrt(self.eps**2 + mindet**2))/2 # $\mu^k$, Eq. (6)
            self.eps = 2*math.sqrt(mu*(mu-mindet)) if mindet<mu else 0 # $\varepsilon^{k+1}$, Eq. (6)
            self.CallCallback()
            if (mindet>0 and F>999e-3*Fprev):
                break # stopping criterion, Alg. 1, line 11
            print(f"{n}, {self.eps:.5f}, {Fprev:.5f}, {mindet:.5f}" )
            n += 1
            if n > 200:
                break

    def ComputePotential(self,U): # compute the energy $F$ and its gradient $\nabla F$ for the map $\vec{u}$, Eq. (5)

        F = 0  # zero out the energy and the gradient
        G = np.zeros( (self.nbNodes,self.pdim), dtype=PBasicFloatType )
        nodes = U.view()
        nodes.shape = self.mesh.nodes.shape

        for name, data , ids in ElementFilter(self.mesh,dimensionality=self.edim):
            ipv, ipw = self.ip[name]
            space = self.space[name]
            sAtIP = self.spaceAtIP[name]
            eps2 =self.eps**2
            for id in ids:
                elConnectivity = data.connectivity[id,:]
                for qc in range(len(ipw)): # for every quad corner
                    J, det, jinv = sAtIP.GetJackAndDetI(qc,nodes[elConnectivity,:])
                    J =J.T
                    chi  = det/2 + math.sqrt(eps2 + det**2)/2    # the penalty function $\chi(\varepsilon^k, \mathrm{det}\ J)$
                    chip = .5 + det/(2*math.sqrt(eps2 + det**2)) # its derivative $\chi'(\varepsilon^k, \mathrm{det}\ J)$
                    #tjtj =
                    tjtj = np.sum(J*J) #♀np.trace(np.transpose(J).dot(J))
                    f = tjtj/(chi**(2/self.edim)) # quad corner shape quality
                    g = (det**2+1)/chi
                    F += (1-self.theta)*f # add the term to the energy
                    F +=    self.theta*g
                    dfdj = (1-self.theta)*(2*J - det*jinv(np.eye(2))*f*chip)/chi #$\frac{\partial f_\varepsilon}{\partial a}$: derivative w.r.t the Jacobian
                    if self.theta != 0:
                        dgdj = (self.theta)*det*jinv( ((2*det-g*chip)/chi) )
                        dfdj += dgdj
                    funcDer = sAtIP.valdphidxi[qc].T
                    dfdu = funcDer.dot(np.transpose(dfdj)) # chain rule for the real variables (Appendix A)
                    G[elConnectivity,:] += dfdu[:,:]

        G[self.nodeMask] = 0
        return F, G.flatten()


    def MinDetJacobian(self):
        res  = np.inf
        for name, data , ids in ElementFilter(self.mesh,dimensionality=self.edim):
            ipv, ipw = self.ip[name]
            sAtIP = self.spaceAtIP[name]
            lenIP = len(ipw)
            for i in ids:
                for qc in range(lenIP): # for every quad corner
                    J, det, jinv = sAtIP.GetJackAndDetI(qc,self.mesh.nodes[data.connectivity[i,:],:])
                    res = min(res,det)
        return res

def CheckIntegrity(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateSquare

    n = 8
    mesh = CreateSquare(dimensions=[n,n],origin=[0,0],spacing=[0.875/(n-1),1.75/(n-1)],ofTriangles=True)

    mask = mesh.nodes[:,1] >.8625
    mesh.nodes[mask,0] += .6
    mesh.nodes[mask,1]  -= .6

    from BasicTools.IO.XdmfWriter import XdmfWriter
    from BasicTools.Helpers.Tests import TestTempDir
    writer = XdmfWriter(fileName= TestTempDir.GetTempPath()+"FoldOverFreeMaps_CheckIntegrity.xdmf")
    print(writer.fileName)
    writer.SetTemporal(True)
    writer.SetBinary()
    writer.SetHdf5(False)
    writer.Open()
    writer.Write(mesh)

    from BasicTools.Containers.Filters import ElementFilter

    FOFM = FoldOverFreeMaps(mesh, boundaryEF=ElementFilter(mesh,dimensionality=-2) )
    FOFM.callback = writer.Write
    FOFM.Start()
    writer.Close()

    if FOFM.MinDetJacobian() < 0 :
        raise Exception("FoldOverFreeMaps unable to unfold the mesh ")

    return "ok"


if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
