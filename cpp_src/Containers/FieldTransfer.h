
#include <string>

template<typename T>
double normsquared(const T& v){
    return v.dot(v);
};

class TransferClass {
    std::string ToStr(){ return "TransferClass"; }
};


template<typename A, typename B, typename C, typename D, typename S1, typename S2, typename S3>
void  ComputeBarycentricCoordinateOnElement(const A& coordAtDofs, const B& localspace, const C& targetPoint, const D& elementType,
                                            S1& inside, S2& xietaphi, S3& xichietaClamped){


const int elemDim = localspace[this->geoSpaceNumber].dimensionality;
const int spaceDim = static_cast<int>(this->nodes->cols());


/*
    linear = ElementNames.linear[elementType]
    spacedim = localspace.GetDimensionality()

    xietaphi = np.array([0.5]*spacedim)
    N = localspace.GetShapeFunc(xietaphi)
    currentPoint = N.dot(coordAtDofs)
    f = targetPoint - currentPoint

    for x in range(10):
        dN = localspace.GetShapeFuncDer(xietaphi)
        df_num = df(f,dN,coordAtDofs)
        H = ddf(f, xietaphi, dN, localspace.GetShapeFuncDerDer, coordAtDofs, linear)
        if spacedim == 2:
            dxietaphi = inv22(H).dot(df_num)
        elif spacedim == 3:
            dxietaphi = hdinv(H).dot(df_num)
        else:
            dxietaphi = df_num/H[0,0]
        xietaphi -= dxietaphi

        # if the cell is linear only one iteration is needed
        if linear :
            break

        N = localspace.GetShapeFunc(xietaphi)
        f = targetPoint - N.dot(coordAtDofs)

        if normsquared(dxietaphi) < 1e-3 and normsquared(f) < 1e-3 :
            break
    else:
        return None, xietaphi,localspace.ClampParamCoorninates(xietaphi)

    xichietaClamped = localspace.ClampParamCoorninates(xietaphi)
    # we treat very closes point as inside
    inside = normsquared(xichietaClamped-xietaphi) < 1e-5
    return inside, xietaphi, xichietaClamped
*/

}


