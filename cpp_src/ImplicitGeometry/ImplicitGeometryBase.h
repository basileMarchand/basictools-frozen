//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//

# pragma once

#include <LinAlg/EigenTypes.h>
#include <Containers/UnstructuredMesh.h>

namespace BasicTools
{

class ImplicitGeometryBase{
    bool insideOut;
    void ApplyInsideOut(MatrixDD1& res) const;
    virtual MatrixDD1 GetDistanceToPointsInternal(const MatrixDDD& pos) const = 0;
public:
    ImplicitGeometryBase();
    //
    virtual MatrixDD1 GetDistanceToPoints(const MatrixDDD& pos) const;
    MatrixDD1 GetDistanceToMesh(UnstructuredMesh& mesh,bool cellCenter=false ) const;
    //
    void SetInsideOutOn();
    void SetInsideOutOff();
    void SetInsideOut(const bool& val);
    bool GetInsideOut()const;
};

} // namespace BasicTools
