//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
#pragma once

#include <Containers/UnstructuredMesh.h>

namespace BasicTools
{
    
MatrixDDD GetElementsCenters(const MapMatrixDDD& nodes, const ElementsContainer& elements);

} // namespace BasicTools
