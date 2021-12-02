//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//

#include "NativeIntegration.h"
#include <Eigen/SparseCore>

#include <iostream>
#include <string>

namespace BasicTools
{

std::ostream& operator <<(std::ostream& stream, const WeakTerm& term) {

if (term.derDegree){
    stream << "Derivative(" << term.fieldName <<"," << term.derCoordName<< ")"  ;
} else {
    stream << term.fieldName ;
}
return stream;
}

std::ostream& operator <<(std::ostream& stream, const WeakMonom& monom) {

    if (monom.prefactor != 0) {
        stream << monom.prefactor << "*";
    }
    for(unsigned int prodn=0; prodn<  monom.prod.size(); ++prodn){
        const WeakTerm& term = monom.prod[prodn];
        if(prodn ) stream << "*";
        stream << term;

    }

    return stream;
}
} // namespace BasicTools
