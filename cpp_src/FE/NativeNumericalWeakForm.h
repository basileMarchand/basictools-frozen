//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
#pragma once

#include <string>
#include <ostream>
#include <vector>

#include <LinAlg/EigenTypes.h>


namespace BasicTools
{

struct WeakTerm {
    WeakTerm ():fieldName("None"),derCoordName("None"),derDegree(0),constant(false),normal(false){};
    std::string fieldName;
    std::string derCoordName;
    int derDegree;
    bool constant;
    bool normal;

    int spaceIndex_;
    int derCoordIndex_;
    int numberingIndex_;
    int valuesIndex_;
    int modeIndex_;
    int internalType;
    friend std::ostream& operator<< (std::ostream& stream, const WeakTerm& term);

};

struct WeakMonom {
    WeakMonom():prefactor(1){};
    double prefactor;
    std::vector<WeakTerm> prod;
    friend std::ostream& operator<< (std::ostream& stream, const WeakMonom& monom);

};

struct WeakForm{
    std::vector<WeakMonom> form;
    CBasicIndexType GetNumberOfTerms(){return static_cast<CBasicIndexType>(this->form.size());}
};

std::ostream& operator <<(std::ostream& stream, const WeakTerm& term);

std::ostream& operator <<(std::ostream& stream, const WeakMonom& monom) ;

} // namespace BasicTools
