//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
#pragma once

#include <string>
#include <ostream>
#include <vector>

#include <LinAlg/EigenTypes.h>

namespace BasicTools {

struct WeakTerm {
    std::string fieldName{"None"};
    std::string derCoordName{"None"};
    int derDegree{0};
    bool constant{false};
    bool normal{false};

    int spaceIndex_;
    int derCoordIndex_;
    int numberingIndex_;
    int valuesIndex_;
    int modeIndex_;
    int internalType;
};

std::ostream& operator<<(std::ostream& stream, const WeakTerm& term);

struct WeakMonom {
    double prefactor{1.0};
    std::vector<WeakTerm> prod;
};

std::ostream& operator<<(std::ostream& stream, const WeakMonom& monom);

struct WeakForm {
    std::vector<WeakMonom> form;
    CBasicIndexType GetNumberOfTerms() const { return static_cast<CBasicIndexType>(this->form.size()); }
};

std::ostream& operator<<(std::ostream& stream, const WeakForm& WeakForm);
}  // namespace BasicTools
