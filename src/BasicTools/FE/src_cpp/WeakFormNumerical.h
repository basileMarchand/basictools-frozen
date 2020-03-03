//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//


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
    int GetNumberOfTerms(){return this->form.size();}
};
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
