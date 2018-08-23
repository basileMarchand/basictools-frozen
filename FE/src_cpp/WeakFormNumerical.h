

struct WeakTerm {
    WeakTerm ():fieldName("Nonne"),derCoordName("None"),derDegree(0),constant(false),normal(false){};
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
};

struct WeakMonom {
    WeakMonom():prefactor(1){};
    double prefactor;
    std::vector<WeakTerm> prod;
};

struct WeakForm{
    std::vector<WeakMonom> form;
    int GetNumberOfTerms(){return this->form.size();}
};