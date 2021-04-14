//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
#pragma once

template<typename T>
std::string ToString(const T& obj){
    std::ostringstream  oss;
    oss << obj ;
    return oss.str();
};

#ifdef NDEBUG
    #define PRINT(args)
    #define PRINTArgs(args)
#else
    #define PRINT(args) std::cout << "file: " << __FILE__ << "::" << __LINE__<< " " << " -> " << args << std::endl;
    #define PRINTArgs(args) std::cout << "file: " << __FILE__ << "::" << __LINE__<< " " << #args << " -> " << args << std::endl;
#endif
