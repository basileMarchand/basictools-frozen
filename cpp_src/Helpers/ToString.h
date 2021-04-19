//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
#pragma once
namespace BasicTools
{

template<typename T>
std::string ToString(const T& obj) {
    std::ostringstream  oss;
    oss << obj ;
    return oss.str();
};

} // BasicTools namespace

#ifdef NDEBUG
#define PRINTDEBUG(args)
#define PRINTDEBUGArgs(args)
#else
#define PRINTDEBUG(args) std::cout << "file: " << __FILE__ << "::" << __LINE__<< " " << " -> " << args << std::endl;
#define PRINTDEBUGArgs(args) std::cout << "file: " << __FILE__ << "::" << __LINE__<< " " << #args << " -> " << args << std::endl;
#endif
