#pragma once
#include<exception>

namespace gftools {

/** Base class for all gftools exceptions. */
class gftools_exception : public std::exception {};

/** Generic error - similar to std::logic_error, but inherited from gftools_exception. */
class ex_generic : public gftools_exception { 
public: 
    ex_generic(std::string s):s_(s){}; 
    virtual const char* what() const throw() { return s_.c_str(); } 
protected:
    std::string s_; 
};

}

