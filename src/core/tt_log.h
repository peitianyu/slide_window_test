#ifndef __LOG_H__
#define __LOG_H__

#include<iostream>

#define NO_PRINT_DEFINE

template <typename... Args>
std::ostream & tt_print(std::ostream &os, Args&&... args) {
    #ifdef NO_PRINT_DEFINE
    #else
    auto a = {(os << std::forward<Args>(args), 0)...};
    (void)a;
    #endif
    return os;
}

#define LOG(...)  tt_print(std::cout, __VA_ARGS__)
#define LOG_FILE(ofs, ...)  tt_print(ofs, __VA_ARGS__)

#endif // __LOG_H__