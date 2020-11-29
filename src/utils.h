#pragma once

#include <string>
#include <iostream>
#include <sys/types.h>
#include <unistd.h>

using std::endl;

std::string getTimestamp();

#define _LOG getTimestamp() << " " << getpid()

#define LOG       std::cout << _LOG << " INFO  : "
#define LOG_ERR   std::cerr << _LOG << " ERROR : "
#define LOG_DEBUG std::cout << _LOG << " DEBUG : "

template<typename T>
inline T max(T one, T other) { return one > other ? one : other; }

template<typename T>
inline T max(const std::vector<T> &vec) {
    T res = T();
    for (int i = 0; i < vec.size(); ++i) {
        res = max(res, vec[i]);
    }
    return res;
}

template<typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &vec) {
    if (vec.empty()) { return out << "[]"; }
    out << "[" << vec[0];
    for (int i = 1; i < vec.size(); ++i) {
        out << ", " << vec[i];
    }
    return out << "]";
}