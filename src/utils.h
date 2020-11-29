#pragma once

#include <string>
#include <iostream>
#include <sys/types.h>
#include <unistd.h>

using std::endl;

std::string getTimestamp();

#define LOG std::cout << getTimestamp() << " " << getpid() << " : "
#define LOG_ERR std::cerr << getTimestamp() << " " << getpid() << " : "

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