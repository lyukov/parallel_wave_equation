#pragma once

#include <string>
#include <iostream>

using std::endl;

std::string getTimestamp();

#define LOG std::cout << getTimestamp() << " : "
#define LOG_ERR std::cerr << getTimestamp() << " : "

template <typename T>
inline T max(T one, T other) { return one > other ? one : other; }