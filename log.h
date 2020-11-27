#pragma once

#include <string>
#include <iostream>

std::string getTimestamp();

#define LOG std::cout << getTimestamp() << " : "
#define LOG_ERR std::cerr << getTimestamp() << " : "
