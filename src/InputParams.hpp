/*
MIT License
Copyright (c) 2019 - present H. Watanabe
The latest version is available at
https://github.com/kaityo256/params

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

// https://qiita.com/kaityo256/items/88a30fe1ecd9abc8449a

#pragma once
#include <fstream>
#include <iostream>
#include <map>
#include <string>

namespace param {

typedef std::map<std::string, std::string> ptype;

class InputParams {
private:
  ptype params;
  bool valid;

  bool contains(std::string &key) const {
    return params.find(key) != params.end();
  }

public:
  // constructer
  InputParams(const std::string filename)
      : valid(true) {
    loadfromfile(filename);
  }

  explicit operator bool() const {
    return valid;
  };


  void loadfromfile(const std::string filename) {
    // class for file input
    std::ifstream is(filename);
    if (is.fail()) {
      std::cerr << "Could not open file " << filename << std::endl;
      valid = false;
    }
    readfromstream(is);
  }


  void readfromstream(std::istream &is) {
    std::string line;
    while (getline(is, line)) {
      if (line.length() > 0 && line[0] == '#') {
        continue;
      }
      size_t index = line.find("=");
      if (std::string::npos != index) {
        std::string key = line.substr(0, index);
        std::string value = line.substr(index + 1, line.length());
        params.insert(ptype::value_type(key, value));
      }
    }
  }


  void check_key(std::string &key) const {
    if (!contains(key)) {
      std::cerr << "No such key: " << key << std::endl;
      std::abort();
    }
  }

  template <class Type>
  Type get(std::string) {
  }

  template <class Type>
  Type get(std::string, Type value) {
    return value;
  }

}; // class InputParams


template <>
bool InputParams::get<bool>(std::string key, bool value) {
  if (!contains(key)) {
    return value;
  }
  if ("yes" == params[key] || "Yes" == params[key]) {
    return true;
  } else {
    return false;
  }
}

template <>
bool InputParams::get(std::string key) {
  check_key(key);
  return get<bool>(key, false);
}

template <>
int InputParams::get(std::string key, int value) {
  if (!contains(key)) {
    return value;
  }
  return std::stoi(params[key]);
}

template <>
int InputParams::get(std::string key) {
  check_key(key);
  return get<int>(key, 0);
}

template <>
double InputParams::get(std::string key, double value) {
  if (!contains(key)) {
    return value;
  }
  return std::stod(params[key]);
}

template <>
double InputParams::get(std::string key) {
  check_key(key);
  return get<double>(key, 0.0);
}

} // namespace param