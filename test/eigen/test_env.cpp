// https://qiita.com/StudyKogaku/items/ce22f1b7e9b9a812debf
// plz define ${CPATH} at .zshrc

#include <stdio.h>
#include <iostream>
#include <boost/version.hpp>
#include <eigen3/Eigen/Core>

int main() {
    printf("Hello World\n");

    std::cout << "boostバージョン:" << BOOST_VERSION << std::endl;
    std::cout << "boostライブラリバージョン:" << BOOST_LIB_VERSION << std::endl;
    return 0;
}
