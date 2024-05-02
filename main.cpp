#include <iostream>
#include <vector>
#include <array>
#include "algebra.hpp"
#include <string>
#include <chrono.hpp>
#include "json.hpp"
using json = nlohmann::json;

int main(int argc, char *argv[]) 
{    
    std::string nomeFile;

    if (argc < 2) {
        std::ifstream f("matrix.json");
        json mtx = json::parse(f);
        nomeFile = mtx.value("file_mtx", "");
    } else {
        nomeFile = argv[1];
    }

    algebra::Matrix<double, algebra::order::row> Mat(nomeFile);

    std::vector<double> vec(Mat.cols(), 1);

    std::vector<double> res(Mat.rows(), 0);
    Timings::Chrono     clock1;
    clock1.start();
    
    res = Mat * vec;
  
    clock1.stop();
    
    std::cout << "The multiplication using the dynamic formatting of the matrix takes: " << clock1 << std::endl;

    Mat.compress();

    clock1.start();
    
    res = Mat * vec;

    clock1.stop();
   
   
    std::cout << "The multiplication using the compressed formatting of the matrix takes: " << clock1 << std::endl;


    return 0;
}
