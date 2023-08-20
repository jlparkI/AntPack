#include "aligners.h"

IMGTAligner::IMGTAligner(
                 py::array_t<double> scoreArray
)
{
    py::buffer_info info = scoreArray.request();
    if (info.shape.size() != 2){
        throw std::runtime_error(std::string("The scoreArray passed to IMGTAligner must be a 2d array"));
    }
    if (info.shape[0] < 1){
        throw std::runtime_error(std::string("The scoreArray passed to IMGTAligner must have >= 1 row"));
    }
    if (info.shape[1] != this->numAAs){
        throw std::runtime_error(std::string("The scoreArray passed to IMGTAligner must have 21 columns (1 per AA)"));
    }
    numPositions = info.shape[0];
    scoreArray = scoreArray;
}

std::vector<std::string> IMGTAligner::align(std::string query_sequence){
    std::vector<std::string> numbering;
    int numElements = query_sequence.length() * this->numPositions;
    double *needleScores = new double[ numElements ];
    delete[] needleScores;
    return numbering;
}