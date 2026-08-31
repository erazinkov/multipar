#ifndef DATA_H
#define DATA_H

#include "structs.h"

#include <map>

const std::map<std::string, ChemResult> data_chem{
    { "sample1.sub", {5.7,  8.5 } },
    { "sample2.sub", {5.8,  9.2 } },
    { "sample3.sub", {8.1,  7.6 } },
    { "sample4.sub", {8.9,  8.4 } },
    { "sample5.sub", {8.8, 13.9 } },
    // { "sample6.sub", {70.1, 0.0 } },
    { "sample7.sub", {9.6, 17.8 } },
    // { "sample8.sub", {73.6,35.8 } },
    { "sample9.sub",  {9.3,  9.4 } },
    { "sample10.sub", {9.2, 10.0 } },
    { "sample11.sub", {70.9, 7.3 } },
    { "sample12.sub", {72.3, 6.1 } },
    { "sample13.sub", {69.2, 8.6 } },
    { "sample14.sub", {69.7, 13.6} },
    { "sample15.sub", {78.0, 7.9 } },
    { "sample16.sub", {81.7, 7.0 } },
    { "sample17.sub", {81.6, 6.2 } },
    { "sample18.sub", {84.6, 5.1 } },
    { "sample19.sub", {83.8, 8.8 } },
    { "sample20.sub", {84.4, 6.2 } },
    { "sample21.sub", {84.1, 6.8 } },
    { "sample22.sub", {14.1, 4.8 } },
    { "sample23.sub", {15.5, 4.3 } },
    { "sample24.sub", {16.9, 5.8 } },
    { "sample25.sub", {19.7, 5.8 } },
    { "sample26.sub", {14.8, 4.4 } },
    { "sample27.sub", {16.2, 4.5 } },
    { "sample28.sub", {16.6, 4.5 } },
    { "sample29.sub", {17.0, 5.4}, },
    { "sample30.sub", {49.9, 5.8 } },
};

#endif // DATA_H
