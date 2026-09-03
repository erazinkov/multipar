#ifndef DATA_H
#define DATA_H

#include "structs.h"

#include <map>

const std::map<std::string, ChemResult> data_chem_cat_1{
    { "sample22.sub", {14.1, 4.8 } },
    { "sample23.sub", {15.5, 4.3 } },
    { "sample24.sub", {16.9, 5.8 } },
    { "sample25.sub", {19.7, 5.8 } },
    { "sample26.sub", {14.8, 4.4 } },
    { "sample27.sub", {16.2, 4.5 } },
    { "sample28.sub", {16.6, 4.5 } },
    { "sample29.sub", {17.0, 5.4}, },
};

const std::map<std::string, ChemResult> data_chem_cat_3{
    { "sample30.sub", {49.9, 5.8 } },
    { "sample31.sub", {53.5, 3.9 } },
    { "sample32.sub", {55.9, 5.7 } },
    { "sample33.sub", {59.5, 5.4 } },
    { "sample34.sub", {59.3, 5.9 } },
    { "sample35.sub", {58.6, 5.6 } },
};

const std::map<std::string, ChemResult> data_chem_cat_4{
    { "sample1.sub", {5.7,  8.5 } },
    { "sample2.sub", {5.8,  9.2 } },
    { "sample3.sub", {8.1,  7.6 } },
    { "sample4.sub", {8.9,  8.4 } },
    { "sample9.sub",  {9.3,  9.4 } },
    { "sample10.sub", {9.2, 10.0 } },
    { "sample37.sub", {8.0, 9.2 } },
    { "sample38.sub", {7.9, 8.0 } },
    { "sample39.sub", {7.7, 8.2 } },
    { "sample40.sub", {9.0, 8.6 } },
    { "sample41.sub", {8.3, 8.3 } },
    { "sample42.sub", {8.7, 8.9 } },
    { "sample43.sub", {9.7, 8.7 } },
    { "sample44.sub", {9.3, 10.0 } },
    { "sample45.sub", {9.4, 9.9 } },
    { "sample46.sub", {9.6, 8.1 } },
    { "sample47.sub", {9.6, 8.8 } },
    { "sample48.sub", {9.6, 7.8 } },
    { "sample49.sub", {10.3,9.4 } },
    { "sample50.sub", {9.3, 9.7 } },
    { "sample51.sub", {9.5, 9.7 } },
    { "sample52.sub", {9.4, 9.1 } },
    { "sample53.sub", {9.4, 8.7 } },
    { "sample76.sub", {9.1,	9.7 } },
    { "sample77.sub", {8.9,	7.0 } },
    { "sample78.sub", {8.6,	9.5 } },
    { "sample79.sub", {9.4,	7.6 } },
    { "sample80.sub", {8.0,	9.2 } },
    { "sample81.sub", {7.3,	9.3 } },
    { "sample82.sub", {7.6,	10.4 } },
    { "sample83.sub", {7.6,	10.7 } },
    { "sample84.sub", {8.3,	10.9 } },
    { "sample85.sub", {8.6,	10.0 } },
};

const std::map<std::string, ChemResult> data_chem_cat_4_grad{
    { "sample1.sub", {5.7,  8.5 } },
    { "sample2.sub", {5.8,  9.2 } },
    { "sample3.sub", {8.1,  7.6 } },
    { "sample4.sub", {8.9,  8.4 } },
    { "sample9.sub",  {9.3,  9.4 } },
    { "sample10.sub", {9.2, 10.0 } },
    { "sample37.sub", {8.0, 9.2 } },
    { "sample38.sub", {7.9, 8.0 } },
    { "sample39.sub", {7.7, 8.2 } },
    { "sample40.sub", {9.0, 8.6 } },
    { "sample41.sub", {8.3, 8.3 } },
    { "sample42.sub", {8.7, 8.9 } },
    { "sample43.sub", {9.7, 8.7 } },
    { "sample44.sub", {9.3, 10.0 } },
    { "sample45.sub", {9.4, 9.9 } },
    { "sample46.sub", {9.6, 8.1 } },
    { "sample47.sub", {9.6, 8.8 } },
};

const std::map<std::string, ChemResult> data_chem_cat_4_check{
    { "sample48.sub", {9.6, 7.8 } },
    { "sample49.sub", {10.3,9.4 } },
    { "sample50.sub", {9.3, 9.7 } },
    { "sample51.sub", {9.5, 9.7 } },
    { "sample52.sub", {9.4, 9.1 } },
    { "sample53.sub", {9.4, 8.7 } },
    { "sample76.sub", {9.1,	9.7 } },
    { "sample77.sub", {8.9,	7.0 } },
    { "sample78.sub", {8.6,	9.5 } },
    { "sample79.sub", {9.4,	7.6 } },
    { "sample80.sub", {8.0,	9.2 } },
    { "sample81.sub", {7.3,	9.3 } },
    { "sample82.sub", {7.6,	10.4 } },
    { "sample83.sub", {7.6,	10.7 } },
    { "sample84.sub", {8.3,	10.9 } },
    { "sample85.sub", {8.6,	10.0 } },
};


const std::map<std::string, ChemResult> data_chem_cat_5{
    { "sample5.sub", {8.8, 13.9 } },
    { "sample7.sub", {9.6, 17.8 } },

    { "sample86.sub", {9.6,11.8 } },
    { "sample87.sub", {8.7, 4.5 } },
    { "sample88.sub", {9.5, 7.3 } },
    { "sample89.sub", {6.6, 17.0 } },
    { "sample90.sub", {7.9, 5.8 } },
    { "sample91.sub", {8.6, 5.5 } },
    { "sample92.sub", {7.0, 11.6 } },
    { "sample93.sub", {6.4, 16.7 } },
    { "sample94.sub", {6.9, 4.7 } },
    { "sample95.sub", {8.3, 7.5 } },
    { "sample96.sub", {6.5, 9.4 } },
    { "sampleX.sub", {9.0, 16.8 } },
    { "sampleXXX.sub", {7.2, 4.9 } },
    { "sample99.sub", {8.1, 8.9 } },

};

const std::map<std::string, ChemResult> data_chem_cat_6{
//    { "sample6.sub", {70.1, 0.0} }, // 0.0
    { "sample8.sub", {73.6, 35.8 } }, // 35.8
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
    { "sample36.sub", {84.1, 6.5 } },
};

#endif // DATA_H
