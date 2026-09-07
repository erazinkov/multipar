#ifndef DATA_H
#define DATA_H

#include "structs.h"

#include <map>

const std::map<std::string, ChemResult> data_chem_cat_1{
    { "sample22.sub", {14.1, 4.8 } },
    { "sample23.sub", {15.5, 4.3 } },
//    { "sample24.sub", {16.9, 5.8 } },
    { "sample25.sub", {19.7, 5.8 } },
    { "sample26.sub", {14.8, 4.4 } },
    { "sample27.sub", {16.2, 4.5 } },
    { "sample28.sub", {16.6, 4.5 } },
    { "sample29.sub", {17.0, 5.4 }, },

    { "sample120.sub", {14.9, 5.6 }, },
    { "sample121.sub", {14.2, 5.0 }, },
    { "sample122.sub", {12.9, 5.5 }, },
    { "sample123.sub", {12.9, 5.0 }, },
    { "sample124.sub", {14.8, 7.3 }, },
    { "sample125.sub", {13.8, 7.3 }, },
    { "sample126.sub", {16.0, 5.9 }, },
    { "sample127.sub", {15.4, 6.4 }, },

    { "sample137.sub", {15.0, 5.7 }, },
    { "sample138.sub", {13.7, 4.7 }, },
    { "sample139.sub", {12.3, 4.7 }, },
    { "sample140.sub", {17.4, 6.4 }, },
    { "sample141.sub", {16.6, 6.2 }, },
    { "sample142.sub", {14.2, 6.3 }, },
    { "sample143.sub", {14.3, 6.1 }, },
    { "sample144.sub", {14.8, 5.9 }, },
    { "sample145.sub", {15.5, 6.5 }, },
    { "sample146.sub", {14.1, 7.4 }, },
    { "sample147.sub", {13.0, 7.0 }, },
    { "sample148.sub", {13.0, 7.0 }, },
    { "sample149.sub", {13.0, 7.0 }, },
    { "sample150.sub", {13.0, 7.0 }, },
    { "sample151.sub", {13.0, 7.0 }, },
    { "sample152.sub", {14.0, 6.8 }, },
    { "sample153.sub", {14.0, 6.8 }, },

    { "sample168.sub", {14.7, 3.9 }, },



};

const std::map<std::string, ChemResult> data_chem_cat_3{
    { "sample30.sub", {49.9, 5.8 } },
    { "sample31.sub", {53.5, 3.9 } },
    { "sample32.sub", {55.9, 5.7 } },
//    { "sample33.sub", {59.5, 5.4 } },
    { "sample34.sub", {59.3, 5.9 } },
    { "sample35.sub", {58.6, 5.6 } },

    { "sample128.sub", {52.4, 6.6 } },
//    { "sample129.sub", {52.3, 6.3 } },
    { "sample130.sub", {55.2, 6.4 } },
    { "sample131.sub", {53.4, 6.0 } },
    { "sample132.sub", {53.4, 6.0 } },
    { "sample133.sub", {53.4, 6.0 } },
    { "sample134.sub", {53.4, 6.0 } },
    { "sample135.sub", {53.4, 6.0 } },

    { "sample154.sub", {53.2, 6.3 } },
    { "sample155.sub", {53.0, 6.2 } },


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
    { "sample58.sub", {9.5, 8.9 } },
    { "sample59.sub", {9.5,	8.9 } },
    { "sample60.sub", {9.5,	8.9 } },
    { "sample61.sub", {9.9,	8.2 } },
    { "sample62.sub", {10.0, 7.5 } },
    { "sample63.sub", {7.7,	7.8 } },
    { "sample64.sub", {8.0,	7.3 } },
    { "sample65.sub", {8.6,	9.5 } },
    { "sample66.sub", {8.5,	7.7 } },
    { "sample67.sub", {8.6,	9.8 } },
    { "sample68.sub", {8.7,	8.1 } },
    { "sample69.sub", {8.9,	9.2 } },
    { "sample72.sub", {8.2,	8.1 } },
    { "sample75.sub", {8.9,	8.1 } },
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

    // --- ! ---

//    { "sample106.sub", {9.2, 9.6 } },
//    { "sample107.sub", {9.7, 8.5 } },
//    { "sample108.sub", {8.8, 9.8 } },
//    { "sample109.sub", {7.3, 10.4 } },
//    { "sample172.sub", {6.0, 9.6 } },
//    { "sample173.sub", {6.2, 8.0 } },
//    { "sample174.sub", {5.5, 8.9 } },
//    { "sample175.sub", {5.8, 8.9 } },

    // --- ! ---
};



const std::map<std::string, ChemResult> data_chem_cat_5{
//    { "sample5.sub", {8.8, 13.9 } },
//    { "sample7.sub", {9.6, 17.8 } },

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
    { "sampleXXX.sub", {9.0, 16.8 } },
    { "sampleX.sub", {7.2, 4.9 } },
    { "sample99.sub", {8.1, 8.9 } },

    { "sample110.sub", {10.2, 8.5 } },
    { "sample111.sub", {10.0, 5.1 } },
    { "sample112.sub", {9.0, 12.1 } },
    { "sample113.sub", {7.4, 18.4 } },
    { "sample114.sub", {9.3, 9.5 } },
    { "sample115.sub", {13.0, 5.2 } },
    { "sample115_2.sub", {13.0, 5.2 } },
    { "sample115_3.sub", {13.0, 5.2 } },
    { "sample115_4.sub", {13.0, 5.2 } },
    { "sample115_5.sub", {13.0, 5.2 } },
    { "sample116.sub", {8.6, 4.9 } },
    { "sample117.sub", {9.1, 6.7 } },
    { "sample118.sub", {7.3, 4.4 } },
    { "sample119.sub", {7.1, 4.8 } },

    { "sample136.sub", {10.3, 8.0} },




};

const std::map<std::string, ChemResult> data_chem_cat_6{
//    { "sample6.sub", {70.1, 0.0} }, // 0.0
    { "sample8.sub", {73.6, 35.8 } }, // 35.8
//    { "sample11.sub", {70.9, 7.3 } },
//    { "sample12.sub", {72.3, 6.1 } },
//    { "sample13.sub", {69.2, 8.6 } },
    { "sample14.sub", {69.7, 13.6} },
//    { "sample15.sub", {78.0, 7.9 } },
    { "sample16.sub", {81.7, 7.0 } },
    { "sample17.sub", {81.6, 6.2 } },
    { "sample18.sub", {84.6, 5.1 } },
    { "sample19.sub", {83.8, 8.8 } },
    { "sample20.sub", {84.4, 6.2 } },
    { "sample21.sub", {84.1, 6.8 } },
    { "sample36.sub", {84.1, 6.5 } },

    { "sample54.sub", {71.4, 12.5 } },
    { "sample55.sub", {73.6, 34.7 } },
    { "sample56.sub", {69.3, 8.1 } },
    { "sample57.sub", {69.9, 13.1 } },

//    { "sample156.sub", {73.7, 7.5 } },
    { "sample157.sub", {74.3, 7.2 } },
    { "sample158.sub", {73.7, 8.2 } },
    { "sample159.sub", {73.7, 8.4 } },
    { "sample160.sub", {79.6, 5.4 } },
    { "sample161.sub", {75.4, 9.0 } },
    { "sample162.sub", {71.5, 10.6 } },
    { "sample163.sub", {83.9, 10.5 } },
    { "sample164.sub", {83.9, 10.5 } },
    { "sample165.sub", {83.9, 10.5 } },
    { "sample166.sub", {83.9, 10.5 } },
    { "sample167.sub", {83.9, 10.5 } },


};

// cat 1
/*
EXT PARAMETER                                   STEP         FIRST
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
   1  p0           1.00193e+00   3.80559e-02   1.74340e-05  -1.70500e-01
   2  p1           5.67864e-01   2.63587e-02   2.16807e-05   1.59611e-01
   3  p2          -4.56241e+00   5.14865e-01   3.21577e-04  -8.01615e-03
   4  p3           8.74908e+01   7.26887e+00   1.52306e-03   3.46246e-04
   5  p4           8.06289e-01   1.82220e-01   2.61688e-04  -6.00637e-03
   6  p5           9.02753e-01   8.96967e-02   2.00929e-05  -2.94535e-02
   7  p6           4.56773e-02   7.16482e-01   1.46032e-03   7.32417e-05
*/

// cat 3
/*
EXT PARAMETER                                   STEP         FIRST
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
   1  p0           9.64211e-01   9.99100e-02   9.05327e-06   1.80228e-01
   2  p1           3.58029e-01   4.21198e-02   6.02678e-06  -2.67769e-01
   3  p2          -9.70307e+00   1.92167e+00   3.26677e-04   4.87244e-03
   4  p3           1.04045e+02   2.26653e+01   6.92634e-03  -7.58536e-05
   5  p4           9.02320e-01   2.47680e+00   1.20030e-03   2.78100e-04
   6  p5           1.11815e+00   5.04027e-01   1.70588e-04   3.63423e-03
   7  p6          -8.81518e-01   2.40003e+01   1.39653e-02   3.78312e-05
*/
// cat 4
/*
EXT PARAMETER                                   STEP         FIRST
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
   1  p0           9.13914e-01   4.67104e-02   2.22216e-05   2.34596e-03
   2  p1           4.25919e-01   4.50505e-02   3.17947e-05  -2.81430e-03
   3  p2          -3.98170e+00   7.48429e-01  -1.20537e-04   2.77128e-04
   4  p3           1.18291e+02   9.54774e+00  -2.55919e-03  -6.11479e-04
   5  p4           1.20214e+00   1.55657e-01  -9.84958e-05   4.31851e-03
   6  p5           1.23644e+00   1.04238e-01  -2.24621e-05   5.02496e-02
   7  p6           1.30710e+00   4.37429e-01   1.71445e-04   2.93510e-04
*/
// cat 5
/*
EXT PARAMETER                                   STEP         FIRST
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
   1  p0           1.14853e+00   2.00570e-02   5.07044e-06   9.50869e-04
   2  p1           5.63544e-01   2.74445e-02   5.22389e-05  -7.60996e-05
   3  p2          -6.87706e+00   3.26017e-01   3.72722e-04   3.86045e-05
   4  p3           1.08000e+02   6.90281e+00   6.05738e-04  -5.45945e-05
   5  p4           9.81035e-01   6.37731e-02   1.02931e-05   6.05055e-04
   6  p5           1.13416e+00   8.10006e-02   4.57107e-06   4.47753e-03
   7  p6           1.18922e+00   3.71947e-01   1.17495e-04   3.05012e-05
*/
// cat 6
/*
EXT PARAMETER                                   STEP         FIRST
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
   1  p0           1.71637e+00   4.57256e-02   1.72171e-05   6.27334e-02
   2  p1           5.42737e-01   1.32650e-02   1.10244e-05  -9.72621e-02
   3  p2          -3.68942e+01   1.81913e+00   8.83660e-04   1.23912e-03
   4  p3           9.89012e+01   5.47881e+00   2.42155e-02   3.69320e-05
   5  p4           2.34177e-01   7.45839e-01   2.13248e-03  -6.13626e-04
   6  p5           1.18567e+00   3.74994e-01   1.50405e-03  -4.87530e-04
   7  p6           7.76461e-01   2.74187e+00   1.39526e-02   3.63202e-05
*/
#endif // DATA_H
