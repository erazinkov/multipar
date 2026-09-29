#ifndef CALCULATIONSAW_H
#define CALCULATIONSAW_H


#include <string>
#include <map>
#include <vector>

class CalculationsAW {
public:
    /**
     * Вычисляет значение A или W для заданного типа.
     *
     * @param concentrations  Карта соответствия символов элементов (например, "O", "C", "N") их значениям.
     * @param value           Какое значение вычислить: "A" или "W".
     * @param type            Ключ набора параметров (например, "cat1"). По умолчанию пустой.
     * @return                Вычисленное значение или 0.0, если тип/значение неизвестно или некорректно.
     */
    static double calculateValueByType(const std::map<std::string, double> &concentrations,
                                       const std::string &value,
                                       const std::string &type = "");

private:
    CalculationsAW() = delete;
    static double getConcentrationByElement(const std::map<std::string, double> &concentrations,
                                            const std::string &element);
    inline static const std::map<std::string, std::vector<double>> typesParameters_{
        { "cat1",
            {
                 8.60112e-01,
                 4.25135e-01,
                -4.10488e+00,
                 9.17170e+01,
                 9.55916e-01,
                 9.47465e-01,
                 5.32751e-02,
            }
        },
        { "cat3",
            {
                4.14040e-01,
                2.10373e-01,
                2.42752e+00,
               -2.36117e+02,
               -1.69721e+02,
                9.70276e+00,
                6.89231e+02,
            }
        },
        { "cat40_43",
            {
                 1.10805e+00,
                 4.42134e-01,
                -7.24584e+00,
                 9.29053e+01,
                 7.87820e-01,
                 9.70413e-01,
                 8.19574e-01,
            }
        },
        { "cat41",
            {
                 9.75127e-01,
                 3.31221e-01,
                -6.69865e+00,
                 8.44098e+01,
                 6.66039e-01,
                 8.78435e-01,
                 1.30395e+00,
            }
        },
        { "cat42",
            {
                 1.12625e+00,
                 4.68697e-01,
                -7.10413e+00,
                 9.21494e+01,
                 6.84809e-01,
                 9.69277e-01,
                 6.58060e-01,
            }
        },
        { "cat44",
            {
                 9.85265e-01,
                 4.60579e-01,
                -4.92370e+00,
                 7.34227e+01,
                 6.65926e-01,
                 7.52364e-01,
                 5.21511e-01,
            }
        },
        { "cat53",
            {
                 9.06990e-01,
                 3.76831e-01,
                -4.00366e+00,
                 1.06400e+02,
                 1.04092e+00,
                 1.11785e+00,
                 8.27098e-01,
            }
        },
    };
};

#endif // CALCULATIONSAW_H
