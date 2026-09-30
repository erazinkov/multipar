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
        { "cat44", // cat44_45
            {
                 1.08732e+00,
                 5.44072e-01,
                -6.13573e+00,
                 1.06681e+02,
                 9.85254e-01,
                 1.11927e+00,
                 1.20606e+00,
            }
        },
        { "cat51",
            {
                5.06248e-01,
                2.24670e-01,
                6.40886e+00,
                1.17249e+02,
                1.44386e+00,
                1.15433e+00,
                8.41011e-01,
            }
        },
        { "cat52",
            {
                 7.76497e-01,
                 3.35370e-01,
                -3.67203e+00,
                 1.08696e+02,
                 1.26253e+00,
                 1.12612e+00,
                 1.38301e+00,
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
        { "cat54",
            {
                 1.18368e+00,
                 4.44667e-01,
                -8.82825e+00,
                 1.19357e+02,
                 9.61402e-01,
                 1.28614e+00,
                 9.05662e-01,
            }
        },
        { "cat6",
            {
                 1.71703e+00,
                 5.31191e-01,
                -3.78993e+01,
                 9.89324e+01,
                 2.24912e-01,
                 1.19709e+00,
                 8.01820e-01,
            }
        },
    };
};

#endif // CALCULATIONSAW_H
