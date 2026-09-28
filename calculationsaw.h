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
                 9.93362e-01,
                 3.67507e-01,
                -1.03398e+01,
                -7.64244e+01,
                -9.54351e+01,
                 5.65644e+00,
                 3.66945e+02,
            }
        },
        { "cat40_43",
            {
                1.95006e-02,
                3.77207e-02,
                3.72222e-01,
                7.69005e+00,
                9.15717e-02,
                8.73175e-02,
                2.52393e-01,
            }
        },
        { "cat41", // +
            {
                 9.96904e-01,
                 3.56129e-01,
                -6.87448e+00,
                 7.90066e+01,
                 6.15664e-01,
                 8.14143e-01,
                 1.27773e+00,
            }
        },
        { "cat42", // +
            {
                 1.09703e+00,
                 4.60467e-01,
                -6.64603e+00,
                 9.06178e+01,
                 6.92772e-01,
                 9.49528e-01,
                 5.99914e-01,
            }
        },
        { "cat44", // +
            {
                 9.71127e-01,
                 3.84690e-01,
                -5.14985e+00,
                 1.07482e+02,
                 1.10473e+00,
                 1.11884e+00,
                 1.03563e+00,
            }
        },
    };
};

#endif // CALCULATIONSAW_H
