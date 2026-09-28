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
    };
};

#endif // CALCULATIONSAW_H
