#include "calculationsaw.h"


double CalculationsAW::calculateValueByType(const std::map<std::string, double>& concentrations,
                                            const std::string& value,
                                            const std::string& type) {
    const auto& table = typesParameters_;
    auto it = table.find(type);
    if (it == table.end() || it->second.empty()) {
        return 0.0;
    }

    const std::vector<double>& p = it->second;
    if (p.size() < 7) {
        return 0.0;
    }

    const double O = getConcentrationByElement(concentrations, "O");
    const double C = getConcentrationByElement(concentrations, "C");
    const double N = getConcentrationByElement(concentrations, "N");

    if (1.0 - p.at(1) * p.at(4) == 0.0) {
        return 0.0;
    }

    const double A = (
        p.at(3)
        - p.at(4) * p.at(0) * O
        - p.at(2) * p.at(4)
        - p.at(5) * C
        - p.at(6) * N
    ) / (1.0 - p.at(1) * p.at(4));

    const double W = (
        p.at(0) * O
        - p.at(1) * A
        + p.at(2)
    );

    if (value == "A") return A;
    if (value == "W") return W;
    return 0.0;
}

double CalculationsAW::getConcentrationByElement(const std::map<std::string, double> &concentrations, const std::string &element)
{
    auto it = concentrations.find(element);
    return (it != concentrations.end()) ? it->second : 0.0;
}
