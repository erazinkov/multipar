#ifndef STRUCTS_H
#define STRUCTS_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <optional>

struct ElementResult {
    std::string name;
    double value;
    double valueError;
    void print() const {
        std::cout << name << " " << value << "\u00B1" << valueError << " ";
    }
};

struct FitResult {
    std::vector<ElementResult> elementResults;
    void print() const {
        for (const auto &er : elementResults) {
            er.print();
        }
    }
    const ElementResult& getElementResultByName(const std::string& name) const {
        auto it = std::find_if(
            elementResults.begin(),
            elementResults.end(),
            [&](const ElementResult& er) { return er.name == name; }
            );
        if (it == elementResults.end()) {
            throw std::runtime_error("Element '" + name + "' not found");
        }
        return *it;
    }
};

struct ChemResult {
    enum class Type {
        A,
        W
    };
    std::optional<double> a;
    std::optional<double> w;
    void print() const {
        auto a{this->a.has_value() ?  this->a.value() : -1};
        auto w{this->w.has_value() ?  this->w.value() : -1};
        std::cout << a << " " << w;
    }
};

struct Data {
    ChemResult chemResult;
    std::vector<FitResult> fitResults;
    void print() const {
        this->chemResult.print();
        std::cout << " ";
        for (const auto& fr : this->fitResults) {
            fr.print();
        }
    }
};

struct Point  {
    std::string sample;
    ChemResult chemResult;
    FitResult fitResult;
    double x;
    double xErr;
    double y;
    double yErr;
};

#endif // STRUCTS_H
