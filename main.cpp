#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TF1.h>
#include <optional>
#include <exception>
#include <TVirtualFitter.h>
#include <TLatex.h>
#include <numeric>
#include <TH2.h>
#include <TLine.h>
#include <TMarker.h>
#include <TStyle.h>

#include <bits/stdc++.h>

#include <regex>

class my_error: public std::exception
{
public:
    my_error(const std::string& message) : _message{message}
    {}
    const char* what() const noexcept override
    {
        return _message.c_str();
    }
private:
    std::string _message;
};


struct FitResult {
    std::string e;
    std::optional<double> value;
    std::optional<double> valueError;
    void print() const {
        auto v{this->value.has_value() ?  this->value.value() : -1};
        auto vE{this->valueError.has_value() ?  this->valueError.value() : -1};
        std::cout << e << " " << v << "\u00B1" << vE << " ";
    }
};

struct ChemResult {
    std::string e;
    std::optional<double> value;
    std::optional<double> valueError;
    void print() const {
        auto v{this->value.has_value() ?  this->value.value() : -1};
        auto vE{this->valueError.has_value() ?  this->valueError.value() : -1};
        std::cout << e << " " << v << "\u00B1" << vE << " ";
    }
};


struct Data {
    std::vector<ChemResult> cr;
    std::vector<std::vector<FitResult>> fr;
    void print() const {
        std::cout << "[ ";
        for (const auto& eItem : this->cr)
        {
             eItem.print();
        }
        std::cout << " ] ";
        for (const auto& eItem : this->fr)
        {
            for (const auto& eItemItem : eItem)
            {
                eItemItem.print();
            }
        }
    }
};

struct Rep {
    std::string e;
    double avg;
    double stdAbs;
    double stdRel;
};

struct GradPar {
    double slope;
    double offset;
};

struct DataGrad {
    std::string e;
    std::vector<double> x;
    std::vector<double> xErr;
    std::vector<double> y;
    std::vector<double> yErr;
};

std::vector<std::string> splitLineToStrs(const std::string &line);

double strToDouble(std::string str);

std::map<std::string, Data> getFitResults(const std::string &fileName,
                   const std::map<int, std::string> &columnElement,
                   const std::map<std::string, std::vector<ChemResult>> &chem,
                   const std::regex &pattern);

void calcRepByElement(const std::map<std::string, Data> &data, const std::string &element, std::map<std::string, std::vector<Rep>> &m, GradPar gP);
void updateDataGradByElement(const std::map<std::string, Data> &data, const std::string &element, DataGrad &dG);

int main()
{

    std::map<std::string, std::vector<ChemResult>> chem1
    {
        {"11979", {{"Zn", 0.77, 0.0}, {"Pb", 0.18, 0.0}, {"Ag", 8.50, 0.0}, {"Cd", 0.0017, 0.0}, {"Fe", 10.1, 0.0}, {"S", 4.570, 0.0}, {"Сu", 0.170, 0.0}}},
        {"12701", {{"Zn", 1.60, 0.0}, {"Pb", 0.39, 0.0}, {"Ag", 20.4, 0.0}, {"Cd", 0.0045, 0.0}, {"Fe", 21.7, 0.0}, {"S", 8.060, 0.0}, {"Сu", 0.044, 0.0}}},
        {"12704", {{"Zn", 2.90, 0.0}, {"Pb", 0.57, 0.0}, {"Ag", 18.7, 0.0}, {"Cd", 0.0063, 0.0}, {"Fe", 9.00, 0.0}, {"S", 10.14, 0.0}, {"Сu", 0.048, 0.0}}},
        {"12708", {{"Zn", 3.70, 0.0}, {"Pb", 0.73, 0.0}, {"Ag", 20.6, 0.0}, {"Cd", 0.0073, 0.0}, {"Fe", 9.30, 0.0}, {"S", 11.18, 0.0}, {"Сu", 0.029, 0.0}}},
        {"12677", {{"Zn", 4.50, 0.0}, {"Pb", 1.10, 0.0}, {"Ag", 39.8, 0.0}, {"Cd", 0.0090, 0.0}, {"Fe", 16.3, 0.0}, {"S", 15.21, 0.0}, {"Сu", 0.042, 0.0}}},
        {"12673", {{"Zn", 5.60, 0.0}, {"Pb", 1.20, 0.0}, {"Ag", 46.1, 0.0}, {"Cd", 0.0120, 0.0}, {"Fe", 21.2, 0.0}, {"S", 19.46, 0.0}, {"Сu", 0.031, 0.0}}},
        {"12648", {{"Zn", 6.10, 0.0}, {"Pb", 0.95, 0.0}, {"Ag", 46.3, 0.0}, {"Cd", 0.0130, 0.0}, {"Fe", 37.3, 0.0}, {"S", 33.45, 0.0}, {"Сu", 0.029, 0.0}}},
        {"12518", {{"Zn", 7.20, 0.0}, {"Pb", 1.30, 0.0}, {"Ag", 41.5, 0.0}, {"Cd", 0.0130, 0.0}, {"Fe", 35.0, 0.0}, {"S", 39.81, 0.0}, {"Сu", 0.022, 0.0}}},
        {"12532", {{"Zn", 8.90, 0.0}, {"Pb", 1.50, 0.0}, {"Ag", 27.3, 0.0}, {"Cd", 0.0140, 0.0}, {"Fe", 29.2, 0.0}, {"S", 26.06, 0.0}, {"Сu", 0.010, 0.0}}},
        {"12517", {{"Zn", 13.3, 0.0}, {"Pb", 1.80, 0.0}, {"Ag", 45.8, 0.0}, {"Cd", 0.0200, 0.0}, {"Fe", 23.6, 0.0}, {"S", 25.63, 0.0}, {"Сu", 0.012, 0.0}}},
        {"11591", {{"Zn", 16.5, 0.0}, {"Pb", 3.60, 0.0}, {"Ag", 97.9, 0.0}, {"Cd", 0.0350, 0.0}, {"Fe", 31.1, 0.0}, {"S", 37.39, 0.0}, {"Сu", 0.066, 0.0}}},
        {"Z93"  , {{"Zn", 21.2, 0.0}, {"Pb", 3.30, 0.0}, {"Ag", 54.3, 0.0}, {"Cd", 0.0350, 0.0}, {"Fe", 15.8, 0.0}, {"S", 22.43, 0.0}, {"Сu", 0.010, 0.0}}}
    };

    const std::map<int, std::string> columnElement1
    {
        {1, "Al"},
        {3, "C"},
        {5, "Ca"},
        {7, "Fe"},
        {9, "K"},
        {11, "Mg"},
        {13, "O"},
        {15, "Pb"},
        {17, "S"},
        {19, "Si"},
        {21, "Zn"},
    };

    std::string fileName{"rea.elts.txt"};
//    std::regex mGrad{"^(?!.*povtor).*$"};
    std::regex mGrad{"ozer"};

    auto dataGrad{getFitResults(fileName, columnElement1, chem1, mGrad)};

    std::vector<DataGrad> dG;

    for (const auto &item : columnElement1)
    {
        DataGrad dGtmp;
        updateDataGradByElement(dataGrad, item.second, dGtmp);
        if (!dGtmp.x.empty())
        {
            dG.push_back(dGtmp);
        }
    }

    std::map<std::string, GradPar> gP;

    for (const auto &item : columnElement1)
    {
        gP[item.second] = {1.0, 0.0};
    }

    const std::string psName{"output.ps"};
    std::unique_ptr<TCanvas> c{new TCanvas("c", "c", 1024, 960)};
    auto cd{static_cast<int>(std::ceil(std::sqrt(dG.size())))};
    c.get()->Divide(cd ,cd);
    c.get()->Print((psName + '[').c_str());
    for (size_t i{0}; i < dG.size(); ++i)
    {
        c.get()->cd(static_cast<int>(i + 1));
        gStyle->SetOptFit(1111);
        std::unique_ptr<TGraphErrors> gr{new TGraphErrors(dG.at(i).x.size(), &dG.at(i).x[0], &dG.at(i).y[0], &dG.at(i).xErr[0], &dG.at(i).yErr[0])};
        gr.get()->SetTitle(dG.at(i).e.c_str());
        std::unique_ptr<TF1> f{new TF1("f", "pol1(0)", *std::min_element(dG.at(i).x.cbegin(), dG.at(i).x.cend()), *std::max_element(dG.at(i).x.cbegin(), dG.at(i).x.cend()))};
        gr.get()->Fit(f.get(), "R");
        gr.get()->DrawClone("AP");
        gP.at(dG.at(i).e).slope = f.get()->GetParameter(1);
        gP.at(dG.at(i).e).offset = f.get()->GetParameter(0);
    }
    c.get()->Print(psName.c_str());
    c.get()->Print((psName + ']').c_str());

    std::map<std::string, std::vector<Rep>> rep;

    std::regex m1{"_povtor_\\d+"};
    auto data1{getFitResults(fileName, columnElement1, chem1, m1)};

    for (const auto &item : columnElement1)
    {
        calcRepByElement(data1, item.second, rep, gP.at(item.second));
    }
    for (const auto &item : rep)
    {
        std::cout << item.first << std::endl;
        for (const auto &entry : item.second)
        {
            std::cout << entry.e << " " << entry.avg << " " << entry.stdAbs << " " << entry.stdRel << std::endl;
        }

    }

    return 0;
}

std::map<std::string, Data> getFitResults(const std::string &fileName,
                   const std::map<int, std::string> &columnElement,
                   const std::map<std::string, std::vector<ChemResult>> &chem,
                   const std::regex &pattern)
{
    std::ifstream ifs(fileName);
    if (!ifs.is_open())
    {
        throw my_error("Can't open file \"" + fileName + "\"");
    }
    std::string line;

    std::map<std::string, Data> data;

    while (getline(ifs, line))
    {
        auto strs{splitLineToStrs(line)};
        try
        {
            auto it = std::find_if(chem.begin(), chem.end(), [&strs, &pattern] (std::pair<std::string, std::vector<ChemResult>> chemItem){
                return strs.at(0).find(chemItem.first) != std::string::npos && std::regex_search(strs.at(0), pattern);
            });

            if (it != chem.end())
            {
                std::cout << strs.at(0) << std::endl;
                std::vector<FitResult> fR;
                for (const auto &item : columnElement)
                {
                    fR.push_back({ item.second, strToDouble(strs.at(static_cast<unsigned int>(item.first))),
                                   strToDouble(strs.at(static_cast<unsigned int>(item.first + 1))) });
                }
                data[(*it).first].cr = it->second;
                data[(*it).first].fr.push_back(fR);
            }
        }  catch (...) {
            std::cout << "Error adding e to data";
        }
    }
    ifs.close();
    return data;
}

std::vector<std::string>splitLineToStrs(const std::string &line)
{
    std::stringstream ss(line);
    std::string str;
    std::vector<std::string> strs;
    while (ss >> str)
    {
        strs.push_back(str);
    }
    return strs;

}

double strToDouble(std::string str)
{
    double d;
    std::stringstream ss(str);
    ss >> d;
    if (ss.fail())
    {
        throw my_error("Can\'t convert: " + str);
    }
    return d;
}

void calcRepByElement(const std::map<std::string, Data> &data, const std::string &element, std::map<std::string, std::vector<Rep>> &m, GradPar gP)
{
    for (auto it{data.begin()}; it != data.end(); ++it)
    {
            std::vector<double> r;
            for (size_t i{0}; i < it->second.fr.size(); ++i)
            {
                auto it1 = std::find_if(it->second.fr.at(i).begin(), it->second.fr.at(i).end(), [&element] (FitResult item){
                    return item.e == element;
                });

                if (it1 != it->second.fr.at(i).end())
                {
                    auto v{(*it1).value};
                    if (v.has_value())
                    {
                        r.push_back(v.value() * gP.slope + gP.offset);
                    }
                }
            }
            auto avg{std::accumulate(r.begin(), r.end(), 0.0) / r.size()};
            auto stdAbs{TMath::RMS(r.begin(), r.end())};
            m[(*it).first].push_back({element, avg, stdAbs, 100.0 * stdAbs / avg});
    }
}

void updateDataGradByElement(const std::map<std::string, Data> &data, const std::string &element, DataGrad &dG)
{
    std::vector<double> x, xErr;
    std::vector<double> y, yErr;
    for (auto it{data.begin()}; it != data.end(); ++it)
    {
        auto it1 = std::find_if(it->second.cr.begin(), it->second.cr.end(), [&element] (ChemResult item){
            return item.e == element;
        });
        if (it1 != it->second.cr.end())
        {

            std::vector<double> r, rErr;
            for (size_t i{0}; i < it->second.fr.size(); ++i)
            {
                auto it2 = std::find_if(it->second.fr.at(i).begin(), it->second.fr.at(i).end(), [&element] (FitResult item){
                    return item.e == element;
                });

                if (it2 != it->second.fr.at(i).end())
                {
                    auto v{(*it2).value};
                    auto vErr{(*it2).valueError};
                    if (v.has_value() && vErr.has_value())
                    {
                        r.push_back(v.value());
                        rErr.push_back(vErr.value());
                    }
                }
            }
            if (!r.empty() && !rErr.empty())
            {
                auto rAvg{std::accumulate(r.begin(), r.end(), 0.0) / r.size()};
                auto rErrAvg{std::accumulate(rErr.begin(), rErr.end(), 0.0) / rErr.size()};
                auto v{(*it1).value};
                auto vErr{(*it1).valueError};
                if (v.has_value() && vErr.has_value())
                {
                    x.push_back(rAvg);
                    xErr.push_back(rErrAvg);
                    y.push_back(v.value());
                    yErr.push_back(vErr.value() == 0.0 ? 0.1 : vErr.value());
                }

            }
        }
        else
        {
            return;
        }
    }
    if (!x.empty() && !y.empty())
    {
        dG.e = element;
        dG.x = x;
        dG.xErr = xErr;
        dG.y = y;
        dG.yErr = yErr;
    }
}
