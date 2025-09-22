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
    double value;
    double valueError;
    void print() const {
        std::cout << e << " " << value << "\u00B1" << valueError << " ";
    }
};

struct ChemResult {
    enum class Value {
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

struct Data1 {
    enum class Value {
        A,
        W
    };
    ChemResult chem;
    std::vector<std::vector<FitResult>> fr;
    void print() const {
        auto a{this->chem.a.has_value() ?  this->chem.a.value() : -1};
        auto w{this->chem.w.has_value() ?  this->chem.w.value() : -1};
        std::cout << a << " " << w << " ";
        for (const auto& eItem : this->fr)
        {
            for (const auto& eItemItem : eItem)
            {
                eItemItem.print();
            }
        }
    }
};

class FitFunction_2
{
public:
    FitFunction_2(std::map<int, std::vector<FitResult>> fR) : _fR{fR} {}

    double operator() (double *x, double *par)
    {
        double arg{x[0]};
        int idx{ std::min(static_cast<int>(std::round(arg)), static_cast<int>(_fR.size() - 1)) };
        const auto nPar{_fR.begin()->second.size()};
        double val{0.0};
        for (size_t i{0}; i < nPar; ++i)
        {
            val += par[i] * _fR.at(idx).at(i).value;
        }
        val += par[nPar];
        return val;
    }
private:
    std::map<int, std::vector<FitResult>> _fR;
};

struct Points {
    std::vector<std::string> l;
    std::vector<double> x;
    std::vector<double> xErr;
    std::vector<double> y;
    std::vector<double> yErr;
};

std::vector<std::string> splitLineToStrs(const std::string &line);

double strToDouble(std::string str);



std::map<std::string, Data1> getFitResults(const std::string &fileName,
                                           const std::map<int, std::string> &columnElement,
                                           const std::map<std::string, ChemResult> &chem,
                                           const std::regex &pattern);

std::map<int, std::vector<FitResult>> getFitResultsByValue(const std::map<std::string, Data1> &data,
                                                           const Data1::Value value);

void addPointsByValue(const std::map<std::string, Data1> &data,
                      Points &points,
                      const Data1::Value value);

void calcRep(const std::map<std::string, Data1> &data,
             const std::unique_ptr<TF1> &f);

void calcConv(const std::map<std::string, Data1> &data,
              const std::unique_ptr<TF1> &f,
              const Data1::Value value);

int main()
{
    std::map<std::string, ChemResult> chemBlind
    {
        { "tochka_5_t", {3.18, std::nullopt} },
        { "tochka_2_t", {1.43, std::nullopt} },
        { "tochka_3_t", {1.75, std::nullopt} },
        { "tochka_4_t", {2.23, std::nullopt} },
        { "tochka_1_t", {1.27, std::nullopt} },
        { "tochka_6_t", {1.35, std::nullopt} },
        { "tochka_7_t", {1.46, std::nullopt} },
        { "tochka_8_t", {1.46, std::nullopt} },
        { "tochka_9_t", {1.01, std::nullopt} },

    };
    //
    const std::map<int, std::string> columnElement
    {
        {1, "Al"},
        {3, "C"},
        {5, "Ca"},
        {7, "Fe"},
        {9, "O"},
        {11, "Si"},
    };
//    //wo_MgCaFeS
//    const std::map<int, std::string> columnElement
//        {
//         {1, "Al"},
//         {3, "C"},
//         {5, "N"},
//         {7, "O"},
//         {9, "Si"},
//         };
//    //wo_MgS
//    const std::map<int, std::string> columnElement
//        {
//         {1,  "Al"},
//         {3,  "C"},
//         {5,  "Ca"},
//         {7,  "Fe"},
//         {9,  "N"},
//         {11, "O"},
//         {13, "Si"},
//         };
//    //all
//    const std::map<int, std::string> columnElement
//    {
//        {1,  "Al"},
//        {3,  "C"},
//        {5,  "Ca"},
//        {7,  "Fe"},
//        {9,  "Mg"},
//        {11, "N"},
//        {13, "O"},
//        {15, "S"},
//        {17, "Si"},
//    };
//    std::string fileName{"rea.elts.txt.12.wo_MgS" }; // wo_MgCaFeS 12
//    std::string fileName{"rea.elts.txt.12.wo_MgCaFeS.all" }; // wo_MgCaFeS 12
//     std::string fileName{"rea.elts.txt.bereza_wo_MgCaFeS.all" }; // wo_MgCaFeS bereza

//    std::string fileName{"rea.elts.txt.12_w_bereza_wo_MgCaFeS.all"}; // wo_MgCaFeS 12+bereza
//    std::string fileName{"rea.elts.txt.12_w_bereza_w_barz_wo_MgCaFeS.grad_w_blind.all"};
//    std::string fileNameBlind{"rea.elts.txt.12_w_bereza_w_barz_wo_MgCaFeS.blind"}; // wo_MgCaFeS barz+12+bereza
//    std::string fileName{"rea.elts.txt.12_w_bereza_w_barz_all.grad_w_blind.all"};

    std::string fileName{"rea.elts.txt.stavropol_t"};

    std::map<std::string, ChemResult> chem
    {
        { "tochka_1_s", {4.74, std::nullopt} },
        { "tochka_2_s", {4.66, std::nullopt} },
        { "tochka_3_s", {4.99, std::nullopt} },
        { "tochka_4_s", {5.14, std::nullopt} },
        { "tochka_5_s", {4.74, std::nullopt} },
        { "tochka_6_s", {4.45, std::nullopt} },
        { "tochka_7_s", {4.78, std::nullopt} },
        { "tochka_8_s", {4.80, std::nullopt} },
        { "tochka_9_s", {3.83, std::nullopt} },
    };

    chem.insert(chemBlind.begin(), chemBlind.end());

    try
    {

        // std::regex p{"_povtor_\\d+\\."};
//        std::regex m{"\\d+_\\d+\\."};
        std::regex m{"\\d+_(s|t)"};
        auto data1{getFitResults(fileName, columnElement, chem, m)};

        Points points;

        auto value{Data1::Value::A};

        addPointsByValue(data1, points, value);

        auto fitResultsByValue{getFitResultsByValue(data1, value)};

        std::cout << fitResultsByValue.size() << std::endl;

        std::unique_ptr<TGraphErrors> gr{new TGraphErrors(static_cast<int>(points.x.size()), &points.x[0], &points.y[0], &points.xErr[0], &points.yErr[0])};
        gr.get()->SetMarkerSize(1.5);
        gr.get()->SetMarkerStyle(21);
        std::string grTitle{";N_probe;"};
        grTitle.append(value == Data1::Value::A ? "A" : "W");
        gr.get()->SetTitle(grTitle.c_str());

        std::vector<TLatex> labels;

        for (size_t i{0}; i < points.x.size(); ++i)
        {
            auto pos{points.l.at(i).find_last_of("_")};
            auto text{points.l.at(i)};
            if (pos != std::string::npos && pos == points.l.at(i).length() - 1)
            {
                text = text.substr(0, pos);
            }
            TLatex l(points.x.at(i), points.y.at(i) + 1.25 * points.yErr.at(i), text.c_str());
            l.SetTextAngle(90);
            l.SetTextAlign(12);
            l.SetTextSize(0.02);
            labels.push_back(l);
        }

        FitFunction_2 fObj(fitResultsByValue);
        std::unique_ptr<TF1> f{new TF1("f", fObj, points.x.front(), points.x.back(), static_cast<int>(columnElement.size() + 1))};

//        f.get()->SetParLimits(1, -5.0, 0.0);
//        f.get()->SetParLimits(f->GetNpar() - 1, 50.0, 150.0);
        f.get()->SetNpx(10 * static_cast<int>(points.x.size()));

        gr.get()->Fit(f.get(), "R");

        const std::string psName{"output.ps"};
        std::unique_ptr<TCanvas> c{new TCanvas("c", "c", 1024, 960)};
        c.get()->Print((psName + '[').c_str());
        gr.get()->Draw("APL");

        auto useSub{true};
        if (useSub)
        {
            std::map<std::pair<std::string, Color_t>, Points> subPoints{
                { std::make_pair("_s", kRed), Points() },
                { std::make_pair("_t", kBlue), Points() },
                { std::make_pair("bereza_blind", kGreen), Points() },
                { std::make_pair("other", kBlack), Points() },
            };


            for (size_t i{0}; i < points.x.size(); ++i)
            {
                auto isOther{false};
                for (auto &item : subPoints)
                {
                    if (points.l.at(i).find(item.first.first) != std::string::npos)
                    {
                        TMarker m{points.x.at(i), points.y.at(i), 21};
                        m.SetMarkerSize(1.5);
                        m.SetMarkerColor(item.first.second);
                        m.DrawClone("SAME");
                        item.second.l.push_back(points.l.at(i));
                        item.second.x.push_back(points.x.at(i));
                        item.second.y.push_back(points.y.at(i));
                        item.second.xErr.push_back(0.1);
                        item.second.yErr.push_back(0.5);

                        isOther = true;
                    }
                }
                if (!isOther)
                {
                    subPoints.at({"other", kBlack}).l.push_back(points.l.at(i));
                    subPoints.at({"other", kBlack}).x.push_back(points.x.at(i));
                    subPoints.at({"other", kBlack}).y.push_back(points.y.at(i));
                    subPoints.at({"other", kBlack}).xErr.push_back(0.1);
                    subPoints.at({"other", kBlack}).yErr.push_back(0.5);
                }
            }
        }

        for (const auto &item : labels)
        {
//            item.DrawClone("SAME");
        }
        c.get()->Print(psName.c_str());
        c.get()->Print((psName + ']').c_str());
        c.get()->Close();

        std::regex s{"\\d+_s"};
//         std::regex s{"\\d+_\\d+\\."};
        auto data1Sum{getFitResults(fileName, columnElement, chem, s)};
//        calcConv(data1Sum, f, value);
//        std::regex p{"_povtor_\\d+\\."};
//        auto data1P{getFitResults(fileName, columnElement, chem, p)};
//        calcRep(data1P, f);

//        std::string fileNameBlind{"rea.elts.txt.12_w_bereza_w_barz_wo_MgCaFeS.blind"}; // wo_MgCaFeS barz+12+bereza
        std::regex t{"\\d+_t"};
        auto dataBlindSum{getFitResults(fileName, columnElement, chemBlind, t)};
        dataBlindSum.insert(data1Sum.begin(), data1Sum.end());
        calcConv(dataBlindSum, f, value);
    }
    catch (const my_error& err)
    {
        std::cout << "Error: " << err.what() << std::endl;
    }
    catch (const std::exception& err)
    {
        std::cout << "Error: " << err.what() << std::endl;
    }
    return 0;
}

void calcConv(const std::map<std::string, Data1> &data,
              const std::unique_ptr<TF1> &f,
              const Data1::Value value)
{
    Points points;
    for (auto it{data.begin()}; it != data.end(); ++it)
    {
        for (size_t i{0}; i < it->second.fr.size(); ++i)
        {
            std::optional<double> v;
            switch (value) {
            case Data1::Value::A:
                v = (*it).second.chem.a;
                break;
            case Data1::Value::W:
                v = (*it).second.chem.w;
                break;
            }
            if (v.has_value())
            {
                auto nPar{f.get()->GetNpar()};
                auto res{0.0};
                for (auto pIdx{0}; pIdx < nPar - 1; ++pIdx)
                {
                    res += f->GetParameter(pIdx) * it->second.fr.at(i).at(static_cast<unsigned int>(pIdx)).value;
                }
                res += f->GetParameter(nPar - 1);
                points.l.push_back(it->first);
                points.x.push_back(res);
                points.y.push_back(v.value());
                points.xErr.push_back(0.1);
                points.yErr.push_back(0.5);
            }
        }
    }
    std::vector<double> d2;
    for (size_t i{0}; i < points.x.size(); ++i)
    {
        d2.push_back(std::pow(points.y.at(i) - points.x.at(i), 2));
    }
    auto stdAbs{std::sqrt(std::accumulate(d2.begin(), d2.end(), 0.0) / d2.size())};
    auto avg{std::accumulate(points.x.begin(), points.x.end(), 0.0) / points.x.size()};
    std::cout << "convergence: " << "avg = " << avg << " stdAbs = " << stdAbs << std::endl;

    std::unique_ptr<TGraphErrors> gr{new TGraphErrors(static_cast<int>(points.x.size()), &points.x[0], &points.y[0], &points.xErr[0], &points.yErr[0])};
    gr.get()->SetMarkerSize(1.5);
    gr.get()->SetMarkerStyle(21);

    std::vector<TLatex> labels;

    for (size_t i{0}; i < points.x.size(); ++i)
    {
        auto pos{points.l.at(i).find_last_of("_")};
        auto text{points.l.at(i)};
        if (pos != std::string::npos && pos == points.l.at(i).length() - 1)
        {
            text = text.substr(0, pos);
        }
        TLatex l(points.x.at(i), points.y.at(i) + 1.25 * points.yErr.at(i), text.c_str());
        l.SetTextAngle(90);
        l.SetTextAlign(12);
        l.SetTextSize(0.02);
        labels.push_back(l);
    }

    std::stringstream ss;
    ss.str("");ss.clear();
    ss << (value == Data1::Value::A ? "Ad" : "Wr")
       << ": avg = "
       << avg << "%"
       << ", stdAbs = "
       << stdAbs << "%"
       << ";AGP-K, %;Chem, %";

    std::unique_ptr<TH2D> h2dConv{new TH2D("h2dConv",
                                           ss.str().c_str(),
                                           static_cast<int>(points.y.size()),
                                           0.75 * (*std::min_element(points.y.begin(), points.y.end())),
                                           1.25 * (*std::max_element(points.y.begin(), points.y.end())),
                                           static_cast<int>(points.y.size()),
                                           0.75 * (*std::min_element(points.y.begin(), points.y.end())),
                                           1.25 * (*std::max_element(points.y.begin(), points.y.end())))};
    h2dConv->SetStats(0);
    std::unique_ptr<TLine> lConv{new TLine(0.75 * (*std::min_element(points.y.begin(), points.y.end())),
                                           0.75 * (*std::min_element(points.y.begin(), points.y.end())),
                                           1.25 * (*std::max_element(points.y.begin(), points.y.end())),
                                           1.25 * (*std::max_element(points.y.begin(), points.y.end())))};
    const std::string psName{"output_conv.ps"};
    std::unique_ptr<TCanvas> c{new TCanvas("c", "c", 1024, 960)};
    c.get()->Print((psName + '[').c_str());
    h2dConv.get()->Draw();
    gr.get()->Draw("P");

    auto useSub{true};
    if (useSub)
    {
        std::map<std::pair<std::string, Color_t>, Points> subPoints{
            { std::make_pair("_s", kRed), Points() },
            { std::make_pair("_t", kBlue), Points() },
            { std::make_pair("bereza_blind", kGreen), Points() },
            { std::make_pair("other", kMagenta), Points() },
        };


        for (size_t i{0}; i < points.x.size(); ++i)
        {
            auto isOther{false};
            for (auto &item : subPoints)
            {
                if (points.l.at(i).find(item.first.first) != std::string::npos)
                {
                    TMarker m{points.x.at(i), points.y.at(i), 21};
                    m.SetMarkerSize(1.5);
                    m.SetMarkerColor(item.first.second);
                    m.DrawClone("SAME");
                    item.second.l.push_back(points.l.at(i));
                    item.second.x.push_back(points.x.at(i));
                    item.second.y.push_back(points.y.at(i));
                    item.second.xErr.push_back(0.1);
                    item.second.yErr.push_back(0.5);

                    isOther = true;
                }
            }
//            if (!isOther)
//            {
                subPoints.at({"other", kMagenta}).l.push_back(points.l.at(i));
                subPoints.at({"other", kMagenta}).x.push_back(points.x.at(i));
                subPoints.at({"other", kMagenta}).y.push_back(points.y.at(i));
                subPoints.at({"other", kMagenta}).xErr.push_back(0.1);
                subPoints.at({"other", kMagenta}).yErr.push_back(0.5);
//            }
        }
        std::stringstream ss;
        ss.str("");ss.clear();
        ss << (value == Data1::Value::A ? "Ad" : "Wr");
        ss << ": stdAbs=" << std::setprecision(3);
        auto stdAbs1 = [](Points points){
            std::vector<double> d2;
            for (size_t i{0}; i < points.x.size(); ++i)
            {
                d2.push_back(std::pow(points.y.at(i) - points.x.at(i), 2));
            }
            return std::sqrt(std::accumulate(d2.begin(), d2.end(), 0.0) / static_cast<double>(d2.size()));
        };
//        auto avg1 = [](Points points){
//            return std::accumulate(points.x.begin(), points.x.end(), 0.0) / points.x.size();
//        };
        for (auto item : subPoints)
        {
            ss << "[#color[" << static_cast<int>(item.first.second) << "]{" << stdAbs1(item.second) << "}] ";
        }

        ss << ";AGP-K, %;Chem, %";
        h2dConv.get()->SetTitle(ss.str().c_str());
    }

    lConv.get()->Draw("SAME");
    for (const auto &item : labels)
    {
        item.DrawClone("SAME");
    }
    c.get()->Print(psName.c_str());
    c.get()->Print((psName + ']').c_str());
    c.get()->Close();
}

void calcRep(const std::map<std::string, Data1> &data,
             const std::unique_ptr<TF1> &f)
{

    std::vector<double> r;
    for (auto it{data.begin()}; it != data.end(); ++it)
    {
        std::cout << (*it).first << " " << it->second.fr.size() << std::endl;
        for (size_t i{0}; i < it->second.fr.size(); ++i)
        {
            auto nPar{f.get()->GetNpar()};
            auto res{0.0};
            for (auto pIdx{0}; pIdx < nPar - 1; ++pIdx)
            {
                res += f->GetParameter(pIdx) * it->second.fr.at(i).at(static_cast<unsigned int>(pIdx)).value;
            }
            res += f->GetParameter(nPar - 1);
            r.push_back(res);
        }
    }
    auto avg{std::accumulate(r.begin(), r.end(), 0.0) / r.size()};
    auto stdAbs{TMath::RMS(r.begin(), r.end())};
    std::cout << "repeatability: " << "avg = " << avg << " stdAbs = " << stdAbs << std::endl;
}

std::map<std::string, Data1> getFitResults(const std::string &fileName,
                   const std::map<int, std::string> &columnElement,
                   const std::map<std::string, ChemResult> &chem,
                   const std::regex &pattern)
{
    std::ifstream ifs(fileName);
    if (!ifs.is_open())
    {
        throw my_error("Can't open file \"" + fileName + "\"");
    }
    std::string line;

    std::map<std::string, Data1> data;

    while (getline(ifs, line))
    {
        auto strs{splitLineToStrs(line)};
        try
        {
            auto it = std::find_if(chem.begin(), chem.end(), [&strs, &pattern] (std::pair<std::string, ChemResult> chemItem){
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
                data[(*it).first].chem.a = it->second.a;
                data[(*it).first].chem.w = it->second.w;
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

std::map<int, std::vector<FitResult>> getFitResultsByValue(const std::map<std::string, Data1> &data,
                                                           const Data1::Value value)
{
    std::map<int, std::vector<FitResult>> fitResultsByValue;
    int xx{0};
    for (auto it{data.begin()}; it != data.end(); ++it)
    {
        for (size_t i{0}; i < it->second.fr.size(); ++i)
        {
            std::string label{it->first + "_" + i};
            std::optional<double> v;
            switch (value) {
            case Data1::Value::A:
                v = (*it).second.chem.a;
                break;
            case Data1::Value::W:
                v = (*it).second.chem.w;
                break;
            }
            if (v.has_value())
            {
                fitResultsByValue[xx++].assign(it->second.fr.at(i).begin(), it->second.fr.at(i).end());
            }
        }
    }
    return fitResultsByValue;
}

void addPointsByValue(const std::map<std::string, Data1> &data,
                        Points &points,
                        const Data1::Value value)
{
    int xx{0};
    for (auto it{data.begin()}; it != data.end(); ++it)
    {
        for (size_t i{0}; i < it->second.fr.size(); ++i)
        {
            std::string label{it->first + "_" + i};
            std::optional<double> v;
            switch (value) {
            case Data1::Value::A:
                v = (*it).second.chem.a;
                break;
            case Data1::Value::W:
                v = (*it).second.chem.w;
                break;
            }
            if (v.has_value())
            {
                points.l.push_back(label);
                points.x.push_back(xx++);
                points.xErr.push_back(0.01);
                points.y.push_back(v.value());
                points.yErr.push_back(0.5);
            }
        }

    }
}
