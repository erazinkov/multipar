#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <optional>
#include <exception>
#include <cmath>
#include <numeric>

#include <TVirtualFitter.h>
#include <TLatex.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TF1.h>

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
    FitFunction_2(const std::map<int, std::vector<double>> &d, const int &aNumber)
        : _d{d}, _aNumber{aNumber} {}

    double operator() (double *x, double *par)
    {
        double arg{x[0]};
        int idx{ std::min(static_cast<int>(std::round(arg)), static_cast<int>(_d.size() - 1)) };
        double val{0.0};
        auto A{
           (par[0] - par[1] * _d.at(idx).at(0) - par[5] * (par[2] * _d.at(idx).at(1) + par[4]))
                   / (1.0 - par[5] * par[3])
        };

        if (idx < _aNumber)
        {
            val = A;
        }
        else
        {
            val = par[2] * _d.at(idx).at(1) - par[3] * A + par[4];
        }
        return val;
   }
private:
    const std::map<int, std::vector<double>> _d;
    const int _aNumber;
};

struct Data {
    enum class Value {
        A,
        W
    };
    std::optional<double> a;
    std::optional<double> w;
};


struct Points {
    std::vector<std::string> l;
    std::vector<double> x;
    std::vector<double> xErr;
    std::vector<double> y;
    std::vector<double> yErr;
};

std::map<std::string, Data1> getFitResults(const std::string &fileName,
                                           const std::map<int, std::string> &columnElement,
                                           const std::map<std::string, ChemResult> &chem,
                                           const std::regex &pattern);

std::vector<std::string> splitLineToStrs(const std::string &line);
std::map<std::string, std::vector<std::string>> getStrMapByChem(const std::string &fileName,
                                                                 const size_t column,
                                                                 std::map<std::string, Data> chem);
void addPointsByValue(const std::map<std::string, Data1> &data,
                        Points &points,
                        const Data1::Value value);
//void addPointsByValue(const std::map<std::string, Data> &chem,
//                      Points &points,
//                      const Data::Value value);
double strToDouble(std::string str);

void calcConv(const std::shared_ptr<TF1> &f,
              const std::map<std::string, Data> &chem,
              const std::map<std::string, std::vector<std::string>> &vC,
              const std::map<std::string, std::vector<std::string>> &vO)
{
    auto p0{f.get()->GetParameter(0)};
    auto p1{f.get()->GetParameter(1)};
    auto p2{f.get()->GetParameter(2)};
    auto p3{f.get()->GetParameter(3)};
    auto p4{f.get()->GetParameter(4)};
    auto p5{f.get()->GetParameter(5)};

    auto sumA2{0.0};
    auto sumA{0.0};
    auto sumW2{0.0};
    auto sumW{0.0};
    for (const auto &item : chem)
    {
        if (item.second.a.has_value())
        {
            auto C{strToDouble(vC.at(item.first).front())};
            auto O{strToDouble(vO.at(item.first).front())};
            auto a{
                (p0 - p1 * C - p5 * (p2 * O + p4)) / (1 - p5 * p3)
            };
            auto w{
                p2 * O - p3 * a + p4
            };
            sumA2 += (item.second.a.value() - a) * (item.second.a.value() - a);
            sumA += a;
            if (item.second.w.has_value())
            {
                sumW2 += (item.second.w.value() - w) * (item.second.w.value() - w);
                sumW += w;
            }
        }
    }
    auto stdAabs{std::sqrt(sumA2 / chem.size())};
    auto stdWabs{std::sqrt(sumW2 / chem.size())};
    std::cout << "convergence" << std::endl;
    std::cout << "avgA = " << sumA / chem.size() << std::endl;
    std::cout << "stdAabs = " << stdAabs << std::endl;
    std::cout << "avgW = " << sumW / chem.size() << std::endl;
    std::cout << "stdWabs = " << stdWabs << std::endl;
}


void calcRep(const std::shared_ptr<TF1> &f,
             const std::map<std::string, Data> &chem,
             const std::map<std::string, std::vector<std::string>> &vC,
             const std::map<std::string, std::vector<std::string>> &vO)
{
    auto p0{f.get()->GetParameter(0)};
    auto p1{f.get()->GetParameter(1)};
    auto p2{f.get()->GetParameter(2)};
    auto p3{f.get()->GetParameter(3)};
    auto p4{f.get()->GetParameter(4)};
    auto p5{f.get()->GetParameter(5)};

    std::vector<double> A;
    std::vector<double> W;
    for (const auto &item : chem)
    {
        if (item.second.a.has_value())
        {
            // std::cout << item.first << " " << vC.at(item.first).front() << std::endl;
            auto C{strToDouble(vC.at(item.first).front())};
            auto O{strToDouble(vO.at(item.first).front())};
            auto a{
                (p0 - p1 * C - p5 * (p2 * O + p4)) / (1 - p5 * p3)
            };
            auto w{
                p2 * O - p3 * a + p4
            };
            A.push_back(a);
            if (item.second.w.has_value())
            {
                W.push_back(w);
            }
        }
    }

    auto avgA{std::accumulate(A.begin(), A.end(), 0.0)};
    avgA /= A.size();
    auto avgW{std::accumulate(W.begin(), W.end(), 0.0)};
    avgW /= W.size();

    auto stdAabs{TMath::RMS(A.begin(), A.end())};
    auto stdWabs{TMath::RMS(W.begin(), W.end())};
    std::cout << "repeatability" << std::endl;
    std::cout << "avgA = " << avgA << std::endl;
    std::cout << "stdAabs = " << stdAabs << std::endl;
    std::cout << "avgW = " << avgW << std::endl;
    std::cout << "stdWabs = " << stdWabs << std::endl;
}

void addMmnByValue(const std::map<std::string, Data1> &data,
                        std::map<int, std::vector<double>> &mmn,
                        const Data1::Value value)
{


    int xx{static_cast<int>(mmn.size())};
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
                mmn[xx].push_back(it->second.fr.at(i).at(static_cast<size_t>(0)).value);
                mmn[xx].push_back(it->second.fr.at(i).at(static_cast<size_t>(1)).value);
                xx++;
            }
        }

    }
}

void calcConv(const std::map<std::string, Data1> &data,
              const std::unique_ptr<TF1> &f,
              const Data1::Value value)
{
    auto p0{f.get()->GetParameter(0)};
    auto p1{f.get()->GetParameter(1)};
    auto p2{f.get()->GetParameter(2)};
    auto p3{f.get()->GetParameter(3)};
    auto p4{f.get()->GetParameter(4)};
    auto p5{f.get()->GetParameter(5)};

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
                auto C{it->second.fr.at(i).at(static_cast<size_t>(0)).value};
                auto O{it->second.fr.at(i).at(static_cast<size_t>(1)).value};

                auto res{0.0};

                auto a{
                    (p0 - p1 * C - p5 * (p2 * O + p4)) / (1 - p5 * p3)
                };
                auto w{
                    p2 * O - p3 * a + p4
                };

                switch (value) {
                case Data1::Value::A:
                    res = a;
                    break;
                case Data1::Value::W:
                    res = w;
                    break;
                }

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
            { std::make_pair("coal_blind", kRed), Points() },
            { std::make_pair("barz_blind", kBlue), Points() },
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
            if (!isOther)
            {
                subPoints.at({"other", kMagenta}).l.push_back(points.l.at(i));
                subPoints.at({"other", kMagenta}).x.push_back(points.x.at(i));
                subPoints.at({"other", kMagenta}).y.push_back(points.y.at(i));
                subPoints.at({"other", kMagenta}).xErr.push_back(0.1);
                subPoints.at({"other", kMagenta}).yErr.push_back(0.5);
            }
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
//        item.DrawClone("SAME");
    }
    c.get()->Print(psName.c_str());
    c.get()->Print((psName + ']').c_str());
    c.get()->Close();
}


int main()
{
    // TVirtualFitter::SetDefaultFitter("Minuit");
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
    std::map<std::string, ChemResult> chemBlind
    {
        { "tochka_5_t", {3.18, 33.48} },
        { "tochka_2_t", {1.43, 19.25} },
        { "tochka_3_t", {1.75, 23.89} },
        { "tochka_4_t", {2.23, 27.53} },
        { "tochka_1_t", {1.27, 19.80} },
        { "tochka_6_t", {1.35, 19.47} },
        { "tochka_7_t", {1.46, 19.03} },
        { "tochka_8_t", {1.46, 19.31} },
        { "tochka_9_t", {1.01, 16.30} },

    };


    const std::map<int, std::string> columnElement
    {
        {1, "Al"},
        {3, "C"},
        {5, "Ca"},
        {7, "Fe"},
        {9, "O"},
        {11, "Si"},
    };

//    const std::map<std::string, size_t> elementColumn
//    {
//        // all
//        // {"C", 3},
//        // {"O", 13},
//        //wo Mg
//        // {"C", 3},
//        // {"O", 11},
//        //wo MgCa
//        // {"C", 3},
//        // {"O", 9},
//        //wo MgCaFe
//        // {"C", 3},
//        // {"O", 7},
//        //wo MgCaFeS
//        // {"C", 3},
//        // {"O", 7},
//        //wo MgCaFeSN
//        // {"C", 3},
//        // {"O", 5},
//        //wo MgCaFeSNAl
//        {"C", 1},
//        {"O", 3},
//    };

    // const auto fileName{"rea.elts.txt_shahta12_all"};
    // const auto fileName{"rea.elts.txt_shahta12_wo_Mg"};
    // const auto fileName{"rea.elts.txt_shahta12_wo_MgCa"};
    // const auto fileName{"rea.elts.txt_shahta12_wo_MgCaFe"};
    // const auto fileName{"rea.elts.txt_shahta12_wo_MgCaFeS"};
    // const auto fileName{"rea.elts.txt_shahta12_wo_MgCaFeSN"};
    // const auto fileName{"rea.elts.txt_shahta12_wo_MgCaFeSNAl"};
//    const auto fileName{"rea.elts.txt_bereza_wo_MgCaFeSNAl"};
//    const auto fileName{"rea.elts.txt.12_w_bereza_w_barz_wo_MgCaFeSNAl.grad_w_blind.all"};
    const auto fileName{"rea.elts.txt.stavropol_t"};
    std::cout << fileName << std::endl;

//    chem.insert(chemBlind.begin(), chemBlind.end());

    try
    {
//        std::regex m{"\\d+_\\d\\."};
        std::regex m{"\\d+_(s|t)"};
        auto data1{getFitResults(fileName, columnElement, chem, m)};

        Points points;

        auto value{Data1::Value::A};

        addPointsByValue(data1, points, Data1::Value::A);
        auto aNumber{points.x.size()};
        addPointsByValue(data1, points, Data1::Value::W);
//        auto wNumber{points.x.size() - aNumber};
        std::map<int, std::vector<double>> mmn;
//        size_t idx{ static_cast<size_t>(std::round(points.x.front())) };
        addMmnByValue(data1, mmn, Data1::Value::A);
        addMmnByValue(data1, mmn, Data1::Value::W);
        std::cout << points.x.size() << " " << mmn.size() << std::endl;

        std::unique_ptr<TGraphErrors> gr{new TGraphErrors(static_cast<int>(points.x.size()), &points.x[0], &points.y[0], &points.xErr[0], &points.yErr[0])};
        gr.get()->SetMarkerSize(1.5);
        gr.get()->SetMarkerStyle(21);
        gr.get()->SetTitle(";N_{probe};[...A, ...W]");

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

        FitFunction_2 fObj(mmn, static_cast<int>(aNumber));
        std::unique_ptr<TF1> f{new TF1("f", fObj, points.x.front(), points.x.back(), 6)};
//        f.get()->SetParameter(0, 100.0);
//        f.get()->SetParameter(1, 1.0);
//        f.get()->SetParameter(2, 1.0);
//        f.get()->SetParameter(3, 0.5);
//        f.get()->SetParameter(4, -1.0);
//        f.get()->SetParameter(5, 1.0);
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
                { std::make_pair("coal_blind", kRed), Points() },
                { std::make_pair("barz_blind", kBlue), Points() },
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

//        std::regex s{"sum"};
        std::regex s{"\\d+_s"};
        auto data1Sum{getFitResults(fileName, columnElement, chem, s)};
        calcConv(data1Sum, f, value);

//        std::regex t{"\\d+_t"};
//        auto dataBlindSum{getFitResults(fileName, columnElement, chemBlind, t)};
//        dataBlindSum.insert(data1Sum.begin(), data1Sum.end());
//        calcConv(dataBlindSum, f, value);

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
    int xx{static_cast<int>(points.x.size())};
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

//void addPointsByValue(const std::map<std::string, Data> &chem,
//                        Points &points,
//                        const Data::Value value)
//{

//    int xx{static_cast<int>(points.x.size())};
//    for (auto it{chem.begin()}; it != chem.end(); ++it)
//    {
//        std::optional<double> v;
//        switch (value) {
//            case Data::Value::A:
//                v = (*it).second.a;
//            break;
//            case Data::Value::W:
//                v = (*it).second.w;
//            break;
//        }
//        if (v.has_value())
//        {
//            points.l.push_back((*it).first);
//            points.x.push_back(xx++);
//            points.xErr.push_back(0.01);
//            points.y.push_back(v.value());
//            points.yErr.push_back(0.5);
//        }
//    }

//}


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

std::map<std::string, std::vector<std::string>> getStrMapByChem(const std::string &fileName,
                                                                const size_t column,
                                                                std::map<std::string, Data> chem)
{
    std::ifstream ifs(fileName);
    if (!ifs.is_open())
    {
        throw my_error("Can't open file \"" + fileName + "\"");
    }
    std::map<std::string, std::vector<std::string>> map;
    std::string line;

    while (getline(ifs, line))
    {
        auto strs{splitLineToStrs(line)};
        try
        {
            auto it = std::find_if(chem.begin(), chem.end(), [&strs] (std::pair<std::string, Data> chemItem){
                    return strs.at(0).find(chemItem.first) != std::string::npos;
            });

            if (it != chem.end())
            {
                map[(*it).first].push_back(strs.at(column));
            }
        }  catch (...) {
            std::cout << "Incorrect column " << column;
        }
    }
    return map;
}
