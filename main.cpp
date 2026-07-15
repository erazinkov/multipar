#include <iostream>
#include <sstream>
#include <fstream>
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
#include <TStyle.h>
#include <TPaveText.h>

#include <bits/stdc++.h>

#include <regex>

std::map<std::string, std::pair<double, double>> preGrad{
    { "C", {-1.70467, 1.074623} },
    { "Al", {0.785556, 1.088333} },
    { "N", {-0.33024, 1.281688} },
    { "O", {-6.21944, 1.068012} },
    { "Si", {0.090895, 0.992175} },
    { "UNKNOWN", {0.0, 1.0}},
};

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

//        mmn[xx].push_back(it->second.fr.at(i).at(static_cast<size_t>(0)).value); // Al
//        mmn[xx].push_back(it->second.fr.at(i).at(static_cast<size_t>(1)).value); // C
//        mmn[xx].push_back(it->second.fr.at(i).at(static_cast<size_t>(2)).value); // N
//        mmn[xx].push_back(it->second.fr.at(i).at(static_cast<size_t>(3)).value); // O
//        mmn[xx].push_back(it->second.fr.at(i).at(static_cast<size_t>(4)).value); // Si

        auto W{
            par[3] * _d.at(idx).at(3) // O
            - par[4] * _d.at(idx).at(4) // Si
            - par[5] * _d.at(idx).at(0) // Al
            - par[6] * _d.at(idx).at(2) // N
            + par[7]
        };

        auto A{
           (par[1] + par[0] * _d.at(idx).at(1)) // C
            / (100.0 - par[2] * W) * 100.0
        };

        if (idx < _aNumber)
        {
            val = A;
        }
        else
        {
            val = W;
        }
        return val;
   }
private:
    const std::map<int, std::vector<double>> _d;
    const int _aNumber;
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
void addPointsByValue(const std::map<std::string, Data1> &data,
                        Points &points,
                        const Data1::Value value);
double strToDouble(std::string str);
void addMmnByValue(const std::map<std::string, Data1> &data,
                        std::map<int, std::vector<double>> &mmn,
                        const Data1::Value value);
void calcConv(const std::map<std::string, Data1> &data,
              const std::unique_ptr<TF1> &f,
              const Data1::Value value);

void writePointsToFile(const std::string fileName, const Points &points);


double calculateStdAbs(const Points& points) {
    if (points.x.empty()) return 0.0;

    double sumSquaredDiff = 0.0;
    for (size_t i = 0; i < points.x.size(); ++i) {
        const double diff = points.y[i] - points.x[i];
        std::cout << points.y[i] << " " << points.x[i] << std::endl;
        sumSquaredDiff += diff * diff;
    }
    return std::sqrt(sumSquaredDiff / points.x.size());
}

Points getPredicatedPoints(const Points& points) {
    Points predicatedPoints;
    std::unique_ptr<TGraphErrors> gr{new TGraphErrors(static_cast<int>(points.x.size()), &points.x[0], &points.y[0], &points.xErr[0], &points.yErr[0])};
    auto min = std::min((*std::min_element(points.x.begin(), points.x.end())), (*std::min_element(points.y.begin(), points.y.end())));
    auto max = std::max((*std::max_element(points.x.begin(), points.x.end())), (*std::max_element(points.y.begin(), points.y.end())));
    std::unique_ptr<TF1> f{std::make_unique<TF1>("f", "pol1", 0.75 * min, 1.25 * max)};
    f.get()->SetParLimits(1, 0.0, 100.0);
    gr.get()->Fit(f.get(), "RQ0");
    for (size_t i{0}; i < points.x.size(); ++i) {
        predicatedPoints.l.push_back(points.l.at(i));
        predicatedPoints.x.push_back(f.get()->GetParameter(0) + f.get()->GetParameter(1) * points.x.at(i));
        predicatedPoints.y.push_back(points.y.at(i));
        predicatedPoints.xErr.push_back(points.xErr.at(i));
        predicatedPoints.yErr.push_back(points.yErr.at(i));
    }
    return predicatedPoints;
}

double getOptimalPar(const std::map<std::string, Data1> &data) {

    std::vector<double> pars;
    std::vector<double> stds;
    for (auto par{-0.1}; par > -1.0; par -= 0.01) {
        Points points;
        for (const auto &i : data) {
            points.l.push_back(i.first);
            points.x.push_back(i.second.fr.at(0).at(3).value + par * i.second.chem.a.value());
            points.y.push_back(i.second.chem.w.value());
            points.xErr.push_back(0.0);
            points.yErr.push_back(0.0);
        }
        pars.push_back(par);
        Points predicatedPoints{getPredicatedPoints(points)};
        stds.push_back(calculateStdAbs(predicatedPoints));
    }
    if (!stds.empty() && !pars.empty() && stds.size() == pars.size()) {
        std::unique_ptr<TGraphErrors> gr{new TGraphErrors(static_cast<int>(pars.size()), &pars[0], &stds[0], nullptr, nullptr)};
        gr.get()->SetMarkerSize(1.05);
        gr.get()->SetMarkerStyle(21);
        std::stringstream ss;
        ss << std::setprecision(2);
        auto it{std::min_element(stds.begin(), stds.end())};
        size_t index = std::distance(stds.begin(), it);
        ss << stds.at(index) << " " << pars.at(index) << ";par;stdAbs";
        gr.get()->SetTitle(ss.str().c_str());
        const std::string psName{"output_optimal.ps"};
        std::unique_ptr<TCanvas> c{new TCanvas("c", "c", 1024, 960)};
        c.get()->SetGrid();
        gStyle->SetOptFit(11111);
        c.get()->Print((psName + '[').c_str());

        gr.get()->Draw("APL");

        c.get()->Print(psName.c_str());
        c.get()->Print((psName + ']').c_str());

        return pars.at(index);
    }
    return 0.0;
}

void drawCorrGraphWr(const std::map<std::string, Data1> &data, double par = 0.0) {

    Points points;

    std::cout << "Ad,%" << " " << "Wr,%" << " " << "Oxy,%" << std::endl;
    for (const auto &i : data) {
//        std::cout << i.second.chem.a.value() << " " << i.second.chem.w.value() << " " << i.second.fr.at(0).at(3).value << std::endl;
        std::stringstream ss;
        ss << std::setprecision(3) << i.second.chem.a.value();
        points.l.push_back(i.first);
        points.x.push_back(i.second.fr.at(0).at(3).value + par * i.second.chem.a.value());
        points.y.push_back(i.second.chem.w.value());
        points.xErr.push_back(0.0);
        points.yErr.push_back(0.0);
    }

    const std::string psName{"output_corr.ps"};
    std::unique_ptr<TCanvas> c{new TCanvas("c", "c", 1024, 960)};
    c.get()->SetGrid();
    gStyle->SetOptFit(11111);
    c.get()->Print((psName + '[').c_str());

    auto min = std::min((*std::min_element(points.x.begin(), points.x.end())), (*std::min_element(points.y.begin(), points.y.end())));
    auto max = std::max((*std::max_element(points.x.begin(), points.x.end())), (*std::max_element(points.y.begin(), points.y.end())));



    std::unique_ptr<TH2D> h2dCorr{new TH2D("h2dCorr",
                                           "h2dCorr",
                                           static_cast<int>(points.y.size()),
                                           0.75 * min,
                                           1.25 * max,
                                           static_cast<int>(points.y.size()),
                                           0.75 * min,
                                           1.25 * max)};
    h2dCorr.get()->SetStats(0);
    h2dCorr.get()->Draw();
    std::unique_ptr<TGraphErrors> gr{new TGraphErrors(static_cast<int>(points.x.size()), &points.x[0], &points.y[0], &points.xErr[0], &points.yErr[0])};
    gr.get()->SetMarkerSize(1.05);
    gr.get()->SetMarkerStyle(21);
    std::stringstream ss;
    ss << std::setprecision(2);
    ss << "mOxy'=mOxy+k#bulletAd," << " k=" << par << ";mOxy,%;Wr,%";
    h2dCorr.get()->SetTitle(ss.str().c_str());

    std::unique_ptr<TF1> f{std::make_unique<TF1>("f", "pol1", 0.75 * min, 1.25 * max)};
//    gr.get()->Fit(f.get(), "R");
    gr.get()->Draw("SAME P");

    std::unique_ptr<TLine> l{new TLine(0.75 * min,
                                       0.75 * min,
                                       1.25 * max,
                                       1.25 * max)};
    l.get()->Draw("SAME");


    std::vector<TLatex> labels;
    std::map<std::pair<std::string, Color_t>, Points> subPoints{
        { std::make_pair("N12", kRed), Points() },
        { std::make_pair("berez_blind_", kBlue), Points() },
        { std::make_pair("berez_11_w", kGreen), Points() },
        { std::make_pair("N12_blind", kOrange), Points() },
        { std::make_pair("barz_blind", kMagenta), Points() },
        { std::make_pair("berez_", kYellow), Points() },
        { std::make_pair("kuz_", kCyan), Points() },
        { std::make_pair("other", kBlack), Points() },
    };


    for (size_t i{0}; i < points.x.size(); ++i)
    {
        for (auto &item : subPoints)
        {
            if (points.l.at(i).find(item.first.first) != std::string::npos)
            {
                TMarker m{points.x.at(i), points.y.at(i), 21};
                m.SetMarkerSize(1.05);
                m.SetMarkerColor(item.first.second);
                m.DrawClone("SAME");
                item.second.l.push_back(points.l.at(i));
                item.second.x.push_back(points.x.at(i));
                item.second.y.push_back(points.y.at(i));
//                item.second.xErr.push_back(0.1);
//                item.second.yErr.push_back(0.5);

            }
        }
    }

    std::unique_ptr<TPaveText> pt{new TPaveText(0.1, 0.7, 0.3, 0.9, "NDC")};
    pt.get()->SetFillColor(0);
    pt.get()->SetBorderSize(1);

    // Add entries from map keys
    for (const auto& entry : subPoints) {
        const std::string& key = entry.first.first;
        Color_t color = entry.first.second;

        TText *text = pt.get()->AddText(key.c_str());
        text->SetTextColor(color);
    }
    pt.get()->Draw("SAME");
    c.get()->Print(psName.c_str());
    c.get()->Print((psName + ']').c_str());
}

void drawCorrGraphWr1(const std::map<std::string, Data1> &data, double par = 0.0) {

    Points points1;

    std::cout << "Ad,%" << " " << "Wr,%" << " " << "Oxy,%" << std::endl;
    for (const auto &i : data) {
//        std::cout << i.second.chem.a.value() << " " << i.second.chem.w.value() << " " << i.second.fr.at(0).at(3).value << std::endl;
        std::stringstream ss;
        ss << std::setprecision(3) << i.second.chem.a.value();
        points1.l.push_back(i.first);
        points1.x.push_back(i.second.fr.at(0).at(3).value + par * i.second.chem.a.value());
        points1.y.push_back(i.second.chem.w.value());
        points1.xErr.push_back(0.0);
        points1.yErr.push_back(0.0);
    }

    Points points{getPredicatedPoints(points1)};

    const std::string psName{"output_corr_grad.ps"};
    std::unique_ptr<TCanvas> c{new TCanvas("c", "c", 1024, 960)};
    c.get()->SetGrid();
    gStyle->SetOptFit(11111);
    c.get()->Print((psName + '[').c_str());

    auto min = std::min((*std::min_element(points.x.begin(), points.x.end())), (*std::min_element(points.y.begin(), points.y.end())));
    auto max = std::max((*std::max_element(points.x.begin(), points.x.end())), (*std::max_element(points.y.begin(), points.y.end())));



    std::unique_ptr<TH2D> h2dCorr{new TH2D("h2dCorr",
                                           "h2dCorr",
                                           static_cast<int>(points.y.size()),
                                           0.75 * min,
                                           1.25 * max,
                                           static_cast<int>(points.y.size()),
                                           0.75 * min,
                                           1.25 * max)};
    h2dCorr.get()->SetStats(0);
    h2dCorr.get()->Draw();
    std::unique_ptr<TGraphErrors> gr{new TGraphErrors(static_cast<int>(points.x.size()), &points.x[0], &points.y[0], &points.xErr[0], &points.yErr[0])};
    gr.get()->SetMarkerSize(1.05);
    gr.get()->SetMarkerStyle(21);
    std::stringstream ss;
    ss << std::setprecision(2);
    ss << "stdAbs=" << calculateStdAbs(points) << ", " << "mOxy'=mOxy+k#bulletAd," << " k=" << par << ";mOxyGrad,%;Wr,%";
    h2dCorr.get()->SetTitle(ss.str().c_str());

    std::unique_ptr<TF1> f{std::make_unique<TF1>("f", "pol1", 0.75 * min, 1.25 * max)};
//    gr.get()->Fit(f.get(), "R");
    gr.get()->Draw("SAME P");

    std::unique_ptr<TLine> l{new TLine(0.75 * min,
                                       0.75 * min,
                                       1.25 * max,
                                       1.25 * max)};
    l.get()->Draw("SAME");


    std::vector<TLatex> labels;
    std::map<std::pair<std::string, Color_t>, Points> subPoints{
        { std::make_pair("N12", kRed), Points() },
        { std::make_pair("berez_blind_", kBlue), Points() },
        { std::make_pair("berez_11_w", kGreen), Points() },
        { std::make_pair("N12_blind", kOrange), Points() },
        { std::make_pair("barz_blind", kMagenta), Points() },
//        { std::make_pair("berez_", kYellow), Points() },
//        { std::make_pair("kuz_", kCyan), Points() },
        { std::make_pair("other", kBlack), Points() },
    };


    for (size_t i{0}; i < points.x.size(); ++i)
    {
        for (auto &item : subPoints)
        {
            if (points.l.at(i).find(item.first.first) != std::string::npos)
            {
                TMarker m{points.x.at(i), points.y.at(i), 21};
                m.SetMarkerSize(1.05);
                m.SetMarkerColor(item.first.second);
                m.DrawClone("SAME");
                item.second.l.push_back(points.l.at(i));
                item.second.x.push_back(points.x.at(i));
                item.second.y.push_back(points.y.at(i));
//                item.second.xErr.push_back(0.1);
//                item.second.yErr.push_back(0.5);

            }
        }
    }

    std::unique_ptr<TPaveText> pt{new TPaveText(0.1, 0.7, 0.3, 0.9, "NDC")};
    pt.get()->SetFillColor(0);
    pt.get()->SetBorderSize(1);

    // Add entries from map keys
    for (const auto& entry : subPoints) {
        const std::string& key = entry.first.first;
        Color_t color = entry.first.second;

        TText *text = pt.get()->AddText(key.c_str());
        text->SetTextColor(color);
    }
    pt.get()->Draw("SAME");
    c.get()->Print(psName.c_str());
    c.get()->Print((psName + ']').c_str());
}


int main()
{
     TVirtualFitter::SetDefaultFitter("Minuit");
    std::map<std::string, ChemResult> chem
    {

        { "pulp_rot_N12_1_", { 7.8, 4.2 } },
        { "pulp_rot_N12_2_", { 9.6, 5.5 } },
        { "pulp_rot_N12_3_", { 11.2, 6.2 } },
        { "pulp_rot_N12_4_", { 11.8, 3.9 } },
        { "pulp_rot_N12_5_", { 15.1, 7.9 } },
        { "pulp_rot_N12_6_", { 18.2, 4.9 } },
        { "pulp_rot_N12_7_", { 20.7, 6.7 } },
        { "pulp_rot_N12_8_", { 27.6, 8.0 } },
        { "pulp_rot_N12_9_", { 28.3, 7.8 } },
        { "pulp_rot_N12_10_", { 30.4, 8.2 } },
        { "pulp_rot_N12_11_", { 32.9, 8.1 } },

        {"pulp_rot_berez_2_", {15.6, 0.9}},
        {"pulp_rot_berez_7_", {19.8, 0.8}},
//        {"pulp_rot_berez_11_", {24.2, 0.8}},
//        {"pulp_rot_berez_5_", {17.5, std::nullopt}},
//        {"pulp_rot_berez_10_", {22.3, std::nullopt}},
        { "pulp_rot_barz_blind_a21p9_", {21.9, 1.6} },
        { "pulp_rot_barz_blind_a6p8_", {6.8, 1.4} },
        { "pulp_rot_berez_blind_a11p7_", {11.7, 3.0} },
        { "pulp_rot_berez_blind_a17p9_", {17.9, 2.8} },
        { "pulp_rot_N12_blind_a8p6_", {8.6, 2.4} },
        { "pulp_rot_N12_blind_a13p6_", {13.6, 2.1} },



//{ "pulp_rot_berez_6_w0p8_", {15.5, 0.8} },
//{ "pulp_rot_berez_6_w5_", {15.5, 5.0} },
//{ "pulp_rot_berez_6_w10_", {15.5, 10.0} },
//{ "pulp_rot_berez_6_w15_", {15.5, 15.0} },


//        {"pulp_rot_berez_11_w5_", {24.2, 5.0}},
//        {"pulp_rot_berez_11_w10_", {24.2, 10.0}},
//        {"pulp_rot_berez_11_w15_", {24.2, 15.0}},

        { "pulp_rot_berez_6_w0p8_", {15.5, 0.8 + 3.0} },
        { "pulp_rot_berez_6_w5_", {15.5, 5.0 + 3.0} },
        { "pulp_rot_berez_6_w10_", {15.5, 10.0 + 3.0} },
        { "pulp_rot_berez_6_w15_", {15.5, 15.0 + 3.0} },

        {"pulp_rot_berez_11_w5_", {24.2, 5.0 + 1.5}},
        {"pulp_rot_berez_11_w10_", {24.2, 10.0 + 1.5}},
        {"pulp_rot_berez_11_w15_", {24.2, 15.0 + 1.5}},

                { "pulp_rot_kuz_1_a7p85_", { 7.85, 0.5 } },
                { "pulp_rot_kuz_6_a18p48_", { 18.48, 0.5 } },
                { "pulp_rot_kuz_7_a27p12_", { 27.12, 0.5 } },
                { "pulp_rot_kuz_8_a5p89_", { 5.89, 0.5 } },
                { "pulp_rot_kuz_9_a5p76_", { 5.76, 0.5 } },

    };

    const std::map<int, std::string> columnElement
        {
         {1, "Al"},
         {3, "C"},
         {5, "N"},
         {7, "O"},
         {9, "Si"},
         };

//    const std::map<int, std::string> columnElement
//    {
//        {1, "Al"}, //0
//        {3, "C"},  //1
//        {5, "Ca"}, //2
//        {7, "Fe"}, //3
//        {9, "O"},  //4
//        {11, "Si"},//5
//    };


    const auto fileName{"rea.elts.stroy.txt"};
    std::cout << fileName << std::endl;

//    chem.insert(chemBlind.begin(), chemBlind.end());

    try
    {
//        std::regex m{"\\d+_\\d\\."};
        // std::regex m{"\\d+_(s|t)"};
        // std::regex m{"\\d+"};
//        std::regex m{"(pulp_rot_berez_6_w\\d+_\\d+)"};
//        std::regex m{"(check_bereza_8)"};
//        std::regex m{"(pulp_rot_N12_\\d+_\\d+)"};
//        std::regex m{R"(pulp_rot_berez_6_w\d+_sum|pulp_rot_berez_6_w\d+p\d+_sum)"};
//        std::regex m{R"(pulp_rot_berez_11_w\d+_sum|pulp_rot_berez_11_w\d+p\d+_sum)"};
//        std::regex m{R"(pulp_rot_berez_6_w\d+_sum|pulp_rot_berez_6_w\d+p\d+_sum|pulp_rot_berez_11_w\d+_sum|pulp_rot_berez_11_w\d+p\d+_sum)"};
//        std::regex m{R"(pulp_rot_berez_6_w\d+_sum|pulp_rot_berez_6_w\d+p\d+_sum|pulp_rot_berez_11_w\d+_sum|pulp_rot_berez_11_w\d+p\d+_sum|pulp_rot_N12_\d+_sum)"};
//        std::regex m{R"(pulp_rot_berez_11_w\d+_sum|pulp_rot_berez_11_w\d+p\d+_sum|pulp_rot_N12_\d+_sum)"}; // used to get optimal par
        std::regex m{R"(pulp_rot_berez_11_w\d+_sum|pulp_rot_berez_11_w\d+p\d+_sum|pulp_rot_N12_\d+_sum|pulp_rot_berez_6_w\d+_sum|pulp_rot_berez_6_w\d+p\d+_sum)"};
//        std::regex m{"(pulp_rot_berez_\\d+_\\d+)"};
//        std::regex m{"(pulp_rot_berez_111_w5_|pulp_rot_berez_111_w10_|pulp_rot_berez_111_w15_)"};
//        std::regex m{"(pulp_rot_N12_\\d+_\\d+|pulp_rot_berez_\\d+_\\d+|_blind_a\\d+p\\d+_\\d+)"};
//        std::regex m{"(pulp_rot_N12_\\d+_\\d+|pulp_rot_berez_\\d+_\\d+)"};
//        std::regex m{"(pulp_rot_N12_\\d+_\\d+)"};
//        std::regex m{"(pulp_rot_N12_\\d*[02468]_\\d+)"};

//        std::regex m{R"(pulp_rot_berez_11_w\d+_sum|pulp_rot_berez_11_w\d+p\d+_sum|pulp_rot_N12_\d+_sum|_blind_a\d+p\d+_sum|pulp_rot_berez_[27]_sum|pulp_rot_kuz_\d+_a\d+p\d+_sum)"};

//        std::regex m{R"(pulp_rot_berez_11_w\d+_sum|pulp_rot_berez_11_w\d+p\d+_sum|pulp_rot_N12_\d+_sum|_blind_a\d+p\d+_sum)"};
//        std::regex m{R"(pulp_rot_berez_11_w\d+_\d+)"};


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

        std::cout << data1.size() << std::endl;

        auto par{getOptimalPar(data1)};

        drawCorrGraphWr(data1, par);
//        drawCorrGraphWr1(data1, -0.39);

//        drawCorrGraphs(data1);


        return 0;

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
        std::unique_ptr<TF1> f{new TF1("f", fObj, points.x.front(), points.x.back(), 8)};


        f.get()->SetParameter(0, -1.0);
        f.get()->SetParLimits(0, -2.0, 0.5);
        f.get()->SetParameter(1, 100.0);
        f.get()->SetParLimits(1, 50.0, 150.0);
        f.get()->FixParameter(2, 1.0);


//        f.get()->SetParameter(2, 1.0);
//        f.get()->SetParLimits(2, 0.0, 2.0);
        for (auto pIdx{3}; pIdx < 7; ++pIdx) {
            f.get()->SetParameter(pIdx, 1.0);
            f.get()->SetParLimits(pIdx, 0.0, 100.0);
        }
        f.get()->SetParameter(7, 0.0);
        f.get()->SetParLimits(7, -50.0, 50.0);

//        f.get()->SetParameter(0, 1.0);
//        f.get()->SetParameter(1, 0.0);
//        f.get()->SetParameter(2, 1.0);
//        f.get()->SetParameter(3, 1.0);
//        f.get()->SetParameter(4, 1.0);
//        f.get()->SetParameter(5, 2.0);
//        f.get()->SetParameter(6, 1.0);
//        f.get()->SetParameter(7, 1.0);
//        f.get()->SetParameter(8, 0.0);

//        f.get()->FixParameter(2, 1.0);

//        f.get()->SetParLimits(3, 0.0, 10.0);
//        f.get()->SetParLimits(3, 0.3, 2.0);
//        f.get()->SetParLimits(4, 0.0, 20.0);

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
                { std::make_pair("tochka_", kRed), Points() },
                { std::make_pair("sector_", kBlue), Points() },
                { std::make_pair("field_", kGreen), Points() },
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
            item.DrawClone("SAME");
        }

        c.get()->Print(psName.c_str());
        c.get()->Print((psName + ']').c_str());
        c.get()->Close();


        const auto fileName_1{"rea.elts.stroy.txt"};
//        std::regex s{"sum"};
        // std::regex s{"\\d+_s"};
//        std::regex s{"check_bereza"};
//        std::regex s{"(pulp_rot_N12_\\d+_\\d+)"};
//        std::regex s{"(pulp_rot_N12_\\d+_sum)"};
//        std::regex s{"(pulp_rot_N12_[13579]_sum|pulp_rot_N12_[1-9]\\d*[13579]_sum)"};
//        std::regex s{"(pulp_rot_N12_\\d+_sum|pulp_rot_berez_\\d+_sum|_blind_a\\d+p\\d+_sum)"};
//        std::regex s{"(pulp_rot_berez_111_w5_|pulp_rot_berez_111_w10_|pulp_rot_berez_111_w15_)"};
//        std::regex s{"(pulp_rot_berez_6_w\\d+_\\d+)"};
//        std::regex s{"(pulp_rot_N12_\\d*[02468]_\\d+)"};

        std::regex s{R"(pulp_rot_berez_11_w\d+_sum|pulp_rot_berez_11_w\d+p\d+_sum|pulp_rot_N12_\d+_sum|_blind_a\d+p\d+_sum)"};
//        std::regex s{R"(pulp_rot_berez_11_w\d+_sum|pulp_rot_berez_11_w\d+p\d+_sum|pulp_rot_N12_\d+_sum|_blind_a\d+p\d+_sum)"};

        auto data1Sum{getFitResults(fileName_1, columnElement, chem, s)};
        calcConv(data1Sum, f, value);

//        std::regex t{"_\\d+_"};
//        auto dataBlindSum{getFitResults(fileName_1, columnElement, chemBlind, t)};
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

//                for (auto &i : fR) {
//                    auto itPreGrad = preGrad.find(i.e);
//                    if (itPreGrad != preGrad.end()) {
//                        i.value = itPreGrad->second.first +  itPreGrad->second.second * i.value;
//                    }
//                }


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
                points.yErr.push_back(0.25 * v.value());
//                auto yE{0.0};
//                if (value == Data1::Value::A) {
//                    yE = 0.1 * v.value();
//                }
//                if (value == Data1::Value::W) {
//                    yE = 0.03 * v.value();
//                }
//                points.yErr.push_back(yE);
            }
        }

    }
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
                mmn[xx].push_back(it->second.fr.at(i).at(static_cast<size_t>(0)).value); // Al
                mmn[xx].push_back(it->second.fr.at(i).at(static_cast<size_t>(1)).value); // C
                mmn[xx].push_back(it->second.fr.at(i).at(static_cast<size_t>(2)).value); // N
                mmn[xx].push_back(it->second.fr.at(i).at(static_cast<size_t>(3)).value); // O
                mmn[xx].push_back(it->second.fr.at(i).at(static_cast<size_t>(4)).value); // Si
                xx++;
            }
        }
    }
}

void useSub(const Points &points,
            std::string &str,
            const Data1::Value value,
            const bool isLabels = false)
{
    std::map<std::pair<std::string, Color_t>, Points> subPoints{
        { std::make_pair("N12", kRed), Points() },
        { std::make_pair("berez", kBlue), Points() },
        { std::make_pair("p43", kGreen), Points() },
        { std::make_pair("blind", kOrange), Points() },
        { std::make_pair("raspad", kCyan), Points() },
        { std::make_pair("other", kBlack), Points() },
    };

    for (size_t i{0}; i < points.x.size(); ++i)
    {
        auto isOther{true};
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

//                isOther = false;
            }
        }
        if (isOther)
        {
            TMarker m{points.x.at(i), points.y.at(i), 21};
            m.SetMarkerSize(1.5);
            m.SetMarkerColor(kBlack);
//            m.DrawClone("SAME");
            subPoints.at({"other", kBlack}).l.push_back(points.l.at(i));
            subPoints.at({"other", kBlack}).x.push_back(points.x.at(i));
            subPoints.at({"other", kBlack}).y.push_back(points.y.at(i));
            subPoints.at({"other", kBlack}).xErr.push_back(0.1);
            subPoints.at({"other", kBlack}).yErr.push_back(0.5);
        }
    }
    std::stringstream ss;
    ss.str("");ss.clear();
    ss << (value == Data1::Value::A ? "Ad" : "Wr");
    ss << ": stdAbs=" << std::setprecision(3);
    auto stdAbs1 = [](const Points &points){
        std::vector<double> d2;
        for (size_t i{0}; i < points.x.size(); ++i)
        {
            d2.push_back(std::pow(points.y.at(i) - points.x.at(i), 2));
        }
        return std::sqrt(std::accumulate(d2.begin(), d2.end(), 0.0) / static_cast<double>(d2.size()));
    };
    for (auto item : subPoints)
    {
        if (!std::isnan(stdAbs1(item.second)))
        {
            ss << "[#color[" << static_cast<int>(item.first.second) << "]{" << stdAbs1(item.second) << "}] ";
        }
    }

    ss << ";AGP-K, %;Chem, %";
    str = ss.str();

    if (isLabels)
    {
        std::vector<TLatex> labels;
        for (size_t i{0}; i < points.x.size(); ++i)
        {
            TLatex l(points.x.at(i), points.y.at(i) + 1.25 * points.yErr.at(i), points.l.at(i).c_str());
            l.SetTextAngle(90);
            l.SetTextAlign(12);
            l.SetTextSize(0.02);
            labels.push_back(l);
        }
        for (const auto &item : labels)
        {
            item.DrawClone("SAME");
        }
    }
}

void addSubPoints(const Points &points, std::map<std::pair<std::string, Color_t>, Points> &subPoints)
{
    for (size_t i{0}; i < points.x.size(); ++i)
    {
        auto isOther{true};
        for (auto &item : subPoints)
        {
            if (points.l.at(i).find(item.first.first) != std::string::npos)
            {
//                TMarker m{points.x.at(i), points.y.at(i), 21};
//                m.SetMarkerSize(1.5);
//                m.SetMarkerColor(item.first.second);
//                m.DrawClone("SAME");
                item.second.l.push_back(points.l.at(i));
                item.second.x.push_back(points.x.at(i));
                item.second.y.push_back(points.y.at(i));
                item.second.xErr.push_back(0.1);
                item.second.yErr.push_back(0.5);

                isOther = false;
            }
        }
        if (isOther)
        {
//            TMarker m{points.x.at(i), points.y.at(i), 21};
//            m.SetMarkerSize(1.5);
//            m.SetMarkerColor(kBlack);
//            m.DrawClone("SAME");
            subPoints.at({"other", kBlack}).l.push_back(points.l.at(i));
            subPoints.at({"other", kBlack}).x.push_back(points.x.at(i));
            subPoints.at({"other", kBlack}).y.push_back(points.y.at(i));
            subPoints.at({"other", kBlack}).xErr.push_back(0.1);
            subPoints.at({"other", kBlack}).yErr.push_back(0.5);
        }
    }
}

void useSub1(const Points &points,
            std::string &str,
            const Data1::Value value,
            const bool isLabels = false)
{


    std::map<std::pair<std::string, Color_t>, Points> subPoints{
        { std::make_pair("N12", kRed), Points() },
        { std::make_pair("berez", kBlue), Points() },
        { std::make_pair("p43", kGreen), Points() },
        { std::make_pair("blind", kOrange), Points() },
        { std::make_pair("raspad", kCyan), Points() },
        { std::make_pair("other", kBlack), Points() },
    };

    addSubPoints(points, subPoints);
    std::stringstream ss;
    ss.str("");ss.clear();
    ss << (value == Data1::Value::A ? "Ad" : "Wr");
    ss << ": stdAbs=" << std::setprecision(3);
    auto stdAbs1 = [](const Points &points){
        std::vector<double> d2;
        for (size_t i{0}; i < points.x.size(); ++i)
        {
            d2.push_back(std::pow(points.y.at(i) - points.x.at(i), 2));
        }
        return std::sqrt(std::accumulate(d2.begin(), d2.end(), 0.0) / static_cast<double>(d2.size()));
    };
    for (auto item : subPoints)
    {
//        for (size_t i{0}; i < item.second.x.size(); ++i)
//        {
//            TMarker m{item.second.x.at(i), item.second.y.at(i), 21};
//            m.SetMarkerSize(1.5);
//            m.SetMarkerColor(item.first.second);
//            m.DrawClone("SAME");
//        }
        if (!std::isnan(stdAbs1(item.second)))
        {
            ss << "[#color[" << static_cast<int>(item.first.second) << "]{" << stdAbs1(item.second) << "}] ";
        }
    }

    ss << ";Number;Chem/AGP-K, %";
    str = ss.str();

    std::map<std::pair<std::string, Color_t>, Points> subPoints1{
        { std::make_pair("N12", kRed), Points() },
        { std::make_pair("berez", kBlue), Points() },
        { std::make_pair("p43", kGreen), Points() },
        { std::make_pair("111", kOrange), Points() },
        { std::make_pair("raspad", kCyan), Points() },
        { std::make_pair("other", kBlack), Points() },
    };



    Points pointsC{points};
    Points pointsR{points};
    for (size_t i{0}; i < points.x.size(); ++i)
    {
        pointsC.x.at(i) = static_cast<double>(i + 1);
        pointsC.yErr.at(i) = 0.0;
        pointsR.x.at(i) = static_cast<double>(i + 1);
        pointsR.y.at(i) = points.x.at(i);
        pointsR.yErr.at(i) = points.xErr.at(i);
    }

    addSubPoints(pointsR, subPoints1);

    writePointsToFile("outputR.txt", pointsR);
    writePointsToFile("outputC.txt", pointsC);

    std::unique_ptr<TGraphErrors> grC{new TGraphErrors(static_cast<int>(pointsC.x.size()),
                                                       &pointsC.x[0],
                                                       &pointsC.y[0],
                                                       &pointsC.xErr[0],
                                                       &pointsC.yErr[0])};
    grC.get()->SetLineColor(kRed);
    grC.get()->SetLineWidth(3);
    grC.get()->SetMarkerSize(1.5);
    grC.get()->SetMarkerStyle(21);

    std::unique_ptr<TGraphErrors> grR{new TGraphErrors(static_cast<int>(pointsR.x.size()),
                                                       &pointsR.x[0],
                                                       &pointsR.y[0],
                                                       &pointsR.xErr[0],
                                                       &pointsR.yErr[0])};
    grR.get()->SetMarkerSize(1.5);
    grR.get()->SetMarkerStyle(21);


    grC.get()->DrawClone("SAME L");
    grR.get()->DrawClone("SAME P");

    for (auto item : subPoints1)
    {
        for (size_t i{0}; i < item.second.x.size(); ++i)
        {
            TMarker m{item.second.x.at(i), item.second.y.at(i), 21};
            m.SetMarkerSize(1.5);
            m.SetMarkerColor(item.first.second);
            m.DrawClone("SAME");
        }
    }

    if (isLabels)
    {
        std::vector<TLatex> labels;
        for (size_t i{0}; i < pointsR.x.size(); ++i)
        {
            TLatex l(pointsR.x.at(i), pointsR.y.at(i) + 1.25 * pointsR.yErr.at(i), pointsR.l.at(i).c_str());
            l.SetTextAngle(90);
            l.SetTextAlign(12);
            l.SetTextSize(0.02);
            labels.push_back(l);
        }
        for (const auto &item : labels)
        {
            item.DrawClone("SAME");
        }
    }
}

void writePointsToFile(const std::string fileName, const Points &points)
{
    std::ofstream ofs;
    ofs.open(fileName, std::ios::out);
    for (size_t i{0}; i < points.l.size(); ++i)
    {
        ofs << points.x.at(i) << " " << points.y.at(i)  << " " << points.l.at(i) << std::endl;
    }
    ofs.close();
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
    auto p6{f.get()->GetParameter(6)};
    auto p7{f.get()->GetParameter(7)};

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
                auto Al{it->second.fr.at(i).at(static_cast<size_t>(0)).value};
                auto C{it->second.fr.at(i).at(static_cast<size_t>(1)).value};
                auto N{it->second.fr.at(i).at(static_cast<size_t>(2)).value};
                auto O{it->second.fr.at(i).at(static_cast<size_t>(3)).value};
                auto Si{it->second.fr.at(i).at(static_cast<size_t>(4)).value};

                auto res{0.0};
                auto w{
                    p3 * O
                    - p4 * Si
                    - p5 * Al
                    - p6 * N
                    + p7
                };
                auto a{
                    (p1 + p0 * C)
                    / (100.0 - p2 * w) * 100.0
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

    std::ofstream ofs("output_data.csv", std::ios_base::app);

    std::vector<TLatex> labels;
    for (size_t i{0}; i < points.x.size(); ++i)
    {

        auto pos{points.l.at(i).find_last_of("_")};
        auto text{points.l.at(i)};
        if (pos != std::string::npos && pos == points.l.at(i).length() - 1)
        {
            text = text.substr(0, pos);
        }
        std::cout << i << " " << text << std::endl;
        TLatex l(points.x.at(i), points.y.at(i) + 1.25 * points.yErr.at(i), text.c_str());
        l.SetTextAngle(90);
        l.SetTextAlign(12);
        l.SetTextSize(0.02);
        labels.push_back(l);

        if (ofs.is_open()) {
            ofs << text << " " << points.x.at(i) << " " << points.xErr.at(i) << " "
                << points.y.at(i) << " " << points.yErr.at(i) << std::endl;
        }
    }

    ofs.close();


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
//    gr.get()->Draw("P");

    lConv.get()->Draw("SAME");
    std::vector<std::string> uniqueL{points.l};
    std::sort(uniqueL.begin(), uniqueL.end());
    auto uniqueIt{std::unique(uniqueL.begin(), uniqueL.end())};
    uniqueL.erase(uniqueIt, uniqueL.end());
    Points avgPoints;
    avgPoints.l = uniqueL;



    for (size_t i{0}; i < avgPoints.l.size(); ++i)
    {
        std::vector<double> tmpX;
        std::vector<double> tmpXerr;
        std::vector<double> tmpY;
        std::vector<double> tmpYerr;
        auto it{points.l.begin()};
        while ((it = std::find(it, points.l.end(), avgPoints.l.at(i))) != points.l.end())
        {
            auto idx{std::distance(points.l.begin(), it)};
            tmpX.push_back(points.x.at(static_cast<size_t>(idx)));
            tmpXerr.push_back(points.xErr.at(static_cast<size_t>(idx)));
            tmpY.push_back(points.y.at(static_cast<size_t>(idx)));
            tmpYerr.push_back(points.yErr.at(static_cast<size_t>(idx)));
            it++;
        }
        avgPoints.x.push_back(std::accumulate(tmpX.begin(), tmpX.end(), 0.0) / static_cast<double>(tmpX.size()));
        avgPoints.xErr.push_back(std::accumulate(tmpXerr.begin(), tmpXerr.end(), 0.0) / static_cast<double>(tmpXerr.size()));
        avgPoints.y.push_back(std::accumulate(tmpY.begin(), tmpY.end(), 0.0) / static_cast<double>(tmpY.size()));
        avgPoints.yErr.push_back(std::accumulate(tmpYerr.begin(), tmpYerr.end(), 0.0) / static_cast<double>(tmpYerr.size()));
    }

    writePointsToFile("output.txt", avgPoints);

    std::string str;
    useSub(avgPoints, str, value, true);
    h2dConv.get()->SetTitle(str.c_str());

    c.get()->Print(psName.c_str());

    std::unique_ptr<TH2D> h2dConv1{new TH2D("h2dConv1",
                                           ss.str().c_str(),
                                           static_cast<int>(avgPoints.x.size()),
                                           0.5,
                                           0.5 + static_cast<int>(avgPoints.x.size()),
                                           static_cast<int>(points.y.size()),
                                           0.75 * (*std::min_element(points.y.begin(), points.y.end())),
                                           1.25 * (*std::max_element(points.y.begin(), points.y.end())))};
    h2dConv1->SetStats(0);
    h2dConv1.get()->Draw();
    std::string str1;
    useSub1(avgPoints, str1, value, true);
    h2dConv1.get()->SetTitle(str1.c_str());
    c.get()->Print(psName.c_str());
    c.get()->Print((psName + ']').c_str());
    c.get()->Close();



//    auto useSub{true};
//    if (useSub)
//    {

//        std::map<std::pair<std::string, Color_t>, Points> subPoints{
//            { std::make_pair("b8", kRed), Points() },
//            { std::make_pair("check_w", kBlue), Points() },
//            { std::make_pair("p43", kGreen), Points() },
//            { std::make_pair("bereza_8", kOrange), Points() },
//            { std::make_pair("raspad", kCyan), Points() },
//            { std::make_pair("other", kBlack), Points() },
//        };


//        for (size_t i{0}; i < points.x.size(); ++i)
//        {
//            auto isOther{false};
//            for (auto &item : subPoints)
//            {
//                if (points.l.at(i).find(item.first.first) != std::string::npos)
//                {
//                    TMarker m{points.x.at(i), points.y.at(i), 21};
//                    m.SetMarkerSize(1.5);
//                    m.SetMarkerColor(item.first.second);
//                    m.DrawClone("SAME");
//                    item.second.l.push_back(points.l.at(i));
//                    item.second.x.push_back(points.x.at(i));
//                    item.second.y.push_back(points.y.at(i));
//                    item.second.xErr.push_back(0.1);
//                    item.second.yErr.push_back(0.5);

//                   isOther = true;
//                }
//            }
//            if (!isOther)
//            {
//                subPoints.at({"other", kMagenta}).l.push_back(points.l.at(i));
//                subPoints.at({"other", kMagenta}).x.push_back(points.x.at(i));
//                subPoints.at({"other", kMagenta}).y.push_back(points.y.at(i));
//                subPoints.at({"other", kMagenta}).xErr.push_back(0.1);
//                subPoints.at({"other", kMagenta}).yErr.push_back(0.5);
//            }
//        }
//        std::stringstream ss;
//        ss.str("");ss.clear();
//        ss << (value == Data1::Value::A ? "Ad" : "Wr");
//        ss << ": stdAbs=" << std::setprecision(3);
//        auto stdAbs1 = [](Points points){
//            std::vector<double> d2;
//            for (size_t i{0}; i < points.x.size(); ++i)
//            {
//                d2.push_back(std::pow(points.y.at(i) - points.x.at(i), 2));
//            }
//            return std::sqrt(std::accumulate(d2.begin(), d2.end(), 0.0) / static_cast<double>(d2.size()));
//        };
////        auto avg1 = [](Points points){
////            return std::accumulate(points.x.begin(), points.x.end(), 0.0) / points.x.size();
////        };
//        for (auto item : subPoints)
//        {
//            ss << "[#color[" << static_cast<int>(item.first.second) << "]{" << stdAbs1(item.second) << "}] ";
//        }

//        ss << ";AGP-K, %;Chem, %";
//        h2dConv.get()->SetTitle(ss.str().c_str());
//    }

//    lConv.get()->Draw("SAME");
//    for (const auto &item : labels)
//    {
//        item.DrawClone("SAME");
//    }
//    c.get()->Print(psName.c_str());
//    c.get()->Print((psName + ']').c_str());
//    c.get()->Close();
}

struct Data { // TODO need this for calcRep
    enum class Value {
        A,
        W
    };
    std::optional<double> a;
    std::optional<double> w;
};


void calcRep(const std::shared_ptr<TF1> &f,
             const std::map<std::string, Data> &chem,
             const std::map<std::string, std::vector<std::string>> &vC,
             const std::map<std::string, std::vector<std::string>> &vO) // TODO make to work with Data1
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
