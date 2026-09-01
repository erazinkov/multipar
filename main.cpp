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
#include <TStyle.h>
#include <TPaveText.h>

#include <bits/stdc++.h>

#include <regex>

#include "structs.h"

#include "data.h"

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


std::map<std::string, Data> getData(const std::string &fileName,
                          const std::map<int, std::string> &columnElement,
                          const std::map<std::string, ChemResult> &chem,
                          const std::regex &pattern);


class FitFunction
{
public:
    FitFunction(const std::map<int, std::vector<double>> &d, const std::map<int, double> &dA)
        : d_{d}, dA_{dA} {}

    double operator() (double *x, double *par)
    {
        double arg{x[0]};
        int idx{ std::min(static_cast<int>(std::round(arg)), static_cast<int>(d_.size() - 1)) };
        auto val{0.0};
        std::map<std::string, double> dm{
            {"Al", d_.at(idx).at(0)},
            {"C", d_.at(idx).at(1)},
            {"N", d_.at(idx).at(2)},
            {"O", d_.at(idx).at(3)},
            {"Si", d_.at(idx).at(4)},
//            {"A", dA_.at(idx - dA_.size())}
        };

        if (idx < static_cast<int>(dA_.size())) {
            val = ( par[3]
                  - par[4] * par[0] * dm.at("O")
                  - par[2] * par[4]
                  - par[5] * dm.at("C")
                  - par[6] * dm.at("N") )
                  / ( 1.0 - par[1] * par[4] );
        } else {
            val = ( par[0] * dm.at("O")
                   - par[1] * dA_.at(idx - dA_.size())
                   + par[2] );
        }
        return val;
   }
private:
    const std::map<int, std::vector<double>> d_;
    const std::map<int, double> dA_;
};

class FitFunction1
{
public:
    FitFunction1(const std::vector<Point> &points, const size_t s)
        : points_{points}, s_{s} {}

    double operator() (double *x, double *par)
    {
        double arg{x[0]};
        int idx{ std::min(static_cast<int>(std::round(arg)), static_cast<int>(points_.size() - 1)) };
        auto val{0.0};
        if (idx < static_cast<int>(s_)) {
            val = ( par[3]
                  - par[4] * par[0] * getElementResultValue_("O", points_.at(idx)).value_or(0.0)
                  - par[2] * par[4]
                  - par[5] * getElementResultValue_("C", points_.at(idx)).value_or(0.0)
                  - par[6] * getElementResultValue_("N", points_.at(idx)).value_or(0.0) )
                  / ( 1.0 - par[1] * par[4] );
        } else {
            val = ( par[0] * getElementResultValue_("O", points_.at(idx)).value_or(0.0)
                   - par[1] * points_.at(idx - s_).chemResult.a.value()
                   + par[2] );
        }
        return val;
   }

private:
    std::optional<double> getElementResultValue_(const std::string &element, const Point &point) {
        for (const auto& e : point.fitResult.elementResults) {
            if (e.name == element) {
                return e.value;
            }
        }
        return std::nullopt;
    }
    const std::vector<Point> points_;
    const size_t s_;
};



struct Points {
    std::vector<std::string> l;
    std::vector<double> x;
    std::vector<double> xErr;
    std::vector<double> y;
    std::vector<double> yErr;
    std::vector<double> d;
};

void writePointsToFile(const std::string fileName, const Points &points);


double calculateStdAbsCon(const std::vector<Point> &points) {
    if (points.empty()) {
        return 0.0;
    }
    double sumSquaredDiff = 0.0;
    for (size_t i = 0; i < points.size(); ++i) {
        const double diff = points.at(i).y - points.at(i).x;
        sumSquaredDiff += diff * diff;
    }
    return std::sqrt(sumSquaredDiff / points.size());
}

double calculateStdAbsCon(const Points& points) {
    if (points.x.empty()) return 0.0;

    double sumSquaredDiff = 0.0;
    for (size_t i = 0; i < points.x.size(); ++i) {
        const double diff = points.y[i] - points.x[i];
        sumSquaredDiff += diff * diff;
    }
    return std::sqrt(sumSquaredDiff / points.x.size());
}

double calculateAvg(const std::vector<double>& values) {
    if (values.empty()) return 0.0;
    return std::accumulate(values.begin(), values.end(), 0.0) / values.size();
}

double calculateStdAbsRep(const std::vector<double>& values) {
    if (values.empty()) return 0.0;

    double sumSquaredDiff = 0.0;
    auto avg = calculateAvg(values);
    for (size_t i = 0; i < values.size() - 1; ++i) {
        const double diff = values[i + 1] - avg;
        sumSquaredDiff += diff * diff;
    }
    return std::sqrt(sumSquaredDiff / values.size());
}

/**
* @brief Извлекает точки данных (Point) из мапы с данными по указанному типу химического результата.
*
* Функция перебирает все записи в data, для каждой записи обходит все fitResults,
* и для каждого fitResult проверяет наличие значения указанного типа (A или W).
* Если значение найдено, создаётся Point с соответствующими данными.
*
* @param data Ссылка на константную мапу, где ключ - идентификатор образца (std::string),
* значение - структура Data, содержащая fitResults и chemResult.
* @param type Тип химического результата (ChemResult::Type), определяющий,
* какое поле использовать: A (зольность) или W (влажность).
*
* @return std::vector<Point> Вектор объектов Point, каждый из которых содержит:
* - sample: идентификатор образца (ключ из data)
* - chemResult: полная структура химического результата
* - x: порядковый номер точки в векторе (начиная с 0)
* - xErr: всегда 0.0 (погрешность по оси X отсутствует)
* - y: значение выбранного параметра (A или W)
* - yErr: погрешность по оси Y (рассчитывается как err * value)
* - fitResult: соответствующий результат фитирования
*
* @note Погрешность по Y вычисляется по формуле: yErr = err * value,
* где err - константа для данного типа (0.1 для A, 0.03 для W).
* @note Порядковый номер точки (x) определяется как текущий размер вектора,
* поэтому нумерация начинается с 0 и идёт по порядку добавления.
* @warning Если значение для указанного типа отсутствует (std::nullopt),
* точка не создаётся и пропускается.
* @see Data, Point, ChemResult::Type
*/
std::vector<Point> getPointsByType(const std::map<std::string, Data> &data, const ChemResult::Type &type) {
    std::vector<Point> points;
    for (const auto &[key, value] : data) {
        for (const auto &fr : value.fitResults) {
            std::optional<double> v;
            double err{0.0};
            switch (type) {
            case (ChemResult::Type::A):
                v = value.chemResult.a;
                err = 0.1;
                break;
            case (ChemResult::Type::W):
                v = value.chemResult.w;
                err = 0.03;
                break;
            }
            if (v.has_value()) {
                Point point;
                point.sample = key;
                point.chemResult = value.chemResult;
                point.x = points.size();
                point.xErr = 0.0;
                point.y = v.value();
                point.yErr = err * v.value();
                point.fitResult = fr;
                points.push_back(point);
            }
        }
    }
    return points;
};

int main()
{
    TVirtualFitter::SetDefaultFitter("Minuit");

    std::map<std::string, ChemResult> chem = data_chem;

    const std::map<int, std::string> columnElement
    {
         {1, "Al"},
         {3, "C"},
         {5, "N"},
         {7, "O"},
         {9, "Si"},
    };


    const auto fileName{"rea.elts.stroy.work.flow.txt"};
    std::cout << fileName << std::endl;

    try {
        std::regex m{R"(\bsample([1-9]|[12][0-9]|30)\b)"};

        auto data{getData(fileName, columnElement, chem, m)};
        for (const auto &[key, value] : data) {
            std::cout << key << " ";
            value.print();
            std::cout << std::endl;
        }

        std::vector<Point> points_a{getPointsByType(data, ChemResult::Type::A)};
        std::vector<Point> points_w{getPointsByType(data, ChemResult::Type::W)};
        for (auto &p : points_w) {
            p.x  = p.x + points_a.size();
        }

        std::vector<Point> points;
        points.insert(points.end(), points_a.cbegin(), points_a.cend());
        points.insert(points.end(), points_w.cbegin(), points_w.cend());

        // std::vector<double> xx, xxErr;
        // std::vector<double> yy, yyErr;

        for (const auto &p : points) {
            // xx.push_back(p.x);
            // yy.push_back(p.y);
            // xxErr.push_back(p.xErr);
            // yyErr.push_back(p.yErr);
            std::cout << p.sample << " " << p.x << " " << p.y;
            p.fitResult.print();
            std::cout << std::endl;
        }

        FitFunction1 fObj(points, points_a.size());
        std::unique_ptr<TF1> f{new TF1("f", fObj, points.front().x, points.back().x, 7)};
        const std::vector<double> parameters = {
              1.55468e+00,
              7.55269e-01,
             -1.22767e+01,
              1.05060e+02,
              6.69100e-01,
              1.17662e+00,
              0.00000e+00
        };

        auto setInitialParameters = [&parameters](TF1 *f){
            for (auto it{parameters.begin()}; it != parameters.end(); it++) {
                f->SetParameter(std::distance(parameters.begin(), it), *it);
            }
        };
        setInitialParameters(f.get());

        std::unique_ptr<TGraphErrors> gr{new TGraphErrors(points.size())};


        for (size_t i{0}; i < points.size(); i++) {
            gr.get()->SetPoint(i, points.at(i).x, points.at(i).y);
            gr.get()->SetPointError(i, points.at(i).xErr, points.at(i).yErr);
        }

        gr.get()->SetMarkerSize(1.5);
        gr.get()->SetMarkerStyle(21);
        gr.get()->SetTitle(";N_{probe};[...A, ...W]");

        f.get()->SetNpx(10 * static_cast<int>(points.size()));

        gr.get()->Fit(f.get(), "R");

        TMarker bM{points_a.back().x + 0.5, TMath::MinElement(gr.get()->GetN(), gr.get()->GetY()) + 5.0, 20};
        bM.SetMarkerSize(1.5);
        bM.SetMarkerColor(kGreen);
        const std::string psName{"output.ps"};
        std::unique_ptr<TCanvas> c{new TCanvas("c", "c", 1024, 960)};
        gPad->SetGrid();
        c.get()->Print((psName + '[').c_str());
        gr.get()->Draw("APL");

        bM.DrawClone("SAME");
        c.get()->Print(psName.c_str());
        c.get()->Print((psName + ']').c_str());
        c.get()->Close();

        auto getPredicatedValueByType = [](const FitResult &fr, const ChemResult::Type &type, TF1 *f){
            auto val_a{0.0};
            auto val_w{0.0};
            val_a = ( f->GetParameter(3)
                     - f->GetParameter(4) * f->GetParameter(0) * fr.getElementResultByName("O").value
                     - f->GetParameter(2) * f->GetParameter(4)
                     - f->GetParameter(5) * fr.getElementResultByName("C").value
                     - f->GetParameter(6) * fr.getElementResultByName("N").value )
                    / ( 1.0 - f->GetParameter(1) * f->GetParameter(4) );
            val_w = ( f->GetParameter(0) * fr.getElementResultByName("O").value
                     - f->GetParameter(1) * val_a
                     + f->GetParameter(2) );

            auto val{0.0};
            switch (type) {
            case (ChemResult::Type::A):
                val = val_a;
                break;
            case (ChemResult::Type::W):
                val = val_w;
                break;
            }
            return val;
        };


        auto getPredicatedPointsByType = [&](const std::map<std::string, Data> &data, const ChemResult::Type &type, TF1 *f) {
            std::vector<Point> points;
            for (const auto &[key, value] : data) {
                for (const auto &fr : value.fitResults) {
                    std::optional<double> v;
                    double err{0.0};
                    switch (type) {
                    case (ChemResult::Type::A):
                        v = value.chemResult.a;
                        err = 0.1;
                        break;
                    case (ChemResult::Type::W):
                        v = value.chemResult.w;
                        err = 0.03;
                        break;
                    }
                    if (v.has_value()) {
                        Point point;
                        point.sample = key;
                        point.chemResult = value.chemResult;
                        point.x = getPredicatedValueByType(fr, type, f);
                        point.xErr = getPredicatedValueByType(fr, type, f) * 0.01;
                        point.y = v.value();
                        point.yErr = err * v.value();
                        point.fitResult = fr;
                        points.push_back(point);
                    }
                }
            }
            return points;
        };

        // std::regex m_a{R"(pulp_rot_berez_7_w\d+_sum|pulp_rot_berez_2_w\d+_sum|pulp_rot_berez_11_w\d+_sum|pulp_rot_berez_7_w\d+p\d+_sum)"};
        std::regex m_a{R"(sample\d+)"};

        auto data_a{getData(fileName, columnElement, chem, m_a)};

        std::vector<Point> points_p_a{getPredicatedPointsByType(data_a, ChemResult::Type::A, f.get())};
        std::vector<Point> points_p_w{getPredicatedPointsByType(data_a, ChemResult::Type::W, f.get())};

        std::vector<Point> points_p;
        // points_p.insert(points_p.end(), points_p_a.cbegin(), points_p_a.cend());
        points_p.insert(points_p.end(), points_p_w.cbegin(), points_p_w.cend());

        std::unique_ptr<TGraphErrors> gr_p{new TGraphErrors(points_p.size())};


        for (size_t i{0}; i < points_p.size(); i++) {
            gr_p.get()->SetPoint(i, points_p.at(i).x, points_p.at(i).y);
            gr_p.get()->SetPointError(i, points_p.at(i).xErr, points_p.at(i).yErr);
        }

        gr_p.get()->SetMarkerSize(1.5);
        gr_p.get()->SetMarkerStyle(21);

        auto itMin = std::max_element(points_p.begin(), points_p.end(),
                                      [](const Point& a, const Point& b) {
                                    return a.y > b.y;
                                      });
        auto itMax = std::max_element(points_p.begin(), points_p.end(),
                                   [](const Point& a, const Point& b) {
                                    return a.y < b.y;
                                   });
        auto min{(*itMin).y};
        auto max{(*itMax).y};

        std::unique_ptr<TH2D> h2d_p{new TH2D("h2d_p",
                                               "h2d_p",
                                               100,
                                               0.25 * min,
                                               1.25 * max,
                                               100,
                                               0.25 * min,
                                               1.25 * max)};
        h2d_p.get()->SetStats(0);
        std::ostringstream ss;
        ss.str("");ss.clear();
        ss << "stdAbsCon=" << calculateStdAbsCon(points_p) << ";A_{m}', %;A_{c}, %";
        h2d_p.get()->SetTitle(ss.str().c_str());

        const std::string psName_p{"output_p.ps"};
        std::unique_ptr<TCanvas> c_p{new TCanvas("c_p", "c_p", 1024, 960)};
        gPad->SetGrid();
        c_p.get()->Print((psName_p + '[').c_str());
        h2d_p.get()->Draw();
        std::unique_ptr<TLine> dLine{new TLine(0.25 * min, 0.25 * min, 1.25 * max, 1.25 * max)};
        dLine.get()->Draw("SAME");
        gr_p.get()->Draw("P");




        std::map<std::pair<std::string, Color_t>, Points> subPoints{
            { std::make_pair(R"(\bsample([1-9]|[12][0-9]|30)\b)", kOrange), Points() },
        };

        for (size_t i{0}; i < points_p.size(); ++i) {
            for (auto &item : subPoints) {
                std::regex pattern(item.first.first);
                if (std::regex_search(points_p.at(i).sample, pattern)) {
                    TMarker m{points_p.at(i).x, points_p.at(i).y, 21};
                    m.SetMarkerSize(1.5);
                    m.SetMarkerColor(item.first.second);
                    m.DrawClone("SAME");
                }
            }
        }

        // std::map<std::pair<std::string, Color_t>, Points> subPoints{
        // { std::make_pair("berez_7_w", kOrange), Points() },
        // { std::make_pair("berez_11", kRed), Points() },
        // { std::make_pair("berez_2", kBlue), Points() },
        // { std::make_pair("berez_7_w0_", kGreen), Points() },
        // { std::make_pair("berez_7_w5_", kGreen), Points() },
        // { std::make_pair("berez_7_w10_", kGreen), Points() },
        // { std::make_pair("berez_7_w15_", kGreen), Points() },
        // };


        // for (size_t i{0}; i < points_p.size(); ++i) {
        //     for (auto &item : subPoints) {
        //         if (points_p.at(i).sample.find(item.first.first) != std::string::npos) {
        //             TMarker m{points_p.at(i).x, points_p.at(i).y, 21};
        //             m.SetMarkerSize(1.5);
        //             m.SetMarkerColor(item.first.second);
        //             m.DrawClone("SAME");
        //         }
        //     }
        // }

        // std::unique_ptr<TPaveText> pt{new TPaveText(0.1, 0.65, 0.6, 0.9, "NDC")};
        // pt.get()->SetFillColor(0);
        // pt.get()->SetBorderSize(1);

        // auto doubleToString = [](double value, int precision = 1) {
        //     std::ostringstream oss;
        //     oss << std::fixed << std::setprecision(precision) << value;
        //     return oss.str();
        // };

        // // Add entries from map keys
        // for (const auto& entry : subPoints) {
        //     std::string key = entry.first.first;
        //     Color_t color = entry.first.second;

        //     for (size_t i{0}; i < entry.second.l.size(); i++) {
        //         key.append("(");
        //         key.append(doubleToString(entry.second.x.at(i)));
        //         key.append(",");
        //         key.append(doubleToString(entry.second.y.at(i)));
        //         key.append(")");
        //     }

        //     TText *text = pt.get()->AddText(key.c_str());
        //     text->SetTextColor(color);
        // }
        // pt.get()->Draw("SAME");

        c_p.get()->Print(psName_p.c_str());
        c_p.get()->Print((psName_p + ']').c_str());
        c_p.get()->Close();

        // std::regex m_a_c{R"(pulp_rot_berez_7_w\d+_\d+|pulp_rot_berez_11_1500g_w\d+_\d+|pulp_rot_berez_7_1500g_w\d+_\d+|N12_\d+_\d+)"};//3
        // auto data_c{getData(fileName, columnElement, chem, m_a_c)};
        // std::vector<Point> points_a_c{getPointsByType(data, ChemResult::Type::A)};


               // check 3
        //          std::regex m_a{R"(pulp_rot_berez_7_w\d+_sum|pulp_rot_berez_11_1500g_w\d+_sum|pulp_rot_berez_7_1500g_w\d+_sum|N12)"};//sum
                  // std::regex m_a{R"(pulp_rot_berez_7_w\d+_\d+|pulp_rot_berez_11_1500g_w\d+_\d+|pulp_rot_berez_7_1500g_w\d+_\d+|N12_\d+_\d+)"};//3
        //std::regex m_a{R"(pulp_rot_berez_7_w\d+_sum|pulp_rot_berez_2_w\d+_sum|pulp_rot_berez_11_w\d+_sum|d+_sum|barz_blind|berez_blind|N12|pulp_rot_berez_6_w\d+_sum)"};//!
               // auto data2{getFitResults(fileName, columnElement, chem, m_a)};
               // drawCorrGraphWr3_a(data2, f);
        //        // check 3
        ////        std::regex m_w{R"(pulp_rot_berez_7_w\d+_sum|pulp_rot_berez_11_1500g_w\d+_sum|pulp_rot_berez_7_1500g_w\d+_sum)"};//sum
        //        std::regex m_w{R"(pulp_rot_berez_7_w\d+_sum|pulp_rot_berez_11_1500g_w\d+_\d+|pulp_rot_berez_7_1500g_w\d+_\d+)"};//3
        ////        std::regex m_w{R"(pulp_rot_berez_7_w\d+_sum|pulp_rot_berez_2_w\d+_sum|pulp_rot_berez_11_w\d+_sum|pulp_rot_berez_6_w\d+_sum)"};//!
        //        auto data3{getFitResults(fileName, columnElement, chem, m_w)};
        //        drawCorrGraphWr3_w(data3, f);

        //        return 0;

        //++++
//        std::vector<TLatex> labels;

//        for (size_t i{0}; i < points.size(); ++i) {
//            auto pos{points.l.at(i).find_last_of("_")};
//            auto text{points.l.at(i)};
//            if (pos != std::string::npos && pos == points.l.at(i).length() - 1)
//            {
//                text = text.substr(0, pos);
//            }
//            TLatex l(points.x.at(i), points.y.at(i) + 1.25 * points.yErr.at(i), text.c_str());
//            l.SetTextAngle(90);
//            l.SetTextAlign(12);
//            l.SetTextSize(0.02);
//            labels.push_back(l);
//        }



//            std::vector<Point> points_w;

//            for (const auto &[key, value] : data) {
//                for (const auto &fr : value.fitResults) {
//                    const auto a{value.chemResult.a};
//                    const auto w{value.chemResult.w};
//                    if (a.has_value()) {
//                        Point point;
//                        point.sample = key;
//                        point.chemResult = value.chemResult;
//                        point.x = points_a.size();
//                        point.xErr = 0.1;
//                        point.y = a.value();
//                        point.yErr = 0.1 * a.value();
//                        point.fitResult = fr;
//                        points_a.push_back(point);
//                    }
//                    if (w.has_value()) {
//                        Point point;
//                        point.sample = key;
//                        point.chemResult = value.chemResult;
//                        point.x = points_w.size();
//                        point.xErr = 0.1;
//                        point.y = w.value();
//                        point.yErr = 0.03 * w.value();
//                        point.fitResult = fr;
//                        points_w.push_back(point);
//                    }
//                }
//            }

//            for (auto &p : points_w) {
//                p.x  = p.x + points_a.size();
//            }
//            std::vector<Point> points;
//            points.insert(points.end(), points_a.cbegin(), points_a.cend());
//            points.insert(points.end(), points_w.cbegin(), points_w.cend());
//            return points;
//        };

//        auto points{getPoints(data)};

//        for (const auto &p : points) {
//            std::cout << p.sample << " " << p.x << " " << p.y;
//            p.fitResult.print();
//            std::cout << std::endl;
//        }

//        FitFunction fObj(mmn, dA);
//        std::unique_ptr<TF1> f{new TF1("f", fObj, points.x.front(), points.x.back(), 7)};

//        auto data1{getFitResults(fileName, columnElement, chem, m)};

//        Points points;

//        addPointsByValue(data1, points, Data1::Value::A);
//        addPointsByValue(data1, points, Data1::Value::W);

//        std::map<int, std::vector<double>> mmn;

//        addMmnByValue(data1, mmn, Data1::Value::A);
//        addMmnByValue(data1, mmn, Data1::Value::W);

//        std::cout << data1.size() << std::endl;

//        std::cout << points.x.size() << " " << mmn.size() << std::endl;

//        std::unique_ptr<TGraphErrors> gr{new TGraphErrors(static_cast<int>(points.x.size()), &points.x[0], &points.y[0], &points.xErr[0], &points.yErr[0])};
//        gr.get()->SetMarkerSize(1.5);
//        gr.get()->SetMarkerStyle(21);
//        gr.get()->SetTitle(";N_{probe};[...A, ...W]");

//        std::vector<TLatex> labels;

//        for (size_t i{0}; i < points.x.size(); ++i)
//        {
//            auto pos{points.l.at(i).find_last_of("_")};
//            auto text{points.l.at(i)};
//            if (pos != std::string::npos && pos == points.l.at(i).length() - 1)
//            {
//                text = text.substr(0, pos);
//            }
//            TLatex l(points.x.at(i), points.y.at(i) + 1.25 * points.yErr.at(i), text.c_str());
//            l.SetTextAngle(90);
//            l.SetTextAlign(12);
//            l.SetTextSize(0.02);
//            labels.push_back(l);
//        }

//        std::map<int, double> dA;

//        int xx{0};
//        for (auto it{data1.begin()}; it != data1.end(); ++it) {
//            for (size_t i{0}; i < it->second.fr.size(); ++i) {
//                std::optional<double> v;
//                v = (*it).second.chem.a;
//                if (v.has_value()) {
//                    dA[xx] = v.value();
//                    xx++;
//                }
//            }
//        }

//        FitFunction_3 fObj(mmn, dA);
//        std::unique_ptr<TF1> f{new TF1("f", fObj, points.x.front(), points.x.back(), 7)};

//// par6 free new
////        1  p0           1.55468e+00   2.91982e-02   5.12994e-06   4.88038e-02
////        2  p1           7.55271e-01   1.62498e-02   4.57209e-06  -5.43493e-02
////        3  p2          -1.22767e+01   6.36519e-01   1.07210e-04   2.25340e-03
////        4  p3           1.00997e+02   9.86184e+00   9.84859e-04  -1.88231e-04
////        5  p4           6.33645e-01   1.05930e-01   8.90741e-05  -2.98749e-03
////        6  p5           1.06754e+00   1.85230e-01   1.44951e-05   1.24332e-02
////        7  p6           2.73195e+00   3.61006e+00   7.23964e-04   1.59404e-04


//        const std::vector<double> parameters = {

//            1.55468e+00,
//            7.55271e-01,
//            -1.22767e+01,
//            1.00997e+02,
//            6.33645e-01,
//            1.06754e+00,
//            2.73195e+00

////            1.55468e+00,
////            7.55269e-01,
////           -1.22767e+01,
////            1.05060e+02,
////            6.69100e-01,
////            1.17662e+00,
////            0.00000e+00
//        };

//        auto setInitialParameters = [&parameters](TF1 *f){
//            for (auto it{parameters.begin()}; it != parameters.end(); it++) {
//                f->SetParameter(std::distance(parameters.begin(), it), *it);
//            }
//        };

//        auto setParametersLimits = [&parameters](TF1 *f, double min = 0.8, double max = 1.2){
//            for (auto it{parameters.begin()}; it != parameters.end(); it++) {
//                auto minPar{min * (*it)};
//                auto maxPar{max * (*it)};
//                if (*it < 0) {
//                    maxPar = min * (*it);
//                    minPar = max * (*it);
//                }
//                f->SetParLimits(std::distance(parameters.begin(), it), minPar, maxPar);
//            }
//        };

//        setInitialParameters(f.get());
////        setParametersLimits(f.get());

//        f.get()->SetNpx(10 * static_cast<int>(points.x.size()));

//        gr.get()->Fit(f.get(), "R");

//        const std::string psName{"output.ps"};
//        std::unique_ptr<TCanvas> c{new TCanvas("c", "c", 1024, 960)};
//        c.get()->Print((psName + '[').c_str());
//        gr.get()->Draw("APL");



//        auto useSub{true};
//        if (useSub)
//        {
//            std::map<std::pair<std::string, Color_t>, Points> subPoints{
//                { std::make_pair("tochka_", kRed), Points() },
//                { std::make_pair("sector_", kBlue), Points() },
//                { std::make_pair("field_", kGreen), Points() },
//                { std::make_pair("other", kBlack), Points() },
//            };


//            for (size_t i{0}; i < points.x.size(); ++i)
//            {
//                auto isOther{false};
//                for (auto &item : subPoints)
//                {
//                    if (points.l.at(i).find(item.first.first) != std::string::npos)
//                    {
//                        TMarker m{points.x.at(i), points.y.at(i), 21};
//                        m.SetMarkerSize(1.5);
//                        m.SetMarkerColor(item.first.second);
//                        m.DrawClone("SAME");
//                        item.second.l.push_back(points.l.at(i));
//                        item.second.x.push_back(points.x.at(i));
//                        item.second.y.push_back(points.y.at(i));
//                        item.second.xErr.push_back(0.1);
//                        item.second.yErr.push_back(0.5);

//                        isOther = true;
//                    }
//                }
//                if (!isOther)
//                {
//                    subPoints.at({"other", kBlack}).l.push_back(points.l.at(i));
//                    subPoints.at({"other", kBlack}).x.push_back(points.x.at(i));
//                    subPoints.at({"other", kBlack}).y.push_back(points.y.at(i));
//                    subPoints.at({"other", kBlack}).xErr.push_back(0.1);
//                    subPoints.at({"other", kBlack}).yErr.push_back(0.5);
//                }
//            }
//        }

//        for (const auto &item : labels)
//        {
//            item.DrawClone("SAME");
//        }

//        c.get()->Print(psName.c_str());
//        c.get()->Print((psName + ']').c_str());
//        c.get()->Close();

//        // check 3
////          std::regex m_a{R"(pulp_rot_berez_7_w\d+_sum|pulp_rot_berez_11_1500g_w\d+_sum|pulp_rot_berez_7_1500g_w\d+_sum|N12)"};//sum
//           std::regex m_a{R"(pulp_rot_berez_7_w\d+_\d+|pulp_rot_berez_11_1500g_w\d+_\d+|pulp_rot_berez_7_1500g_w\d+_\d+|N12_\d+_\d+)"};//3
////std::regex m_a{R"(pulp_rot_berez_7_w\d+_sum|pulp_rot_berez_2_w\d+_sum|pulp_rot_berez_11_w\d+_sum|d+_sum|barz_blind|berez_blind|N12|pulp_rot_berez_6_w\d+_sum)"};//!
//        auto data2{getFitResults(fileName, columnElement, chem, m_a)};
//        drawCorrGraphWr3_a(data2, f);
//        // check 3
////        std::regex m_w{R"(pulp_rot_berez_7_w\d+_sum|pulp_rot_berez_11_1500g_w\d+_sum|pulp_rot_berez_7_1500g_w\d+_sum)"};//sum
//        std::regex m_w{R"(pulp_rot_berez_7_w\d+_sum|pulp_rot_berez_11_1500g_w\d+_\d+|pulp_rot_berez_7_1500g_w\d+_\d+)"};//3
////        std::regex m_w{R"(pulp_rot_berez_7_w\d+_sum|pulp_rot_berez_2_w\d+_sum|pulp_rot_berez_11_w\d+_sum|pulp_rot_berez_6_w\d+_sum)"};//!
//        auto data3{getFitResults(fileName, columnElement, chem, m_w)};
//        drawCorrGraphWr3_w(data3, f);

//        return 0;


//        const auto fileName_1{"rea.elts.stroy.txt"};
////        std::regex s{"sum"};
//        // std::regex s{"\\d+_s"};
////        std::regex s{"check_bereza"};
////        std::regex s{"(pulp_rot_N12_\\d+_\\d+)"};
////        std::regex s{"(pulp_rot_N12_\\d+_sum)"};
////        std::regex s{"(pulp_rot_N12_[13579]_sum|pulp_rot_N12_[1-9]\\d*[13579]_sum)"};
////        std::regex s{"(pulp_rot_N12_\\d+_sum|pulp_rot_berez_\\d+_sum|_blind_a\\d+p\\d+_sum)"};
////        std::regex s{"(pulp_rot_berez_111_w5_|pulp_rot_berez_111_w10_|pulp_rot_berez_111_w15_)"};
////        std::regex s{"(pulp_rot_berez_6_w\\d+_\\d+)"};
////        std::regex s{"(pulp_rot_N12_\\d*[02468]_\\d+)"};

////        std::regex s{R"(pulp_rot_berez_11_w\d+_sum|pulp_rot_berez_11_w\d+p\d+_sum|pulp_rot_N12_\d+_sum|_blind_a\d+p\d+_sum)"};
////        std::regex s{R"(pulp_rot_berez_11_w\d+_sum|pulp_rot_berez_11_w\d+p\d+_sum|pulp_rot_N12_\d+_sum|_blind_a\d+p\d+_sum)"};
////        std::regex s{R"(pulp_rot_berez_11_w\d+_sum|pulp_rot_berez_11_w\d+p\d+_sum|pulp_rot_berez_7_w\d+_sum|pulp_rot_berez_2_sum)"};//!
////        auto data1Sum{getFitResults(fileName_1, columnElement, chem, s)};
////        calcConv(data1Sum, f, value);

////        std::regex t{"_\\d+_"};
////        auto dataBlindSum{getFitResults(fileName_1, columnElement, chemBlind, t)};
////        dataBlindSum.insert(data1Sum.begin(), data1Sum.end());
////        calcConv(dataBlindSum, f, value);

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
    if (ss.fail()) {
        throw my_error("Can\'t convert: " + str);
    }
    return d;
}

std::vector<std::string>splitLineToStrs(const std::string &line)
{
    std::stringstream ss(line);
    std::string str;
    std::vector<std::string> strs;
    while (ss >> str) {
        strs.push_back(str);
    }
    return strs;
}

std::map<std::string, Data> getData(const std::string &fileName,
                          const std::map<int, std::string> &columnElement,
                          const std::map<std::string, ChemResult> &chem,
                          const std::regex &pattern)
{
    std::ifstream ifs(fileName);
    if (!ifs.is_open()) {
        throw my_error("Can't open file \"" + fileName + "\"");
    }
    std::string line;

    std::map<std::string, Data> d;

    while (getline(ifs, line)) {
        auto strs{splitLineToStrs(line)};
        try {
            auto it = std::find_if(chem.begin(), chem.end(), [&strs, &pattern] (std::pair<std::string, ChemResult> chemItem){
                return strs.at(0).find(chemItem.first) != std::string::npos && std::regex_search(strs.at(0), pattern);
            });

            if (it != chem.end()) {
                std::cout << strs.at(0) << std::endl;
                FitResult fR;
                for (const auto &[key, value] : columnElement)
                {
                    fR.elementResults.push_back({value,
                                            strToDouble(strs.at(static_cast<unsigned int>(key))),
                                            strToDouble(strs.at(static_cast<unsigned int>(key + 1)))
                                           });
                }
                const auto &[key, value] = *it;
                d[key].chemResult.a = value.a;
                d[key].chemResult.w = value.w;
                d[key].fitResults.push_back(fR);
            }
        }  catch (...) {
            std::cout << "Error adding e to data";
        }
    }
    ifs.close();
    return d;
}




//std::map<std::string, Data1> getFitResults(const std::string &fileName,
//                   const std::map<int, std::string> &columnElement,
//                   const std::map<std::string, ChemResult> &chem,
//                   const std::regex &pattern)
//{
//    std::ifstream ifs(fileName);
//    if (!ifs.is_open())
//    {
//        throw my_error("Can't open file \"" + fileName + "\"");
//    }
//    std::string line;

//    std::map<std::string, Data1> data;

//    while (getline(ifs, line))
//    {
//        auto strs{splitLineToStrs(line)};
//        try
//        {
//            auto it = std::find_if(chem.begin(), chem.end(), [&strs, &pattern] (std::pair<std::string, ChemResult> chemItem){
//                return strs.at(0).find(chemItem.first) != std::string::npos && std::regex_search(strs.at(0), pattern);
//            });

//            if (it != chem.end())
//            {
//                std::cout << strs.at(0) << std::endl;
//                std::vector<FitResult> fR;
//                for (const auto &item : columnElement)
//                {
//                    fR.push_back({ item.second, strToDouble(strs.at(static_cast<unsigned int>(item.first))),
//                                   strToDouble(strs.at(static_cast<unsigned int>(item.first + 1))) });
//                }

////                for (auto &i : fR) {
////                    auto itPreGrad = preGrad.find(i.e);
////                    if (itPreGrad != preGrad.end()) {
////                        i.value = itPreGrad->second.first +  itPreGrad->second.second * i.value;
////                    }
////                }


//                data[(*it).first].chem.a = it->second.a;
//                data[(*it).first].chem.w = it->second.w;
//                data[(*it).first].fr.push_back(fR);
//            }
//        }  catch (...) {
//            std::cout << "Error adding e to data";
//        }
//    }
//    ifs.close();
//    return data;
//}

//void addPointsByValue(const std::map<std::string, Data1> &data,
//                        Points &points,
//                        const Data1::Value value)
//{
//    int xx{static_cast<int>(points.x.size())};
//    for (auto it{data.begin()}; it != data.end(); ++it)
//    {
//        for (size_t i{0}; i < it->second.fr.size(); ++i)
//        {
//            std::string label{it->first + "_" + i};
//            std::optional<double> v;
//            switch (value) {
//            case Data1::Value::A:
//                v = (*it).second.chem.a;
//                break;
//            case Data1::Value::W:
//                v = (*it).second.chem.w;
//                break;
//            }
//            if (v.has_value())
//            {
//                points.l.push_back(label);
//                points.x.push_back(xx++);
//                points.xErr.push_back(0.01);
//                points.y.push_back(v.value());
////                points.yErr.push_back(0.25 * v.value());
//                auto yE{0.0};
//                if (value == Data1::Value::A) {
//                    yE = 0.1 * v.value();
//                }
//                if (value == Data1::Value::W) {
//                    yE = 0.03 * v.value();
//                }
//                points.yErr.push_back(yE);
//            }
//        }

//    }
//}



//void addMmnByValue(const std::map<std::string, Data1> &data,
//                        std::map<int, std::vector<double>> &mmn,
//                        const Data1::Value value)
//{


//    int xx{static_cast<int>(mmn.size())};
//    for (auto it{data.begin()}; it != data.end(); ++it)
//    {
//        for (size_t i{0}; i < it->second.fr.size(); ++i)
//        {
//            std::string label{it->first + "_" + i};
//            std::optional<double> v;
//            switch (value) {
//            case Data1::Value::A:
//                v = (*it).second.chem.a;
//                break;
//            case Data1::Value::W:
//                v = (*it).second.chem.w;
//                break;
//            }
//            if (v.has_value())
//            {
//                mmn[xx].push_back(it->second.fr.at(i).at(static_cast<size_t>(0)).value); // Al
//                mmn[xx].push_back(it->second.fr.at(i).at(static_cast<size_t>(1)).value); // C
//                mmn[xx].push_back(it->second.fr.at(i).at(static_cast<size_t>(2)).value); // N
//                mmn[xx].push_back(it->second.fr.at(i).at(static_cast<size_t>(3)).value); // O
//                mmn[xx].push_back(it->second.fr.at(i).at(static_cast<size_t>(4)).value); // Si
//                xx++;
//            }
//        }
//    }
//}

//void useSub(const Points &points,
//            std::string &str,
//            const Data1::Value value,
//            const bool isLabels = false)
//{
//    std::map<std::pair<std::string, Color_t>, Points> subPoints{
//        { std::make_pair("N12", kRed), Points() },
//        { std::make_pair("berez", kBlue), Points() },
//        { std::make_pair("p43", kGreen), Points() },
//        { std::make_pair("blind", kOrange), Points() },
//        { std::make_pair("raspad", kCyan), Points() },
//        { std::make_pair("other", kBlack), Points() },
//    };

//    for (size_t i{0}; i < points.x.size(); ++i)
//    {
//        auto isOther{true};
//        for (auto &item : subPoints)
//        {
//            if (points.l.at(i).find(item.first.first) != std::string::npos)
//            {
//                TMarker m{points.x.at(i), points.y.at(i), 21};
//                m.SetMarkerSize(1.5);
//                m.SetMarkerColor(item.first.second);
//                m.DrawClone("SAME");
//                item.second.l.push_back(points.l.at(i));
//                item.second.x.push_back(points.x.at(i));
//                item.second.y.push_back(points.y.at(i));
//                item.second.xErr.push_back(0.1);
//                item.second.yErr.push_back(0.5);

////                isOther = false;
//            }
//        }
//        if (isOther)
//        {
//            TMarker m{points.x.at(i), points.y.at(i), 21};
//            m.SetMarkerSize(1.5);
//            m.SetMarkerColor(kBlack);
////            m.DrawClone("SAME");
//            subPoints.at({"other", kBlack}).l.push_back(points.l.at(i));
//            subPoints.at({"other", kBlack}).x.push_back(points.x.at(i));
//            subPoints.at({"other", kBlack}).y.push_back(points.y.at(i));
//            subPoints.at({"other", kBlack}).xErr.push_back(0.1);
//            subPoints.at({"other", kBlack}).yErr.push_back(0.5);
//        }
//    }
//    std::stringstream ss;
//    ss.str("");ss.clear();
//    ss << (value == Data1::Value::A ? "Ad" : "Wr");
//    ss << ": stdAbs=" << std::setprecision(3);
//    auto stdAbs1 = [](const Points &points){
//        std::vector<double> d2;
//        for (size_t i{0}; i < points.x.size(); ++i)
//        {
//            d2.push_back(std::pow(points.y.at(i) - points.x.at(i), 2));
//        }
//        return std::sqrt(std::accumulate(d2.begin(), d2.end(), 0.0) / static_cast<double>(d2.size()));
//    };
//    for (auto item : subPoints)
//    {
//        if (!std::isnan(stdAbs1(item.second)))
//        {
//            ss << "[#color[" << static_cast<int>(item.first.second) << "]{" << stdAbs1(item.second) << "}] ";
//        }
//    }

//    ss << ";AGP-K, %;Chem, %";
//    str = ss.str();

//    if (isLabels)
//    {
//        std::vector<TLatex> labels;
//        for (size_t i{0}; i < points.x.size(); ++i)
//        {
//            TLatex l(points.x.at(i), points.y.at(i) + 1.25 * points.yErr.at(i), points.l.at(i).c_str());
//            l.SetTextAngle(90);
//            l.SetTextAlign(12);
//            l.SetTextSize(0.02);
//            labels.push_back(l);
//        }
//        for (const auto &item : labels)
//        {
//            item.DrawClone("SAME");
//        }
//    }
//}

//void addSubPoints(const Points &points, std::map<std::pair<std::string, Color_t>, Points> &subPoints)
//{
//    for (size_t i{0}; i < points.x.size(); ++i)
//    {
//        auto isOther{true};
//        for (auto &item : subPoints)
//        {
//            if (points.l.at(i).find(item.first.first) != std::string::npos)
//            {
////                TMarker m{points.x.at(i), points.y.at(i), 21};
////                m.SetMarkerSize(1.5);
////                m.SetMarkerColor(item.first.second);
////                m.DrawClone("SAME");
//                item.second.l.push_back(points.l.at(i));
//                item.second.x.push_back(points.x.at(i));
//                item.second.y.push_back(points.y.at(i));
//                item.second.xErr.push_back(0.1);
//                item.second.yErr.push_back(0.5);

//                isOther = false;
//            }
//        }
//        if (isOther)
//        {
////            TMarker m{points.x.at(i), points.y.at(i), 21};
////            m.SetMarkerSize(1.5);
////            m.SetMarkerColor(kBlack);
////            m.DrawClone("SAME");
//            subPoints.at({"other", kBlack}).l.push_back(points.l.at(i));
//            subPoints.at({"other", kBlack}).x.push_back(points.x.at(i));
//            subPoints.at({"other", kBlack}).y.push_back(points.y.at(i));
//            subPoints.at({"other", kBlack}).xErr.push_back(0.1);
//            subPoints.at({"other", kBlack}).yErr.push_back(0.5);
//        }
//    }
//}

//void useSub1(const Points &points,
//            std::string &str,
//            const Data1::Value value,
//            const bool isLabels = false)
//{


//    std::map<std::pair<std::string, Color_t>, Points> subPoints{
//        { std::make_pair("N12", kRed), Points() },
//        { std::make_pair("berez", kBlue), Points() },
//        { std::make_pair("p43", kGreen), Points() },
//        { std::make_pair("blind", kOrange), Points() },
//        { std::make_pair("raspad", kCyan), Points() },
//        { std::make_pair("other", kBlack), Points() },
//    };

//    addSubPoints(points, subPoints);
//    std::stringstream ss;
//    ss.str("");ss.clear();
//    ss << (value == Data1::Value::A ? "Ad" : "Wr");
//    ss << ": stdAbs=" << std::setprecision(3);
//    auto stdAbs1 = [](const Points &points){
//        std::vector<double> d2;
//        for (size_t i{0}; i < points.x.size(); ++i)
//        {
//            d2.push_back(std::pow(points.y.at(i) - points.x.at(i), 2));
//        }
//        return std::sqrt(std::accumulate(d2.begin(), d2.end(), 0.0) / static_cast<double>(d2.size()));
//    };
//    for (auto item : subPoints)
//    {
////        for (size_t i{0}; i < item.second.x.size(); ++i)
////        {
////            TMarker m{item.second.x.at(i), item.second.y.at(i), 21};
////            m.SetMarkerSize(1.5);
////            m.SetMarkerColor(item.first.second);
////            m.DrawClone("SAME");
////        }
//        if (!std::isnan(stdAbs1(item.second)))
//        {
//            ss << "[#color[" << static_cast<int>(item.first.second) << "]{" << stdAbs1(item.second) << "}] ";
//        }
//    }

//    ss << ";Number;Chem/AGP-K, %";
//    str = ss.str();

//    std::map<std::pair<std::string, Color_t>, Points> subPoints1{
//        { std::make_pair("N12", kRed), Points() },
//        { std::make_pair("berez", kBlue), Points() },
//        { std::make_pair("p43", kGreen), Points() },
//        { std::make_pair("111", kOrange), Points() },
//        { std::make_pair("raspad", kCyan), Points() },
//        { std::make_pair("other", kBlack), Points() },
//    };



//    Points pointsC{points};
//    Points pointsR{points};
//    for (size_t i{0}; i < points.x.size(); ++i)
//    {
//        pointsC.x.at(i) = static_cast<double>(i + 1);
//        pointsC.yErr.at(i) = 0.0;
//        pointsR.x.at(i) = static_cast<double>(i + 1);
//        pointsR.y.at(i) = points.x.at(i);
//        pointsR.yErr.at(i) = points.xErr.at(i);
//    }

//    addSubPoints(pointsR, subPoints1);

//    writePointsToFile("outputR.txt", pointsR);
//    writePointsToFile("outputC.txt", pointsC);

//    std::unique_ptr<TGraphErrors> grC{new TGraphErrors(static_cast<int>(pointsC.x.size()),
//                                                       &pointsC.x[0],
//                                                       &pointsC.y[0],
//                                                       &pointsC.xErr[0],
//                                                       &pointsC.yErr[0])};
//    grC.get()->SetLineColor(kRed);
//    grC.get()->SetLineWidth(3);
//    grC.get()->SetMarkerSize(1.5);
//    grC.get()->SetMarkerStyle(21);

//    std::unique_ptr<TGraphErrors> grR{new TGraphErrors(static_cast<int>(pointsR.x.size()),
//                                                       &pointsR.x[0],
//                                                       &pointsR.y[0],
//                                                       &pointsR.xErr[0],
//                                                       &pointsR.yErr[0])};
//    grR.get()->SetMarkerSize(1.5);
//    grR.get()->SetMarkerStyle(21);


//    grC.get()->DrawClone("SAME L");
//    grR.get()->DrawClone("SAME P");

//    for (auto item : subPoints1)
//    {
//        for (size_t i{0}; i < item.second.x.size(); ++i)
//        {
//            TMarker m{item.second.x.at(i), item.second.y.at(i), 21};
//            m.SetMarkerSize(1.5);
//            m.SetMarkerColor(item.first.second);
//            m.DrawClone("SAME");
//        }
//    }

//    if (isLabels)
//    {
//        std::vector<TLatex> labels;
//        for (size_t i{0}; i < pointsR.x.size(); ++i)
//        {
//            TLatex l(pointsR.x.at(i), pointsR.y.at(i) + 1.25 * pointsR.yErr.at(i), pointsR.l.at(i).c_str());
//            l.SetTextAngle(90);
//            l.SetTextAlign(12);
//            l.SetTextSize(0.02);
//            labels.push_back(l);
//        }
//        for (const auto &item : labels)
//        {
//            item.DrawClone("SAME");
//        }
//    }
//}

//void writePointsToFile(const std::string fileName, const Points &points)
//{
//    std::ofstream ofs;
//    ofs.open(fileName, std::ios::out);
//    for (size_t i{0}; i < points.l.size(); ++i)
//    {
//        ofs << points.x.at(i) << " " << points.y.at(i)  << " " << points.l.at(i) << std::endl;
//    }
//    ofs.close();
//}

//void calcConv(const std::map<std::string, Data1> &data,
//              const std::unique_ptr<TF1> &f,
//              const Data1::Value value)
//{
//    auto p0{f.get()->GetParameter(0)};
//    auto p1{f.get()->GetParameter(1)};
//    auto p2{f.get()->GetParameter(2)};
//    auto p3{f.get()->GetParameter(3)};
//    auto p4{f.get()->GetParameter(4)};
//    auto p5{f.get()->GetParameter(5)};
//    auto p6{f.get()->GetParameter(6)};
//    auto p7{f.get()->GetParameter(7)};

//    Points points;
//    for (auto it{data.begin()}; it != data.end(); ++it)
//    {
//        for (size_t i{0}; i < it->second.fr.size(); ++i)
//        {
//            std::optional<double> v;
//            switch (value) {
//            case Data1::Value::A:
//                v = (*it).second.chem.a;
//                break;
//            case Data1::Value::W:
//                v = (*it).second.chem.w;
//                break;
//            }
//            if (v.has_value())
//            {
//                auto Al{it->second.fr.at(i).at(static_cast<size_t>(0)).value};
//                auto C{it->second.fr.at(i).at(static_cast<size_t>(1)).value};
//                auto N{it->second.fr.at(i).at(static_cast<size_t>(2)).value};
//                auto O{it->second.fr.at(i).at(static_cast<size_t>(3)).value};
//                auto Si{it->second.fr.at(i).at(static_cast<size_t>(4)).value};

//                auto res{0.0};
//                auto w{
//                    p3 * O
//                    - p4 * Si
//                    - p5 * Al
//                    - p6 * N
//                    + p7
//                };
//                auto a{
//                    (p1 + p0 * C)
//                    / (100.0 - p2 * w) * 100.0
//                };

//                switch (value) {
//                case Data1::Value::A:
//                    res = a;
//                    break;
//                case Data1::Value::W:
//                    res = w;
//                    break;
//                }

//                points.l.push_back(it->first);
//                points.x.push_back(res);
//                points.y.push_back(v.value());
//                points.xErr.push_back(0.1);
//                points.yErr.push_back(0.5);
//            }
//        }
//    }

//    std::vector<double> d2;
//    for (size_t i{0}; i < points.x.size(); ++i)
//    {
//        d2.push_back(std::pow(points.y.at(i) - points.x.at(i), 2));
//    }



//    auto stdAbs{std::sqrt(std::accumulate(d2.begin(), d2.end(), 0.0) / d2.size())};
//    auto avg{std::accumulate(points.x.begin(), points.x.end(), 0.0) / points.x.size()};
//    std::cout << "convergence: " << "avg = " << avg << " stdAbs = " << stdAbs << std::endl;

//    std::unique_ptr<TGraphErrors> gr{new TGraphErrors(static_cast<int>(points.x.size()), &points.x[0], &points.y[0], &points.xErr[0], &points.yErr[0])};
//    gr.get()->SetMarkerSize(1.5);
//    gr.get()->SetMarkerStyle(21);

//    std::ofstream ofs("output_data.csv", std::ios_base::app);

//    std::vector<TLatex> labels;
//    for (size_t i{0}; i < points.x.size(); ++i)
//    {

//        auto pos{points.l.at(i).find_last_of("_")};
//        auto text{points.l.at(i)};
//        if (pos != std::string::npos && pos == points.l.at(i).length() - 1)
//        {
//            text = text.substr(0, pos);
//        }
//        std::cout << i << " " << text << std::endl;
//        TLatex l(points.x.at(i), points.y.at(i) + 1.25 * points.yErr.at(i), text.c_str());
//        l.SetTextAngle(90);
//        l.SetTextAlign(12);
//        l.SetTextSize(0.02);
//        labels.push_back(l);

//        if (ofs.is_open()) {
//            ofs << text << " " << points.x.at(i) << " " << points.xErr.at(i) << " "
//                << points.y.at(i) << " " << points.yErr.at(i) << std::endl;
//        }
//    }

//    ofs.close();


//    std::stringstream ss;
//    ss.str("");ss.clear();
//    ss << (value == Data1::Value::A ? "Ad" : "Wr")
//       << ": avg = "
//       << avg << "%"
//       << ", stdAbs = "
//       << stdAbs << "%"
//       << ";AGP-K, %;Chem, %";

//    std::unique_ptr<TH2D> h2dConv{new TH2D("h2dConv",
//                                           ss.str().c_str(),
//                                           static_cast<int>(points.y.size()),
//                                           0.75 * (*std::min_element(points.y.begin(), points.y.end())),
//                                           1.25 * (*std::max_element(points.y.begin(), points.y.end())),
//                                           static_cast<int>(points.y.size()),
//                                           0.75 * (*std::min_element(points.y.begin(), points.y.end())),
//                                           1.25 * (*std::max_element(points.y.begin(), points.y.end())))};
//    h2dConv->SetStats(0);
//    std::unique_ptr<TLine> lConv{new TLine(0.75 * (*std::min_element(points.y.begin(), points.y.end())),
//                                           0.75 * (*std::min_element(points.y.begin(), points.y.end())),
//                                           1.25 * (*std::max_element(points.y.begin(), points.y.end())),
//                                           1.25 * (*std::max_element(points.y.begin(), points.y.end())))};
//    const std::string psName{"output_conv.ps"};
//    std::unique_ptr<TCanvas> c{new TCanvas("c", "c", 1024, 960)};
//    c.get()->Print((psName + '[').c_str());
//    h2dConv.get()->Draw();
////    gr.get()->Draw("P");

//    lConv.get()->Draw("SAME");
//    std::vector<std::string> uniqueL{points.l};
//    std::sort(uniqueL.begin(), uniqueL.end());
//    auto uniqueIt{std::unique(uniqueL.begin(), uniqueL.end())};
//    uniqueL.erase(uniqueIt, uniqueL.end());
//    Points avgPoints;
//    avgPoints.l = uniqueL;



//    for (size_t i{0}; i < avgPoints.l.size(); ++i)
//    {
//        std::vector<double> tmpX;
//        std::vector<double> tmpXerr;
//        std::vector<double> tmpY;
//        std::vector<double> tmpYerr;
//        auto it{points.l.begin()};
//        while ((it = std::find(it, points.l.end(), avgPoints.l.at(i))) != points.l.end())
//        {
//            auto idx{std::distance(points.l.begin(), it)};
//            tmpX.push_back(points.x.at(static_cast<size_t>(idx)));
//            tmpXerr.push_back(points.xErr.at(static_cast<size_t>(idx)));
//            tmpY.push_back(points.y.at(static_cast<size_t>(idx)));
//            tmpYerr.push_back(points.yErr.at(static_cast<size_t>(idx)));
//            it++;
//        }
//        avgPoints.x.push_back(std::accumulate(tmpX.begin(), tmpX.end(), 0.0) / static_cast<double>(tmpX.size()));
//        avgPoints.xErr.push_back(std::accumulate(tmpXerr.begin(), tmpXerr.end(), 0.0) / static_cast<double>(tmpXerr.size()));
//        avgPoints.y.push_back(std::accumulate(tmpY.begin(), tmpY.end(), 0.0) / static_cast<double>(tmpY.size()));
//        avgPoints.yErr.push_back(std::accumulate(tmpYerr.begin(), tmpYerr.end(), 0.0) / static_cast<double>(tmpYerr.size()));
//    }

//    writePointsToFile("output.txt", avgPoints);

//    std::string str;
//    useSub(avgPoints, str, value, true);
//    h2dConv.get()->SetTitle(str.c_str());

//    c.get()->Print(psName.c_str());

//    std::unique_ptr<TH2D> h2dConv1{new TH2D("h2dConv1",
//                                           ss.str().c_str(),
//                                           static_cast<int>(avgPoints.x.size()),
//                                           0.5,
//                                           0.5 + static_cast<int>(avgPoints.x.size()),
//                                           static_cast<int>(points.y.size()),
//                                           0.75 * (*std::min_element(points.y.begin(), points.y.end())),
//                                           1.25 * (*std::max_element(points.y.begin(), points.y.end())))};
//    h2dConv1->SetStats(0);
//    h2dConv1.get()->Draw();
//    std::string str1;
//    useSub1(avgPoints, str1, value, true);
//    h2dConv1.get()->SetTitle(str1.c_str());
//    c.get()->Print(psName.c_str());
//    c.get()->Print((psName + ']').c_str());
//    c.get()->Close();



////    auto useSub{true};
////    if (useSub)
////    {

////        std::map<std::pair<std::string, Color_t>, Points> subPoints{
////            { std::make_pair("b8", kRed), Points() },
////            { std::make_pair("check_w", kBlue), Points() },
////            { std::make_pair("p43", kGreen), Points() },
////            { std::make_pair("bereza_8", kOrange), Points() },
////            { std::make_pair("raspad", kCyan), Points() },
////            { std::make_pair("other", kBlack), Points() },
////        };


////        for (size_t i{0}; i < points.x.size(); ++i)
////        {
////            auto isOther{false};
////            for (auto &item : subPoints)
////            {
////                if (points.l.at(i).find(item.first.first) != std::string::npos)
////                {
////                    TMarker m{points.x.at(i), points.y.at(i), 21};
////                    m.SetMarkerSize(1.5);
////                    m.SetMarkerColor(item.first.second);
////                    m.DrawClone("SAME");
////                    item.second.l.push_back(points.l.at(i));
////                    item.second.x.push_back(points.x.at(i));
////                    item.second.y.push_back(points.y.at(i));
////                    item.second.xErr.push_back(0.1);
////                    item.second.yErr.push_back(0.5);

////                   isOther = true;
////                }
////            }
////            if (!isOther)
////            {
////                subPoints.at({"other", kMagenta}).l.push_back(points.l.at(i));
////                subPoints.at({"other", kMagenta}).x.push_back(points.x.at(i));
////                subPoints.at({"other", kMagenta}).y.push_back(points.y.at(i));
////                subPoints.at({"other", kMagenta}).xErr.push_back(0.1);
////                subPoints.at({"other", kMagenta}).yErr.push_back(0.5);
////            }
////        }
////        std::stringstream ss;
////        ss.str("");ss.clear();
////        ss << (value == Data1::Value::A ? "Ad" : "Wr");
////        ss << ": stdAbs=" << std::setprecision(3);
////        auto stdAbs1 = [](Points points){
////            std::vector<double> d2;
////            for (size_t i{0}; i < points.x.size(); ++i)
////            {
////                d2.push_back(std::pow(points.y.at(i) - points.x.at(i), 2));
////            }
////            return std::sqrt(std::accumulate(d2.begin(), d2.end(), 0.0) / static_cast<double>(d2.size()));
////        };
//////        auto avg1 = [](Points points){
//////            return std::accumulate(points.x.begin(), points.x.end(), 0.0) / points.x.size();
//////        };
////        for (auto item : subPoints)
////        {
////            ss << "[#color[" << static_cast<int>(item.first.second) << "]{" << stdAbs1(item.second) << "}] ";
////        }

////        ss << ";AGP-K, %;Chem, %";
////        h2dConv.get()->SetTitle(ss.str().c_str());
////    }

////    lConv.get()->Draw("SAME");
////    for (const auto &item : labels)
////    {
////        item.DrawClone("SAME");
////    }
////    c.get()->Print(psName.c_str());
////    c.get()->Print((psName + ']').c_str());
////    c.get()->Close();
//}

//struct Data { // TODO need this for calcRep
//    enum class Value {
//        A,
//        W
//    };
//    std::optional<double> a;
//    std::optional<double> w;
//};


//void calcRep(const std::shared_ptr<TF1> &f,
//             const std::map<std::string, Data> &chem,
//             const std::map<std::string, std::vector<std::string>> &vC,
//             const std::map<std::string, std::vector<std::string>> &vO) // TODO make to work with Data1
//{
//    auto p0{f.get()->GetParameter(0)};
//    auto p1{f.get()->GetParameter(1)};
//    auto p2{f.get()->GetParameter(2)};
//    auto p3{f.get()->GetParameter(3)};
//    auto p4{f.get()->GetParameter(4)};
//    auto p5{f.get()->GetParameter(5)};

//    std::vector<double> A;
//    std::vector<double> W;
//    for (const auto &item : chem)
//    {
//        if (item.second.a.has_value())
//        {
//            // std::cout << item.first << " " << vC.at(item.first).front() << std::endl;
//            auto C{strToDouble(vC.at(item.first).front())};
//            auto O{strToDouble(vO.at(item.first).front())};
//            auto a{
//                (p0 - p1 * C - p5 * (p2 * O + p4)) / (1 - p5 * p3)
//            };
//            auto w{
//                p2 * O - p3 * a + p4
//            };
//            A.push_back(a);
//            if (item.second.w.has_value())
//            {
//                W.push_back(w);
//            }
//        }
//    }

//    auto avgA{std::accumulate(A.begin(), A.end(), 0.0)};
//    avgA /= A.size();
//    auto avgW{std::accumulate(W.begin(), W.end(), 0.0)};
//    avgW /= W.size();

//    auto stdAabs{TMath::RMS(A.begin(), A.end())};
//    auto stdWabs{TMath::RMS(W.begin(), W.end())};
//    std::cout << "repeatability" << std::endl;
//    std::cout << "avgA = " << avgA << std::endl;
//    std::cout << "stdAabs = " << stdAabs << std::endl;
//    std::cout << "avgW = " << avgW << std::endl;
//    std::cout << "stdWabs = " << stdWabs << std::endl;
//}
