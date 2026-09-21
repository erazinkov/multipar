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

void findExcludedPoints(const std::vector<Point> &points);
void process(const std::vector<Point> &points, const ChemResult::Type &value);
double getPredicatedValueByType(const FitResult &fr, const ChemResult::Type &type, TF1 *f);
std::vector<Point> getPredicatedPointsByType(const std::map<std::string, Data> &data, const ChemResult::Type &type, TF1 *f);



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

double calculateCorr(const std::vector<Point>& points) {
    const size_t n{points.size()};
    if (n < 2) {
        throw std::invalid_argument("Need at least 2 points to calculate correlation");
    }

    double sumX{0.0};
    double sumY{0.0};
    for (const auto& p : points) {
        sumX += p.x;
        sumY += p.y;
    }
    const double meanX = sumX / n;
    const double meanY = sumY / n;

    double cov{0.0};
    double varX{0.0};
    double varY{0.0};
    for (const auto& p : points) {
        const double dx = p.x - meanX;
        const double dy = p.y - meanY;
        cov  += dx * dy;
        varX += dx * dx;
        varY += dy * dy;
    }

    const double denom{std::sqrt(varX * varY)};
    if (denom == 0.0) {
        throw std::runtime_error("Cannot compute correlation: zero variance");
    }
    return cov / denom;
}

std::pair<double, double> calculateAvgXY(const std::vector<Point>& points) {
    if (points.empty()) return {0.0, 0.0};
    auto avgX{0.0};
    auto avgY{0.0};
    for (const auto &item : points) {
        avgX += item.x;
        avgY += item.y;
    }
    avgX /= points.size();
    avgY /= points.size();
    return {avgX, avgY};
}

double calculateStdAbsRep(const std::vector<Point>& points) {
    if (points.empty()) return 0.0;
    auto avg = 0.0;
    for (const auto &item : points) {
        avg += item.x;
    }
    avg /= points.size();

    double sumSquaredDiff = 0.0;

    for (const auto &item : points) {
        const double diff = item.x - avg;
        sumSquaredDiff += diff * diff;
    }
    return std::sqrt(sumSquaredDiff / points.size());
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


void findExcludedPoints(const std::vector<Point> &points, const double limit = 1.0) {
    struct ExcludedPoint : Point {
        double d;
    };

    std::vector<ExcludedPoint> excludedPoints{};
    for (size_t i{0}; i < points.size(); ++i) {
        excludedPoints.push_back(ExcludedPoint{points.at(i), std::abs(points.at(i).x - points.at(i).y)});
    }
    excludedPoints.erase(std::remove_if(excludedPoints.begin(), excludedPoints.end(),
                                 [&limit](const ExcludedPoint& a) { return a.d <= limit; }), excludedPoints.end());
    std::cout << "It is proposed to exclude..." << std::endl;
    for (const auto &item : excludedPoints) {
        std::cout << item.sample << " " << std::endl;
    }
    std::cout << std::endl;
}

void process(const std::vector<Point> &points, const ChemResult::Type &value) {
    std::vector<Point> points_p;
    const std::string valueStr{value == ChemResult::Type::A ? "A" : "W"};
    points_p.insert(points_p.end(), points.cbegin(), points.cend());
    std::unique_ptr<TGraphErrors> gr_p{new TGraphErrors(points_p.size())};
    for (size_t i{0}; i < points_p.size(); i++) {
        gr_p.get()->SetPoint(i, points_p.at(i).x, points_p.at(i).y);
        gr_p.get()->SetPointError(i, points_p.at(i).xErr, points_p.at(i).yErr);
    }

    gr_p.get()->SetMarkerSize(1.0);
    gr_p.get()->SetMarkerStyle(21);

    auto itMin = std::max_element(points_p.begin(), points_p.end(),
                                  [](const Point& a, const Point& b) {
                                return a.y > b.y;
                                  });
    auto itMax = std::max_element(points_p.begin(), points_p.end(),
                               [](const Point& a, const Point& b) {
                                return a.y < b.y;
                               });
    auto min{(*itMin).y * 0.95};
    auto max{(*itMax).y * 1.05};

    std::unique_ptr<TH2D> h2d_p{new TH2D("h2d_p", "h2d_p", 100, min, max, 100, min, max)};
    h2d_p.get()->SetStats(0);

    auto doubleToString = [](double value, int precision = 2) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(precision) << value;
        return oss.str();
    };

    std::ostringstream ss;
    ss.str("");ss.clear();
    ss << "stdAbs=" << doubleToString(calculateStdAbsCon(points_p)) << ";" + valueStr + "_{m}', %;" + valueStr + "_{c}, %";
    h2d_p.get()->SetTitle(ss.str().c_str());

    const std::string psName_p{"output_p_" + valueStr + ".pdf"};
    std::unique_ptr<TCanvas> c_p{new TCanvas("c_p", "c_p", 1024, 960)};
    gPad->SetGrid();
    c_p.get()->Print((psName_p + '[').c_str());
    h2d_p.get()->Draw();
    std::unique_ptr<TLine> dLine{new TLine(min, min, max, max)};
    dLine.get()->Draw("SAME");
    gr_p.get()->Draw("P");

    std::vector<double> xAvg_p;
    for (const auto &item : points_p) {
        xAvg_p.push_back(item.x);
    }

    std::cout << "Avgs: " << calculateAvg(xAvg_p) << " " << "Number: " << xAvg_p.size() << std::endl;

//    findExcludedPoints(points_p, 1.0);

    std::map<std::pair<std::string, Color_t>, std::vector<Point>> subRanks{
        { std::make_pair(R"(ГБФ)", kBlue), {} }, // ГБФ
        { std::make_pair(R"(СКВП)", kOrange), {} }, // СКВП
        { std::make_pair(R"(ТИМАКС)", kGreen), {} }, // ТИМАКС
        { std::make_pair(R"(ТСГЦ)", kRed), {} }, // ТСГЦ
    };

    for (size_t i{0}; i < points_p.size(); ++i) {
        auto it = data_sample_rank.find(points_p.at(i).sample);
        if (it != data_sample_rank.end()) {
            for (auto &item : subRanks) {
                std::regex pattern(item.first.first);
                if (std::regex_search(it->second, pattern)) {
                    TMarker m{points_p.at(i).x, points_p.at(i).y, 21};
                    m.SetMarkerSize(1.5);
                    m.SetMarkerColor(item.first.second);
//                    m.DrawClone("SAME");
                    item.second.push_back(points_p.at(i));
                }
            }
        }
    }



    auto saveRanksToFile = [&](){
        const auto fileName{"ranks.txt"};
        std::ofstream ofs(fileName, std::ios::out);
        if (ofs.is_open()) {
            for (const auto &[key, value] : subRanks) {
                ofs << key.first << " ";
                for (const auto &item : value) {
                    ofs << item.sample << " ";
                }
                ofs << std::endl;
            }
            ofs.close();
        }
    };

//    saveRanksToFile();


    for (const auto &[key, value] : subRanks) {
        std::cout << key.first << " " << value.size() << std::endl;
    }




    std::map<std::pair<std::string, Color_t>, std::vector<Point>> subPoints{

        { std::make_pair(R"(sample(175|176_1|17[7-9]|18[0-7]|22[4-9]|23[0-2]|230_1|23[7-9]|24[0-3]|25[2-5]|264|267)\.)", kGreen), {} }, // grad
        { std::make_pair(R"(sample(27[4-9]|28[0-3]|287|290|293|29[6-9]|30[0-5]|308|311|313|31[5-6]|318|32[1-3]|32[5-9]|33[0-1]|334|336|33[8-9]|34[1-4]|34[8-9]|35[0-4]|357|363|36[5-6]|369|371|38[1-9]|39[0-5]|39[8-9]|40[0-1])\.)", kRed), {} },
    };

     auto sampleToLabel = [](const std::string &sample){
         auto pos{sample.find_first_of(".")};
         auto label{sample};
         if (pos != std::string::npos) {
             label.erase(pos, 1);
         }
         const std::string subStr{"sample"};
         pos = label.find(subStr);
         if (pos != std::string::npos) {
             label = label.erase(pos, subStr.length());
         }
         return label;
     };

     for (size_t i{0}; i < points_p.size(); ++i) {
//         auto it = data_sample_rank.find(points_p.at(i).sample);
//         if (it != data_sample_rank.end()) {
//             std::cout << it->second << std::endl;
//         }
         for (auto &item : subPoints) {
             std::regex pattern(item.first.first);
             if (std::regex_search(points_p.at(i).sample, pattern)) {
                 TLatex l(points_p.at(i).x, points_p.at(i).y + 1.25 * points_p.at(i).xErr, sampleToLabel(points_p.at(i).sample).c_str());
                 l.SetTextAngle(90);
                 l.SetTextAlign(12);
                 l.SetTextSize(0.02);
//                 l.DrawClone("SAME");
                 TMarker m{points_p.at(i).x, points_p.at(i).y, 21};
                 m.SetMarkerSize(1.5);
                 m.SetMarkerColor(item.first.second);
                 m.DrawClone("SAME");
                 item.second.push_back(points_p.at(i));
             }
         }
     }

    // choose points by idx
    std::map<std::pair<std::string, Color_t>, std::vector<Point>> subPoints_{
        {std::make_pair("grad", kGreen), {} },
        {std::make_pair("check", kRed), {} }
    };

    for (size_t i{0}; i < points_p.size(); ++i) {
        auto color{kBlack};
        if (points_p.at(i).idx < points_p.size() / 2) {
            subPoints_.at({"grad", kGreen}).push_back(points_p.at(i));
            color = kGreen;
        } else {
            subPoints_.at({"check", kRed}).push_back(points_p.at(i));
            color = kRed;
        }
        TLatex l(points_p.at(i).x, points_p.at(i).y + 1.25 * points_p.at(i).xErr, sampleToLabel(points_p.at(i).sample).c_str());
        l.SetTextAngle(90);
        l.SetTextAlign(12);
        l.SetTextSize(0.02);
//        l.DrawClone("SAME");
        TMarker m{points_p.at(i).x, points_p.at(i).y, 21};
        m.SetMarkerSize(1.5);
        m.SetMarkerColor(color);
        m.DrawClone("SAME");
    }


    std::map<std::pair<std::string, Color_t>, std::map<std::string, double>> subStats;
    for (auto &item : subPoints_) {
        subStats[item.first] = {};
        subStats.at(item.first).insert({"stdAbs", calculateStdAbsCon(item.second)});
        auto avgXY{calculateAvgXY(item.second)};
        subStats.at(item.first).insert({"avgX", avgXY.first});
        subStats.at(item.first).insert({"avgY", avgXY.second});
        subStats.at(item.first).insert({"r", calculateCorr(item.second)});
    }
    if (!subStats.empty()) {
        double yRow = 0.85;
        double xCol = 0.15;
        double colWidth = 0.1;
        auto it = subStats.begin();
        for (const auto &statsItem : (*it).second) {
            TLatex *t = new TLatex(xCol, yRow, statsItem.first.c_str());
            t->SetTextAlign(22);   // 2 = center horizontally, 2 = center vertically
            t->SetTextSize(0.03);
            t->SetNDC();
            t->Draw();
            xCol += colWidth;
        }
        yRow -= 0.05;
        for (const auto& item : subStats) {
            double xCol = 0.15;
            Color_t color = item.first.second;
            std::map<std::string, double> stats = item.second;
            for (const auto &statsItem : stats) {
                TLatex *t = new TLatex(xCol, yRow, doubleToString(statsItem.second).c_str());
                t->SetTextColor(color);
                t->SetTextAlign(22);   // 2 = center horizontally, 2 = center vertically
                t->SetTextSize(0.03);
                t->SetNDC();
                t->Draw();
                xCol += colWidth;
            }
            yRow -= 0.05;
        }
    }


//    std::unique_ptr<TPaveText> pt{new TPaveText(0.1, 0.65, 0.5, 0.9, "NDC")};
//    pt.get()->SetFillColor(0);
//    pt.get()->SetBorderSize(1);

//    if (!subStats.empty()) {
//        std::ostringstream oss;
//        auto it = subStats.begin();
//        for (const auto &statsItem : (*it).second) {
//            oss << std::left << std::setw(8) << statsItem.first;
//        }
//        pt->AddText(oss.str().c_str());
//    }

//        for (const auto& item : subStats) {
//            std::string key = "";
//            Color_t color = item.first.second;
//            std::map<std::string, double> stats = item.second;
//            for (const auto &statsItem : stats) {
//    //            key.append(statsItem.first);
//    //            key.append("=");
//                key.append(doubleToString(statsItem.second));
//                key.append(" ");
//            }

//            size_t pos = key.find_last_not_of(" \t\n\r\f\v");
//            if (pos != std::string::npos) {
//                key.erase(pos + 1);
//            } else {
//                key.clear();
//            }
//            TText *text = pt.get()->AddText(key.c_str());
//            text->SetTextColor(color);
//        }

//    if (!subStats.empty()) {
//        std::string key = "";
//        auto it = subStats.begin();
//        std::map<std::string, double> stats = (*it).second;
//        for (const auto &statsItem : stats) {
//            key.append(statsItem.first);
//            key.append(";");
//        }
//        pt.get()->AddText(key.c_str());
//    }


//    for (const auto& item : subStats) {
//        std::string key = "";
//        Color_t color = item.first.second;
//        std::map<std::string, double> stats = item.second;
//        for (const auto &statsItem : stats) {
////            key.append(statsItem.first);
////            key.append("=");
//            key.append(doubleToString(statsItem.second));
//            key.append(" ");
//        }

//        size_t pos = key.find_last_not_of(" \t\n\r\f\v");
//        if (pos != std::string::npos) {
//            key.erase(pos + 1);
//        } else {
//            key.clear();
//        }
//        TText *text = pt.get()->AddText(key.c_str());
//        text->SetTextColor(color);
//    }
//    pt.get()->Draw("SAME");

    c_p.get()->Print(psName_p.c_str());
    c_p.get()->Print((psName_p + ']').c_str());
    c_p.get()->Close();
}

double getPredicatedValueByType(const FitResult &fr, const ChemResult::Type &type, TF1 *f) {
    auto val_a{0.0};
    auto val_w{0.0};
    val_a = ( f->GetParameter(3)
             - f->GetParameter(4) * f->GetParameter(0) * fr.getElementResultByName("O").value
             - f->GetParameter(2) * f->GetParameter(4)
             - f->GetParameter(5) * fr.getElementResultByName("C").value
             - f->GetParameter(6) * fr.getElementResultByName("N").value
            )
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


std::vector<Point> getPredicatedPointsByType(const std::map<std::string, Data> &data, const ChemResult::Type &type, TF1 *f) {
    std::vector<Point> points;
    for (const auto &[key, value] : data) {
        for (const auto &fr : value.fitResults) {
            std::optional<double> v;
            double err{0.0};
            switch (type) {
            case (ChemResult::Type::A):
                v = value.chemResult.a;
                err = 0.01;
                break;
            case (ChemResult::Type::W):
                v = value.chemResult.w;
                err = 0.01;
                break;
            }
            if (v.has_value()) {
                Point point;
                point.idx = value.idx;
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

std::vector<std::string>splitLineToStrs(const std::string &line);
double strToDouble(std::string str);

int main()
{
    TVirtualFitter::SetDefaultFitter("Minuit");

    std::map<std::string, ChemResult> chem{};
//    chem.insert(data_chem_cat_1.begin(), data_chem_cat_1.end());
//    chem.insert(data_chem_cat_3.begin(), data_chem_cat_3.end());
    chem.insert(data_chem_cat_4.begin(), data_chem_cat_4.end());
//    chem.insert(data_chem_cat_5.begin(), data_chem_cat_5.end());
//    chem.insert(data_chem_cat_6.begin(), data_chem_cat_6.end());

    const std::map<int, std::string> columnElement
    {
         {1, "Al"},
         {2, "C"},
         {3, "N"},
         {4, "O"},
         {5, "Si"},
    };

//    const std::map<int, std::string> columnElement
//    {
//         {1, "Al"},
//         {3, "C"},
//         {5, "O"},
//         {7, "Si"},
//    };


    const auto fileName{"OF_data.subcat.csv"};

    std::cout << fileName << std::endl;

    auto splitLineToStrs_ = [](const std::string& line) {
        std::vector<std::string> result;
        std::stringstream ss(line);
        std::string field;
        while (std::getline(ss, field, '\t')) {
            result.push_back(field);
        }

        if (!line.empty() && line.back() == '\t') {
            result.push_back("");
        }

        return result;
    };

    // exclude
    std::vector<std::string> excludeSamples{
//        R"(sample28(0|1|2)\.)",
//        R"(sample447_\d\.)",
    };

    auto isExclude = [](const std::string &sample, std::vector<std::string> &excludeSamples){
        bool exclude{false};
        for (auto &eS : excludeSamples) {
            std::regex pattern(eS);
            if (std::regex_search(sample, pattern)) {
                exclude = true;
            }
        }
        return exclude;
    };


    auto getData_ = [&](const std::string &fileName, const std::map<int, std::string> &columnElement) {
        std::ifstream ifs(fileName);
        if (!ifs.is_open()) {
            throw my_error("Can't open file \"" + fileName + "\"");
        }
        std::string line;

        std::vector<Data> data;
        try {
            getline(ifs, line);
            while (getline(ifs, line)) {
                auto strs{splitLineToStrs_(line)};
                if (isExclude(strs.at(0), excludeSamples)) {
                    continue;
                }
                Data d;
                FitResult fR;
                for (const auto &[key, value] : columnElement)
                {
                    fR.elementResults.push_back({value,
                                            strToDouble(strs.at(static_cast<unsigned int>(key))),
                                            strToDouble(strs.at(static_cast<unsigned int>(key + 1)))
                                           });
                }
                d.sample = strs.at(0);
                d.category = static_cast<uint>(std::stoul(strs.at(9)));
                auto str_a{strs.at(10)};
                d.chemResult.a = !str_a.empty() ? std::optional<double>(strToDouble(str_a)) : std::nullopt;
                auto str_w{strs.at(11)};
                d.chemResult.w = !str_w.empty() ? std::optional<double>(strToDouble(str_w)) : std::nullopt;
                d.subCategory = static_cast<uint>(std::stoul(strs.at(12)));
                d.fitResults.push_back(fR);
                data.push_back(d);
            }

        }  catch (...) {
            std::cout << "Error reading data from " << fileName;
        }
        ifs.close();
        return data;
    };


    auto data_{getData_(fileName, columnElement)};

    // category
    std::map<std::string, Data> data;
    uint idx{0};
    for (const auto& d : data_) {
        if (d.subCategory == 41) {
            data[d.sample] = d;
            data[d.sample].idx = idx;
            idx++;
        }
    }
    // choose grad by idx
    std::map<std::string, Data> data_grad;
    for (const auto& [key, value] : data) {
        if (value.idx < data.size() / 2) {
            data_grad[key] = value;
        }
    }

//    try {
//        std::regex m{R"(\bsample([1-9]|[12][0-9]|30)\b)"}; //30
//        std::regex m{R"((sample(?:[1-4]|9|10|3[7-9]|4[0-9]|5[0-3]|5[8-9]|6[0-9]|72|7[5-9]|8[0-5])\.))"}; // grad
//        std::regex m(R"((sample(?:10[6-9]|17[2-9]|18[0-7])(?:_1)?\.))"); // check
//        std::regex m{R"(sample\d+\.)"};
//        std::regex m{R"(sample(1(7[6-9]|8[0-9]|9[0-9])|2([0-5][0-9]|6[0-7]))\.)"};
        std::regex m{R"(sample(175|176_1|17[7-9]|18[0-7]|22[4-9]|23[0-2]|230_1|23[7-9]|24[0-3]|25[2-5]|264|267)\.)"};
//        auto data{getData(fileName, columnElement, chem, m)};
//        for (const auto &[key, value] : data_grad) {
//            std::cout << key << " ";
//            value.print();
//            std::cout << std::endl;
//        }

        std::vector<Point> points_a{getPointsByType(data_grad, ChemResult::Type::A)};
        std::vector<Point> points_w{getPointsByType(data_grad, ChemResult::Type::W)};
        for (auto &p : points_w) {
            p.x  = p.x + points_a.size();
        }

        std::vector<Point> points;
        points.insert(points.end(), points_a.cbegin(), points_a.cend());
        points.insert(points.end(), points_w.cbegin(), points_w.cend());

        for (const auto &p : points) {
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

//        f.get()->FixParameter(6, 0.0);

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

        // std::regex m_a{R"(pulp_rot_berez_7_w\d+_sum|pulp_rot_berez_2_w\d+_sum|pulp_rot_berez_11_w\d+_sum|pulp_rot_berez_7_w\d+p\d+_sum)"};
//        std::regex m_a{R"((sample(?:[1-4]|9|10|3[7-9]|4[0-9]|5[0-3]|5[8-9]|6[0-9]|72|7[5-9]|8[0-5])\.))"};
        std::regex m_a{R"(sample\d+)"};
//        std::regex m_a{R"(sample(1(7[6-9]|8[0-9]|9[0-9])|2([0-9][0-9])|3([0-3][0-9]|4[0-4]))\.)"};


        std::vector<Point> points_p_a{getPredicatedPointsByType(data, ChemResult::Type::A, f.get())};
        std::vector<Point> points_p_w{getPredicatedPointsByType(data, ChemResult::Type::W, f.get())};

//        auto data_a{getData(fileName, columnElement, chem, m_a)};

//        std::vector<Point> points_p_a{getPredicatedPointsByType(data_a, ChemResult::Type::A, f.get())};
//        std::vector<Point> points_p_w{getPredicatedPointsByType(data_a, ChemResult::Type::W, f.get())};

        process(points_p_a, ChemResult::Type::A);
        process(points_p_w, ChemResult::Type::W);

//        std::regex m_r{R"(sample(115|115_\d+))"};

//        auto data_r{getData(fileName, columnElement, data_chem_cat_5, m_r)};

//        for (const auto &[key, value] : data_r) {
//            std::cout << key << " ";
//            value.print();
//            std::cout << std::endl;
//        }

//        std::vector<Point> points_r_a{getPredicatedPointsByType(data_r, ChemResult::Type::A, f.get())};
//        std::vector<Point> points_r_w{getPredicatedPointsByType(data_r, ChemResult::Type::W, f.get())};

//        std::cout << "RepA=" << calculateStdAbsRep(points_r_a) << ::std::endl;
//        std::cout << "RepW=" << calculateStdAbsRep(points_r_w) << ::std::endl;

//    }
//    catch (const my_error& err)
//    {
//        std::cout << "Error: " << err.what() << std::endl;
//    }
//    catch (const std::exception& err)
//    {
//        std::cout << "Error: " << err.what() << std::endl;
//    }
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
