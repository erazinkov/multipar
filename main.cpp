#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <regex>
#include <memory>
#include <cmath>
#include <optional>

#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLatex.h>
#include <TH2.h>
#include <TLine.h>
#include <TMarker.h>
#include <TMath.h>
#include <TVirtualFitter.h>

// ============================================================================
// Exception Classes
// ============================================================================

class DataProcessingError : public std::exception {
public:
    explicit DataProcessingError(const std::string& message) : message_(message) {}
    const char* what() const noexcept override { return message_.c_str(); }
private:
    std::string message_;
};

// ============================================================================
// Data Structures
// ============================================================================

struct FitResult {
    std::string element;
    double value;
    double error;

    void print() const {
        std::cout << element << " " << value << "\u00B1" << error << " ";
    }
};

struct ChemicalComposition {
    std::optional<double> a;  // Some parameter A
    std::optional<double> w;  // Some parameter W

    void print() const {
        std::cout << a.value_or(-1) << " " << w.value_or(-1);
    }
};

struct SampleData {
    ChemicalComposition chem;
    std::vector<std::vector<FitResult>> fitResults;

    void print() const {
        chem.print();
        std::cout << " ";
        for (const auto& resultSet : fitResults) {
            for (const auto& result : resultSet) {
                result.print();
            }
        }
        std::cout << std::endl;
    }
};

struct DataPoints {
    std::vector<std::string> labels;
    std::vector<double> x;
    std::vector<double> xError;
    std::vector<double> y;
    std::vector<double> yError;
};

// ============================================================================
// Fit Function Class
// ============================================================================

class MultiElementFitFunction {
public:
    explicit MultiElementFitFunction(std::map<int, std::vector<FitResult>> fitData)
        : fitData_(std::move(fitData)) {}

    double operator()(double* x, double* parameters) const {
        const int index = std::min(static_cast<int>(std::round(x[0])),
                                   static_cast<int>(fitData_.size() - 1));
        const size_t numParams = fitData_.begin()->second.size();

        double result = 0.0;
        for (size_t i = 0; i < numParams; ++i) {
            result += parameters[i] * fitData_.at(index).at(i).value;
        }
        result += parameters[numParams];
        return result;
    }

private:
    std::map<int, std::vector<FitResult>> fitData_;
};

// ============================================================================
// Utility Functions
// ============================================================================

std::vector<std::string> splitLine(const std::string& line) {
    std::stringstream ss(line);
    std::string token;
    std::vector<std::string> tokens;
    while (ss >> token) {
        tokens.push_back(token);
    }
    return tokens;
}

double stringToDouble(const std::string& str) {
    try {
        return std::stod(str);
    } catch (const std::exception&) {
        throw DataProcessingError("Cannot convert to double: " + str);
    }
}

double calculateStdAbs(const DataPoints& points) {
    if (points.x.empty()) return 0.0;

    double sumSquaredDiff = 0.0;
    for (size_t i = 0; i < points.x.size(); ++i) {
        const double diff = points.y[i] - points.x[i];
        std::cout << points.y[i] << " " << points.x[i] << std::endl;
        sumSquaredDiff += diff * diff;
    }
    return std::sqrt(sumSquaredDiff / points.x.size());
}

double calculateAverage(const std::vector<double>& values) {
    if (values.empty()) return 0.0;
    return std::accumulate(values.begin(), values.end(), 0.0) / values.size();
}

// ============================================================================
// Data Processing Functions
// ============================================================================

std::map<std::string, SampleData> loadFitResults(
    const std::string& filename,
    const std::map<int, std::string>& columnElements,
    const std::map<std::string, ChemicalComposition>& chemicalData,
    const std::regex& pattern)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw DataProcessingError("Cannot open file: " + filename);
    }

    std::map<std::string, SampleData> results;
    std::string line;

    while (std::getline(file, line)) {
        const auto tokens = splitLine(line);

        auto match = std::find_if(chemicalData.begin(), chemicalData.end(),
            [&tokens, &pattern](const auto& entry) {
                return tokens[0].find(entry.first) != std::string::npos &&
                       std::regex_search(tokens[0], pattern);
            });

        if (match == chemicalData.end()) continue;

        std::vector<FitResult> fitResults;
        for (const auto& [colIndex, elementName] : columnElements) {
            fitResults.push_back({
                elementName,
                stringToDouble(tokens[colIndex]),
                stringToDouble(tokens[colIndex + 1])
            });
        }

        const auto& sampleName = match->first;
        results[sampleName].chem = match->second;
        results[sampleName].fitResults.push_back(fitResults);
    }

    return results;
}

std::map<int, std::vector<FitResult>> extractFitResultsByValue(
    const std::map<std::string, SampleData>& data,
    bool useAValue)
{
    std::map<int, std::vector<FitResult>> extractedResults;
    int index = 0;

    for (const auto& [sampleName, sampleData] : data) {
        for (const auto& resultSet : sampleData.fitResults) {
            std::optional<double> value = useAValue ? sampleData.chem.a : sampleData.chem.w;
            if (value.has_value()) {
                extractedResults[index++] = resultSet;
            }
        }
    }

    return extractedResults;
}

void addPointsByValue(
    const std::map<std::string, SampleData>& data,
    DataPoints& points,
    bool useAValue)
{
    int index = 0;
    for (const auto& [sampleName, sampleData] : data) {
        for (size_t i = 0; i < sampleData.fitResults.size(); ++i) {
            std::optional<double> value = useAValue ? sampleData.chem.a : sampleData.chem.w;
            if (value.has_value()) {
                points.labels.push_back(sampleName + "_" + std::to_string(i));
                points.x.push_back(index++);
                points.xError.push_back(0.01);
                points.y.push_back(value.value());
                points.yError.push_back(0.5);
            }
        }
    }
}

void writePointsToFile(const std::string& filename, const DataPoints& points) {
    std::ofstream file(filename);
    if (!file.is_open()) return;

    for (size_t i = 0; i < points.labels.size(); ++i) {
        file << points.x[i] << " " << points.y[i] << " " << points.labels[i] << std::endl;
    }
}

// ============================================================================
// Visualization Functions
// ============================================================================

std::string createSubGroupLegend(const DataPoints& points) {
    std::map<std::pair<std::string, Color_t>, DataPoints> subgroups = {
        {{"N12", kRed}, DataPoints()},
        {{"berez", kBlue}, DataPoints()},
        {{"kuz", kGreen}, DataPoints()},
        {{"bereza_8", kOrange}, DataPoints()},
        {{"raspad", kCyan}, DataPoints()},
        {{"other", kBlack}, DataPoints()}
    };

    // Categorize points into subgroups
    for (size_t i = 0; i < points.x.size(); ++i) {
        bool found = false;
        for (auto& [key, subgroup] : subgroups) {
            if (points.labels[i].find(key.first) != std::string::npos) {
                subgroup.labels.push_back(points.labels[i]);
                subgroup.x.push_back(points.x[i]);
                subgroup.y.push_back(points.y[i]);
                subgroup.xError.push_back(0.1);
                subgroup.yError.push_back(0.5);
                found = true;
                break;
            }
        }
        if (!found) {
            subgroups[{"other", kBlack}].labels.push_back(points.labels[i]);
            subgroups[{"other", kBlack}].x.push_back(points.x[i]);
            subgroups[{"other", kBlack}].y.push_back(points.y[i]);
            subgroups[{"other", kBlack}].xError.push_back(0.1);
            subgroups[{"other", kBlack}].yError.push_back(0.5);
        }
    }

    // Build legend string
    std::stringstream legend;
    for (const auto& [key, subgroup] : subgroups) {
        const double stdAbs = calculateStdAbs(subgroup);
        if (!std::isnan(stdAbs)) {
            legend << "[#color[" << static_cast<int>(key.second) << "]{"
                   << stdAbs << "}] ";
        }
    }
    return legend.str();
}

void drawLabels(const DataPoints& points, double offset = 1.25) {
    for (size_t i = 0; i < points.x.size(); ++i) {
        TLatex label(points.x[i], points.y[i] + offset * points.yError[i],
                     points.labels[i].c_str());
        label.SetTextAngle(90);
        label.SetTextAlign(12);
        label.SetTextSize(0.02);
        label.DrawClone("SAME");
    }
}

// ============================================================================
// Calculation Functions
// ============================================================================

void calculateConvergence(
    const std::map<std::string, SampleData>& data,
    const std::unique_ptr<TF1>& fitFunction,
    bool useAValue)
{
    DataPoints points;

    // Calculate predicted values
    for (const auto& [sampleName, sampleData] : data) {
        for (const auto& resultSet : sampleData.fitResults) {
            std::optional<double> chemValue = useAValue ? sampleData.chem.a : sampleData.chem.w;
            if (!chemValue.has_value()) continue;

            const int numParams = fitFunction->GetNpar();
            double predicted = 0.0;
            for (int p = 0; p < numParams - 1; ++p) {
                predicted += fitFunction->GetParameter(p) * resultSet[p].value;
            }
            predicted += fitFunction->GetParameter(numParams - 1);

            points.labels.push_back(sampleName);
            points.x.push_back(predicted);
            points.y.push_back(chemValue.value());
            points.xError.push_back(0.1);
            points.yError.push_back(0.5);
        }
    }

    // Calculate statistics
    const double stdAbs = calculateStdAbs(points);
    const double avgPredicted = calculateAverage(points.x);

    std::cout << "Convergence: avg = " << avgPredicted
              << ", stdAbs = " << stdAbs << std::endl;

    // Create convergence plot
    const auto minY = 0.5 * *std::min_element(points.y.begin(), points.y.end());
    const auto maxY = 1.5 * *std::max_element(points.y.begin(), points.y.end());

    auto h2dConv = std::make_unique<TH2D>("h2dConv",
        (std::string(useAValue ? "Ad" : "Wr") +
         ": avg = " + std::to_string(avgPredicted) + "%" +
         ", stdAbs = " + std::to_string(stdAbs) + "%" +
         ";AGP-K, %;Chem, %").c_str(),
        points.y.size(), minY, maxY,
        points.y.size(), minY, maxY);

    h2dConv->SetStats(0);

    auto line = std::make_unique<TLine>(minY, minY, maxY, maxY);

    TCanvas canvas("convCanvas", "Convergence Plot", 1024, 960);
    canvas.Print("output_conv.ps[");

    h2dConv->Draw();
    line->Draw("SAME");

    // Draw subgroup markers
    for (size_t i = 0; i < points.x.size(); ++i) {
        TMarker marker(points.x[i], points.y[i], 21);
        marker.SetMarkerSize(1.5);

        // Color by subgroup
        if (points.labels[i].find("N12") != std::string::npos)
            marker.SetMarkerColor(kRed);
        else if (points.labels[i].find("berez") != std::string::npos)
            marker.SetMarkerColor(kBlue);
        else if (points.labels[i].find("kuz") != std::string::npos)
            marker.SetMarkerColor(kGreen);
        else if (points.labels[i].find("bereza_8") != std::string::npos)
            marker.SetMarkerColor(kOrange);
        else if (points.labels[i].find("raspad") != std::string::npos)
            marker.SetMarkerColor(kCyan);
        else
            marker.SetMarkerColor(kBlack);

        marker.DrawClone("SAME");
    }

    drawLabels(points);

    canvas.Print("output_conv.ps");
    canvas.Print("output_conv.ps]");
    canvas.Close();
}

void calculateRepeatability(
    const std::map<std::string, SampleData>& data,
    const std::unique_ptr<TF1>& fitFunction)
{
    std::vector<double> results;

    for (const auto& [sampleName, sampleData] : data) {
        std::cout << sampleName << " " << sampleData.fitResults.size() << std::endl;

        for (const auto& resultSet : sampleData.fitResults) {
            const int numParams = fitFunction->GetNpar();
            double predicted = 0.0;
            for (int p = 0; p < numParams - 1; ++p) {
                predicted += fitFunction->GetParameter(p) * resultSet[p].value;
            }
            predicted += fitFunction->GetParameter(numParams - 1);
            results.push_back(predicted);
        }
    }

    const double avg = calculateAverage(results);
    const double stdAbs = TMath::RMS(results.begin(), results.end());

    std::cout << "Repeatability: avg = " << avg
              << ", stdAbs = " << stdAbs << std::endl;
}

DataPoints getPredictedPoints(
    const std::map<std::string, SampleData>& data,
    const std::unique_ptr<TF1>& fitFunction,
    bool useAValue);

// ============================================================================
// Main Application
// ============================================================================

int main() {
    try {
        // Configuration
        const std::map<int, std::string> columnElements = {
            {1, "Al"},
            {3, "C"},
            {5, "N"},
            {7, "O"},
            {9, "Si"}
        };

        const std::string inputFile = "rea.elts.stroy.txt";

        // Chemical composition data for different sample types
        std::map<std::string, ChemicalComposition> chemicalData = {
            // Standard coal samples
            {"std_coal_proba_1_", {8.2, 5.2}},
            {"std_coal_proba_2_", {8.0, 5.3}},
            {"std_coal_proba_4_", {9.0, 5.5}},
            {"std_coal_proba_5_", {6.9, 4.9}},
            {"std_coal_proba_7_", {8.6, 6.1}},
            {"std_coal_proba_8_", {9.2, 6.5}},
            {"std_coal_proba_9_", {7.5, 5.0}},
            {"std_coal_proba_10_", {6.9, 5.6}},
            {"std_coal_proba_11_", {7.3, 5.7}},
            {"std_coal_proba_12_", {7.35, 5.5}},
            {"std_coal_proba_13_", {7.55, 5.3}},
            {"std_coal_proba_14_", {7.5, 4.9}},
            {"std_coal_proba_15_", {8.6, 7.2}},
            {"std_coal_proba_16_", {10.1, 5.6}},

            // Coal check samples with varying water content
            {"coal_check_w_1p35_", {7.0, 1.35}},
            {"coal_check_w_5p0_", {7.3, 5.0}},
            {"coal_check_w_10p0_", {7.6, 10.0}},
            {"coal_check_w_15p0_", {8.6, 15.0}},
            {"coal_check_w_20p0_", {9.2, 20.0}},

            // Bereza samples
            {"coal_check_bereza_8_w_0p8", {24.5, 0.8}},
            {"coal_check_bereza_8_w_5p0_", {24.5, 5.0}},
            {"coal_check_bereza_8_w_10p0_", {24.5, 10.0}},
            {"coal_check_bereza_8_w_25p0_", {24.5, 25.0}},

            {"coal_check_bereza_9_w_0p0_", {27.0, 0.9}},
            {"coal_check_bereza_9_w_5p0_", {27.0, 5.0}},
            {"coal_check_bereza_9_w_10p0_", {27.0, 10.0}},
            {"coal_check_bereza_9_w_15p0_", {27.0, 15.0}},

            // Pulp rotor samples
            {"pulp_rot_N12_1_", {7.8, 4.2}},
            {"pulp_rot_N12_2_", {9.6, 5.5}},
            {"pulp_rot_N12_3_", {11.2, 6.2}},
            {"pulp_rot_N12_4_", {11.8, 3.9}},
            {"pulp_rot_N12_5_", {15.1, 7.9}},
            {"pulp_rot_N12_6_", {18.2, 4.9}},
            {"pulp_rot_N12_7_", {20.7, 6.7}},
            {"pulp_rot_N12_8_", {27.6, 8.0}},
            {"pulp_rot_N12_9_", {28.3, 7.8}},
            {"pulp_rot_N12_10_", {30.4, 8.2}},
            {"pulp_rot_N12_11_", {32.9, 8.1}},

            {"pulp_rot_berez_2_", {15.6, 0.9}},
            {"pulp_rot_berez_7_", {19.8, 0.8}},
            {"pulp_rot_berez_11_", {24.2, 0.8}},

            {"pulp_rot_kuz_1_a7p85_", {7.85, 0.5}},
            {"pulp_rot_kuz_6_a18p48_", {18.48, 0.5}},
            {"pulp_rot_kuz_7_a27p12_", {27.12, 0.5}},
            {"pulp_rot_kuz_8_a5p89_", {5.89, 0.5}},
            {"pulp_rot_kuz_9_a5p76_", {5.76, 0.5}},

            {"pulp_rot_berez_6_w5_", {15.5, 5.0}},
            {"pulp_rot_berez_6_w10_", {15.5, 10.0}},
            {"pulp_rot_berez_6_w15_", {15.5, 15.0}},
            {"pulp_rot_berez_6_w20_", {15.5, 20.0}},

            {"pulp_rot_berez_11_w5_", {24.2, 5.0}},
            {"pulp_rot_berez_11_w10_", {24.2, 10.0}},
            {"pulp_rot_berez_11_w15_", {24.2, 15.0}},
        };

        // Load data
//        std::regex pattern(R"((pulp_rot_N12_\d+_\d+))");
        std::regex pattern(R"((pulp_rot_berez_(11|6)_w\d+_\d+))");
        auto fitData = loadFitResults(inputFile, columnElements, chemicalData, pattern);

        // Build points and fit
        DataPoints points;
        const bool useAValue = true;  // Use 'A' values
        addPointsByValue(fitData, points, useAValue);

        auto fitResultsByValue = extractFitResultsByValue(fitData, useAValue);
        std::cout << "Number of fit results: " << fitResultsByValue.size() << std::endl;

        // Create graph
        auto graph = std::make_unique<TGraphErrors>(
            points.x.size(), &points.x[0], &points.y[0],
            &points.xError[0], &points.yError[0]);
        graph->SetMarkerSize(1.5);
        graph->SetMarkerStyle(21);
        graph->SetTitle((";N_probe;" + std::string(useAValue ? "A" : "W")).c_str());

        // Create fit function
        MultiElementFitFunction fitObj(fitResultsByValue);
        auto fitFunction = std::make_unique<TF1>("fit", fitObj,
            points.x.front(), points.x.back(),
            columnElements.size() + 1);

        // Set parameter limits for W value
        if (!useAValue) {
            fitFunction->SetParLimits(0, -10.0, 0.0);
            fitFunction->SetParLimits(1, -10.0, 0.0);
            fitFunction->SetParLimits(2, -10.0, 0.0);
            fitFunction->SetParLimits(3, 0.0, 10.0);
            fitFunction->SetParLimits(4, -10.0, 0.0);
        }

        fitFunction->SetNpx(10 * points.x.size());
        graph->Fit(fitFunction.get(), "R");

        // Draw plot
        TCanvas canvas("mainCanvas", "Main Plot", 1024, 960);
        canvas.SetGrid();
        canvas.Print("output.ps[");
        graph->Draw("APL");

        // Draw labels
        drawLabels(points);

        // Draw subgroups with different colors
        const std::vector<std::pair<std::string, Color_t>> subgroups = {
            {"N12", kRed},
            {"berez", kBlue},
            {"kuz", kGreen},
            {"bereza_8", kOrange},
            {"raspad", kCyan}
        };

        for (size_t i = 0; i < points.x.size(); ++i) {
            Color_t color = kBlack;
            for (const auto& [name, col] : subgroups) {
                if (points.labels[i].find(name) != std::string::npos) {
                    color = col;
                    break;
                }
            }

            TMarker marker(points.x[i], points.y[i], 21);
            marker.SetMarkerSize(1.5);
            marker.SetMarkerColor(color);
            marker.DrawClone("SAME");
        }

        canvas.Print("output.ps");
        canvas.Print("output.ps]");
        canvas.Close();

        // Calculate convergence for summed data
//        std::regex sumPattern(R"((pulp_rot_N12_\d+_sum|pulp_rot_berez_\d+_sum|pulp_rot_kuz_\d+_a\d+p\d+_sum))");
//        std::regex sumPattern(R"((pulp_rot_N12_\d+_sum))");
//        auto sumData = loadFitResults(inputFile, columnElements, chemicalData, sumPattern);
//        calculateConvergence(sumData, fitFunction, useAValue);

        struct Results {
            std::string title;
            DataPoints predictedPoints;
            double stdAbs;
            Color_t color;
        };

        auto getResults = [&inputFile, &columnElements, &chemicalData, &fitFunction](const std::string &title,
                const std::regex &pattern,
                const bool useAValue, const Color_t color){
            Results r;
            r.title = title;
            auto fitResults = loadFitResults(inputFile, columnElements, chemicalData, pattern);
            for (const auto &[key, value] : fitResults) {
                std::cout << key << " ";
                value.print();
            }
            r.predictedPoints = getPredictedPoints(fitResults, fitFunction, useAValue);
            r.stdAbs = calculateStdAbs(r.predictedPoints);
            r.color = color;
            return r;
        };

        std::vector<std::tuple<std::string, std::regex, Color_t>> preResults{
//            {"N12", std::regex(R"((pulp_rot_N12_\d+_sum))"), kRed},
//            {"berez", std::regex(R"((pulp_rot_berez_\d+_sum))"), kBlue},
//            {"kuz", std::regex(R"((pulp_rot_kuz_\d+_a\d+p\d+_sum))"), kGreen},
            {"berez_11_w", std::regex(R"((pulp_rot_berez_11_w\d+_\d+))"), kOrange},
            {"berez_6_w", std::regex(R"((pulp_rot_berez_6_w\d+_\d+))"), kMagenta},
        };
        std::vector<Results> results;
        for (const auto &r : preResults) {
            auto [title, pattern, color] = r;
            auto rr = getResults(title, pattern, useAValue, color);
            results.push_back(rr);
        }

        const auto minY = 0.5 * *std::min_element(points.y.begin(), points.y.end());
        const auto maxY = 1.5 * *std::max_element(points.y.begin(), points.y.end());

        std::string h2dConvTitle = useAValue ? "Ad" : "Wr";
        auto h2dConv = std::make_unique<TH2D>("h2dConv", h2dConvTitle.c_str(), points.y.size(), minY, maxY, points.y.size(), minY, maxY);

        h2dConv->SetStats(0);

        auto line = std::make_unique<TLine>(minY, minY, maxY, maxY);

        TCanvas convCanvas("convCanvas", "Convergence Plot", 1024, 960);
        convCanvas.Print("output_conv.ps[");
        convCanvas.SetGrid();
        h2dConv->Draw();
        line->Draw("SAME");

        for (const auto &r : results) {
            for (size_t i{0}; i < r.predictedPoints.x.size(); ++i) {
                std::cout << r.predictedPoints.x.at(i) << " " <<  r.predictedPoints.y.at(i) << std::endl;
                TMarker marker(r.predictedPoints.x.at(i), r.predictedPoints.y.at(i), 21);
                marker.SetMarkerSize(1.5);
                marker.SetMarkerColor(r.color);
                marker.DrawClone("SAME");
                drawLabels(r.predictedPoints);
            }
            std::stringstream ss;
            ss.precision(2);
            ss << " " << "#color[" << r.color << "]{" << r.stdAbs << "}";
            h2dConvTitle.append(ss.str());
        }
        h2dConvTitle.append("; AGP-K, %; Chem, %");
        h2dConv->SetTitle(h2dConvTitle.c_str());
        convCanvas.Print("output_conv.ps");
        convCanvas.Print("output_conv.ps]");
        convCanvas.Close();


    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

DataPoints getPredictedPoints(
    const std::map<std::string, SampleData>& data,
    const std::unique_ptr<TF1>& fitFunction,
    bool useAValue)
{
    DataPoints points;

    // Calculate predicted values
    for (const auto& [sampleName, sampleData] : data) {
        for (const auto& resultSet : sampleData.fitResults) {
            std::optional<double> chemValue = useAValue ? sampleData.chem.a : sampleData.chem.w;
            if (!chemValue.has_value()) continue;

            const int numParams = fitFunction->GetNpar();
            double predicted = 0.0;
            for (int p = 0; p < numParams - 1; ++p) {
                predicted += fitFunction->GetParameter(p) * resultSet[p].value;
            }
            predicted += fitFunction->GetParameter(numParams - 1);

            points.labels.push_back(sampleName);
            points.x.push_back(predicted);
            points.y.push_back(chemValue.value());
            points.xError.push_back(0.1);
            points.yError.push_back(0.5);
        }
    }
    return points;
}
