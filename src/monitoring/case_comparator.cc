#include "monitoring/case_comparator.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <nlohmann/json.hpp>

namespace mpm {

// Constructor
CaseComparator::CaseComparator(const ComparisonConfig& config) : config_(config) {
  console_ = spdlog::get("mpm_comparator")
             ? spdlog::get("mpm_comparator")
             : spdlog::stdout_color_mt("mpm_comparator");
}

// Load case data
bool CaseComparator::load_case(const std::string& case_name, 
                              const std::string& case_path) {
  try {
    auto case_data = std::make_shared<CaseData>();
    case_data->case_name = case_name;
    case_data->case_path = case_path;
    case_data->created_at = std::chrono::system_clock::now();
    
    // Load JSON data file
    std::ifstream file(case_path);
    if (!file.is_open()) {
      console_->error("Failed to open case file: {}", case_path);
      return false;
    }
    
    nlohmann::json json_data;
    file >> json_data;
    file.close();
    
    // Parse time series data
    if (json_data.contains("time_series")) {
      for (auto& [key, series_json] : json_data["time_series"].items()) {
        auto series = std::make_shared<TimeSeriesData>();
        series->name = series_json["name"];
        series->type = static_cast<MonitorDataType>(series_json["type"]);
        
        // Parse timestamps
        for (const auto& ts : series_json["timestamps"]) {
          series->timestamps.push_back(ts);
        }
        
        // Parse step numbers
        for (const auto& step : series_json["step_numbers"]) {
          series->step_numbers.push_back(step);
        }
        
        // Parse values
        for (const auto& value_json : series_json["values"]) {
          Eigen::VectorXd value(value_json.size());
          for (size_t i = 0; i < value_json.size(); ++i) {
            value(i) = value_json[i];
          }
          series->values.push_back(value);
        }
        
        case_data->time_series[key] = series;
      }
    }
    
    // Parse metadata
    if (json_data.contains("metadata")) {
      for (auto& [key, value] : json_data["metadata"].items()) {
        case_data->metadata[key] = value;
      }
    }
    
    {
      std::lock_guard<std::mutex> lock(cases_mutex_);
      cases_[case_name] = case_data;
    }
    
    console_->info("Loaded case: {} from {}", case_name, case_path);
    return true;
    
  } catch (const std::exception& e) {
    console_->error("Error loading case {}: {}", case_name, e.what());
    return false;
  }
}

// Load case from monitor data
bool CaseComparator::load_case_from_monitor(const std::string& case_name,
                                           std::shared_ptr<DataMonitor<2>> monitor) {
  auto case_data = std::make_shared<CaseData>();
  case_data->case_name = case_name;
  case_data->created_at = std::chrono::system_clock::now();
  
  // Export monitor data to temporary file and load it
  std::string temp_file = "/tmp/" + case_name + "_monitor_data.json";
  monitor->export_data(temp_file);
  
  return load_case(case_name, temp_file);
}

bool CaseComparator::load_case_from_monitor(const std::string& case_name,
                                           std::shared_ptr<DataMonitor<3>> monitor) {
  auto case_data = std::make_shared<CaseData>();
  case_data->case_name = case_name;
  case_data->created_at = std::chrono::system_clock::now();
  
  // Export monitor data to temporary file and load it
  std::string temp_file = "/tmp/" + case_name + "_monitor_data.json";
  monitor->export_data(temp_file);
  
  return load_case(case_name, temp_file);
}

// Compare two cases
ComparisonResult CaseComparator::compare_cases(const std::string& case1_name,
                                              const std::string& case2_name,
                                              ComparisonMetric metric) {
  ComparisonResult result;
  result.case1_name = case1_name;
  result.case2_name = case2_name;
  result.metric = metric;
  
  // Get case data
  auto case1 = get_case_data(case1_name);
  auto case2 = get_case_data(case2_name);
  
  if (!case1 || !case2) {
    result.value = std::numeric_limits<double>::quiet_NaN();
    result.interpretation = "Error: Case not found";
    return result;
  }
  
  // Find common time series to compare
  std::vector<std::string> common_series;
  for (const auto& [key, series] : case1->time_series) {
    if (case2->time_series.find(key) != case2->time_series.end()) {
      common_series.push_back(key);
    }
  }
  
  if (common_series.empty()) {
    result.value = std::numeric_limits<double>::quiet_NaN();
    result.interpretation = "Error: No common data series found";
    return result;
  }
  
  // Compare each common series
  std::vector<double> metric_values;
  for (const auto& series_name : common_series) {
    auto& data1 = case1->time_series[series_name];
    auto& data2 = case2->time_series[series_name];
    
    // Interpolate to common timeline if needed
    if (config_.interpolate_missing) {
      interpolate_to_common_timeline(*data1, *data2);
    }
    
    // Normalize if needed
    if (config_.normalize_data) {
      normalize_data(*data1);
      normalize_data(*data2);
    }
    
    double metric_value = compute_metric(*data1, *data2, metric);
    metric_values.push_back(metric_value);
    result.detailed_results[series_name] = metric_value;
    result.error_distribution.push_back(metric_value);
  }
  
  // Compute overall metric value (average)
  result.value = std::accumulate(metric_values.begin(), metric_values.end(), 0.0) 
                  / metric_values.size();
  
  // Interpret result
  result.interpretation = interpret_metric(metric, result.value);
  
  console_->info("Compared {} vs {} using {}: {:.6e} ({})",
                  case1_name, case2_name, 
                  static_cast<int>(metric), result.value, result.interpretation);
  
  return result;
}

// Batch compare multiple cases
BatchComparisonResults CaseComparator::batch_compare(
    const std::vector<std::string>& case_names,
    const std::vector<ComparisonMetric>& metrics) {
  
  BatchComparisonResults results;
  
  size_t total_comparisons = case_names.size() * (case_names.size() - 1) / 2 * metrics.size();
  size_t completed_comparisons = 0;
  
  // Compare all pairs of cases
  for (size_t i = 0; i < case_names.size(); ++i) {
    for (size_t j = i + 1; j < case_names.size(); ++j) {
      for (const auto& metric : metrics) {
        auto result = compare_cases(case_names[i], case_names[j], metric);
        results.individual_comparisons.push_back(result);
        
        completed_comparisons++;
        
        // Update progress
        if (progress_callback_) {
          double progress = static_cast<double>(completed_comparisons) / total_comparisons;
          progress_callback_(progress);
        }
      }
    }
  }
  
  // Compute summary statistics
  for (const auto& metric : metrics) {
    std::vector<double> values;
    for (const auto& result : results.individual_comparisons) {
      if (result.metric == metric && std::isfinite(result.value)) {
        values.push_back(result.value);
      }
    }
    
    if (!values.empty()) {
      std::string metric_name = "avg_" + std::to_string(static_cast<int>(metric));
      results.summary_statistics[metric_name] = 
          std::accumulate(values.begin(), values.end(), 0.0) / values.size();
      
      // Find min and max
      auto [min_it, max_it] = std::minmax_element(values.begin(), values.end());
      results.summary_statistics["min_" + std::to_string(static_cast<int>(metric))] = *min_it;
      results.summary_statistics["max_" + std::to_string(static_cast<int>(metric))] = *max_it;
    }
  }
  
  results.generated_at = std::chrono::system_clock::now();
  results.report_html = generate_html_report(results);
  
  return results;
}

// Compute metric between two datasets
double CaseComparator::compute_metric(const TimeSeriesData& data1,
                                    const TimeSeriesData& data2,
                                    ComparisonMetric metric) {
  switch (metric) {
    case ComparisonMetric::STRESS_DIFFERENCE:
      return compute_stress_difference(data1, data2);
    case ComparisonMetric::STRAIN_DIFFERENCE:
      return compute_strain_difference(data1, data2);
    case ComparisonMetric::DISPLACEMENT_DIFFERENCE:
      return compute_displacement_difference(data1, data2);
    case ComparisonMetric::CORRELATION_COEFFICIENT:
      return compute_correlation(data1, data2);
    case ComparisonMetric::RMS_ERROR:
      return compute_rms_error(data1, data2);
    case ComparisonMetric::MAX_ABSOLUTE_ERROR:
      return compute_max_absolute_error(data1, data2);
    case ComparisonMetric::RELATIVE_ERROR:
      return compute_relative_error(data1, data2);
    default:
      return std::numeric_limits<double>::quiet_NaN();
  }
}

// Compute correlation coefficient
double CaseComparator::compute_correlation(const TimeSeriesData& data1,
                                          const TimeSeriesData& data2) {
  if (data1.size() != data2.size() || data1.size() == 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  
  // For simplicity, we'll compute correlation on the first component
  std::vector<double> values1, values2;
  for (size_t i = 0; i < data1.size(); ++i) {
    values1.push_back(data1.values[i](0));  // First component
    values2.push_back(data2.values[i](0));  // First component
  }
  
  // Compute means
  double mean1 = std::accumulate(values1.begin(), values1.end(), 0.0) / values1.size();
  double mean2 = std::accumulate(values2.begin(), values2.end(), 0.0) / values2.size();
  
  // Compute correlation
  double numerator = 0.0, denom1 = 0.0, denom2 = 0.0;
  for (size_t i = 0; i < values1.size(); ++i) {
    double diff1 = values1[i] - mean1;
    double diff2 = values2[i] - mean2;
    numerator += diff1 * diff2;
    denom1 += diff1 * diff1;
    denom2 += diff2 * diff2;
  }
  
  if (denom1 == 0.0 || denom2 == 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  
  return numerator / std::sqrt(denom1 * denom2);
}

// Compute RMS error
double CaseComparator::compute_rms_error(const TimeSeriesData& data1,
                                        const TimeSeriesData& data2) {
  if (data1.size() != data2.size() || data1.size() == 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  
  double sum_squared_error = 0.0;
  size_t total_points = 0;
  
  for (size_t i = 0; i < data1.size(); ++i) {
    const auto& vec1 = data1.values[i];
    const auto& vec2 = data2.values[i];
    
    if (vec1.size() == vec2.size()) {
      for (int j = 0; j < vec1.size(); ++j) {
        double diff = vec1(j) - vec2(j);
        sum_squared_error += diff * diff;
        total_points++;
      }
    }
  }
  
  if (total_points == 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  
  return std::sqrt(sum_squared_error / total_points);
}

// Interpret metric value
std::string CaseComparator::interpret_metric(ComparisonMetric metric, double value) {
  switch (metric) {
    case ComparisonMetric::CORRELATION_COEFFICIENT:
      if (value > 0.95) return "Excellent";
      else if (value > 0.80) return "Good";
      else if (value > 0.60) return "Fair";
      else return "Poor";
      
    case ComparisonMetric::RMS_ERROR:
    case ComparisonMetric::MAX_ABSOLUTE_ERROR:
    case ComparisonMetric::RELATIVE_ERROR:
      if (value < 0.01) return "Excellent";
      else if (value < 0.05) return "Good";
      else if (value < 0.10) return "Fair";
      else return "Poor";
      
    default:
      return "Unknown";
  }
}

// Get loaded cases
std::vector<std::string> CaseComparator::get_loaded_cases() const {
  std::lock_guard<std::mutex> lock(cases_mutex_);
  std::vector<std::string> case_names;
  for (const auto& [name, data] : cases_) {
    case_names.push_back(name);
  }
  return case_names;
}

// Get case data
std::shared_ptr<CaseData> CaseComparator::get_case_data(const std::string& case_name) const {
  std::lock_guard<std::mutex> lock(cases_mutex_);
  auto it = cases_.find(case_name);
  if (it != cases_.end()) {
    return it->second;
  }
  return nullptr;
}

// Clear all cases
void CaseComparator::clear_cases() {
  std::lock_guard<std::mutex> lock(cases_mutex_);
  cases_.clear();
  console_->info("All cases cleared");
}

// Generate comparison report
std::string CaseComparator::generate_report(const BatchComparisonResults& results,
                                           const std::string& report_type) {
  if (report_type == "html") {
    return generate_html_report(results);
  } else if (report_type == "csv") {
    return generate_csv_report(results);
  } else if (report_type == "json") {
    return generate_json_report(results);
  } else {
    return "Unsupported report type: " + report_type;
  }
}

// Generate HTML report
std::string CaseComparator::generate_html_report(const BatchComparisonResults& results) {
  std::stringstream html;
  
  html << "<!DOCTYPE html>\n";
  html << "<html>\n<head>\n";
  html << "<title>MPM Case Comparison Report</title>\n";
  html << "<style>\n";
  html << "body { font-family: Arial, sans-serif; margin: 20px; }\n";
  html << "table { border-collapse: collapse; width: 100%; margin: 20px 0; }\n";
  html << "th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }\n";
  html << "th { background-color: #f2f2f2; }\n";
  html << ".excellent { background-color: #d4edda; }\n";
  html << ".good { background-color: #d1ecf1; }\n";
  html << ".fair { background-color: #fff3cd; }\n";
  html << ".poor { background-color: #f8d7da; }\n";
  html << "</style>\n</head>\n<body>\n";
  
  html << "<h1>MPM Case Comparison Report</h1>\n";
  html << "<p>Generated: " << std::chrono::system_clock::to_time_t(results.generated_at) << "</p>\n";
  
  // Summary statistics
  html << "<h2>Summary Statistics</h2>\n";
  html << "<table>\n";
  html << "<tr><th>Metric</th><th>Value</th></tr>\n";
  for (const auto& [key, value] : results.summary_statistics) {
    html << "<tr><td>" << key << "</td><td>" << value << "</td></tr>\n";
  }
  html << "</table>\n";
  
  // Individual comparisons
  html << "<h2>Individual Comparisons</h2>\n";
  html << "<table>\n";
  html << "<tr><th>Case 1</th><th>Case 2</th><th>Metric</th><th>Value</th><th>Interpretation</th></tr>\n";
  
  for (const auto& result : results.individual_comparisons) {
    std::string css_class = result.interpretation;
    std::transform(css_class.begin(), css_class.end(), css_class.begin(), ::tolower);
    
    html << "<tr class=\"" << css_class << "\">";
    html << "<td>" << result.case1_name << "</td>";
    html << "<td>" << result.case2_name << "</td>";
    html << "<td>" << static_cast<int>(result.metric) << "</td>";
    html << "<td>" << result.value << "</td>";
    html << "<td>" << result.interpretation << "</td>";
    html << "</tr>\n";
  }
  
  html << "</table>\n";
  html << "</body>\n</html>";
  
  return html.str();
}

}  // namespace mpm