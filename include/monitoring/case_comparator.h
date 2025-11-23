#ifndef MPM_CASE_COMPARATOR_H_
#define MPM_CASE_COMPARATOR_H_

#include <string>
#include <vector>
#include <map>
#include <memory>
#include <functional>

#include "Eigen/Dense"
#include "data_monitor.h"

namespace mpm {

//! Comparison metric types
enum class ComparisonMetric {
  STRESS_DIFFERENCE,      //!< Stress difference
  STRAIN_DIFFERENCE,      //!< Strain difference
  DISPLACEMENT_DIFFERENCE,//!< Displacement difference
  VELOCITY_DIFFERENCE,    //!< Velocity difference
  CORRELATION_COEFFICIENT,//!< Correlation coefficient
  RMS_ERROR,              //!< Root mean square error
  MAX_ABSOLUTE_ERROR,     //!< Maximum absolute error
  RELATIVE_ERROR          //!< Relative error
};

//! Case comparison configuration
struct ComparisonConfig {
  std::vector<ComparisonMetric> metrics;
  double tolerance = 1e-6;
  bool normalize_data = true;
  bool interpolate_missing = true;
  std::string output_format = "json";  // json, csv, html
  bool generate_plots = true;
  std::string plot_format = "png";
};

//! Case data structure
struct CaseData {
  std::string case_name;
  std::string case_path;
  std::map<std::string, std::shared_ptr<TimeSeriesData>> time_series;
  std::map<std::string, double> metadata;
  std::chrono::system_clock::time_point created_at;
  std::string description;
};

//! Comparison result
struct ComparisonResult {
  std::string case1_name;
  std::string case2_name;
  ComparisonMetric metric;
  double value;
  std::string interpretation;  // "Good", "Fair", "Poor", etc.
  std::map<std::string, double> detailed_results;
  std::vector<double> error_distribution;
};

//! Batch comparison results
struct BatchComparisonResults {
  std::vector<ComparisonResult> individual_comparisons;
  std::map<std::string, double> summary_statistics;
  std::string report_html;
  std::chrono::system_clock::time_point generated_at;
};

//! Case comparator for multi-case analysis
//! \brief Compares multiple simulation cases and generates reports
class CaseComparator {
 public:
  //! Constructor
  explicit CaseComparator(const ComparisonConfig& config);
  
  //! Load case data
  bool load_case(const std::string& case_name, const std::string& case_path);
  
  //! Load case from monitor data
  bool load_case_from_monitor(const std::string& case_name, 
                               std::shared_ptr<DataMonitor<2>> monitor);
  bool load_case_from_monitor(const std::string& case_name,
                               std::shared_ptr<DataMonitor<3>> monitor);
  
  //! Compare two cases
  ComparisonResult compare_cases(const std::string& case1_name,
                                const std::string& case2_name,
                                ComparisonMetric metric);
  
  //! Batch compare multiple cases
  BatchComparisonResults batch_compare(
      const std::vector<std::string>& case_names,
      const std::vector<ComparisonMetric>& metrics);
  
  //! Generate comparison report
  std::string generate_report(const BatchComparisonResults& results,
                               const std::string& report_type = "html");
  
  //! Export comparison data
  bool export_comparison(const ComparisonResult& result,
                        const std::string& filename);
  
  //! Get loaded cases
  std::vector<std::string> get_loaded_cases() const;
  
  //! Get case data
  std::shared_ptr<CaseData> get_case_data(const std::string& case_name) const;
  
  //! Clear all cases
  void clear_cases();
  
  //! Set progress callback
  void set_progress_callback(std::function<void(double)> callback) {
    progress_callback_ = callback;
  }
  
 private:
  //! Compute metric between two datasets
  double compute_metric(const TimeSeriesData& data1,
                       const TimeSeriesData& data2,
                       ComparisonMetric metric);
  
  //! Compute stress difference
  double compute_stress_difference(const TimeSeriesData& data1,
                                  const TimeSeriesData& data2);
  
  //! Compute strain difference
  double compute_strain_difference(const TimeSeriesData& data1,
                                  const TimeSeriesData& data2);
  
  //! Compute displacement difference
  double compute_displacement_difference(const TimeSeriesData& data1,
                                        const TimeSeriesData& data2);
  
  //! Compute correlation coefficient
  double compute_correlation(const TimeSeriesData& data1,
                            const TimeSeriesData& data2);
  
  //! Compute RMS error
  double compute_rms_error(const TimeSeriesData& data1,
                          const TimeSeriesData& data2);
  
  //! Compute maximum absolute error
  double compute_max_absolute_error(const TimeSeriesData& data1,
                                     const TimeSeriesData& data2);
  
  //! Compute relative error
  double compute_relative_error(const TimeSeriesData& data1,
                               const TimeSeriesData& data2);
  
  //! Interpolate data to common time points
  void interpolate_to_common_timeline(TimeSeriesData& data1,
                                     TimeSeriesData& data2);
  
  //! Normalize data
  void normalize_data(TimeSeriesData& data);
  
  //! Generate HTML report
  std::string generate_html_report(const BatchComparisonResults& results);
  
  //! Generate CSV report
  std::string generate_csv_report(const BatchComparisonResults& results);
  
  //! Generate JSON report
  std::string generate_json_report(const BatchComparisonResults& results);
  
  //! Interpret metric value
  std::string interpret_metric(ComparisonMetric metric, double value);
  
  //! Configuration
  ComparisonConfig config_;
  
  //! Loaded cases
  std::map<std::string, std::shared_ptr<CaseData>> cases_;
  mutable std::mutex cases_mutex_;
  
  //! Progress callback
  std::function<void(double)> progress_callback_;
  
  //! Logger
  std::shared_ptr<spdlog::logger> console_;
};

}  // namespace mpm

#endif  // MPM_CASE_COMPARATOR_H_