#ifndef MPM_DATA_MONITOR_H_
#define MPM_DATA_MONITOR_H_

#include <memory>
#include <string>
#include <vector>
#include <map>
#include <mutex>
#include <thread>
#include <atomic>
#include <chrono>

#include "Eigen/Dense"
#include "spdlog/spdlog.h"

namespace mpm {

//! Data types for monitoring
enum class MonitorDataType {
  STRESS,      //!< Stress tensor
  STRAIN,      //!< Strain tensor
  VELOCITY,    //!< Velocity vector
  DISPLACEMENT,//!< Displacement vector
  PRESSURE,    //!< Pressure
  MASS,        //!< Mass
  VOLUME,      //!< Volume
  STRAIN_RATE, //!< Strain rate
  TEMPERATURE  //!< Temperature
};

//! Monitoring configuration
struct MonitorConfig {
  bool enabled = true;
  std::chrono::milliseconds update_interval{100};  //!< Update interval
  std::string output_directory = "monitoring_data";
  bool real_time_output = true;
  bool batch_mode = false;
  std::vector<MonitorDataType> data_types;
  double stress_threshold = 0.0;  //!< Stress threshold for alerts
  double strain_threshold = 0.0;  //!< Strain threshold for alerts
};

//! Real-time data point
struct DataPoint {
  size_t step;
  double time;
  MonitorDataType type;
  Eigen::VectorXd data;
  size_t particle_id;
  std::chrono::system_clock::time_point timestamp;
};

//! Time series data for a specific metric
struct TimeSeriesData {
  std::string name;
  MonitorDataType type;
  std::vector<double> timestamps;
  std::vector<Eigen::VectorXd> values;
  std::vector<size_t> step_numbers;
  
  void add_point(double time, const Eigen::VectorXd& value, size_t step) {
    timestamps.push_back(time);
    values.push_back(value);
    step_numbers.push_back(step);
  }
  
  size_t size() const { return timestamps.size(); }
  
  void clear() {
    timestamps.clear();
    values.clear();
    step_numbers.clear();
  }
};

//! Data monitor class for real-time monitoring
//! \brief Monitors and extracts mechanical data during simulation
//! \details Provides real-time extraction of stress, strain, velocity, etc.
template <unsigned Tdim>
class DataMonitor {
 public:
  //! Constructor
  explicit DataMonitor(const MonitorConfig& config);
  
  //! Destructor
  ~DataMonitor();
  
  //! Initialize monitoring
  void initialize(const std::string& case_name);
  
  //! Start monitoring
  void start_monitoring();
  
  //! Stop monitoring
  void stop_monitoring();
  
  //! Update monitoring data at current step
  void update_step(size_t step, double time);
  
  //! Add particle data
  void add_particle_data(size_t particle_id, MonitorDataType type, 
                        const Eigen::VectorXd& data);
  
  //! Add nodal data
  void add_nodal_data(size_t node_id, MonitorDataType type,
                     const Eigen::VectorXd& data);
  
  //! Get time series data
  std::shared_ptr<TimeSeriesData> get_time_series(
      const std::string& name, MonitorDataType type) const;
  
  //! Get current step data
  std::vector<DataPoint> get_current_step_data() const;
  
  //! Check if monitoring is active
  bool is_monitoring() const { return is_monitoring_.load(); }
  
  //! Get configuration
  const MonitorConfig& config() const { return config_; }
  
  //! Export data to file
  void export_data(const std::string& filename) const;
  
  //! Get statistics
  std::map<std::string, double> get_statistics() const;
  
  //! Set alert callback
  void set_alert_callback(std::function<void(const std::string&)> callback) {
    alert_callback_ = callback;
  }
  
 private:
  //! Monitoring thread function
  void monitoring_loop();
  
  //! Process data batch
  void process_data_batch();
  
  //! Check thresholds and trigger alerts
  void check_thresholds(const DataPoint& point);
  
  //! Write real-time data
  void write_realtime_data();
  
  //! Configuration
  MonitorConfig config_;
  
  //! Monitoring state
  std::atomic<bool> is_monitoring_{false};
  std::atomic<bool> should_stop_{false};
  
  //! Thread for monitoring
  std::unique_ptr<std::thread> monitoring_thread_;
  
  //! Data storage
  mutable std::mutex data_mutex_;
  std::map<std::string, std::shared_ptr<TimeSeriesData>> time_series_data_;
  std::vector<DataPoint> current_step_data_;
  std::vector<DataPoint> data_buffer_;
  
  //! Case information
  std::string case_name_;
  size_t current_step_ = 0;
  double current_time_ = 0.0;
  
  //! Logger
  std::shared_ptr<spdlog::logger> console_;
  
  //! Alert callback
  std::function<void(const std::string&)> alert_callback_;
  
  //! Statistics
  mutable std::mutex stats_mutex_;
  std::map<std::string, double> statistics_;
};

}  // namespace mpm

#endif  // MPM_DATA_MONITOR_H_