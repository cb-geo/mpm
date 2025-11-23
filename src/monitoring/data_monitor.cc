#include "monitoring/data_monitor.h"
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>

namespace mpm {

// Constructor
template <unsigned Tdim>
DataMonitor<Tdim>::DataMonitor(const MonitorConfig& config) : config_(config) {
  console_ = spdlog::get("mpm_monitor") 
             ? spdlog::get("mpm_monitor")
             : spdlog::stdout_color_mt("mpm_monitor");
}

// Destructor
template <unsigned Tdim>
DataMonitor<Tdim>::~DataMonitor() {
  stop_monitoring();
}

// Initialize monitoring
template <unsigned Tdim>
void DataMonitor<Tdim>::initialize(const std::string& case_name) {
  case_name_ = case_name;
  
  // Create output directory
  std::string full_path = config_.output_directory + "/" + case_name;
  std::filesystem::create_directories(full_path);
  
  console_->info("Data monitor initialized for case: {}", case_name);
}

// Start monitoring
template <unsigned Tdim>
void DataMonitor<Tdim>::start_monitoring() {
  if (is_monitoring_.load()) {
    console_->warn("Monitoring is already running");
    return;
  }
  
  is_monitoring_.store(true);
  should_stop_.store(false);
  
  monitoring_thread_ = std::make_unique<std::thread>(
      &DataMonitor<Tdim>::monitoring_loop, this);
  
  console_->info("Data monitoring started");
}

// Stop monitoring
template <unsigned Tdim>
void DataMonitor<Tdim>::stop_monitoring() {
  if (!is_monitoring_.load()) return;
  
  should_stop_.store(true);
  
  if (monitoring_thread_ && monitoring_thread_->joinable()) {
    monitoring_thread_->join();
  }
  
  is_monitoring_.store(false);
  
  // Write final data
  if (!data_buffer_.empty()) {
    process_data_batch();
  }
  
  console_->info("Data monitoring stopped");
}

// Update monitoring data at current step
template <unsigned Tdim>
void DataMonitor<Tdim>::update_step(size_t step, double time) {
  std::lock_guard<std::mutex> lock(data_mutex_);
  current_step_ = step;
  current_time_ = time;
  
  // Process current step data
  if (!current_step_data_.empty()) {
    data_buffer_.insert(data_buffer_.end(), 
                         current_step_data_.begin(), 
                         current_step_data_.end());
    current_step_data_.clear();
  }
  
  // Update statistics
  {
    std::lock_guard<std::mutex> stats_lock(stats_mutex_);
    statistics_["total_steps"] = static_cast<double>(step);
    statistics_["current_time"] = time;
    statistics_["data_points"] = static_cast<double>(data_buffer_.size());
  }
}

// Add particle data
template <unsigned Tdim>
void DataMonitor<Tdim>::add_particle_data(size_t particle_id, 
                                          MonitorDataType type,
                                          const Eigen::VectorXd& data) {
  if (!is_monitoring_.load()) return;
  
  DataPoint point;
  point.step = current_step_;
  point.time = current_time_;
  point.type = type;
  point.data = data;
  point.particle_id = particle_id;
  point.timestamp = std::chrono::system_clock::now();
  
  {
    std::lock_guard<std::mutex> lock(data_mutex_);
    current_step_data_.push_back(point);
  }
  
  // Check thresholds
  check_thresholds(point);
}

// Add nodal data
template <unsigned Tdim>
void DataMonitor<Tdim>::add_nodal_data(size_t node_id,
                                       MonitorDataType type,
                                       const Eigen::VectorXd& data) {
  if (!is_monitoring_.load()) return;
  
  DataPoint point;
  point.step = current_step_;
  point.time = current_time_;
  point.type = type;
  point.data = data;
  point.particle_id = node_id;  // Using particle_id field for node_id
  point.timestamp = std::chrono::system_clock::now();
  
  {
    std::lock_guard<std::mutex> lock(data_mutex_);
    current_step_data_.push_back(point);
  }
  
  check_thresholds(point);
}

// Get time series data
template <unsigned Tdim>
std::shared_ptr<TimeSeriesData> DataMonitor<Tdim>::get_time_series(
    const std::string& name, MonitorDataType type) const {
  std::lock_guard<std::mutex> lock(data_mutex_);
  
  auto it = time_series_data_.find(name);
  if (it != time_series_data_.end() && it->second->type == type) {
    return it->second;
  }
  
  return nullptr;
}

// Get current step data
template <unsigned Tdim>
std::vector<DataPoint> DataMonitor<Tdim>::get_current_step_data() const {
  std::lock_guard<std::mutex> lock(data_mutex_);
  return current_step_data_;
}

// Export data to file
template <unsigned Tdim>
void DataMonitor<Tdim>::export_data(const std::string& filename) const {
  std::lock_guard<std::mutex> lock(data_mutex_);
  
  nlohmann::json json_data;
  json_data["case_name"] = case_name_;
  json_data["export_time"] = std::chrono::system_clock::now().time_since_epoch().count();
  
  // Export time series data
  nlohmann::json time_series_json;
  for (const auto& [name, data] : time_series_data_) {
    nlohmann::json series;
    series["name"] = name;
    series["type"] = static_cast<int>(data->type);
    series["timestamps"] = data->timestamps;
    series["step_numbers"] = data->step_numbers;
    
    // Convert Eigen vectors to arrays
    nlohmann::json values_array;
    for (const auto& value : data->values) {
      nlohmann::json vec;
      for (int i = 0; i < value.size(); ++i) {
        vec.push_back(value(i));
      }
      values_array.push_back(vec);
    }
    series["values"] = values_array;
    
    time_series_json[name] = series;
  }
  json_data["time_series"] = time_series_json;
  
  // Write to file
  std::ofstream file(filename);
  if (file.is_open()) {
    file << json_data.dump(2);
    file.close();
    console_->info("Data exported to: {}", filename);
  } else {
    console_->error("Failed to export data to: {}", filename);
  }
}

// Get statistics
template <unsigned Tdim>
std::map<std::string, double> DataMonitor<Tdim>::get_statistics() const {
  std::lock_guard<std::mutex> lock(stats_mutex_);
  return statistics_;
}

// Monitoring thread function
template <unsigned Tdim>
void DataMonitor<Tdim>::monitoring_loop() {
  console_->info("Monitoring thread started");
  
  while (!should_stop_.load()) {
    auto start_time = std::chrono::steady_clock::now();
    
    // Process data batch
    if (!data_buffer_.empty()) {
      process_data_batch();
    }
    
    // Write real-time data if enabled
    if (config_.real_time_output) {
      write_realtime_data();
    }
    
    // Sleep for the remaining time
    auto end_time = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        end_time - start_time);
    
    if (elapsed < config_.update_interval) {
      std::this_thread::sleep_for(config_.update_interval - elapsed);
    }
  }
  
  console_->info("Monitoring thread stopped");
}

// Process data batch
template <unsigned Tdim>
void DataMonitor<Tdim>::process_data_batch() {
  std::vector<DataPoint> batch;
  
  {
    std::lock_guard<std::mutex> lock(data_mutex_);
    batch = std::move(data_buffer_);
    data_buffer_.clear();
  }
  
  // Process each data point
  for (const auto& point : batch) {
    // Create time series key
    std::string key = std::to_string(static_cast<int>(point.type)) + "_" +
                     std::to_string(point.particle_id);
    
    std::lock_guard<std::mutex> lock(data_mutex_);
    
    // Create or get time series
    if (time_series_data_.find(key) == time_series_data_.end()) {
      auto series = std::make_shared<TimeSeriesData>();
      series->name = key;
      series->type = point.type;
      time_series_data_[key] = series;
    }
    
    // Add data point
    time_series_data_[key]->add_point(point.time, point.data, point.step);
  }
  
  // Update statistics
  {
    std::lock_guard<std::mutex> stats_lock(stats_mutex_);
    statistics_["processed_points"] += static_cast<double>(batch.size());
  }
}

// Check thresholds and trigger alerts
template <unsigned Tdim>
void DataMonitor<Tdim>::check_thresholds(const DataPoint& point) {
  bool should_alert = false;
  std::string alert_message;
  
  switch (point.type) {
    case MonitorDataType::STRESS:
      if (config_.stress_threshold > 0) {
        double max_stress = point.data.maxCoeff();
        if (max_stress > config_.stress_threshold) {
          should_alert = true;
          alert_message = fmt::format(
              "High stress alert: {:.3e} at step {} (threshold: {:.3e})",
              max_stress, point.step, config_.stress_threshold);
        }
      }
      break;
      
    case MonitorDataType::STRAIN:
      if (config_.strain_threshold > 0) {
        double max_strain = point.data.maxCoeff();
        if (max_strain > config_.strain_threshold) {
          should_alert = true;
          alert_message = fmt::format(
              "High strain alert: {:.3e} at step {} (threshold: {:.3e})",
              max_strain, point.step, config_.strain_threshold);
        }
      }
      break;
      
    default:
      break;
  }
  
  if (should_alert && alert_callback_) {
    alert_callback_(alert_message);
    console_->warn("Alert triggered: {}", alert_message);
  }
}

// Write real-time data
template <unsigned Tdim>
void DataMonitor<Tdim>::write_realtime_data() {
  // This would typically write to a real-time database or message queue
  // For now, we'll just log the current state
  std::lock_guard<std::mutex> lock(data_mutex_);
  
  console_->debug("Real-time update - Step: {}, Time: {:.3e}, Data points: {}",
                  current_step_, current_time_, data_buffer_.size());
}

// Explicit template instantiation
template class DataMonitor<2>;
template class DataMonitor<3>;

}  // namespace mpm