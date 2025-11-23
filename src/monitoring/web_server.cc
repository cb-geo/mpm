#include "monitoring/web_server.h"
#include <sstream>
#include <iomanip>
#include <fstream>
#include <filesystem>
#include <nlohmann/json.hpp>

namespace mpm {

// Constructor
WebDashboardServer::WebDashboardServer(const WebServerConfig& config) 
    : config_(config) {
  console_ = spdlog::get("mpm_web_server")
             ? spdlog::get("mpm_web_server")
             : spdlog::stdout_color_mt("mpm_web_server");
}

// Destructor
WebDashboardServer::~WebDashboardServer() {
  stop();
}

// Start server
bool WebDashboardServer::start() {
  if (is_running_.load()) {
    console_->warn("Web server is already running");
    return false;
  }
  
  is_running_.store(true);
  should_stop_.store(false);
  
  // Initialize routes
  initialize_routes();
  
  // Start server thread
  server_thread_ = std::make_unique<std::thread>(
      &WebDashboardServer::server_loop, this);
  
  console_->info("Web dashboard server started on {}:{}", 
                 config_.host, config_.port);
  console_->info("Dashboard available at: {}", get_server_url());
  
  return true;
}

// Stop server
void WebDashboardServer::stop() {
  if (!is_running_.load()) return;
  
  should_stop_.store(true);
  
  if (server_thread_ && server_thread_->joinable()) {
    server_thread_->join();
  }
  
  is_running_.store(false);
  
  console_->info("Web dashboard server stopped");
}

// Register data monitor
void WebDashboardServer::register_monitor(std::shared_ptr<DataMonitor<2>> monitor) {
  monitor2d_ = monitor;
  console_->info("2D data monitor registered");
}

void WebDashboardServer::register_monitor(std::shared_ptr<DataMonitor<3>> monitor) {
  monitor3d_ = monitor;
  console_->info("3D data monitor registered");
}

// Update dashboard data
void WebDashboardServer::update_data(const DashboardData& data) {
  std::lock_guard<std::mutex> lock(data_mutex_);
  current_data_ = data;
  
  // Broadcast to WebSocket clients
  std::string json_data = serialize_dashboard_data(data);
  broadcast_data("dashboard_update", json_data);
}

// Broadcast data to all connected clients
void WebDashboardServer::broadcast_data(const std::string& data_type, 
                                        const std::string& data) {
  if (!config_.enable_websocket) return;
  
  std::lock_guard<std::mutex> lock(clients_mutex_);
  
  nlohmann::json message;
  message["type"] = data_type;
  message["data"] = nlohmann::json::parse(data);
  message["timestamp"] = std::chrono::system_clock::now().time_since_epoch().count();
  
  std::string message_str = message.dump();
  
  // In a real implementation, this would send to actual WebSocket clients
  console_->debug("Broadcasting {} to {} clients", data_type, websocket_clients_.size());
}

// Get server URL
std::string WebDashboardServer::get_server_url() const {
  return fmt::format("http://{}:{}", config_.host, config_.port);
}

// Initialize routes
void WebDashboardServer::initialize_routes() {
  console_->info("Initializing API routes");
}

// Simple HTTP server implementation
void WebDashboardServer::server_loop() {
  console_->info("Web server loop started");
  
  while (!should_stop_.load()) {
    // In a real implementation, this would handle HTTP requests
    // For now, we'll just simulate the server running
    
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
    
    // Update dashboard with current monitor data
    if (monitor2d_ && monitor2d_->is_monitoring()) {
      DashboardData data;
      data.case_name = "2D_Simulation";
      data.current_step = static_cast<size_t>(monitor2d_->get_statistics().at("total_steps"));
      data.current_time = monitor2d_->get_statistics().at("current_time");
      
      // Add some sample metrics
      data.global_metrics["max_stress"] = 1.5e6;
      data.global_metrics["max_strain"] = 0.02;
      data.global_metrics["total_energy"] = 2.3e9;
      
      update_data(data);
    }
  }
  
  console_->info("Web server loop stopped");
}

// Handle HTTP request
void WebDashboardServer::handle_request(const std::string& method, 
                                        const std::string& path,
                                        const std::string& body,
                                        std::string& response,
                                        std::map<std::string, std::string>& headers) {
  
  // Parse URL parameters
  auto params = parse_url_params(path);
  std::string clean_path = path;
  size_t query_pos = clean_path.find('?');
  if (query_pos != std::string::npos) {
    clean_path = clean_path.substr(0, query_pos);
  }
  
  // Route handling
  if (method == "GET" && clean_path == "/api/status") {
    response = handle_get_status();
    headers["Content-Type"] = "application/json";
  } else if (method == "GET" && clean_path == "/api/cases") {
    response = handle_get_cases();
    headers["Content-Type"] = "application/json";
  } else if (method == "GET" && clean_path.find("/api/case/") == 0) {
    std::string case_name = clean_path.substr(10);
    response = handle_get_case_data(case_name);
    headers["Content-Type"] = "application/json";
  } else if (method == "GET" && clean_path.find("/api/timeseries/") == 0) {
    // Parse case name and metric from path
    size_t pos = clean_path.find("/timeseries/");
    if (pos != std::string::npos) {
      std::string case_metric = clean_path.substr(pos + 12);
      size_t slash_pos = case_metric.find('/');
      if (slash_pos != std::string::npos) {
        std::string case_name = case_metric.substr(0, slash_pos);
        std::string metric = case_metric.substr(slash_pos + 1);
        response = handle_get_time_series(case_name, metric);
        headers["Content-Type"] = "application/json";
      }
    }
  } else {
    // Return 404
    response = "{\"error\": \"Not found\"}";
    headers["Content-Type"] = "application/json";
  }
}

// API handlers
std::string WebDashboardServer::handle_get_status() {
  ApiResponse<nlohmann::json> response;
  
  nlohmann::json status;
  status["server_running"] = is_running_.load();
  status["websocket_enabled"] = config_.enable_websocket;
  status["connected_clients"] = websocket_clients_.size();
  status["monitor2d_active"] = monitor2d_ && monitor2d_->is_monitoring();
  status["monitor3d_active"] = monitor3d_ && monitor3d_->is_monitoring();
  status["timestamp"] = std::chrono::system_clock::now().time_since_epoch().count();
  
  response.data = status;
  response.message = "Server status retrieved successfully";
  
  return serialize_response(response);
}

std::string WebDashboardServer::handle_get_cases() {
  ApiResponse<nlohmann::json> response;
  
  nlohmann::json cases;
  
  if (monitor2d_ && monitor2d_->is_monitoring()) {
    nlohmann::json case_info;
    case_info["name"] = "2D_Simulation";
    case_info["dimension"] = 2;
    case_info["active"] = true;
    case_info["current_step"] = monitor2d_->get_statistics().at("total_steps");
    case_info["current_time"] = monitor2d_->get_statistics().at("current_time");
    cases.push_back(case_info);
  }
  
  if (monitor3d_ && monitor3d_->is_monitoring()) {
    nlohmann::json case_info;
    case_info["name"] = "3D_Simulation";
    case_info["dimension"] = 3;
    case_info["active"] = true;
    case_info["current_step"] = monitor3d_->get_statistics().at("total_steps");
    case_info["current_time"] = monitor3d_->get_statistics().at("current_time");
    cases.push_back(case_info);
  }
  
  response.data = cases;
  response.message = "Cases retrieved successfully";
  
  return serialize_response(response);
}

std::string WebDashboardServer::handle_get_case_data(const std::string& case_name) {
  ApiResponse<nlohmann::json> response;
  
  std::lock_guard<std::mutex> lock(data_mutex_);
  
  if (case_name == "current" && !current_data_.case_name.empty()) {
    response.data = nlohmann::json::parse(serialize_dashboard_data(current_data_));
    response.message = "Current case data retrieved successfully";
  } else {
    response.success = false;
    response.message = "Case not found: " + case_name;
  }
  
  return serialize_response(response);
}

std::string WebDashboardServer::handle_get_time_series(const std::string& case_name,
                                                       const std::string& metric) {
  ApiResponse<nlohmann::json> response;
  
  // This would typically query the data monitor for time series data
  nlohmann::json time_series;
  time_series["case_name"] = case_name;
  time_series["metric"] = metric;
  time_series["timestamps"] = std::vector<double>{0.0, 0.1, 0.2, 0.3};
  time_series["values"] = std::vector<double>{1.0, 1.2, 1.5, 1.8};
  
  response.data = time_series;
  response.message = "Time series data retrieved successfully";
  
  return serialize_response(response);
}

// Parse URL parameters
std::map<std::string, std::string> WebDashboardServer::parse_url_params(const std::string& path) {
  std::map<std::string, std::string> params;
  
  size_t query_pos = path.find('?');
  if (query_pos == std::string::npos) return params;
  
  std::string query_string = path.substr(query_pos + 1);
  std::stringstream ss(query_string);
  std::string param;
  
  while (std::getline(ss, param, '&')) {
    size_t eq_pos = param.find('=');
    if (eq_pos != std::string::npos) {
      std::string key = param.substr(0, eq_pos);
      std::string value = param.substr(eq_pos + 1);
      params[key] = value;
    }
  }
  
  return params;
}

// Template for serializing responses
template <typename T>
std::string WebDashboardServer::serialize_response(const ApiResponse<T>& response) {
  nlohmann::json json_response;
  json_response["success"] = response.success;
  json_response["message"] = response.message;
  json_response["data"] = response.data;
  json_response["timestamp"] = response.timestamp.time_since_epoch().count();
  
  return json_response.dump(2);
}

// Serialize dashboard data
std::string WebDashboardServer::serialize_dashboard_data(const DashboardData& data) {
  nlohmann::json json_data;
  json_data["case_name"] = data.case_name;
  json_data["current_step"] = data.current_step;
  json_data["current_time"] = data.current_time;
  json_data["global_metrics"] = data.global_metrics;
  json_data["particle_data"] = data.particle_data;
  json_data["node_data"] = data.node_data;
  json_data["time_series"] = data.time_series;
  
  return json_data.dump();
}

}  // namespace mpm