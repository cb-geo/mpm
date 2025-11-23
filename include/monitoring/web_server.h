#ifndef MPM_WEB_SERVER_H_
#define MPM_WEB_SERVER_H_

#include <string>
#include <vector>
#include <map>
#include <memory>
#include <functional>
#include <thread>
#include <atomic>

#include "data_monitor.h"
#include "spdlog/spdlog.h"

namespace mpm {

//! Web server configuration
struct WebServerConfig {
  int port = 8080;
  std::string host = "0.0.0.0";
  bool enable_cors = true;
  std::string static_files_path = "./web_dashboard";
  int update_interval_ms = 100;
  bool enable_websocket = true;
};

//! API response structure
template <typename T>
struct ApiResponse {
  bool success;
  std::string message;
  T data;
  std::chrono::system_clock::time_point timestamp;
  
  ApiResponse() : success(true), timestamp(std::chrono::system_clock::now()) {}
  ApiResponse(bool s, const std::string& m, const T& d) 
      : success(s), message(m), data(d), timestamp(std::chrono::system_clock::now()) {}
};

//! Dashboard data structure
struct DashboardData {
  std::string case_name;
  size_t current_step;
  double current_time;
  std::map<std::string, double> global_metrics;
  std::vector<std::map<std::string, double>> particle_data;
  std::vector<std::map<std::string, double>> node_data;
  std::map<std::string, std::vector<double>> time_series;
};

//! Web dashboard server
//! \brief Provides REST API and WebSocket for real-time data visualization
class WebDashboardServer {
 public:
  //! Constructor
  explicit WebDashboardServer(const WebServerConfig& config);
  
  //! Destructor
  ~WebDashboardServer();
  
  //! Start server
  bool start();
  
  //! Stop server
  void stop();
  
  //! Register data monitor
  void register_monitor(std::shared_ptr<DataMonitor<2>> monitor2d);
  void register_monitor(std::shared_ptr<DataMonitor<3>> monitor3d);
  
  //! Update dashboard data
  void update_data(const DashboardData& data);
  
  //! Broadcast data to all connected clients
  void broadcast_data(const std::string& data_type, const std::string& data);
  
  //! Check if server is running
  bool is_running() const { return is_running_.load(); }
  
  //! Get server URL
  std::string get_server_url() const;
  
 private:
  //! Initialize routes
  void initialize_routes();
  
  //! API handlers
  std::string handle_get_status();
  std::string handle_get_cases();
  std::string handle_get_case_data(const std::string& case_name);
  std::string handle_get_time_series(const std::string& case_name, 
                                     const std::string& metric);
  std::string handle_get_particle_data(const std::string& case_name, 
                                       size_t step);
  std::string handle_get_statistics(const std::string& case_name);
  
  //! WebSocket handlers
  void handle_websocket_connect(int client_id);
  void handle_websocket_disconnect(int client_id);
  void handle_websocket_message(int client_id, const std::string& message);
  
  //! Data conversion helpers
  template <typename T>
  std::string serialize_response(const ApiResponse<T>& response);
  
  std::string serialize_dashboard_data(const DashboardData& data);
  
  //! Simple HTTP server implementation
  void server_loop();
  
  //! Handle HTTP request
  void handle_request(const std::string& method, const std::string& path,
                     const std::string& body, std::string& response,
                     std::map<std::string, std::string>& headers);
  
  //! Parse URL parameters
  std::map<std::string, std::string> parse_url_params(const std::string& path);
  
  //! Configuration
  WebServerConfig config_;
  
  //! Server state
  std::atomic<bool> is_running_{false};
  std::atomic<bool> should_stop_{false};
  
  //! Server thread
  std::unique_ptr<std::thread> server_thread_;
  
  //! Data monitors
  std::shared_ptr<DataMonitor<2>> monitor2d_;
  std::shared_ptr<DataMonitor<3>> monitor3d_;
  
  //! Current dashboard data
  DashboardData current_data_;
  std::mutex data_mutex_;
  
  //! WebSocket clients
  std::map<int, bool> websocket_clients_;
  std::mutex clients_mutex_;
  
  //! Logger
  std::shared_ptr<spdlog::logger> console_;
  
  //! Simple socket for demo (in production, use proper HTTP library)
  int server_socket_ = -1;
};

}  // namespace mpm

#endif  // MPM_WEB_SERVER_H_