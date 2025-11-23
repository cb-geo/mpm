#include <iostream>
#include <memory>
#include <thread>
#include <chrono>
#include <csignal>
#include <nlohmann/json.hpp>

#include "mpm.h"
#include "monitoring/data_monitor.h"
#include "monitoring/web_server.h"
#include "monitoring/case_comparator.h"

using namespace mpm;

// Global variables for signal handling
std::unique_ptr<WebDashboardServer> web_server;
std::unique_ptr<DataMonitor<3>> monitor;
std::atomic<bool> running{true};

void signal_handler(int signal) {
  std::cout << "Received signal " << signal << ", shutting down..." << std::endl;
  running = false;
  
  if (web_server) {
    web_server->stop();
  }
  if (monitor) {
    monitor->stop_monitoring();
  }
}

// Example usage of the monitoring system
int main(int argc, char* argv[]) {
  try {
    // Set up signal handling
    std::signal(SIGINT, signal_handler);
    std::signal(SIGTERM, signal_handler);
    
    std::cout << "MPM Real-time Monitoring System Demo" << std::endl;
    std::cout << "=====================================" << std::endl;
    
    // Initialize monitoring
    auto monitor_config = std::make_shared<MonitorConfig>();
    monitor_config->update_interval = 0.1;  // Update every 0.1 seconds
    monitor_config->max_history_size = 1000;  // Keep last 1000 data points
    monitor_config->export_interval = 10.0;   // Export every 10 seconds
    monitor_config->export_format = "json";
    
    monitor = std::make_unique<DataMonitor<3>>(monitor_config);
    
    // Add monitoring targets
    monitor->add_particle_monitor("stress", MonitorDataType::STRESS);
    monitor->add_particle_monitor("strain", MonitorDataType::STRAIN);
    monitor->add_particle_monitor("velocity", MonitorDataType::VELOCITY);
    monitor->add_node_monitor("displacement", MonitorDataType::DISPLACEMENT);
    monitor->add_node_monitor("force", MonitorDataType::FORCE);
    
    // Initialize web server
    WebServerConfig server_config;
    server_config.port = 8080;
    server_config.enable_cors = true;
    server_config.enable_websocket = true;
    server_config.websocket_port = 8081;
    
    web_server = std::make_unique<WebDashboardServer>(server_config);
    
    // Register monitor with web server
    web_server->register_monitor("main_simulation", monitor);
    
    // Start monitoring and web server
    monitor->start_monitoring();
    web_server->start();
    
    std::cout << "Web server started at http://localhost:8080" << std::endl;
    std::cout << "WebSocket server started at ws://localhost:8081" << std::endl;
    std::cout << "Press Ctrl+C to stop..." << std::endl;
    
    // Simulate some data generation
    std::cout << "Generating sample data..." << std::endl;
    
    double time = 0.0;
    int step = 0;
    
    while (running && step < 100) {
      // Simulate particle data
      Eigen::Vector3d particle_stress(100.0 + 10.0 * std::sin(time), 
                                       50.0 + 5.0 * std::cos(time), 
                                       25.0 + 2.0 * std::sin(2.0 * time));
      
      Eigen::Vector3d particle_strain(0.1 + 0.01 * std::sin(time),
                                       0.05 + 0.005 * std::cos(time),
                                       0.025 + 0.002 * std::sin(2.0 * time));
      
      Eigen::Vector3d particle_velocity(1.0 + 0.1 * std::sin(time),
                                        0.5 + 0.05 * std::cos(time),
                                        0.25 + 0.02 * std::sin(2.0 * time));
      
      // Simulate node data
      Eigen::Vector3d node_displacement(0.1 * time,
                                         0.05 * time,
                                         0.025 * time);
      
      Eigen::Vector3d node_force(100.0 * std::sin(time),
                                  50.0 * std::cos(time),
                                  25.0 * std::sin(2.0 * time));
      
      // Update monitor data
      monitor->update_data("stress", step, time, particle_stress);
      monitor->update_data("strain", step, time, particle_strain);
      monitor->update_data("velocity", step, time, particle_velocity);
      monitor->update_data("displacement", step, time, node_displacement);
      monitor->update_data("force", step, time, node_force);
      
      // Simulate multiple particles and nodes
      for (int i = 0; i < 10; ++i) {
        Eigen::Vector3d particle_pos(0.1 * i, 0.2 * i, 0.3 * i);
        
        Eigen::Vector3d stress = particle_stress + Eigen::Vector3d(5.0 * i, 2.5 * i, 1.0 * i);
        Eigen::Vector3d strain = particle_strain + Eigen::Vector3d(0.005 * i, 0.0025 * i, 0.001 * i);
        Eigen::Vector3d velocity = particle_velocity + Eigen::Vector3d(0.1 * i, 0.05 * i, 0.02 * i);
        
        monitor->update_particle_data(i, particle_pos, stress, strain, velocity);
      }
      
      for (int i = 0; i < 5; ++i) {
        Eigen::Vector3d node_pos(0.2 * i, 0.3 * i, 0.4 * i);
        
        Eigen::Vector3d displacement = node_displacement + Eigen::Vector3d(0.01 * i, 0.005 * i, 0.002 * i);
        Eigen::Vector3d force = node_force + Eigen::Vector3d(10.0 * i, 5.0 * i, 2.0 * i);
        
        monitor->update_node_data(i, node_pos, displacement, force);
      }
      
      time += 0.1;
      step++;
      
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
    
    std::cout << "Sample data generation completed." << std::endl;
    
    // Example case comparison
    std::cout << "Performing case comparison..." << std::endl;
    
    CaseComparator comparator(ComparisonConfig{});
    
    // Save current monitor data as a case
    monitor->export_data("case1.json");
    comparator.load_case("case1", "case1.json");
    
    // Simulate a slightly different case
    step = 0;
    time = 0.0;
    
    while (running && step < 100) {
      Eigen::Vector3d particle_stress(95.0 + 8.0 * std::sin(time), 
                                       45.0 + 4.0 * std::cos(time), 
                                       20.0 + 1.5 * std::sin(2.0 * time));
      
      monitor->update_data("stress", step, time, particle_stress);
      
      time += 0.1;
      step++;
      
      std::this_thread::sleep_for(std::chrono::milliseconds(50));
    }
    
    monitor->export_data("case2.json");
    comparator.load_case("case2", "case2.json");
    
    // Perform comparison
    auto result = comparator.compare_cases("case1", "case2", 
                                          ComparisonMetric::CORRELATION_COEFFICIENT);
    
    std::cout << "Comparison result: " << result.value 
              << " (" << result.interpretation << ")" << std::endl;
    
    // Generate report
    std::vector<std::string> cases = {"case1", "case2"};
    std::vector<ComparisonMetric> metrics = {
      ComparisonMetric::CORRELATION_COEFFICIENT,
      ComparisonMetric::RMS_ERROR,
      ComparisonMetric::MAX_ABSOLUTE_ERROR
    };
    
    auto batch_results = comparator.batch_compare(cases, metrics);
    std::string html_report = comparator.generate_report(batch_results, "html");
    
    std::ofstream report_file("comparison_report.html");
    report_file << html_report;
    report_file.close();
    
    std::cout << "Comparison report saved to: comparison_report.html" << std::endl;
    
    // Keep server running
    std::cout << "Web server is still running. Press Ctrl+C to stop..." << std::endl;
    
    while (running) {
      std::this_thread::sleep_for(std::chrono::seconds(1));
    }
    
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
  
  return 0;
}