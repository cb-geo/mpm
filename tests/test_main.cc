#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

int main(int argc, char* argv[]) {
  Catch::Session session;  // There must be exactly one instance
  std::cout << "TESTING CATCH MAIN\n";
  // Let Catch (using Clara) parse the command line
  int returnCode = session.applyCommandLine(argc, argv);
  if (returnCode != 0)  // Indicates a command line error
    return returnCode;
  std::cout << "Finishing CATCH MAIN\n";
  int result = session.run();
  std::cout << "Finishing Result\t" << result << "\n";
  return result;
}
