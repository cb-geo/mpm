#ifndef _FACTORY_H_
#define _FACTORY_H_

#include <functional>
#include <map>
#include <memory>
#include <string>

//! \brief Singleton factory implementation
template <typename Tbaseclass, typename... Targs>
class Factory {
 public:
  // Get the single instance of the factory
  static Factory* instance() {
    static Factory factory;
    return &factory;
  }

  // register a factory function to create an instance of classname
  template <typename Tderivedclass>
  void register_factory(const std::string& key) {
    registry[key].reset(new Creator<Tderivedclass>);
  }

  // create an instance of a registered class
  std::shared_ptr<Tbaseclass> create(const std::string& key, Targs&&... args) {
    return registry.at(key)->create(std::forward<Targs>(args)...);
  }

  // List registered elements
  std::vector<std::string> list() const {
    std::vector<std::string> factory_items;
    for (const auto& keyvalue : registry)
      factory_items.push_back(keyvalue.first);
    return factory_items;
  }

 private:
  // Private constructor
  Factory() = default;

  struct CreatorBase {
    virtual std::shared_ptr<Tbaseclass> create(Targs&&...) = 0;
  };
  template <typename Tderivedclass>
  struct Creator : public CreatorBase {
    std::shared_ptr<Tbaseclass> create(Targs&&... args) override {
      return std::make_shared<Tderivedclass>(std::forward<Targs>(args)...);
    }
  };
  // Register of factory functions
  std::map<std::string, std::shared_ptr<CreatorBase>> registry;
};

// A helper class to register a factory function
template <typename Tbaseclass, typename Tderivedclass, typename... Targs>
class Register {
 public:
  explicit Register(const std::string& key) {
    // register the class factory function
    Factory<Tbaseclass, Targs...>::instance()
        ->template register_factory<Tderivedclass>(key);
  }
};

#endif  // _FACTORY_H_
