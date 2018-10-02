#ifndef _FACTORY_H_
#define _FACTORY_H_

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

//! \brief Singleton factory implementation
//! \tparam Tbaseclass Base class
//! \tparam Targs variadic template arguments
template <typename Tbaseclass, typename... Targs>
class Factory {
 public:
  //! Get the single instance of the factory
  static Factory* instance() {
    static Factory factory;
    return &factory;
  }

  //! Register a factory function to create an instance of classname
  //! \param[in] key Register key
  //! \tparam Tderivedclass Derived class
  template <typename Tderivedclass>
  void register_factory(const std::string& key) {
    registry[key].reset(new Creator<Tderivedclass>);
  }

  //! Create an instance of a registered class
  //! \param[in] key key to item in registry
  //! \param[in] args Variadic template arguments
  //! \retval shared_ptr<Tbaseclass> Shared pointer to a base class
  std::shared_ptr<Tbaseclass> create(const std::string& key, Targs&&... args) {
    return registry.at(key)->create(std::forward<Targs>(args)...);
  }

  //! List registered elements
  //! \retval factory_items Return list of items in the registry
  std::vector<std::string> list() const {
    std::vector<std::string> factory_items;
    for (const auto& keyvalue : registry)
      factory_items.push_back(keyvalue.first);
    return factory_items;
  }

  //! Check if an element is registered
  //! \param[in] item Item to be checked in registry
  //! \retval status Return true if element is registered or false otherwise
  bool check(const std::string& item) const {
    bool status = false;
    for (const auto& keyvalue : registry)
      if (keyvalue.first == item) status = true;
    return status;
  }

 private:
  // Private constructor
  Factory() = default;

  //! A base class creator struct
  struct CreatorBase {
    //! A virtual create function
    virtual std::shared_ptr<Tbaseclass> create(Targs&&...) = 0;
  };

  //! Creator class
  //! \tparam Tderivedclass Derived class
  template <typename Tderivedclass>
  struct Creator : public CreatorBase {
    //! Create instance of object
    std::shared_ptr<Tbaseclass> create(Targs&&... args) override {
      return std::make_shared<Tderivedclass>(std::forward<Targs>(args)...);
    }
  };
  // Register of factory functions
  std::map<std::string, std::shared_ptr<CreatorBase>> registry;
};

//! A helper class to register a factory function
//! \tparam Tbaseclass Base class
//! \tparam Tderivedclass Derived class
//! \tparam Targs variadic template arguments
template <typename Tbaseclass, typename Tderivedclass, typename... Targs>
class Register {
 public:
  //! Register with a given key
  //! \param[in] key Key to item in registry
  explicit Register(const std::string& key) {
    // register the class factory function
    Factory<Tbaseclass, Targs...>::instance()
        ->template register_factory<Tderivedclass>(key);
  }
};

#endif  // _FACTORY_H_
