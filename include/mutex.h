#ifndef _MUTEX_
#define _MUTEX_

#include <atomic>

//! Hybrid SpinMutex class
//! \brief SpinMutex class that offers a fast locking and unlocking mechanism
class SpinMutex {

 public:
  //! Attempt locking
  bool try_lock() { return !lock_.test_and_set(std::memory_order_acquire); }

  //! Call lock 
  void lock() {
    std::size_t spin_count{0};
    while (!try_lock()) {
      ++spin_count;

      if (spin_count < spin_pred_ * 2) continue;

      std::this_thread::sleep_for(std::chrono::nanoseconds(1));
    }

    spin_pred_ += (spin_count - spin_pred_) / 8;
  }

  //! Call unlock
  void unlock() { lock_.clear(std::memory_order_release); }

 private:
  //! Lock variable
  std::atomic_flag lock_ = ATOMIC_FLAG_INIT;
  //! Spin prediction
  std::atomic<std::size_t> spin_pred_{0};
};
#endif  // _MUTEX_
