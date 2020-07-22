#ifndef _MUTEX_
#define _MUTEX_

#include <atomic>

class SpinMutex {
 private:
  std::atomic_flag _lock = ATOMIC_FLAG_INIT;
  std::atomic<std::size_t> _spin_pred{0};

 public:
  bool try_lock() { return !_lock.test_and_set(std::memory_order_acquire); }

  void lock() {
    std::size_t spin_count{0};
    while (!try_lock()) {
      ++spin_count;

      if (spin_count < _spin_pred * 2) continue;

      // REVISIT (fbrereto) : Iff sleep_for is guaranteed to block even
      // when the duration is 0, then std::chrono::nanoseconds(0) will
      // suffice here.
      std::this_thread::sleep_for(std::chrono::nanoseconds(1));
    }

    _spin_pred += (spin_count - _spin_pred) / 8;
  }

  void unlock() { _lock.clear(std::memory_order_release); }
};
#endif
