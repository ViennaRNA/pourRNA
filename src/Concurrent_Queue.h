//
// Copyright (c) 2013 Juan Palacios juan.palacios.puyana@gmail.com
// Subject to the BSD 2-Clause License
// - see < http://opensource.org/licenses/BSD-2-Clause>
//

#ifndef CONCURRENT_QUEUE_
#define CONCURRENT_QUEUE_

#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>

template <typename T>
class Concurrent_Queue
{
public:

T
pop()
{
  std::unique_lock<std::mutex>  mlock(mutex_);

  while (queue_.empty())
    cond_.wait(mlock);
  auto                          val = queue_.front();
  queue_.pop();
  return val;
}


void
pop(T& item)
{
  std::unique_lock<std::mutex> mlock(mutex_);

  while (queue_.empty())
    cond_.wait(mlock);
  item = queue_.front();
  queue_.pop();
}


void
push(const T& item)
{
  std::unique_lock<std::mutex> mlock(mutex_);

  queue_.push(item);
  mlock.unlock();
  cond_.notify_one();
}


bool
empty()
{
  std::unique_lock<std::mutex>  mlock(mutex_);
  bool                          status = queue_.empty();

  mlock.unlock();
  cond_.notify_one();
  return status;
}


void
clear()
{
  std::unique_lock<std::mutex> mlock(mutex_);

  while (!queue_.empty())
    queue_.pop();
  mlock.unlock();
  cond_.notify_one();
}


Concurrent_Queue()                        = default;
Concurrent_Queue(const Concurrent_Queue&) = delete;                   // disable copying
Concurrent_Queue&
operator=(const Concurrent_Queue&) = delete;                          // disable assignment


private:
std::queue<T> queue_;
std::mutex mutex_;
std::condition_variable cond_;
};

#endif
