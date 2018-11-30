//
// Copyright (c) 2013 Juan Palacios juan.palacios.puyana@gmail.com
// Subject to the BSD 2-Clause License
// - see < http://opensource.org/licenses/BSD-2-Clause>
//

#ifndef CONCURRENT_PAIR_HASH_MAP_
#define CONCURRENT_PAIR_HASH_MAP_

#include <unordered_map>
#include <thread>
#include <mutex>
#include <condition_variable>

template <typename T, typename U, typename V, typename W>
class Concurrent_Pair_Hash_Map
{
public:
U&
at(T& k)
{
  std::unique_lock<std::mutex>  mlock(mutex_);
  U&                            res = hashmap.at(k);
  mlock.unlock();
  cond_.notify_one();
  return res;
}


typename
std::unordered_map<T, U, V, W>::iterator
find(T& k)
{
  std::unique_lock<std::mutex> mlock(mutex_);
  typename std::unordered_map<T, U, V, W>::iterator res = hashmap.find(k);
  mlock.unlock();
  cond_.notify_one();
  return res;
}


typename
std::unordered_map<T, U, V, W>::iterator
begin() noexcept
{
  std::unique_lock<std::mutex> mlock(mutex_);
  typename std::unordered_map<T, U, V, W>::iterator res = hashmap.begin();
  mlock.unlock();
  cond_.notify_one();
  return res;
}


typename
std::unordered_map<T, U, V, W>::iterator
end() noexcept
{
  std::unique_lock<std::mutex> mlock(mutex_);
  typename std::unordered_map<T, U, V, W>::iterator res = hashmap.end();
  mlock.unlock();
  cond_.notify_one();
  return res;
};

void
clear() noexcept
{
  std::unique_lock<std::mutex> mlock(mutex_);
  hashmap.clear();
  mlock.unlock();
  cond_.notify_one();
}


U&
operator[](T& k)
{
  std::unique_lock<std::mutex>  mlock(mutex_);
  U&                            res = hashmap[k];
  mlock.unlock();
  cond_.notify_one();
  return res;
}


typename
std::unordered_map<T, U, V, W>&
operator=(const std::unordered_map<T, U, V, W>& ump)
{
  return ump;
}


private:
std::unordered_map<T, U, V, W>  hashmap;
std::mutex                      mutex_;
std::condition_variable         cond_;
};


#endif
