//===----------------------------------------------------------------------===//
//
//                         BusTub
//
// count_min_sketch.cpp
//
// Identification: src/primer/count_min_sketch.cpp
//
// Copyright (c) 2015-2025, Carnegie Mellon University Database Group
//
//===----------------------------------------------------------------------===//

#include "primer/count_min_sketch.h"

#include <queue>
#include <stdexcept>
#include <string>

namespace bustub {

/**
 * Constructor for the count-min sketch.
 *
 * @param width The width of the sketch matrix.
 * @param depth The depth of the sketch matrix.
 * @throws std::invalid_argument if width or depth are zero.
 */
template <typename KeyType>
CountMinSketch<KeyType>::CountMinSketch(uint32_t width, uint32_t depth) : width_(width), depth_(depth) {
  /** @TODO(student) Implement this function! */
  if (width <= 0 || depth <= 0) {
    throw std::invalid_argument("Invalid argument value for CountMinSketch, width and depth must be greater than 0");
  }
  data_ = std::make_unique<std::atomic<uint32_t>[]>(depth * width);

  /** @fall2025 PLEASE DO NOT MODIFY THE FOLLOWING */
  // Initialize seeded hash functions
  hash_functions_.reserve(depth_);
  for (size_t i = 0; i < depth_; i++) {
    hash_functions_.push_back(this->HashFunction(i));
  }
}

template <typename KeyType>
CountMinSketch<KeyType>::CountMinSketch(CountMinSketch &&other) noexcept : width_(other.width_), depth_(other.depth_) {
  /** @TODO(student) Implement this function! */
  data_ = std::move(other.data_);
  hash_functions_ = std::move(other.hash_functions_);
}

template <typename KeyType>
auto CountMinSketch<KeyType>::operator=(CountMinSketch &&other) noexcept -> CountMinSketch & {
  /** @TODO(student) Implement this function! */
  if (this != &other) {
    data_ = std::move(other.data_);
    hash_functions_ = std::move(other.hash_functions_);
  }
  return *this;
}

template <typename KeyType>
void CountMinSketch<KeyType>::Insert(const KeyType &item) {
  /** @TODO(student) Implement this function! */
  for (uint32_t row = 0; row < this->depth_; ++row) {
    auto col = this->hash_functions_[row](item);
    data_[row * width_ + col].fetch_add(1, std::memory_order_relaxed);
  }
}

template <typename KeyType>
void CountMinSketch<KeyType>::Merge(const CountMinSketch<KeyType> &other) {
  if (width_ != other.width_ || depth_ != other.depth_) {
    throw std::invalid_argument("Incompatible CountMinSketch dimensions for merge.");
  }
  /** @TODO(student) Implement this function! */
  for (uint32_t row = 0; row < this->depth_; ++row) {
    for (uint32_t col = 0; col < this->width_; col++) {
      auto index = row * width_ + col;
      auto other_val = other.data_[index].load();
      data_[row * width_ + col].fetch_add(other_val, std::memory_order_relaxed);
    }
  }
}

template <typename KeyType>
auto CountMinSketch<KeyType>::Count(const KeyType &item) const -> uint32_t {
  uint32_t count = UINT32_MAX;
  for (uint32_t row = 0; row < this->depth_; ++row) {
    auto col = this->hash_functions_[row](item);
    count = std::min(data_[row * width_ + col].load(), count);
  }
  return count;
}

template <typename KeyType>
void CountMinSketch<KeyType>::Clear() {
  /** @TODO(student) Implement this function! */
  for (uint32_t row = 0; row < depth_ * width_; row++) {
    data_[row].store(0, std::memory_order_relaxed);
  }
}

template <typename KeyType>
auto CountMinSketch<KeyType>::TopK(uint16_t k, const std::vector<KeyType> &candidates)
    -> std::vector<std::pair<KeyType, uint32_t>> {
  /** @TODO(student) Implement this function! */
  auto cmp = [](const std::pair<KeyType, uint32_t> &a, const std::pair<KeyType, uint32_t> &b) {
    return a.second < b.second;
  };
  std::priority_queue<std::pair<KeyType, uint32_t>, std::vector<std::pair<KeyType, uint32_t>>, decltype(cmp)> pq(cmp);
  std::vector<std::pair<KeyType, uint32_t>> res(k);
  for (auto &candidate : candidates) {
    pq.push(std::make_pair(candidate, Count(candidate)));
  }
  for (int i = 0; i < k; ++i) {
    res[i] = pq.top();
    pq.pop();
  }
  return res;
}

// Explicit instantiations for all types used in tests
template class CountMinSketch<std::string>;
template class CountMinSketch<int64_t>;  // For int64_t tests
template class CountMinSketch<int>;      // This covers both int and int32_t
}  // namespace bustub
