#include <stdexcept>
#include <vector>

template<typename T>
class Deque {
 private:
  struct position {
    int index;
    int block;

    position& operator+=(int value) {
      if (value < 0) {
        return operator-=(-value);
      }
      index += value;
      block += index / kSize;
      index %= kSize;
      return *this;
    }

    position& operator-=(int value) {
      if (value < 0) {
        return operator+=(-value);
      }
      index = block * kSize + index - value;
      block = index / kSize;
      index %= kSize;
      return *this;
    }

    position& operator++() {
      return operator+=(1);
    }

    position& operator--() {
      return operator-=(1);
    }

    position operator++(int) {
      position copy = *this;
      operator++();
      return copy;
    }

    position operator--(int) {
      position copy = *this;
      operator--();
      return copy;
    }

    int operator-(const position& other) const {
      int difference =
          index + block * kSize - other.index - other.block * kSize;
      return difference;
    }

    position operator+(int value) const {
      position copy = *this;
      copy += value;
      return copy;
    }

    position operator-(int value) const {
      position copy = *this;
      copy -= value;
      return copy;
    }

    std::strong_ordering operator<=>(const position& other) const {
      if (block < other.block) {
        return std::strong_ordering::less;
      }
      if (block > other.block) {
        return std::strong_ordering::greater;
      }
      if (index < other.index) {
        return std::strong_ordering::less;
      }
      if (index > other.index) {
        return std::strong_ordering::greater;
      }
      return std::strong_ordering::equivalent;
    }

    bool operator==(const position& other) const {
      return block == other.block && index == other.index;
    }

    bool operator!=(const position& other) const {
      return !operator==(other);
    }
  };


 public:
  template<bool IsConst>
  class base_iterator {
   public:
    using value_type = std::conditional_t<IsConst, const T, T>;
    using reference_type = std::conditional_t<IsConst, const T&, T&>;
    using pointer_type = std::conditional_t<IsConst, const T*, T*>;
    using iterator_category = std::random_access_iterator_tag;
    using difference_type = std::ptrdiff_t;

    base_iterator(T** data, position pos) : data_(data), it(pos) {}

    reference_type operator*() {
      return data_[it.block][it.index];
    }

    const T& operator*() const {
      return data_[it.block][it.index];
    };

    base_iterator& operator++() {
      it++;
      return *this;
    }

    base_iterator operator++(int) {
      base_iterator copy = *this;
      it++;
      return copy;
    }

    base_iterator& operator--() {
      it--;
      return *this;
    }

    base_iterator operator--(int) {
      base_iterator copy = *this;
      it--;
      return copy;
    }

    base_iterator& operator+=(difference_type value) {
      it += value;
      return *this;
    }

    base_iterator& operator-=(difference_type value) {
      it -= value;
      return *this;
    }

    base_iterator operator+(difference_type value) const {
      base_iterator copy = *this;
      copy += value;
      return copy;
    }

    base_iterator operator-(difference_type value) const {
      base_iterator copy = *this;
      copy -= value;
      return copy;
    }

    difference_type operator-(const base_iterator& other) const {
      return it - other.it;
    }

    bool operator<(const base_iterator& other) const {
      if (other.it.block != it.block) {
        return it.block < other.it.block;
      }
      return it.index < other.it.index;
    }

    std::strong_ordering operator<=>(const base_iterator& other) const {
      return it <=> other.it;
    }

    bool operator==(const base_iterator& other) const {
      return it == other.it;
    }

    bool operator!=(const base_iterator& other) const {
      return it != other.it;
    }

    pointer_type operator->() {
      return data_[it.block] + it.index;
    }

    const T* operator->() const {
      return data_[it.block] + it.index;
    }

    base_iterator& operator=(const base_iterator<false>& other) {
      data_ = other.data_;
      it = other.it;
      return *this;
    }

    base_iterator(const base_iterator<false>& other) : data_(other.data_),
        it(other.it) {}


    friend class Deque;

   private:
    T** data_;
    position it;
  };

  Deque() : data_(nullptr), size_(0), block_count(0), begin_{0, 0},
      end_{0, 0} {}

  Deque(int n) : data_(allocate((n + kSize - 1) / kSize * 3)), size_(n),
      block_count((n + kSize - 1) / kSize * 3),
      begin_{.index = 0, .block = block_count / 3},
      end_(begin_ + size_) {
    fill(data_, size_, begin_);
  }

  Deque(int n, const T& value) : data_(allocate((n + kSize - 1) / kSize * 3)),
      size_(n),
      block_count((n + kSize - 1) / kSize * 3),
      begin_{.index = 0, .block = block_count / 3},
      end_(begin_ + size_) {
    fill(data_, size_, begin_, value);
  }

  Deque(const Deque& other) : data_(allocate(other.block_count)),
      size_(other.size_), block_count(other.block_count),
      begin_(other.begin_), end_(other.end_) {
    copy(other.data_, data_, begin_, size_);
  }

  Deque& operator=(Deque other) {
    std::swap(data_, other.data_);
    std::swap(size_, other.size_);
    std::swap(block_count, other.block_count);
    std::swap(begin_, other.begin_);
    std::swap(end_, other.end_);
    return *this;
  }

  int size() const {
    return size_;
  }

  T& operator[](int index) {
    position pos = begin_ + index;
    return data_[pos.block][pos.index];
  }

  const T& operator[](int index) const {
    position pos = begin_ + index;
    return data_[pos.block][pos.index];
  }

  T& at(int index) {
    if (index >= size_ || index < 0) {
      throw std::out_of_range("index out of range");
    }
    position pos = begin_ + index;
    return data_[pos.block][pos.index];
  }

  const T& at(int index) const {
    position pos = begin_ + index;
    return data_[pos.block][pos.index];
  }

  void push_back(const T& value) {
    if (end_.block == block_count || data_ == nullptr) {
      data_ = transfer();
    }
    new(data_[end_.block] + end_.index) T(value);
    ++end_;
    size_++;
  }

  void pop_back() {
    --end_;
    data_[end_.block][end_.index].~T();
    size_--;
  }

  void push_front(const T& value) {
    if (begin_.block == 0 && begin_.index == 0) {
      data_ = transfer();
    }
    begin_--;
    try {
      new(data_[begin_.block] + begin_.index) T(value);
    } catch (...) {
      begin_++;
      throw;
    }
    size_++;
  }

  void pop_front() {
    data_[begin_.block][begin_.index].~T();
    ++begin_;
    size_--;
  }

  using iterator = base_iterator<false>;
  using const_iterator = base_iterator<true>;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  iterator begin() {
    return iterator{data_, begin_};
  }

  const_iterator begin() const {
    return const_iterator{data_, begin_};
  }

  iterator end() {
    return iterator{data_, end_};
  }

  const_iterator end() const {
    return const_iterator{data_, end_};
  }

  const_iterator cbegin() const {
    return const_iterator{data_, begin_};
  }

  const_iterator cend() const {
    return const_iterator{data_, end_};
  }

  reverse_iterator rbegin() const {
    return reverse_iterator(iterator{data_, end_});
  }

  reverse_iterator rend() const {
    return reverse_iterator(iterator{data_, begin_});
  }

  const_reverse_iterator crbegin() const {
    return const_reverse_iterator(const_iterator(data_, end_));
  }

  const_reverse_iterator crend() const {
    return const_reverse_iterator(const_iterator(data_, begin_));
  }

  void insert(const const_iterator& iter, const T& value) {
    int idx = iter.it - begin_;
    int cur = size_;
    push_back(value);
    position pos = end_ - 1;
    while (idx != cur) {
      position prev = pos - 1;
      try {
        std::swap(data_[pos.block][pos.index], data_[prev.block][prev.index]);
      } catch (...) {

      }
      cur--;
      pos = prev;
    }
  }

  void erase(const const_iterator& iter) {
    position current = iter.it;
    while (current != end_ - 1) {
      position next = current + 1;
      std::swap(data_[next.block][next.index], data_[current.block][current.index]);
      current = next;
    }
    pop_back();
  }

  ~Deque() {
    for (int i = 0; i < block_count; ++i) {
      delete[] reinterpret_cast<char*>(data_[i]);
    }
    delete[] data_;
  }

 private:
  T** data_;
  static const int kSize = 128;
  int size_;
  int block_count;
  position begin_, end_;

  T** allocate(int count) {
    T** data = new T* [count];
    int last = 0;
    try {
      for (; last < count; ++last) {
        data[last] = reinterpret_cast<T*>(new char[kSize * sizeof(T)]);
      }
    } catch (...) {
      for (int i = 0; i < last; ++i) {
        delete[] reinterpret_cast<char*>(data[last]);
      }
      delete[] data;
      throw;
    }
    return data;
  }

  static void fill(T** data, int count, position pos, const T& value = T()) {
    int cnt = 0;
    try {
      for (; cnt < count; ++cnt) {
        new(data[pos.block] + pos.index) T(value);
        ++pos;
      }
    } catch (...) {
      --pos;
      for (int i = 0; i < cnt; ++i) {
        data[pos.block][pos.index].~T();
        --pos;
      }
      throw;
    }
  }

  static void copy(T** from, T** to, position start, int count) {
    int cnt = 0;
    try {
      for (; cnt < count; ++cnt) {
        new(to[start.block] + start.index) T(from[start.block][start.index]);
        ++start;
      }
    } catch (...) {
      --start;
      for (int i = 0; i < cnt; ++i) {
        to[start.block][start.index].~T();
        start--;
      }
      throw;
    }
  }

  T** transfer() {
    if (data_ == nullptr) {
      begin_ = {.index = 0, .block = 1};
      end_ = {.index = 0, .block = 1};
      size_ = 0;
      block_count = 3;
      return allocate(3);
    }

    T** new_data = new T* [block_count * 3];
    int cnt = 0;

    try {
      for (; cnt < block_count; ++cnt) {
        new_data[cnt] = reinterpret_cast<T*>(new char[kSize * sizeof(T)]);
      }
    } catch (...) {
      for (int i = 0; i < cnt; ++i) {
        delete[] reinterpret_cast<char*>(new_data[i]);
      }
      throw;
    }
    for (int i = 0; i < block_count; ++i) {
      new_data[i + block_count] = data_[i];
    }
    cnt = 0;
    try {
      for (; cnt < block_count; ++cnt) {
        new_data[cnt + 2 * block_count] = reinterpret_cast<T*>(new char[kSize *
                                                                        sizeof(T)]);
      }
    } catch (...) {
      for (int i = 0; i < cnt; ++i) {
        delete[] reinterpret_cast<char*>(new_data[i + 2 * block_count]);
      }
      throw;
    }
    delete[] data_;
    block_count *= 3;
    begin_.block += block_count / 3;
    end_.block += block_count / 3;
    return new_data;
  }

};

