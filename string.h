#include <algorithm>
#include <iostream>
#include <cstring>

class String {
 public:
  String() : size_(0), capacity_(1), data_(new char[1]) {
    data_[size_] = '\0';
  }

  String(const char* string) : size_(strlen(string)), capacity_(size_ + 1),
      data_(new char[size_ + 1]) {
    std::copy(string, string + size_, data_);
    data_[size_] = '\0';
  }

  String(char sym) : String(1, sym) {
  }

  String(size_t size, char sym) : size_(size), capacity_(size + 1),
      data_(new char[size + 1]) {
    std::fill(data_, data_ + size_, sym);
    data_[size] = '\0';
  }

  String(const String& other) : size_(other.size_), capacity_(other.capacity_),
      data_(new char[other.capacity_]) {
    std::copy(other.data_, other.data_ + size_, data_);
    data_[size_] = '\0';
  }

  String& operator=(String other) {
    swap(other);
    return *this;
  }

  char& operator[](size_t index) { return data_[index]; }

  const char& operator[](size_t index) const { return data_[index]; }

  size_t size() const { return size_; };

  size_t capacity() const { return capacity_ - 1; }

  size_t length() const { return size_; }

  void push_back(char sym) {
    if (size_ + 1 == capacity_) {
      reallocation(2 * capacity_);
      capacity_ *= 2;
    }
    data_[size_] = sym;
    size_++;
    data_[size_] = '\0';
  }

  void pop_back() {
    size_--;
    data_[size_] = '\0';
  }

  char& front() { return data_[0]; }

  const char& front() const { return data_[0]; }

  char& back() { return data_[size_ - 1]; }

  const char& back() const { return data_[size_ - 1]; }

  String& operator+=(const String& other) {
    if (capacity_ <= size_ + other.size_) {
      reallocation(std::max(2 * capacity_, 2 * (size_ + other.size_ + 1)));
      capacity_ = std::max(2 * capacity_, 2 * (size_ + other.size_ + 1));
    }
    std::copy(other.data_, other.data_ + other.size_, data_ + size_);
    size_ += other.size_;
    data_[size_] = '\0';
    return *this;
  }

  String& operator+=(char sym) {
    push_back(sym);
    return *this;
  }

  size_t find(const String& substring) const {
    if (substring.size() > size_) {
      return size_;
    }
    for (size_t i = 0; i <= size_; ++i) {
      if (check_equal(i, substring)) {
        return i;
      }
    }
    return size_;
  }

  size_t rfind(const String& substring) const {
    if (substring.size() > size_) {
      return size_;
    }
    for (int i = static_cast<int>(size_); i >= 0; --i) {
      if (check_equal(i, substring)) {
        return i;
      }
    }
    return size_;
  }

  String substr(size_t start, size_t count) const {
    if (start + count > size_) {
      count = size_ - start;
    }
    String ret(count, '\0');
    std::copy(data_ + start, data_ + start + count, ret.data_);
    ret.size_ = count;
    ret.capacity_ = count + 1;
    ret.data_[count] = '\0';
    return ret;
  }

  bool empty() const { return size_ == 0; }

  void clear() {
    size_ = 0;
    data_[size_] = '\0';
  }

  void shrink_to_fit() {
    if (size_ + 1 == capacity_) {
      return;
    }
    capacity_ = size_ + 1;
    char* new_data = new char[size_ + 1];

    std::copy(data_, data_ + size_, new_data);
    delete[] data_;
    data_ = new_data;
    data_[size_] = '\0';
  }

  const char* data() const { return data_; }

  char* data() { return data_; }

  ~String() {
    delete[] data_;
  }

 private:
  void swap(String& other) {
    std::swap(data_, other.data_);
    std::swap(size_, other.size_);
    std::swap(capacity_, other.capacity_);
  }

  void reallocation(size_t capacity) {
    char* new_data = new char[capacity];
    std::copy(data_, data_ + size_, new_data);
    delete[] data_;
    data_ = new_data;
  }

  bool check_equal(size_t index, const String& substring) const {
    if (index + substring.size() > size_) {
      return false;
    }
    return memcmp(data_ + index, substring.data_, substring.size_) == 0;
  }

  size_t size_;
  size_t capacity_;
  char* data_;
};

std::istream& operator>>(std::istream& in, String& str) {
  str.clear();
  int cur;
  cur = in.get();
  size_t size = 0;
  while ((!std::isspace(cur) || size == 0) && cur != -1) {
    if (std::isspace(cur) && size != 0) {
      break;
    }
    if (std::isspace(cur)) {
      cur = in.get();
      continue;
    }
    str.push_back(cur);
    ++size;
    cur = in.get();
  }

  return in;
}

std::ostream& operator<<(std::ostream& out, const String& str) {
  for (size_t i = 0; i < str.size(); ++i) {
    out << str[i];
  }
  return out;
}

String operator+(const String& first, const String& second) {
  String copy(first.size() + second.size(), 'a');
  copy.clear();
  copy = first;
  copy += second;
  return copy;
}

bool operator<(const String& first, const String& second) {
  if (first.size() <= second.size()) {
    return memcmp(first.data(), second.data(), first.size()) < 0;
  }
  return memcmp(second.data(), first.data(), second.size()) > 0;
}

bool operator>(const String& first, const String& second) {
  return second < first;
}

bool operator<=(const String& first, const String& second) {
  return !(second > first);
}

bool operator>=(const String& first, const String& second) {
  return !(second < first);
}

bool operator==(const String& first, const String& second) {
  if (first.size() != second.size()) {
    return false;
  }
  return memcmp(first.data(), second.data(), first.size()) == 0;
}

bool operator!=(const String& first, const String& second) {
  return !(first == second);
}
