#include <memory>

template<size_t N>
class StackStorage {
 public:
  StackStorage() : storage_(), end_(storage_), space_(N) {}

  StackStorage(const StackStorage& other) = delete;

  template<typename T>
  void* allocate(size_t size);

  void* get_pointer() {
    return end_;
  }

  size_t get_space() {
    return space_;
  }

  void set_pointer(void* end) {
    end_ = end;
  }

  void set_space(size_t space) {
    space_ = space;
  }

 private:
  alignas(std::max_align_t) char storage_[N];
  void* end_;
  size_t space_;
};

template<size_t N>
template<typename T>
void* StackStorage<N>::allocate(size_t size) {
  if (std::align(alignof(T), size * sizeof(T), end_, space_)) {
    void* ret = end_;
    end_ = reinterpret_cast<void*>(reinterpret_cast<char*>(end_) +
                                   size * sizeof(T));
    space_ -= size * sizeof(T);
    return ret;
  }
  throw std::bad_alloc();
}

template<typename T, size_t N>
class StackAllocator {
 public:
  using value_type = T;

  StackAllocator(StackStorage<N>& storage) : storage_(&storage) {}

  T* allocate(size_t size) {
    void* end = storage_->get_pointer();
    size_t space = storage_->get_space();
    if (std::align(alignof(T), size * sizeof(T), end, space)) {
      void* ret = end;
      end = reinterpret_cast<void*>(reinterpret_cast<char*>(end) +
                                     size * sizeof(T));
      space -= size * sizeof(T);
      storage_->set_space(space);
      storage_->set_pointer(end);
      return reinterpret_cast<T*>(ret);
    }
    throw std::bad_alloc();
  }

  void deallocate(T*, size_t) {}

  template<typename T1, size_t M>
  StackAllocator(const StackAllocator<T1, M>& alloc) : storage_(
      alloc.storage_) {}

  template<typename T1>
  struct rebind {
    using other = StackAllocator<T1, N>;
  };
  StackStorage<N>* storage_;
};

template<typename T, typename Allocator = std::allocator<T>>
class List {
 private:
  struct BaseNode {
    // public:
    BaseNode* next, * prev;

    BaseNode() : next(nullptr), prev(nullptr) {}

    ~BaseNode() = default;

    BaseNode(BaseNode* next, BaseNode* prev) : next(next), prev(prev) {}
  };

  struct Node : BaseNode {
    // public:
    T value;

    Node() = default;

    ~Node() = default;

    Node(BaseNode* next, BaseNode* prev, const T& value) : BaseNode(next, prev),
        value(value) {}

    Node(const T& value) : BaseNode(), value(value) {}
  };

  [[no_unique_address]] typename std::allocator_traits<Allocator>::template rebind_alloc<Node> alloc_;
  using AllocTraits = std::allocator_traits<typename std::allocator_traits<Allocator>::template rebind_alloc<Node>>;

  //  using AllocTraits = typename std::allocator_traits<Allocator>::template rebind_traits<Node>;
  //  [[no_unique_address]] AllocTraits::allocator_type alloc_;
  BaseNode fake_node_;
  size_t size_;


  template<bool IsConst>
  class base_iterator {
   public:
    using value_type = std::conditional_t<IsConst, const T, T>;
    using pointer = std::conditional_t<IsConst, const T*, T*>;
    using reference = std::conditional_t<IsConst, const T&, T&>;
    using iterator_category = std::bidirectional_iterator_tag;
    using difference_type = std::ptrdiff_t;

    friend class List<T, Allocator>;

    base_iterator& operator=(const base_iterator<false>& other) {
      node = other.node;
      return *this;
    }

    base_iterator(const base_iterator<false>& other) : node(other.node) {}

    bool operator==(const base_iterator& other) {
      return other.node == node;
    }

    bool operator!=(const base_iterator& other) {
      return !(*this == other);
    }

    /*explicit operator base_iterator<true>() const {
      return base_iterator<true>(node);
    }*/

    base_iterator& operator++() {
      node = node->next;
      return *this;
    }

    base_iterator& operator--() {
      node = node->prev;
      return *this;
    }

    base_iterator operator++(int) {
      base_iterator copy = *this;
      ++(*this);
      return copy;
    }

    base_iterator operator--(int) {
      base_iterator copy = *this;
      --(*this);
      return copy;
    }

    reference operator*() {
      return static_cast<Node*>(node)->value;
    }

    const reference operator*() const {
      return static_cast<Node*>(node)->value;
    }

    pointer operator->() {
      return node;
    }

    const pointer operator->() const {
      return node;
    }

   private:
    base_iterator(BaseNode* node) : node(node) {}

    BaseNode* node;
  };


 public:
  using iterator = base_iterator<false>;
  using const_iterator = base_iterator<true>;
  using reverse_iterator = std::reverse_iterator<base_iterator<false>>;
  using const_reverse_iterator = std::reverse_iterator<base_iterator<true>>;

  List() : fake_node_(&fake_node_, &fake_node_), size_(0) {}

  List(const Allocator& alloc) : alloc_(alloc),
      fake_node_(&fake_node_, &fake_node_), size_(0) {}

  List(size_t size_);

  List(size_t size, const Allocator& alloc)
      : alloc_(alloc), fake_node_(&fake_node_, &fake_node_), size_(0) {
    size_t i = 0;
    try {
      for (; i < size; ++i) {
        insert_default(end());
      }
    } catch (...) {
      for (size_t j = 0; j < i; ++i) {
        pop_back();
      }
      throw;
    }
  }

  List(size_t size, const T& value);

  List(size_t size, const T& value, const Allocator& alloc) : alloc_(alloc),
      List(size, value) {}

  List(const List<T, Allocator>& other)
      : alloc_(
      AllocTraits::select_on_container_copy_construction(other.alloc_)),
      fake_node_(&fake_node_, &fake_node_), size_(0) {
    size_t i = 0;
    auto cur_it = other.begin();
    try {
      for (; i < other.size_; ++i) {
        push_back(*cur_it);
        cur_it++;
      }
    } catch (...) {
      for (size_t j = 0; j < i; ++j) {
        pop_back();
      }
      throw;
    }
  }

  List& operator=(const List<T, Allocator>& other);

  void push_back(const T& value);

  void pop_back();

  AllocTraits::allocator_type get_allocator();

  size_t size() const;

  void push_front(const T& value);

  void pop_front();

  void insert(const_iterator it, const T& value);

  void erase(const_iterator it);

  // iterators

  iterator begin();

  iterator end();

  const_iterator begin() const;

  const_iterator end() const;

  const_iterator cbegin() const;

  const_iterator cend() const;

  reverse_iterator rbegin();

  reverse_iterator rend();

  const_reverse_iterator rbegin() const;

  const_reverse_iterator rend() const;

  const_reverse_iterator crbegin() const;

  const_reverse_iterator crend() const;

  ~List();

 private:

  void set_pointers(BaseNode* node, Node* ins_node);

  void insert_default(const_iterator it);

};

// private

template<typename T, typename Allocator>
void List<T, Allocator>::set_pointers(List::BaseNode* node,
                                      List::Node* ins_node) {
  ins_node->prev = node->prev;
  ins_node->next = node;
  node->prev->next = ins_node;
  node->prev = ins_node;
}

template<typename T, typename Allocator>
void List<T, Allocator>::insert_default(const_iterator it) {
  BaseNode* node = it.node;
  Node* ins_node = AllocTraits::allocate(alloc_, 1);
  try {
    AllocTraits::construct(alloc_, ins_node);
  } catch (...) {
    AllocTraits::deallocate(alloc_, ins_node, 1);
    throw;
  }
  set_pointers(node, ins_node);
  ++size_;
}

// insert etc

template<typename T, typename Allocator>
void List<T, Allocator>::insert(const_iterator it, const T& value) {
  BaseNode* node = it.node;
  Node* ins_node = AllocTraits::allocate(alloc_, 1);
  try {
    AllocTraits::construct(alloc_, ins_node, value);
  } catch (...) {
    AllocTraits::deallocate(alloc_, ins_node, 1);
    throw;
  }
  ++size_;
  set_pointers(node, ins_node);
}

template<typename T, typename Allocator>
void List<T, Allocator>::push_back(const T& value) {
  insert(end(), value);
}

template<typename T, typename Allocator>
void List<T, Allocator>::pop_back() {
  iterator it = end();
  it--;
  erase(it);
}

template<typename T, typename Allocator>
void List<T, Allocator>::push_front(const T& value) {
  insert(begin(), value);
}

template<typename T, typename Allocator>
void List<T, Allocator>::pop_front() {
  erase(begin());
}

template<typename T, typename Allocator>
void List<T, Allocator>::erase(const_iterator it) {
  BaseNode* cur = it.node;
  BaseNode* prev = cur->prev;
  BaseNode* next = cur->next;
  if (cur == &fake_node_) {
    return;
  }
  AllocTraits::destroy(alloc_, static_cast<Node*>(it.node));
  AllocTraits::deallocate(alloc_, static_cast<Node*>(it.node), 1);
  prev->next = next;
  next->prev = prev;
  --size_;
}


// constructors

template<typename T, typename Allocator>
List<T, Allocator>::List(size_t size, const T& value) {
  fake_node_.prev = &fake_node_;
  fake_node_.next = &fake_node_;
  int cnt = 0;
  try {
    for (; cnt < size; ++cnt) {
      push_back(value);
    }
  } catch (...) {
    for (int i = 0; i < cnt; ++i) {
      pop_back();
    }
    throw;
  }
}

template<typename T, typename Allocator>
List<T, Allocator>::List(size_t size) {
  size_t i = 0;
  size_ = 0;
  fake_node_.prev = &fake_node_;
  fake_node_.next = &fake_node_;
  try {
    for (; i < size; ++i) {
      insert_default(end());
    }
  } catch (...) {
    for (size_t j = 0; j < i; ++j) {
      pop_back();
    }
    throw;
  }
}

template<typename T, typename Allocator>
List<T, Allocator>& List<T, Allocator>::operator=(
    const List<T, Allocator>& other) {
  typename AllocTraits::allocator_type alloc = alloc_;
  if (AllocTraits::propagate_on_container_copy_assignment::value) {
    while (size_ != 0) {
      pop_back();
    }
    alloc_ = other.alloc_;
  }
  size_t i = 0;
  size_t size = size_;
  try {
    auto cur_it = other.begin();
    for (; i < other.size_; ++i) {
      push_back(*cur_it);
      ++cur_it;
    }
  } catch (...) {
    for (size_t j = 0; j < i; ++j) {
      pop_back();
    }
    alloc_ = alloc;
    throw;
  }
  if (!AllocTraits::propagate_on_container_copy_assignment::value) {
    for (size_t j = 0; j < size; ++j) {
      pop_front();
    }
  }
  size_ = other.size_;
  return *this;
}

// smth
template<typename T, typename Allocator>
List<T, Allocator>::AllocTraits::allocator_type
List<T, Allocator>::get_allocator() {
  return alloc_;
}

template<typename T, typename Allocator>
size_t List<T, Allocator>::size() const {
  return size_;
}

// iterators

template<typename T, typename Allocator>
List<T, Allocator>::const_iterator List<T, Allocator>::begin() const {
  return const_iterator(const_cast<BaseNode*>(fake_node_.next));
}

template<typename T, typename Allocator>
List<T, Allocator>::iterator List<T, Allocator>::begin() {
  return iterator(fake_node_.next);
}

template<typename T, typename Allocator>
List<T, Allocator>::const_iterator List<T, Allocator>::cbegin() const {
  return const_iterator(const_cast<BaseNode*>(fake_node_.next));
}

template<typename T, typename Allocator>
List<T, Allocator>::const_reverse_iterator List<T, Allocator>::rbegin() const {
  return reverse_iterator(iterator(const_cast<BaseNode*>(&fake_node_)));
}

template<typename T, typename Allocator>
List<T, Allocator>::reverse_iterator List<T, Allocator>::rbegin() {
  return reverse_iterator(iterator(&fake_node_));
}

template<typename T, typename Allocator>
List<T, Allocator>::const_reverse_iterator List<T, Allocator>::crbegin() const {
  return reverse_iterator(iterator(const_cast<BaseNode*>(&fake_node_)));
}

template<typename T, typename Allocator>
List<T, Allocator>::const_iterator List<T, Allocator>::end() const {
  return const_iterator(const_cast<BaseNode*>(&fake_node_));
}

template<typename T, typename Allocator>
List<T, Allocator>::iterator List<T, Allocator>::end() {
  return iterator(&fake_node_);
}

template<typename T, typename Allocator>
List<T, Allocator>::const_iterator List<T, Allocator>::cend() const {
  return const_iterator(const_cast<BaseNode*>(&fake_node_));
}

template<typename T, typename Allocator>
List<T, Allocator>::const_reverse_iterator List<T, Allocator>::rend() const {
  return reverse_iterator(iterator(const_cast<BaseNode*>(fake_node_.prev)));
}

template<typename T, typename Allocator>
List<T, Allocator>::reverse_iterator List<T, Allocator>::rend() {
  return reverse_iterator(iterator(fake_node_.prev));
}

template<typename T, typename Allocator>
List<T, Allocator>::const_reverse_iterator List<T, Allocator>::crend() const {
  return reverse_iterator(iterator(const_cast<BaseNode*>(fake_node_.prev)));
}

template<typename T, typename Allocator>
List<T, Allocator>::~List() {
  auto cur = begin();
  size_t size = size_;
  for (size_t i = 0; i < size; ++i) {
    auto next = cur;
    std::advance(next, 1);
    erase(cur);
    cur = next;
  }
}
