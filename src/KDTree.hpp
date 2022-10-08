// Copyright

#ifndef SRC_KDTREE_HPP_
#define SRC_KDTREE_HPP_

#include <cmath>
#include <iostream>
#include <map>    // <---- para calcular el nn mas comun
#include <stdexcept>
#include <utility>
#include <vector>
#include "Point.hpp"

template <size_t N, typename ElemType>
class KDTree {
 public:
  typedef std::pair<Point<N>, ElemType> value_type;

  KDTree();

  ~KDTree();

  KDTree(const KDTree &rhs);
  KDTree &operator=(const KDTree &rhs);

  size_t dimension() const;

  size_t size() const;
  bool empty() const;

  bool contains(const Point<N> &pt) const;

  void insert(const Point<N> &pt, const ElemType &value);

  ElemType &operator[](const Point<N> &pt);

  ElemType &at(const Point<N> &pt);
  const ElemType &at(const Point<N> &pt) const;

  ElemType knn_value(const Point<N> &key, size_t k) const;

  std::vector<ElemType> knn_query(const Point<N> &key, size_t k) const;

 private:
  size_t dimension_;
  size_t size_;

  struct Node {
    Node(const value_type &value);
    value_type val;
    Node *children[2];
  };

  typedef std::pair<Node*, double> knn_node;

  Node *root;

  bool find(const Point<N> &x, Node* const* &p) const;
  void knn(const Point<N> &key, Node* n, int d, std::vector<knn_node> &best, int k) const;
};

template <size_t N, typename ElemType>
KDTree<N, ElemType>::Node::Node(const value_type &value)
{
  val = value;
  children[0] = nullptr;
  children[1] = nullptr;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree() {
  size_ = 0;
  dimension_ = N;
  root = nullptr;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::~KDTree() {
  if (!empty())
  {
    std::vector<Node *> BFS;
    BFS.reserve(size());
    BFS.push_back(root);

    for (int i = 0; i < BFS.size(); i++){
      if (BFS[i]->children[0])
        BFS.push_back(BFS[i]->children[0]);
      if (BFS[i]->children[1])
        BFS.push_back(BFS[i]->children[1]);
    }

    for (int i = 1; i < BFS.size(); i++)
      delete BFS[i];
    delete root;
  }
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree(const KDTree &rhs) {
  size_ = 0;
  dimension_ = rhs.dimension();
  root = nullptr;
  if (!rhs.empty())
  {
    // almacenar un puntero a cada nodo en rhs
    std::vector<Node *> BFS;
    BFS.reserve(size());
    BFS.push_back(rhs.root);
    for (int i = 0; i < BFS.size(); i++){
      if (BFS[i]->children[0])
        BFS.push_back(BFS[i]->children[0]);
      if (BFS[i]->children[1])
        BFS.push_back(BFS[i]->children[1]);
    }

    // insertar uno a uno el valor de dichos nodos
    for (int i = 0; i < BFS.size(); i++){
      insert(BFS[i]->val.first, BFS[i]->val.second);
    }
  }
}

template <size_t N, typename ElemType>
KDTree<N, ElemType> &KDTree<N, ElemType>::operator=(const KDTree &rhs) {
  size_ = 0;
  dimension_ = rhs.dimension();
  root = nullptr;
  if (!rhs.empty())
  {
    std::vector<Node *> BFS;
    BFS.reserve(size());
    BFS.push_back(rhs.root);
    for (int i = 0; i < BFS.size(); i++){
      if (BFS[i]->children[0])
        BFS.push_back(BFS[i]->children[0]);
      if (BFS[i]->children[1])
        BFS.push_back(BFS[i]->children[1]);
    }

    for (int i = 0; i < BFS.size(); i++){
      insert(BFS[i]->val.first, BFS[i]->val.second);
    }
  }
  return *this;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::dimension() const {
  return dimension_;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::size() const {
  return size_;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::empty() const {
  return (size() == 0);
}

// Funcion auxiliar que busca un punto en el arbol y
// de encontrarlo, retorna un puntero al puntero que
// apunta al nodo que contiene dicho punto.
template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::find(const Point<N> &x, Node* const* &p) const{
  int d = 0;
  for (p = &root; *p && (*p)->val.first != x; d = (d + 1) % N)
  {
    p = &((*p)->children[(*p)->val.first[d] < x[d]]);
  }
  return *p && (*p)->val.first == x;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::contains(const Point<N> &pt) const {
  Node * const *ptr;
  bool b = find(pt,ptr);
  return b;
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::insert(const Point<N> &pt, const ElemType &value) {
  if (empty()) {
    root = new Node(std::make_pair(pt,value));
    size_++;
  }
  else
  {
    Node * const *ptr;
    if (find(pt, ptr)) {
      (*ptr)->val.second = value; // sobre-escrbir valor
    }
    else
    {
       *(const_cast<Node **>(ptr)) = new Node(std::make_pair(pt,value));
      size_++;
    }
  }
}

template <size_t N, typename ElemType>
ElemType &KDTree<N, ElemType>::operator[](const Point<N> &pt) {
  Node * const *ptr;
  if (!find(pt,ptr))
    insert(pt,ElemType());
  return (*ptr)->val.second;
}

template <size_t N, typename ElemType>
ElemType &KDTree<N, ElemType>::at(const Point<N> &pt) {
  Node * const *ptr;
  if (find(pt,ptr))
    return (*ptr)->val.second;
  else
    throw std::out_of_range("Valor no encontrado");
}

template <size_t N, typename ElemType>
const ElemType &KDTree<N, ElemType>::at(const Point<N> &pt) const {
  Node * const *ptr;
  if (find(pt,ptr))
    return (*ptr)->val.second;
  else
    throw std::out_of_range("Valor no encontrado");
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::knn(const Point<N> &key, Node* n, int d,
                              std::vector<knn_node> &best, int k) const{
  if (!n)
    return;
  
  // enqueue
  double n_dist = distance(key, n->val.first);
  if (best.size() < k) {
    best.push_back(std::make_pair(n, n_dist));
    for (int i = best.size() - 1; i > 0; i--)
      if (best[i].second < best[i - 1].second)
        std::swap(best[i], best[i - 1]);
  }
  else
  {
    if (best.back().second > n_dist){
      best[k - 1] = std::make_pair(n, n_dist);
      for (int i = best.size() - 1; i > 0; i--)
        if (best[i].second < best[i - 1].second)
          std::swap(best[i], best[i - 1]);
    }
  }

  // buscar para abajo
  int axis = d % N;
  bool right = false;
  if (n->val.first[d] < key[d]){
    right = true;
    knn(key, n->children[1], d+1, best, k);
  }
  else
  {
    knn(key, n->children[0], d+1, best, k);
  }

  // buscar en la otra rama si es necesario
  if (best.size() < k || fabs(n->val.first[axis] - key[axis]) < best.back().second){
    if (right)
      knn(key, n->children[0], d+1, best, k);
    else
      knn(key, n->children[1], d+1, best, k);
  }
}

template <size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::knn_value(const Point<N> &key, size_t k) const {
  std::vector<ElemType> elems = knn_query(key, k);
  std::map<ElemType, int> my_map;
  for(int i = 0; i < elems.size(); i++)
    my_map[elems[i]]++;

  ElemType max_element;
  int max = 0;
  typename std::map<ElemType,int>::iterator it;
  for (it = my_map.begin(); it != my_map.end(); it++){
    if (it->second > max){
      max = it->second;
      max_element = it->first;
    }
  }
  
  return max_element;
}

template <size_t N, typename ElemType>
std::vector<ElemType> KDTree<N, ElemType>::knn_query(const Point<N> &key,
                                                     size_t k) const {
  std::vector<ElemType> values;
  std::vector<knn_node> best;
  values.reserve(k);
  best.reserve(k);
  knn(key, root, 0, best, k);

  for (int i = 0; i < best.size(); i++)
    values.push_back(best[i].first->val.second);

  return values;
}

#endif  // SRC_KDTREE_HPP_
