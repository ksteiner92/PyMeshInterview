//
// Created by klaus on 19.01.19.
//

#ifndef LBM_ATTRIBUTES_H
#define LBM_ATTRIBUTES_H

#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <valarray>
#include <vector>

// #include "expression.hh"

namespace mesh {

class MeshBase;

class SystemBase;

class SegmentBase;

enum StorageLocation { VERTEX, EDGE, FACE, CELL };

class AttributeBase {
public:
  virtual std::size_t size() const noexcept = 0;
};

template <typename T> class Attribute : public AttributeBase {
public:
  Attribute(size_t size) : values(size) {}

  T &get(size_t loc) { return values[loc]; }

  T &operator[](size_t loc) { return values[loc]; }

  void set(size_t loc, const T &value) { values[loc] = value; }

  std::size_t size() const noexcept override { return values.size(); }

  //  Attribute &operator[](std::size_t dim) {
  //    if (idim >= iextents.size()) {
  //      std::stringstream ss;
  //      ss << "Extent " << idim << " does not exist";
  //      throw std::range_error(ss.str());
  //    }
  //    iextents[idim] = dim;
  //    ++idim;
  //    return *this;
  //  }
  //
  //  Attribute &operator()(const SelectorBase<Location> &selector) {
  //    mask &= selector(mesh);
  //    return *this;
  //  }
  //
  //  Attribute &operator=(const T &value) {
  //    values[mask] = value;
  //    return *this;
  //  }
  //
  //  std::string getName() const { return name; }
  //
  //  AttributeExtent &getExtents() override { return extents; }
  //
  //  const AttributeExtent &getExtents() const override { return extents; }

private:
  std::valarray<T> values;
};

} // namespace mesh

#endif // LBM_ATTRIBUTES_H
