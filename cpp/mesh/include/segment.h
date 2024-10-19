//
// Created by klaus on 2019-08-25.
//

#ifndef ULB_SEGMENT_H
#define ULB_SEGMENT_H

#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "attribute.h"
#include "mesh.h"

namespace mesh {

class SegmentBase {
public:
  virtual MeshBase *mesh() const = 0;

  virtual ID getID() const noexcept = 0;

  virtual std::string getName() const noexcept = 0;
};

template <uint Dim, uint TopDim> class Interface;

template <uint Dim, uint TopDim> class Surface;

template <uint Dim, uint TopDim> class System;

template <uint Dim, uint TopDim> class Mesh;

template <uint Dim, uint TopDim> class Segment : public SegmentBase {
  friend System<Dim, TopDim>;

public:
  explicit Segment(Mesh<Dim, TopDim> *system, ID id, std::string name,
                   bool isHole);

  Segment() = delete;

  Segment(Segment &) = delete;

  Segment(const Segment &) = delete;

  Segment &operator=(Segment &) = delete;

  Segment &operator=(const Segment &) = delete;

  Segment(Segment &&seg) noexcept = default;

  Segment &operator=(Segment &&seg) noexcept = default;

  Mesh<Dim, TopDim> *mesh();

  MeshBase *mesh() const override;

  Surface<Dim, TopDim> &surface() const;

  bool isHole() const { return is_hole_; }

  const std::vector<Interface<Dim, TopDim> *> &interfaces() const;

  ID getID() const noexcept override;

  std::string getName() const noexcept override;

  template <typename T>
  Attribute<T> *getOrCreateAttribute(StorageLocation location,
                                     std::string name) {
    const auto it = attributes_.find(name);
    if (it != attributes_.end()) {
      return dynamic_cast<Attribute<T> *>(it->second.get());
    }
    if (location == StorageLocation::VERTEX) {
      return static_cast<Attribute<T> *>(
          attributes_
              .emplace(name, std::make_unique<Attribute<T>>(
                                 mesh()->vertices().size()))
              .first->second.get());
    } else {
      throw std::invalid_argument("Unsupported storage location");
    }
  }

  template <uint BuilderDim, uint BuilderTopDim> class Builder {
  public:
    Builder(System<Dim, TopDim> *system, const std::string &name,
            bool is_hole = false);

    Builder &&addEdges(const EigenDRef<const MatrixXid> &indices);

    Builder addSegment(std::string name);

    Builder addHole(std::string name);

  private:
    System<Dim, TopDim> *system_;
    Segment *segment_;
  };

private:
  std::unordered_map<std::string, std::unique_ptr<AttributeBase>> attributes_;
  std::unique_ptr<Mesh<Dim, TopDim>> _mesh;
  std::vector<Interface<Dim, TopDim> *> _interfaces;
  Surface<Dim, TopDim> *surface_;
  bool is_hole_;
  ID id;
  std::string name;
};

template <uint Dim, uint TopDim> class Hole : public Segment<Dim, TopDim> {
public:
  using Segment<Dim, TopDim>::Segment;

  explicit Hole(Mesh<Dim, TopDim> *system, ID id, std::string name)
      : Segment<Dim, TopDim>(system, id, name, true) {}

  Hole() = delete;

  Hole &operator=(Hole &) = delete;

  Hole &operator=(const Hole &) = delete;

  Hole &operator=(Hole &&seg) noexcept = default;
};

template <uint Dim, uint TopDim>
class Surface : public Segment<Dim, TopDim - 1> {
public:
  using Segment<Dim, TopDim - 1>::Segment;

  explicit Surface(Mesh<Dim, TopDim> *system, ID id, std::string name)
      : Segment<Dim, TopDim - 1>(system, id, name, false) {}

  Surface() = delete;

  Surface &operator=(Surface &) = delete;

  Surface &operator=(const Surface &) = delete;

  Surface &operator=(Surface &&seg) noexcept = default;

  Segment<Dim, TopDim> *segment() const;

private:
  Segment<Dim, TopDim> *seg;
};

template <uint Dim, uint TopDim>
class Interface : public Segment<Dim, TopDim - 1> {
public:
  using Segment<Dim, TopDim - 1>::Segment;

  explicit Interface(Mesh<Dim, TopDim> *system, ID id, std::string name)
      : Segment<Dim, TopDim - 1>(system, id, name, false) {}

  Interface() = delete;

  Interface &operator=(Interface &) = delete;

  Interface &operator=(const Interface &) = delete;

  Interface &operator=(Interface &&seg) noexcept = default;

  std::pair<Segment<Dim, TopDim> *, Segment<Dim, TopDim> *> segments() const;

private:
  Segment<Dim, TopDim> *seg1;
  Segment<Dim, TopDim> *seg2;
};

} // namespace mesh

#endif // ULB_SEGMENT_H
