//
// Created by klaus on 2019-08-25.
//

#include "segment.h"
#include "system.h"
#include <optional>
#include <stdexcept>

namespace mesh {

template <uint Dim, uint TopDim>
Segment<Dim, TopDim>::Segment(Mesh<Dim, TopDim> *root_mesh, ID id,
                              std::string name, bool is_hole)
    : id(id), name(std::move(name)), is_hole_(is_hole) {
  _mesh = std::make_unique<Mesh<Dim, TopDim>>(root_mesh);
}

template <uint Dim, uint TopDim>
Mesh<Dim, TopDim> *Segment<Dim, TopDim>::mesh() {
  return _mesh.get();
}

template <uint Dim, uint TopDim> MeshBase *Segment<Dim, TopDim>::mesh() const {
  return mesh();
}

template <uint Dim, uint TopDim>
ID Segment<Dim, TopDim>::getID() const noexcept {
  return id;
}

template <uint Dim, uint TopDim>
std::string Segment<Dim, TopDim>::getName() const noexcept {
  return name;
}

template <uint Dim, uint TopDim>
const std::vector<Interface<Dim, TopDim> *> &
Segment<Dim, TopDim>::interfaces() const {
  return _interfaces;
}

template <uint Dim, uint TopDim>
std::pair<Segment<Dim, TopDim> *, Segment<Dim, TopDim> *>
Interface<Dim, TopDim>::segments() const {
  return std::make_pair(seg1, seg2);
}

template <uint Dim, uint TopDim>
Segment<Dim, TopDim> *Surface<Dim, TopDim>::segment() const {
  return seg;
}

template <uint Dim, uint TopDim>
template <uint BuilderDim, uint BuilderTopDim>
Segment<Dim, TopDim>::Builder<BuilderDim, BuilderTopDim>::Builder(
    System<Dim, TopDim> *system, const std::string &name, bool is_hole)
    : system_(system), segment_(system->getOrCreateSegment(name, is_hole)) {}

template <uint Dim, uint TopDim>
template <uint BuilderDim, uint BuilderTopDim>
typename Segment<Dim, TopDim>::template Builder<BuilderDim, BuilderTopDim> &&
Segment<Dim, TopDim>::Builder<BuilderDim, BuilderTopDim>::addEdges(
    const EigenDRef<const MatrixXid> &edges) {
  if constexpr (TopDim > 0) {
    segment_->mesh()->edges().add(edges);
  } else {
    throw new std::invalid_argument("Zero dimensions don't have edges");
  }
  return std::move(*this);
}

template <uint Dim, uint TopDim>
template <uint BuilderDim, uint BuilderTopDim>
typename Segment<Dim, TopDim>::template Builder<BuilderDim, BuilderTopDim>
Segment<Dim, TopDim>::Builder<BuilderDim, BuilderTopDim>::addSegment(
    std::string name) {
  return std::move(Builder(system_, name));
}

template <uint Dim, uint TopDim>
template <uint BuilderDim, uint BuilderTopDim>
typename Segment<Dim, TopDim>::template Builder<BuilderDim, BuilderTopDim>
Segment<Dim, TopDim>::Builder<BuilderDim, BuilderTopDim>::addHole(std::string name) {
  return std::move(Builder(system_, name, true));
}

template class Segment<1, 0>;
template class Segment<1, 1>;
template class Segment<2, 0>;
template class Segment<2, 1>;
template class Segment<2, 2>;
template class Segment<3, 0>;
template class Segment<3, 1>;
template class Segment<3, 2>;
template class Segment<3, 3>;

template class Segment<1, 1>::Builder<1, 1>;
template class Segment<2, 1>::Builder<2, 1>;
template class Segment<2, 2>::Builder<2, 2>;
template class Segment<3, 1>::Builder<3, 1>;
template class Segment<3, 2>::Builder<3, 2>;
template class Segment<3, 3>::Builder<3, 3>;

template class Interface<1, 1>;
template class Interface<2, 1>;
template class Interface<3, 1>;
template class Interface<2, 2>;
template class Interface<3, 2>;
template class Interface<3, 3>;

template class Hole<1, 1>;
template class Hole<2, 1>;
template class Hole<3, 1>;
template class Hole<2, 2>;
template class Hole<3, 2>;
template class Hole<3, 3>;

} // namespace mesh
