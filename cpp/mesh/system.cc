//
// Created by klaus on 2019-11-11.
//

#include "system.h"
#include "expression.hh"
#include "segment.h"
#include <memory>

#define VOID void
#define REAL double

#include "triangle.h"

namespace mesh {

template <uint Dim, uint TopDim> System<Dim, TopDim>::System() {
  _mesh = std::make_unique<Mesh<Dim, TopDim>>();
  _voronoi = std::make_unique<Mesh<Dim, TopDim>>();
}

template <uint Dim, uint TopDim>
System<Dim, TopDim> &System<Dim, TopDim>::operator=(System &&sys) noexcept {
  segment_names = std::move(sys.segment_names);
  segments = std::move(sys.segments);
  inthash2idx = std::move(sys.inthash2idx);
  interfaces = std::move(sys.interfaces);
  _mesh = std::move(sys._mesh);
  _voronoi = std::move(sys._voronoi);
  // attributes = std::move(sys.attributes);
  return *this;
}

template <uint Dim, uint TopDim> MeshBase *System<Dim, TopDim>::mesh() const {
  return _mesh.get();
}

template <uint Dim, uint TopDim>
Mesh<Dim, TopDim> *System<Dim, TopDim>::mesh() {
  return _mesh.get();
}

template <uint Dim, uint TopDim>
Mesh<Dim, TopDim> *System<Dim, TopDim>::voronoi() {
  return _voronoi.get();
}

template <uint Dim, uint TopDim>
Segment<Dim, TopDim> *
System<Dim, TopDim>::getOrCreateSegment(const std::string &name, bool is_hole) {
  const auto it = find(segment_names.begin(), segment_names.end(), name);
  if (it != segment_names.end()) {
    return segments[distance(segment_names.begin(), it)].get();
  }
  segment_names.push_back(name);
  surfaces.emplace_back(
      std::make_unique<Surface<Dim, TopDim>>(mesh(), surfaces.size(), name));
  segments.emplace_back(std::make_unique<Segment<Dim, TopDim>>(
      mesh(), segments.size() + interfaces.size(), name, is_hole));
  segments.back()->surface_ = surfaces.back().get();
  return segments.back().get();
}

template <uint Dim, uint TopDim>
Segment<Dim, TopDim> *System<Dim, TopDim>::segment(const std::string &name) {
  const auto it = find(segment_names.begin(), segment_names.end(), name);
  if (it != segment_names.end())
    return segments[distance(segment_names.begin(), it)].get();
  return nullptr;
}

template <uint Dim, uint TopDim>
SegmentBase *System<Dim, TopDim>::segment(const std::string &name) const {
  return static_cast<const SystemBase *>(this)->segment(name);
}

template <uint Dim, uint TopDim>
Segment<Dim, TopDim> *System<Dim, TopDim>::segment(ID id) {
  if (id < segments.size())
    return segments[id].get();
  return nullptr;
}

template <uint Dim, uint TopDim>
SegmentBase *System<Dim, TopDim>::segment(ID id) const {
  return static_cast<const SystemBase *>(this)->segment(id);
}

template <uint Dim, uint TopDim>
Interface<Dim, TopDim> *
System<Dim, TopDim>::interface(const std::string &seg1_name,
                               const std::string &seg2_name) {
  Segment<Dim, TopDim> *seg1 = segment(seg1_name);
  if (seg1 == nullptr)
    return nullptr;
  Segment<Dim, TopDim> *seg2 = segment(seg2_name);
  if (seg2 == nullptr)
    return nullptr;
  return interface(seg1->getID(), seg2->getID());
}

template <uint Dim, uint TopDim>
SegmentBase *
System<Dim, TopDim>::interface(const std::string &seg1_name,
                               const std::string &seg2_name) const {
  return static_cast<const SystemBase *>(this)->interface(seg1_name, seg2_name);
}

template <uint Dim, uint TopDim>
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> System<Dim, TopDim>::laplace(
    const EigenDRef<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> &val)
    const {
  if (val.rows() != laplace_op_.cols() || val.rows() != laplace_op_.rows()) {
    throw std::invalid_argument("Input dimension invalid");
  }
  return laplace_op_ * val;
}

template <uint Dim, uint TopDim>
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> System<Dim, TopDim>::gradient(
    const EigenDRef<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> &val)
    const {
  if (val.rows() != laplace_op_.cols() || val.rows() != laplace_op_.rows()) {
    throw std::invalid_argument("Input dimension invalid");
  }
  return gradient_op_ * val;
}

template <uint Dim, uint TopDim>
Eigen::SparseMatrix<double, Eigen::RowMajor>
System<Dim, TopDim>::getLaplacianMatrix() const {
  return laplace_op_;
}

template <uint Dim, uint TopDim>
Eigen::SparseMatrix<double, Eigen::RowMajor>
System<Dim, TopDim>::getGradientMatrix() const {
  return gradient_op_;
}

template <uint Dim, uint TopDim>
Eigen::SparseMatrix<double, Eigen::RowMajor>
System<Dim, TopDim>::getMassMatrix() const {
  return mass_op_;
}

template <uint Dim, uint TopDim>
void System<Dim, TopDim>::assembleMassOperator() {
  const size_t nvertices = _mesh->vertices().size();
  mass_op_.resize(nvertices, nvertices);
  std::vector<Eigen::Triplet<double>> tripletList;
  for (size_t ib = 0; ib < _mesh->bodies().size(); ib++) {
    const mesh::MeshElement *body = _mesh->bodies()[ib];
    const Eigen::Vector3d v0 = body->getPoint(0);
    const Eigen::Vector3d v1 = body->getPoint(1);
    const Eigen::Vector3d v2 = body->getPoint(2);
    const Eigen::Vector3d e0 = v1 - v0;
    const Eigen::Vector3d e1 = v2 - v1;
    const double area = 0.5 * (e0.cross(e1)).norm();
    // Each vertex of the triangle gets 1/3 of the area
    for (size_t i = 0; i < Dim; ++i) {
      tripletList.push_back(
          Eigen::Triplet<double>((*body)[i], (*body)[i], area / 3.0));
    }
  }
  mass_op_.setFromTriplets(tripletList.begin(), tripletList.end());
}

template <uint Dim, uint TopDim>
void System<Dim, TopDim>::assembleLaplaceOperator() {
  const size_t numVertices = _mesh->vertices().size();
  laplace_op_.resize(numVertices,
                     numVertices); // Laplacian is a square matrix of size
  std::vector<Eigen::Triplet<double>> tripletList;
  for (size_t ib = 0; ib < _mesh->bodies().size(); ib++) {
    mesh::MeshElement *body = _mesh->bodies()[ib];
    const Eigen::Vector3d v0 = body->getPoint(0);
    const Eigen::Vector3d v1 = body->getPoint(1);
    const Eigen::Vector3d v2 = body->getPoint(2);
    const Eigen::Vector3d e0 = v1 - v0;
    const Eigen::Vector3d e1 = v2 - v1;
    const Eigen::Vector3d e2 = v0 - v2;
    // Compute area of the triangle
    const Eigen::Vector3d normal = e0.cross(e1);
    const double area = 0.5 * normal.norm();
    // Compute cotangents of the angles for each vertex
    const double cotan0 = e0.dot(e2) / normal.norm();
    const double cotan1 = e1.dot(e0) / normal.norm();
    const double cotan2 = e2.dot(e1) / normal.norm();
    // Add contributions to the Laplacian matrix
    tripletList.push_back(
        Eigen::Triplet<double>((*body)[0], (*body)[1], cotan2 / (2 * area)));
    tripletList.push_back(
        Eigen::Triplet<double>((*body)[1], (*body)[0], cotan2 / (2 * area)));
    tripletList.push_back(
        Eigen::Triplet<double>((*body)[1], (*body)[2], cotan0 / (2 * area)));
    tripletList.push_back(
        Eigen::Triplet<double>((*body)[2], (*body)[1], cotan0 / (2 * area)));
    tripletList.push_back(
        Eigen::Triplet<double>((*body)[2], (*body)[0], cotan1 / (2 * area)));
    tripletList.push_back(
        Eigen::Triplet<double>((*body)[0], (*body)[2], cotan1 / (2 * area)));
    // Diagonal entries (sf cotaEigen::Triplet<double>(ngents for each vertex)
    tripletList.push_back(Eigen::Triplet<double>(
        (*body)[0], (*body)[0], -(cotan1 + cotan2) / (2 * area)));
    tripletList.push_back(Eigen::Triplet<double>(
        (*body)[1], (*body)[1], -(cotan0 + cotan2) / (2 * area)));
    tripletList.push_back(Eigen::Triplet<double>(
        (*body)[2], (*body)[2], -(cotan0 + cotan1) / (2 * area)));
  }
  laplace_op_.setFromTriplets(tripletList.begin(), tripletList.end());
}

template <uint Dim, uint TopDim>
void System<Dim, TopDim>::assembleGradientOperator() {
  const size_t nvertices = _mesh->vertices().size();
  const size_t nbodies = _mesh->bodies().size();
  // Gradient matrix: maps scalar field on vertices to vector field on edges
  gradient_op_.resize(nvertices * 3,
                      nbodies); // 3 rows per triangle (for x, y, z components)
  std::vector<Eigen::Triplet<double>> tripletList;

  for (size_t ib = 0; ib < nbodies; ++ib) {
    mesh::MeshElement *body = _mesh->bodies()[ib];
    const Eigen::Vector3d v0 = body->getPoint(0);
    const Eigen::Vector3d v1 = body->getPoint(1);
    const Eigen::Vector3d v2 = body->getPoint(2);
    // Compute the edge vectors
    const Eigen::Vector3d e0 = v1 - v0;
    const Eigen::Vector3d e1 = v2 - v1;
    const Eigen::Vector3d e2 = v0 - v2;
    // Compute the area of the triangle (half the cross product magnitude)
    Eigen::Vector3d normal = e0.cross(e1);
    double area = 0.5 * normal.norm();
    for (size_t i = 0; i < 3; ++i) {
      tripletList.push_back(
          Eigen::Triplet<double>(3 * ib, (*body)[i], (e0[i] / (2 * area))));
      tripletList.push_back(
          Eigen::Triplet<double>(3 * ib + 1, (*body)[i], (e1[i] / (2 * area))));
      tripletList.push_back(
          Eigen::Triplet<double>(3 * ib + 2, (*body)[i], (e2[i] / (2 * area))));
    }
  }
  gradient_op_.setFromTriplets(tripletList.begin(), tripletList.end());
}

template <uint Dim, uint TopDim>
Interface<Dim, TopDim> *System<Dim, TopDim>::interface(ID seg1_id, ID seg2_id) {
  SimplexHash<1> hash;
  const ullong key = hash(seg1_id, seg2_id);
  const auto it = inthash2idx.find(key);
  if (it != inthash2idx.end())
    return interfaces[it->second].get();
  const size_t idx = interfaces.size();
  inthash2idx[key] = idx;
  interfaces.emplace_back(std::make_unique<Interface<Dim, TopDim>>(
      mesh(), idx, segment_names[seg1_id] + "_" + segment_names[seg2_id]));
  Interface<Dim, TopDim> *intf = interfaces.back().get();
  for (size_t iseg = 0; iseg < segments.size(); ++iseg) {
    if (segments[iseg]->getID() == seg1_id ||
        segments[iseg]->getID() == seg2_id)
      segments[iseg]->_interfaces.push_back(intf);
  }
  return interfaces.back().get();
}

template <uint Dim, uint TopDim>
SegmentBase *System<Dim, TopDim>::interface(ID seg1_id, ID seg2_id) const {
  return static_cast<const SystemBase *>(this)->interface(seg1_id, seg2_id);
}

template <uint Dim, uint TopDim>
Surface<Dim, TopDim> *System<Dim, TopDim>::surface(ID id) {
  if (id < surfaces.size()) {
    return surfaces[id].get();
  }
  return nullptr;
}

template <uint Dim, uint TopDim>
Surface<Dim, TopDim> *System<Dim, TopDim>::surface(const std::string &name) {
  Segment<Dim, TopDim> *seg = segment(name);
  if (seg != nullptr) {
    return surfaces[seg->getID()].get();
  }
  return nullptr;
}

template <uint Dim, uint TopDim>
SegmentBase *System<Dim, TopDim>::surface(ID id) const {
  return static_cast<const SystemBase *>(this)->surface(id);
}

template <uint Dim, uint TopDim>
SegmentBase *System<Dim, TopDim>::surface(const std::string &name) const {
  return static_cast<const SystemBase *>(this)->surface(name);
}

template <uint Dim, uint TopDim>
System<Dim, TopDim>::Builder::Builder() : vertices_proxy_(system_.mesh()) {}

template <uint Dim, uint TopDim>
typename System<Dim, TopDim>::Builder &
System<Dim, TopDim>::Builder::addVertices(
    const EigenDRef<const Eigen::MatrixXd> &points) {
  vertices_proxy_.add(points);
  return *this;
}

template <uint Dim, uint TopDim>
typename Segment<Dim, TopDim>::template Builder<Dim, TopDim>
System<Dim, TopDim>::Builder::addSegment(std::string name) {
  return typename Segment<Dim, TopDim>::template Builder<Dim, TopDim>(&system_,
                                                                      name);
}

template <uint Dim, uint TopDim>
System<Dim, TopDim>::Builder::operator System<Dim, TopDim> &&() {
  return std::move(system_);
}

template <uint Dim, uint TopDim>
Mesh<Dim, TopDim> *System<Dim, TopDim>::Builder::getMesh() {
  return system_.mesh();
}

template <uint Dim, uint TopDim>
std::unique_ptr<System<Dim, TopDim>>
System<Dim, TopDim>::Builder::create(double area, double min_angle) {
  throw std::runtime_error("Not yet implemented");
}

template class System<1, 1>;
template class System<2, 1>;
template class System<2, 2>;
template class System<3, 1>;
template class System<3, 2>;
template class System<3, 3>;

} // namespace mesh
