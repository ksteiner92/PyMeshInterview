//
// Created by klaus on 2019-11-11.
//

#ifndef PYULB_SYSTEM_H
#define PYULB_SYSTEM_H

#include "attribute.h"
#include "expression.hh"
#include "mesh.h"
#include "segment.h"
#include <memory>

namespace mesh {

class SegmentBase;

template <uint Dim, uint TopDim> class Segment;

template <uint Dim, uint TopDim> class Interface;

class SystemBase {
public:
  virtual MeshBase *mesh() const = 0;

  virtual SegmentBase *segment(const std::string &name) const = 0;

  virtual SegmentBase *segment(ID id) const = 0;

  virtual SegmentBase *interface(ID seg1_id, ID seg2_id) const = 0;

  virtual SegmentBase *interface(const std::string &seg1,
                                 const std::string &seg2) const = 0;

  virtual SegmentBase *surface(ID seg) const = 0;

  virtual SegmentBase *surface(const std::string &seg) const = 0;
};

template <uint Dim, uint TopDim = Dim> class System : public SystemBase {
  friend Segment<Dim, TopDim>;
  friend std::unique_ptr<System<Dim, TopDim>>
  std::make_unique<System<Dim, TopDim>>();

private:
  System();

public:
  class Builder {
  public:
    Builder();

    Builder &addVertices(const EigenDRef<const Eigen::MatrixXd> &points);

    typename Segment<Dim, TopDim>::template Builder<Dim, TopDim>
    addSegment(std::string name);

    operator System &&();

    Mesh<Dim, TopDim> *getMesh();

    std::unique_ptr<System> create(double area = 0.0, double min_angle = 0.0);

  private:
    System system_;
    typename Mesh<Dim, TopDim>::VerticesProxy vertices_proxy_;
  };

  System(System &&) noexcept = default;

  System(System &) = delete;

  System(const System &) = delete;

  System &operator=(System &&) noexcept;

  System &operator=(System &) = delete;

  System &operator=(const System &) = delete;

  Mesh<Dim, TopDim> *mesh();

  MeshBase *mesh() const override;

  Mesh<Dim, TopDim> *voronoi();

  Segment<Dim, TopDim> *segment(const std::string &name);

  Segment<Dim, TopDim> *segment(ID id);

  SegmentBase *segment(const std::string &name) const override;

  SegmentBase *segment(ID id) const override;

  Surface<Dim, TopDim> *surface(ID seg);

  Surface<Dim, TopDim> *surface(const std::string &seg);

  SegmentBase *surface(ID seg) const override;

  SegmentBase *surface(const std::string &seg) const override;

  Interface<Dim, TopDim> *interface(ID seg1_id, ID seg2_id);

  Interface<Dim, TopDim> *interface(const std::string &seg1,
                                    const std::string &seg2);

  SegmentBase *interface(ID seg1_id, ID seg2_id) const override;

  SegmentBase *interface(const std::string &seg1,
                         const std::string &seg2) const override;

  template <typename T>
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> laplace(
      const EigenDRef<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> &val)
      const;

  template <typename T>
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> gradient(
      const EigenDRef<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> &val)
      const;

  Eigen::SparseMatrix<double, Eigen::RowMajor> getLaplacianMatrix() const;

  Eigen::SparseMatrix<double, Eigen::RowMajor> getGradientMatrix() const;

  Eigen::SparseMatrix<double, Eigen::RowMajor> getMassMatrix() const;

  static Builder create() { return Builder(); }

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

private:
  std::unordered_map<std::string, std::unique_ptr<AttributeBase>> attributes_;
  std::vector<std::string> segment_names;
  std::vector<std::unique_ptr<Segment<Dim, TopDim>>> segments;
  std::unordered_map<ullong, size_t> inthash2idx;
  std::vector<std::unique_ptr<Interface<Dim, TopDim>>> interfaces;
  std::vector<std::unique_ptr<Surface<Dim, TopDim>>> surfaces;
  std::unique_ptr<Mesh<Dim, TopDim>> _mesh;
  std::unique_ptr<Mesh<Dim, TopDim>> _voronoi;
  Eigen::SparseMatrix<double, Eigen::RowMajor> laplace_op_;
  Eigen::SparseMatrix<double, Eigen::RowMajor> gradient_op_;
  Eigen::SparseMatrix<double, Eigen::RowMajor> mass_op_;
  // std::unordered_map<std::string, std::unique_ptr<AttributeBase>> attributes;

  Segment<Dim, TopDim> *getOrCreateSegment(const std::string &name, bool isHole = false);

  void assembleMassOperator();

  void assembleLaplaceOperator();

  void assembleGradientOperator();
};

} // namespace mesh

#endif // PYULB_SYSTEM_H
