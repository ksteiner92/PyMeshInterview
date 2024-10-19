
#include "diffoperators.h"

namespace math {

template <uint Dim, uint TopDim>
DifferentialOperators<Dim, TopDim>::DifferentialOperators(
    mesh::System<Dim, TopDim> *system)
    : system_(system) {
  mesh::Mesh<Dim, TopDim> *mesh = system_->mesh();
  const size_t nvertices = mesh->vertices().size();
  M.resize(nvertices, nvertices);
  std::vector<Eigen::Triplet<double>> tripletList;
  for (size_t ib = 0; ib < mesh->bodies().size(); ib++) {
    const mesh::MeshElement *body = mesh->bodies()[ib];
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
  M.setFromTriplets(tripletList.begin(), tripletList.end());
}

template <uint Dim, uint TopDim>
Eigen::MatrixXd DifferentialOperators<Dim, TopDim>::laplace(
    const EigenDRef<const Eigen::MatrixXd> &f) const {
  const size_t numVertices = system_->mesh()->vertices().size();
  Eigen::SparseMatrix<double, Eigen::RowMajor> L(
      numVertices, numVertices); // Laplacian is a square matrix of size
  std::vector<Eigen::Triplet<double>> tripletList;

  mesh::Mesh<Dim, TopDim> *mesh = system_->mesh();
  for (size_t ib = 0; ib < mesh->bodies().size(); ib++) {
    mesh::MeshElement *body = mesh->bodies()[ib];
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
  L.setFromTriplets(tripletList.begin(), tripletList.end());
  return std::move(L * f);
}

template class DifferentialOperators<1, 1>;
template class DifferentialOperators<2, 1>;
template class DifferentialOperators<2, 2>;
template class DifferentialOperators<3, 1>;
template class DifferentialOperators<3, 2>;
template class DifferentialOperators<3, 3>;

} // namespace math
