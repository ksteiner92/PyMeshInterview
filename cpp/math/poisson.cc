//
// Created by klaus on 02.03.19.
//

#include <array>
#include <limits>
#include <vector>

#include "attribute.h"
#include "logger.h"
#include "poisson.h"

namespace math {

template <uint Dim, uint TopDim>
Poisson<Dim, TopDim>::Poisson(mesh::System<Dim, TopDim> *system)
    : system_(system) {
  LOG_T(INFO) << "[Poisson] Pre-processing ..." << io::LogFlags::ENDL;
  io::ProgressBar progress;
  mesh::Mesh<Dim, TopDim> *mesh = system->mesh();
  const size_t nvertices = mesh->vertices().size();
  K.resize(nvertices, nvertices);
  std::vector<Eigen::Triplet<double>> kij;
  progress.start(mesh->bodies().size());
  ID max_k = 0;
  ID max_l = 0;
  ID min_k = std::numeric_limits<ID>::max();
  ID min_l = std::numeric_limits<ID>::max();
  for (size_t ic = 0; ic < mesh->bodies().size(); ic++) {
    mesh::MeshElement *body = mesh->bodies()[ic];
    std::vector<Eigen::Matrix<double, Dim, 1>> p(body->getNumVertices());
    for (uint8_t i = 0; i < body->getNumVertices(); i++) {
      p[i] = body->getPoint(i);
    }
    std::array<double, 3> b;
    std::array<double, 3> c;
    for (uint8_t i = 0; i < 3; i++) {
      b[i] = p[(i + 1) % 3](1) - p[(i + 2) % 3](1);
      c[i] = p[(i + 2) % 3](0) - p[(i + 1) % 3](0);
    }
    const double Omega = 2 * std::fabs(p[1](0) * p[2](1) + p[0](0) * p[1](1) +
                                       p[0](1) * p[2](0) - p[1](0) * p[0](1) -
                                       p[2](0) * p[1](1) - p[2](1) * p[0](0));
    for (uint8_t i = 0; i < 3; i++) {
      for (uint8_t j = i; j < 3; j++) {
        const double Kij = (b[j] * b[i] + c[j] * c[i]) / Omega;
        const ID k = std::min((*body)[i], (*body)[j]);
        const ID l = std::max((*body)[i], (*body)[j]);
        max_k = std::max(k, max_k);
        max_l = std::max(l, max_l);
        min_k = std::min(k, min_k);
        min_l = std::min(l, min_l);
        kij.push_back(Eigen::Triplet<double>(k, l, Kij));
      }
    }
    progress.update(ic);
  }
  std::cout << "nertices: " << nvertices << std::endl;
  std::cout << "kij: " << kij.size() << std::endl;
  std::cout << "max_k: " << max_k << std::endl;
  std::cout << "max_l: " << max_l << std::endl;
  std::cout << "min_k: " << min_k << std::endl;
  std::cout << "min_l: " << min_l << std::endl;
  K.setFromTriplets(kij.begin(), kij.end());
  progress.stop();
}

template <uint Dim, uint TopDim>
void Poisson<Dim, TopDim>::grad(const std::string &u,
                                const std::string &w) const {
  mesh::Mesh<Dim, TopDim> *mesh = system_->mesh();
  auto phi = *(system_->template getOrCreateAttribute<double>(
      mesh::StorageLocation::VERTEX, u));
  // if (phi == nullptr) {
  //   std::stringstream ss;
  //   ss << "Attribute '" << u << "' does not exist on vertices";
  //   throw std::invalid_argument(ss.str());
  // }
  auto E = *reinterpret_cast<mesh::Attribute<Eigen::Matrix<double, Dim, 1>> *>(
      system_->template getOrCreateAttribute<Eigen::Matrix<double, Dim, 1>>(
          mesh::StorageLocation::FACE, w));
  // if (E == nullptr) {
  //   std::stringstream ss;
  //   ss << "Attribute '" << w << "' does not exist on mesh bodies";
  //   throw std::invalid_argument(ss.str());
  // }
  LOG_T(INFO) << "[Poisson] Calculating gradient..." << io::LogFlags::ENDL;
  io::ProgressBar progress;
  progress.start(mesh->bodies().size());
  // vector<double> vvol(mesh->getNumVertices(), 0.0);
  // vector<Matrix<double, Dim, 1>> gradient(mesh1D->getNumVertices(),
  // Matrix<double, Dim, 1>::Zero());
  for (size_t ib = 0; ib < mesh->bodies().size(); ib++) {
    mesh::MeshElement *body = mesh->bodies()[ib];
    const size_t nvertices = body->getNumVertices();
    std::vector<Eigen::Matrix<double, Dim, 1>> p(nvertices);
    for (uint8_t i = 0; i < nvertices; i++) {
      p[i] = body->getPoint(i);
    }
    const double Omega =
        0.5 * (p[1](0) * p[2](1) + p[0](0) * p[1](1) + p[0](1) * p[2](0) -
               p[1](0) * p[0](1) - p[2](0) * p[1](1) - p[2](1) * p[0](0));
    E[ib] = Eigen::Matrix<double, Dim, 1>::Zero();
    for (int i = 0; i < nvertices; i++) {
      E[ib](0) += (p[(i + 1) % 3](1) - p[(i + 2) % 3](1)) * phi[(*body)[i]] /
                  (2.0 * Omega);
      E[ib](1) += (p[(i + 2) % 3](0) - p[(i + 1) % 3](0)) * phi[(*body)[i]] /
                  (2.0 * Omega);
    }
    progress.update(ib);
  }
  progress.stop();
}

template <uint Dim, uint TopDim>
void Poisson<Dim, TopDim>::solve(const std::string &phi_str,
                                 const std::string &rho_str,
                                 const std::vector<ID> &tmp) {

  // auto phi_attr = system_->template getOrCreateAttribute<double>(
  //     mesh::StorageLocation::VERTEX, phi_str);
  // Eigen::SparseMatrix<double, Eigen::RowMajor> L =
  //     system_->getLaplacianMatrix();
  // Eigen::SparseMatrix<double, Eigen::RowMajor> M = system_->getMassMatrix();

  // LOG_T(INFO) << "[Poisson] Assembling ..." << io::LogFlags::ENDL;
  // mesh::Mesh<Dim, TopDim> *mesh = system_->mesh();
  // const size_t nvertices = mesh->vertices().size();
  // Eigen::VectorXd F = Eigen::VectorXd::Zero(nvertices);
  // for (int idx : tmp) {
  //   for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(L,
  //   idx);
  //        it; ++it) {
  //     it.valueRef() = it.row() == idx ? 1.0 : 0.0;
  //   }
  //   for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(M,
  //   idx);
  //        it; ++it) {
  //     it.valueRef() = it.row() == idx ? 1.0 : 0.0;
  //   }
  // }
  // LOG_T(INFO) << "[Poisson] Solving ..." << io::LogFlags::ENDL;
  // Eigen::ConjugateGradient<Eigen::SparseMatrix<double, Eigen::RowMajor>,
  //                          Eigen::Upper>
  //     solver;
  // Eigen::VectorXd phi = solver.compute(L).solve(F);
  // if (solver.info() != Eigen::Success) {
  //   throw std::runtime_error("Could not solve Poisson system");
  // }
  // for (size_t i = 0; i < nvertices; i++) {
  //   (*phi_attr)[i] = phi(i);
  // }
  LOG_T(INFO) << "[Poisson] Assembling ..." << io::LogFlags::ENDL;
  // auto rho_attr = mesh1D->template
  // getOrCreateAttributeOnVertex<double>(rho_str);

  mesh::Mesh<Dim, TopDim> *mesh = system_->mesh();
  auto phi_attr = system_->template getOrCreateAttribute<double>(
      mesh::StorageLocation::VERTEX, phi_str);
  std::vector<ID> dirichlet(tmp.size());
  std::copy(tmp.begin(), tmp.end(), dirichlet.begin());
  std::sort(dirichlet.begin(), dirichlet.end());
  const size_t nvertices = mesh->vertices().size();
  Eigen::VectorXd F = Eigen::VectorXd::Zero(nvertices);
  Eigen::SparseMatrix<double, Eigen::RowMajor> A = K;
  for (size_t i = 0; i < dirichlet.size(); i++) {
    const ID dv = dirichlet[i];
    for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(A, dv);
         it; ++it) {
      it.valueRef() = it.col() == dv ? 1.0 : 0.0;
    }
    const double phi = phi_attr->get(dv);
    F(dv) = phi;
    const auto eov =
        static_cast<typename mesh::Mesh<Dim, 1>::EdgesProxy &>(mesh->edges())
            .edgesOfVertex(dv);
    for (const auto edge : *eov) {
      const ID neighbor = (*edge)[0] == dv ? (*edge)[1] : (*edge)[0];
      const ID k = std::min(neighbor, dv);
      const ID l = std::max(neighbor, dv);
      A.coeffRef(k, l) = 0.0;
      if (!std::binary_search(dirichlet.begin(), dirichlet.end(), neighbor)) {
        F(neighbor) -= K.coeff(k, l) * phi;
      }
    }
  }
  LOG_T(INFO) << "[Poisson] Solving ..." << io::LogFlags::ENDL;
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double, Eigen::RowMajor>,
                           Eigen::Upper>
      solver;
  Eigen::VectorXd phi = solver.compute(A).solve(F);
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("Could not solve Poisson system");
  }
  for (size_t i = 0; i < nvertices; i++) {
    (*phi_attr)[i] = phi(i);
  }
}

template class Poisson<1, 1>;
template class Poisson<2, 1>;
template class Poisson<2, 2>;
template class Poisson<3, 1>;
template class Poisson<3, 2>;
template class Poisson<3, 3>;

} // namespace math
