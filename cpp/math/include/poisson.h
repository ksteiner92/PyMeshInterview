//
// Created by klaus on 02.03.19.
//

#ifndef ULB_POISSON_H
#define ULB_POISSON_H

#include <vector>

#include "eigen.h"
#include "segment.h"
#include "system.h"

namespace math {

template <uint Dim, uint TopDim> class Poisson {
  static_assert(Dim >= TopDim && TopDim >= 1,
                "Dimension combination not supported");

public:
  Poisson(mesh::System<Dim, TopDim> *system);

  void solve(const std::string &phi, const std::string &rho,
             const std::vector<ID> &dirichlet);

  void grad(const std::string &u, const std::string &E) const;

private:
  mesh::System<Dim, TopDim> *system_;
  Eigen::SparseMatrix<double, Eigen::RowMajor> K;
  std::vector<long long int> interior;
};

} // namespace math

#endif // ULB_POISSON_H
