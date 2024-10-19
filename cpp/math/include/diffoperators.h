#ifndef PYMESH_DIFFOPERATORS_H
#define PYMESH_DIFFOPERATORS_H

#include <vector>

#include "system.h"

namespace math {

template <uint Dim, uint TopDim> class DifferentialOperators {
  static_assert(Dim >= TopDim && TopDim >= 1,
                "Dimension combination not supported");

public:
  DifferentialOperators(mesh::System<Dim, TopDim> *system);

  Eigen::MatrixXd laplace(const EigenDRef<const Eigen::MatrixXd> &f) const;
private:
  mesh::System<Dim, TopDim> *system_;
  Eigen::SparseMatrix<double, Eigen::RowMajor> K;
  Eigen::SparseMatrix<double, Eigen::RowMajor> M;
  std::vector<long long int> interior;
};

} // namespace math

#endif // PYMESH_DIFFOPERATORS_H
