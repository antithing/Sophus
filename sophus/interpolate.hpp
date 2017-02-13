#ifndef SOPHUS_INTERPOLATE_HPP
#define SOPHUS_INTERPOLATE_HPP

#include <Eigen/Eigenvalues>

#include "interpolate_details.hpp"

namespace Sophus {

// This function interpolates between two Lie group elements ``foo_T_bar``
// and ``foo_T_baz`` with an interpolation factor of ``alpha`` in [0, 1].
//
// It returns a pose ``foo_T_quiz`` with ``quiz`` being a frame between ``bar``
// and ``baz``. If ``alpha=0`` it returns ``foo_T_bar``. If it is 1, it returns
// ``foo_T_bar``.
//
// (Since interpolation on Lie groups is inverse-invariant, we can equivalently
// think of the input arguments as being ``bar_T_foo``, ``baz_T_foo`` and the
// return value being ``quiz_T_foo``.)
//
// Precondition: ``alpha`` must be in [0, 1].
//
template <class G, class Scalar2>
std::enable_if_t<interp_details::Traits<G>::supported, G> interpolate(
    G const& foo_T_bar, G const& foo_T_baz, Scalar2 alpha = Scalar2(0.5f)) {
  using Scalar = typename G::Scalar;
  SOPHUS_ENSURE(alpha >= Scalar2(0) && alpha <= Scalar2(1),
                "alpha (%) must in [0, 1].");
  return foo_T_bar *
         G::exp(Scalar(alpha) * (foo_T_bar.inverse() * foo_T_baz).log());
}

}  // namespace Sophus

#endif  // SOPHUS_INTERPOLATE_HPP
