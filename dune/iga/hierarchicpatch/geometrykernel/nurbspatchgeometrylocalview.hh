// SPDX-FileCopyrightText: 2022 Alexander Müller mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

// #include <dune/iga/hierarchicpatch/geometrykernel/geohelper.hh>
#include <dune/common/diagonalmatrix.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/iga/hierarchicpatch/enums.hh>
#include <dune/iga/hierarchicpatch/geometrykernel/higherorderalgorithms.hh>
#include <dune/iga/hierarchicpatch/splines/nurbsalgorithms.hh>

namespace Dune::IGANEW::GeometryKernel {

  namespace Impl {
    struct IndexPair {
      int row;
      int col;
    };

    template <int dim>
    auto& voigtIndices() {
      if constexpr (dim == 1) {
        static std::array<IndexPair, dim*(dim + 1) / 2> voigt1D = {{0, 0}};
        return voigt1D;
      } else if constexpr (dim == 2) {
        static std::array<IndexPair, dim*(dim + 1) / 2> voigt2D = {{{0, 0}, {1, 1}, {0, 1}}};
        return voigt2D;
      } else if constexpr (dim == 3) {
        static std::array<IndexPair, dim*(dim + 1) / 2> voigt3D = {{{0, 0}, {1, 1}, {2, 2}, {1, 2}, {0, 2}, {0, 1}}};
        return voigt3D;
      }
    }
  }  // namespace Impl

  template <int codim, typename PatchGeometry, Trimming trim_>
  struct PatchGeometryLocalView {
    using ctype                                         = typename PatchGeometry::ctype;
    static constexpr int gridDimension                  = PatchGeometry::mydimension;
    static constexpr int mydimension                    = gridDimension - codim;
    static constexpr int numberOfSecondDerivatives      = mydimension * (mydimension + 1) / 2;
    static constexpr int patchNumberOfSecondDerivatives = gridDimension * (gridDimension + 1) / 2;
    static constexpr Trimming trim                      = trim_;

    static constexpr std::integral auto worlddimension = PatchGeometry::worlddimension;

    using LocalCoordinate         = FieldVector<ctype, mydimension>;
    using GlobalCoordinate        = typename PatchGeometry::GlobalCoordinate;
    using JacobianTransposed      = FieldMatrix<ctype, mydimension, worlddimension>;
    using PatchJacobianTransposed = typename PatchGeometry::JacobianTransposed;
    using PatchHessian            = typename PatchGeometry::Hessian;
    // TODO trim ParameterSpaceGeometry
    using ParameterSpaceGeometry = typename PatchGeometry::template ParameterSpaceGeometry<codim>;

    //! if we have codim==0, then the Jacobian in the parameter space of the grid entity itself is a DiagonalMatrix, and
    // Coordinates in a single knot span differ from coordinates on the B-spline patch
    // by an affine transformation.  This transformation is stored in the diagonal entries.
    // If trimming is disabled the Jacobian in the parameter space of subentities (edges or surfaces) is a
    // DiagonalMatrixBlock but is treated as a FieldMatrix, since is also just cubes and an axis-aligned geometry, if
    // trimming is enabled the Jacobian of the parameterspace for subentities is a potentially fully populated
    // FieldMatrix, Think about an arbitrary curve in the axis-aligned cube parameter grid
    //                      ------->
    // `____________________C_______D_______`
    // ::````````````````::``````````'|````::
    // ::                ::          `|    ::
    // ::              ^F::          `|    ::
    // ::              | ::          `|    ::
    // ::              | ::          `|    ::
    // ::              |E::          `|    ::
    // ::                ::          .|  x ::
    // ::                ::          _:    ::
    // ::................::..........|:....::
    // ::''''''''''''''''::''''''''''/'''''::
    // ::                ::         _:B    ::
    // ::                ::        .|      ::
    // ::                ::      `::A      ::
    // ::           '::::||::::::'         ::
    // ::         ::`    ::          x     ::
    // ::        |'      ::                ::
    // ::        |   x   ::                ::
    // ::````````|```````::````````````````::
    // `____________________________________`
    // at point A the tangent of the curve (the JacobianInParameterSpace) points into the oblique direction (to B)
    // for the axis-aligned geometry these tangents are from E to F or from C to D, depending on the orientation of the
    // subentity, which are aligned with the axis
    using JacobianTransposedInParameterSpace = typename ParameterSpaceGeometry::JacobianTransposed;
    using GlobalInParameterSpace             = typename ParameterSpaceGeometry::GlobalCoordinate;

    using Hessian                   = FieldMatrix<ctype, numberOfSecondDerivatives, worlddimension>;
    using Jacobian                  = FieldMatrix<ctype, worlddimension, mydimension>;
    using JacobianInverseTransposed = FieldMatrix<ctype, worlddimension, mydimension>;
    using JacobianInverse           = FieldMatrix<ctype, mydimension, worlddimension>;
    using MatrixHelper              = typename PatchGeometry::MatrixHelper;
    using Volume                    = ctype;

    using Nurbs          = Splines::Nurbs<gridDimension, ctype>;
    using NurbsLocalView = typename Nurbs::LocalView;

    using ControlPointCoordinateNetType = typename PatchGeometry::ControlPointCoordinateNetType;

    PatchGeometryLocalView() = default;
    explicit PatchGeometryLocalView(const PatchGeometry& patchGeometry) : patchGeometry_{&patchGeometry} {}

    [[nodiscard]] GlobalCoordinate center() const { return global(LocalCoordinate(0.5)); }

    void bind(const ParameterSpaceGeometry& lGeo) {
      parameterSpaceGeometry = std::make_optional<ParameterSpaceGeometry>(lGeo);
      std::tie(nurbsLocalView_, localControlPointNet, spanIndices_)
          = patchGeometry_->calculateNurbsAndControlPointNet(lGeo.center());
    }

    /**
     * \brief evaluates the geometric position
     * \param[in] local coordinate
     * \return position in world space
     */
    [[nodiscard]] GlobalCoordinate global(const LocalCoordinate& u) const {
      checkState();
      return GeometryKernel::position(globalInParameterSpace(u), nurbsLocalView_, localControlPointNet);
    }

    [[nodiscard]] JacobianTransposed jacobianTransposed(const LocalCoordinate& u) const {
      checkState();
      const PatchJacobianTransposed dfdg
          = GeometryKernel::jacobianTransposed(globalInParameterSpace(u), nurbsLocalView_, localControlPointNet);
      const JacobianTransposedInParameterSpace dgdt = jacobianTransposedInParameterSpace(u);
      return dgdt * dfdg;
    }

    [[nodiscard]] Jacobian jacobian(const LocalCoordinate& local) const { return transpose(jacobianTransposed(local)); }

    [[nodiscard]] ctype integrationElement(const LocalCoordinate& local) const {
      auto jT = jacobianTransposed(local);
      return MatrixHelper::template sqrtDetAAT<mydimension, worlddimension>(jT);
    }

    [[nodiscard]] double volume(int scaleOrder = 1) const {
      if constexpr (trim == Trimming::Enabled) {
        // if constexpr (mydimension == 2)
        // if (subgrid_) {
        //   Dune::QuadratureRule<double, mydimension> rule;
        //   fillQuadratureRuleImpl(rule, *subgrid_.get(), (*std::ranges::max_element(patchData().degree)));
        //   Volume vol = 0.0;
        //   for (auto& gp : rule)
        //     vol += integrationElement(gp.position()) * gp.weight();
        //   return vol;
        // }
        // TODO here the integration of trimmed quantities has to happen and also the new edge geometries
      }

      const auto& rule = QuadratureRules<ctype, mydimension>::rule(
          GeometryTypes::cube(mydimension), (*std::ranges::max_element(patchGeometry_->patchData_.degree)));
      Volume vol = 0.0;
      for (auto& gp : rule)
        vol += integrationElement(gp.position()) * gp.weight();
      return vol;
    }

    /** \brief Type of the type of the parameter space element */
    [[nodiscard]] GeometryType type() const { return GeometryTypes::cube(mydimension); }

    /** \brief Return the number of corners of the element */
    [[nodiscard]] int corners() const { return 1 << mydimension; }

    /** \brief Return world coordinates of the k-th corner of the element */
    [[nodiscard]] GlobalCoordinate corner(int k) const {
      LocalCoordinate localcorner;
      for (size_t i = 0; i < mydimension; i++)
        localcorner[i] = (k & (1 << i)) ? 1 : 0;

      return global(localcorner);
    }

    /** \brief Inverse of global this function gets a point defined in the world space and return
     * the closest point in local coordinates, i.e. in [0,1] domain for each grid dimension
     *
     *  \param global global coordinates for the point where the local coordinates are searched for
     */
    [[nodiscard]] LocalCoordinate local(const GlobalCoordinate& global) const {
      return GeometryKernel::computeParameterSpaceCoordinate(*this, global);
    }

    [[nodiscard]] JacobianInverseTransposed jacobianInverseTransposed(const LocalCoordinate& local) const {
      JacobianTransposed Jt = jacobianTransposed(local);
      JacobianInverseTransposed jacobianInverseTransposed;
      MatrixHelper::template rightInvA<mydimension, worlddimension>(Jt, jacobianInverseTransposed);
      return jacobianInverseTransposed;
    }

    std::tuple<GlobalCoordinate, JacobianTransposed, Hessian> zeroFirstAndSecondDerivativeOfPosition(
        const LocalCoordinate& u) const {
      checkState();

      return std::make_tuple(global(u), jacobianTransposed(u), hessian(u));
    }

    /**
     * @brief Computes the second derivatives of the local view.
     *
     * @param[in] u The parameter value of the local view on the geometry
     * @return A tuple of the global position, the Jacobian matrix and Hessian matrix
     * @details
     * **Example curve on surface**
     *
     * The local view contains a `parameterSpaceGeometry`, which describes the  composition function
     *
     * \f[g: \begin{cases}\mathbb{R} \rightarrow \mathbb{R}^2 \\ t \mapsto g(t) \end{cases},\f]
     *
     * which describes the curve's local parametrization in terms of the patch geometries parametrization
     * The patch geometry is described by
     *
     * \f[f: \begin{cases}\mathbb{R}^2 \rightarrow \mathbb{R}^3 \\ (u,v) \mapsto * f(u,v) \end{cases},\f]
     *
     * which describes a surface in 3D.
     * Then the curve's realization on the surface in the 3D space can be written by
     *
     * \f[h: \begin{cases}\mathbb{R} \rightarrow \mathbb{R}^3 \\ t \mapsto f(g_1(t),g_2(t)) \end{cases},\f]
     *
     * This function calculates the Jacobian and Hessian matrices for the composition function \f$h(g_1(t), g_2(t))\f$
     * where \f$g(t)\f$ is the original curve and \f$h(u, v)\f$  is the patch geometry function. The computation
     * involves chain and product rules between \f$g(t)\f$ and \f$h(u, v)\f$.
     *     *
     * The Jacobian matrix J is given by
     * \f[
     * J =
     * \frac{\partial f }{\partial g_1}  \frac{\partial g_1}{\partial t} + \frac{\partial f}{\partial g_2}
     * \frac{\partial g_2}{\partial t} \f] and has dimensions \f$3\times 1\f$ or in general
     * \f$\verb+worlddimension+\times \verb+mydimension+\f$.
     *
     * The Hessian matrix H is given by:
     * \f[
     * H =
     * \frac{\partial^2 f}{\partial g_1\partial g_1}  \left(\frac{\partial g_1}{\partial t}\right)^2 +\frac{\partial^2
     * f}{\partial g_2\partial g_2}  \left(\frac{\partial g_2}{\partial t}\right)^2 + 2\frac{\partial^2 f}{\partial
     * g_1\partial g_2}  \frac{\partial g_1}{\partial t} \frac{\partial g_2}{\partial t} + \frac{\partial f}{\partial
     * g_1}\frac{\partial^2 g_1}{\partial t^2}+ \frac{\partial f}{\partial g_2}\frac{\partial^2 g_2}{\partial t^2} \f]
     *
     * and has dimensions \f$1\times 3\f$ or in general \f$\verb|mydimension*(mydimension+1)/2| \times
     * \verb+worlddimension+\f$.
     */
    Hessian hessian(const LocalCoordinate& u) const {
      checkState();
      const auto ouInPatch = globalInParameterSpace(u);

      const auto [p, dfdg, dfdgdg] = patchGeometry_->zeroFirstAndSecondDerivativeOfPositionImpl(
          ouInPatch, nurbsLocalView_, localControlPointNet);

      const JacobianTransposedInParameterSpace dgdt = jacobianTransposedInParameterSpace(u);
      static_assert(JacobianTransposedInParameterSpace::rows == mydimension);
      static_assert(JacobianTransposedInParameterSpace::cols == gridDimension);
      // dgdt (1x2, curve on surface), (2x2 surface in surface), (2x3 surface in 3D), (3x3 volume in 3D)
      FieldMatrix<ctype, patchNumberOfSecondDerivatives, numberOfSecondDerivatives> dgdtSquared(0);

      for (int patchIndex = 0; auto [patchRow, patchCol] : Impl::voigtIndices<gridDimension>()) {
        for (int myIndex = 0; auto [myRow, myCol] : Impl::voigtIndices<mydimension>()) {
          if constexpr (not std::is_same_v<JacobianTransposedInParameterSpace, DiagonalMatrix<ctype, gridDimension>>)
            dgdtSquared[patchIndex++][myIndex++]
                = (dgdt[myRow][patchRow] * dgdt[myCol][patchCol] + dgdt[myCol][patchRow] * dgdt[myRow][patchCol])
                  * ((patchRow == patchCol) ? 0.5 : 1.0);
          else {
            const auto dgdtmyRowpatchRow = (myRow == patchRow) ? dgdt[myRow][patchRow] : 0;
            const auto dgdtmyColpatchCol = (myCol == patchCol) ? dgdt[myCol][patchCol] : 0;
            const auto dgdtmyColpatchRow = (myRow == patchCol) ? dgdt[myCol][patchRow] : 0;
            const auto dgdtmyRowpatchCol = (myCol == patchRow) ? dgdt[myRow][patchCol] : 0;
            dgdtSquared[patchIndex++][myIndex++]
                = (dgdtmyRowpatchRow * dgdtmyColpatchCol + dgdtmyColpatchRow * dgdtmyRowpatchCol)
                  * ((patchRow == patchCol) ? 0.5 : 1.0);
          }
        }
      }

      Hessian h;
      // h = mydimension * (mydimension + 1) / 2 X worlddimension = 1x2 (curve in 2D), 1x3 (curve in 3D), (3x2 surface
      // in 2D), (3x3 surface in 3D), (6x3 volume in 3D)
      //
      // dfdgdg = gridDimension*(gridDimension + 1) / 2 X worlddimension
      // = 3x2 (curve on surface flat),3x3 (curve on surface in 3D), 6x3 (curve in volume), 6x3 (surface in volume),(3x3
      // surface in surface in 3D)
      //
      // dgdtSquared = gridDimension * (gridDimension + 1) / 2 X mydimension * (mydimension +
      // 1) / 2 = 3x1 (curve on surface), 6x1 (curve in volume), 6x3 (surface in volume), 3x3 ( surface in surface), 6x6
      // ( volume in volume)
      h = transposedView(dgdtSquared) * dfdgdg;

      /* if trimming is enabled the parameter space geometry is potentially non-linear,
       * the resutling Hessian has another contribution due to chain-rule, namely the second derivative of g */
      if constexpr (trim == Trimming::Enabled) {
        const auto dgdtdt = parameterSpaceGeometry.value().hessian(ouInPatch);

        assert(trim == Trimming::Disabled
               && "This can not be checked yet. Check if this works with trimming and then remove assert");
        h += transpose(transposedView(dfdg) * dgdtdt);
      }
      return h;
    }

    [[nodiscard]] LocalCoordinate domainMidPoint() const {
      auto dom = domain();
      LocalCoordinate result{};
      for (int i = 0; i < mydimension; ++i)
        result[i] = dom[i].center();

      return result;
    }

    [[nodiscard]] std::array<Utilities::Domain<double>, mydimension> domain() const { return {}; }

    [[nodiscard]] bool affine() const { return false; }

   private:
    void checkState() const { assert(parameterSpaceGeometry && "Bind the local view first!"); }
    GlobalInParameterSpace globalInParameterSpace(const LocalCoordinate& local) const {
      checkState();
      return parameterSpaceGeometry.value().global(local);
    }

    [[nodiscard]] JacobianTransposedInParameterSpace jacobianTransposedInParameterSpace(
        const LocalCoordinate& u) const {
      checkState();
      return parameterSpaceGeometry.value().jacobianTransposed(u);
    }

    ControlPointCoordinateNetType localControlPointNet;
    NurbsLocalView nurbsLocalView_;
    std::array<int, gridDimension> spanIndices_;
    const PatchGeometry* patchGeometry_;
    std::optional<ParameterSpaceGeometry> parameterSpaceGeometry;
  };

}  // namespace Dune::IGANEW::GeometryKernel
