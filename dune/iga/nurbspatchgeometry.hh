// SPDX-FileCopyrightText: 2022 Alexander Müller mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "dune/iga/geometry/closestpointprojection.hh"
#include "dune/iga/geometry/geohelper.hh"
#include "dune/iga/nurbsalgorithms.hh"
#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/quadraturerules.hh>

namespace Dune::IGA {

  template <std::integral auto dim, std::integral auto dimworld, typename ScalarType = double>
  class NURBSPatchGeometry {
   public:
    static constexpr std::integral auto mydimension = dim;

    static constexpr std::integral auto coorddimension = dimworld;

    using ctype                     = ScalarType;
    using LocalCoordinate           = Dune::FieldVector<ScalarType, mydimension>;
    using GlobalCoordinate          = Dune::FieldVector<ScalarType, coorddimension>;
    using JacobianTransposed        = Dune::FieldMatrix<ScalarType, mydimension, coorddimension>;
    using JacobianInverseTransposed = Dune::FieldMatrix<ScalarType, coorddimension, mydimension>;
    using Jacobian                  = Dune::FieldMatrix<ctype, coorddimension, mydimension>;

    using Volume                    = ctype;

    using ControlPointType = typename NURBSPatchData<mydimension, dimworld, ScalarType>::ControlPointType;
    using ControlPointNetType = typename NURBSPatchData<mydimension, dimworld, ScalarType>::ControlPointNetType;
    using Nurbs = Dune::IGA::Nurbs<mydimension, ScalarType>;
    using NurbsLocalView = typename Nurbs::LocalView;

   private:
    /* Helper class to compute a matrix pseudo inverse */
    using MatrixHelper = typename MultiLinearGeometryTraits<ctype>::MatrixHelper;

   public:
    NURBSPatchGeometry() = default;

    struct LocalView {

      explicit LocalView(const NURBSPatchGeometry& patchGeometry):patchGeometry_{patchGeometry}{}
      [[nodiscard]] GlobalCoordinate center() const {
        LocalCoordinate localcenter(0.5);
        return global(localcenter);
      }

      void bind(const LocalCoordinate& u) {
        std::tie(nurbsLocalView_,localControlPointNet)= patchGeometry_->calculateNurbsAndControlPointNet(u);
      }

      /** \brief evaluates the geometric position
 *
 *  \param[in] u coordinates for each dimension in the [knotSpan.front(), knotSpan.back() ] domain
 */
      [[nodiscard]] FieldVector<ctype, dimworld> global(const LocalCoordinate& u) const {
        assert(subNetStart && "Bind the local view first!");
        return this->global(u, nurbsLocalView_,localControlPointNet);
      }

      [[nodiscard]] JacobianTransposed jacobianTransposed(const LocalCoordinate& u) const {
        return patchGeometry_->jacobianTransposedImpl(u, nurbsLocalView_,localControlPointNet);
      }

      [[nodiscard]] Jacobian jacobian(const LocalCoordinate& local) const {
        return transpose(jacobianTransposed(local));
      }

      [[nodiscard]] ctype integrationElement(const LocalCoordinate& local) const {
        auto j = jacobianTransposed(local);
        return MatrixHelper::template sqrtDetAAT<mydimension, coorddimension>(j);
      }

      [[nodiscard]] double volume(int scaleOrder = 1) const {
        const auto rule = QuadratureRules<ctype, mydimension>::rule(
            this->type(), scaleOrder * mydimension * (*std::ranges::max_element(patchGeometry_->patchData_.degree)));
        ctype vol = 0.0;
        for (auto& gp : rule)
          vol += integrationElement(gp.position()) * gp.weight();
        return vol;
      }

      [[nodiscard]] GeometryType type () const {
        return GeometryTypes::cube(mydimension);
      }


    private:
      std::optional<std::array<int, dim>> subNetStart;
      ControlPointNetType localControlPointNet;
      typename Nurbs::LocalView nurbsLocalView_;
      const NURBSPatchGeometry* patchGeometry_;
    };

    explicit NURBSPatchGeometry(const NURBSPatchData<mydimension, dimworld, ScalarType>& patchData)
        : patchData_(patchData), nurbs_{patchData} {}

    /** \brief Map the center of the element to the geometry */
    [[nodiscard]] GlobalCoordinate center() const {
      LocalCoordinate localcenter(0.5);
      return global(localcenter);
    }


    /** \brief Computes the volume of the element with an integration rule for order max(order)*elementdim */
    [[nodiscard]] double volume(int scaleOrder = 1) const {
      const auto rule = Dune::QuadratureRules<ctype, mydimension>::rule(
          this->type(), scaleOrder * mydimension * (*std::ranges::max_element(patchData_.degree)));
      ctype vol = 0.0;
      for (auto& gp : rule)
        vol += integrationElement(gp.position()) * gp.weight();
      return vol;
    }

    [[nodiscard]] bool affine() const { return false; }

    /** \brief Return the number of corners of the element */
    [[nodiscard]] int corners() const { return 1 << mydimension; }

    /** \brief Return world coordinates of the k-th corner of the element */
    [[nodiscard]] GlobalCoordinate corner(int k) const {
      LocalCoordinate localcorner;
      for (size_t i = 0; i < mydimension; i++) {
        localcorner[i] = (k & (1 << i)) ? 1 : 0;
      }
      return global(localcorner);
    }

    FieldVector<ctype, dimworld> operator()(const LocalCoordinate& local) const { return global(local); }


    /** \brief evaluates the geometric position
     *
     *  \param[in] u coordinates for each dimension in the [knotSpan.front(), knotSpan.back() ] domain
     */
    [[nodiscard]] FieldVector<ctype, dimworld> global(const LocalCoordinate& u) const {
      auto [nurbsLocalView,cpNet] = calculateNurbsAndControlPointNet(u);
      return this->global(u, nurbsLocalView,cpNet);
    }

    auto zeroFirstAndSecondDerivativeOfPosition(const LocalCoordinate& u) const {
      FieldVector<ctype, dimworld> pos;
      JacobianTransposed J;
      FieldMatrix<ctype, dim*(dim + 1) / 2, coorddimension> H;
      std::array<unsigned int, dim> subDirs;
      for (int subI = 0, i = 0; i < mydimension; ++i)
        subDirs[subI++] = i;

      const auto basisFunctionDerivatives = nurbs_.basisFunctionDerivatives(u, 2);
      auto cpCoordinateNet                = netOfSpan(u, patchData_.knotSpans, patchData_.degree,
                                                      extractControlCoordinates(patchData_.controlPoints));
      std::array<unsigned int, mydimension> ithVecZero{};
      pos = Dune::IGA::dot(basisFunctionDerivatives.get(ithVecZero), cpCoordinateNet);

      for (int dir = 0; dir < mydimension; ++dir) {
        std::array<unsigned int, mydimension> ithVec{};
        ithVec[subDirs[dir]] = 1;
        J[dir]               = dot(basisFunctionDerivatives.get(ithVec), cpCoordinateNet);
      }
      for (int dir = 0; dir < dim; ++dir) {
        std::array<unsigned int, mydimension> ithVec{};
        ithVec[subDirs[dir]] = 2;  // second derivative in dir direction
        H[dir]               = dot(basisFunctionDerivatives.get(ithVec), cpCoordinateNet);
      }
      if constexpr (dim > 1) {
        std::array<int, dim> mixeDerivs;
        std::ranges::fill_n(mixeDerivs.begin(), 2, 1);  // first mixed derivatives
        int mixedDireCounter = dim;
        do {
          H[mixedDireCounter++] = dot(basisFunctionDerivatives.get(mixeDerivs), cpCoordinateNet);
        } while (std::ranges::next_permutation(mixeDerivs, std::greater()).found);
      }
      return std::make_tuple(pos, J, H);
    }

    /** \brief Inverse of global this function gets a point defined in the world space and return
     * the closest point in local coordinates, i.e. in [0,1] domain for each grid dimension
     *
     *  \param global global coordinates for the point where the local coordinates are searched for
     */
    [[nodiscard]] LocalCoordinate local(const GlobalCoordinate& global) const {
      if constexpr (dim == 0)
        return {};
      else if constexpr (dim != dimworld) {
        auto [u, Ru, fu, gap] = Dune::IGA::closestPointProjectionByTrustRegion(*this, global);
        return u;
      } else {
        const ctype tolerance = ctype(16) * std::numeric_limits<ctype>::epsilon();
        LocalCoordinate x     = domainMidPoint();

        LocalCoordinate dx{};
        do {  // from multilinearGeometry
          const GlobalCoordinate dglobal = (*this).global(x) - global;
          MatrixHelper::template xTRightInvA<mydimension, coorddimension>(jacobianTransposed(x), dglobal, dx);
          const bool invertible
              = MatrixHelper::template xTRightInvA<mydimension, coorddimension>(jacobianTransposed(x), dglobal, dx);

          if (!invertible) return LocalCoordinate(std::numeric_limits<ctype>::max());
          x -= dx;
          // if local is outside the maximum knot vector span bound, thus we clamp it to it and hope for convergence
          for (int i = 0; i < mydimension; ++i) {
            if (Dune::FloatCmp::gt(x[i], patchData_.knotSpans[i].back())) x[i] = patchData_.knotSpans[i].back();
            if (Dune::FloatCmp::lt(x[i], patchData_.knotSpans[i].front())) x[i] = patchData_.knotSpans[i].front();
          }

        } while (dx.two_norm2() > tolerance);
        return x;
      }
    }

    /** \brief compute the Jacobian transposed matrix
     *
     *  \param local coordinates for each dimension
     */
    [[nodiscard]] JacobianTransposed jacobianTransposed(const LocalCoordinate& u) const {
      auto [nurbsLocalView,cpNet] = calculateNurbsAndControlPointNet(u);
      return jacobianTransposedImpl(u, nurbsLocalView,cpNet);
    }

    [[nodiscard]] Jacobian jacobian(const LocalCoordinate& local) const {
      return transpose(jacobianTransposed(local));
    }


    [[nodiscard]] ctype integrationElement(const LocalCoordinate& local) const {
      auto j = jacobianTransposed(local);
      return MatrixHelper::template sqrtDetAAT<mydimension, coorddimension>(j);
    }

    /** \brief Type of the element: a hypercube of the correct dimension */
    [[nodiscard]] GeometryType type() const { return GeometryTypes::cube(mydimension); }

    // The following are functions are not part of the Geometry Interface
    [[nodiscard]] std::array<Utilities::Domain<double>, mydimension> domain() const {
      std::array<Utilities::Domain<double>, mydimension> result{};
      for (int i = 0; i < mydimension; ++i)
        result[i] = {patchData_.knotSpans[i].front(), patchData_.knotSpans[i].back()};

      return result;
    }

    [[nodiscard]] std::array<int, mydimension> degree() const { return patchData_.degree; }

    [[nodiscard]] FieldVector<ctype, mydimension> domainMidPoint() const {
      auto dom = domain();
      Dune::FieldVector<ctype, mydimension> result{};
      for (int i = 0; i < mydimension; ++i)
        result[i] = dom[i].center();

      return result;
    }

  private:
    auto calculateNurbsAndControlPointNet(const LocalCoordinate& u) const {
      auto subNetStart = findSpanCorrected(patchData_.degree, u, patchData_.knotSpans);

      auto cpCoordinateNet = netOfSpan(subNetStart,patchData_.degree,
                                       extractControlCoordinates(patchData_.controlPoints));
      auto nurbsLocalView = nurbs_.localView();
      nurbsLocalView.bind(subNetStart);
      return std::make_pair(nurbsLocalView,cpCoordinateNet);
    }
    [[nodiscard]] static FieldVector<ctype, dimworld> global(const LocalCoordinate& u,const NurbsLocalView& nurbsLocalView,
      const MultiDimensionNet<dim, ControlPoint<FieldVector<ctype,dimworld>>>& localControlPointNet)  {
      auto basis           = nurbsLocalView.basisFunctions(u);

      return IGA::dot(basis, localControlPointNet);
    }

    [[nodiscard]] static JacobianTransposed jacobianTransposedImpl(const LocalCoordinate& u,const NurbsLocalView& nurbsLocalView,
      const MultiDimensionNet<dim, ControlPoint<FieldVector<ctype,dimworld>>>& localControlPointNet)  {

      JacobianTransposed result;

      const auto basisFunctionDerivatives = nurbsLocalView.basisFunctionDerivatives(u, 1);

      for (int dir = 0; dir < mydimension; ++dir) {
        std::array<unsigned int, mydimension> ithVec{};
        ithVec[dir] = 1;
        result[dir] = dot(basisFunctionDerivatives.get(ithVec), localControlPointNet);
      }
      return result;
    }

   public:
    NURBSPatchData<mydimension, dimworld, ScalarType> patchData_;

   private:
    Dune::IGA::Nurbs<mydimension, ScalarType> nurbs_;
  };

}  // namespace Dune::IGA
