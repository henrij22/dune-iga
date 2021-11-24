//
// Created by lex on 21.10.21.
//

#pragma once

#include <concepts>
#include <numbers>
#include <ranges>

#include <dune/iga/bsplinealgorithms.hh>
#include <dune/iga/multidimensionNet.hh>
#include <dune/iga/nurbspatchdata.hh>
namespace Dune::IGA {

  template <Vector VectorType>
  requires(VectorType::dimension == 3) VectorType cross(const VectorType& a, const VectorType& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
  }

  template <std::integral auto dim>
  auto ordersFromDegrees(const std::array<int, dim>& degree) {
    std::array<int, dim> order;
    std::ranges::transform(degree, order.begin(), [](const auto& p) { return p + 1; });
    return order;
  }

  template <std::floating_point ScalarType, std::integral auto dim, typename NetValueType>
  auto netOfSpan(const std::array<ScalarType, dim> u, const std::array<std::vector<ScalarType>, dim>& knots,
                 const std::array<int, dim>& degree, const MultiDimensionNet<dim, NetValueType>& net) {
    std::array<int, dim> order = ordersFromDegrees(degree);
    auto subNetStart           = findSpan(degree, u, knots);
    for (std::size_t i = 0; i < dim; ++i)
      subNetStart[i] -= degree[i];

    return net.subNet(subNetStart, order);
  }

  template <std::integral auto dim>
  auto createPartialSubDerivativPermutations(const FieldVector<int, dim>& v) {
    std::vector<FieldVector<int, dim>> perm;
    for (int i = 1; i < std::pow(2, v.size()); ++i) {
      FieldVector<int, dim> x{};
      for (int j = 0; j < v.size(); ++j)
        if ((i & (1 << j)) != 0) {
          x[j] = v[j];
          if (v[j] == 0) goto outer;
        }
      perm.push_back(x);
    outer:;
    }
    return perm;
  }

  template <Vector VectorType>
  int binom(const VectorType& n, const VectorType& k) {
    return std::inner_product(n.begin(), n.end(), k.begin(), 1, std::multiplies{},
                              [](auto& ni, auto& ki) { return Dune::binomial(ni, ki); });
  }

  template <std::integral auto dim, ControlPointConcept ValueType>
  auto extractWeights(const MultiDimensionNet<dim, ValueType>& cpsandWeight) {
    auto viewOverWeights = std::ranges::transform_view(cpsandWeight.directGetAll(), [](auto& cp) { return cp.w; });
    return MultiDimensionNet<dim, typename ValueType::VectorType::value_type>(cpsandWeight.size(), viewOverWeights);
  }

  template <std::integral auto dim, ControlPointConcept ValueType>
  auto extractControlpoints(const MultiDimensionNet<dim, ValueType>& cpsandWeight) {
    auto viewOverCps = std::ranges::transform_view(cpsandWeight.directGetAll(), [](auto& cp) { return cp.p; });
    return MultiDimensionNet<dim, typename ValueType::VectorType>(cpsandWeight.size(), viewOverCps);
  }

  template <std::floating_point ScalarType, std::size_t dim>
  class Nurbs {
  public:
    template <std::integral auto dimworld>
    Nurbs(const Dune::IGA::NURBSPatchData<dim, dimworld>& data)
        : knots_{data.knotSpans}, degree_{data.order}, weights_{extractWeights(data.controlPoints)} {}

    Nurbs(const std::array<std::vector<ScalarType>, dim>& knots, const std::array<int, dim>& degree, const std::vector<ScalarType>& weights)
        : knots_{knots}, degree_{degree}, weights_{weights} {}

    auto operator()(const std::array<ScalarType, dim>& u) { return basisFunctions(u, knots_, degree_, weights_).directGetAll(); }

    auto basisFunctionNet(const std::array<ScalarType, dim>& u) const { return basisFunctions(u, knots_, degree_, weights_); }

    static auto basisFunctions(const std::array<ScalarType, dim> u, const std::array<std::vector<ScalarType>, dim>& knots,
                               const std::array<int, dim>& degree, const MultiDimensionNet<dim, ScalarType>& weights) {
      const std::array<int, dim> order = ordersFromDegrees(degree);
      std::array<std::vector<ScalarType>, dim> bSplines;

      for (std::size_t i = 0; i < dim; ++i)
        bSplines[i] = Bspline<ScalarType>::basisFunctions(u[i], knots[i], degree[i]);

      auto Nnet = MultiDimensionNet<dim, ScalarType>(bSplines);

      const auto subNetWeights = netOfSpan<ScalarType, dim>(u, knots, degree, weights);

      const ScalarType invSumWeight = dot(Nnet, subNetWeights);
      Nnet *= subNetWeights;
      Nnet /= invSumWeight;
      return Nnet;
    }

    // \*brief This function return the basis function and the corresponding derivatives
    // it generalizes the formula (4.20) of Piegl and Tiller 97 for a basis of dimension dim and up to the derivative order derivativeOrder
    static auto basisFunctionDerivatives(const std::array<ScalarType, dim> u, const std::array<std::vector<ScalarType>, dim>& knots,
                                         const std::array<int, dim>& degree, const MultiDimensionNet<dim, double>& weights,
                                         const int derivativeOrder, const bool triangleDerivatives = false) {
      std::array<DynamicMatrix<ScalarType>, dim> bSplineDerivatives;
      for (int i = 0; i < dim; ++i)
        bSplineDerivatives[i] = Bspline<ScalarType>::basisFunctionDerivatives(u[i], knots[i], degree[i], derivativeOrder);

      std::array<std::vector<ScalarType>, dim> dimArrayOfVectors;
      FieldVector<int, dim> dimSize(derivativeOrder + 1);
      MultiDimensionNet<dim, MultiDimensionNet<dim, ScalarType>> netsOfDerivativeNets(dimSize);

      for (int j = 0; auto& derivNet : netsOfDerivativeNets.directGetAll()) {
        auto multiIndex = netsOfDerivativeNets.directToMultiIndex(j++);
        for (int i = 0; i < multiIndex.size(); ++i)
          dimArrayOfVectors[i] = bSplineDerivatives[i][multiIndex[i]].container();
        derivNet = MultiDimensionNet<dim, ScalarType>(dimArrayOfVectors);
      }

      MultiDimensionNet<dim, MultiDimensionNet<dim, ScalarType>> R = netsOfDerivativeNets;
      MultiDimensionNet<dim, ScalarType> netsOfWeightfunctions(dimSize);
      const auto subNetWeights = netOfSpan<ScalarType, dim>(u, knots, degree, weights);

      for (int j = 0; j < R.directSize(); ++j) {
        R.directGet(j) *= subNetWeights;
        netsOfWeightfunctions.directGet(j) = dot(netsOfDerivativeNets.directGet(j), subNetWeights);
      }

      for (int j = 0; j < R.directSize(); ++j) {
        const auto derivOrders = R.template directToMultiIndex<FieldVector<int, dim>>(j);
        if (triangleDerivatives)
          if (std::accumulate(derivOrders.begin(), derivOrders.end(), derivativeOrder)) continue;
        const auto perms = createPartialSubDerivativPermutations(derivOrders);

        for (const auto& perm : perms) {
          MultiDimensionNet<dim, int> kNet(perm + FieldVector<int, dim>(1));
          auto startMultiIndex = perm;
          std::ranges::transform(startMultiIndex, startMultiIndex.begin(), [](auto& v) { return (v != 0); });

          for (int kk = kNet.index(startMultiIndex); kk < kNet.directSize(); ++kk) {
            const auto k = kNet.template directToMultiIndex<FieldVector<int, dim>>(kk);
            R.directGet(j) -= binom(perm, k) * (R.get(derivOrders - k) * netsOfWeightfunctions.get(k));  // generalized Piegl Tiller (4.20)
          }
        }
        R.directGet(j) /= netsOfWeightfunctions.directGet(0);
      }
      return R;
    }

    auto basisFunctionDerivatives(const std::array<ScalarType, dim>& u, const int derivativeOrder) const {
      return basisFunctionDerivatives(u, knots_, degree_, weights_, derivativeOrder);
    }

  private:
    std::array<std::vector<ScalarType>, dim> knots_;
    std::array<int, dim> degree_;
    MultiDimensionNet<dim, ScalarType> weights_;
  };

  template <std::integral auto dim, std::integral auto dimworld, NurbsGridLinearAlgebra NurbsGridLinearAlgebraTraitsImpl>
  auto knotRefinement(const NURBSPatchData<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl>& oldData, const std::vector<double>& newKnots,
                      const int refinementDirection) {
    using NurbsPatchData = NURBSPatchData<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl>;
    using ControlPoint   = typename NURBSPatchData<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl>::ControlPointType;
    using namespace std::ranges;
    std::array<std::vector<double>, dim> newKnotsArray;
    std::array<int, dim - 1> otherDirections;
    for (int counter = 0, i = 0; i < dim; ++i) {
      if (i == refinementDirection) continue;
      newKnotsArray[i]           = oldData.knotSpans[i];
      otherDirections[counter++] = i;
    }

    typename NurbsPatchData::ControlPointNetType newCPv(oldData.controlPoints.size());

    auto oldCPv     = oldData.controlPoints;
    auto oldKnotVec = oldData.knotSpans[refinementDirection];

    auto& newKnotVec = newKnotsArray[refinementDirection];

    newKnotVec.reserve(oldKnotVec.size() + newKnots.size());
    merge(oldKnotVec, newKnots, std::back_inserter(newKnotVec));

    const auto newKSize                   = newKnots.size();
    std::array<int, dim> numberOfCPperDir = oldCPv.size();
    numberOfCPperDir[refinementDirection] += newKSize;

    newCPv.resize(numberOfCPperDir);

    auto scaleCPWithW = [](const auto& cp) -> ControlPoint { return {.p = cp.w * cp.p, .w = cp.w}; };

    const int p = oldData.order[refinementDirection];
    const int a = findSpan(p, newKnots.front(), oldKnotVec);
    const int b = findSpan(p, newKnots.back(), oldKnotVec);

    auto newSurfI = newCPv.hyperSurfBegin(otherDirections);
    auto oldSurfI = oldCPv.hyperSurfBegin(otherDirections);

    for (auto hSOld = oldSurfI, hSNew = newSurfI; hSOld != oldCPv.hyperSurfBegin(otherDirections) + (a - p + 1); ++hSOld, ++hSNew)
      std::ranges::transform(*hSOld, (*hSNew).begin(), scaleCPWithW);
    for (auto hSOld = oldSurfI + b, hSNew = newSurfI + b; hSOld != oldCPv.hyperSurfEnd(otherDirections); ++hSOld, ++hSNew)
      std::ranges::transform(*hSOld, (*hSNew).begin(), scaleCPWithW);

    newSurfI = newCPv.hyperSurfEnd(otherDirections) - 1;
    oldSurfI = oldCPv.hyperSurfEnd(otherDirections) - 1;

    auto currentNewKnot = newKnotVec.end();
    auto currentOldKnot = oldKnotVec.end();

    for (auto const& currentAdditionalKnot : reverse_view(newKnots)) {
      while (currentAdditionalKnot <= *(currentOldKnot - 1) && std::distance(oldKnotVec.begin(), currentOldKnot) > a + 1) {
        std::ranges::transform(*oldSurfI, (*newSurfI).begin(), scaleCPWithW);
        --currentNewKnot, --currentOldKnot;
        --newSurfI, --oldSurfI;
      }

      std::ranges::copy(*(newSurfI + 1), (*newSurfI).begin());
      ++newSurfI;
      currentOldKnot -= p;
      for ([[maybe_unused]] int _ : iota_view{0, p}) {
        ++newSurfI;

        auto alpha = *currentNewKnot - currentAdditionalKnot;
        using std::abs;
        if (abs(alpha) < 1e-7)
          std::ranges::copy(*newSurfI, (*(newSurfI - 1)).begin());
        else {
          alpha /= (*currentNewKnot - *(currentOldKnot));
          std::ranges::transform(*newSurfI, *(newSurfI - 1), (*(newSurfI - 1)).begin(),
                                 [&alpha](auto& cp, auto& cpL) { return alpha * cpL + cp * (1.0 - alpha); });
        }
        ++currentOldKnot, ++currentNewKnot;
      }
      currentNewKnot -= p + 1;
      newSurfI -= p + 2;
    }

    std::ranges::transform(newCPv.directGetAll(), newCPv.directGetAll().begin(),
                           [](auto& cp) -> ControlPoint { return {.p = cp.p / cp.w, .w = cp.w}; });

    return NURBSPatchData<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl>(newKnotsArray, newCPv, oldData.order);
  }

  namespace Impl {
    template <Vector VectorType>
    auto projectPointOntoLine(const VectorType basepoint, const VectorType revolutionaxis, const VectorType point) {
      const VectorType e1 = basepoint + revolutionaxis - basepoint;
      const VectorType e2 = point - basepoint;
      using Dune::dot;
      const typename VectorType::value_type angle = dot(e1, e2);
      const typename VectorType::value_type len2  = e1.two_norm();

      VectorType p = basepoint + angle * e1 / len2;
      return p;
    }

    template <Vector VectorType>
    auto Intersect3DLines(const VectorType basepoint1, const VectorType direction1, const VectorType basepoint2,
                          const VectorType direction2) {
      using ScalarType                     = typename VectorType::value_type;
      const ScalarType tol                 = 1e-8;
      const VectorType basePointDifference = basepoint1 - basepoint2;
      using Dune::dot;
      const ScalarType dir1squaredLength = dot(direction1, direction1);
      const ScalarType angle             = dot(direction1, direction2);
      const ScalarType dir2squaredLength = dot(direction2, direction2);
      const ScalarType dir1BasePdiff     = dot(direction1, basePointDifference);
      const ScalarType dir2BasePdiff     = dot(direction2, basePointDifference);
      const ScalarType D                 = dir1squaredLength * dir2squaredLength - angle * angle;

      if (std::abs(D) < tol) throw std::logic_error("The two lines are almost parallel.");

      const ScalarType parameter1 = (angle * dir2BasePdiff - dir2squaredLength * dir1BasePdiff) / D;
      const ScalarType parameter2 = (dir1squaredLength * dir2BasePdiff - angle * dir1BasePdiff) / D;

      VectorType c1 = basepoint1 + parameter1 * direction1;
      assert(dot(c1 - (basepoint2 + parameter2 * direction2), c1 - (basepoint2 + parameter2 * direction2)) < tol
             && "Both calculated points do not coincide");
      return c1;
    }
  }  // namespace Impl

  // Piegl and Tiller Algo A7.1
  template <NurbsGridLinearAlgebra NurbsGridLinearAlgebraTraitsImpl = Dune::IGA::LinearAlgebraTraits<double, 2UL, 3UL>>
  auto makeCircularArc(const typename NurbsGridLinearAlgebraTraitsImpl::GlobalCoordinateType::value_type radius     = 1.0,
                       const typename NurbsGridLinearAlgebraTraitsImpl::GlobalCoordinateType::value_type startAngle = 0.0,
                       typename NurbsGridLinearAlgebraTraitsImpl::GlobalCoordinateType::value_type endAngle         = 360.0,
                       const typename NurbsGridLinearAlgebraTraitsImpl::GlobalCoordinateType origin                 = {0, 0, 0},
                       const typename NurbsGridLinearAlgebraTraitsImpl::GlobalCoordinateType X                      = {1, 0, 0},
                       const typename NurbsGridLinearAlgebraTraitsImpl::GlobalCoordinateType Y                      = {0, 1, 0}) {
    using ScalarType           = typename NurbsGridLinearAlgebraTraitsImpl::GlobalCoordinateType::value_type;
    using GlobalCoordinateType = typename NURBSPatchData<1UL, 3UL, NurbsGridLinearAlgebraTraitsImpl>::GlobalCoordinateType;
    const auto pi              = std::numbers::pi_v<typename NurbsGridLinearAlgebraTraitsImpl::GlobalCoordinateType::value_type>;

    if (endAngle < startAngle) endAngle += 360.0;
    const ScalarType theta = endAngle - startAngle;
    const int narcs        = std::ceil(theta / 90);

    typename NURBSPatchData<1UL, 3UL, NurbsGridLinearAlgebraTraitsImpl>::ControlPointNetType circleCPs(2 * narcs + 1);
    const ScalarType dtheta  = theta / narcs * pi / 180;
    const int n              = 2 * narcs;
    const ScalarType w1      = cos(dtheta / 2.0);
    GlobalCoordinateType PO  = origin + radius * cos(startAngle) * X + radius * sin(startAngle) * Y;
    GlobalCoordinateType TO  = -sin(startAngle) * X + cos(startAngle) * Y;
    circleCPs.directGet(0).p = PO;
    ScalarType angle         = startAngle;
    for (int index = 0, i = 0; i < narcs; ++i) {
      angle += dtheta;
      const GlobalCoordinateType P2    = origin + radius * cos(angle) * X + radius * sin(angle) * Y;
      circleCPs.directGet(index + 2).p = P2;
      const GlobalCoordinateType T2    = -sin(angle) * X + cos(angle) * Y;
      const GlobalCoordinateType P1    = Impl::Intersect3DLines(PO, TO, P2, T2);
      circleCPs.directGet(index + 1)   = {.p = P1, .w = w1};
      index += 2;
      if (i < narcs - 1) {
        PO = P2;
        TO = T2;
      }
    }

    auto knotVec = std::array<std::vector<double>, 1>{};
    auto& U      = knotVec[0];
    U.resize(2 * (narcs + 2));
    if (narcs == 1) {
    } else if (narcs == 2) {
      U[3] = U[4] = 1.0 / 2.0;
    } else if (narcs == 3) {
      U[3] = U[4] = 1.0 / 3.0;
      U[5] = U[6] = 2.0 / 3.0;
    } else {
      U[3] = U[4] = 1.0 / 4.0;
      U[5] = U[6] = 2.0 / 4.0;
      U[7] = U[8] = 3.0 / 4.0;
    }

    std::ranges::fill_n(U.begin(), 3, 0.0);
    std::ranges::fill_n(std::ranges::reverse_view(U).begin(), 3, 1.0);
    return NURBSPatchData<1UL, 3UL, NurbsGridLinearAlgebraTraitsImpl>(knotVec, circleCPs, {2});
  }


  // Algo 8.1
  template <NurbsGridLinearAlgebra NurbsGridLinearAlgebraTraitsImpl = Dune::IGA::LinearAlgebraTraits<double, 2UL, 3UL>>
  auto makeSurfaceOfRevolution(const NURBSPatchData<1UL, 3UL, NurbsGridLinearAlgebraTraitsImpl>& generatrix,
                               const typename NurbsGridLinearAlgebraTraitsImpl::GlobalCoordinateType point,
                               const typename NurbsGridLinearAlgebraTraitsImpl::GlobalCoordinateType revolutionaxisI,
                               const typename NurbsGridLinearAlgebraTraitsImpl::GlobalCoordinateType::value_type revolutionAngle = 360.0) {
    const auto& genCP          = generatrix.controlPoints;
    using ScalarType           = typename NurbsGridLinearAlgebraTraitsImpl::GlobalCoordinateType::value_type;
    using ControlPoint         = typename NURBSPatchData<2UL, 3UL, NurbsGridLinearAlgebraTraitsImpl>::ControlPointType;
    using GlobalCoordinateType = typename NURBSPatchData<2UL, 3UL, NurbsGridLinearAlgebraTraitsImpl>::GlobalCoordinateType;
    const auto pi              = std::numbers::pi_v<typename NurbsGridLinearAlgebraTraitsImpl::GlobalCoordinateType::value_type>;

    const auto revolutionaxis = revolutionaxisI / revolutionaxisI.two_norm();
    auto newKnotsArray        = std::array<std::vector<double>, 2UL>();
    newKnotsArray[1]          = generatrix.knotSpans[0];
    auto& U                   = newKnotsArray[0];

    const int narcs = std::ceil(revolutionAngle / 90);
    U.resize(2 * (narcs + 2));
    if (revolutionAngle <= 90.0) {
    } else if (revolutionAngle <= 180.0) {
      U[3] = U[4] = 0.5;
    } else if (revolutionAngle <= 270.0) {
      U[3] = U[4] = 1.0 / 3.0;
      U[5] = U[6] = 2.0 / 3.0;
    } else {
      U[3] = U[4] = 0.25;
      U[5] = U[6] = 0.5;
      U[7] = U[8] = 0.75;
    }
    const ScalarType dtheta = revolutionAngle / narcs * pi / 180;
    std::ranges::fill_n(U.begin(), 3, 0.0);
    std::ranges::fill_n(std::ranges::reverse_view(U).begin(), 3, 1.0);

    typename NURBSPatchData<2UL, 3UL, NurbsGridLinearAlgebraTraitsImpl>::ControlPointNetType surfaceCP(2 * narcs + 1,
                                                                                                       generatrix.controlPoints.size()[0]);
    using std::cos;
    using std::sin;
    const ScalarType wm = cos(dtheta / 2.0);
    std::vector<ScalarType> cosines(narcs);
    std::vector<ScalarType> sines(narcs);
    ScalarType angle = 0.0;
    for (int i = 0; i < narcs; i++) {
      angle += dtheta;
      cosines[i] = cos(angle);
      sines[i]   = sin(angle);
    }
    ControlPoint PO = genCP.directGet(0);
    for (int j = 0; j < genCP.size()[0]; j++) {
      const GlobalCoordinateType Om = Impl::projectPointOntoLine(point, revolutionaxis, genCP.directGet(j).p);
      GlobalCoordinateType X        = genCP.directGet(j).p - Om;
      const ScalarType r            = X.two_norm();
      X /= r;
      const GlobalCoordinateType Y = cross(revolutionaxis, X);
      surfaceCP(0, j) = PO = genCP.directGet(j);
      GlobalCoordinateType TO    = Y;
      for (int index = 0, i = 0; i < narcs; ++i) {
        const GlobalCoordinateType P2 = Om + r * cosines[i] * X + r * sines[i] * Y;
        surfaceCP(index + 2, j) = {.p = P2, .w = genCP.directGet(j).w};

        const GlobalCoordinateType T2   = -sines[i] * X + cosines[i] * Y;
        surfaceCP(index + 1, j).p = Impl::Intersect3DLines(PO.p, TO, P2, T2);
        surfaceCP(index + 1, j).w = wm * genCP.directGet(j).w;
        index += 2;
        if (i < narcs - 1) {
          PO.p = P2;
          TO   = T2;
        }
      }
    }
    return NURBSPatchData<2UL, 3UL, NurbsGridLinearAlgebraTraitsImpl>(newKnotsArray, surfaceCP, {2, generatrix.order[0]});
  }
}  // namespace Dune::IGA