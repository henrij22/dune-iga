// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/localfefunctions/expressions/linearStrainsExpr.hh>
#include <dune/localfefunctions/impl/standardLocalFunction.hh>
#include <dune/localfefunctions/manifolds/realTuple.hh>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <ikarus/assembler/simpleAssemblers.hh>
#include <ikarus/controlRoutines/loadControl.hh>
#include <ikarus/finiteElements/feBases/autodiffFE.hh>
#include <ikarus/finiteElements/feBases/powerBasisFE.hh>
#include <ikarus/finiteElements/feBases/scalarFE.hh>
#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/solver/nonLinearSolver/newtonRaphson.hh>
#include <ikarus/utils/algorithms.hh>
#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/drawing/griddrawer.hh>
#include <ikarus/utils/eigenDuneTransformations.hh>
#include <ikarus/utils/observer/controlVTKWriter.hh>
#include <ikarus/utils/observer/nonLinearSolverLogger.hh>

namespace Ikarus {

  template <typename Basis>
  class ReissnerMindlinPlate : public PowerBasisFE<typename Basis::FlatBasis> {
  public:
    using FlatBasis              = typename Basis::FlatBasis;
    using BaseDisp               = PowerBasisFE<FlatBasis>;  // Handles globalIndices function
    using GlobalIndex            = typename PowerBasisFE<FlatBasis>::GlobalIndex;
    using FERequirementType      = FErequirements<>;
    using ResultRequirementsType = ResultRequirements<FERequirementType>;
    using LocalView              = typename FlatBasis::LocalView;
    using GridView               = typename FlatBasis::GridView;

    template <typename VolumeLoad = std::nullptr_t>
    ReissnerMindlinPlate(Basis& basis, const typename LocalView::Element& element, const double p_Emodul,
                         const double p_nu, const double p_thickness, const VolumeLoad& p_volumeLoad = nullptr)
        : BaseDisp(basis.flat(), element),
          Emodul{p_Emodul},
          nu{p_nu},
          thickness{p_thickness},
          volumeLoad_(p_volumeLoad) {
      this->localView().bind(element);
      const int order = 2 * (this->localView().tree().child(0).finiteElement().localBasis().order());
      element.impl().fillQuadratureRule(rule, order);
    }

    static Eigen::Matrix<double, 5, 5> constitutiveMatrix(double Emod, double p_nu, double p_thickness) {
      const double factor = Emod * Dune::power(p_thickness, 3) / (12.0 * (1.0 - p_nu * p_nu));
      Eigen::Matrix<double, 5, 5> D;
      D.setZero();
      D(0, 0) = D(1, 1) = 1;
      D(0, 1) = D(1, 0) = p_nu;
      D(2, 2)           = (1 - p_nu) / 2.0;
      D *= factor;
      const double shear_term = (5.0 / 6.0) * p_thickness * Emod / (2.0 * (1.0 + p_nu));
      D(3, 3) = D(4, 4) = shear_term;
      return D;
    }

    using Traits = TraitsFromLocalView<LocalView>;

  public:
    double calculateScalar(const FERequirementType& par) const {
      std::cerr << "Returning zero energy via calculateScalar(feRequirement)" << std::endl;
      return 0.0;
    }

    void calculateVector(const FERequirementType& par, typename Traits::VectorType& g) const {
      const auto& lambda = par.getParameter(Ikarus::FEParameter::loadfactor);
      const auto D       = constitutiveMatrix(Emodul, nu, thickness);
      using namespace Dune::DerivativeDirections;
      auto& ele = this->localView().element();
      auto& fe  = this->localView().tree().child(0).finiteElement();

      const auto& localBasis = fe.localBasis();
      const auto geo         = this->localView().element().geometry();

      g.template setZero(this->localView().size());
      if (volumeLoad_) {
        for (const auto& gp : rule) {
          Eigen::Vector<double, 3> fext = (*volumeLoad_)(geo.global(gp.position()), lambda);
          const auto Jinv               = toEigen(geo.jacobianInverseTransposed(gp.position())).transpose().eval();
          const double intElement       = geo.integrationElement(gp.position()) * gp.weight();
          std::vector<Dune::FieldVector<double, 1>> shapeFunctionValues;
          localBasis.evaluateFunction(gp.position(), shapeFunctionValues);
          Eigen::MatrixXd N;
          N.setZero(3, this->localView().size());
          for (size_t nn = 0; nn < shapeFunctionValues.size(); ++nn) {
            N(0, 3 * nn)     = shapeFunctionValues[nn];
            N(1, 3 * nn + 1) = shapeFunctionValues[nn];
            N(2, 3 * nn + 2) = shapeFunctionValues[nn];
          }
          g -= N.transpose() * fext * intElement;
        }
      }
    }

    void calculateMatrix(const FERequirementType& par, typename Traits::MatrixType& h) const {
      const auto& lambda = par.getParameter(Ikarus::FEParameter::loadfactor);
      const auto D       = constitutiveMatrix(Emodul, nu, thickness);
      using namespace Dune::DerivativeDirections;
      auto& ele = this->localView().element();
      auto& fe  = this->localView().tree().child(0).finiteElement();

      const auto& localBasis = fe.localBasis();
      const auto geo         = this->localView().element().geometry();
      h.template setZero(this->localView().size(), this->localView().size());

      for (const auto& gp : rule) {
        const auto Jinv         = geo.jacobianInverseTransposed(gp.position());
        const double intElement = geo.integrationElement(gp.position()) * gp.weight();
        std::vector<Dune::FieldMatrix<double, 1, 2>> referenceGradients;
        localBasis.evaluateJacobian(gp.position(), referenceGradients);
        std::vector<Dune::FieldVector<double, 2>> gradients(referenceGradients.size());

        for (size_t i = 0; i < gradients.size(); i++)
          Jinv.mv(referenceGradients[i][0], gradients[i]);

        std::vector<Dune::FieldVector<double, 1>> shapeFunctionValues;
        localBasis.evaluateFunction(gp.position(), shapeFunctionValues);

        Eigen::VectorXd dNdx = Eigen::VectorXd::Zero(shapeFunctionValues.size());
        Eigen::VectorXd dNdy = Eigen::VectorXd::Zero(shapeFunctionValues.size());
        for (size_t i = 0; i < shapeFunctionValues.size(); i++) {
          dNdx[i] = gradients[i][0];
          dNdy[i] = gradients[i][1];
        }

        Eigen::MatrixXd bop;
        bop.setZero(5, this->localView().size());
        for (auto i = 0U; i < shapeFunctionValues.size(); ++i) {
          bop(0, 3 * i + 2) = dNdx(i);

          bop(1, 3 * i + 1) = -dNdy(i);

          bop(2, 3 * i + 2) = dNdy(i);
          bop(2, 3 * i + 1) = -dNdx(i);

          bop(3, 3 * i)     = dNdx(i);
          bop(3, 3 * i + 2) = shapeFunctionValues[i];

          bop(4, 3 * i)     = dNdy(i);
          bop(4, 3 * i + 1) = -shapeFunctionValues[i];
        }
        h += bop.transpose() * D * bop * intElement;
      }
    }

    void calculateAt(const ResultRequirementsType& req, const Dune::FieldVector<double, Traits::mydim>& local,
                     ResultTypeMap<double>& result) const {
      using namespace Dune::Indices;
      const auto& disp       = req.getGlobalSolution(Ikarus::FESolutions::displacement);
      const auto D           = constitutiveMatrix(Emodul, nu, thickness);
      auto& fe               = this->localView().tree().child(0).finiteElement();
      const auto& localBasis = fe.localBasis();
      const auto geo         = this->localView().element().geometry();
      Eigen::VectorXd local_disp;
      local_disp.setZero(this->localView().size());

      int disp_counter = 0;
      for (size_t i = 0; i < fe.size(); ++i)
        for (size_t j = 0; j < 3; ++j) {
          auto globalIndex         = this->localView().index(this->localView().tree().child(j).localIndex(i));
          local_disp[disp_counter] = disp[globalIndex];
          disp_counter++;
        }

      const auto Jinv = geo.jacobianInverseTransposed(local);
      std::vector<Dune::FieldMatrix<double, 1, 2>> referenceGradients;
      localBasis.evaluateJacobian(local, referenceGradients);
      std::vector<Dune::FieldVector<double, 2>> gradients(referenceGradients.size());

      for (size_t i = 0; i < gradients.size(); i++)
        Jinv.mv(referenceGradients[i][0], gradients[i]);

      std::vector<Dune::FieldVector<double, 1>> shapeFunctionValues;
      localBasis.evaluateFunction(local, shapeFunctionValues);

      Eigen::VectorXd dNdx = Eigen::VectorXd::Zero(shapeFunctionValues.size());
      Eigen::VectorXd dNdy = Eigen::VectorXd::Zero(shapeFunctionValues.size());
      for (size_t i = 0; i < shapeFunctionValues.size(); i++) {
        dNdx[i] = gradients[i][0];
        dNdy[i] = gradients[i][1];
      }

      Eigen::MatrixXd bop;
      bop.setZero(5, this->localView().size());

      for (auto i = 0U; i < shapeFunctionValues.size(); ++i) {
        bop(0, 3 * i + 2) = dNdx(i);

        bop(1, 3 * i + 1) = -dNdy(i);

        bop(2, 3 * i + 2) = dNdy(i);
        bop(2, 3 * i + 1) = -dNdx(i);

        bop(3, 3 * i)     = dNdx(i);
        bop(3, 3 * i + 2) = shapeFunctionValues[i];

        bop(4, 3 * i)     = dNdy(i);
        bop(4, 3 * i + 1) = -shapeFunctionValues[i];
      }

      Eigen::Vector<double, 5> req_res;
      req_res.setZero();
      req_res = D * bop * local_disp;

      Eigen::Matrix<double, 3, 3> sent_res;
      sent_res(0, 0) = req_res[0];                              // moment m_xx
      sent_res(1, 1) = req_res[1];                              // moment m_yy
      sent_res(0, 1) = sent_res(1, 0) = req_res[2];   // moment m_xy
      sent_res(0, 2) = sent_res(2, 0) = req_res[3];   // shear force vx
      sent_res(1, 2) = sent_res(2, 1) = req_res[4];   // shear force vy

      typename ResultTypeMap<double>::ResultArray resv;
      if (req.isResultRequested(ResultType::cauchyStress)) {
        resv.resize(3, 3);
        resv = sent_res;
        result.insertOrAssignResult(ResultType::cauchyStress, resv);
      }
    }

    Dune::QuadratureRule<double, 2> rule;
    double Emodul;
    double nu;
    double thickness;
    std::optional<
        std::function<Eigen::Vector<double, 3>(const Dune::FieldVector<double, Traits::worlddim>&, const double&)>>
        volumeLoad_;
  };

}  // namespace Ikarus
