
#pragma once



    namespace Dune::IGANEW::DefaultTrim {

      template <typename GridImp>
      struct PatchTrimDataImpl {
        using TrimmingCurve = typename GridImp::GridFamily::TrimmerTraits::TrimmingCurve;

        void insertTrimCurve(const TrimmingCurve& curve) { curves_.push_back(curve); }
        const auto& curves() const { return curves_; }
        std::vector<TrimmingCurve> curves_;

        //???
      };

    } // namespace Dune::IGANEW::DefaultTrim


