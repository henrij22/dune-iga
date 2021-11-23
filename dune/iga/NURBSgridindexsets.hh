// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#pragma once

#include <dune/grid/common/indexidset.hh>
#include <dune/iga/NURBSgridentity.hh>
#include <span>
#include <map>

namespace Dune::IGA
  {
    template<class GridViewImp>
    class NURBSGridLeafIndexSet
    {
    public:
      using Types = std::array<GeometryType ,1>;
      using IndexType = unsigned int;

      //! constructor
      explicit NURBSGridLeafIndexSet (const GridViewImp& g) : gridView_(&g)
      {
          for(int id= 0;auto& entTypes: types_) {
              entTypes = {GeometryTypes::cube(GridViewImp::dimension - id)};
              ++id;
          }
      }

      template<class Entity>
      bool contains(const Entity& e) const
      {
        return gridView_->contains(e);
      }

      //! get index of an entity, need to change the type to a property from GridImp
      template<int codim>
      int index (const NURBSGridEntity<codim, GridViewImp> & e) const
      {
        return e.getIndex();
      }

        template<int codimElement>
        int subIndex (const NURBSGridEntity<codimElement, GridViewImp> & e,int i, unsigned int codim) const
        {
          if (codimElement==0 && NURBSGridEntity<codimElement, GridViewImp>::mydim==codim)
            return gridView_->NURBSpatch_->getGlobalVertexIndexFromElementIndex(e.getIndex(),i);
          else if (i==0 && codim == 0 && codimElement ==0)
            return this->index(e);
          else if (codim == 1 && codimElement ==0 || (codim == 2 && codimElement ==0 && GridViewImp::dimension==3))
            return gridView_->NURBSpatch_->getGlobalEdgeIndexFromElementIndex(e.getIndex(),i);
          else
            throw std::logic_error("subIndex only defined from element to vertices and edges");
        }



      auto& types(int codim) const
      {
          return types_[codim];
      }


      auto size(int codim) const
      {
          return gridView_->size(codim);
      }

        auto size(const GeometryType& gt) const
        {
          const int codim = GridViewImp::dimension - gt.dim();
            return gridView_->size(codim);
        }



    private:
      std::array<std::array<GeometryType ,1>,GridViewImp::dimension+1> types_;
      std::array<std::map<int ,int>,GridViewImp::dimension+1> indices;
      const GridViewImp* gridView_;
    };
  }
