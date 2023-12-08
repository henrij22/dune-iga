
#pragma once

namespace Dune {
  namespace IGANEW {
    namespace DefaultTrim {
      // : TrimPatchEntitiy
      template <int codim_, int dim, class GridImp>
      class TrimmedParameterSpaceGridEntity {
        using ctype = typename GridImp::ctype;

        static constexpr int mydimension = dim;
        // [[nodiscard]] GeometryType type() const { return GeometryTypes::cube(mydimension); }
        using LocalCoordinate = FieldVector<ctype, mydimension>;

        using Trimmer           = typename GridImp::Trimmer;
        using GlobalIdSetIdType = typename Trimmer::TrimmerTraits::GlobalIdSetId;
        using EntityInfo        = typename Trimmer::TrimmerTraits::template Codim<codim_>::EntityInfo;
        using ElementTrimData   = typename Trimmer::ElementTrimData;
        using HostParameterSpaceGridEntity =
            typename Trimmer::TrimmerTraits::template Codim<codim_>::HostParameterSpaceGridEntity;
        using UntrimmedParameterSpaceGeometry =
            typename GridImp::Trimmer::TrimmerTraits::template Codim<codim_>::UntrimmedParameterSpaceGeometry;
        using TrimmedParameterSpaceGeometry =
            typename GridImp::Trimmer::TrimmerTraits::template Codim<codim_>::TrimmedParameterSpaceGeometry;
        // using LocalParameterSpaceGeometry=typename GridImp::Trimmer::TrimmerTraits::template
        // Codim<codim_>::LocalParameterSpaceGeometry;
        using LocalParameterSpaceGeometry =
            typename GridImp::Trimmer::TrimmerTraits::template Codim<codim_>::LocalParameterSpaceGeometry;
        using ParameterSpaceGridEntitySeed =
            typename GridImp::Trimmer::TrimmerTraits::template Codim<codim_>::ParameterSpaceGridEntitySeed;
        // using LocalParameterSpaceGeometry= typename Trimmer::TrimmerTraits::template
        // Codim<codim_>::LocalParameterSpaceGeometry;
       public:
        TrimmedParameterSpaceGridEntity()                                                      = default;
        TrimmedParameterSpaceGridEntity(const TrimmedParameterSpaceGridEntity& other) noexcept = default;
        TrimmedParameterSpaceGridEntity(TrimmedParameterSpaceGridEntity&& other) noexcept      = default;
        TrimmedParameterSpaceGridEntity& operator=(const TrimmedParameterSpaceGridEntity& other) = default;
        TrimmedParameterSpaceGridEntity& operator=(TrimmedParameterSpaceGridEntity&& other) noexcept = default;

        // Entity with codim 0 but trimmed thus needs trimdata
        template <typename = void>
        requires(codim_ == 0)
            TrimmedParameterSpaceGridEntity(const GridImp* grid, const HostParameterSpaceGridEntity& untrimmedElement,
                                            EntityInfo entInfo, const ElementTrimData& trimData)
            : grid_{grid},
              hostEntity_{untrimmedElement},
              trimmedlocalGeometry_{},
              entityInfo_{entInfo},
              trimData_{trimData} {
          assert(entityInfo_.lvl == untrimmedElement.level());
          isTrimmed = true;
          if (isTrimmed) trimmedlocalGeometry_ = std::make_optional<TrimmedParameterSpaceGeometry>();
        }

        // Entity untrimmed does not need trimdata but untrimmedElement
        TrimmedParameterSpaceGridEntity(const GridImp* grid, const HostParameterSpaceGridEntity& untrimmedElement,
                                        EntityInfo entInfo)
            : grid_{grid}, hostEntity_{untrimmedElement}, entityInfo_{entInfo} {
          assert(entityInfo_.lvl == untrimmedElement.level());
          isTrimmed = false;

          // DUNE_THROW(NotImplemented,"This constructor should accept a geometry object");
        }

        // Entity with codim!=0 but trimmed does need trimdata but no untrimmedElement
        template <typename = void>
        requires(codim_ != 0)
            TrimmedParameterSpaceGridEntity(const GridImp* grid, const ElementTrimData& trimData, EntityInfo entInfo)
            : grid_{grid}, trimData_{trimData}, entityInfo_{entInfo} {
          isTrimmed             = true;
          trimmedlocalGeometry_ = std::make_optional<TrimmedParameterSpaceGeometry>();
          // DUNE_THROW(NotImplemented,"This constructor should accept a geometry object");
        }

        auto& id() const { return entityInfo_.id; }
        template <typename = void>
        requires(codim_ == 0) auto& subId(int i, int codim) const {
          return grid_->trimmer().entityContainer_.subId(entityInfo_.id, i, codim);
        }

        HostParameterSpaceGridEntity getHostEntity() const {
          if (isTrimmed and codim_ != 0)
            DUNE_THROW(NotImplemented, "getHostEntity");
          else
            return hostEntity_;
        }

       private:
        EntityInfo entityInfo_;
        struct Empty {};
        HostParameterSpaceGridEntity hostEntity_;
        // The optional is only here since geometries are not default constructable
        std::optional<TrimmedParameterSpaceGeometry> trimmedlocalGeometry_;

        std::optional<std::reference_wrapper<const ElementTrimData>> trimData_;

       public:
        [[nodiscard]] bool operator==(const TrimmedParameterSpaceGridEntity& other) const {
          if constexpr (codim_ == 0)
            return hostEntity_ == other.hostEntity_;
          else
            return entityInfo_.id == other.entityInfo_.id;
        }

        //! returns true if father entity exists
        template <typename T = void>
        requires(codim_ == 0) [[nodiscard]] bool hasFather() const {
          if constexpr (codim_ == 0)
            if (entityInfo_.id.elementState == GlobalIdSetIdType::ElementState::full) return hostEntity_.hasFather();
          //@todo Trim this is crasy
          DUNE_THROW(NotImplemented, " hasFather");

          // return hostEntity_.hasFather();
        }

        //! Create EntitySeed
        [[nodiscard]] ParameterSpaceGridEntitySeed seed() const {
          DUNE_THROW(NotImplemented, " seed");
          if constexpr (codim_ == 0) return hostEntity_.seed();
          return {};
        }

        //! Level of this element
        [[nodiscard]] int level() const { return entityInfo_.lvl; }

        /** @brief The partition type for parallel computing */
        [[nodiscard]] PartitionType partitionType() const {
          //@todo Trim this is crasy
          if constexpr (codim_ == 0)
            if (entityInfo_.id.elementState == GlobalIdSetIdType::ElementState::full)
              return hostEntity_.partitionType();
          DUNE_THROW(NotImplemented, "partitionType not implemented for codim!=0 objects");
        }

        //! Geometry of this entity
        [[nodiscard]] LocalParameterSpaceGeometry geometry() const {
          //@todo Trim this is crasy
          // if(trimData_)
          //   return trimData_.template geometry<codim_>(localId_);
          // if constexpr (codim_==0) {
          //   if(id_.elementState==GlobalIdSetIdType::ElementState::full)
          //     return hostEntity_.geometry();
          //   else {
          //     DUNE_THROW(NotImplemented,"geometry not implemented for trimmed codim==0 objects");
          //
          //     return localGeometry_.value();
          //   }
          // }else {
          if (not isTrimmed) return hostEntity_.geometry();
          DUNE_THROW(NotImplemented, "geometry not implemented for trimmed codim!=0 objects");
          return trimmedlocalGeometry_.value();
          // }
        }

        /** @brief Return the number of subEntities of codimension codim.
         */
        [[nodiscard]] unsigned int subEntities(unsigned int codim) const {
          //@todo Trim this is crasy
          // if(trimData_)
          //   return trimData_. subEntities(codim,localId_);
          if constexpr (codim_ == 0) {
            if (entityInfo_.id.elementState == GlobalIdSetIdType::ElementState::full)
              return hostEntity_.subEntities(codim);
            else {
              DUNE_THROW(NotImplemented, "subEntities not implemented for codim==0 objects");
              return {};
            }
          } else {
            DUNE_THROW(NotImplemented, "subEntities not implemented for trimmed codim!=0 objects");
            return {};
          }
        }

        /** @brief Provide access to sub entity i of given codimension. Entities
         *  are numbered 0 ... subEntities(cc)-1
         */
        template <int cc>
        requires(codim_ == 0)
            [[nodiscard]] TrimmedParameterSpaceGridEntity<cc, mydimension, GridImp> subEntity(int i) const {
          // if(trimData_)
          //   return trimData_.template subEntity<codim_,cc>(i,localId_);
          // auto id = grid_->entityContainer().subId(id_,i,cc);
          if constexpr (cc == 0) return *this;
          return grid_->trimmer().entityContainer_.template entity<cc>(subId(i, cc));
          // {
          //   if constexpr (cc==0)
          //   {
          //     if (isTrimmed)
          //       return TrimmedParameterSpaceGridEntity<cc, mydimension, GridImp>(grid_, hostEntity_,
          //       entityInfo_,ElementTrimData()); // trimmed element with trimdata
          //     else {
          //       return TrimmedParameterSpaceGridEntity<cc, mydimension, GridImp>(grid_, hostEntity_, entityInfo_); //
          //       untrimmed element without trimdata
          //     }
          //   } else if ( not isTrimmed) {
          //     return TrimmedParameterSpaceGridEntity<cc, mydimension, GridImp>(grid_, hostEntity_.template
          //     subEntity<cc>(i), grid_->trimmer().entityContainer_.idToSubEntityInfoMap[cc+1].at( subId(i,cc))); //
          //     untrimmed subentity without trimdata
          //   }else
          //   {
          //     DUNE_THROW(Dune::NotImplemented, "trimmed subEntity can not be requested");
          //     return TrimmedParameterSpaceGridEntity<cc, mydimension, GridImp>(grid_,
          //     ElementTrimData(),grid_->trimmer().entityContainer_.idToSubEntityInfoMap[cc+1].at( subId(i,cc))); //
          //     Trimmed subentity without trimdata
          //
          //   }
          // }
        }

        //! First level intersection
        template <typename = void>
        requires(codim_ == 0) [[nodiscard]] decltype(auto) ilevelbegin() const {
          // if(trimData_)
          //   return trimData_.template ilevelbegin<codim_>(localId_);
          return hostEntity_.ilevelbegin();
        }

        //! Reference to one past the last neighbor
        template <typename = void>
        requires(codim_ == 0) decltype(auto) ilevelend() const { return hostEntity_.ilevelend(); }

        //! First leaf intersection
        template <typename = void>
        requires(codim_ == 0) decltype(auto) ileafbegin() const { return hostEntity_.ileafbegin(); }

        //! Reference to one past the last leaf intersection
        template <typename = void>
        requires(codim_ == 0) decltype(auto) ileafend() const { return hostEntity_.ileafend(); }

        //! returns true if Entity has NO children
        template <typename = void>
        requires(codim_ == 0) bool isLeaf() const { return hostEntity_.isLeaf(); }

        //! Inter-level access to father element on coarser grid.
        //! Assumes that meshes are nested.
        template <typename = void>
        requires(codim_ == 0) decltype(auto) father() const {
          assert(entityInfo_.fatherId.has_value());
          return grid_->trimmer().entityContainer_.template entity<0>(entityInfo_.fatherId.value());
          // return TrimmedParameterSpaceGridEntity(grid_, hostEntity_.father(),
          // grid_->trimmer().entityContainer_.idToElementInfoMap.at( entityInfo_.fatherId.value()));
        }

        /** @brief Location of this element relative to the reference element element of the father.
         * This is sufficient to interpolate all dofs in conforming case.
         * Nonconforming may require access to neighbors of father and
         * computations with local coordinates.
         * On the fly case is somewhat inefficient since dofs  are visited several times.
         * If we store interpolation matrices, this is tolerable. We assume that on-the-fly
         * implementation of numerical algorithms is only done for simple discretizations.
         * Assumes that meshes are nested.
         */
        template <typename = void>
        requires(codim_ == 0) decltype(auto) geometryInFather() const { return hostEntity_.geometryInFather(); }

        /** @brief Inter-level access to son elements on higher levels<=maxlevel.
         * This is provided for sparsely stored nested unstructured meshes.
         * Returns iterator to first son.
         */
        template <typename = void>
        requires(codim_ == 0) decltype(auto) hbegin(int maxLevel) const { return hostEntity_.hbegin(maxLevel); }

        //! Returns iterator to one past the last son
        template <typename = void>
        requires(codim_ == 0) decltype(auto) hend(int maxLevel) const { return hostEntity_.hend(maxLevel); }

        //! @todo Please doc me !
        template <typename = void>
        requires(codim_ == 0) bool wasRefined() const { return hostEntity_.wasRefined(); }

        //! @todo Please doc me !
        template <typename = void>
        requires(codim_ == 0)

            bool mightBeCoarsened() const {
          return hostEntity_.mightBeCoarsened();
        }

        const auto& hostEntity() const { return hostEntity_; }

        const GridImp* grid_;
        bool isTrimmed{};
      };

    }  // namespace DefaultTrim
  }    // namespace IGANEW
}  // namespace Dune
