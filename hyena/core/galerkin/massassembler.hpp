// Copyright (C) 2009-2010 Matthias Messner, Michael Messner, Franz
// Rammerstorfer, Peter Urthaler
// 
// This file is part of HyENA - a C++ boundary element methods library.
// 
// HyENA is free software: you can redistribute it and/or modify it under the
// terms of the GNU Lesser Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
// 
// HyENA is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU Lesser Public License for more
// details.
// 
// You should have received a copy of the GNU Lesser Public License along with
// HyENA. If not, see <http://www.gnu.org/licenses/>.
#ifndef massassembler_hpp
#define massassembler_hpp

// stl includes 
#include <vector>

// boost
#include <boost/utility.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

// own includes
#include "hyena/core/common/storageadapter.hpp"
#include "hyena/core/common/tags.H"
#include "hyena/core/galerkin/identitykernel.hpp"
#include "hyena/core/galerkin/singleintegrator3d.hpp"


namespace hyena
{

template<typename MATRIX>
class StorageAdapterLHS;
  
template<typename MATRIX>
class StorageAdapterSparseLHS;

template<typename VECTOR>
class StorageAdapterRHS;

  template<typename DTX, 
           typename DTY,
           typename X_ASSEMBLE_TAG = IDX_tag,
           typename Y_ASSEMBLE_TAG = IDX_tag,
           typename STORAGE_TAG = DENSE_tag>
  class MassAssembler 
    : boost::noncopyable
  {

    // compile time
    typedef typename DTX::gdof_type x_gdof_type;
    typedef typename DTY::gdof_type y_gdof_type;

    static const ELEMENT_SHAPE shape     = DTX::shape;
    static const PROBLEM       problem   = DTX::problem;
    static const APPROXIMATION approx_x  = DTX::field_approximation;
    static const APPROXIMATION approx_y  = DTY::field_approximation;
 
    typedef ReferenceElement<shape,approx_x> x_type;
    typedef ReferenceElement<shape,approx_y> y_type;
 
    typedef IdentityKernel<DTX, DTY>                       kernel_type;
    typedef SingleIntegrator3d<kernel_type>            integrator_type;
    typedef typename kernel_type::kernel_block_type  kernel_block_type;


    typedef typename DTX::const_ldof_array_type x_ldof_array_type;
    typedef typename DTY::const_ldof_array_type y_ldof_array_type;
  
    typedef typename DTX::superelement_type                      x_se_type;
    typedef typename DTY::superelement_type                      y_se_type;
  
    typedef typename DTX::superelement_extractor_type  x_se_extractor_type;
    typedef typename DTY::superelement_extractor_type  y_se_extractor_type;
  
    typedef typename DTX::const_superelement_array_type  x_se_array_type;
    typedef typename DTY::const_superelement_array_type  y_se_array_type;

    static const unsigned int  num_ldofs_per_gdof_ = 
      ProblemTraits<problem>::num_ldofs_per_gdof;
    static const unsigned int x_num_gdofs_per_element_ = 
      x_type::num_gdofs_per_element;
    static const unsigned int y_num_gdofs_per_element_ = 
      y_type::num_gdofs_per_element;

    const x_se_extractor_type x_se_extractor_;
    const y_se_extractor_type y_se_extractor_;
    const kernel_type kernel_;
    const integrator_type integrate_;
    const X_ASSEMBLE_TAG x_tag_;
    const Y_ASSEMBLE_TAG y_tag_;
    const STORAGE_TAG storage_tag_;



  public:
  
    // ctor
    MassAssembler(const unsigned int iorder)
      : x_se_extractor_(),
        y_se_extractor_(),
        kernel_(),
        integrate_(kernel_, iorder),
        x_tag_(X_ASSEMBLE_TAG() ),
        y_tag_(Y_ASSEMBLE_TAG() ),
        storage_tag_(STORAGE_TAG() )
    {}
    

  
    //! assemble LHS 
    template<typename MATRIX>
    void operator()(const x_ldof_array_type& x_ldofs, 
                    const y_ldof_array_type& y_ldofs,
                    const double scale,
                    MATRIX& matrix) const
    {
      assemble_dispatch(x_ldofs, y_ldofs, scale, matrix, storage_tag_);
    }
    


    //! assemble LHS symmetric
    template<typename MATRIX>
    void operator()(const x_ldof_array_type& ldofs, 
                    const double scale,
                    MATRIX& matrix) const 
    {
      assemble_dispatch(ldofs, ldofs, scale, matrix, storage_tag_);
    }


    
    //! assemble RHS
    template<typename VECTOR> 
    void operator()(const x_ldof_array_type& x_ldofs, 
                    const y_ldof_array_type& y_ldofs,
                    const VECTOR& data, 
                    const double scale,
                    VECTOR& result) const
    {
      // functor for storage
      StorageAdapterRHS<VECTOR> store(data, 1., result);
      assemble(x_ldofs, y_ldofs, scale, store);
      // finalize result vector, needed for parallel computation.
      store.finalize();
    }

    
  private:
  
    template<typename DATA>
    void assemble_dispatch (const x_ldof_array_type& x_ldofs, 
                            const y_ldof_array_type& y_ldofs,
                            const double scale,
                            DATA& dat,
                            const DENSE_tag&) const
    {
      // functor for storage
      StorageAdapterLHS<DATA> store(dat);
      assemble(x_ldofs, y_ldofs, scale, store);
    }

    template<typename DATA>
    void assemble_dispatch (const x_ldof_array_type& x_ldofs, 
                            const y_ldof_array_type& y_ldofs,
                            const double scale,
                            DATA& dat,
                            const SPARSE_tag&) const
    {
      // functor for storage
      StorageAdapterSparseLHS<typename DATA::value_type> store(dat);
      assemble(x_ldofs, y_ldofs, scale, store);
    }



    template<typename STORAGE>
    void assemble(const x_ldof_array_type& x_ldofs, 
                  const y_ldof_array_type& y_ldofs,
                  const double scale,
                  STORAGE& store) const
    {
      // abort in case the boundary does not exist
      if (x_ldofs.size()==0 or y_ldofs.size()==0)
        return;

      // get flags
      const LDOF_TYPE x_ldof_flag = x_ldofs[0]->getType();
      const LDOF_TYPE y_ldof_flag = y_ldofs[0]->getType();
    
      // extract needed superelement arrays
      x_se_array_type x_se;
      y_se_array_type y_se;
      x_se_extractor_(x_ldofs, x_se);
      y_se_extractor_(y_ldofs, y_se);
    
#ifdef _MSC_VER && ASSEMBLE_GALERKIN_OMP
#pragma omp parallel for schedule(guided)
      for (long it_p=0; it_p<x_se.size(); ++it_p) { typename x_se_array_type::const_iterator iter_x=x_se.begin()+it_p;
#else
      // start loop over row support
      typename x_se_array_type::const_iterator iter_x;
#ifdef ASSEMBLE_GALERKIN_OMP
#pragma omp parallel for schedule(guided)
#endif
      for (iter_x=x_se.begin(); iter_x<x_se.end(); ++iter_x) {
#endif
        // instantiate block to be reused by the integrator
        kernel_block_type block;
        // start loop over column support
        typename y_se_array_type::const_iterator iter_y;
        for (iter_y=y_se.begin();iter_y<y_se.end(); ++iter_y) {
          // skip the loop if elements do not coincide
          if( (*iter_x)->getID() != (*iter_y)->getID() ) 
            continue;
          // set all block entries to zero
          for(unsigned int i=0;i<block.getNumEntries(); ++i) 
            block[i].zeros();
          // integrate
          integrate_(*iter_x, *iter_y, block);
          // store
          storeBlock(*iter_x, *iter_y, x_ldof_flag, y_ldof_flag, 
                     scale, block, store, x_tag_, y_tag_);
        }//se_y
      }//se_x
    }



  private:

    // store block with respect to ldof_type aka. boundary using idx
    template<typename BLOCK,
             typename STORAGE>
    void storeBlock(const x_se_type *const x_se,
                    const y_se_type *const y_se,
                    const LDOF_TYPE x_ldof_flag,
                    const LDOF_TYPE y_ldof_flag,
                    const double scale,
                    const BLOCK& block,
                    STORAGE& store,
                    const IDX_tag&,
                    const IDX_tag&) const
    {
      // loop over all shape_functions on tau_x
      for(unsigned int ngdx=0;ngdx<x_num_gdofs_per_element_; ++ngdx){
        const x_gdof_type *const x_gdof = x_se->getGDof(ngdx);
        // loop over all shape_functions on tau_y
        for(unsigned int ngdy=0; ngdy<y_num_gdofs_per_element_; ++ngdy){
          const y_gdof_type *const y_gdof = y_se->getGDof(ngdy);
          // loop over all FUNDSOL::num_ldofs_per_gdof on tau_x
          for(unsigned int nldx=0; nldx<num_ldofs_per_gdof_; ++nldx){
            // get index i with boundary distinction
            if (x_gdof->getLDof(nldx)->getType()!=x_ldof_flag)
              continue;
            const unsigned int i = x_gdof->getLDof(nldx)->getIDX();
            // loop over all FUNDSOL::num_ldofs_per_gdof on tau_y
            for(unsigned int nldy=0; nldy<num_ldofs_per_gdof_; ++nldy){
              // get index j with boundary distinction
              if (y_gdof->getLDof(nldy)->getType()!=y_ldof_flag)
                continue;
              const unsigned int j = y_gdof->getLDof(nldy)->getIDX();
              // store
              store( i, j,block(ngdx, ngdy)(nldx, nldy)*scale );
            }//y_ldof
          }//x_lodf
        }//y_gdof
      }//x_gdof
    }


    template<typename BLOCK,
             typename STORAGE>
    void storeBlock(const x_se_type *const x_se,
                    const y_se_type *const y_se,
                    const LDOF_TYPE x_ldof_flag,
                    const LDOF_TYPE y_ldof_flag,
                    const double scale,
                    const BLOCK& block,
                    STORAGE& store,
                    const IDX_tag&,
                    const ID_tag&) const
    {
      // loop over all shape_functions on tau_x
      for(unsigned int ngdx=0;ngdx<x_num_gdofs_per_element_; ++ngdx){
        const x_gdof_type *const x_gdof = x_se->getGDof(ngdx);
        // loop over all shape_functions on tau_y
        for(unsigned int ngdy=0; ngdy<y_num_gdofs_per_element_; ++ngdy){
          const y_gdof_type *const y_gdof = y_se->getGDof(ngdy);
          // loop over all FUNDSOL::num_ldofs_per_gdof on tau_x
          for(unsigned int nldx=0; nldx<num_ldofs_per_gdof_; ++nldx){
            // get index i with boundary distinction
            if (x_gdof->getLDof(nldx)->getType()!=x_ldof_flag)
              continue;
            const unsigned int i = x_gdof->getLDof(nldx)->getIDX();
            // loop over all FUNDSOL::num_ldofs_per_gdof on tau_y
            for(unsigned int nldy=0; nldy<num_ldofs_per_gdof_; ++nldy){
              // get index j without
              const unsigned int j = y_gdof->getLDof(nldy)->getID();
              // store
              store( i, j,block(ngdx, ngdy)(nldx, nldy)*scale );
            }//y_ldof
          }//x_lodf
        }//y_gdof
      }//x_gdof
    }





    template<typename BLOCK,
             typename STORAGE>
    void storeBlock(const x_se_type *const x_se,
                    const y_se_type *const y_se,
                    const LDOF_TYPE x_ldof_flag,
                    const LDOF_TYPE y_ldof_flag,
                    const double scale,
                    const BLOCK& block,
                    STORAGE& store,
                    const ID_tag&,
                    const IDX_tag&) const
    {
      // loop over all shape_functions on tau_x
      for(unsigned int ngdx=0;ngdx<x_num_gdofs_per_element_; ++ngdx){
        const x_gdof_type *const x_gdof = x_se->getGDof(ngdx);
        // loop over all shape_functions on tau_y
        for(unsigned int ngdy=0; ngdy<y_num_gdofs_per_element_; ++ngdy){
          const y_gdof_type *const y_gdof = y_se->getGDof(ngdy);
          // loop over all FUNDSOL::num_ldofs_per_gdof on tau_x
          for(unsigned int nldx=0; nldx<num_ldofs_per_gdof_; ++nldx){
            // get index i without
            const unsigned int i = x_gdof->getLDof(nldx)->getID();
            // loop over all FUNDSOL::num_ldofs_per_gdof on tau_y
            for(unsigned int nldy=0; nldy<num_ldofs_per_gdof_; ++nldy){
              // get index j with boundary distinction
              if (y_gdof->getLDof(nldy)->getType()!=y_ldof_flag)
                continue;
              const unsigned int j = y_gdof->getLDof(nldy)->getIDX();
              // store
              store( i, j,block(ngdx, ngdy)(nldx, nldy)*scale );
            }//y_ldof
          }//x_lodf
        }//y_gdof
      }//x_gdof
    }







    // store block irrespectible of ldof_type  using id (whole boundary)
    template<typename BLOCK,
             typename STORAGE>
    void storeBlock(const x_se_type *const x_se,
                    const y_se_type *const y_se,
                    const LDOF_TYPE x_ldof_flag,
                    const LDOF_TYPE y_ldof_flag,
                    const double scale, 
                    const BLOCK& block,
                    STORAGE& store,
                    const ID_tag&,
                    const ID_tag&) const
    {
      // loop over all shape_functions on tau_x
      for(unsigned int ngdx=0;ngdx<x_num_gdofs_per_element_; ++ngdx){
        const x_gdof_type *const x_gdof = x_se->getGDof(ngdx);
        // loop over all shape_functions on tau_y
        for(unsigned int ngdy=0; ngdy<y_num_gdofs_per_element_; ++ngdy){
          const y_gdof_type *const y_gdof = y_se->getGDof(ngdy);
          // loop over all FUNDSOL::num_ldofs_per_gdof on tau_x
          for(unsigned int nldx=0; nldx<num_ldofs_per_gdof_; ++nldx){
            // get index i without boundary distinction
            const unsigned int i = x_gdof->getLDof(nldx)->getID();
            // loop over all FUNDSOL::num_ldofs_per_gdof on tau_y
            for(unsigned int nldy=0; nldy<num_ldofs_per_gdof_; ++nldy){
              // get index j without boundary distinction
              const unsigned int j = y_gdof->getLDof(nldy)->getID();
              // store
              store( i, j, block(ngdx, ngdy)(nldx, nldy)*scale );
            }//y_ldof
          }//x_lodf
        }//y_gdof
      }//x_gdof
    }


    

    // store block with respect to ldof_type aka. boundary using vidx
    template<typename BLOCK,
             typename STORAGE>
    void storeBlock(const x_se_type *const x_se,
                    const y_se_type *const y_se,
                    const LDOF_TYPE x_ldof_flag,
                    const LDOF_TYPE y_ldof_flag,
                    const double scale,
                    const BLOCK& block,
                    STORAGE& store,
                    const VIDX_tag&,
                    const VIDX_tag&) const
    {
      // loop over all shape_functions on tau_x
      for(unsigned int ngdx=0;ngdx<x_num_gdofs_per_element_; ++ngdx){
        const x_gdof_type *const x_gdof = x_se->getGDof(ngdx);
        // loop over all shape_functions on tau_y
        for(unsigned int ngdy=0; ngdy<y_num_gdofs_per_element_; ++ngdy){
          const y_gdof_type *const y_gdof = y_se->getGDof(ngdy);
          // loop over all FUNDSOL::num_ldofs_per_gdof on tau_x
          for(unsigned int nldx=0; nldx<num_ldofs_per_gdof_; ++nldx){
            // get index i with boundary distinction
            if (x_gdof->getLDof(nldx)->getType()!=x_ldof_flag)
              continue;
            const unsigned int i = x_gdof->getLDof(nldx)->vidx;
            // loop over all FUNDSOL::num_ldofs_per_gdof on tau_y
            for(unsigned int nldy=0; nldy<num_ldofs_per_gdof_; ++nldy){
              // get index j with boundary distinction
              if (y_gdof->getLDof(nldy)->getType()!=y_ldof_flag)
                continue;
              const unsigned int j = y_gdof->getLDof(nldy)->vidx;
              // store
              store( i, j,block(ngdx, ngdy)(nldx, nldy)*scale );
            }//y_ldof
          }//x_lodf
        }//y_gdof
      }//x_gdof
    }




  };

} // end namespace hyena

#endif

               
