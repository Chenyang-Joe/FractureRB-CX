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


/**
 * @file    galerkinassembler.hpp
 * @author  MiM, MaM, PU, TT, Rf
 * @date    created:     16.08.11
 *          last change: 16.08.11
 */
#ifndef galerkinassembler_hpp
#define galerkinassembler_hpp

// stl includes 
#include <vector>
#include <algorithm>

// boost
#include <boost/utility.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

// own includes
#include "hyena/core/common/storageadapter.hpp"
#include "hyena/core/utilities/weightedifftblocks.hpp"
#include "hyena/core/common/tags.H"


namespace hyena
{

template<typename MATRIX>
class StorageAdapterLHS;

template<typename VECTOR>
class StorageAdapterRHS;

template<typename MATRIX>
class StorageAdapterArrayLHS;

template<typename VECTOR>
class StorageAdapterConvolutionRHS;

template<typename VECTOR>
class StorageAdapterArrayRHS;

template<typename IN_BLOCK>
class WeightedIFFTBlocks;

// TODO FIX ME
  /**
   * This template function gets the row and column support and uses the
   * integrator to assemble the matrix. The row and column support is provided
   * by superelement iterators that have the ++ operator defined. The instance
   * of the integrator type must have the () operator defined. The matrix must
   * have the access to be defined in the way /a matrix(row index, column
   * index).
   *
   * @ingroup galerkin
   *
   * @tparam TAU_X      type of superelement iterator providing row support
   * @tparam TAU_Y      type of superelement iterator providing column support
   * @tparam INTEGRATOR type of integrator
   * @tparam MATRIX     type of matrix to be assembled
   *
   * @param[in]  start_x    start iterator providing row support
   * @param[in]  end_iter_x end iterator providing row support
   * @param[in]  start_y    start iterator providing column support
   * @param[in]  end_iter_y end iterator providing column support
   * @param[in]  integrate  functor that integrates over computational support
   * @param[out] matrix     to be assembled
   * @param[in]  WITH_JUMP  true is with and false without 1/2 identity
   */
  template<typename KERNEL,
           typename X_ASSEMBLE_TAG = IDX_tag,
           typename Y_ASSEMBLE_TAG = IDX_tag>
  class GalerkinAssembler
  {
  private:
    typedef KERNEL                                kernel_type;
    typedef GalerkinIntegrator3d<kernel_type> integrator_type;

    typedef typename kernel_type::dof_traits_x_type           dtx_type;
    typedef typename kernel_type::dof_traits_y_type           dty_type;
  
    typedef typename kernel_type::value_type                value_type;
    typedef typename kernel_type::kernel_block_type  kernel_block_type;
    typedef typename kernel_type::x_type                        x_type;
    typedef typename kernel_type::y_type                        y_type;
  
    typedef typename dtx_type::superelement_extractor_type  x_se_extractor_type;
    typedef typename dtx_type::gdof_type                            x_gdof_type;
    typedef typename dtx_type::superelement_type                      x_se_type;
    typedef typename dtx_type::const_superelement_array_type    x_se_array_type;
    typedef typename dtx_type::const_ldof_array_type          x_ldof_array_type;
 
    typedef typename dty_type::superelement_extractor_type  y_se_extractor_type;
    typedef typename dty_type::gdof_type                            y_gdof_type;
    typedef typename dty_type::superelement_type                      y_se_type;
    typedef typename dty_type::const_superelement_array_type    y_se_array_type;
    typedef typename dty_type::const_ldof_array_type          y_ldof_array_type;
  
  private:
  
    static const unsigned int num_ldofs_per_gdof_ = 
      ProblemTraits<kernel_type::problem>::num_ldofs_per_gdof;
    static const unsigned int x_num_gdofs_per_element_ = 
      x_type::num_gdofs_per_element;
    static const unsigned int y_num_gdofs_per_element_ = 
      y_type::num_gdofs_per_element;
  
    const integrator_type integrate_;
    const x_se_extractor_type x_se_extractor_;
    const y_se_extractor_type y_se_extractor_;
    const X_ASSEMBLE_TAG x_tag_;
    const Y_ASSEMBLE_TAG y_tag_;


  public:

    // ctor
    template<typename FSOL>
    GalerkinAssembler(FSOL& fsol, const unsigned int iorder)
      : integrate_(fsol, iorder),
        x_se_extractor_(),
        y_se_extractor_(), 
        x_tag_(X_ASSEMBLE_TAG() ), 
        y_tag_(Y_ASSEMBLE_TAG() )
    {}
  


    //! assemble LHS
    template<typename MATRIX>
    void operator()(const x_ldof_array_type& x_ldofs, 
                    const y_ldof_array_type& y_ldofs,
                    MATRIX& matrix, int x_store_offset=0, int y_store_offset=0,
                    int x_crop=0, int y_crop=0,int x_clip=-1, int y_clip=-1) const // mod. DH, 16.12.2013: offset position in storage matrix
    {
      // functor for storage
      StorageAdapterLHS<MATRIX> store(matrix,x_store_offset,y_store_offset,x_crop,y_crop,x_clip,y_clip); // mod. DH, 16.12.2013: offset position in storage matrix
      assemble(x_ldofs, y_ldofs, store);
    }



    //! assemble LHS symmetric
    template<typename MATRIX>
    void operator()(const x_ldof_array_type& ldofs, 
                    MATRIX& matrix) const
    {
      // functor for storage
      StorageAdapterLHS<MATRIX> store(matrix);
      assembleSymmetric(ldofs, store);
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
      StorageAdapterRHS<VECTOR> store(data, scale, result);
      assemble(x_ldofs, y_ldofs, store);
      // finalize result vector, needed for parallel computation.
      store.finalize();
    }


    //! assemble array LHS
    template<typename MATRIX>
    void operator()(const x_ldof_array_type& x_ldofs, 
                    const y_ldof_array_type& y_ldofs,
                    const std::vector<std::complex<double> >& laplace_params,
                    const unsigned int relevant_steps,
                    std::vector<MATRIX>& matrices) const
    {
      // functor for storage
      StorageAdapterArrayLHS<MATRIX> store(matrices);
      // functor for wifft
      WeightedIFFTBlocks<kernel_block_type> transform(laplace_params.size() );
      assemble(x_ldofs, y_ldofs, laplace_params, relevant_steps, 
               store, transform);
    }


    //! assemble array LHS symmetric
    template<typename MATRIX>
    void operator()(const x_ldof_array_type& ldofs,
                    const std::vector<std::complex<double> >& laplace_params,
                    const unsigned int relevant_steps,
                    std::vector<MATRIX>& matrices) const
    {
      // functor for storage
      StorageAdapterArrayLHS<MATRIX> store(matrices);
      // functor for wifft
      WeightedIFFTBlocks<kernel_block_type> transform(laplace_params.size() );
      assembleSymmetric(ldofs,laplace_params,relevant_steps,store,transform);
    }
  

    //! assemble array RHS with convolution
    template<typename VECTOR>
    void operator()(const x_ldof_array_type& x_ldofs, 
                    const y_ldof_array_type& y_ldofs,
                    const std::vector<std::complex<double> >& laplace_params,
                    const unsigned int relevant_steps,
                    const VECTOR& const_data, 
                    const std::vector<double>& time_func, 
                    const double scale,
                    std::vector<VECTOR>& result) const
    {
      // functor for storage
      StorageAdapterConvolutionRHS<VECTOR> 
        store(const_data, time_func, scale, result);
      // functor for wifft
      WeightedIFFTBlocks<kernel_block_type> transform(laplace_params.size() );
      assemble(x_ldofs,y_ldofs,laplace_params,relevant_steps,store,transform);
    }


    //! assemble array RHS with constant time function
    template<typename VECTOR>
    void operator()(const x_ldof_array_type& x_ldofs, 
                    const y_ldof_array_type& y_ldofs,
                    const std::vector<std::complex<double> >& laplace_params,
                    const unsigned int relevant_steps,
                    const VECTOR& const_data, 
                    const double scale,
                    std::vector<VECTOR>& result) const
    {
      // functor for storage
      StorageAdapterArrayRHS<VECTOR> store(const_data, scale, result);
      // functor for wifft
      WeightedIFFTBlocks<kernel_block_type> transform(laplace_params.size() );
      assemble(x_ldofs,y_ldofs,laplace_params,relevant_steps,store,transform);
    }



  private:
  
    template<typename STORAGE>
    void assemble(const x_ldof_array_type& x_ldofs,
                  const y_ldof_array_type& y_ldofs,
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
      for (long it_p=0; it_p<x_se.size(); ++it_p) {typename x_se_array_type::const_iterator iter_x=x_se.begin()+it_p;
#else
      // start loop over row support
      typename x_se_array_type::const_iterator iter_x;
      // possibly parallel x-loop
#ifdef ASSEMBLE_GALERKIN_OMP
#pragma omp parallel for schedule(guided)
#endif
      for (iter_x=x_se.begin();iter_x<x_se.end(); ++iter_x) {
#endif

        // instantiate block to be reused by the integrator
        kernel_block_type block;
 
        // start loop over column support
        typename y_se_array_type::const_iterator iter_y = y_se.begin();
        for (; iter_y<y_se.end(); ++iter_y) {
        
          // set all block entries to zero
          for(unsigned int i=0;i<block.getNumEntries(); ++i)
            block[i].zeros();
     
          // integrate
          integrate_(*iter_x, *iter_y, block);
          // store
          storeBlock(*iter_x, *iter_y, x_ldof_flag, y_ldof_flag, 
                     block, store, x_tag_, y_tag_);
        } // end y_se_it
      } // end x_se_it
    }


    // symmetric
    template<typename STORAGE>
    void assembleSymmetric(const x_ldof_array_type& ldofs,
                           STORAGE& store) const
    {
      // abort in case the boundary does not exist
      if (ldofs.size()==0)
        return;
    
      // get flags
      const LDOF_TYPE ldof_flag = ldofs[0]->getType();
 
      // extract needed superelement arrays
      x_se_array_type se;
      x_se_extractor_(ldofs, se);
    
      // possibly parallel x-loop
      // start loop over row support
#ifdef _MSC_VER && ASSEMBLE_GALERKIN_OMP
#pragma omp parallel for schedule(runtime)
	  for (long it_p=0; it_p<se.size(); ++it_p) {typename x_se_array_type::const_iterator iter_x=se.begin()+it_p;
#else
      typename x_se_array_type::const_iterator iter_x;
#ifdef ASSEMBLE_GALERKIN_OMP
#pragma omp parallel for schedule(runtime)
#endif
      for (iter_x=se.begin(); iter_x<se.end(); ++iter_x) {
#endif
      
        // instantiate block to be reused by the integrator
        kernel_block_type block;
  
        typename x_se_array_type::const_iterator iter_y=iter_x;
        // start loop over column support
        for (; iter_y<se.end(); ++iter_y) {
        
          // set all block entries to zero
          for(unsigned int i=0;i<block.getNumEntries(); ++i) 
            block[i].zeros();
          // integrate
          integrate_(*iter_x, *iter_y, block);
          // store
          storeBlockSymmetric(*iter_x, *iter_y, ldof_flag, ldof_flag, 
                              block, store, x_tag_, y_tag_);
        } // end y_se_it
      } // end x_se_it
    }





    template<typename STORAGE,
             typename TRANSFORM_BLOCK>
    void assemble(const x_ldof_array_type& x_ldofs, 
                  const y_ldof_array_type& y_ldofs,
                  const std::vector<std::complex<double> >& laplace_params,
                  const unsigned int relevant_steps,
                  STORAGE& store, 
                  const TRANSFORM_BLOCK& transform) const
    {
      // abort in case the boundary does not exist
      if (x_ldofs.size()==0 or y_ldofs.size()==0)
        return;
    
      // get number of laplace parameters
      const unsigned int nsteps = laplace_params.size();
    
      // use half of laplace_params for integration
      // 2nd half is conjugate complex
      std::vector<std::complex<double> > half_laplace_params(nsteps/2+1);
      std::copy(laplace_params.begin(),
                laplace_params.begin()+nsteps/2+1,
                half_laplace_params.begin() );  
      
      // get flags
      const LDOF_TYPE x_ldof_flag = x_ldofs[0]->getType();
      const LDOF_TYPE y_ldof_flag = y_ldofs[0]->getType();

      // extract needed superelement arrays
      x_se_array_type x_se;
      y_se_array_type y_se;
      x_se_extractor_(x_ldofs, x_se);
      y_se_extractor_(y_ldofs, y_se);
   
      // copy functors for OpenMP
#ifdef ASSEMBLE_GALERKIN_OMP
      boost::ptr_vector<TRANSFORM_BLOCK> transform_omp;
      boost::ptr_vector<STORAGE> store_omp;
      // allocate transform functor for each thread
      for (unsigned int i=0; i<omp_get_max_threads(); ++i) {
        store_omp.push_back(new STORAGE(store) );
        transform_omp.push_back(new TRANSFORM_BLOCK(transform) );
      }
#endif

      // start loop over row support
      typename x_se_array_type::const_iterator iter_x;
 
      // possibly parallel x-loop
#ifdef ASSEMBLE_GALERKIN_OMP
#pragma omp parallel for schedule(guided)
#endif
      for (iter_x=x_se.begin();iter_x<x_se.end(); ++iter_x) {

        // create kernel blocks
        std::vector<kernel_block_type> in_blocks(nsteps);
        std::vector<typename TRANSFORM_BLOCK::out_block_type> 
          out_blocks(nsteps);

        // start loop over column support
        typename y_se_array_type::const_iterator iter_y = y_se.begin();
        for (; iter_y<y_se.end(); ++iter_y) {
        
          // set blocks to zero
          for (unsigned int i=0; i<nsteps; ++i)
            for(unsigned int j=0;j<in_blocks[i].getNumEntries(); ++j)
              in_blocks[i][j].zeros();
        
          // integrate for half of the laplace parameters
          integrate_(*iter_x, *iter_y, half_laplace_params, in_blocks);
 
          // fill rest with conjugate complex values
          for (unsigned int s=nsteps/2+1; s<nsteps; s++)
            for(unsigned int ngdx=0;ngdx<x_num_gdofs_per_element_; ++ngdx)
              for(unsigned int ngdy=0; ngdy<y_num_gdofs_per_element_; ++ngdy)
                for(unsigned int nldx=0; nldx<num_ldofs_per_gdof_; ++nldx)
                  for(unsigned int nldy=0; nldy<num_ldofs_per_gdof_; ++nldy)
                    in_blocks[s](ngdx, ngdy)(nldx, nldy) = 
                      NumberTraits<value_type>::conj( in_blocks[nsteps-s](ngdx, ngdy)(nldx, nldy) );
        
#ifndef ASSEMBLE_GALERKIN_OMP // sequential
          // weighted ifft
          transform(in_blocks, out_blocks);
          for (unsigned int i=0; i<relevant_steps; ++i) {
            // set internal storage position
            store.setCurrentStep(i);
            storeBlock(*iter_x, *iter_y, x_ldof_flag, y_ldof_flag, 
                       out_blocks[i], store, x_tag_, y_tag_);
          }
#else // parallel
          // weighted ifft
          transform_omp[omp_get_thread_num()](in_blocks, out_blocks);
          for (unsigned int i=0; i<relevant_steps; ++i) {
            // set internal storage position
            store_omp[omp_get_thread_num()].setCurrentStep(i);
            storeBlock(*iter_x, *iter_y, x_ldof_flag, y_ldof_flag, 
                       out_blocks[i], store_omp[omp_get_thread_num()], x_tag_, y_tag_);
          }
#endif
        } // end y_se_it
      } // end x_se_it
    }
  




    template<typename STORAGE,
             typename TRANSFORM_BLOCK>
    void assembleSymmetric(const x_ldof_array_type& ldofs,
                           const std::vector<std::complex<double> >& laplace_params,
                           const unsigned int relevant_steps,
                           STORAGE& store, 
                           const TRANSFORM_BLOCK& transform) const
    {
      // abort in case the boundary does not exist
      if (ldofs.size()==0)
        return;
    
      // get number of laplace parameters
      const unsigned int nsteps = laplace_params.size();
    
      // use half of laplace_params for integration
      // 2nd half is conjugate complex
      std::vector<std::complex<double> > half_laplace_params(nsteps/2+1);
      std::copy(laplace_params.begin(),
                laplace_params.begin()+nsteps/2+1,
                half_laplace_params.begin() );  
      
      // get flags
      const LDOF_TYPE ldof_flag = ldofs[0]->getType();
 
      // extract needed superelement arrays
      x_se_array_type se;
      x_se_extractor_(ldofs, se);
 
      // copy functors for OpenMP
#ifdef ASSEMBLE_GALERKIN_OMP
      boost::ptr_vector<TRANSFORM_BLOCK> transform_omp;
      boost::ptr_vector<STORAGE> store_omp;
      // allocate transform functor for each thread
      for (unsigned int i=0; i<omp_get_max_threads(); ++i) {
        store_omp.push_back(new STORAGE(store) );
        transform_omp.push_back(new TRANSFORM_BLOCK(transform) );
      }
#endif

      // start loop over row support
      typename x_se_array_type::const_iterator iter_x;
 
      // possibly parallel x-loop
#ifdef ASSEMBLE_GALERKIN_OMP
#pragma omp parallel for schedule(guided)
#endif
      for (iter_x=se.begin(); iter_x<se.end(); ++iter_x) {
 
        // create kernel blocks
        std::vector<kernel_block_type> in_blocks(nsteps);
        std::vector<typename TRANSFORM_BLOCK::out_block_type> 
          out_blocks(nsteps);

        // start loop over column support
        typename x_se_array_type::const_iterator iter_y=iter_x;
        for (; iter_y<se.end(); ++iter_y) {
        
          // set blocks to zero
          for (unsigned int i=0; i<nsteps; ++i)
            for(unsigned int j=0;j<in_blocks[i].getNumEntries(); ++j)
              in_blocks[i][j].zeros();
        
          // integrate for half of the laplace parameters
          integrate_(*iter_x, *iter_y, half_laplace_params, in_blocks);
 
          // fill rest with conjugate complex values
          for (unsigned int s=nsteps/2+1; s<nsteps; s++)
            for(unsigned int ngdx=0;ngdx<x_num_gdofs_per_element_; ++ngdx)
              for(unsigned int ngdy=0; ngdy<y_num_gdofs_per_element_; ++ngdy)
                for(unsigned int nldx=0; nldx<num_ldofs_per_gdof_; ++nldx)
                  for(unsigned int nldy=0; nldy<num_ldofs_per_gdof_; ++nldy)
                    in_blocks[s](ngdx, ngdy)(nldx, nldy) = 
                      NumberTraits<value_type>::conj( in_blocks[nsteps-s](ngdx, ngdy)(nldx, nldy) );
        
#ifndef ASSEMBLE_GALERKIN_OMP // sequential
          // weighted ifft
          transform(in_blocks, out_blocks);
          for (unsigned int i=0; i<relevant_steps; ++i) {
            // set internal storage position
            store.setCurrentStep(i);
            storeBlockSymmetric(*iter_x, *iter_y, ldof_flag, ldof_flag, 
                                out_blocks[i], store, x_tag_, y_tag_);
          }
#else // parallel
          // weighted ifft
          transform_omp[omp_get_thread_num()](in_blocks, out_blocks);
          for (unsigned int i=0; i<relevant_steps; ++i) {
            // set internal storage position
            store_omp[omp_get_thread_num()].setCurrentStep(i);
            storeBlockSymmetric(*iter_x, *iter_y, ldof_flag, ldof_flag, 
                                out_blocks[i],
                                store_omp[omp_get_thread_num()], 
                                x_tag_, y_tag_);
          }
#endif
        } // end y_se_it
      } // end x_se_it
    }
  




    // store block with respect to ldof_type aka. boundary using idx
    template<typename BLOCK,
             typename STORAGE>
    void storeBlock(const x_se_type *const x_se,
                    const y_se_type *const y_se,
                    const LDOF_TYPE x_ldof_flag,
                    const LDOF_TYPE y_ldof_flag,
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
              store( i, j,block(ngdx, ngdy)(nldx, nldy) );
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
              store( i, j,block(ngdx, ngdy)(nldx, nldy) );
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
              store( i, j,block(ngdx, ngdy)(nldx, nldy) );
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
              store( i, j, block(ngdx, ngdy)(nldx, nldy) );
            }//y_ldof
          }//x_lodf
        }//y_gdof
      }//x_gdof
    }





    // store block with respect to ldof_type aka. boundary using idx
    template<typename BLOCK,
             typename STORAGE>
    void storeBlockSymmetric(const x_se_type *const x_se,
                             const y_se_type *const y_se,
                             const LDOF_TYPE x_ldof_flag,
                             const LDOF_TYPE y_ldof_flag,
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
            if( x_gdof==NULL ) HYENA_ERROR(); if (x_gdof->getLDof(nldx)->getType()!=x_ldof_flag)
              continue;
            const unsigned int i = x_gdof->getLDof(nldx)->getIDX();
            // loop over all FUNDSOL::num_ldofs_per_gdof on tau_y
            for(unsigned int nldy=0; nldy<num_ldofs_per_gdof_; ++nldy){
              // get index j with boundary distinction
              if( y_gdof==NULL ) HYENA_ERROR(); if (y_gdof->getLDof(nldy)->getType()!=y_ldof_flag)
                continue;
              const unsigned int j = y_gdof->getLDof(nldy)->getIDX();
              // store
              store( i, j, block(ngdx, ngdy)(nldx, nldy) );
              if (x_se != y_se)
                store( j, i, block(ngdx, ngdy)(nldx, nldy) );
            }//y_ldof
          }//x_lodf
        }//y_gdof
      }//x_gdof
    }




    template<typename BLOCK,
             typename STORAGE>
    void storeBlockSymmetric(const x_se_type *const x_se,
                             const y_se_type *const y_se,
                             const LDOF_TYPE x_ldof_flag,
                             const LDOF_TYPE y_ldof_flag,
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
              store( i, j, block(ngdx, ngdy)(nldx, nldy) );
              if (x_se != y_se)
                store( j, i, block(ngdx, ngdy)(nldx, nldy) );
            }//y_ldof
          }//x_lodf
        }//y_gdof
      }//x_gdof
    }




    template<typename BLOCK,
             typename STORAGE>
    void storeBlockSymmetric(const x_se_type *const x_se,
                             const y_se_type *const y_se,
                             const LDOF_TYPE x_ldof_flag,
                             const LDOF_TYPE y_ldof_flag,
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
              store( i, j, block(ngdx, ngdy)(nldx, nldy) );
              if (x_se != y_se)
                store( j, i, block(ngdx, ngdy)(nldx, nldy) );
            }//y_ldof
          }//x_lodf
        }//y_gdof
      }//x_gdof
    }




    // store block irrespectible of ldof_type  using id (whole boundary)
    template<typename BLOCK,
             typename STORAGE>
    void storeBlockSymmetric(const x_se_type *const x_se,
                             const y_se_type *const y_se,
                             const LDOF_TYPE x_ldof_flag,
                             const LDOF_TYPE y_ldof_flag,
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
              store( i, j, block(ngdx, ngdy)(nldx, nldy) );
              if (x_se != y_se)
                store( j, i, block(ngdx, ngdy)(nldx, nldy) );
            }//y_ldof
          }//x_lodf
        }//y_gdof
      }//x_gdof
    }


    
  };

} // end namespace hyena

#endif
 
