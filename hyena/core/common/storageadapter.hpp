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
 * @file    storageadapter.hpp
 *          Provide store functors for assemble methods 
 *          to unify LHS and RHS assembling routines.
 *          Matrices and vectors are stored so the call structure is identical
 *          \code operator(i,j,entry) \endcode
 * @author  TT, Rf
 * @date    created:     05.08.2011
 *          last change: 10.08.2011
 */
#ifndef storageadapter_hpp
#define storageadapter_hpp

// stl includes 
#include <vector>

// boost
#include <boost/utility.hpp>
#include <boost/ptr_container/ptr_vector.hpp>


namespace hyena
{

  //------------------------------------------------------------------------------
  /**
   * store entry in matrix.
   * matrix assembling.
   */
  template<typename MATRIX>
  class StorageAdapterLHS
    : boost::noncopyable
  {
  private:
    MATRIX& mat_;
	int x_store_offset_, y_store_offset_, x_crop_, y_crop_, x_clip_, y_clip_; // mod. DH, 16.12.2013: offset position in storage matrix

  public:
    // ctor
    StorageAdapterLHS(MATRIX& mat, int x_store_offset=0, int y_store_offset=0,
            int x_crop=0, int y_crop=0,int x_clip=-1, int y_clip=-1) // mod. DH, 16.12.2013: offset position in storage matrix
      : mat_(mat)
    {
        x_store_offset_=x_store_offset; y_store_offset_=y_store_offset;
        x_crop_=x_crop; y_crop_=y_crop;
        x_clip_=(x_clip<0)?mat_.rows():x_clip;
        y_clip_=(y_clip<0)?mat_.cols():y_clip;
    } // mod. DH, 16.12.2013: offset position in storage matrix

    // store
    template <typename T>
    void operator() (const unsigned int i, 
                     const unsigned int j, 
                     const T& entry) 
    {
      // attention: mass matrix assembling collo not working for parallel!
      // just for galerkin parallelization critical
      // superelement vs. collocation point assembling 
		int i_off=i+x_store_offset_, j_off=j+y_store_offset_; // mod. DH, 16.12.2013: offset position in storage matrix
		if(i_off<x_clip_ && j_off<y_clip_ && i_off>=x_crop_ && j_off>=y_crop_) // mod. DH, 16.12.2013: guard against writing beyond matrix size & offset position in storage matrix
#ifdef ASSEMBLE_GALERKIN_OMP
#pragma omp critical
#endif
        {
            mat_(i_off,j_off) += entry; // mod. DH, 16.12.2013: guard against writing beyond matrix size & offset position in storage matrix
        }
    }

  };





  //---------------------------------------------------------------------------
  /**
   * store entry for sparse matrix construction.
   * here just fill up a std::vector with triplets (i,j,val)
   * @fixme parallel, efficiency
   */
  template<typename TRIPLET>
  class StorageAdapterSparseLHS
    : boost::noncopyable
  {
  private:
    std::vector<TRIPLET>& dat_;
    
  public:
    // ctor
    StorageAdapterSparseLHS(std::vector<TRIPLET>& dat)
      : dat_(dat)
    {
      // first reserve entries in triplet list!!
    }
  
    // store
    template <typename T>
    void operator() (const unsigned int i,
                     const unsigned int j, 
                     const T& entry) 
    {
      // just for galerkin parallelization critical
      // superelement vs. collocation point assembling 
#ifdef ASSEMBLE_GALERKIN_OMP
#pragma omp critical
#endif
      dat_.push_back( TRIPLET(i, j, entry) );
    }

  };





  //---------------------------------------------------------------------------
  //#ifndef ASSEMBLE_GALERKIN_OMP // sequential
  /**
   * compute matrix-vector product, and store entry in vector.
   * 'on-the-fly' assembling.
   * sequential or with OpenMP parallelization
   */
  template<typename VECTOR>
  class StorageAdapterRHS
    : boost::noncopyable
  {
  private:
    const VECTOR& data_;
    const double scale_;
    VECTOR& result_;
  
#ifdef ASSEMBLE_GALERKIN_OMP
    boost::ptr_vector<VECTOR> part_vector_;
#endif
  
  public:
  
    // ctor
    StorageAdapterRHS(const VECTOR& data, const double scale, VECTOR& result)
      : data_(data),
        scale_(scale),
        result_(result)
    {
#ifdef ASSEMBLE_GALERKIN_OMP
      // allocate part vectors for each thread
      for (unsigned int i=0; i<omp_get_max_threads(); ++i) {
        part_vector_.push_back(new VECTOR);
        part_vector_.back().resize(result_.size() );
        part_vector_.back().setZero();
      }
#endif
    }
  
    //! matrix (M_{ij}=entry) - vector (data) - product
    //! \f$ res_i = M_{ij} \; v_j \f$
    template <typename T>
    void operator() (const unsigned int i, 
                     const unsigned int j, 
                     const T& entry)
    { 

#ifndef ASSEMBLE_GALERKIN_OMP
      result_(i) += scale_ * entry * data_(j);
#else    
      part_vector_[omp_get_thread_num()](i) += scale_ * entry * data_(j);
#endif
    }
  
    void finalize()
    {
#ifdef ASSEMBLE_GALERKIN_OMP
      for (unsigned int i=0; i<omp_get_max_threads(); ++i)
        result_ += part_vector_[i];
#endif
    }
  
  };





  //---------------------------------------------------------------------------
  /**
   * store entry in array of matrices.
   * matrix assembling.
   * Used in assembling routine for Lubich CQM
   * copyable for parallel
   */
  template<typename MATRIX>
  class StorageAdapterArrayLHS
  {
  private:
    std::vector<MATRIX>& mat_;
    unsigned int curr_step_;
  
  public:
    // ctor
    StorageAdapterArrayLHS(std::vector<MATRIX>& mat)
      : mat_(mat),
        curr_step_(0)
    {}
  
    void setCurrentStep(const unsigned int curr_step)
    {
      HYENA_ASSERT( curr_step < mat_.size() );
      curr_step_ = curr_step;
    }
  
    // store
    template <typename T>
    void operator() (const unsigned int i, 
                     const unsigned int j, 
                     const T& entry) 
    {
#ifdef ASSEMBLE_GALERKIN_OMP
#pragma omp critical
#endif
      mat_[curr_step_](i,j) += entry;
    }
  
  };





  //---------------------------------------------------------------------------
  /**
   * compute convolution with time function applied on const data.
   * matrix-vector product, store entry in array of vectors.
   * 'on-the-fly' assembling. Used in assembling routine for Lubich CQM
   * copyable for parallel
   */
  template<typename VECTOR>
  class StorageAdapterConvolutionRHS
  {
  private:
    const VECTOR& const_data_;
    const std::vector<double>& time_func_;
    const double scale_;
    std::vector<VECTOR>& result_;
    const unsigned int nsteps_;
    unsigned int curr_step_;
  
  public:
    // ctor
    StorageAdapterConvolutionRHS(const VECTOR& const_data, 
                                 const std::vector<double>& time_func, 
                                 const double scale, 
                                 std::vector<VECTOR>& result)
      : const_data_(const_data),
        time_func_(time_func),
        scale_(scale),
        result_(result),
        nsteps_(time_func.size() ),
        curr_step_(0)
    {
      HYENA_ASSERT(time_func_.size() == result_.size() );
    }
  
    void setCurrentStep(const unsigned int curr_step)
    {
      curr_step_ = curr_step;
      HYENA_ASSERT( curr_step_ < nsteps_ );
    }
  
    /**
     * convolution of matrix(M)-vector(v) product \f$ M_{ij} * v_j \f$
     *  \f$ res_{i,k} = \sum_{m=0}^{k-1} M_{ij,k-m} * data_j\; timefun_m  \f$
     */  
    template <typename T>
    void operator() (const unsigned int i, 
                     const unsigned int j, 
                     const T& entry) 
    {
      for (unsigned int step=0; step<nsteps_-curr_step_; ++step) {
        result_[step+curr_step_](i) += scale_ * entry * const_data_(j)
          *time_func_[step];
      }
    }  
  
  };





  //---------------------------------------------------------------------------
  /**
   * compute  matrix-vector product, store entry in array of vectors.
   * 'on-the-fly' assembling. Used in assembling routine for Lubich CQM
   */
  template<typename VECTOR>
  class StorageAdapterArrayRHS
  {
  private:
    const VECTOR& data_;
    const double scale_;
    std::vector<VECTOR>& result_;
    unsigned int curr_step_;
  
  public:
    // ctor
    StorageAdapterArrayRHS(const VECTOR& data, 
                           const double scale, 
                           std::vector<VECTOR>& result)
      : data_(data),
        scale_(scale),
        result_(result)
    { }
  
    void setCurrentStep(const unsigned int curr_step)
    {
      curr_step_ = curr_step;
      HYENA_ASSERT( curr_step_ < result_.size());
    }
  
    /**
     * matrix(M)-vector(v) product \f$ M_{ij} \; v_j \f$
     */  
    template <typename T>
    void operator() (const unsigned int i, 
                     const unsigned int j, 
                     const T& entry) 
    {
      result_[curr_step_](i) += scale_ * entry * data_(j); 
    }  

  };

} // end namespace hyena

#endif // include guard

