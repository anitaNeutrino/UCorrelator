#ifndef _lazylookup_hh
#define _lazylookup_hh

/* lazylookup.hh 
 *
 * Some lazy lookup helpers (maps, vectors, arrays, linear interpolating grid)
 * meant for caching expensive functions "as you go." 
 *
 *
 *  Cosmin Deaconu
 *  <cozzyd@kicp.uchicago.edu> 
 *
 **/ 

#include <map> 
#include <vector> 
#include <array> 
#include <bitset> 
#include <stdio.h> 
#include <functional> 


/* You can use the robin hood unordered_map instead of the STL one
 * by including robin hood before this header */ 
#ifdef ROBIN_HOOD_H_INCLUDED
#define lazylookup_umap robin_hood::unordered_map
#else
#include <unordered_map> 
#define lazylookup_umap std::unordered_map
#endif 


namespace lazylookup
{


  /** A lazy generic map (allowing choice of anything that has a mappish
   * interface for storage). 
   *
   * To avoid two lookups, an insertion with a dummy val is always attempted.
   * You may define what that is in the constructor if it's somehow helpful
   * (i.e. if your value type might be expensive to create in some cases).  If
   * that's a problem, using a pointer type might make sense. 
   *
   **/

  template<class M, typename K, typename V> 
  class generic_map 
  {
    private: 
      V dummyVal;
      M m;
      std::function<V(K)> Fn; 
      std::size_t nfilled; 

    public:
      generic_map(std::function<V(K)> fn, V dummy = V()) 
        : dummyVal(dummy), Fn(fn), nfilled(0) { ; }

      V at(K key) 
      {
        auto res = m.insert(std::pair<K,V>(key,dummyVal)); 
        if (res.second)
        {
          //now actually evaluate! 
          res.first->second = Fn(key); 
          nfilled++; 
        }
        return res.first->second; 
      }
      V operator[](K key) { return at(key); } 

      std::size_t get_n_filled() const { return nfilled; } 
  }; 

  
  /** dict is an alias for using an unordered_map (either
   * robin_hood::unordered_map if robin_hood has been included before this, or
   * std::unordered_map otherwise) */
  template <typename K, typename V> 
  using dict  = generic_map<lazylookup_umap<K,V>,K,V>; 

  /** map is an alias for using an ordered map (std::map). This is often a bad choice.  */
  template <typename K, typename V> 
  using map  = generic_map<std::map<K,V>,K,V>; 



  /** A lazy vector.  This is suitable for when the function takes an unsigned integer
   * and will be densly packed. 
   * 
   * THIS WILL BE AUTOMATICALLY EXTENDED. USE ARR FOR A STATICLY SIZED VERSION. 
   *
   * It is possible to set an upper limit to how large it might be extended by
   * using set_grow_limit(). In that case Fn will be called for all invocations
   * above the max size. 
   *
   */ 

  template<typename T> 
  class vec
  {
    private: 
      std::vector<T> v; 
      std::vector<bool> filled; 
      std::size_t nfilled; 
      std::size_t max_size; 
      std::function<T(std::size_t)> Fn;
    public: 
      vec(std::function<T(std::size_t)> fn) : nfilled(0), max_size(0), Fn(fn) {; }
      vec(std::function<T(std::size_t)> fn, std::size_t N, std::size_t grow_limit = 0) : 
        v(N),filled(N), nfilled(0), max_size(grow_limit), Fn(fn) { ; }
      
      void resize(std::size_t N) 
      {
        v.resize(N); 
        filled.resize(N); 
      }

      // this won't truncate if we're already past max size.
      void set_grow_limit(std::size_t grow_limit) { max_size = grow_limit; } 
      T at(std::size_t i) 
      {

        bool must_resize = i >=v.size(); 
        if (must_resize)
        {
          //too big. we can't do it. 
          if (max_size &&  i+1 > max_size) 
          {
            return Fn(i); 
          }

          v.resize(i+1); 
          filled.resize(i+1); 
        }

        if (must_resize || !filled[i]) 
        {
          filled[i] = true; 
          nfilled++; 
          v[i] = Fn(i); 
        }

        return v[i];
      }

      T operator[](std::size_t i) { return at(i); } 

      std::size_t get_n_filled() { return nfilled; } 
  };

  /** A lazy array
   *  Same as a lazy vector but fixed size. 
   */
  template<std::size_t N, typename T>
  class arr 
  {
    private: 
      std::array<T,N> a; 
      std::bitset<N> filled; 
      std::size_t nfilled; 
      std::function<T(std::size_t)> Fn; 
    public: 

      arr(std::function<T(std::size_t)> fn)
        : nfilled(0), Fn(fn)  {;} 

      T at(std::size_t i)
      {
        // if out of bounds, just evaluate function. 
        if (i >= N)
          return Fn(i); 

        if (!filled[i]) 
        {
          filled[i] = true; 
          a[i] = Fn(i); 
          nfilled++; 
        }
        return a[i]; 
      }

      T operator[](std::size_t i) { return at(i); } 
      std::size_t get_n_filled() { return nfilled; } 
  }; 


  //forward declaration for grid interpolator storage types 
  namespace grid_interpolator_storage
  {
    class dense;
    class sparse;
  }

  /** A lazy grid interpolator. 
   *
   *  This requires a function pointer taking a std::array<double,Ndims>. 
   *  If the point requested is in the bounds, the function will be evaluated
   *  at the appropriate grid points. A linearly-interpolated value will be returned. 
   *  Out of bounds calls will just call the function directly. 
   *
   *  For now just linear interpolation is implemented. 
   *
   *  By default this uses dense grid storage (memory is allocated for the entire grid) for Ndim < 3, sparse otherwise
   * 
   * */ 

  template <std::size_t Ndims, typename Storage =grid_interpolator_storage::dense>
  class grid_interpolator
  {
    private: 
      std::array<std::size_t,Ndims > Ns; 
      std::array<double,Ndims > mins; 
      std::array<double,Ndims > maxes; 

      std::array<double,Ndims > deltas; 
      std::array<std::size_t,Ndims> prods;  //products of dimensions (used in calculating glboal bin number); 
      
      std::size_t Nbins; 
      Storage storage; 
      double eps; 
      std::function<double(const std::array<double, Ndims>&)>  Fn; 


      std::size_t global_bin(const std::array<std::size_t, Ndims> & bins) 
      {
        std::size_t bin = 0; 
        for (std::size_t i = 0; i < Ndims; i++) 
        {
          bin += prods[i] * bins[i]; 
        }
        return bin; 
      }
      double bin_val(const std::array<std::size_t, Ndims> & bins) 
      {
         size_t gbin = global_bin(bins); 
         if (!storage.have(gbin)) 
         {
           std::array<double, Ndims> X; 
           for (std::size_t i = 0; i < Ndims; i++)
           {
             X[i] = mins[i] + bins[i]*deltas[i]; 
           }
           storage.fill(gbin, Fn(X));
         }
         return storage.val(gbin);
      }

      //TODO see about removing some of these branches if possible 
      bool get_bins(const std::array<double, Ndims> & vals, 
                   std::array<std::size_t,Ndims> & low_bins, 
                   std::array<std::size_t,Ndims> & high_bins, 
                   std::array<double,Ndims> & fracs)
      {
        for (std::size_t i = 0; i < Ndims; i++) 
        {
          if (vals[i] > maxes[i]) return false; 
          if (vals[i] < mins[i]) return false; 

          double fbin =  (vals[i] - mins[i]) / deltas[i]; 
          int ibin = fbin;
          low_bins[i] = ibin; 
          high_bins[i] = ibin+1; 
          fracs[i] = fbin-ibin; 

          if (fracs[i] >= 0 && fracs[i] < eps)
          {
            fracs[i] = 0; 
            high_bins[i] = ibin; 
          }
          else if (fracs[i] < 0 && fracs[i] > -eps) 
          {
            low_bins[i]=high_bins[i]; 
            fracs[i] = 0; 
          }
        }

        return true; 
      }



    public: 
      /** Initialize the grid interpolator. 
       *  @param dim_Ns the number of bins in each dimension (number of grid points will be this plus 1)
       *  @param dim_mins the minimum grid point in each dimension
       *  @param dim_maxes the maximum grid point in each timension 
       *  @param bin_eps the fractional part of a bin that will be considered to be exactly the grid point
       *                 (this is useful because it can avoid evaluating the neighboring grid points) 
       */ 
      grid_interpolator( std::function<double(const std::array<double,Ndims>&)>  f, 
                        const std::array<std::size_t, Ndims> & dim_Ns, 
                        const std::array<double, Ndims> & dim_mins, 
                        const std::array<double, Ndims> & dim_maxes,
                        double bin_eps = 1e-4)  
        : Ns(dim_Ns), mins(dim_mins), maxes(dim_maxes) , Nbins(1), eps(bin_eps), Fn(f) 
      {

        for (std::size_t i = 0; i < Ndims; i++) 
        {
          prods[i] = Nbins; 
          Nbins *= (Ns[i]+1); 
          deltas[i] = (maxes[i]-mins[i])/(Ns[i]); 
        }

        storage.make_room(Nbins); 
      }

      double val( const std::array<double, Ndims> & X ) 
      {
         //figure out what bin we are in 
          std::array<std::size_t, Ndims > low_bins; 
          std::array<std::size_t, Ndims > high_bins; 
          std::array<double,Ndims> fracs; 

          //out of bounds, just evaluate the function 
          if (!get_bins(X, low_bins, high_bins, fracs))
          {
            return Fn(X); 
          }


          //specialization for 1D 
          if (Ndims == 1) 
          {

            return !fracs[0] ? bin_val(low_bins) : 
                   fracs[0] ==1 ? bin_val(high_bins) : 
                   fracs[0] * bin_val(high_bins) + (1-fracs[0]) * bin_val(low_bins);
          }


          //specialization for 2D
          if (Ndims ==2) 
          {
            double sum = 0; 
            double f00 = (1-fracs[0]) * (1-fracs[1]); 
            double f11 = fracs[0] * fracs[1]; 
            double f01 = (1-fracs[0]) * fracs[1]; 
            double f10 = (1-fracs[1]) * fracs[0]; 

            if (f00) sum+= f00 * bin_val(low_bins); 
            if (f11) sum+= f11 * bin_val(high_bins); 
            if (f10) sum += f10 * bin_val({high_bins[0], low_bins[1]}); 
            if (f01) sum += f01 * bin_val({low_bins[0], high_bins[1]}); 

            return sum; 
          }

          double sum = 0; 
          std::array<std::size_t,Ndims> bins;

          //loop over all vertices. 
          //Each vertex is either lower (0) or upper (1); 
          //There is probably some template magic to express this at compile time
          // TODO: specialize for low numbers of dimensions 
          for (std::size_t V = 0; V < (1 << Ndims); V++) 
          {
            double prod = 1;
            for (std::size_t i = 0; i < Ndims; i++) 
            {
              bool high = V & (1 << i); 
              prod *= high ? fracs[i]: 1-fracs[i]; 
              if (!prod) break; //this coefficient will be zero. no need to calculate. 
              bins[i] = high ? high_bins[i] : low_bins[i]; 
            }
            if (!prod) continue; 

            sum += prod*bin_val(bins); 
          }
          
          return sum; 
      }

      double operator[](std::array<double,Ndims> where ) { return val(where); } 

  };

  namespace grid_interpolator_storage
  {
    //These are only meant to be used by the grid_interpolator!!!! 
    class dense
    {
        template <std::size_t Ndims,typename S> 
        friend class lazylookup::grid_interpolator; 

        private: 
          void make_room(std::size_t N) 
          {
            vals.resize(N); 
            filled.resize(N);
          }


          bool have(std::size_t bin)
          {
            return filled[bin]; 
          }

          void fill(std::size_t bin, double val) 
          {
            vals[bin] = val; 
            filled[bin] = true; 
          }

          double val(std::size_t bin) 
          {
            return vals[bin]; 
          }


        private: 
          std::vector<double> vals; 
          std::vector<bool> filled; 
    };


    class sparse
    {

      template <std::size_t Ndims,typename S> 
      friend class lazylookup::grid_interpolator; 
      private: 
   
        void make_room(std::size_t N) 
        {
          (void) N; 
        }

        bool have(std::size_t bin) 
        {
          //try to insert a dummy value, so we have the cached version of this afterwards
          cached_insert = m.insert(std::pair<std::size_t,double> (bin,0));
          return !cached_insert.second; 
        }

        //We always check if we have it first, so we can use our cached insert return value
        double val(std::size_t bin) 
        {
          (void) bin; 
          return cached_insert.first->second; 
        }

        void fill(std::size_t bin, double val) 
        {
          (void) bin; 
          cached_insert.first->second = val; 
        }

      private: 
        lazylookup_umap<std::size_t, double> m; 
        std::pair<lazylookup_umap<std::size_t,double>::iterator,bool> cached_insert;
    }; 

  }

  
}
#endif
