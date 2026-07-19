/* **
deltanull-Poisson-Dirichlet(alpha,0) coalescent
   copyright (C) 2026 bjarki eldon
Sun May 31 10:22:12 AM CEST 2026
   use write_mergers_to_files.cpp for writing out merger sizes
lldb --debug --arch x86_64 ./a.out 1234567
 valgrind -v --leak-check=full --leak-resolution=high --num-callers=40 --vgdb=full <call>
**** */
#include <iostream>
#include <cstdlib>
#include <iterator>
#include <limits>
#include <random>
#include <fstream>
#include <iomanip>
#include <vector>
#include <numeric>
#include <functional>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <assert.h>
#include <float.h>
#include <fenv.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>

 std::random_device randomseed;
 std::mt19937_64 rng(randomseed());

gsl_rng * rngtype {nullptr};
static void setup_rng( const  unsigned long int s )
{
        const gsl_rng_type *T ; 
        gsl_rng_env_setup(); 
        T = gsl_rng_default ;
        rngtype = gsl_rng_alloc(T);
        gsl_rng_set( rngtype,  s) ;
}


struct STIKAR
{
  double ALPHA {};
  double KAPPA {};
  double CEPS {};

  STIKAR (
	  const double ialpha__,
	  const double ikappa__,
	  const double icc__) :
    ALPHA (ialpha__), KAPPA (ikappa__), CEPS (icc__) {}
} ;

struct FASTAR
{
  double MINF {} ;
  double BKAPPA {};
  double AKAPPA {} ;
  double CKAPPA {} ;

  FASTAR () {} ;
};


static double descending_factorial(const double x, const unsigned int m)
{

  double p = (m == 1 ? x : 1 );
 switch(static_cast<unsigned int>( (m > 1 ? 1 : 0)))
   {
   case 0 : break ;
   case 1 : {
  for( unsigned int  i = 0; i < m; ++i){
    p *= (x - static_cast<double>(i)); };
  break ; };
   default : break ; };

  return p ;
}


static double lambdanks(const double n,
                        const std::vector<unsigned int>& v_k,
			const STIKAR & s__,
			const FASTAR & f__ )
{

  /* n is current number of blocks */
  /* v_k is the merger sizes for which to compute merger rate */
  auto veldi = []( const double x, const double y )
  {
     feclearexcept(FE_ALL_EXCEPT);
     const double d = pow(x,y);
     return( fetestexcept(FE_UNDERFLOW) ? 0. : (fetestexcept(FE_OVERFLOW) ? FLT_MAX : d) ) ;   
  };
  
  assert(v_k.size()>0);
  assert( std::all_of(std::cbegin(v_k), std::cend(v_k), [](const auto &x){return x > 1;}));
 
  double d {};
  double k {} ;
  double f {1} ;
  const double r = static_cast<double>(v_k.size());
  
  std::unordered_map<unsigned int, unsigned int> counts {} ;
  
  for( std::size_t i = 0; i < v_k.size(); ++i){
    f *= descending_factorial( static_cast<double>(v_k[i]) - 1. - s__.ALPHA,
                               v_k[i] - 1 );
    /* get number of  occurrences  of each  merger size */
    ++counts[v_k[i]];
    k += static_cast<double>(v_k[i]) ;
    d += lgamma(static_cast<double>(v_k[i] + 1)) ; }

  assert( k < n + 1 ); 
  const double s = n - k;
  
  const double p = static_cast<double>(std::accumulate(counts.begin(), counts.end(), 0.,
                                                       [](double a, const auto &x)
						       {return a + lgamma((double)x.second + 1.);}));
  
  const double l = ((v_k.size() < 2 ? (v_k[0] < 3 ? 1. : 0) : 0)*f__.CKAPPA)
                    + (s__.CEPS * veldi(s__.ALPHA, r+s-1) * tgamma(r+s) * f/tgamma(n));

  assert(isnormal(l) != 0);
  
  return ( veldi(exp(1), (lgamma( n + 1.) - d)  - lgamma( n - k + 1 ) - p)
           * l / (f__.CKAPPA + (s__.CEPS*(1-s__.ALPHA))));
}


void lesaskra(const unsigned int n_,
              std::vector< std::vector<unsigned int>>& a__sizes )
{
  a__sizes.clear();
  assert(a__sizes.size() < 1 );
  
  auto summa_sizes = [](const std::vector<unsigned int>& z)
   {
     return static_cast<unsigned int>(std::accumulate( std::cbegin(z), std::cend(z),0));
   } ;

  auto lesa = [&a__sizes,&summa_sizes](const std::string & lina)
   {
     std::stringstream ss(lina) ;
     a__sizes.push_back(std::vector<unsigned int>(std::istream_iterator<unsigned int>(ss), {}));
     assert(a__sizes.back().size() > 0);
     return summa_sizes(a__sizes.back());
     /* 
     return static_cast<unsigned int>(std::accumulate(std::cbegin(a__sizes.back()),
                                                      std::cend(a__sizes.back() ),0));
     */
   };
  /* the merger sizes are stored in file merger_sizes.txt */ 
   std::ifstream f("merger_sizes.txt",  std::ios_base::in) ;
   assert(f.is_open());
   std::string line {} ;
   line.clear();
   for(unsigned int i = 0; (std::getline(f,line) && (lesa(line) < n_ + 1)); ++i){};
   f.close();
   assert(a__sizes.size() > 0);
}


void compute_rates(std::vector<double>& v__totalrate,
                   std::vector<std::vector<double>>& a__rates_for_sorting,
		   const unsigned int m,
		   const std::vector<std::vector<unsigned int>>& a__, 
		   const STIKAR & s_, const FASTAR & f_)
{
  assert(m > 1) ;
  assert(a__rates_for_sorting.size() > m) ;
  a__rates_for_sorting[m].clear();
  auto summa = [](const std::vector<unsigned int>& v)
  {
    return (static_cast<unsigned int>(std::accumulate(std::cbegin(v), std::cend(v),0)));
  };
  double x {} ;
  for( unsigned int i = 0; summa(a__[i]) <= m; ++i)
    {
      assert(a__[i].size() > 0);
      assert(std::all_of(std::cbegin(a__[i]), std::cend(a__[i]), [](const auto &x){return x > 1;}) );
       x = lambdanks(m, a__[i], s_, f_);
       assert(isnormal(x) != 0);
       /* get total rate for m blocks - use for sampling time */
      v__totalrate[m] += x;
      /* record rates for sorting with n as number of blocks */
      a__rates_for_sorting[m].push_back(x);
    };
}


static void updatetree (std::vector<unsigned int> & tre, const std::vector<unsigned int>& mergersizes)
{
  assert(tre.size() > 1);
  assert( mergersizes.size() > 0) ;
  assert(std::all_of(std::cbegin(mergersizes), std::cend(mergersizes), [](const auto &x){return x > 1;}));

  std::vector< unsigned int> newblocks {};
  newblocks.clear() ;
  newblocks.reserve(tre.size());
  std::shuffle( tre.begin(), tre.end(), rng) ;
  std::size_t s = tre.size();
  for(const auto &k : mergersizes)
    {
      assert( k > 1);
      assert( k <= s) ;
      s -= k ;
      /* record the size of the merging blocks */
      newblocks.push_back( static_cast<unsigned int>(std::accumulate( tre.rbegin(), tre.rbegin() + k, 0)));
      assert( newblocks.back() > 1 );
      /* remove the blocks that merged */
      tre.resize(s);
    }
  assert( newblocks.size() > 0) ;
  // assert( std::all_of(std::cbegin(newblocks), std::cend(newblocks), [](const auto &x){return x > 1;}));
  
  tre.insert( tre.end(), std::cbegin(newblocks), std::cend(newblocks));

  newblocks.clear();
  newblocks.resize(0);
  std::vector<unsigned int>().swap(newblocks);
}


static void onexp(const unsigned int n__,
                  const std::vector<double>& lambdan,
		  const std::vector<std::vector<unsigned int>>& a__indx,
		  const std::vector<std::vector<unsigned int>>& a__sizes,
		  const std::vector<std::vector<double>>& a__rates,
		  std::vector<double>& v__ri)
{
  std::vector<unsigned int> tree__ (n__, 1) ;
  std::vector<double> v__l (n__, 0.) ;

  auto get__line = [](const std::vector<double>& v__rates,
                      const std::vector<unsigned int>& indx,
		      const double s)
  {
    unsigned int j {};
    const double u = gsl_rng_uniform_pos(rngtype);
    double f = v__rates[indx[j]] ;
    while(u > f/s){++j ; f += v__rates[indx[j]]; };
    return indx[j] ;
  };
  double t {};
  unsigned int ell {} ;
  while( tree__.size() > 1 )
    {
      t = gsl_ran_exponential(rngtype, 1./lambdan[ tree__.size()]);
      std::for_each(std::cbegin(tree__),
                    std::cend(tree__),
		    [&v__l,&t](const auto &b){v__l[0] +=t;
		                              v__l[b] += t;});
      /* get merger sizes for transition */
      ell = get__line(a__rates[tree__.size()],
                      a__indx[tree__.size()],
		      lambdan[tree__.size()] );
      
      /* update tree */
      updatetree(tree__,  a__sizes[ell] );
    };
    /* record relative lengths */
  assert(v__l[0] > std::numeric_limits<double>::epsilon());
  const double d = v__l[0];
  std::transform(std::cbegin(v__l),
                 std::cend(v__l),
		 std::begin(v__ri),
		 std::begin(v__ri),
		 [&d](const auto &x, const auto &y){return y + (x/d);});
  tree__.clear(); tree__.resize(0);
  std::vector<unsigned int>().swap(tree__);
  v__l.clear();  v__l.resize(0);
  std::vector<double>().swap(v__l);
}


static void approximate(const unsigned int sample_size,
              		const std::vector<double>& v__lambdan,
			const std::vector<std::vector<unsigned int>>& a__indx,
			const std::vector<std::vector<unsigned int>>& a__sizes,
			const std::vector<std::vector<double>>& a__rates)
{
    std::vector<double> v__eri (sample_size, 0.);

    const int EXPERIMENTS = 25 ;
    int r = EXPERIMENTS + 1;
    while(--r > 0)
    {onexp(sample_size,
           v__lambdan,
	   a__indx,
	   a__sizes,
	   a__rates,
	   v__eri);};

   std::for_each(std::cbegin(v__eri),
                 std::cend(v__eri),
		 [](const auto &x){std::cout << x << '\n';});

   v__eri.clear();
   v__eri.resize(0);
   std::vector<double>().swap(v__eri);
}


/* sort the rates in descending order and record the index giving the sorted rates */
static void generate_sorted_index(const std::vector<double>& v__rates,
                                  std::vector<unsigned int>& indx)
{
  assert( v__rates.size() > 0);
  indx.clear();
  indx.resize (v__rates.size(), 0);
  std::iota( indx.begin(), indx.end(), 0);
  std::stable_sort(indx.begin(), indx.end(),
                   [&v__rates](const unsigned int x, const unsigned int y)
		                 {return v__rates[x] > v__rates[y];});
}


static void prufa(const unsigned int sample_size,
           const double ialpha__,
	   const double ikappa__,
	   const double icc__  )
{
  const STIKAR s__ (ialpha__, ikappa__, icc__ ) ;
  
  FASTAR f__ {} ;
  f__.MINF = (2. + (1. + pow(2., 1. - s__.KAPPA))/(s__.KAPPA - 1.))/2. ;
  f__.BKAPPA = pow(2., s__.KAPPA)*(s__.KAPPA - 2.)*(s__.KAPPA - 1.);
  f__.AKAPPA = ((s__.KAPPA+2.) + pow(s__.KAPPA,2.))/2.;
  f__.CKAPPA = ((s__.KAPPA > 2 ? f__.AKAPPA / f__.BKAPPA : 1.) * 2.) / pow(f__.MINF, 2.) ;
  
  std::vector<std::vector<unsigned int> > a__merger__sizes {} ;
  a__merger__sizes.clear();

  std::vector< std::vector<double >> a__rates (sample_size + 1, std::vector<double>()) ;

  std::vector< std::vector<unsigned int>> a__indx (sample_size + 1, std::vector<unsigned int>());

  std::vector<double> v__lambdan (sample_size + 1, 0.) ;
  
  lesaskra(sample_size + 1, a__merger__sizes);

  for(unsigned int n_ = 2; n_ <= sample_size; ++n_ )
    { compute_rates(v__lambdan, a__rates,  n_, a__merger__sizes, s__, f__ );
      generate_sorted_index(a__rates[n_], a__indx[n_] );};
  assert(a__rates.size() > 0);
  assert(a__indx.size() > 0);
  
  approximate(sample_size,
              v__lambdan,
	      a__indx,
	      a__merger__sizes,
	      a__rates);
  
  std::for_each(std::begin(a__merger__sizes),
                std::end(a__merger__sizes),
		[](auto &x){x.clear();
		            x.resize(0);
			    std::vector<unsigned int>().swap(x);});
 
  a__merger__sizes.clear();
  a__merger__sizes.resize(0);
  std::vector<std::vector<unsigned int>>().swap(a__merger__sizes);

  std::for_each(std::begin(a__rates),
                std::end(a__rates),
		[](auto &x){x.clear() ;
		            x.resize(0);
			    std::vector<double>().swap(x);});
  a__rates.clear();
  a__rates.resize(0);
  std::vector<std::vector<double>>().swap(a__rates);
}


int main(int argc, const char * argv[])
{
  /* 
  setup_rng(1234567);
  prufa(static_cast<unsigned int>(atoi(argv[1])), 0.5, 2, 2);
  */
    /*  SAMPLE_SIZE -- ALPHA -- KAPPA -- CEPS */
    setup_rng(static_cast<unsigned long int>(atoi(argv[argc - 1])));
  prufa(static_cast<unsigned int>(atoi(argv[1])),
	atof(argv[2]),
	atof(argv[3]),
	atof(argv[4]));
    /*
  */
  gsl_rng_free(rngtype);
  return GSL_SUCCESS ; 
}
