/*
  copyright (C) 2026 bjarki eldon
  prepare merger size files for sampling from the deltanull_Poisson_Dirichlet(alpha,0)-coalescent
 */
#include <iostream>
#include <cstdlib>
#include <iterator>
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

/*
GenPartitions  taken from 
https://stackoverflow.com/questions/22988674/partition-of-an-integer-number-of-partitions
 */
static void  GenPartitions(const unsigned int m,  const unsigned int myInt,
                   const unsigned int PartitionSize,
                   const unsigned int MinVal,
			 const  unsigned int MaxVal )
{

    std::vector<unsigned int> partition(PartitionSize);
    unsigned int idx_Last = PartitionSize - 1;
    unsigned int idx_Dec = idx_Last; 
    unsigned int idx_Spill = 0;         
    unsigned int idx_SpillPrev; 

    auto writefile = [&partition,&m](void)
    {
      std::ofstream f ;
      f.open("hh_" + std::to_string(m) + "_.txt", std::ios::app);
      for(const auto &p:partition){f << p << ' ';};
      f << '\n';
      f.close();
    };

    
    unsigned int LeftRemain = myInt - MaxVal - (idx_Dec - 1)*MinVal; 
    partition[idx_Dec] = MaxVal + 1; 

    do {
        unsigned int val_Dec = partition[idx_Dec] - 1; 
        partition[idx_Dec] = val_Dec;

        idx_SpillPrev = idx_Spill; 
        idx_Spill = idx_Dec - 1; 

        while (LeftRemain > val_Dec) 
        {
            partition[idx_Spill--] = val_Dec;
            LeftRemain -= val_Dec - MinVal; 
	}

        partition[idx_Spill] = LeftRemain; 

        char a = (idx_Spill) ? ~((-3 >> (LeftRemain - MinVal)) << 2) : 11; 
        char b = (-3 >> (val_Dec - LeftRemain));

        switch (a & b) 
        {
            case 1:
            case 2:
            case 3: idx_Dec = idx_Spill;
                    LeftRemain = 1 + (idx_Spill - idx_Dec + 1)*MinVal; 
                    break;

            case 5: for (++idx_Dec, LeftRemain = (idx_Dec - idx_Spill)*val_Dec; (idx_Dec <= idx_Last) && (partition[idx_Dec] <= MinVal); idx_Dec++)
                        LeftRemain += partition[idx_Dec];

                    LeftRemain += 1 + (idx_Spill - idx_Dec + 1)*MinVal;
                    break;

            case 6:
            case 7:
            case 11:idx_Dec = idx_Spill + 1;
                    LeftRemain += 1 + (idx_Spill - idx_Dec + 1)*MinVal;
                    break;


            case 9: for (++idx_Dec, LeftRemain = idx_Dec * val_Dec; (idx_Dec <= idx_Last) && (partition[idx_Dec] <= (val_Dec + 1)); idx_Dec++)
                        LeftRemain += partition[idx_Dec];

                    LeftRemain += 1 - (idx_Dec - 1)*MinVal;
                    break;

            case 10:for (LeftRemain += idx_Spill * MinVal + (idx_Dec - idx_Spill)*val_Dec + 1, ++idx_Dec; (idx_Dec <= idx_Last) && (partition[idx_Dec] <= (val_Dec - 1)); idx_Dec++)
                        LeftRemain += partition[idx_Dec];

                    LeftRemain -= (idx_Dec - 1)*MinVal;
                    break;
        }

        while (idx_Spill > idx_SpillPrev) 
            partition[--idx_Spill] = MinVal;
	
        assert( static_cast<unsigned int>(std::accumulate( partition.begin(), partition.end(), 0)) == myInt);

	/* record partition in file */
	writefile();

    } while (idx_Dec <= idx_Last);

    partition.clear();
    partition.resize(0);
    std::vector<unsigned int>().swap(partition);
}



static void  allmergers_sum_m( const unsigned int n)
{

  /* n is the number of blocks;
     the partitions sum to m */

  auto writefile = [&n](const unsigned int k)
  {
    std::ofstream f ;
      f.open("hh_" + std::to_string(n) + "_.txt", std::ios::app);
      assert(f.is_open());
      f << k << '\n';
      f.close();
  };
  
  for( unsigned int m = 2; m <= n; ++m){
    writefile(m);
    if( m > 3){
      for( unsigned int s = 2; s <= m/2; ++s){
	
	assert(m > 2*(s-1) );
	GenPartitions(n, m, s, 2, m - (2*(s-1)) ); }
    }
  }
}


int main(int argc, const char * argv[])
{
  /* the file hh_<sample_size>_.txt  contains all mergers from 2 up and
     possible when sample_size  blocks */
  allmergers_sum_m(static_cast<unsigned int>(atoi(argv[argc  -  1])));
  return 0 ;
}
