%$ NAFN=delta_poisson_dirichlet_wexpgrowth_usingrates
%$ echo $NAFN
%$ echo "constexpr unsigned int SAMPLE_SIZE = 50;" > headerfile.hpp
%$ echo "constexpr double ALPHA = 0.24;" >> headerfile.hpp
%$ echo "constexpr double KAPPA = 2.0 ;" >> headerfile.hpp
%$ echo "constexpr double CEPS = 100. ;" >> headerfile.hpp
%$ echo "constexpr double RHO = 25. ;" >> headerfile.hpp
%$ echo "constexpr double MINF = (2. + (1 + pow(2., 1-KAPPA)/(KAPPA -1)))/2. ;" >> headerfile.hpp
%$ echo "constexpr double BKAPPA = pow(2,KAPPA)*(KAPPA - 2)*(KAPPA - 1);">> headerfile.hpp
%$ echo "constexpr double AKAPPA = ((KAPPA+2) + pow(KAPPA,2))/2.;" >> headerfile.hpp
%$ echo "constexpr double CKAPPA = (KAPPA > 2 ? AKAPPA / BKAPPA : 1.) * 2/pow(MINF,2.) ;">> headerfile.hpp
%$ echo "constexpr int EXPERIMENTS = 1e4;" >> headerfile.hpp
%$ ctangle $NAFN
%$ echo $NAFN
%$ g++ -std=c++26 -m64 -march=native -O3 -x c++ -DNDEBUG $NAFN.c -lm -lgsl -lgslcblas
%$ rm -f gg_*_.txt
%$ ./a.out $(shuf -i 433433-292201013 -n1) > logitresout
%$ seq 49 | awk '{n=50;print log($1/n) - log(1 - ($1/n))}' > nlogits
%$ paste -d',' nlogits logitresout > forplottingfile1
%$ NAFN=delta_poisson_dirichlet_wexpgrowth_usingrates
%$ cweave $NAFN
%$ echo $NAFN
%$ tail -n4 deltanull_incbeta_expgrowth.tex > endi
%$ for i in $(seq 5); do $(sed -i '$d' $NAFN.tex) ; done
%$ cat endi >> $NAFN.tex
%$ emacs --version | head -n1 > innleggemacs
%$ g++ --version | head -n1 > innleggcpp
%$ lualatex --version | head -n2 > innlegglualatex
%$ cweave  --version | head -n1 > innleggcweave
%$ ctangle  --version | head -n1  > innleggctangle
%$ uname  --kernel-release -o  > innleggop
%$ bash --version | head -n1 > innleggbash
%$ sed -i 's/x86/x86\\/g' innleggbash
%$ gsl-config --version > innlegggsl
%$ spix --version > innleggspix
%$ awk --version | head -n1 > innleggawk
%$ lualatex $NAFN
%$ lualatex $NAFN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Copyright (C) 2025 aber der voll idiot
\documentclass[a4paper,10pt,final]{cweb}
\usepackage{fancyvrb,  graphicx, url}
\usepackage[svgnames]{xcolor}
\usepackage[inputlexerlinenos=true]{minted}
\setminted[cpp]{breaklines}
\usepackage{piton}
%% width=min
%% width=12cm
\PitonOptions{width=min,line-numbers,break-lines,indent-broken-lines,background-color=gray!15}
%% a4 paper size 210 x 297 millimeters
%% for xelatex
%\usepackage{xunicode}
\usepackage{fontspec}
\usepackage{xltxtra}
\usepackage{lineno}
\usepackage{amsfonts, amsmath, amssymb}
\usepackage{fullpage}
\usepackage{marvosym}
\usepackage{bm}
\usepackage{bbm}
    \usepackage[]{hyperref}
\hypersetup{
    colorlinks,
    linkcolor={gray!50!black},
    citecolor={blue!50!black},
    urlcolor={blue!80!black}
}
\usepackage{color}
\usepackage{a4wide,fullpage}
\usepackage{setspace}
\usepackage{hyperref}
\usepackage{enumerate}
\usepackage{dsfont}
\usepackage{tabto}
\usepackage{lipsum}
\usepackage{orcidlink}
\usepackage{abstract}
\renewcommand{\abstractname}{}
\usepackage[en-GB,showdow]{datetime2}
\DTMlangsetup[en-GB]{ord=raise,monthyearsep={,\space}}
\usepackage[rightcaption]{sidecap}
\sidecaptionvpos{figure}{m}
\usepackage{siunitx}
\usepackage{tikz}
\usepackage{pgfplots}
\pgfplotsset{compat=newest}
\usepgfplotslibrary{units}
\sisetup{
  round-mode          = places,
  round-precision     = 3,
}
\setstretch{1}
\newcommand{\one}[1]{\ensuremath{\mathds{1}_{\left\{ #1 \right\}}}}
\newcommand{\EE}[1]{\ensuremath{\mathds{E}\left[ #1 \right]}}%
\newcommand{\set}[1]{\ensuremath{\left\{ #1 \right\}}}
\newcommand{\svigi}[1]{\ensuremath{\left( #1 \right)}}
\newcommand{\prb}[1]{\ensuremath{\mathds{P}\left( #1 \right) } }%
\newcommand{\h}[1]{\ensuremath{\uptheta_{ #1 } } }%
\newcommand{\VV}[1]{\ensuremath{ \mathbb{V}\left( #1 \right)}}%
\newcommand{\hp}{\ensuremath{\theta_1}}
\newcommand{\hs}{\ensuremath{\theta_2}}
\newcommand{\D}{\ensuremath{\mathbb{D}}}
\newcommand{\IN}{\ensuremath{\mathds{N}}}
\newcommand{\F}{\ensuremath{\mathbb{F}} }
\newcommand{\bt}[1]{\textcolor{blue}{\tt #1}}
\newcommand{\aths}[1]{\textcolor{violet}{\textbf{\small #1}}}
\newcommand{\s}{{\textcolor{blue}{\ensuremath\bullet}}}
\newcommand{\sa}[1]{\mathrel{\reflectbox{\rotatebox[origin=c]{#1}{\scalebox{0.5}{$\uparrow$}}}}}
\newcommand{\tf}{\ensuremath{\tau_{\text{\sf fix}}}}
\newcommand{\tl}{\ensuremath{\tau_{\text{\sf loss}}}}
\newcommand{\x}{{\circ}}
\newtheorem{thm}{Theorem}[section]
\newtheorem{propn}[thm]{Proposition}%
\newtheorem{lemma}[thm]{Lemma}%
\newtheorem{eg}[thm]{Example}
\newtheorem{defn}[thm]{Definition}
\newtheorem{remark}[thm]{Remark}
\newtheorem{notn}[thm]{Notation}
\newtheorem{corollary}[thm]{Corollary}
\newtheorem{conjecture}[thm]{Conjecture}
\newtheorem{assumption}[thm]{Assumption}
\makeatletter
\renewcommand{\maketitle}{\bgroup\setlength{\parindent}{0pt}
\begin{flushleft}
  {\textsf{\textbf{\@@title}}}
\end{flushleft}

\begin{center}
\textsc{\@@author}
\end{center}
\egroup
}
\makeatother
\title{Gene genealogies in a  haploid population evolving
according to  sweepstakes reproduction  \\ --  approximating
$\EE{R_{i}(n)}$ for the time-changed $\delta_{0}$-Poisson-Dirichlet$(\alpha,0)$-coalescent}
 \author{Bjarki Eldon\footnote{\href{beldon11@@gmail.com}{beldon11@@gmail.com}}\footnote{ compiled @@
{\DTMcurrenttime} on  {\today}; using   \\ kernel  \input{innleggop}
\\
  \input{innleggbash} \\
\input{innleggctangle} \\
 \input{innleggcpp} \\  
GSL \input{innlegggsl} \\ \input{innleggcweave} \\  
\input{innlegglualatex}  \\ \input{innleggspix} \\ \input{innleggawk} \\  \input{innleggemacs}}
\orcidlink{0000-0001-9354-2391}}
\begin{document}
\maketitle
\renewcommand{\abstractname}{\vspace{-\baselineskip}}


\begin{abstract}
Fix $0 < \alpha < 2$ and $0 <\gamma \le 1$.  Let
$\{\xi^{n}\} \equiv \{ \xi^{n}(t) : t \ge 0 \}$ be the continuous-time
time-changed $\delta_{0}$-Poisson-Dirichlet$(\alpha,0)$-coalescent
with time-change function $G(t) = \int_{0}^{t}e^{\rho s}ds$.  Write
$\# A$ for the number of elements in a given finite set $A$,
$L_{i}(n) \equiv \int_{0}^{\tau(n) } \# \left\{ \xi \in \xi^{n}(t) :
\#\xi = i \right\}dt $ and
$L(n) \equiv \int_{0}^{\tau(n)} \# \xi^{n}(t)dt $ and
$\tau(n) \equiv \inf \left\{ t \ge 0 : \# \xi^{n}(t) = 1 \right\} $
for $i \in \{1, 2, \ldots, n-1\}$.  We then have
$L(n) = L_{1}(n) + \cdots + L_{n-1}(n)$. Define
$R_{i}(n) \equiv L_{i}(n)/\sum_{j}L_{j}(n) $ for $i=1,2,\ldots, n-1$.
Interpreting $\{\xi^{n}\}$ as `trees' we may view $L_{i}(n)$ as the
random total length of branches supporting
$i \in \{1, 2, \ldots, n-1\}$ leaves, with the length measured in
coalescent time units, and $n$ sample size.  With this C++ code we
sample the time-changed $\delta_{0}$-Poisson-Dirichlet$(\alpha,0)$-coalescent with
$0 < \alpha < 2$ and exponential growth parameter $\rho \ge 0$ fixed,
and approximate $\EE{R_{i}(n)}$
\end{abstract}

\tableofcontents



@* {\bf Copyright}.
\label{sec:copy}

Copyright {\textcopyright} {\the\year} Bjarki Eldon \newline 

This document and any source code it contains is distributed under the
terms of the GNU General Public Licence (version $\ge 3$).  You should
have received a copy of the licence along with this file (see file
COPYING).


The source codes described in this document are free software: you can
redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation, either
version 3 of the License, or (at your option) any later version.

This document and the code it contains is distributed in the hope that
it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
the GNU General Public License for more details.  You should have
received a copy of the GNU General Public License along with this file
(see COPYING).  If not, see \url{http://www.gnu.org/licenses/}.



@* {\bf compilation and output}. 
\label{sec:compi}


This CWEB \cite{knuth1994cweb} document (the {\tt .w} file) can be
compiled with {\tt cweave} to generate a {\tt .tex} file, and with
{\tt ctangle} to generate a {\tt .c} \cite{kernighan1988c} C++ code
file. 

Use the shell tool {\tt spix} on the  script appearing before  the
preamble (the lines starting with  \%\$); simply 

{\tt spix /path/to/the/sourcefile}

where {\tt sourcefile} is the {\tt .w} file

One   may also  copy the script into a file and run  {\tt parallel}
\cite{tange11:_gnu_paral} :

{\tt parallel --gnu -j1 :::: /path/to/scriptfile}


@* {\bf intro}.
\label{sec:intro}


We consider the continuous-time time-changed
$\delta_{0}$-Poisson-Dirichlet$(\alpha,0)$-coalescent,   with $A \equiv
1/\left( \prod_{j=2}^{n} \left( \sum_{i}\one{k_{i} = j} \right)! \right) $, 
transition rates
\begin{equation}
\label{eq:1}
\begin{split}
\lambda_{n;k_{1},\ldots, k_{r};s} &  = A\binom{n }{k_{1}\ldots k_{r}\, s} \left( \one{r=1,k_{1}=2} \frac{C_{\kappa}}{C_{\kappa} + c(1-\alpha) }    +   \frac{c p_{n;k_{1},\ldots,k_{r};s} }{C_{\kappa} + c(1-\alpha)}  \right)
\end{split}
\end{equation}
where $0 < \alpha < 1$, $\kappa \ge 2$, $c \ge 0$ all fixed, and 
\begin{equation}
\label{eq:2}
\begin{split}
 p_{n;k_{1},\ldots,k_{r};s}  & =  \frac{\alpha^{r + s - 1} \Gamma(r+s)  }{\Gamma(n) } \prod_{i=1}^{r}(k_{i} - \alpha - 1)_{k_{i} - 1} \\
 C_{\kappa} & =   \one{\kappa = 2} \frac{2}{m_{\infty}^{2}} + \one{\kappa > 2} \frac{c_{\kappa}}{2^{\kappa}(\kappa - 2)(\kappa - 1) }
\end{split}
\end{equation}
and $\kappa + 2 < c_{\kappa} < \kappa^{2} $ for $\kappa > 2$. 

@* {\bf code}.
\label{sec:code}



@  {\bf includes}.

the included libraries

@<includes@>=@#
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
#include "headerfile.hpp"



@  {\bf the random number generators}.

the random number generators

@<rngs@>=@#
  std::random_device randomseed;
  /* \newline \aths{ Standard mersenne twister  random number engine} */
  std::mt19937_64 rng(randomseed());

gsl_rng * rngtype ;
static void setup_rng( unsigned long int s )
{
        const gsl_rng_type *T ; 
        gsl_rng_env_setup(); 
        T = gsl_rng_default ;
        rngtype = gsl_rng_alloc(T);
        gsl_rng_set( rngtype,  s) ;
}



@*1 {\bf the coalescent rate}.
\label{sec:lambdanks}


the coalescent rate  $\lambda_{n;k_{1},\ldots,k_{r};s}$


@<$\lambda_{n;k_{1},\ldots,k_{r};s}$@>=@#
static double lambdanks( const double n, const std::vector<unsigned int>& v_k )
{
@#
  assert(v_k[0] > 1);

@#
  
  double d {};
  double k {} ;
  double f {1} ;
  const double r = static_cast<double>(v_k.size());
  
  std::unordered_map<unsigned int, unsigned int> counts {} ;

@#
  const auto descending_factorial = [](const double x, const double m)
  {
  double p = 1;
@#
  for( double i = 0; i < m; ++i){
    p *= (x - i); }
@#
  return p ; };

@#
  const auto veldi = [](const double x, const double y)
  {
    feclearexcept(FE_ALL_EXCEPT);
    const double d = pow(x,y);
   @#
    return( fetestexcept(FE_UNDERFLOW) ? 0. : (fetestexcept(FE_OVERFLOW) ? FLT_MAX : d) ) ;
  };


@#

  for( std::size_t i = 0; i < v_k.size(); ++i){
    f *= descending_factorial( static_cast<double>(v_k[i]) - 1. - ALPHA,  static_cast<double>(v_k[i]) -1);
    /* \newline \aths{ count occurrence of each  merger size} */
    ++counts[v_k[i]];
    k += static_cast<double>(v_k[i]) ;
    d += lgamma(static_cast<double>(v_k[i] + 1)) ; }

@#
  assert( k < n + 1 );

@#

  const double s = n - k;

@#

  const double p = static_cast<double>(std::accumulate( counts.begin(), counts.end(), 0, [](double a, const auto &x){ return a +  lgamma( (double)x.second + 1);}));

@#

  const double l = ((v_k.size() < 2 ? (v_k[0] < 3 ? 1. : 0) : 0)*CKAPPA) + (CEPS * veldi(ALPHA, r+s-1) * tgamma(r+s)   *f/tgamma(n) ) ;

@#

  return ( veldi(exp(1), (lgamma( n + 1.) - d)  - lgamma( n - k + 1 ) - p) * l / (CKAPPA + (CEPS*(1-ALPHA))) ) ;
}


@*1 {\bf partitions summing to a given number}. 
\label{sec:GenPartitions}


generate the  partitions of a given size  summing to a given number |myInt|

@<GenPartitions@>=@#
static double  GenPartitions( const unsigned int m,   const unsigned int myInt,
                   const unsigned int PartitionSize,
                   unsigned int MinVal,
			  unsigned int MaxVal,
			  std::vector< std::pair< double, std::vector<unsigned int> > > & v_l_k,
			  std::vector<double>& lrates_sorting)
{

  /*\newline \aths{ |m| is the given number of blocks; |PartitionSize|
  is the size of the partitions; 
     the partitions sum to  |myInt|} */ 
  

  double lrate {} ;
  double sumrates {} ;
  
    std::vector<unsigned int> partition(PartitionSize);
    unsigned int idx_Last = PartitionSize - 1;
    unsigned int idx_Dec = idx_Last;   
    unsigned int idx_Spill = 0;        
    unsigned int idx_SpillPrev;        

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
        

        const char a = (idx_Spill) ? ~((-3 >> (LeftRemain - MinVal)) << 2) : 11;  
        const char b = (-3 >> (val_Dec - LeftRemain));

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
            @#
        assert( static_cast<unsigned int>(std::accumulate( partition.begin(), partition.end(), 0)) == myInt);
           /* \newline \aths{\S~\ref{sec:lambdanks}} */
        lrate =  lambdanks( static_cast<double>(m), partition) ;

	assert( lrate >= 0) ;
	
	v_l_k.push_back( std::make_pair( lrate, partition) ) ;

	lrates_sorting.push_back( lrate );

	sumrates += lrate;

    } while (idx_Dec <= idx_Last);

    assert( sumrates >= 0) ;
    return sumrates ;
}


@*1 {\bf all partitions summing to a given number}.
\label{sec:allmergers_sum_m}

get all partitions (merger sizes)  summing to a given number 

@<allmergers\_sum\_m@>=@#
static double allmergers_sum_m( const unsigned int n,  const unsigned int m, std::vector< std::pair< double, std::vector<unsigned int> > >& v__l_k,  std::vector<double>&  v_lrates_sort )
{

  /* \newline \aths{ |n| is the number of blocks;
     the partitions sum to |m|} */
  
  const std::vector<unsigned int> v__m {m};
  /* \newline \aths{\S~\ref{sec:lambdanks}} */
  double sumr = lambdanks( static_cast<double>(n), v__m ) ;

  v__l_k.push_back( std::make_pair( sumr, v__m ) ) ;

  v_lrates_sort.push_back( sumr );

  if( m > 3){
    for( unsigned int s = 2; s <= m/2; ++s){
      
      assert(m > 2*(s-1) );
      /* \newline \aths{\S~\ref{sec:GenPartitions}} */
      sumr += GenPartitions(n, m, s, 2, m - (2*(s-1)), v__l_k, v_lrates_sort ); }
  }
  @#
  assert( sumr >= 0) ;
  return sumr ;
}


@  {\bf write mergers to file}.


@<ratesmergersfile@>=@#
static void ratesmergersfile( const unsigned int n, const std::vector<unsigned int>& v__indx, const std::vector< std::pair< double, std::vector<unsigned int> > > & vlk, const double s,   std::vector< std::vector< double > > & a__cmf )
{
@#
  assert( s > 0);
  double cmf {} ;
  std::ofstream f ;
  f.open("gg_" + std::to_string(n) + "_.txt", std::ios::app);
@#
  a__cmf[n].clear() ;
  @#
  for( const auto &i : v__indx)
    {
    @#
      cmf += (vlk[ i ].first) / s ;
      @#
      assert( cmf >= 0) ;
      @#
      a__cmf[n].push_back(cmf) ;
      @#
      assert((vlk[ i ].second).size() > 0 );
      @#
      for( const auto &x : vlk[ i ].second){
      @#
	f << x << ' ' ;}
      f << '\n' ; }
      @#
  f.close();
@#
  assert( abs(cmf - 1.) < 0.999999 );
}


@*1 {\bf get all partitions}.
\label{sec:allmergers_when_n_blocks}


@<allmergers\_when\_n\_blocks@>=@#
static void allmergers_when_n_blocks( const unsigned int n, std::vector<double> & v__lambdan, std::vector< std::vector< double > > & a__cmf )
{
  std::vector< std::pair< double, std::vector<unsigned int> > > vlk {} ;
  std::vector< double > ratetosort {} ;
  ratetosort.clear() ;

  double lambdan {} ;
  vlk.clear() ;
  assert( n > 1);
  @#
  for( unsigned int k = 2 ; k  <= n ; ++k){
    /*\newline \aths{ the partition sums to |k|; the number of blocks is |n|} */
    lambdan += allmergers_sum_m(n, k, vlk, ratetosort); }
  /* \newline \aths{ record the total rate when |n| blocks} */
  /* \newline \aths{ use for sampling time} */
  assert( lambdan > 0);
  @#
  v__lambdan[n] = lambdan ;

  @#
  std::vector<unsigned int> indx (ratetosort.size());
  std::iota( indx.begin(), indx.end(), 0);

@#
  std::stable_sort( indx.begin(), indx.end(), [&ratetosort](const unsigned int x, const unsigned int y){return ratetosort[x] > ratetosort[y];});


  /*\newline \aths{ merger rates sorted in descending order;
   print the cmf and rates to file} */

  ratesmergersfile(n, indx, vlk,  v__lambdan[n], a__cmf);
}



@  {\bf get all partitions}. 


 generate all mergers up to sample size

@<all possible partitions@>=@#
static void allmergers( std::vector<double>& vlmn, std::vector< std::vector<double> > & acmf  )
{
@#   
  for( unsigned int tmpn = 2; tmpn <= SAMPLE_SIZE; ++tmpn )
    {
  @#
  /* \newline \aths{\S~\ref{sec:allmergers_when_n_blocks}} */
      allmergers_when_n_blocks( tmpn, vlmn, acmf );
    }
}


@*1 {\bf tree modules}.
\label{sec:updatetree}



@<updatetree@>=@#
static void  updatetree( std::vector<unsigned int> & tre,   const std::vector<unsigned int>& mergersizes )
{
@#
  assert( mergersizes.size() > 0) ;

  std::vector< unsigned int> newblocks {};
  newblocks.clear() ;
  std::shuffle( tre.begin(), tre.end(), rng) ;
  std::size_t s = tre.size();
  for(const auto &k : mergersizes)
    {
      assert( k > 1);
      assert( k <= s) ;
      s -= k ;
      /* \newline \aths{ record the size of the merging blocks} */
      newblocks.push_back( std::accumulate( tre.rbegin(), tre.rbegin() + k, 0));
      /* \newline \aths{ remove the blocks that merged} */
      tre.resize(s);
    }
    @#
  assert( newblocks.size() > 0) ;
  @#
  assert( static_cast<unsigned int>(std::accumulate( newblocks.begin(), newblocks.end(), 0)) <= SAMPLE_SIZE) ;
@#
  tre.insert( tre.end(), newblocks.begin(), newblocks.end() ) ;
@#
  assert( static_cast<unsigned int>(std::accumulate( tre.begin(), tre.end(), 0)) == SAMPLE_SIZE) ;
}


@  {\bf read partitions}.



@<readmergersizes@>=@#
static void readmergersizes(const unsigned int n,  const unsigned int j, std::vector<unsigned int> & v__mergers )
{
@#
  std::ifstream f("gg_" + std::to_string(n) + "_.txt") ;

@#

  std::string line {} ;

@#

  v__mergers.clear(); 

@#

  for( unsigned int i = 0; std::getline(f, line) && i < j; ++i ){
    if(i >= j-1 ){
  @#    
      std::stringstream ss(line) ;

@#
      v__mergers = std::vector<unsigned int>( std::istream_iterator<unsigned int>(ss), {} );
    } }
  @#
  assert(v__mergers.size() > 0) ;
  @#
  assert(v__mergers[0] > 1);
  @#
  assert(v__mergers.back() > 1);
  @#
  f.close();
}



@*1 {\bf one experiment}.
\label{sec:onexperiment}


@<onexperiment@>=@#
static void onexperiment( std::vector<double>& v__ri,  std::vector<double>& vl, const std::vector<double>& v__lambda__n, const std::vector< std::vector< double > > & a__cmf )
{@#
  std::vector<unsigned int> v__tre(SAMPLE_SIZE, 1) ;
  std::fill( vl.begin(), vl.end(), 0);
  unsigned int lina {} ;
  double tau {};
  double new_time {};

  std::vector<unsigned int> v_merger_sizes (SAMPLE_SIZE/2);
  v_merger_sizes.reserve(SAMPLE_SIZE/2) ;


  auto newtime = [](const double lambdam,  const double tau)
  {
    return  (RHO > DBL_EPSILON ? log1p(-(RHO * exp(-RHO * tau)/lambdam) * log(1. - gsl_rng_uniform_pos(rngtype)))/RHO : gsl_ran_exponential(rngtype, 1./lambdam)) ;
  };

  auto updatelengths = [&v__tre, &vl](const double t)
  {@#
  for( const auto &b: v__tre)
    {@#
      assert( b > 0) ;
      assert( b < SAMPLE_SIZE) ;
      vl[0] += t ;
      vl[b] += t ; }};

@#

  auto samplemerger = [](const unsigned int n, const std::vector< double > & v__cmf )
  {@#
  unsigned int j {} ;
  const double u = gsl_rng_uniform( rngtype);
   @#
  while( u > v__cmf[j]){ ++j;}
@#
  return j ; };

@#

  while( v__tre.size() > 1)
    {@#
      new_time = newtime( v__lambda__n[ v__tre.size()],  tau);
      @#
      tau += new_time;
      @#
      updatelengths( new_time );
      @#
      lina = samplemerger( v__tre.size(),  a__cmf[v__tre.size()] );
      @#
      readmergersizes( v__tre.size(), 1 + lina, v_merger_sizes) ;
      /* \newline \aths{\S~\ref{sec:updatetree}} */
      updatetree( v__tre, v_merger_sizes );
    }
  assert( v__tre.back() == SAMPLE_SIZE) ;
@#
  assert( vl[0] > DBL_EPSILON);
  const double d =  vl[0];
  std::transform(vl.cbegin(), vl.cend(), v__ri.begin(), v__ri.begin(), [&d](const auto &x, const auto &y){return y + (x/d);});
}



@  {\bf go ahead -- approximate $\EE{R_{i}(n)}$}. 


 approximate $\EE{R_{i}(n)}$


@<approximate@>=@#
static void approximate()
{@#
  std::vector<double> vri (SAMPLE_SIZE) ;
  vri.reserve(SAMPLE_SIZE) ;
  std::vector<double> v__l (SAMPLE_SIZE);
  v__l.reserve(SAMPLE_SIZE) ;
  std::vector<double> v_l_n (SAMPLE_SIZE + 1) ;
  v_l_n.reserve(SAMPLE_SIZE + 1) ;
  std::vector< std::vector< double > > a__cmfs (SAMPLE_SIZE + 1, std::vector<double> {} ) ;
@#
  allmergers(v_l_n, a__cmfs );  
  
  int r = EXPERIMENTS + 1 ;

  while( --r > 0)
    {
      /* \newline \aths{\S~\ref{sec:onexperiment}} */
      onexperiment( vri, v__l, v_l_n, a__cmfs);
    }
@#
double x {};
  for( unsigned int  i = 1; i < SAMPLE_SIZE; ++i){
  @#
  x = vri[i]/static_cast<double>(EXPERIMENTS);
  @#
     std::cout << log(x) - log1p( -x) << '\n';}
}






@  {\bf main}.

the |main| module


@C


@<includes@>@#
@<rngs@>@#
/* \newline \aths{\S~\ref{sec:lambdanks}} */
@<$\lambda_{n;k_{1},\ldots,k_{r};s}$@>@#
/* \newline \aths{\S~\ref{sec:GenPartitions}} */
@<GenPartitions@>@#
/* \newline \aths{\S~\ref{sec:allmergers_sum_m}} */
@<allmergers\_sum\_m@>@#
@<ratesmergersfile@>@#
/* \newline \aths{\S~\ref{sec:allmergers_when_n_blocks}} */
@<allmergers\_when\_n\_blocks@>@#
@<all possible partitions@>@#
/* \newline \aths{\S~\ref{sec:updatetree}} */
@<updatetree@>@#
@<readmergersizes@>@#
/* \newline \aths{\S~\ref{sec:onexperiment}} */
@<onexperiment@>@#
@<approximate@>@#

int main(int argc, const char * argv[])
{

@#
setup_rng( static_cast<unsigned long>( atoi(argv[1]))) ;

/* \newline \aths{\S~\ref{sec:onexperiment}} */
 approximate();

@#
gsl_rng_free(rngtype);

@#
return GSL_SUCCESS ;
}



@* {\bf conclusions and bibliography}.
\label{sec:concl}


We approximate $\EE{R_{i}(n)}$ when the coalescent is the time-changed 
$\delta_{0}$-Poisson-Dirichlet$(\alpha,0)$-coalescent. An example in 
 Figure~\ref{fig:graph}



\begin{SCfigure}[0.8][htb]
    \begin{tikzpicture}
      \begin{axis}[
        xlabel = $\log(i/n) - \log(1 - i/n)$,
        axis line style = {draw=none},
        tick style = {draw=none},
        xticklabels={draw=none},
        yticklabels={draw=none},
        legend pos=north east]
        \addplot table[col sep=comma] {forplottingfile1};
        %\addplot table[col sep=comma] {forplottingfile2};
        %\addplot table[col sep=comma] {forplottingfile3};
        %\addlegendentry{$\alpha = 0.05$}
        %\addlegendentry{$\alpha = 1$}
        %\addlegendentry{$\alpha = 1.5$}
       \end{axis}
       \end{tikzpicture}
       \caption{\textcolor{violet}{ \it  An example approximation of $\EE{R_{i}(n)}$
 for the  given  parameter values and graphed as logits against
 $\log(i/n) - \log(1 - i/n)$ for $i = 1,2,\ldots, n-1$ where $n$ is
 sample size}}
    \label{fig:graph}
       \end{SCfigure}



\bibliographystyle{plain}
\begin{thebibliography}{99}

\bibitem[D2024]{D2024}
 Diamantidis,  Dimitrios and Fan,  Wai-Tong (Louis) and Birkner,
 Matthias and Wakeley,  John.  {Bursts of coalescence within
 population pedigrees whenever big families occur}.  Genetics Volume 227,  February
  2024.
 \\
 \url{https://dx.doi.org/10.1093/genetics/iyae030}.


\bibitem[CDEH25]{chetwyn-diggle_beta}
JA~Chetwyn-Diggle and  Bjarki Eldon
\newblock Beta-coalescents when sample size is large.
\newblock 2026
\newblock bioRxiv
\newblock \url{https://doi.org/10.64898/2025.12.30.697022}


\bibitem[DK99]{DK99}
P~Donnelly and T~G Kurtz.
\newblock Particle representations for measure-valued population models.
\newblock {\em Ann Probab}, 27:166--205, 1999.
\\
\url{https://dx.doi.org/10.1214/aop/1022677258}

\bibitem[KL94]{knuth1994cweb}
Donald~Ervin Knuth and Silvio Levy.
\newblock {\em The CWEB system of structured documentation: version 3.0}.
\newblock Addison-Wesley Longman Publishing Co., Inc., Reading, Massachusetts,
  1994.


\bibitem[KR88]{kernighan1988c}
Brian~W Kernighan and Dennis~M Ritchie.
\newblock The {C} programming language, 1988.



\bibitem[Pit99]{P99}
J~Pitman.
\newblock Coalescents with multiple collisions.
\newblock {\em Ann Probab}, 27:1870--1902, 1999.
\\
\url{https://dx.doi.org/10.1214/aop/1022874819}


\bibitem[Sag99]{S99}
S~Sagitov.
\newblock The general coalescent with asynchronous mergers of ancestral lines.
\newblock {\em J Appl Probab}, 36:1116--1125, 1999.
\\
\url{https://doi.org/10.1239/jap/1032374759}


\bibitem[Sch03]{schweinsberg03}
J~Schweinsberg.
\newblock Coalescent processes obtained from supercritical {G}alton-{W}atson
  processes.
\newblock {\em Stoch Proc Appl}, 106:107--139, 2003.
\\
\url{https://doi.org/10.1016/S0304-4149(03)00028-0}


\bibitem[S00]{S00}
J~Schweinsberg.
\newblock Coalescents with simultaneous multiple collisions.
\newblock {\em Electronic Journal of Probability}, 5:1--50, 2000.
\\
\url{https://dx.doi.org/10.1214/EJP.v5-68}



\bibitem[Tan11]{tange11:_gnu_paral}
O~Tange.
\newblock {GNU} parallel -- the command-line power tool.
\newblock The USENIX Magazine, 2011.

\end{thebibliography}


@
\end{document}
