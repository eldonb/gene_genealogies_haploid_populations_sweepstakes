%% Copyright (C) 2025 Bjarki Eldon
%$ NAFN=delta_poisson_dirichlet_usingrates_wplotting
%$ echo 'constexpr unsigned int SAMPLE_SIZE = 20;' > $NAFN.hpp
%$ echo 'constexpr double ALPHA = 0.05;' >> $NAFN.hpp
%$ echo 'constexpr double KAPPA = 2.0 ;' >> $NAFN.hpp
%$ echo 'constexpr double CEPS = 1. ;' >> $NAFN.hpp
%$ echo 'constexpr double MINF = (2. + (1 + pow(2., 1-KAPPA)/(KAPPA - 1)))/2. ;' >> $NAFN.hpp
%$ echo 'constexpr double BKAPPA = pow(2,KAPPA)*(KAPPA - 2)*(KAPPA - 1);' >> $NAFN.hpp
%$ echo 'constexpr double AKAPPA = ((KAPPA+2) + pow(KAPPA,2))/2.;' >> $NAFN.hpp
%$ echo 'constexpr double CKAPPA = (KAPPA > 2 ? AKAPPA / BKAPPA : 1.) * 2/pow(MINF,2.) ;' >> $NAFN.hpp
%$ echo 'constexpr int EXPERIMENTS = 10000 ;' >> $NAFN.hpp
%$ ctangle $NAFN.w
%$ g++ -std=c++26 -m64 -march=native -O3 -x c++ $NAFN.c -lm -lgsl -lgslcblas
%$ rm -f gg_*_.txt
%$ ./a.out $(shuf -i 43484-2392022 -n1) > logitresout
%$ seq 19 | awk '{S=20;print log($1/S) - log(1 - ($1/S))}' > nlogits
%$ paste -d',' nlogits logitresout > forplottingfile1
%$ sed -i 's/ALPHA = 0.05/ALPHA = 0.95/g' $NAFN.hpp
%$ ctangle $NAFN.w
%$ g++ -std=c++26 -m64 -march=native -O3 -DNDEBUG -x c++ $NAFN.c -lm -lgsl -lgslcblas
%$ rm -f gg_*_.txt
%$ ./a.out $(shuf -i 43484-2392022 -n1) > logitresout
%$ paste -d',' nlogits logitresout > forplottingfile2
%$ emacs --version | head -n1 > innleggemacs
%$ g++ --version | head -n1 > innleggcpp
%$ xelatex --version | head -n1  > innleggxelatex
%$ lualatex --version | head -n2 > innlegglualatex
%$ ctangle  --version | head -n1 > innleggctangle
%$ cweave  --version | head -n1 > innleggcweave
%$ uname  --kernel-release -o  > innleggop
%$ bash --version | head -n1 > innleggbash
%$ spix --version > innleggspix
%$ bash --version | head -n1 | sed -z 's/x86/x86\\/g' > innleggbash
%$ NAFN=delta_poisson_dirichlet_usingrates_wplotting
%$ cweave $NAFN.w
%$ tail -n4 $NAFN.tex > endi
%$ for i in $(seq 5); do $(sed -i '$d' $NAFN.tex) ; done
%$ cat endi >> $NAFN.tex
%$ xelatex $NAFN
%$ xelatex $NAFN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% installing gprof2dot :
%% 
%% python3 -m pip install gprof2dot
%% using gprof :
%% compile with -pg -g option
%% run to generate a gmon.out file 
%% ./a.out (shuf 44-333 -n1)
%% gprof --no-static  ./ |gprof2dot.py | dot -Tpdf -o callgraph.pdf
%% for generating callgraph :
%% compile code with -g option 
%% valgrind --tool=callgrind ./a.out $(shuf -i fjf-kfd -n1)
%% use valgrind  callgrind to generate  callgrindfile.out
%% $ ./gprof2dot.py --format=callgrind --output=out.dot /path/to/callgrindfile.out
%% $ dot -Tpdf out.dot -o graph.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\documentclass[a4paper,10pt,final]{cweb}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{fancyvrb,  graphicx, url}
\usepackage[svgnames]{xcolor}
\usepackage{minted}
%% a4 paper size 210 x 297 millimeters
%% for xelatex
\usepackage{fontspec}
\usepackage{xunicode}
\usepackage{fontspec}
\usepackage{xltxtra}
\usepackage{lineno}
\usepackage[all]{xy}
\usepackage{hyperref}
\hypersetup{
    colorlinks,
    linkcolor={gray!50!black},
    citecolor={blue!50!black},
    urlcolor={blue!80!black}
}
%%\usepackage[bigdelims,vvarbb]{newtxmath}
\usepackage{amsfonts, amsmath, amssymb}
\usepackage{fullpage}
\usepackage{marvosym}
\usepackage{bm}
%%\usepackage[round,numbers,super]{natbib}
\usepackage{natbib}
\usepackage{color}
\usepackage{a4wide,fullpage}
\usepackage{setspace}
\usepackage{hyperref}
%\hypersetup{colorlinks=true}
\hypersetup{
    colorlinks,
    linkcolor={gray!50!black},
    citecolor={blue!50!black},
    urlcolor={blue!80!black}
}
\usepackage{enumerate}
\usepackage{dsfont}
\usepackage{tabto}
\usepackage{lipsum}
\usepackage{orcidlink}
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
\usepackage[en-GB,showdow]{datetime2}
\DTMlangsetup[en-GB]{ord=raise,monthyearsep={,\space}}
\setstretch{1}
\newcommand{\one}[1]{\ensuremath{\mathds{1}_{\left\{ #1 \right\}}}}
\newcommand{\EE}[1]{\ensuremath{\mathds{E}\left[ #1 \right]}}%
\newcommand{\set}[1]{\ensuremath{\left\{ #1 \right\}}}
\newcommand{\prb}[1]{\ensuremath{\mathds{P}\left( #1 \right) } }%
\newcommand{\h}[1]{\ensuremath{\uptheta_{ #1 } } }%
\newcommand{\VV}[1]{\ensuremath{ \mathbb{V}\left( #1 \right)}}%
\newcommand{\D}{\ensuremath{\mathds{D}}}
\newcommand{\F}{\ensuremath{\mathds{F}} }
\newcommand{\IN}{\ensuremath{\mathds{N}} }
\newcommand{\bt}[1]{\textcolor{blue}{\tt #1}}
\newcommand{\aths}[1]{\textcolor{violet}{\textbf{\small #1}}}
\newcommand{\svigi}[1]{\ensuremath{\left( #1 \right)}}
\makeatletter
\renewcommand{\maketitle}{\bgroup\setlength{\parindent}{0pt}
\begin{flushleft}
  {\textbf{\textsf{\@@title}}}
\end{flushleft}

\begin{center}
\textsc{\@@author}
\end{center}
\egroup
}
\makeatother
%% https://www.latextemplates.com/svgnames-colors
%%\pagecolor{PaleGoldenrod}
%%  {\XeLaTeX}  \input{innleggxelatex}
\title{Gene genealogies in haploid populations evolving according to
sweepstakes reproduction \\  -- approximating $\EE{R_i(n)}$ for the
$\delta_{0}$-Poisson-Dirichlet$(\alpha,0)$ coalescent}
\author{Bjarki
Eldon\footnote{\href{beldon11@@gmail.com}{beldon11@@gmail.com} \\
compiled @@ {\DTMcurrenttime} on  {\today}; using   \\ \input{innleggop} \\ \input{innleggbash}  \\  \input{innleggcpp} \\ \input{innleggctangle} \\  \input{innleggcweave}  \\ {\XeLaTeX}  \input{innleggxelatex}  \\ \input{innleggspix} \\ \input{innleggawk}  \\  \input{innleggemacs} }  \orcidlink{0000-0001-9354-2391}}
\begin{document}


\maketitle
\renewcommand{\abstractname}{\vspace{-\baselineskip}}

\begin{abstract} Let  $[n] \equiv \set{1,2,\ldots, n} $ for any $n\in
\IN \equiv \set{1,2,\ldots} $. Write $\set{\xi^{n}} \equiv
\set{\xi^{n}(t); t \ge 0}$  for the
$\delta_{0}$-Poisson-Dirichlet$(\alpha,0)$ coalescent. Write  $\# A$
for the number of elements in a given (finite) set $A$,   $\tau(n)
\equiv  \inf\set{t \ge 0 :  \#\xi^{n}(t) = 1}$,    $L_{i}(n)
\equiv   \int_{0}^{\tau(n)} \# \set{ \xi \in \xi^{n}(t) :  \#\xi = i
}dt $ for all $i\in [n-1]$,   $L(n) \equiv  \int_{0}^{\tau(n)} \#\xi^{n}(t) dt $. Then $L(n)
= L_{1}(n) + \cdots + L_{n-1}(n)$;  write $R_{i}(n) \equiv
L_{i}(n)/L(n)$.  With this C++ simulation code we approximate
$\EE{R_{i}(n)}$.   
\end{abstract}


\tableofcontents


@* {\bf Copyright}.
\label{sec:copy}

Copyright {\textcopyright} {\the\year}  Bjarki Eldon \newline

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


This CWEB \citep{knuth1994cweb} document (the {\tt .w} file) can be
compiled with {\tt cweave} to generate a {\tt .tex} file, and with
{\tt ctangle} to generate a {\tt .c} \citep{kernighan1988c} C++ code
file. 


Use the shell tool {\tt spix} on the  script appearing before  the preamble (the lines starting with  \%\$)  
\begin{Verbatim}[frame=leftline,numbers=left,fontsize=\small,fontshape=it,formatcom=\color{teal}]
%$ NAFN=delta_poisson_dirichlet_usingrates
%$ echo 'constexpr unsigned int SAMPLE_SIZE = 20;' > $NAFN.hpp
%$ echo 'constexpr double ALPHA = 0.05;' >> $NAFN.hpp
%$ echo 'constexpr double KAPPA = 2.0 ;' >> $NAFN.hpp
%$ echo 'constexpr double CEPS = 1. ;' >> $NAFN.hpp
%$ echo 'constexpr double MINF = (2. + (1 + pow(2., 1-KAPPA)/(KAPPA - 1)))/2. ;' >> $NAFN.hpp
%$ echo 'constexpr double CKAPPA = 2/pow(MINF,2.) ;' >> $NAFN.hpp
%$ echo 'constexpr int EXPERIMENTS = 10000 ;' >> $NAFN.hpp
%$ ctangle $NAFN.w
%$ g++ -std=c++26 -m64 -march=native -O3 -x c++ $NAFN.c -lm -lgsl -lgslcblas
%$ rm -f gg_*_.txt
%$ ./a.out $(shuf -i 43484-2392022 -n1) > logitresout
%$ seq 19 P awk '{S=20;print log($1/S) - log(1 - ($1/S))}' > nlogits
%$ paste -d',' nlogits logitresout > forplottingfile1
%$ sed -i 's/ALPHA = 0.05/ALPHA = 0.95/g' $NAFN.hpp
%$ ctangle $NAFN.w
%$ g++ -std=c++26 -m64 -march=native -O3 -DNDEBUG -x c++ $NAFN.c -lm -lgsl -lgslcblas
%$ rm -f gg_*_.txt
%$ ./a.out $(shuf -i 43484-2392022 -n1) > logitresout
%$ paste -d',' nlogits logitresout > forplottingfile2
%$ emacs --version P head -n1 > innleggemacs
%$ g++ --version P head -n1 > innleggcpp
%$ xelatex --version P head -n1  > innleggxelatex
%$ cweave  --version P head -n1 > innleggcweave
%$ uname  --kernel-release -o  > innleggop
%$ bash --version P head -n1 > innleggbash
%$ sed -i 's/x86/x86\\/g' innleggbash
%$ NAFN=delta_poisson_dirichlet_usingrates
%$ cweave $NAFN.w
%$ tail -n4 $NAFN.tex > endi
%$ for i in $(seq 5); do $(sed -i '$d' $NAFN.tex) ; done
%$ cat endi >> $NAFN.tex
%$ xelatex $NAFN.tex
\end{Verbatim}
where {\tt P} is the system pipe operator.   Figure~\ref{fig:graph} records an example
of estimates of $\EE{R_{i}(n)}$.  



One   may also  copy the script into a file and run  {\tt parallel}
\cite{tange11:_gnu_paral} :

{\tt parallel --gnu -j1 :::: /path/to/scriptfile}


@* {\bf intro}.
\label{sec:intro}

Let $n, r, k_{1},\ldots, k_{r},\in \IN$, $n,k_{1},\ldots, k_{r} \ge
2$,  $\sum_{i}k_{i} \le n$, $s = n - \sum_{i} k_{i}$.   
The $\delta_{0}$-Poisson-Dirichlet$(\alpha,0)$ coalescent
$\set{\xi^{n}}$   has
transition rate
\begin{equation}
\label{eq:1}
\begin{split}
\lambda_{n;k_{1},\ldots, k_{r};s} &  = \one{r=1,k_{1}=2} \binom{n}{2} \frac{C_{\kappa}}{C_{\kappa} + c(1-\alpha) }  \\
 & +   \binom{n }{k_{1}\ldots k_{r}\, s} \frac{1}{\prod_{j=2}^{n} \left( \sum_{i}\one{k_{i} = j} \right)!  }   \frac{c}{C_{\kappa} + c(1-\alpha)} p_{n;k_{1},\ldots,k_{r};s}
\end{split}
\end{equation}
where $0 < \alpha < 1$, $\kappa \ge 2$, $c \ge 0$ all fixed, and 
\begin{equation}
\label{eq:2}
\begin{split}
 p_{n;k_{1},\ldots,k_{r};s}  & =  \frac{\alpha^{r + s - 1} \Gamma(r+s)  }{\Gamma(n) } \prod_{i=1}^{r}(k_{i} - \alpha - 1)_{k_{i} - 1} \\
 C_{\kappa} & =   \one{\kappa = 2} \frac{2}{m_{\infty}^{2}} + \one{\kappa > 2} c_{\kappa}
\end{split}
\end{equation}
and $\kappa + 2 < c_{\kappa} < \kappa^{2} $ for $\kappa > 2$.    The
$\delta_{0}$-Poisson-Dirichlet$(\alpha,0)$ coalescent is an example of
a simultaneous multiple-merger coalescent. It allows 
simultaneous mergers in   up to $\lfloor n/2 \rfloor $ groups for any
given $n$.   



In \S~\ref{sec:code} we summarise the algorithm, the code follows in
\S~\ref{sec:includes} -- \S~\ref{sec:main}, we
conclude in \S~\ref{sec:concl}.  Comments within the code are
in \aths{this font and colour}



@* {\bf code}.
\label{sec:code}


we summarise the algorithm; $\lambda_{m}$ denotes the total transition
rate
\begin{equation}
\label{eq:3}
\lambda_{m} = 
 \sum_{\substack{ 2 \le k_{1} \le \ldots \le k_{r} \le  m \\ k_{1} + \cdots + k_{r}   \le m  }}
\lambda_{m;k_{1},\ldots, k_{r}; s}
\end{equation}
\begin{enumerate}
\item generate all possible (simultaneous)  mergers for $m= 2,3,\ldots, n$ blocks
\item for every merger(s) sizes   compute and record the transition
rate \eqref{eq:1}
\item record the total transition rate $\lambda_{m}$ \eqref{eq:3}  for
$m=2,3,\ldots, m$
\item $(r_{1},\ldots, r_{n-1}) \leftarrow (0,\ldots, 0)$
\item for each of $M$ experiments :
\begin{enumerate}
\item $(\ell_{1}, \ldots, \ell_{n-1}) \leftarrow (0,\ldots, 0) $
\item $m \leftarrow n$
\item $(\xi_{1},\ldots, \xi_{n}) \leftarrow (1,\ldots, 1) $
\item {\bf while}   $m > 1$ :
\begin{enumerate}
\item $t \leftarrow $ Exp($\lambda_{m})$ 
\item $\ell_{\xi} \leftarrow  \ell_{\xi} + t $ for $\xi = \xi_{1},
\ldots, \xi_{m}$
\item sample merger sizes $k_{1}, \ldots, k_{r}$ using the transition rates
\item given merger sizes merge blocks
\item $m \leftarrow m - \sum_{i}k_{i} + r $
\end{enumerate}
\item $r_{i} \leftarrow  \ell_{i}/\sum_{j}\ell_{j} $ for  $i =
1,2,\ldots, n-1$
\end{enumerate}
\item return $r_{i}/M$ as an approximation of $\EE{R_{i}(n)}$ for $i=1,2,\ldots,n-1$
\end{enumerate}


@*1 {\bf includes}.
\label{sec:includes}


@<includes@>=@#
#include <iostream>
#include <cstdlib>
#include <iterator>
#include <random>
#include <fstream>
#include <iomanip>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <assert.h>
#include <float.h>
#include <fenv.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include "./delta_poisson_dirichlet_usingrates_wplotting.hpp"




@*1 {\bf the random number generators}.
\label{sec:rngs}

the STL and GSL random number generators


@<rngs@>=@#
/* \newline \aths{ obtain a seed out of thin air for the random number engine} */
  std::random_device randomseed;
  /* \newline \aths{ Standard mersenne twister  random number engine
  seeded with $|rng()|$ } */
  std::mt19937_64 rng(randomseed());



gsl_rng * rngtype ;


static void setup_rng( unsigned long int s )
{

@#
        const gsl_rng_type *T ; 
        gsl_rng_env_setup(); 
        T = gsl_rng_default ;
        rngtype = gsl_rng_alloc(T);
        gsl_rng_set( rngtype,  s) ;
}



@*1 {\bf descending factorial}.
\label{sec:descfact}

the descending factorial  $(x)_{m} \equiv x(x-1)\cdots (x- m+1)$ with
$(x)_{0} \equiv 1$; recall \eqref{eq:2}


@<descending factorial@>=@#
static double descending_factorial(const double x, const double m)
{

@#

  double p = 1;

  for( double i = 0; i < m; ++i){
    p *= (x - i); }

  return p ;
}



@*1 {\bf $x^{y}$}. 
\label{sec:power}


compute $x^{y}$ guarding against over- and underflow 

@<guarded $x^{y}$ module@>=@#
static double veldi( const double x, const double y )
  {

@#

    feclearexcept(FE_ALL_EXCEPT);
    const double d = pow(x,y);

    return( fetestexcept(FE_UNDERFLOW) ? 0. : (fetestexcept(FE_OVERFLOW) ? FLT_MAX : d) ) ;   
  }


@*1 {\bf $\lambda_{n;k_{1},\ldots,k_{r};s}$}.
\label{sec:lambdanks}



compute   $\lambda_{n;k_{1},\ldots,k_{r};s}$ \eqref{eq:1}; see the
header file for value of  |ALPHA|, |KAPPA|, |CEPS|, and   |CKAPPA| as
in \eqref{eq:2}



@<compute  $\lambda_{n;k_{1},\ldots,k_{r};s}$@>=@#
static double lambdanks( const double n, const std::vector<unsigned int>& v_k )
{

@#

 /* \newline \aths{$n$ is the given number of blocks;  |v_k| records
 the given merger sizes $k_{1},\ldots, k_{r}$ } */


  assert(v_k[0] > 1);

  
  double d {};
  double k {} ;
  double f {1} ;
  const double r = static_cast<double>(v_k.size());
  
  std::unordered_map<unsigned int, unsigned int> counts {} ;
  
  for( std::size_t i = 0; i < v_k.size(); ++i){
    f *= descending_factorial( static_cast<double>(v_k[i]) - 1. - ALPHA,  static_cast<double>(v_k[i]) -1);
    /* \newline  \aths{ count occurrence of each  merger size} */
    ++counts[v_k[i]];
    k += static_cast<double>(v_k[i]) ;
    d += lgamma(static_cast<double>(v_k[i] + 1)) ; }

  assert( k < n + 1);

  const double s = n - k;


@#
  /* \newline \aths{$p = \sum_{j} \log\Gamma(1 +  \sum_{i} \one{k_{i}=j })$ } */
  const double p = static_cast<double>(std::accumulate( counts.begin(), counts.end(), 0, [](double a, const auto &x){ return a +  lgamma(x.second + 1);}));


@#

  /* \newline \aths{|veldi|  \S~\ref{sec:power}} */
  const double l = ((v_k.size() < 2 ? (v_k[0] < 3 ? 1. : 0) :0)*CKAPPA) + (CEPS * veldi(ALPHA, r+s-1)*tgamma(r+s) *f/tgamma(n) ) ;

  /* \newline \aths{$\ell =   \one{r=1,k_{1} = 2}C_{\kappa}  +  c
  \alpha^{r + s - 1} (\Gamma(r+s)/\Gamma(n)) * \prod_{i=1}^{r}(k_{i} - 1 -
  \alpha)_{k_{i} - 1} $}   */
  
  return ( veldi(exp(1), (lgamma( n + 1.) - d)   - lgamma(n-k +1 ) - p) * l / (CKAPPA + (CEPS*(1-ALPHA))) ) ;
}




@*1 {\bf $r$ mergers summing to $m$}.
\label{sec:msumm}

generate all ordered $r$ merger sizes  $2 \le k_{1} \le \ldots \le k_{r}
\le n$ where $\sum_{i} k_{i} = |myInt|$ and $r$ is fixed;  for example for
$|myInt|=6$ and  $r = 2$  we should get $(2,4)$ and $(3,3)$. 



@<$r$-merger summing to $m$@>=@#
static double  GenPartitions(const unsigned int m,  const unsigned int myInt,
                   const unsigned int PartitionSize,
                   unsigned int MinVal,
			  unsigned int MaxVal,
			  std::vector< std::pair< double, std::vector<unsigned int> > > & v_l_k,
			  std::vector<double>& lrates_sorting)
{

@#

  /* \newline \aths{compute all |PartitionSize| ordered  merger sizes summing
  to |myInt| } */


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
        /* Value AFTER decrementing */
        unsigned int val_Dec = partition[idx_Dec] - 1;
          /* Decrement at the Decrement Point */
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

       @#

        /* \newline \aths{compute $\lambda_{m;k_{1}, \ldots, k_{r};s}$
        \eqref{eq:1} for the given merger size(s)   |lambdanks| \S~\ref{sec:lambdanks} } */
        lrate =  lambdanks( static_cast<double>(m), partition) ;

	assert( lrate >= 0) ;
	
	v_l_k.push_back( std::make_pair( lrate, partition) ) ;

	lrates_sorting.push_back( lrate );

	sumrates += lrate ; 
	

    } while (idx_Dec <= idx_Last);

    assert( sumrates >= 0) ;

@#

   /* \newline \aths{ return the sum of the   |PartitionSize|  merger sizes summing
   to |myInt| } */

    return sumrates ;
}


@*1 {\bf merger sizes summing to $m$}.
\label{sec:mergerstom}


generate all  ordered merger sizes  summing to a given  $m \le n$  when $n$ is the  given number of blocks

@<merger sizes suming to $m$@>=@#
static double allmergers_sum_m( const unsigned int n,   const unsigned int m, std::vector< std::pair< double, std::vector<unsigned int> > >& v__l_k,  std::vector<double>&  v_lrates_sort )
{

 @#


  /* \newline \aths{generate all  ordered merger sizes summing to |m|
  } */

  const std::vector<unsigned int> v__m {m};

  /* \newline \aths{\S~\ref{sec:lambdanks}} */
  double sumr = lambdanks( static_cast<double>(n), v__m ) ;

  v__l_k.push_back( std::make_pair( sumr, v__m ) ) ;

  v_lrates_sort.push_back( sumr );

  if( m > 3){
    for( unsigned int s = 2; s <= m/2; ++s){
      
      assert(m > 2*(s-1) );
      /* \newline \aths{\S~\ref{sec:msumm} } */
      sumr += GenPartitions(n, m, s, 2, m - (2*(s-1)), v__l_k, v_lrates_sort ); }
  }
  assert( sumr >= 0) ;
  return sumr ;
}



@*1 {\bf merger sizes to file}.
\label{sec:mergerstofile}


print merger sizes to file

@<print merger sizes to file@>=@#
static void ratesmergersfile( const unsigned int n, const std::vector<unsigned int>& v__indx, const std::vector< std::pair< double, std::vector<unsigned int> > > & vlk, const double s,   std::vector< std::vector< double > > & a__cmf )
{

@#

  assert( s > 0);
  double cmf {} ;
  std::ofstream f ;
  f.open("gg_" + std::to_string(n) + "_.txt", std::ios::app);

  a__cmf[n].clear() ;
  for( const auto &i : v__indx)
    {
    @#
    /* \newline \aths{record the corresponding cmf} */
      cmf += (vlk[ i ].first) / s ;
      assert( cmf >= 0) ;
      /* \newline \aths{ possibly write also the corresponding  cmf to
      file;  f << cmf << ' ' ;} */
      a__cmf[n].push_back(cmf) ;
      assert((vlk[ i ].second).size() > 0 );
      for( const auto &x : vlk[ i ].second){
	f << x << ' ' ;}
      f << '\n' ; }

@#

  f.close();
  
  
  assert( abs(cmf - 1.) < 0.999999 );
}




@*1 {\bf all merger sizes when $m$ blocks}. 
\label{sec:allmergersm}


generate all possible  ordered  merger sizes when  $m$ blocks 

@<possible merger sizes when $m$ blocks@>=@#
static void allmergers_when_n_blocks( const unsigned int n, std::vector<double> & v__lambdan, std::vector< std::vector< double > > & a__cmf )
{

@#

  std::vector< std::pair< double, std::vector<unsigned int> > > vlk {} ;
  std::vector< double > ratetosort {} ;
  ratetosort.clear() ;

  double lambdan {} ;
  vlk.clear() ;
  assert( n > 1);
  for( unsigned int k = 2 ; k  <= n ; ++k){
  /* \newline \aths{\S~\ref{sec:mergerstom}} */
    lambdan += allmergers_sum_m(n, k, vlk, ratetosort); }
  /* \newline \aths{ record the total rate for each given number of
  blocks; 
 use for sampling time between mergers } */
  assert( lambdan > 0);
  v__lambdan[n] = lambdan ;

@#
  
  std::vector<unsigned int> indx (ratetosort.size());


@#

  std::iota( indx.begin(), indx.end(), 0);

@#

  std::stable_sort( indx.begin(), indx.end(), [&ratetosort](const unsigned int x, const unsigned int y){return ratetosort[x] > ratetosort[y];});


  /* \newline \aths{\S~\ref{sec:mergerstofile} } */
  ratesmergersfile(n, indx, vlk,  v__lambdan[n], a__cmf);
}



@*1 {\bf all possible merger sizes given sample size}.
\label{sec:allmergers}


generate all possible  merger sizes given  sample size  |SAMPLE_SIZE| 


@<allmergers@>=@#
static void allmergers( std::vector<double>& vlmn, std::vector< std::vector<double> > & acmf  )
{


@#


  for( unsigned int tmpn = 2; tmpn <= SAMPLE_SIZE; ++tmpn )
    {
      /* \newline \aths{\S~\ref{sec:allmergersm}} */
      allmergers_when_n_blocks( tmpn, vlmn, acmf );
    }
}




@*1 {\bf update tree}.
\label{sec:updatetree}


merge blocks and update the current  tree  configuration, the current 
block sizes. given merger sizes  $k_{1}, \ldots, k_{r}$, $k =
\sum_{i}k_{i}$,   and the
current  tree with block sizes  $(\xi_{1}, \ldots, \xi_{m})$ where
$\xi_{j}\in [n]$ and $\sum_{j}\xi_{j} = n$    with $m$ blocks the goal
here is to
\begin{displaymath}
\begin{split}
\text{randomize the blocks: } & \svigi{\xi_{\sigma(1)},\ldots, \xi_{\sigma(m)} } \\
\text{merge blocks: } & \svigi{\xi_{\sigma(1)}, \ldots,
\xi_{\sigma(m-k)}, \underset{ \xi_{r}^{\prime} \text{ is sum of $k_{r}$ blocks} }{\underbrace{\xi, \ldots, \xi}}, \ldots,  \underset{\xi_{1}^{\prime}   \text{ is sum of $k_{1}$ blocks}  }{\underbrace{\xi, \ldots, \xi_{\sigma(m) }}}   }  \\
\text{return the new tree: } & \svigi{\xi_{\sigma(1)}, \ldots, \xi_{\sigma(m-k)}, \xi_{r}^{\prime}, \ldots, \xi_{1}^{\prime}  }
\end{split}
\end{displaymath}


@<update the tree configuration@>=@#
static void  updatetree( std::vector<unsigned int> & tre,   const std::vector<unsigned int>& mergersizes )
{

 @#

  assert( mergersizes.size() > 0) ;

@#

  std::vector< unsigned int> newblocks {};
  newblocks.clear() ;
  std::shuffle( tre.begin(), tre.end(), rng) ;
  std::size_t s = tre.size();

@#

  for(const auto &k : mergersizes)
    {
      assert( k > 1);
      assert( k <= s) ;
      s -= k ;
      /* \newline \aths{ record the size of the merging blocks} */
      newblocks.push_back( std::accumulate( tre.rbegin(), tre.rbegin() + k, 0));
      /* \newline \aths{ remove the blocks that merged}  */
      tre.resize(s);
    }

@#

  assert( newblocks.size() > 0) ;

@#

  tre.insert( tre.end(), newblocks.begin(), newblocks.end() ) ;

@#


  assert( static_cast<unsigned int>(std::accumulate( tre.begin(), tre.end(), 0)) == SAMPLE_SIZE) ;
}



@*1 {\bf update lengths}.
\label{sec:updatelengths}


given current block sizes $\xi_{1},\ldots, \xi_{m}$ of the current $m$
blocks the goal here is
to
\begin{enumerate}
\item sample a random exponential  time  with rate $\lambda_{m}$
\item $\ell_{\xi} \leftarrow \ell_{\xi} + t$ for  $\xi =
\xi_{1},\ldots, \xi_{m}$
\end{enumerate}
where $\ell_{\xi}$ a realisation of $L_{\xi}(n)$ 

@<update $\ell_{i}(n)$@>=@#
static void updatelengths( const std::vector<unsigned int>& tre,  std::vector<double>& v__lengths,  const std::vector<double>& v_lambdan)
{

@#

  const double t = gsl_ran_exponential(rngtype,  1./v_lambdan[tre.size()] );


@#


  for( const auto &b: tre)
    {
      assert( b > 0) ;
      assert( b < SAMPLE_SIZE) ;
      v__lengths[0] += t ;
      v__lengths[b] += t ;
    }
}




@*1 {\bf update $r_{i}$}. 
\label{sec:updateri}

update the approximation of $\EE{R_{i}(n)}$; given lengths $\ell_{1},
\ldots, \ell_{n-1}$ the goal here is to
\begin{displaymath}
r_{i} \leftarrow r_{i} +   \frac{\ell_{i}}{\ell_{1} + \cdots + \ell_{n-1}}
\end{displaymath}
for $i = 1,2,\ldots, n-1$


@<update $r_{i}$@>=@#
static void updateri( const std::vector<double>& v__l, std::vector<double>& v_ri )
{

@#

/* \newline \aths{ |v__l[0]| is the sum $\sum_{i}\ell_{i}$ } */
  assert( v__l[0] > 0);

@#


  const double d = v__l[0];


@#

  std::transform(v__l.begin(), v__l.end(), v_ri.begin(), v_ri.begin(), [&d](const auto &x, const auto &y){return y + (x/d);});
}



@*1 {\bf sample merger sizes}.
\label{sec:samplemerger}

sample merger size(s) $(k_{1}, \ldots, k_{r})$ by  returning
$\inf\set{j \in \set{0,1,2} : F_{j} > u }$ where $u$ is a random
uniform and  $F_{j}$ a cumulative mass function generated by ordering
the merger size(s) in descending order by the rate $\lambda_{n;k_{1},
\ldots, k_{r};s}$ \eqref{eq:1} 

@<get merger size(s)@>=@#
static unsigned int samplemerger(const unsigned int n, const std::vector< double > & v__cmf )
{

@#

/* \newline \aths{|n| the current number of blocks; |v__cmf| the
cumulative mass function  }  */


  unsigned int j {} ;
  const double u = gsl_rng_uniform( rngtype);

@#

  while( u > v__cmf[j]){ ++j;}

  return j ;
}



@*1 {\bf read merger size(s) from file}.
\label{sec:readmergersizes}


read the merger sizes(s) corresponding to the index sampled in \S~\ref{sec:samplemerger}


@<read in merger size(s) from file@>=@#
static void readmergersizes(const unsigned int n,  const unsigned int j, std::vector<unsigned int> & v__mergers )
{

@#

/* \newline \aths{|n| the current number of blocks; |j| the index
sampled using  \S~\ref{sec:samplemerger}; |v__mergers| will record the
merger size(s)} */

  std::ifstream f("gg_" + std::to_string(n) + "_.txt") ;

@#

  std::string line {} ;

@#

  v__mergers.clear(); 


@#

  for( unsigned int i = 0; std::getline(f, line) && i < j; ++i ){


@#
    if(i >= j-1 ){

@#

      std::stringstream ss(line) ;

@#

      v__mergers = std::vector<unsigned int>( std::istream_iterator<unsigned int>(ss), {} );
    } }

@#


  assert(v__mergers.size() > 0) ;
  assert(v__mergers[0] > 1);
  assert(v__mergers.back() > 1);
  f.close();
}



@*1 {\bf one realisation of $\svigi{L_{1}(n),\ldots, L_{n-1}(n)}$}. 
\label{sec:onexperiment}


get one realisation of   $\svigi{L_{1}(n),\ldots, L_{n-1}(n)}$; given
$\svigi{\ell_{1}, \ldots, \ell_{n-1}}$ update $\svigi{r_{1},\ldots, r_{n-1}}$


the goal here is
\begin{enumerate}
\item {\bf while} at least two blocks :
\begin{enumerate}
\item sample time until merger and update $\ell_{i} \leftarrow
\ell_{i} + t$  \S~\ref{sec:updatelengths}
\item  sample merger size(s) $k_{1},\ldots, k_{r}$ \S~\ref{sec:samplemerger}
\item merge blocks and update tree   \S~\ref{sec:updatetree}
\end{enumerate}
\item update $r_{i} \leftarrow r_{i} +  \ell_{i}/\sum_{j}\ell_{j}$ \S~\ref{sec:updateri}
\end{enumerate}


@<get $\svigi{\ell_{1}, \ldots, \ell_{n-1}}$@>=@#
static void onexperiment( std::vector<double>& v__ri,  std::vector<double>& vl, const std::vector<double>& v__lambda__n, const std::vector< std::vector< double > > & a__cmf )
{

@#

/* \newline \aths{ |v__ri| holds the approximations $r_{1},\ldots,r_{n-1}$;  |vl| holds $\ell_{1},\ldots, \ell_{n-1}$;  |v__lambda__n|
the total rates $\lambda_{n}$;  |a__cmf| the cumulative mass functions
for sampling merger size(s)} */


/* \newline \aths{initialise the block sizes} */

  std::vector<unsigned int> v__tre(SAMPLE_SIZE, 1) ;
  std::fill( vl.begin(), vl.end(), 0);
  unsigned int lina {} ;

@#

  std::vector<unsigned int> v_merger_sizes (SAMPLE_SIZE/2);
  v_merger_sizes.reserve(SAMPLE_SIZE/2) ;

@#

  while( v__tre.size() > 1)
    {
    /* \newline \aths{\S~\ref{sec:updatelengths} } */
      updatelengths( v__tre, vl, v__lambda__n );
      /* \newline \aths{\S~\ref{sec:samplemerger} } */
      lina = samplemerger( v__tre.size(),  a__cmf[v__tre.size()] );
      /* \newline \aths{\S~\ref{sec:readmergersizes}} */
      readmergersizes( v__tre.size(), 1 + lina, v_merger_sizes) ;
      /* \newline \aths{\S~\ref{sec:updatetree}} */
      updatetree( v__tre, v_merger_sizes );
    }

@#

  assert( v__tre.size() < 2) ;
  assert( v__tre.back() == SAMPLE_SIZE) ;

@#


/* \newline \aths{\S~\ref{sec:updateri}} */
  updateri( vl,  v__ri );
}



@*1 {\bf approximate $\EE{R_{i}(n)}$}.
\label{sec:approximate}


approximate  $\EE{R_{i}(n)}$ as predicted by the
$\delta_{0}$-Poisson-Dirichlet$(\alpha,0)$ coalescent with transition
rates \eqref{eq:1} given sample size $n$;  the goal here is
\begin{enumerate}
\item given sample size generate all possible merger sizes and
corresponding transition rates   \S~\ref{sec:allmergers}
\item for each of $M = |EXPERIMENTS|$     experiments :
\begin{enumerate}
\item sample a realisation $\ell_{1},\ldots, \ell_{n-1}$   of $L_{1}(n), \ldots, L_{n-1}(n)$
\S~\ref{sec:onexperiment}
\item update $r_{i} \leftarrow r_{i} + \ell_{i}/\sum_{j}\ell_{j}$ \S~\ref{sec:updateri}
\end{enumerate}
\item return $r_{i}/M$ as approximation of  $\EE{R_{i}(n)}$ for
$i=1,2,\ldots, n-1$
\end{enumerate}


@<go ahead -- approximate $\EE{R_{i}(n)}$@>=@#
static void approximate()
{

@#


 /* \newline \aths{ |vri| holds $r_{1},\ldots, r_{n-1}$;  |v__l| holds
 $\ell_{1},\ldots, \ell_{n-1}$;  |v_l_n| holds  $\lambda_{n}$ for
 $n=2,\ldots, |SAMPLE_SIZE|$; |a__cmfs| holds the cumulative mass
 functions for sampling merger size(s) } */

  std::vector<double> vri (SAMPLE_SIZE) ;


@#

  vri.reserve(SAMPLE_SIZE) ;

@#

  std::vector<double> v__l (SAMPLE_SIZE);


@#

  v__l.reserve(SAMPLE_SIZE) ;

@#

  std::vector<double> v_l_n (SAMPLE_SIZE + 1) ;

@#

  v_l_n.reserve(SAMPLE_SIZE + 1) ;


@#

  std::vector< std::vector< double > > a__cmfs (SAMPLE_SIZE + 1, std::vector<double> {} ) ;

 /* \newline \aths{\S~\ref{sec:allmergers}} */
  allmergers(v_l_n, a__cmfs );  
  
  int r = EXPERIMENTS + 1 ;

  while( --r > 0)
    {
    /* \newline \aths{\S~\ref{sec:onexperiment}} */
      onexperiment( vri, v__l, v_l_n, a__cmfs);
    }

  for( unsigned int i = 1 ; i < SAMPLE_SIZE; ++i)
   { std::cout << log( vri[i] / static_cast<double>(EXPERIMENTS)) - log(1. - (vri[i]/static_cast<double>(EXPERIMENTS))) << '\n';}
}





@*1 {\bf main}.
\label{sec:main}


the |main| module; here we
\begin{enumerate}
\item initialise the GSL  random number generator \S~\ref{sec:rngs}
\item approximate   $\EE{R_{i}(n)}$  \S~\ref{sec:approximate}
\end{enumerate}


@C

/* \newline \aths{\S~\ref{sec:includes}} */
@<includes@>@#
/* \newline \aths{\S~\ref{sec:rngs}} */
@<rngs@>@#
/* \newline \aths{\S~\ref{sec:descfact}} */
@<descending factorial@>@#
/* \newline \aths{\S~\ref{sec:power}} */
@<guarded $x^{y}$ module@>@#
/* \newline \aths{\S~\ref{sec:lambdanks}} */
@<compute  $\lambda_{n;k_{1},\ldots,k_{r};s}$@>@#
/* \newline \aths{\S~\ref{sec:msumm}} */
@<$r$-merger summing to $m$@>@#
/* \newline \aths{\S~\ref{sec:mergerstom}} */
@<merger sizes suming to $m$@>@#
/* \newline \aths{\S~\ref{sec:mergerstofile}} */
@<print merger sizes to file@>@#
/* \newline \aths{\S~\ref{sec:allmergersm}} */
@<possible merger sizes when $m$ blocks@>@#
/* \newline \aths{\S~\ref{sec:allmergers}} */
@<allmergers@>@#
/* \newline \aths{\S~\ref{sec:updatetree}} */
@<update the tree configuration@>@#
/* \newline \aths{\S~\ref{sec:updatelengths}} */
@<update $\ell_{i}(n)$@>@#
/* \newline \aths{\S~\ref{sec:updateri}} */
@<update $r_{i}$@>@#
/* \newline \aths{\S~\ref{sec:samplemerger}} */
@<get merger size(s)@>@#
/* \newline \aths{\S~\ref{sec:readmergersizes}} */
@<read in merger size(s) from file@>@#
/* \newline \aths{\S~\ref{sec:onexperiment}} */
@<get $\svigi{\ell_{1}, \ldots, \ell_{n-1}}$@>@#
/* \newline \aths{\S~\ref{sec:approximate}} */
@<go ahead -- approximate $\EE{R_{i}(n)}$@>@#

int main(int argc, char *argv[])
{

@#



/* \newline \aths{\S \ref{sec:rngs}} */
setup_rng(  static_cast<unsigned long int>(atoi(argv[1])) );

@#

/* \newline \aths{\S~\ref{sec:approximate}} */
approximate();

@#

gsl_rng_free( rngtype ); 

return 0 ;

}



@* {\bf conclusions and bibliography}.
\label{sec:concl}


We approximate $\EE{R_{i}(n)}$  when associated with the
$\delta_{0}$-Poisson-Dirichlet$(\alpha,0)$ coalescent with transition
rates \eqref{eq:1}; see \S~\ref{sec:code} for a summary of the
algorithm.  Figure \ref{fig:graph} records an example.






\begin{SCfigure}[0.8][htb]
    \begin{tikzpicture}
      \begin{axis}[
        xlabel = $\log(i/n) - \log(1 - i/n)$,
        axis line style = {draw = none},
        tick style = {draw=none},
        xticklabels={draw=none},
        yticklabels={draw=none},
        legend pos=north east]
        \addplot table[col sep=comma] {forplottingfile1};
        \addplot table[col sep=comma] {forplottingfile2};
        \addlegendentry{$\alpha = 0.05$};
        \addlegendentry{$\alpha = 0.95$};
       \end{axis}
       \end{tikzpicture}
       \caption{\it  An example approximation of $\EE{R_{i}(n)}$
graphed as logits  as a function of $\log(i/n) - \log(1 - i/n)$ for 
$i = 1,2,\ldots, n-1$ where $n$ is sample size }
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
JA~Chetwyn-Diggle, Bjarki Eldon.
\newblock Beta-coalescents when sample size is large.
\newblock \url{https://doi.org/10.64898/2025.12.30.697022}

\bibitem[DK99]{DK99}
P~Donnelly and T~G Kurtz.
\newblock Particle representations for measure-valued population models.
\newblock {\em Ann Probab}, 27:166--205, 1999.
\\
\url{https://doi.org/10.1214/aop/1022677258}

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
