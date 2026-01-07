%$ NAFN=beta_poisson_dirichlet_usingrates
%$ echo $NAFN
%$ echo "constexpr unsigned int SAMPLE_SIZE = 50;" > $NAFN.hpp
%$ echo "constexpr double ALPHA = 0.99;" >> $NAFN.hpp
%$ echo "constexpr double BETA = 1.5 ;" >> $NAFN.hpp
%$ echo "constexpr double CEPS = 1.0e4 ;" >> $NAFN.hpp
%$ echo "constexpr double MINF = (2. + (1 + pow(2., 1-BETA)/(BETA-1)))/2.;" >> $NAFN.hpp
%$ echo "constexpr double CBETA =BETA*tgamma(2.-BETA)*tgamma(BETA)/pow(MINF,BETA);" >> $NAFN.hpp
%$ echo "constexpr int EXPERIMENTS = 1e5 ;" >> $NAFN.hpp
%$ ctangle beta_poisson_dirichlet_usingrates.w
%$ g++ -std=c++26 -m64 -march=native -O3 -x c++ -DNDEBUG $NAFN.c -lm -lgsl -lgslcblas
%$ ./a.out $(shuf -i 43434-2929101 -n1) > tmpresout
%$ sed '1d' tmpresout | awk -F\n '{M=1e5; for (i=1;i<=NF;i++) print log($i/M) - log(1-($i/M))}' > logitresout
%$ seq 49 | awk '{n=50;print log($1/n) - log(1 - ($1/n))}' > nlogits
%$ paste -d',' nlogits logitresout > forplottingfile1
%$ cweave beta_poisson_dirichlet_usingrates.w
%$ tail -n4 beta_poisson_dirichlet_usingrates.tex > endi
%$ for i in $(seq 5); do $(sed -i '$d' beta_poisson_dirichlet_usingrates.tex) ; done
%$ cat endi >> beta_poisson_dirichlet_usingrates.tex
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
%$ lualatex beta_poisson_dirichlet_usingrates
%$ lualatex beta_poisson_dirichlet_usingrates
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
%%\usepackage{xunicode}
\usepackage{fontspec}
\usepackage{xltxtra}
\usepackage{lineno}
\usepackage{amsfonts, amsmath, amssymb}
\usepackage{fullpage}
\usepackage{marvosym}
\usepackage{bm}
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
  {\textsf{\large \textbf{\@@title}}}
\end{flushleft}

\begin{center}
\textsc{\large \@@author}
\end{center}
\egroup
}
\makeatother
\title{Gene genealogies in a  haploid panmictic population evolving
according to  sweepstakes reproduction -- \\ sampling the Beta$(2-\beta,\beta)$-Poisson-Dirichlet$(\alpha,0)$-coalescent}
 \author{Bjarki Eldon\footnote{\href{beldon11@@gmail.com}{beldon11@@gmail.com}}\footnote{ compiled @@
{\DTMcurrenttime} on  {\today}  \\
\input{innleggctangle} \\
 \input{innleggcpp} \\ kernel  \input{innleggop} \\  \input{innleggbash} \\
GSL \input{innlegggsl} \\ \input{innleggcweave} \\  
\input{innlegglualatex}  \\ \input{innleggspix} \\  written using  \input{innleggemacs}}
\orcidlink{0000-0001-9354-2391}}
\begin{document}
\maketitle
\renewcommand{\abstractname}{\vspace{-\baselineskip}}


\begin{abstract}
 Let $\{\xi^{n}\} \equiv \{ \xi^{n}(t) : t \ge 0 \}$ be 
the Beta$(2-\beta,\beta)$-Poisson-Dirichlet$(\alpha,0)$-coalescent.
 Write  $\# A$ for the number of elements
in a given finite set $A$, 
$L_{i}(n) \equiv \int_{0}^{\tau(n) } \# \left\{ \xi \in \xi^{n}(t)
: \#\xi = i \right\}dt $ and   $L(n) \equiv \int_{0}^{\tau(n)} \# \xi^{n}(t)dt $  and 
$\tau(n) \equiv \inf \left\{ t \ge 0 : \# \xi^{n}(t) = 1 \right\} $
for $i \in \{1, 2, \ldots, n-1\}$.   We then have
$L(n) = L_{1}(n) + \cdots + L_{n-1}(n)$. Define 
$R_{i}(n) \equiv L_{i}(n)/\sum_{j}L_{j}(n) $ for $i=1,2,\ldots, n-1$.
Interpreting $\{\xi^{n}\}$ as `trees'  we may view   $L_{i}(n)$  as the random total length of
branches supporting $i \in \{1, 2, \ldots, n-1\}$ leaves, with the
length measured in coalescent time units,  and $n$ sample size.
 With this C++ code we sample the
Beta$(2-\beta,\beta)$-Poisson-Dirichlet$(\alpha,0)$-coalescent with $0
 < \alpha < 1$ and $1 < \beta < 2$ fixed, and approximate $\EE{R_{i}(n)}$
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



Use the shell tool {\tt spix} on the  script appearing before  the
preamble (the lines starting with  \%\$); simply 

{\tt spix /path/to/the/sourcefile}

where {\tt sourcefile} is the {\tt {.w}} file

One   may also  copy the script into a file and run  {\tt parallel}
~\cite{tange11:_gnu_paral}:

{\tt parallel --gnu -j1 :::: /path/to/scriptfile}



@* {\bf intro}.
\label{sec:intro}

Write $\IN \equiv \set{1,2,\ldots}$ and
$\IN_{0} \equiv \set{0,1,2,\ldots}$ and
$[n] \equiv \set{1,2,\ldots,n}$ for all $n\in \IN$.  Fix
$0 < \alpha < 1 $ and $1 < \beta < 2$.  The
Beta$(2-\beta,\beta)$-Poisson-Dirichlet$(\alpha,0)$-coalescent
coalescent is a Markov-chain on the partitions on $[n]$ with
transition rates
\begin{subequations}
\begin{align}
C_{\beta} & = \frac{\beta  }{m_{\infty}^{\beta}}B(2-\beta,\beta) \label{eq:2}  \\
\lambda_{n;k_{1}, \ldots, k_{r}; s} & = \one{r=1}\frac{ C_{\beta}}{C_{\beta} + c(1-\alpha) }B(k-\beta,n-k+\beta) + \frac{cp_{n;k_{1},\ldots, k_{r};s} }{C_{\beta} + c(1-\alpha)}   \label{eq:3}
\end{align}
\end{subequations}
where the transition is the merging of blocks in $r$ groups of size(s)
$k_{1},\ldots, k_{r}$. We approximate $m_{\infty}$ in \eqref{eq:2}
with $\svigi{2 + \svigi{1 + 2^{1-\beta}}/(\beta - 1)}/2$. In
\eqref{eq:3} $p_{n;k_{1},\ldots, k_{r};s} $ is given by
\begin{equation}
\label{eq:4}
p_{n;k_{1},\ldots, k_{r};s} =  \frac{\alpha^{r+s-1}\svigi{r+s-1}!}{(n-1)!}\prod_{i=1}^{r}\svigi{k_{i} - 1 - \alpha}_{k_{i}-1}
\end{equation}
where ${(x)}_{m} \equiv x(x-1)\cdots (x-m+1)$ and ${(x)}_{0} \equiv 1$ for
$x$ real and $m \in \IN_{0}$.   


The background for the
Beta$(2-\beta,\beta)$-Poisson-Dirichlet$(\alpha,0)$-coalescent is
this \cite{e2025+}. 
Fix $0 < \alpha < 1$ and $1 < \beta < 2$.  Consider a haploid  population of
constant finite  size $N$ 
evolving in non-overlapping generations.  In each generation the
current individuals independently produce potential offspring
according to
\begin{equation}
\label{eq:1}
\prb{X=k} = C\svigi{k^{-a} - {(k+1)}^{-a}}
\end{equation}
for all $k = 1,2,\ldots, \zeta(N)$ 
where $a > 0$ and $C$ is such that $\prb{1 \le X \le \zeta(N)} = 1$,
and  and such that  $\zeta(N)/N^{1/\alpha} \to \infty$
as $N\to \infty$.
Let $X_{1},\ldots,X_{N}$ be the iid random numbers of potential
offspring produced by the current $N$ individuals.  Write
$X_{1},\ldots, X_{N} \vartriangleright \mathds L(a,\zeta(N))$ when
have law \eqref{eq:1}.       Take
$\svigi{\varepsilon_{N}}_{N}$ a positive sequence where $0  <
\varepsilon_{N}  < 1$. It may hold that $\varepsilon_{N}\to 0$ as
$N\to \infty$.  With probability  $\varepsilon_{N}$ it holds that
$a=\alpha$, and with probability $1 - \varepsilon_{N}$ it holds that
$a=\beta$, where $a$ is the one in \eqref{eq:1},
\begin{displaymath}
X_{1},\ldots, X_{N} \vartriangleright \mathds L\svigi{ \one{E}\alpha + \one{E^{\sf c}}\beta,  \zeta(N)  }
\end{displaymath}
where $E$ is the event that $a=\alpha$.  
Taking     $\varepsilon_{N} \overset{c}\sim c_{N}$ where $N^{\beta - 1}c_{N}
\overset c \sim 1$ as $N\to \infty$   leads to the 
coalescent with transition rates \eqref{eq:3}. 


The algorithm is summarised in \S, the code follows in \S~\ref{sec:includes}--\S~\ref{sec:main}, we
conclude in \S~\ref{sec:concl}. Comments within the code in \aths{this
font and color}



@* {\bf code}. 
\label{sec:code}


we briefly summarise the algorithm, let $\lambda_{m}$ denote the
total rate of merging $m$ blocks  obtained by summing the rates
\eqref{eq:3} over all possible  ordered merger sizes given $m$ blocks 

\begin{enumerate}
\item $r_{i}\leftarrow 0$ for $i = 1,2,\ldots, n-1$
\item for each of $M$ experiments: \S~\ref{sec:approximate}
\begin{enumerate}
\item $m \leftarrow n$
\item $\ell_{i}\leftarrow 0$ for all  $i \in [n-1]$
\item the current block sizes  $b_{i}\leftarrow 1$ for all $i \in [n]$
\item {\bf while} $m > 1$ : \S~\ref{sec:onexperiment}
\begin{enumerate}
\item sample time until next merger $t \leftarrow $ Exp$(\lambda_{m})$
\item update $\ell_{b} \leftarrow \ell_{b}  + t$ for all $b\in
\set{b_{1},\ldots, b_{m}}$ \S~\ref{sec:updatelengths}
\item sample merger size(s) \S~\ref{sec:samplemerger}
$k_{1},\ldots,k_{r}$ and merge blocks \S~\ref{sec:updatetree}
\item $m\leftarrow m - k_{1} - \cdots - k_{r} + r$ 
 \end{enumerate}
 \item $r_{i} \leftarrow l_{i}/\sum_{j}\ell_{j}$ for all $i\in [n-1]$ \S~\ref{sec:updateri}
\end{enumerate}
\item return $\overline \varrho_{i}(n)$ as   $r_{i}/M$  for all $i\in [n-1]$
\end{enumerate}



@*1 {\bf the includes}.
\label{sec:includes}


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
#include <gsl/gsl_sf.h>
#include "beta_poisson_dirichlet_usingrates.hpp"



@*1 {\bf random number generators}.
\label{sec:rngs}


defining the random number generators

\begin{Piton}
gsl_rng * rngtype ;
static void setup_rng(unsigned long int s)
{
        const gsl_rng_type *T ; 
        gsl_rng_env_setup(); 
        T = gsl_rng_default ;
        rngtype = gsl_rng_alloc(T);
        gsl_rng_set( rngtype,  s) ;
}
\end{Piton}


@<rngs@>=@#
  std::random_device randomseed;
  /* \newline \aths{ Standard mersenne twister  random number engine} */
  std::mt19937_64 rng(randomseed());

gsl_rng * rngtype ;
static void setup_rng(unsigned long int s)
{
        const gsl_rng_type *T ; 
        gsl_rng_env_setup(); 
        T = gsl_rng_default ;
        rngtype = gsl_rng_alloc(T);
        gsl_rng_set( rngtype,  s) ;
}



@*1 {\bf descending factorial}.
\label{sec:descending}


the descending factorial ${(x)}_{m}$

\begin{Piton}
static double descending_factorial(const double x, const double m)
{
  double p = 1;

  for( double i = 0; i < m; ++i){
    p *= (x - i); }

  return p ;
}
\end{Piton}


@<descending@>=@#
static double descending_factorial(const double x, const double m)
{
  double p = 1;

  for( double i = 0; i < m; ++i){
    p *= (x - i); }

  return p ;
}


@*1 {\bf veldi}. 
\label{sec:veldi}


compute $x^{y}$ checking for underflow and overflow

\begin{Piton}
static double veldi (const double x, const double y)
  {
    feclearexcept (FE_ALL_EXCEPT);
    const double d = pow (x,y);

    return (fetestexcept (FE_UNDERFLOW) ? 0. : (fetestexcept (FE_OVERFLOW) ? FLT_MAX : d));
  }
\end{Piton}



@<veldi@>=@#
static double veldi (const double x, const double y)
  {
    feclearexcept(FE_ALL_EXCEPT);
    const double d = pow(x,y);

    return( fetestexcept(FE_UNDERFLOW) ? 0. : (fetestexcept(FE_OVERFLOW) ? FLT_MAX : d) ) ;   
  }


@*1 {\bf the beta part in $\lambda_{n;k_{1},\ldots, k_{r};s}$}. 
\label{sec:betapart}


compute the beta part in  $\lambda_{n;k_{1},\ldots, k_{r};s}$
\eqref{eq:3}, the
\begin{displaymath}
\frac{ C_{\beta}}{C_{\beta} + c(1-\alpha) }B(k-\beta,n-k+\beta)
\end{displaymath}
part

\begin{Piton}
double betapart (const double n, const double k)
{
  return  (CBETA * gsl_sf_beta (k - BETA,  n - k + BETA)) / (CBETA +  CEPS*(1-ALPHA)) ;
}
\end{Piton}


@<beta part@>=@#
double betapart (const double n, const double k)
{
/* \newline \aths{|CBETA| is $C_{\beta}$; |BETA| is $\beta$; |ALPHA|
is $\alpha$} */ 
  return  (CBETA * gsl_sf_beta(k - BETA,  n - k + BETA)) / (CBETA + (CEPS * (1-ALPHA)));
}


@*1 {\bf $\lambda_{n;k_{1},\ldots, k_{r};s}$}. 
\label{sec:lambdanks}

compute   $\lambda_{n;k_{1},\ldots, k_{r};s}$ ~\eqref{eq:3} and
multiply by
\begin{displaymath}
\binom{n}{k_{1} \ldots k_{r} s }\frac{1}{\prod_{j=1}^{n}\svigi{ \sum_{i}\one{k_{i} = j}}! }
\end{displaymath}

@<lambdanks@>=@#
static double lambdanks (const double n, const std::vector<unsigned int>& v_k)
{

  assert(v_k[0] > 1);

  
  double d {};
  double k {} ;
  double f {1} ;
  const double r = static_cast<double>(v_k.size());
  
  std::unordered_map<unsigned int, unsigned int> counts {} ;
  
  for( std::size_t i = 0; i < v_k.size(); ++i){
  /* \newline \aths{\S~\ref{sec:descending}} */
    f *= descending_factorial( static_cast<double>(v_k[i]) - 1. - ALPHA,  static_cast<double>(v_k[i]) -1);
    /* \newline \aths{ count occurrence of each  merger size}  */
    ++counts[v_k[i]];
    k += static_cast<double>(v_k[i]) ;
    d += lgamma(static_cast<double>(v_k[i] + 1)) ; }

  assert( k < n + 1 ); 
  const double s = n - k;

@#

  const double p = static_cast<double>(std::accumulate( counts.begin(), counts.end(), 0, [](double a, const auto &x){ return a +  lgamma( (double)x.second + 1);}));


@#
/* \newline \aths{|betapart| \S~\ref{sec:betapart}} */
  const double l = ((v_k.size() < 2 ? 1. : 0)*betapart(n, v_k[0])) + (CEPS * veldi(ALPHA, r+s-1) * tgamma(r+s)   *f/tgamma(n) ) ;


@#

  return ( veldi(exp(1), (lgamma( n + 1.) - d)  - lgamma( n - k + 1 ) - p) * l / (CBETA + (CEPS*(1-ALPHA))) ) ;
}



@*1 {\bf merger sizes}.
\label{sec:GenPartitions}


generate merger size(s) of $m$ groups summing to |myInt| 
 
@<merger sizes@>=@#
static double  GenPartitions( const unsigned int m,   const unsigned int myInt,
                   const unsigned int PartitionSize,
                   unsigned int MinVal,
			  unsigned int MaxVal,
			  std::vector< std::pair< double, std::vector<unsigned int> > > & v_l_k,
			  std::vector<double>& lrates_sorting)
{

  /* \newline \aths{ |m| is the given number of blocks; the partitions sum to  |myInt|} */
  
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
@#
            case 5: for (++idx_Dec, LeftRemain = (idx_Dec - idx_Spill)*val_Dec; (idx_Dec <= idx_Last) && (partition[idx_Dec] <= MinVal); idx_Dec++) 
                        LeftRemain += partition[idx_Dec];
@#
                    LeftRemain += 1 + (idx_Spill - idx_Dec + 1)*MinVal;
                    break;
@#
            case 6:
            case 7:
            case 11:idx_Dec = idx_Spill + 1;
                    LeftRemain += 1 + (idx_Spill - idx_Dec + 1)*MinVal;
                    break;

@#
            case 9: for (++idx_Dec, LeftRemain = idx_Dec * val_Dec; (idx_Dec <= idx_Last) && (partition[idx_Dec] <= (val_Dec + 1)); idx_Dec++)
                        LeftRemain += partition[idx_Dec];

                    LeftRemain += 1 - (idx_Dec - 1)*MinVal;
                    break;
@#
            case 10:for (LeftRemain += idx_Spill * MinVal + (idx_Dec - idx_Spill)*val_Dec + 1, ++idx_Dec; (idx_Dec <= idx_Last) && (partition[idx_Dec] <= (val_Dec - 1)); idx_Dec++)   
                        LeftRemain += partition[idx_Dec];

                    LeftRemain -= (idx_Dec - 1)*MinVal;
                    break;
        }
@#
        while (idx_Spill > idx_SpillPrev) 
            partition[--idx_Spill] = MinVal;

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




@*1 {\bf all mergers with fixed sum}. 
\label{sec:allmergers_sum_m}


get all mergers  summing to |m|


@<mergers fixed sum@>=@#
static double allmergers_sum_m( const unsigned int n,  const unsigned int m, std::vector< std::pair< double, std::vector<unsigned int> > >& v__l_k,  std::vector<double>&  v_lrates_sort )
{
  /* \newline \aths{ |n| is the number of blocks;
     the partitions sum to |m|} */
  
  const std::vector<unsigned int> v__m {m};
  /* \newline \aths{\S~\ref{sec:lambdanks}; first record
  $\lambda_{n;m;n-m}$ the rate of a single merger of size $m$} */
  double sumr = lambdanks( static_cast<double>(n), v__m ) ;

  v__l_k.push_back( std::make_pair( sumr, v__m ) ) ;

  v_lrates_sort.push_back( sumr );

  if( m > 3){
    for( unsigned int s = 2; s <= m/2; ++s){
      
      assert(m > 2*(s-1) );
      /* \newline \aths{\S~\ref{sec:GenPartitions}; add the rates of
      simultaneous mergers} */
      sumr += GenPartitions(n, m, s, 2, m - (2*(s-1)), v__l_k, v_lrates_sort ); }
  }
  assert( sumr >= 0) ;
  return sumr ;
}


@*1 {\bf record merger sizes in order}.
\label{sec:ratesmergersfile}


record ordered merger sizes  

@<record merger sizes in order@>=@#
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
      cmf += (vlk[ i ].first) / s ;
      assert( cmf >= 0) ;
@#    
      a__cmf[n].push_back(cmf) ;

@#

      assert((vlk[ i ].second).size() > 0 );


@#

      for( const auto &x : vlk[ i ].second){
	f << x << ' ' ;}
      f << '\n' ; }

@#

  f.close();

@#

  assert( abs(cmf - 1.) < 0.999999 );
}


@*1 {\bf all mergers when $n$ blocks}. 
\label{sec:allmergers_when_n_blocks}


get all mergers when a given number of blocks 


@<mergers when $n$ blocks@>=@#
static void allmergers_when_n_blocks (const unsigned int n, std::vector<double> & v__lambdan, std::vector< std::vector< double > > & a__cmf)
{
  std::vector< std::pair< double, std::vector<unsigned int> > > vlk {} ;
  std::vector< double > ratetosort {} ;
  ratetosort.clear() ;

  double lambdan {} ;
  vlk.clear() ;
  assert( n > 1);
  for( unsigned int k = 2 ; k  <= n ; ++k){
    /* \newline  \aths{the partition sums to |k|; the number of blocks
    is |n|; \S~\ref{sec:allmergers_sum_m}} */
    lambdan += allmergers_sum_m(n, k, vlk, ratetosort); }
  /* \newline \aths{ record the total rate when |n| blocks; use for sampling time} */
  assert( lambdan > 0);
  v__lambdan[n] = lambdan ;

  std::vector<unsigned int> indx (ratetosort.size());
  std::iota(indx.begin(), indx.end(), 0);

  std::stable_sort( indx.begin(), indx.end(), [&ratetosort](const unsigned int x, const unsigned int y){return ratetosort[x] > ratetosort[y];});
  
  /* \newline \aths{ merger rates sorted in descending order;  print
  the cmf and rates to file; \S~\ref{sec:ratesmergersfile}} */

  ratesmergersfile(n, indx, vlk,  v__lambdan[n], a__cmf);
}



@*1 {\bf all mergers up to sample size}.
\label{sec:allmergers}


 generate all mergers up to sample size 


@<all mergers@>=@#
static void allmergers( std::vector<double>& vlmn, std::vector< std::vector<double> > & acmf  )
{
   
  for( unsigned int tmpn = 2; tmpn <= SAMPLE_SIZE; ++tmpn)
    {
      /* \newline \aths{\S~\ref{sec:allmergers_when_n_blocks}} */
      allmergers_when_n_blocks (tmpn, vlmn, acmf );
    }
}



@*1 {\bf update tree}.
\label{sec:updatetree}


update block sizes


@<update block sizes@>=@#
static void  updatetree( std::vector<unsigned int> & tre,   const std::vector<unsigned int>& mergersizes )
{

@#

  assert( mergersizes.size() > 0) ;


@#

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
      /* \newline \aths{  remove the blocks that merged} */
      tre.resize(s);
    }
  assert( newblocks.size() > 0) ;
  assert( static_cast<unsigned int>(std::accumulate( newblocks.begin(), newblocks.end(), 0)) <= SAMPLE_SIZE) ;
  
  tre.insert( tre.end(), newblocks.begin(), newblocks.end() ) ;

  assert( static_cast<unsigned int>(std::accumulate( tre.begin(), tre.end(), 0)) == SAMPLE_SIZE) ;
}



@*1 {\bf update lengths}. 
\label{sec:updatelengths}


update the lengths $\ell_{i}$ given the current block sizes


@<lengths@>=@#
static void updatelengths( const std::vector<unsigned int>& tre,  std::vector<double>& v__lengths,  const std::vector<double>& v_lambdan)
{

/* \newline \aths{get the time until merger} */
  const double t = gsl_ran_exponential(rngtype, 1./v_lambdan[tre.size()] );
  for( const auto &b: tre)
    {
    @#
      assert( b > 0) ;
      @#
      assert( b < SAMPLE_SIZE) ;
      @#
      v__lengths[0] += t ;
      v__lengths[b] += t ;
    }
}



@*1 {\bf update approximations $r_{i}$}.
\label{sec:updateri}

update the approximations $r_{i}$


@<$r_{i}$@>=@#
static void updateri( const std::vector<double>& v__l, std::vector<double>& v_ri )
{
@#
  assert( v__l[0] > 0);
  @#
  const double d = v__l[0];
  std::transform(v__l.begin(), v__l.end(), v_ri.begin(), v_ri.begin(), [&d](const auto &x, const auto &y){return y + (x/d);});
}



@*1 {\bf sample merger sizes}.
\label{sec:samplemerger}

get merger size(s)


@<get merger sizes@>=@#
static unsigned int samplemerger(const unsigned int n, const std::vector< double > & v__cmf )
{
  unsigned int j {} ;
  const double u = gsl_rng_uniform( rngtype);

  while( u > v__cmf[j]){ ++j;}

/* \newline \aths{|j| is the sampled index of the  ordered  merger size(s)} */
  return j ;
}



@*1 {\bf read merger sizes}.
\label{sec:readmergersizes}


@<read in merger sizes@>=@#
static void readmergersizes(const unsigned int n,  const unsigned int j, std::vector<unsigned int> & v__mergers )
{
@#
  std::ifstream f("gg_" + std::to_string(n) + "_.txt") ;
  std::string line {} ;

  v__mergers.clear(); 

@#
  for( unsigned int i = 0; std::getline(f, line) && i < j; ++i ){
    if(i >= j-1 ){
      
      std::stringstream ss(line) ;
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

generate one experiment

@<one experiment@>=@#
static void onexperiment( std::vector<double>& v__ri,  std::vector<double>& vl, const std::vector<double>& v__lambda__n, const std::vector< std::vector< double > > & a__cmf )
{
@#
  std::vector<unsigned int> v__tre(SAMPLE_SIZE, 1) ;
  @#
  std::fill( vl.begin(), vl.end(), 0);
  @#
  unsigned int lina {} ;

@#


  std::vector<unsigned int> v_merger_sizes (SAMPLE_SIZE/2);
  v_merger_sizes.reserve(SAMPLE_SIZE/2) ;
  @#
  while( v__tre.size() > 1)
    {
    /* \newline \aths{\S~\ref{sec:updatelengths}} */
      updatelengths( v__tre, vl, v__lambda__n );
      /* \newline \aths{\S~\ref{sec:samplemerger}} */
      lina = samplemerger( v__tre.size(),  a__cmf[v__tre.size()] );
      /* \newline \aths{\S~\ref{sec:readmergersizes}} */
      readmergersizes( v__tre.size(), 1 + lina, v_merger_sizes) ;
      /* \newline \aths{\S~\ref{sec:updatetree}} */
      updatetree( v__tre, v_merger_sizes );
    }
  assert( v__tre.back() == SAMPLE_SIZE) ;

/* \newline \aths{\S~\ref{sec:updateri}} */
  updateri( vl,  v__ri );
}


@*1 {\bf approximate}.
\label{sec:approximate}



@<go ahead -- get $\overline \varrho_i(n)$@>=@#
static void approximate()
{
@#
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

 
  for(const auto &z: vri){std::cout << z << '\n';}
}






@*1 {\bf main}.
\label{sec:main}

the |main| function


@C


/* \newline \aths{\S~\ref{sec:includes}} */
@<includes@>@#
/* \newline \aths{\S~\ref{sec:rngs}} */
@<rngs@>@#
/* \newline \aths{\S~\ref{sec:descending}} */
@<descending@>@#
/* \newline \aths{\S~\ref{sec:veldi}} */
@<veldi@>@#
/* \newline \aths{\S~\ref{sec:betapart}} */
@<beta part@>@#
/* \newline \aths{\S~\ref{sec:lambdanks}} */
@<lambdanks@>@#
/* \newline \aths{\S~\ref{sec:GenPartitions}} */
@<merger sizes@>@#
/* \newline \aths{\S~\ref{sec:allmergers_sum_m}} */
@<mergers fixed sum@>@#
/* \newline \aths{\S~\ref{sec:ratesmergersfile}} */
@<record merger sizes in order@>@#
/* \newline \aths{\S~\ref{sec:allmergers_when_n_blocks}} */
@<mergers when $n$ blocks@>@#
/* \newline \aths{\S~\ref{sec:allmergers}} */
@<all mergers@>@#
/* \newline \aths{\S~\ref{sec:updatetree}} */
@<update block sizes@>@#
/* \newline \aths{\S~\ref{sec:updatelengths}} */
@<lengths@>@#
/* \newline \aths{\S~\ref{sec:updateri}} */
@<$r_{i}$@>@#
/* \newline \aths{\S~\ref{sec:samplemerger}} */
@<get merger sizes@>@#
/* \newline \aths{\S~\ref{sec:readmergersizes}} */
@<read in merger sizes@>@#
/* \newline \aths{\S~\ref{sec:onexperiment}} */
@<one experiment@>@#
/* \newline \aths{\S~\ref{sec:approximate}} */
@<go ahead -- get $\overline \varrho_i(n)$@>@#

int main(int argc, char * argv[])
{

/* \newline \aths{\S~\ref{sec:rngs}} */
   setup_rng( static_cast<unsigned long>( atoi(argv[1]))) ;
/* \newline \aths{\S~\ref{sec:approximate}} */
approximate();

@#
gsl_rng_free(rngtype);


return 0;
}


@* {\bf conclusions and references}.
\label{sec:concl}


we sample the
Beta$(2-\beta,\beta)$-Poisson-Dirichlet$(\alpha,0)$-coalescent and
approximate  $\EE{R_{i}(n)}$.       An example approximation  in
Figure~\ref{fig:graph1}.  








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
        %\addplot table[col sep=comma] {forplottingfile2};
        %\addplot table[col sep=comma] {forplottingfile3};
        %\addlegendentry{$\gamma = 0.1$}
        %\addlegendentry{$\gamma = 0.5$}
        %\addlegendentry{$\gamma = 1$}
       \end{axis}
       \end{tikzpicture}
       \caption{\textcolor{violet}{ \it  An example approximation of $\EE{R_{i}(n)}$
 for the  given  parameter values and graphed as logits against
 $\log(i/n) - \log(1 - i/n)$ for $i = 1,2,\ldots, n-1$ where $n$ is
 sample size; }}
      \label{fig:graph1}
       \end{SCfigure}







\bibliographystyle{acm}
%%\bibliography{refs}
\begin{thebibliography}{1}

\bibitem{e2025+}
{\sc Eldon, B.}
\newblock Gene genealogies in haploid populations evolving according
to sweepstakes reproduction.
\newblock In preparation, 2025+.


\bibitem{tange11:_gnu_paral}
{\sc Tange, O.}
\newblock {GNU} parallel -- the command-line power tool.
\newblock The USENIX Magazine, 2011.

\end{thebibliography}


@
\end{document}
