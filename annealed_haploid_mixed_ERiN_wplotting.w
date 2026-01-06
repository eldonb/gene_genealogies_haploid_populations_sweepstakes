%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Copyright (C) 2025 Bjarki Eldon
%$ NAFN=annealed_haploid_mixed_ERiN_wplotting
%$ ctangle $NAFN.w
%$ g++ -std=c++26 -m64 -march=native -O3 -x c++ $NAFN.c -lm -lgsl -lgslcblas
%$ rm -f tmp_runs
%$ for i in $(seq 4); do echo "./a.out 50 "  $(shuf -i 39393-282929191 -n1) "> resout"$i >> tmp_runs; done
%$ parallel --gnu -j4 :::: ./tmp_runs
%$ paste resout* -d, | sed '1d' | awk -F, '{s=0; for (i=1;i<=NF;i++) {s+=$i} print log(s/1e4) - log(1-(s/1e4))}' > logitresout
%$ seq 49 | awk '{S=50;print log($1/S) - log(1 - ($1/S))}' > nlogits
%$ paste -d',' nlogits logitresout > forplottingfile1
%$ rm -f tmp_runs
%$ for i in $(seq 4); do echo "./a.out 500 "  $(shuf -i 39393-282929191 -n1) "> resout"$i >> tmp_runs; done
%$ parallel --gnu -j4 :::: ./tmp_runs
%$ paste resout* -d, | sed '1d' | awk -F, '{s=0; for (i=1;i<=NF;i++) {s+=$i} print log(s/1e4) - log(1-(s/1e4))}' > logitresout
%$ seq 499 | awk '{S=500;print log($1/S) - log(1 - ($1/S))}' > nlogits
%$ paste -d',' nlogits logitresout > forplottingfile2
%$ g++ --version | head -n1 > innleggcpp
%$ xelatex --version | head -n1  > innleggxelatex
%$ cweave  --version | head -n1 > innleggcweave
%$ uname  --kernel-release -o  > innleggop
%$ ctangle  --version | head -n1  > innleggctangle
%$ emacs --version | head -n1 > innleggemacs
%$ parallel --version | head -n1 > innleggparallel
%$ awk --version | head -n1 > innleggawk
%$ spix --version > innleggspix
%$ cweave $NAFN.w
%$ tail -n4 $NAFN.tex > endi
%$ for i in $(seq 5); do $(sed -i '$d' $NAFN.tex) ; done
%$ cat endi >> $NAFN.tex
%$ xelatex $NAFN.tex
%$ xelatex $NAFN.tex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\documentclass[a4paper,10pt,final]{cweb}
\usepackage{fancyvrb,  graphicx,  amsmath, amssymb, url,dsfont}
\usepackage[svgnames]{xcolor}
\usepackage{minted}
%% a4 paper size 210 x 297 millimeters
\usepackage[a4paper,text={15.5cm,24.2cm},centering]{geometry}
\usepackage{orcidlink}
%% for xelatex
\usepackage{xunicode}
\usepackage{fontspec}
\usepackage{xltxtra}
\usepackage{lineno}
\usepackage[all]{xy}
\usepackage[]{hyperref}
\hypersetup{
    colorlinks,
    linkcolor={gray!50!black},
    citecolor={blue!50!black},
    urlcolor={blue!80!black}
}
%   \hypersetup{colorlinks=true}
\usepackage{dsfont}
%%\usepackage[round,numbers,super]{natbib}
\usepackage{fancyvrb}
\usepackage{natbib}
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
\newcommand{\set}[1]{\ensuremath{\left\{ #1 \right\}}}
\newcommand{\svigi}[1]{\ensuremath{\left( #1 \right)}}
\newcommand{\EE}[1]{\ensuremath{\mathds {E}\left[ #1 \right]}}
\newcommand{\IN}{\ensuremath{\mathds N}}
\newcommand{\prb}[1]{\ensuremath{\mathds P\left( #1 \right)}}
\newcommand{\one}[1]{\ensuremath{\mathds 1_{\set{#1}} }}
\newcommand{\aths}[1]{\textcolor{violet}{\bf \small  #1}}
\makeatletter
\renewcommand{\maketitle}{\bgroup\setlength{\parindent}{0pt}
\begin{flushleft}
  {\textsf{\textbf{\@@title}}}
\end{flushleft}
\begin{center}
  {\sc \@@author}
\end{center}
\egroup
} \makeatother
\title{Gene genealogies in haploid  populations evolving according to sweepstakes reproduction \\ --- approximating \EE{R_i^N(n)} }
\author{Bjarki Eldon\orcidlink{0000-0001-2391-9354}\footnote{\href{mailto:beldon11@@gmail.com}{beldon11@@gmail.com} \\  compiled @@ {\DTMcurrenttime} on  {\today}  \\ \input{innleggop}\\   \input{innleggctangle}  \\  \input{innleggcpp} \\  \input{innleggcweave} \\    {\XeLaTeX}  \input{innleggxelatex}  \\   \input{innleggparallel}  \\ \input{innleggawk} \\  \input{innleggspix} \\   \input{innleggemacs} } }
\date{}
\begin{document}
\maketitle
\renewcommand{\abstractname}{\vspace{-\baselineskip}}

\begin{abstract}
Write $[n] \equiv \set{1,2,\ldots, n}$ for
$n\in \IN \equiv \set{1,2,\ldots}$,
$\IN_{0} \equiv \set{0,1,2,\ldots} = \IN\cup\set{0}$.  Let
$\set{\xi^{n,N}} \equiv \set{\xi^{n,N}(t), t \ge \IN_{0}}$ on the
partitions of $n$ denote the process keeping track of the random
ancestral relations of a random set of individuals (gene copies)
sampled at an arbitrary time from a haploid panmictic population of
constant finite  size $N$.  The population may be evolving according to  a
model of sweepstakes reproduction (heavy-tailed offspring number
distribution).   Write  $\# A$ for  the number of elements in a given
finite set $A$. Let     $\tau^{N}(n) \equiv \inf \set{t\in \IN_{0} :
\#\xi^{n,N}(t) = 1 } $,    $L_{i}^{N}(n) \equiv
\sum_{t=0}^{\tau^{N}(n)} \#\set{\xi \in \xi^{n,N}(t) : \#\xi = i} $
for all  $i \in [n-1]$, 
$L^{N}(n) \equiv  \sum_{t=0}^{\tau^{N}(n)} \# \xi^{n,N}(t) $. Then
$L^{N}(n) \equiv \sum_{i}L_{i}^{N}(n)$.   Define $R_{i}^{N}(n) \equiv
L_{i}^{N}(n)/L^{N}(n)$ for all $i\in [n-1]$.  Interpreting
$\set{\xi^{n,N}}$ as `trees' we may  view   $R_{i}^{N}(n)$ as the
random relative branch length supporting $i$ leaves (sampled gene
copies).   With this C++ code we use simulations to approximate
$\EE{R_{i}^{N}(n)}$ 
\end{abstract}

\tableofcontents


@* {\bf Copyright}.
\label{sec:copyright}

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
\label{sec:compile}


This CWEB \citep{knuth1994cweb} document (the {\tt .w} file) can be
compiled with {\tt cweave} to generate a {\tt .tex} file, and with
{\tt ctangle} to generate a {\tt .c} \citep{kernighan1988c} C++ code
file. 


Use {\tt valgrind} to check for memory leaks:

{\tt valgrind -v ---leak-check=full ---leak-resolution=high ---num-callers=40 ---vgdb=full <program call>}


Use {\tt cppcheck} to check the code:

{\tt cppchek  ---enable=all ----language=c++ <prefix>.c}


See \S~\ref{sec:compscript} for a compilation script


@*1 {\bf a compilation script}.
\label{sec:compscript}



Since the script below appears before the preamble one may use the
shell tool {\tt spix} on the script below; simply

{\tt spix <name>.w}

\begin{Verbatim}[frame=leftline,numbers=left,fontsize=\small,fontshape=it,formatcom=\color{teal}]
 NAFN=annealed_haploid_mixed_ERiN_wplotting
 ctangle $NAFN.w
 g++ -std=c++26 -m64 -march=native -O3 -x c++ $NAFN.c -lm -lgsl -lgslcblas
 rm -f tmp_runs
 for i in $(seq 4); do echo "./a.out 50 "  $(shuf -i 39393-282929191 -n1) "> resout"$i >> tmp_runs; done
 parallel --gnu -j4 :::: ./tmp_runs
 paste resout* -d, P sed '1d' P awk -F, '{s=0; for (i=1;i<=NF;i++) {s+=$i} print log(s/1e4) - log(1-(s/1e4))}' > logitresout
 seq 49 P awk '{S=50;print log($1/S) - log(1 - ($1/S))}' > nlogits
 paste -d',' nlogits logitresout > forplottingfile1
 rm -f tmp_runs
 for i in $(seq 4); do echo "./a.out 500 "  $(shuf -i 39393-282929191 -n1) "> resout"$i >> tmp_runs; done
 parallel --gnu -j4 :::: ./tmp_runs
 paste resout* -d, P sed '1d' P awk -F, '{s=0; for (i=1;i<=NF;i++) {s+=$i} print log(s/1e4) - log(1-(s/1e4))}' > logitresout
 seq 499 P awk '{S=500;print log($1/S) - log(1 - ($1/S))}' > nlogits
 paste -d',' nlogits logitresout > forplottingfile2
 cweave $NAFN.w
 tail -n4 $NAFN.tex > endi
 for i in $(seq 5); do $(sed -i '$d' $NAFN.tex) ; done
 cat endi >> $NAFN.tex
 xelatex $NAFN.tex
\end{Verbatim}
In the script {\tt P} stands for the shell pipe operator.  For each of
two parameter settings the simulations are run in parallel using {\tt
parallel}, we use {\tt paste}, {\tt sed}, and {\tt awk} to process the
results for plotting using the \LaTeX package {\tt tikz/pgf}; lines
21--23 in the script are just a hack to ensure a smooth compilation
with \XeTeX

One may also save the script in a file (say, compile) and  apply {\tt parallel} \citep{tange11:_gnu_paral}:

{\tt parallel --gnu -j1 :::: ./compile}


Figure~\ref{fig:graph} in \S~\ref{sec:concl}  graphs the results of the simulations generated by the script


@* {\bf intro}. 
\label{sec:intro}

Recalling from the abstract let 
$\set{\xi^{n,N}} \equiv \set{\xi^{n,N}(t), t \ge \IN_{0}}$ on the
partitions of $n$ denote the process keeping track of the random
ancestral relations of a random set of individuals (gene copies)
sampled at an arbitrary time from a haploid panmictic population of
constant size $N$. The population may be evolving according to  a
model of sweepstakes reproduction (heavy-tailed offspring number
distribution). Write  $\# A$ for  the number of elements in a given
finite set $A$. Let     $\tau^{N}(n) \equiv \inf \set{t\in \IN_{0} :
\#\xi^{n,N}(t) = 1 } $,    $L_{i}^{N}(n) \equiv
\sum_{t=0}^{\tau^{N}(n)} \#\set{\xi \in \xi^{n,N}(t) : \#\xi = i} $
for all  $i \in [n-1]$, 
$L^{N}(n) \equiv  \sum_{t=0}^{\tau^{N}(n)} \# \xi^{n,N}(t) $. Then
$L^{N}(n) \equiv \sum_{i}L_{i}^{N}(n)$.   Define $R_{i}^{N}(n) \equiv
L_{i}^{N}(n)/L^{N}(n)$ for all $i\in [n-1]$.  Interpreting
$\set{\xi^{n,N}}$ as `trees' we may  view   $R_{i}^{N}(n)$ as the
random relative branch length supporting $i$ leaves (sampled gene
copies).  We use simulations to estimate $\EE{R_{i}^{N}(n)}$ for $i =
1,2, \ldots, n-1$.   



In each generation the  current $N$ individuals independently produce
{\it potential}  offspring (offspring that may survive to maturity)  according to 
\begin{equation}
\label{eq:1}
\mathds P\left( X = k\right) = C\left( k^{-a} - (1+k)^{-a} \right)
\end{equation}
for $a > 0$ fixed and $C>0$ chosen such that $\prb{X \ge 1} =
1$. We may use a {\it finite} (fixed)  upper bound $\zeta(N)$ on $X$, and
choose $C$ in \eqref{eq:1} such that  $\prb{1 \le X \le \zeta(N)} =
1$.

Write  $X\vartriangleright L(\alpha,\zeta(N))$ when the law of  $X$ is
as in \eqref{eq:1} for some given  $a$ and $\zeta(N)$. Write
\begin{displaymath}
\begin{split}
E & =  \set{X_{i} \vartriangleright L(\alpha,\zeta(N)) \text{ for all } i\in [N]} \\
E^{\sf c} & =    \set{X_{i} \vartriangleright L(\kappa,\zeta(N)) \text{ for all } i\in [N]}
\end{split}
\end{displaymath}
for some fixed $0 < \alpha < 2$ and $\kappa \ge 2$.   Let
$\svigi{\varepsilon_{N}}_{N\in \IN}$  be a fixed positive sequence
where $0 \le \varepsilon_{N} \le  1$ and it may hold that $\varepsilon_N
\to 0$ as $N\to \infty$.   Suppose
\begin{equation}
\label{eq:2}
\prb{E} = \varepsilon_N, \quad \prb{E^{\sf c}} = 1 - \varepsilon_N
\end{equation}
When $0 < \alpha < 1$ we take  $\varepsilon_{N} = c (\log N)/N$.  

In \S~\ref{sec:algo} we summarize the algorithm, the code follows in  \S~\ref{sec:includes} -- \S~\ref{sec:main}, we
conclude in \S \ref{sec:concl}.  Comments within the code are in \aths{this font and colour}


@* {\bf code}.
\label{sec:code}


@*1 {\bf the algorithm}.
\label{sec:algo}


\begin{enumerate}
\item $\svigi{r_1, \ldots, r_{n-1}} \leftarrow \svigi{0,\ldots, 0} $
\item for each of $M$ experiments
\begin{enumerate}
\item $\svigi{\ell_{1}, \ldots, \ell_{n-1}} \leftarrow
\svigi{0,\ldots, 0}$
\item  $\svigi{\xi_{1}, \ldots, \xi_{n}} \leftarrow
\svigi{1,\ldots, 1}$
\item $m\leftarrow n$ the current number of blocks
\item {\bf while } $m > 1$ :
\begin{enumerate}
\item $\ell_{\xi}\leftarrow \ell_{\xi} + 1$ for $\xi = \xi_{1},\ldots,
\xi_{m}$
\item sample $N$ numbers of potential offspring $X_{1},\ldots, X_{N}$ \eqref{eq:1}
\item get merger sizes $k_{1}, \ldots, k_{r}$ by  assigning  blocks to families using  $X_{1},\ldots, X_{N}$ 
\item merge blocks assigned to same family 
\item $m\leftarrow m - \sum_{i}\one{k_{i} \ge 2}k_{i} + \sum_{i}\one{k_{i} \ge 2} $
\end{enumerate}
\item $r_{i} \leftarrow r_{i} +   \ell_{i}/\sum_{j}\ell_{j} $ for $i =
1,2, \ldots, n-1$
\end{enumerate}
\item return $r_{i}/M$ as an approximation of $\EE{R_{i}^{N}(n)}$ for   $i =
1,2, \ldots, n-1$
\end{enumerate}



@*1 {\bf includes}.
\label{sec:includes}

the included libraries


@<includes@>=@#
#include <iostream>
#include <vector>
#include <random>
#include <functional>
#include <memory>
#include <utility>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <list>
#include <string>
#include <fstream>
#include <forward_list>
#include <chrono>
#include <assert.h>
#include <math.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_gamma.h>



@*1 {\bf constants}.
\label{sec:constants}

the global constants

@<constants@>=@#
/* \newline \aths{$\alpha$ in \eqref{eq:2}} */
constexpr double CONST_ALPHA = 0.25 ;
/* \newline \aths{$\kappa$ in \eqref{eq:2}} */
constexpr double KAPPA = 2.0;
/* \newline \aths{population size $N$} */
constexpr size_t CONST_POP_SIZE = 5e2 ;
/* \newline \aths{$\varepsilon_N =  cf(N)$} */
constexpr double Cc = 1.;
/* \newline \aths{$\varepsilon_N$ in \eqref{eq:2}} */
constexpr double EPSILON = Cc*log(static_cast<double>(CONST_POP_SIZE))/static_cast<double>(CONST_POP_SIZE) ;
/* \newline \aths{$\zeta(N)$ in \eqref{eq:1}} */
constexpr double CONST_CUTOFF = floor(pow(static_cast<double>(CONST_POP_SIZE), 0.)*log(static_cast<double>(CONST_POP_SIZE))) ;
/* \newline \aths{number of experiments} */
constexpr double CONST_NUMBER_EXPERIMENTS = 2500. ;
/* \newline \aths{number of experiments for approximating coalescence
probabilities} */
constexpr std::size_t NUMBER_EXPS = 10000 ;



@*1 {\bf the random number generators}.
\label{sec:rngs}

the standard library and GSL random number generators

@<the rngs@>=@#
std::random_device randomseed;
std::mt19937_64 rng(randomseed());

gsl_rng * rngtype ;
static void setup_rng( unsigned long int s )
{
        const gsl_rng_type *T ; 
        gsl_rng_env_setup(); 
        T = gsl_rng_default ;
        rngtype = gsl_rng_alloc(T);
        gsl_rng_set( rngtype, s) ;
}



@*1 {\bf the probability mass function \eqref{eq:1}}. 
\label{sec:pmf}


compute $C \svigi{j^{-a} - (1+j)^{-a}} $ where $C = 1/\svigi{ 1 - (1 +
\zeta(N))^{-a}} $

@<compute \eqref{eq:1}@>=@#
static double massfunction( const double& j)
{

return ( (pow(1.0/j, CONST_ALPHA) -   pow(1.0/(1. + j),
CONST_ALPHA))/( 1. -  pow( 1./(CONST_CUTOFF + 1.), CONST_ALPHA)) ) ;
}


@*1 {\bf generate cumulative mass function of \eqref{eq:1}}.
\label{sec:cmf}


generate the cumulative mass function of \eqref{eq:1};  the mass
function in \eqref{eq:1} is monotone decreasing 
$\prb{X=j+1} \le \prb{X= j}$


@<cmf@>=@#
static void generate_cdf( std::vector<double>& vcdf )
{ @#
  vcdf.clear() ; @#
   vcdf.push_back(0.0) ; @#
  /* \newline  \aths{using |massfunction| in  \S~\ref{sec:pmf} } */
  for( double j = 1; j <= CONST_CUTOFF ; ++j){
    vcdf.push_back( vcdf.back() + massfunction(j) );}

  assert( abs( vcdf.back() - 1.) < 0.999999 );
}


@*1 {\bf sample one random number of potential offspring}.
\label{sec:sampleoneX}

sample the random number of potential offspring produced by one
individual; use when there is an  upper bound $\zeta(N)$,    $\zeta(N) < \infty$ 

@<sample an $X$@>=@#
static size_t sample_juveniles( const std::vector<double>& vcdf  )
{@#
    /* \newline \aths{the CDF is generated in \S~\ref{sec:cmf} } */
  size_t j = 1; @#
  const double u = gsl_rng_uniform( rngtype );  @#
  while( vcdf[j] < u){
    ++j ; }
  return(j) ;
}


@*1 {\bf sample an $X$ when $\zeta(N) = \infty$}.
\label{sec:unboundedX}


sample one random number of potential offspring when $\zeta(N) =
\infty$,  return   $\lfloor U^{-1/a} \rfloor$ where $U$ a standard
random uniform on $(0,1]$ 

@<$X$ when unbounded@>=@#
static std::size_t randomX(const double &a)
{
 
   return static_cast<std::size_t>( floor( 1./pow( gsl_rng_uniform_pos(rngtype), 1./a))) ;
}



@*1 {\bf sample a pool of potential offspring}.
\label{sec:pool}

get a pool of potential offspring produced by the current $N$
individuals


@<pool@>=@#
static std::size_t sample_pool_juveniles( std::vector<size_t>& pool_juvs, const std::vector<double>& v_cdf)
{@#

  const double a = gsl_rng_uniform(rngtype) < EPSILON ? CONST_ALPHA :
  KAPPA ; @#
  /* \newline \aths{|s| will be the total number sampled} */
  std::size_t s = 0; @# 
  /* \newline \aths{use |randomX| \S~\ref{sec:unboundedX}  when
  unbounded;   |x= randomX(a)|; } */
  /* \newline \aths{use |sample_juveniles| \S~\ref{sec:sampleoneX}
  when bounded; |x=  sample_juveniles(v_cdf);|} */
  std::transform(pool_juvs.begin(), pool_juvs.end(), pool_juvs.begin(), [&s, &a](auto &x){x = randomX(a); s += x; return x;}); 

  assert( s >= CONST_POP_SIZE);
  return (s);
}



@*1 {\bf approximate coalescence probabilities}.
\label{sec:estimatecoals}



estimate coalescence probabilities for  speeding up reaching the most
recent common ancestor.  When only two blocks left we can sample a
geometric with success probability the pairwise coalescence
probability.  When only three blocks left can sample between a
pairwise merger and  a triple merger. Given a realisation $x_{1},
\ldots, x_{N}$ of $X_{1},
\ldots, X_{N}$ with $s_{N} := x_{1} + \cdots + x_{N}$  the pairwise coalescence
probability is


\begin{equation}
\label{eq:3}
c_{N} =  \sum_{j=1}^{N} \frac{x_{j}(x_{j} - 1)}{s_{N}(s_{N}-1)},
\end{equation}
a 3-merger when three blocks is
\begin{equation}
\label{eq:4}
c_{N}(3;3) = \sum_{j=1}^{N} \frac{(x_{j})_{3}}{(s_{N})_{3}},
\end{equation}
a 2-merger when three blocks is
\begin{equation}
\label{eq:5}
c_{N}(3;2) =  \sum_{j=1}^{N} \frac{3(x_{j})_{2}(s_{N} - x_{j})}{(s_{N})_{3}}.
\end{equation}



@<coalescence probs@>=@#
static void estimate_coalescence_probabilities( std::vector<double>& v_cN, const std::vector<double>& v_cdf,   std::vector<size_t>& v_pool_jvs)
{
@#
  size_t SN {} ;
 @#
  /* \newline \aths{|NUMBER_EXPS| \S~\ref{sec:constants}} */ 
  for( std::size_t i = 0 ; i < NUMBER_EXPS ; ++i){
 
    SN = sample_pool_juveniles( v_pool_jvs,  v_cdf) ;

@#

    for( size_t j = 0 ; j < CONST_POP_SIZE; ++j){

     /* \newline \aths{\eqref{eq:3}} */
      v_cN[0] += (static_cast<double>( v_pool_jvs[j])/static_cast<double>(SN)) * (static_cast<double>(v_pool_jvs[j] - 1) / static_cast<double>(SN - 1)) ;
      /* \newline \aths{\eqref{eq:4}} */
      v_cN[1] += (static_cast<double>( v_pool_jvs[j])/static_cast<double>(SN)) *  (static_cast<double>( v_pool_jvs[j]-1)/static_cast<double>(SN-1)) *  (static_cast<double>( v_pool_jvs[j]-2)/static_cast<double>(SN-2)) ;
      /* \newline \aths{\eqref{eq:5}} */
      v_cN[2] += 3.*(static_cast<double>(SN - v_pool_jvs[j]) / static_cast<double>(SN)) * (static_cast<double>( v_pool_jvs[j])/static_cast<double>(SN - 1)) *(static_cast<double>(v_pool_jvs[j] - 1) / static_cast<double>(SN-2)) ;
    } }
 
  v_cN[0] /= static_cast<double>(NUMBER_EXPS) ; @#
  v_cN[1] /= static_cast<double>(NUMBER_EXPS) ; @#
  v_cN[2] /= static_cast<double>(NUMBER_EXPS) ; @#
}



@*1 {\bf assign blocks to families}.
\label{sec:blockstofamilies}


assign blocks to families by sampling without replacement using  given
numbers of potential offspring 

@<blocks to families@>=@#
static void rmvhyper( std::vector<size_t>& merger_sizes,  size_t k, const std::vector<size_t>& v_juvs, const size_t SN, gsl_rng *r )
{

@#

  /* \newline \aths{ |k| is the current number of lines} */
  merger_sizes.clear();

@#

  size_t n_others =  SN - v_juvs[0] ;

@#

  /* \newline  \aths{ sample the number of blocks assigned to the first family}  */
  size_t x = gsl_ran_hypergeometric( r, v_juvs[0], n_others, k);

@#

  if( x > 1){
    /* \newline \aths{ only record  merger sizes} */
    merger_sizes.push_back(x ); }
    /* \newline \aths{ update the remaining number of blocks} */
  k -= x ;
  /* \newline \aths{ update new number of lines} */
  
  size_t i =0 ;
    /* \newline  \aths{ we can stop as soon as all lines 
       have been assigned to a family } */
  while( (k > 0) && (i < CONST_POP_SIZE-1) ){
  /* \newline \aths{ set the index to the one being sampled from} */
    ++i ;
    /* \newline \aths{ update |n_others|} */
    n_others -= v_juvs[i] ;
    x =  gsl_ran_hypergeometric( r,  v_juvs[i], n_others, k );
    if( x > 1){
      merger_sizes.push_back( x) ; }
      /* \newline \aths{ update the remaining number of blocks } */
      k -= x ;
  }
  /* \newline \aths{ check if at least two lines assigned to last individual } */
  if( k > 1){
    merger_sizes.push_back( k); }
}



@*1 {\bf merge blocks}.
\label{sec:mergeblocks}


@<block merging@>=@#
static void update_tree( std::vector<size_t>& tree, const std::vector<size_t>& merger_sizes )
{

@#

  std::vector<size_t> new_blocks {} ;

@#

  if( merger_sizes.size() > 0){
    /* \newline \aths{ at least one merger} */
    new_blocks.clear() ;
    /* \newline \aths{ shuffle the tree} */
    std::ranges::shuffle( tree,  rng );
    /* \newline \aths{ loop over the mergers} */
    for( const auto &m: merger_sizes){
      /* \newline \aths{ |m| is number of blocks merging;  |m| is at least two; append
      new block to vector of  new blocks} */
      assert( m > 1) ;
      /* \newline \aths{ record the size of the new block by summing the sizes of the
      merging blocks} */
      new_blocks.push_back( std::accumulate( std::rbegin( tree), std::rbegin(tree) + m, 0) ) ;
      assert( new_blocks.back() > 1) ;
      /* \newline \aths{ remove the rightmost |m| merged  blocks from tree} */
      tree.resize( tree.size() - m) ;
    }
    /*  \newline \aths{ append new blocks to tree} */
    tree.insert( tree.end(),  new_blocks.begin(),  new_blocks.end() ) ;
  }
  /* \newline \aths{ if no mergers then tree is unchanged} */
}




@*1 {\bf update $L_{i}^{N}(n)$}.
\label{sec:updatelengths}


update lengths $\ell_{i}$  given current block sizes


@<update lengths@>=@#
static void update_ebib( const std::vector< std::size_t>& tree,  std::vector<double>& vebib)
{  
   std::for_each(tree.begin(), tree.end(), [&vebib](const std::size_t &i){ ++vebib[0]; ++vebib[i]; });
}



@*1 {\bf update approximations of $\EE{R_{i}^{N}(n)}$}. 
\label{sec:updateri}

given lengths $\ell_{i}$ from one tree  update approximations $r_{i}$ of  $\EE{R_{i}^{N}(n)}$


@<update $r_{i}$@>=@#
static void update_estimate_ebib( const std::vector<double>& v_l,  std::vector<double>& v_ri)
{
  assert( v_l[0] > 0.); 
  const double d = v_l[0]; 

  std::transform( v_l.begin(), v_l.end(), v_ri.begin(), v_ri.begin(), [&d](const auto &x, const auto &y){return y + (x/d);});
}



@*1 {\bf three or two blocks left}.
\label{sec:threeortwo}


when at most three blocks left we use the coalescence probabilites
\eqref{eq:3}, \eqref{eq:4}, \eqref{eq:5} approximated in  \S~\ref{sec:estimatecoals} 


@<2 or 3 blocks@>=@#
static void three_or_two_blocks_left(  std::vector<double>& tmp_bib, const std::vector<double>& v_cN,  std::vector<size_t>& v_tree)
{

@#

  double Tk = 0.;
  double Tkk = 0. ;
  size_t newblock {} ;
  switch( v_tree.size() ){
  case 3 : {
    /* \newline \aths{ three lines left so  sample the two waiting times for a 3-merger and a 2-merger} */
    Tk = static_cast<double>( gsl_ran_geometric(rngtype, v_cN[1] ) ) ;
    Tkk = static_cast<double>( gsl_ran_geometric(rngtype, v_cN[2]));
    if( Tk < Tkk){
      /* \newline \aths{ all three blocks merge;  update the branch lengths }   */

      tmp_bib[0] += (3. * Tk) ;
      tmp_bib[ v_tree[0]] += Tk ;
      tmp_bib[ v_tree[1]] += Tk ;
      tmp_bib[ v_tree[2]] += Tk ;
      /* \newline  \aths{ clear the tree} */
      v_tree.clear() ;
      assert( v_tree.size() < 1);
    }
    else{
      /* \newline \aths{ a 2-merger occurs followed by a merger of the last two
      blocks} */
      tmp_bib[0] += (3. * Tkk) ;
      tmp_bib[ v_tree[0]] += Tkk ;
      tmp_bib[ v_tree[1]] += Tkk ;
      tmp_bib[ v_tree[2]] += Tkk ;
      /* \newline \aths{ shuffle  the tree } */
      std::ranges::shuffle( v_tree,  rng );
      newblock = v_tree[1] + v_tree[2] ;
      v_tree.resize(1) ;
      v_tree.push_back( newblock);
      assert( v_tree.size() == 2 );
      /* \newline \aths{ sample waiting time until merger of last two blocks} */
      Tk = static_cast<double>( gsl_ran_geometric(rngtype, v_cN[0] ) ) ;
      tmp_bib[0] += (2. * Tk) ;
      tmp_bib[ v_tree[0]] += Tk ;
      tmp_bib[ v_tree[1]] += Tk ;
      v_tree.clear() ;
      assert( v_tree.size() < 1); }
      break; }
        case 2 : {
    /* \newline \aths{ two blocks left} */
  Tk = static_cast<double>( gsl_ran_geometric(rngtype, v_cN[0] ) ) ;
    tmp_bib[0] += (2. * Tk) ;
    tmp_bib[ v_tree[0]] += Tk ;
    tmp_bib[ v_tree[1]] += Tk ;
    v_tree.clear() ;
    assert( v_tree.size() < 1);
    break ; }
  default : break ;
  }
}



@*1 {\bf approximate $\EE{R_{i}^{N}(n)}$}.
\label{sec:approximateri}


approximate $\EE{R_{i}^{N}(n)}$ given the settings in
\S~\ref{sec:constants}



@<go ahead --- approximate $\EE{R_{i}^{N}(n)}$@>=@#
static void estimate_ebib( const std::size_t & n_leaves )
{

@#
 /* \newline \aths{|n_leaves| is sample size $n$} */

@#
  
   std::vector< double > v_cdf (static_cast<std::size_t>(CONST_CUTOFF)
   + 1) ;

@#

  /* \newline  \aths{ compute the CDF function for sampling juveniles
  \S \ref{sec:cmf}  | generate_cdf( v_cdf);|} */
  

  std::vector<size_t> v_number_juvs (CONST_POP_SIZE) ;

@#

/* \newline \aths{ the tree; initially all blocks are singletons } */
  std::vector<size_t> v_tree (n_leaves, 1) ;
  std::vector<size_t> v_merger_sizes {};
  v_merger_sizes.reserve( n_leaves );

  std::vector<double> v_tmp_ebib (n_leaves, 0.0) ;
  std::vector<double> v_ebib (n_leaves, 0.0) ;
  std::vector<double> v_coal_probs (3, 0.0) ;

  /* \newline \aths{ estimate the coalescence probs
  \S \ref{sec:estimatecoals}} */
  estimate_coalescence_probabilities( v_coal_probs, v_cdf, v_number_juvs) ;
  
  size_t SN = 0;
  double number_experiments = CONST_NUMBER_EXPERIMENTS + 1.;
  while( --number_experiments  > 0.){
  /* \newline \aths{ initialise the tree as all singletons } */
  v_tree.clear();
  v_tree.assign( n_leaves, 1);
  /* \newline \aths{ initialise the container for the branch length for  the current
  realisation } */
  std::fill(  std::begin( v_tmp_ebib), std::end(  v_tmp_ebib ), 0.0 );
 

   while( v_tree.size() > 1){
    
    /* \newline \aths{update $\ell_{i}$ \S \ref{sec:updatelengths}  for the current tree configuration} */
    update_ebib( v_tree, v_tmp_ebib) ;
    if( v_tree.size() > 3 ){
    /* \newline \aths{ sample pool of juveniles \S~\ref{sec:pool}  } */
    SN = sample_pool_juveniles( v_number_juvs, v_cdf) ;
    /* \newline \aths{ compute the merger sizes \S \ref{sec:blockstofamilies}} */
    rmvhyper( v_merger_sizes, v_tree.size(),  v_number_juvs, SN, rngtype) ;
    /* \newline \aths{ update the tree \S \ref{sec:mergeblocks}} */
    update_tree( v_tree, v_merger_sizes);}
    else{
      /* \newline \aths{ at most three blocks left \S \ref{sec:threeortwo}} */
      three_or_two_blocks_left( v_tmp_ebib, v_coal_probs, v_tree) ;
    }
  }
    /* \newline \aths{ update approximation  of   $\EE{R_{i}^{N}(n)}$ \S \ref{sec:updateri} }  */
  update_estimate_ebib( v_tmp_ebib, v_ebib);
  }
  /* \newline \aths{ print approximation  of $\EE{R_{i}^{N}}$ } */
  for( const auto&r: v_ebib){
    std::cout << r << '\n' ; }
}


@*1 {\bf main}.
\label{sec:main}

the |main| module

@C

/* \newline \aths{\S~\ref{sec:includes}} */
@<includes@>@#
/* \newline \aths{\S~\ref{sec:constants}} */
@<constants@>@#
/* \newline \aths{\S~\ref{sec:rngs}} */
@<the rngs@>@#
/* \newline \aths{\S~\ref{sec:pmf}} */
@<compute \eqref{eq:1}@>@#
/* \newline \aths{\S~\ref{sec:cmf}} */
@<cmf@>@#
/* \newline \aths{\S~\ref{sec:sampleoneX}} */
@<sample an $X$@>@#
/* \newline \aths{\S~\ref{sec:unboundedX}} */
@<$X$ when unbounded@>@#
/* \newline \aths{\S~\ref{sec:pool}} */
@<pool@>@#
/* \newline \aths{\S~\ref{sec:estimatecoals}} */
@<coalescence probs@>@#
/* \newline \aths{\S~\ref{sec:blockstofamilies}} */
@<blocks to families@>@#
/* \newline \aths{\S~\ref{sec:mergeblocks}} */
@<block merging@>@#
/* \newline \aths{\S~\ref{sec:updatelengths}} */
@<update lengths@>@#
/* \newline \aths{\S~\ref{sec:updateri}} */
@<update $r_{i}$@>@#
/* \newline \aths{\S~\ref{sec:threeortwo}} */
@<2 or 3 blocks@>@#
/* \newline \aths{\S~\ref{sec:approximateri}} */
@<go ahead --- approximate $\EE{R_{i}^{N}(n)}$@>@#




int main(int argc, char *argv[])
{

@#


/* \newline \aths{\S \ref{sec:rngs}} */
setup_rng(  static_cast<unsigned long int>(atoi(argv[2])) );


/* \newline \aths{\S \ref{sec:approximateri} } */
 estimate_ebib( static_cast<std::size_t>(atoi(argv[1]))) ;

@#

gsl_rng_free( rngtype ); 

return 0 ;
}


@* {\bf conclusions and references}.
\label{sec:concl}


The idea behind approximating  $\EE{R_{i}^{N}(n)}$ is to  compare
$\EE{R_{i}^{N}(n)}$ to 
$\EE{R_{i}(n)}$, the quantities corresponding to  $\EE{R_{i}^{N}(n)}$
and predicted by $\set{\xi^{n}}$,   the limit of  $\set{\xi^{n,N}}$
as $N\to \infty$.  For example, for $1\le \alpha < 2$ the
(time-rescaled)  limit (in the sense of convergence of
finite-dimensional distributions)   of     $\set{\xi^{n,N}}$ is the
Beta$(2-\alpha,\alpha)$ coalescent \cite{schweinsberg03}. Here we
provide the code for  approximating   $\EE{R_{i}^{N}(n)}$ when the
sample comes from a haploid panmictic population of constant size. The
population may be evolving according to sweepstakes reproduction
(heavy-tailed offspring number distribution; \eqref{eq:1}).    





\begin{SCfigure}[0.8][htb]
    \begin{tikzpicture}
      \begin{axis}[
        xlabel = $\log(i/n) - \log(1 - i/n)$,
        axis line style = {draw = none},
        tick style = {draw=none},
        xticklabels={draw=none},
        yticklabels={draw=none}]
        \addplot table[col sep=comma] {forplottingfile1};
        \addplot table[col sep=comma] {forplottingfile2};
       \end{axis}
       \end{tikzpicture}
       \caption{\it  An example approximation of $\EE{R_{i}^{N}(n)}$
graphed as logits   given the
settings in \S~\ref{sec:constants}; see also  \S~\ref{sec:compscript}  }
       \label{fig:graph}
       \end{SCfigure}



We use simulations to approximate $\EE{R_{i}^{N}(n)}$; an example is
given in Figure~\ref{fig:graph}. 
One  may (for example) investigate the effect of increasing $n$, or
the upper bound $\zeta(N)$, on 
$\EE{R_{i}^{N}(n)}$.



\bibliographystyle{plain}
%%\bibliography{refs}


\begin{thebibliography}{99}

\bibitem[D2024]{D2024}
 Diamantidis,  Dimitrios and Fan,  Wai-Tong (Louis) and Birkner,
 Matthias and Wakeley,  John.  {Bursts of coalescence within
 population pedigrees whenever big families occur}.  Genetics Volume 227,  February
  2024.
 \\
 \url{https://dx.doi.org/10.1093/genetics/iyae030}.


\bibitem[CDEH25]{chetwyn-diggle_beta}
JA~Chetwyn-Diggle, Bjarki Eldon, and Matthias Hammer.
\newblock Beta-coalescents when sample size is large.
\newblock In preparation, 2025+.

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



\bibitem[Tan11]{tange11:_gnu_paral}
O~Tange.
\newblock {GNU} parallel -- the command-line power tool.
\newblock The USENIX Magazine, 2011.

\end{thebibliography}




@
\end{document}
