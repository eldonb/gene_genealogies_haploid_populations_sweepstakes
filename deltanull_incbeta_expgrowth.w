%$ NAFN=deltanull_incbeta_expgrowth
%$ echo $NAFN
%$ echo "const double ALPHA = 1.5;" > $NAFN.hpp
%$ echo "const double KAPPA = 2;" >> $NAFN.hpp
%$ echo "const double CEPS = 1;"  >> $NAFN.hpp
%$ echo "const double GAMMA = 0.1;"  >> $NAFN.hpp
%$ echo "const double RHO =  25. ;"  >> $NAFN.hpp
%$ echo "const double MINF = (2 + (1 + pow(2,1-KAPPA))/(KAPPA  -1))/2;"  >> $NAFN.hpp
%$ echo "const double DK = KAPPA > 2 ? (pow(2,KAPPA+1)*(KAPPA -2)*(KAPPA - 1)) : 1;" >> $NAFN.hpp
%$ echo "const double CK = (KAPPA > 2 ? ((KAPPA + 2) + pow(KAPPA,2))/DK : 1);" >> $NAFN.hpp
%$ echo "const double CKAPPA = 2*CK/pow(MINF,2);" >> $NAFN.hpp
%$ echo "const double BETA = (GAMMA < 1 ? gsl_sf_beta_inc(2-ALPHA,ALPHA,GAMMA) : 1.)*gsl_sf_beta(2-ALPHA,ALPHA);" >> $NAFN.hpp
%$ echo "const double CKAG = CKAPPA + (CEPS*ALPHA*BETA/pow(MINF,ALPHA));" >> $NAFN.hpp
%$ echo "const double SAMPLE_SIZE = 100;" >> $NAFN.hpp
%$ echo "const int EXPERIMENTS = 2000;" >> $NAFN.hpp
%$ ctangle $NAFN
%$ echo $NAFN
%$ g++ -std=c++26 -m64 -march=native -O3 -x c++ -DNDEBUG $NAFN.c -lm -lgsl -lgslcblas
%$ rm -f tmp_runs
%$ for i in $(seq 5); do echo "./a.out  "  $(shuf -i 39393-282929191 -n1) "> resout"$i >> tmp_runs; done
%$ parallel --gnu -j5 :::: ./tmp_runs
%$ paste resout* -d, | sed '1d' | awk -F, '{M=1e4;s=0; for (i=1;i<=NF;i++) {s+=$i} print log(s/M) - log(1-(s/M))}' > logitresout
%$ seq 99 | awk '{n=100;print log($1/n) - log(1 - ($1/n))}' > nlogits
%$ paste -d',' nlogits logitresout > forplottingfile1
%$ NAFN=deltanull_incbeta_expgrowth
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
$\EE{R_{i}(n)}$ for the time-changed $\delta_{0}$-Beta$(\gamma,2-\alpha,\alpha)$-coalescent}
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
time-changed $\delta_{0}$-Beta$(\gamma,2-\alpha,\alpha)$-coalescent
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
sample the time-changed $\delta_{0}$-Beta$(\gamma,2-\alpha,\alpha)$-coalescent with
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
\label{sec:comp}

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





Here we consider the continuous-time
$\delta_{0}$-Beta$(\gamma,2-\alpha,\alpha)$-coalescent with the
time-change corresponding to exponential population growth where the
time-change function $G(t) = \int_{0}^{t}e^{\rho s}ds $ for some fixed
$\rho \ge 0$.


One samples  the
$\delta_{0}$-Beta$(\gamma,2-\alpha,\alpha)$-coalescent by
sampling group sizes $2 \le k_{1},\ldots,k_{r} \le n$ with
$\sum_{i}k_{i}\le n$, $s = n - \sum_{i}k_{i}$,    according to  
\begin{equation}
\label{eq:1}
\begin{split}
\lambda_{n;k;s} & = \frac{1}{C_{\kappa,\alpha,\gamma}}\binom{n}{k}\svigi{\one{r=1,k_{1}=2}C_{\kappa} +  \frac{c \alpha }{\mathbbm m^{\alpha}}B(\gamma,k-\alpha,n-k + \alpha)}\\
C_{\kappa} &  =  \frac{2}{\mathbbm m^{2}}\svigi{\one{\kappa = 2} +  \one{\kappa > 2} \frac{c_{\kappa}}{2^{\kappa}(\kappa - 2)(\kappa - 1) } } \\
C_{\kappa,\alpha,\gamma} & =    C_{\kappa} +  c \alpha \mathbbm m^{-\alpha}B(\gamma,2-\alpha,\alpha) \\
\mathbbm m & =  1 + \frac{1 +  2^{1-\kappa} }{2(\kappa - 1)}
\end{split}
\end{equation}
where $2 + \kappa < c_{\kappa} < \kappa^{2}$ when $\kappa > 2$, and
$B(p,a,b) \equiv \int_{0}^{1}\one{0<t\le p}t^{a-1}(1-t)^{b-1}dt$ for
$a,b>0$ and $0< p \le 1$ fixed.




Let $\# A$ denote the number of elements in a finite set $A$.  For a
given coalescent $\set{\xi^{n}}$ write
$L_{i}(n) \equiv \int_{0}^{\tau(n)} \# \set{ \xi \in \xi^{n}(t) :
\#\xi = i }dt$ and $L(n) \equiv \int_{0}^{\tau(n)} \# \xi^{n}(t)dt $
where $\tau(n) \equiv \inf \set{t \ge 0 : \#\xi^{n}(t) = 1}$. Then
$L(n) = L_{1}(n) + \cdots + L_{n-1}(n)$. Write
$R_{i}(n) \equiv L_{i}(n)/L(n)$ for $i = 1,2,\ldots, n-1$.  With this
C++ code we use 
simulations to approximate $\EE{R_{i}(n)}$ when the coalescent is the
time-changed continuous-time
$\delta_{0}$-Beta$(\gamma,2-\alpha,\alpha)$-coalescent with
time-change function $G(t) = \int_{0}^{t}e^{\rho s}ds$ for a given
fixed $\rho \ge 0$.  See Figure~\ref{fig:graph} in \S~\ref{sec:concl}
for an example


The code follows in \S~\ref{sec:includes}--\S~\ref{sec:main}, we
conclude in \S~\ref{sec:concl}. Comments within the code in \aths{this
font and colour}



@* {\bf code}.
\label{sec:code}


@*1 {\bf includes}. 
\label{sec:includes}


the included libraries

@<includes@>=@#
#include <iostream>
#include <fstream>
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
#include <chrono>
#include <forward_list>
#include <assert.h>
#include <math.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>
#include <boost/math/special_functions/beta.hpp>
#include "deltanull_incbeta_expgrowth.hpp"


@*1 {\bf random number generators}.
\label{sec:rgs}


the random number generators

@<rngs@>=@#
  std::random_device randomseed;
  /* \newline \aths{ Standard Mersenne twister  random number engine} */
  std::mt19937_64 rng(randomseed());


gsl_rng * rngtype ;
void setup_rng( const  unsigned long int s )
{
        const gsl_rng_type *T ; 
        gsl_rng_env_setup(); 
        T = gsl_rng_default ;
        rngtype = gsl_rng_alloc(T);
        gsl_rng_set( rngtype,  s) ;
}



@*1 {\bf the merging rate}. 
\label{sec:rate}


@<$\lambda_{n;k;s}$@>=@#
double rate(  const   double m, const   double k)
{
@#  
  assert( k <= m);

@#
/* \newline \aths{the full (non-normalised) incomplete beta function}*/
    auto incbeta = [](const double a, const double b, const double x)
    {
    return  boost::math::beta( static_cast<long double>(a), static_cast<long double>(b), static_cast<long double>(x));
    };
  @#
   const long double lbinom =   lgammal( m + 1) - lgammal(k+1) - lgammal(m-k+1) ;
@#
 const long double dnull =  (int)k < 3 ? logl((long double)(CKAPPA)) : 0 ;
@#
 const long double bpart =  logl((long double)(CEPS * ALPHA * incbeta(k-ALPHA, m + ALPHA -k, GAMMA)/pow(MINF ,ALPHA)));
@#
 return (static_cast<double>( expl( lbinom + dnull + bpart) ) / CKAG) ;
}



@  {\bf totalrate}

compute the total rate

@<$\lambda_{n}$@>=@#
void totalrate( std::vector<double>& v )
{
@#
  for( double m = 2 ; m <= SAMPLE_SIZE ; ++m){
  @#
    for ( double j = 2 ; j <= m ; ++j){
    @#
      assert( j <= m);
      /* \newline \aths{\S~\ref{sec:rate}} */
      v[m] += rate( m, j); }}
}


@*1 {\bf one experiment}.
\label{sec:genealogy}


@<genealogy@>=@#
static void genealogy( const std::vector<double>& vlambdan,  std::vector<double>& vri )
{
@#
  std::vector<int> t (SAMPLE_SIZE, 1);
  std::size_t ms {} ;
  double timi {} ;
  int newb {} ; 
  std::vector<double> vb( SAMPLE_SIZE) ;
  std::size_t q {} ;
  double otimi {} ;

  auto newtime = [](const double lambdab, const double oldtime)
  {
      return (RHO > DBL_EPSILON ? log1p(  -(RHO * exp(- RHO * oldtime) / lambdab) * log(gsl_rng_uniform_pos(rngtype)))/RHO : gsl_ran_exponential(rngtype, 1./lambdab)) ;
  } ;


@#

  auto getmerger = [&vlambdan](const double m)
  {@#
  const double u = gsl_rng_uniform( rngtype);
  @#
  double j = 2;
  /* \newline \aths{\S~\ref{sec:rate}} */
  double s =  rate( m, j) ;
@#
  while( u > s/vlambdan[ static_cast<int>(m) ]){  ++j; assert( j <= m); s += rate( m, j);  }
@#
  return static_cast<std::size_t>(j); };

@#

  while(t.size() > 1){
  @#
    timi = newtime( vlambdan[t.size()], otimi);
    assert(timi > DBL_EPSILON);
    otimi += timi ;
   @#
    std::for_each( t.cbegin(), t.cend(), [&timi, &vb]( const int t ){ assert(t > 0);  vb[0] += timi; vb[t] += timi; } );
    ms =  getmerger(static_cast<double>(t.size()));
    std::shuffle( t.begin(), t.end(), rng);
    newb = std::accumulate( t.rbegin(), t.rbegin() + ms, 0 );
    q = t.size() ;
    t.resize( q - ms );
    t.push_back( newb); }
@#
  assert( vb[0] > DBL_EPSILON);
  @#
  const double d = vb[0] ;
  @#
  std::transform( vb.cbegin(), vb.cend(), vri.begin(), vri.begin(), [&d]( const auto &x, const auto &y ){return y + (static_cast<double>(x)/d); });
}



@  {\bf approximate $\EE{R_{i}(n)}$}. 


@<go ahead@>=@#
void approximate()
{
@#
  std::vector<double> v_lambdan ( static_cast<int>(SAMPLE_SIZE) + 1 );
  @#
  totalrate( v_lambdan);
@#
  std::vector<double> vri (static_cast<int>(SAMPLE_SIZE) );
  int r = EXPERIMENTS + 1;
  while( --r > 0){
  /* \newline \aths{\S~\ref{sec:genealogy}} */
    genealogy( v_lambdan, vri); }
@#
  std::for_each( vri.begin(), vri.end(), [](const auto &x){ std::cout << x << '\n'; });
}






@*1 {\bf main}.
\label{sec:main}

the |main| module

@C

/* \newline \aths{\S~\ref{sec:includes}} */
@<includes@>@#
/* \newline \aths{\S~\ref{sec:rgs}} */
@<rngs@>@#
/* \newline \aths{\S~\ref{sec:rate}} */
@<$\lambda_{n;k;s}$@>@#
@<$\lambda_{n}$@>@#
/* \newline \aths{\S~\ref{sec:genealogy}} */
@<genealogy@>@#
@<go ahead@>@#


int main(int argc, const char *argv[])
{

/* \newline \aths{\S~\ref{sec:rgs}} */
  setup_rng( static_cast<long unsigned>( atoi(argv[1])));

/* \newline \aths{\S~\ref{sec:genealogy}} */
   approximate( ) ;

@#
  gsl_rng_free( rngtype);

@#
  return GSL_SUCCESS; 
}


@* {\bf conclusions and bibliography}. 
\label{sec:concl}






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
JA~Chetwyn-Diggle and  Bjarki Eldon.
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
