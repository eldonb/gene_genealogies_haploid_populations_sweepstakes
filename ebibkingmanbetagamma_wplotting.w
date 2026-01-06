%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Copyright (C) 2025 Bjarki Eldon
%$ NAFN=ebibkingmanbetagamma_wplotting
%$ ctangle $NAFN
%$ gcc -std=c23 -m64 -march=native -O3 $NAFN.c -lm -lgsl -lgslcblas
%$ ./a.out > forplottingfile1
%$ gcc --version | head -n1 > innleggcpp
%$ gsl-config --version > innleggsl
%$ xelatex --version | head -n1  > innleggxelatex
%$ lualatex --version | head -n2 > innlegglualatex
%$ cweave  --version | head -n1 > innleggcweave
%$ uname  --kernel-release -o  > innleggop
%$ ctangle  --version | head -n1  > innleggctangle
%$ emacs --version | head -n1 > innleggemacs
%$ awk --version | head -n1 > innleggawk
%$ cweave $NAFN
%$ tail -n4 $NAFN.tex > endi
%$ for i in $(seq 5); do $(sed -i '$d' $NAFN.tex) ; done
%$ cat endi >> $NAFN.tex
%$ lualatex $NAFN
%$ lualatex $NAFN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\documentclass[a4paper,10pt]{cweb}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
%%\usepackage[lf]{Baskervaldx}
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
\hypersetup{
    colorlinks,
    linkcolor={gray!50!black},
    citecolor={blue!50!black},
    urlcolor={blue!80!black}
}
\usepackage{enumerate}
\usepackage{dsfont}
\usepackage[right]{lineno}
\usepackage{verbatim}
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
\setstretch{1.}
\newcommand{\uN}{\ensuremath{\zeta(N)}}
\newcommand{\at}{\ensuremath{\kappa}}
\newcommand{\one}[1]{\ensuremath{\mathds{1}_{\left\{ #1 \right\}}}}
\newcommand{\EE}[1]{\ensuremath{\mathds{E}\left[ #1 \right]}}
\newcommand{\im}{\ensuremath{\imath} }
\newcommand{\jm}{\ensuremath{\jmath} }
\newcommand{\set}[1]{\ensuremath{\left\{ #1 \right\}}}
\newcommand{\prob}[1]{\ensuremath{\mathds{P}\left( #1 \right) } }%
\newcommand{\h}[1]{\ensuremath{\uptheta_{ #1 } } }%
\newcommand{\VV}[1]{\ensuremath{ \mathbb{V}\left( #1 \right)}}%
\newcommand{\hp}{\ensuremath{\theta_1}}
\newcommand{\hs}{\ensuremath{\theta_2}}
\newcommand{\D}{\ensuremath{\mathbb{D}}}
\newcommand{\F}{\ensuremath{\mathbb{F}} }
\newcommand{\G}{\ensuremath{\mathbb{G}} }
\newcommand{\bt}[1]{\textcolor{blue}{\tt #1}}
\newcommand{\aths}[1]{\textcolor{violet}{\bf #1}}
%%
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
%%
\makeatletter
\renewcommand{\maketitle}{\bgroup\setlength{\parindent}{0pt}
\begin{flushleft}
  \textsf{\textbf{\@@title}}
\end{flushleft}
\begin{center}
\textsc{\@@author}
\end{center}
\egroup
}
\makeatother
\title{Gene genealogies in  haploid populations evolving according to  sweepstakes reproduction \\  --- exact $\EE{L_{i}(n)}$ for the  $\delta_{0}$-Beta($\gamma, 2-\alpha, \alpha$)-coalescent }
\author{Bjarki Eldon\footnote{\href{mailto:beldon11@@gmail.com}{beldon11@@gmail.com}}
\footnote{compiled @@ {\DTMcurrenttime} on  {\today}  \\ \input{innleggop}\\   \input{innleggctangle}  \\  \input{innleggcpp} \\ GSL \input{innleggsl}\\  \input{innleggcweave} \\  \input{innlegglualatex}   \\ \input{innleggawk}  \\ \input{innleggspix} \\   \input{innleggemacs}}
 \orcidlink{0000-0001-9354-2391}}
\begin{document}
\maketitle
\renewcommand{\abstractname}{\vspace{-\baselineskip}}

%%\rule{\textwidth}{.8pt}

\begin{abstract}
Let $\set{\xi^{n}} \equiv \set{\xi^{n}(t); t \ge 0}$ denote the
$\delta_{0}$-Beta($\gamma, 2-\alpha, \alpha$)-coalescent for
$0 < \gamma \le 1$ and $0 < \alpha < 2$. Write
$L_{i}(n) \equiv \int_{0}^{\tau(n)} \#\set{\xi \in \xi^{n}(t) : \#\xi
= i } dt $ for $i = 1, \ldots, n-1$,
$L(n) \equiv L_{1}(n) + \cdots + L_{n-1}(n)$, and
$\tau(n) \equiv \inf\set{ t \ge 0 : \#\xi^{n}(t) = 1} $.  With this C code we compute
$\EE{L_{i}(n)}$ predicted by the
$\delta_{0}$-Beta($\gamma, 2-\alpha, \alpha$)-coalescent.  This family
of multiple-merger coalescents can be obtained from a model of a
haploid panmictic population of constant size evolving according to
sweepstakes reproduction (heavy right-tailed offspring number
distribution).
\end{abstract}

\tableofcontents


@* {\bf Copyright}. 


Copyright {\textcopyright} {\the\year}  Bjarki Eldon \newline



This document and any source code it contains is distributed under the
terms of the GNU General Public Licence (version $\ge 3$).  You should
have received a copy of the licence along with this file (see file
COPYING).


    The source codes described in this document are free software: you
    can redistribute it and/or modify it under the terms of the GNU
    General Public License as published by the Free Software
    Foundation, either version 3 of the License, or (at your option)
    any later version.

    This document and the code it contains is distributed in the hope
    that it will be useful, but WITHOUT ANY WARRANTY; without even the
    implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
    PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this file (see COPYING).  If not,
    see \url{http://www.gnu.org/licenses/}.


@* {\bf Compilation,  output and execution}. 
\label{compile}

 This CWEB
      \citep{knuth1994cweb} document (the {\tt .w} file) can be
      compiled with {\tt cweave} to generate a {\tt .tex} file, and
      with {\tt ctangle} to generate a {\tt .c} \citep{kernighan1988c}
      file.

One can use {\tt cweave} to generate a {\tt .tex} file, and {\tt
ctangle} to generate a {\tt .c} file. To compile the C++ code (the {\tt
.c} file), one needs the GNU Scientific Library. 
Using a Makefile can be helpful, naming this file {\tt iguana.w}


 {\tt
iguana.pdf : iguana.tex \\
\tab\quad\quad\quad\quad cweave iguana.w \\
\tab\quad\quad\quad\quad        pdflatex iguana \\
\tab\quad\quad\quad\quad        bibtex iguana \\
\tab\quad\quad\quad\quad        pdflatex iguana \\
\tab\quad\quad\quad\quad        pdflatex iguana \\
\tab\quad\quad\quad\quad        ctangle iguana \\
\tab\quad\quad\quad\quad        gcc -std=c23 -O3 -march=native -m64 -x c iguana.c -lm -lgsl -lgslcblas \\
        
       
clean :  \\
\tab\quad\quad\quad\quad        rm -vf iguana.c iguana.tex \\
}


Use {\tt valgrind} to check for memory leaks:

{\tt valgrind -v ---leak-check=full ---show-leak-kinds=all --vgdb=full --leak-resolution=high --num-callers=50 <program call>}

Use {\tt cppcheck} to check the code:

{\tt  cppcheck ---enable=all ---language=c  iguana.c}


Use {\tt splint} to check the code:



@* {\bf intro}.  
\label{sec:intro}

Suppose the population evolves as in Definition~\ref{hschwpop}.  

\begin{defn}[Evolution of a haploid panmictic population]
\label{hschwpop}%
Consider a haploid, i.e.\ each individual carries exactly  one gene
copy, panmictic population of constant size $N$ evolving in discrete
(non-overlapping) generations.  In any given generation each
individual independently produces a random number of potential
offspring according to some given law.  If the total number of
potential offspring produced in this way is at least $N$, then $N$ of
them sampled uniformly at random without replacement form a new set of
reproducing individuals, and the remaining potential offspring perish.
Otherwise we will assume an unchanged population over the generation
(all the potential offspring perish).
\end{defn}%


Let $\{\pi_{N,n}(t); t \ge 0\}$ denote the Markov sequence describing
the random ancestral relations of $n$ gene copies (leaves) sampled
from a population of size $N$ evolving according to  Definition~\ref{hschwpop}.   
If  the random number of potential offspring ($X$)   of an arbitrary
individual is  distributed according to 
\begin{equation}
\label{eq:4}
\lim_{x \to \infty}Cx^{\alpha}\prob{X \ge x}  = 1
\end{equation}
\cite[Eq~11]{schweinsberg03} then
$\{\pi_{N,n}(\lfloor t/c_{N} \rfloor ); t \ge 0\}$ converges as
$N\to \infty$ in the $J_{1}$-Skorohod topology to the
Beta$(2-\alpha,\alpha)$-coalescent provided $1 \le \alpha < 2$
\cite[Thm 4c]{schweinsberg03}.


Suppose $X$ is distributed according to 
\begin{equation}
\label{eq:PXiJ}
g(k) \left( \frac{1}{k^\alpha} - \frac{1}{(1+k)^\alpha} \right)  \leq
\prob{X = k} \leq f(k) \left( \frac{1}{k^\alpha} -
\frac{1}{(1+k)^\alpha} \right), \quad k  \in
[\uN],
\end{equation}
where the functions $g,f$ are independent of $N$ and so that  $\EE{X} > 1$ and  $\prob{X \le \uN} = 1$.   
Furthermore, the population evolves according to Definition~\ref{def:haplrandomalpha} or Definition~\ref{def:alpha-random-one}.
\begin{defn}[The random environment of type $A$]
\label{def:haplrandomalpha}
%%  by Eq~\eqref{eq:PXiJ}.
Suppose a population evolves according to Definition~\ref{hschwpop}
with  the number of potential offspring  produced by each individual distributed
according to   \eqref{eq:4} or   \eqref{eq:PXiJ}. 
Fix  $0 < \alpha < 2 \le \at $.   Write $E $ for
the event  $\{ X_i \sim L(\alpha, \uN) \text{ for all }
 i \in [N] \}$ and $E^{\sf c}$  when $\at$ replaces $\alpha$ in $ E$   for some given $\uN$.  
Suppose, with  $(\varepsilon_{N})_{N}$ a positive sequence with $0 < \varepsilon_{N} < 1$ and $\varepsilon_{N}\to 0$, 
\begin{displaymath}
\prob{ E}  = \varepsilon_N, \quad \prob{ E^{\sf c}}  = 1 - \varepsilon_N
\end{displaymath}
%% Further, it will be assumed that Assumption~\ref{a:mN1plus} holds so that $\EE{X_{1}} > 1$ for all $N$. 
\end{defn}%


\begin{defn}[The random environment of type $B$]
\label{def:alpha-random-one}
Suppose a population evolves according to Definition~\ref{hschwpop} with
the law for the number of potential offspring given by
Eq~\eqref{eq:PXiJ}.  Fix $0 < \alpha < 2 \le \at$.  Write $ E_{1}$ for
the event 
$\{\text{there exists exactly one } i \in [N] : X_i \sim L(\alpha, \uN), X_{j} \sim L(\at,
\uN) \text{ for all } j \in [N]\setminus \{i\} \}$ with $i$ picked
uniformly at random, and $ E _{1}^{\sf c}$ when $\at$ replaces
$\alpha$ in $ E_{1}$ so that $ E_{1}^{\sf c} = E^{\sf c}$ from
Definition~\ref{def:haplrandomalpha}. Suppose
%%write $E_{2}$ for the event  $\{ X_i^{N} \sim L(\at, \uN)\, \forall\, i \in [N] \}$.
\begin{displaymath}
\begin{split}
\prob{ E_{1} } & = \overline \varepsilon_N, \quad \prob{ E_{1}^{\sf c}}  = 1 - \overline \varepsilon_N
\end{split}
\end{displaymath}
where  $(\overline \varepsilon_{N})_{N}$ is  a positive sequence with $0 < \overline \varepsilon_{N} < 1$ for all $N$. 
\end{defn}
Then $\{\pi_{N,n}(\lfloor t/c_{N} \rfloor ); t \ge 0\}$ converges as
$N\to \infty$ in the $J_{1}$-Skorohod topology to the
Kingman-Beta$(\gamma,2-\alpha, \alpha)$-coalescent provided $\liminf_{N\to \infty} \uN/N > 0$.   



Let $L_{i}(n)$ denote the random length of branches supporting
$i \in \{1,2,\ldots, n-1\}$ leaves; write $L(n) =
\sum_{i}L_{i}(n)$. We compute
\begin{equation}
\label{eq:varphi}
\varphi_{i}(n) :=  \frac{ \EE{L_{i}(n)}}{ \EE{L(n)}} 
\end{equation}
by adapting a recursion for computing $\EE{L_{i}(n)}$ for
$\Lambda$-coalescents \cite{BBE2013a}.  We compute $\varphi_{i}(n) $
for the Kingman-Beta($\gamma, 2-\alpha, \alpha$)-coalescent with
transition rates
\begin{subequations}
\begin{align}
\label{eq:ratesarandall}
\lambda_{n,k} & =  \one{k=2}\binom{n}{k}\frac{ C_{\kappa} }{C} +   \binom{n}{k} \frac{ \alpha c_{\alpha}^{\prime}  }{ Cm_{\infty}^{\alpha} } B(\gamma, k-\alpha,n - k+ \alpha)  \\
\label{eq:27}
%% \mathfrak{w}
\gamma & =  \one{\frac{\uN}{N} \to K } \frac{K}{m_{\infty} + K} + \one{ \frac{\uN}{N} \to \infty }  \\
\label{eq:23}
C &  =   C_{\kappa}   +  \frac{ \alpha c }{m_{\infty}^{\alpha}} \int_{\one{\uN/N \to K}\frac{m_{\infty} }{K + m_{\infty} }  }^{1} (1-y)^{1-\alpha}y^{\alpha - 1}dy  \\
\label{eq:1}
C_{\kappa} & =   2 m_{\infty}^{-2}   \one{\kappa = 2}    +  2 m_{\infty}^{-2} \frac{c_\kappa}{2^{\kappa}(\kappa - 2)(\kappa - 1)  } \one{\kappa > 2 }
\end{align}
\end{subequations}
The range of $\alpha$ in Definition~\ref{def:haplrandomalpha} and
Definition~\ref{def:alpha-random-one} intersects at 1 so one must decide
which of the two definitions holds when $\alpha = 1$.  Let where
$0 < c^{\prime} < 1$ and $c > 0$. Write $c_{\alpha}^{\prime} = c$ if
Definition~\ref{def:haplrandomalpha} holds when $\alpha = 1$, and
$c_{\alpha}^{\prime} = \one{\alpha = 1}c^{\prime} + \one{\alpha \in
(0,2)\setminus \{1\} }c $ when Definition~\ref{def:alpha-random-one} holds
when $\alpha = 1$ (and $\kappa > 2$).  We approximate $m_{\infty}$
with $m_{\infty} \approx (2 + (1 + 2^{1 - \kappa })/(\kappa - 1))/2$.  When $\kappa > 2$  
\begin{displaymath}
%%\frac{\kappa 2^{2-\kappa} }{\kappa - 2 } - \frac{3\kappa 2^{1-\kappa} }{\kappa - 1} + 2^{-\kappa} \le c_{\kappa} \le  \frac{\kappa 2^{2-\kappa}  }{\kappa - 2} - \frac{\kappa 2^{1-\kappa} }{\kappa - 1}
\kappa + 2 < c_{\kappa} < \kappa^{2}
\end{displaymath}

@* {\bf code}. 
\label{sec:code}


The code requires the GSL Library.  Comments within the code are in
\aths{this font and colour}; the code follows in
\S~\ref{sec:includes}--\S~\ref{sec:main}, we conclude in
\S~\ref{sec:concl}.

@*1 {\bf the included libraries}.
\label{sec:includes}


the included C libraries 

@<includes@>=@#
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_rng.h> 
#include <gsl/gsl_randist.h> 
#include <gsl/gsl_vector.h> 
#include <gsl/gsl_matrix.h> 
#include <gsl/gsl_sf.h> 
#include <gsl/gsl_errno.h> 
#include <gsl/gsl_fit.h> 
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_sort.h> 
#include <gsl/gsl_statistics_double.h> 
#include <gsl/gsl_integration.h> 
#include <gsl/gsl_errno.h> 


@*1 {\bf the beta function}. 
\label{sec:betafunc}

 return the logarithm of  the (incomplete) beta function  $B(x,a,b) =  \int_{0}^{x} t^{a-1}(1-t)^{b-1}dt$ for $0 < x \le 1$ and $a,b > 0$ using the Gauss hypergeometric function
 $B(x,a,b) =  a^{-1}x^{a}(1-x)^{b}F(a+b, 1, a+1, x)$


@<beta function@>=@#
static long double betafunc( const double x, const double a, const  double b )
{
  /* \newline  \aths{the GSL incomplete beta function is normalised by the complete beta function}  */
  assert( x > 0);
  assert( a > 0);
  assert( b > 0);
  /* \newline \aths{ the standard way would be | gsl_sf_beta_inc( a, b, x) * gsl_sf_beta( a, b) |} */
  /* \newline \aths{ return  the logarithm  of the beta function  as $\log\Gamma(a) + \log\Gamma(b) - \log\Gamma(a+b)$ } */
  const long  double f = (long double)(x < 1 ? gsl_sf_hyperg_2F1(a + b, 1, a+1,x) : 1) ;
  
  assert( f > 0);
  /* \newline  \aths{ return $\one{x < 1}( \log f  +a\log x + b\log(1-x) - \log a) + \one{x = 1} (\log\Gamma(a) + \log\Gamma(b) - \log\Gamma(a+b))$} */
  return( x < 1 ?  (logl( f ) +  ((long double)( (a*log(x)) + (b * log(1-x)) - log(a) ) ) ) : (lgammal( (long double)a) + lgammal( (long double)b) - lgammal( (long double)(a + b) ) ) ); 
}


@*1 {\bf binomial constant}.
\label{SEC:binomialconstant}

compute the log of the  binomial constant $\binom{n}{k}$  

@<binomial@>=@#
static long  double binom( const int n, const int k)
{


assert( k <= n);
assert( k >= 0);
assert( n >= 0);
/* \newline\aths{  $\log \binom{n}{k} =  \log \Gamma(n+1) - \log\Gamma(k+1) - \log\Gamma(n-k + 1)$} */
   return( lgammal( (long double)(n+1)) - lgammal( (long double)(k+1)) - lgammal( (long double)(n-k + 1)) ) ; 
}


@*1 {\bf total beta rate}. 
\label{sec:totalbetarate}

return the total beta rate  $ \binom{n}{k}B(\gamma, k-\alpha, n+\alpha - k)$. 

@<betarate@>=@#
static double lbetank( const int n, const int k, const double a, const double g)
{
/* \newline\aths{ return $\exp( \log \binom{n}{k} +  \log B(\gamma, k-\alpha, n-k + \alpha) )$} */
/* \newline\aths{ \S~\ref{SEC:binomialconstant} and \S~\ref{sec:betafunc}} */
  return(  (double) expl( binom(n, k)  +   betafunc( g, ((double)k) - a, ((double)(n-k)) + a) ) ) ;
}


@*1 {\bf the constant $C$ in  \eqref{eq:ratesarandall}}. 
\label{sec:C}

compute the constant $C$ in  \eqref{eq:ratesarandall} where we approximate   $m_{\infty} \approx  (2 + (1 + 2^{1-\kappa})/(\kappa - 1))/2 $

@<constantC@>=@#
static double constantC( const double a, const double k, const double c, const double cprime, const double g)
{
  /* \newline $a = \alpha$, $k = \kappa$, $c =c$, $cprime = c^{\prime}$ */
  /* \newline $g = \gamma =  \one{\uN/N \to K} \tfrac{K}{m_{\infty} + K} + \one{\uN/N \to \infty}$  and $ \tfrac{m_{\infty}}{K + m_{\infty}} = 1 - \gamma$ */

  const double m = (2. + (1. + pow(2, 1-k))/(k-1))/2. ;

  /* \newline \eqref{eq:1}; first taking $\kappa = 2$ */
  const double Ck = 1/pow(m,2.) ;

/* \newline $c_{\alpha}^{\prime} =  \one{\alpha = 1}c^{\prime} + \one{ \alpha \in (0,2)\setminus\{1\} }c $ */
  const double caprime =   (a > 1 ? c : (a < 1 ? c : cprime) );

  /* \newline the integral in \eqref{eq:23}  $ \int_{\one{\uN/N \to K}\frac{m_{\infty} }{K + m_{\infty} }  }^{1} (1-y)^{1-\alpha}y^{\alpha - 1}dy  =  \int_0^1 \one{0 < t \le \gamma} t^{1 - \alpha} (1- t)^{\alpha - 1}dt $ */
 /* \newline the incomplete beta function  is normalised by the complete beta function in GSL */
const double tegur =  (g < 1 ? (gsl_sf_beta_inc( 2 - a, a, g)*gsl_sf_beta(2-a,a)) : gsl_sf_beta(2-a,a) ) ;

  return (Ck +  (a*caprime*tegur/pow(m, a))) ;
}


@*1 {\bf $\lambda_{n,k}$  \eqref{eq:ratesarandall}}.
\label{sec:lambdank}

compute the transition rate $\lambda_{n,k}$ in   \eqref{eq:ratesarandall}

@<lambdank@>=@#
static  double lambdank( const int n, const int k, const double a, const double kappa, const double c, const double cprime,  const double gm )
{
   assert( k > 1);
   assert(n > 1) ;
   assert( k <= n);
   /* \newline we approximate $m_{\infty} \approx  \tfrac 12 \left(2 +  \tfrac{1 + 2^{1-\kappa}}{\kappa - 1} \right) $ */
    const double m = (2. + (1. + pow(2, 1-kappa))/(kappa-1))/2. ;
    /* \newline  taking $\kappa = 2$ so that $C_\kappa =  2m_\infty^{-2} $  \eqref{eq:1}  */
      const double Ck = 2/pow(m, 2.) ;
      const double caprime =   (a != 1 ? c : cprime);

/*  \newline   return  \eqref{eq:ratesarandall} using |lbetank|  \S~\ref{sec:totalbetarate} and  |constantC| \S~\ref{sec:C}   */
   return( ((k < 3 ?  (Ck * ((double)(k*(k-1)))/2.) : 0) + (a*caprime*lbetank(n, k, a, gm)/pow(m, a))) / constantC(a, kappa,c, cprime, gm) );
}


@*1 {\bf check the functions}. 
\label{sec:checkfunctions}


check the functions in  \S~\ref{sec:betafunc}, \S~\ref{SEC:binomialconstant}, \ref{sec:totalbetarate}


@<check@>=@#
static void checking(const double ix, const int ib)
{

   for( int k = 3; k <= ib; ++k){

printf( "%f %d %d\n", ix, ib, k);
   printf("%Lf \n",   betafunc( ix, (double)ib, (double)k) );
      printf( "%g\n", log( gsl_sf_beta_inc((double)ib, (double)k, ix) * gsl_sf_beta( (double)ib, (double)k) ) ); 

   printf("%Lf %g %g\n", binom( ib, k),  log( (double)gsl_sf_choose( (unsigned)ib, (unsigned)k)), gsl_sf_lnchoose( (unsigned int)ib, (unsigned int)k) );
   }

 }


@*1 {\bf the jump rate}.
\label{sec:jump}

compute the jump rate of going from $i$ to $j$ blocks $\lambda_{i,  i - k + 1}$ with $\lambda_{n,k}$ the rate of merging $k$ blocks; then $k = n-j+1$

@<jumprate@>=@#
static double qij( const int i, const int j, const double a,  const double kpp, const double c,   const double cpr, const double g)
{

/* \newline |kpp| is $\kappa \ge 2$ */
/* \newline |cpr| is $c_{\alpha}^\prime = c$ when $\alpha \neq 1$, otherwise $c^{\prime} \in (0,1)$ when $\alpha = 1 $ and $\kappa > 2$ */ 
/* \newline |g| is $\gamma$ \eqref{eq:27}  */
assert( i > 1);
assert( j > 0);
assert( j < i); 
  /* \newline    $i \to j = i - k + 1$ so $k = i - j + 1$ ; 
    using |lambdank|  \S~\ref{sec:lambdank} */
    return( lambdank( i, i -j + 1, a, kpp, c, cpr, g) );
}

@*1 {\bf QP}. 
\label{sec:QP}

compute the matrices  $Q$ and $P$

@<QP@>=@#
static void QP( const  int n, const  double a, const  double kpp,  const  double cconst,  const double cprime, const double g,    gsl_matrix * Q, gsl_matrix *P  )
{

  /* \newline  compute matrices |Q| qij and |P| pij */
  /* \newline  |a| is $\alpha \in (0,2)$ */
 /* \newline  |kpp| is $\kappa \ge 2$ */
  /* \newline  |cconst| is $c > 0$ */
   /* \newline  |cprime| is $c^{\prime} \in (0,1)$ */
   /* \newline |g| is $\gamma$  \eqref{eq:27}  */

  int i, j;
  double s = 0;
  double x = 0;
  for( i = 2; i <= n ; i++){
    assert( i <= n);
    s = 0. ;
    for( j = 1; j < i ; j++){
    /* \newline using |qij| from  \S~\ref{sec:jump} */
      x = qij( i, j, a, kpp, cconst, cprime, g); 
      s +=   x;
      gsl_matrix_set( Q, i, j, x);
      gsl_matrix_set( P, i, j, x); }
      assert( s > 0);
    gsl_matrix_set( Q, i, i, s);
    for( j = 1; j < i ; j++){
      assert( j < i);
      gsl_matrix_set( P, i, j,  gsl_matrix_get( P, i, j)/s ); } }
}


@*1 {\bf G}. 
\label{sec:Gmatrix}

compute the $G$ matrix 

@<Gmatrix@>=@#
static void gmatrix( const  int n,  gsl_matrix * G , gsl_matrix * Q, gsl_matrix * P)
{

  int i, k, m ;
  double s = 0.0 ;

  /* \newline  initialise the diagonal */
  for ( i = 2 ; i <= n ; i++){
  assert( gsl_matrix_get( Q, i,i) > 0 );
    gsl_matrix_set( G, i, i, 1./gsl_matrix_get( Q, i,i) );}

  for( i = 3; i <= n ; i++){
    assert( i <= n );
    for( m = 2; m < i; m++){
        s = 0. ;
        for( k = m; k < i ; k++){
          assert( i <= n);
          assert( k <= n);
          assert( m <= n);
          s +=   gsl_matrix_get( P, i, k) * gsl_matrix_get( G, k, m) ; }
      gsl_matrix_set( G, i, m, s); } }
}


@*1 {\bf the matrix $ p^{(n)}[k,b]$}.
\label{sec:pnkb}

 compute the matrix  $p^{(n)}[k,b]$

@<pnkb@>=@#
static void pnb( int n, gsl_matrix * lkb, gsl_matrix * G, gsl_matrix * P  )
{

  /* \newline  $j$ is $n^{\prime}$  nprime from  \cite[Prop A1]{BBE2013a}  */
  /* \newline  |lnb| is the matrix  $p^{(nprime)}[k,b]$; used for each fixed k  */
  gsl_matrix * lnb = gsl_matrix_calloc( n+1, n+1);
  
  int  k, b, j, i ;
  double s = 0.0 ;
  gsl_matrix_set( lkb, n, 1, 1.0 ); 
  for ( k = 2 ; k < n ; k++){
    for( i = k ; i <= n ; i++){
    for ( b = 1; b <= i - k + 1 ; b++){
      gsl_matrix_set( lnb, i, b, ( k == i ? ( b == 1 ? 1.0 : 0.0) : 0.0 ) );
      s = 0. ;
      for( j = k ; j < i ; j++){
        gsl_matrix_set( lnb,  i, b, gsl_matrix_get( lnb, i,b) +  (b > i - j ? (((double)(b - i + j)) * gsl_matrix_get( lnb, j, b - i + j) * (gsl_matrix_get( P, i, j) * gsl_matrix_get( G, j, k) / gsl_matrix_get( G, i, k)) / ((double)j)) : 0.0 ) +  (b < j ? ((((double)(j-b))*gsl_matrix_get( lnb, j, b) *  (gsl_matrix_get( P, i, j) * gsl_matrix_get( G, j, k) / gsl_matrix_get( G, i, k)) / ((double)j) )) : 0.0 ) ) ;
      } } }

    for( j = 1 ; j <= (n - k + 1) ; j++){
      gsl_matrix_set( lkb, k, j, gsl_matrix_get( lnb, n, j) ); }
    gsl_matrix_set_zero( lnb ); }

  gsl_matrix_free( lnb ); 
}


@*1 {\bf compute the expected spectrum}.
\label{sec:ebi}


compute the expected spectrum  $\EE{L_{i}(n)}/\EE{L(n)}$


@<ebi@>=@#
static void ebi( const  int n, const  double a, const  double kappa, const double cconst, const double cprime, const double gm )
{

  /* \newline |n| is sample size */
  /* \newline  |a| is $\alpha \in (0,2)$ */
  /* \newline |kappa| is $\kappa \ge 2$ */
  /* \newline |cconst| is $c > 0$ */
  /* \newline |cprime| is $c^{\prime} \in (0,1)$ */
  /* \newline |gm| is $\gamma$   \eqref{eq:27}  */

  gsl_matrix * P = gsl_matrix_calloc( n+1, n+1);
  gsl_matrix * Q = gsl_matrix_calloc( n+1, n+1);
  gsl_matrix * G = gsl_matrix_calloc( n+1, n+1);
  gsl_matrix * Pn = gsl_matrix_calloc( n+1, n+1);

/* \newline |iA| is $\kappa$ */
  const double iA = kappa ; 
  const double minfty =  (2. +  (1. + pow(2.,  1. - iA))/(iA - 1.))/2. ;

  /* \newline  compute the |Q| and |P|  matrices \S~\ref{sec:QP}  */
  QP( n, a, iA, cconst, cprime, gm,    Q, P);
  /* \newline  the  |G| matrix \S~\ref{sec:Gmatrix} */
  gmatrix( n,  G, Q, P);
  /* \newline  the   $p^{(n)}[k,b]$ matrix  \S~\ref{sec:pnkb} */
  pnb( n, Pn, G, P);
  
  int b, k;
  double s = 0.0 ;
  /* \newline $eb$ will store $\EE{L(n)}$ */
  double eb = 0.0 ;
  double * ebi  = (double *)calloc( n, sizeof(double) ); 

  for ( b = 1; b < n ; b++){
    s = 0. ;
    for( k = 2; k <= n-b+1 ; k++){
      s = s +  (gsl_matrix_get( Pn, k, b) * ((double)k) * gsl_matrix_get( G, n, k)); }
    ebi[b] = s ; 
    eb +=  s ;}
  assert( eb > 0);
  double x ;
  for( b = 1 ; b < n ; ++b){
      x = ((double)b)/((double)n) ;
      printf("%g,%g\n", log(x) - log(1-x),  log(ebi[b]/eb) - log(1.0 -  (ebi[b]/eb)) ) ;}

  gsl_matrix_free( P);
  gsl_matrix_free( Q);
  gsl_matrix_free( G);
  gsl_matrix_free( Pn);
  free( ebi);
}


@*1 {\bf the main function}. 
\label{sec:main}

@C

/* \newline \S~\ref{sec:includes} */
@<includes@>@#
/* \newline \S~\ref{sec:betafunc} */
@<beta function@>@#
/* \newline \S~\ref{SEC:binomialconstant} */
@<binomial@>@#
/* \newline \S~\ref{sec:totalbetarate} */
@<betarate@>@#
/* \newline \S~\ref{sec:C} */
@<constantC@>@#
/* \newline \S~\ref{sec:lambdank} */
@<lambdank@>@#
/* \newline \S~\ref{sec:checkfunctions} */
@<check@>@#
/* \newline \S~\ref{sec:jump} */
@<jumprate@>@#
/* \newline \S~\ref{sec:QP} */
@<QP@>@#
/* \newline \S~\ref{sec:Gmatrix} */
@<Gmatrix@>@#
/* \newline \S~\ref{sec:pnkb} */
@<pnkb@>@#
/* \newline \S~\ref{sec:ebi} */
@<ebi@>@#
int main( int argc, char *argv[])
{

/* \newline  |ebi| \S~\ref{sec:ebi} */
/* \newline   ebi( const  int n, const  double a, const  double kappa, const double cconst, const double cprime, const double gm ); */
 /*  ebi( atoi(argv[1]), atof(argv[2]), atof(argv[3]),  atof(argv[4]),  atof(argv[5]),  atof(argv[6]) ); */
ebi( 100, 0.01, 2, 1, 1, 1 );

@#

return GSL_SUCCESS ;
}



@* {\bf conclusions and bibliography}. 
\label{sec:concl}


For $i = 1,2, \ldots, n-1$ where $n$ is sample size   we compute
\begin{equation}
\label{eq:master}
\EE{L_{i}(n)} =  \sum_{k=2}^{n-i+1} k  p^{(n)}[k,i] g(n,k)
\end{equation}
by a recursion;
$g(n,m) = \EE{ \int_{0}^{\infty} \one{ \# \xi^{n}(s) = m }ds :
\xi^{(n)}(0) = n }$ is the expected amount time $\set{\xi^{n}}$ has
$m$ blocks when starting with $n$, and $p^{(n)}[k,b]$ is the
probability a given branch when there are $k$ branches will support
exactly $b$ leaves; see \cite{BBE2013a} for details.
Equation~\eqref{eq:master} can be used to compute $\EE{L_{i}(n)} $ for
any coalescent provided $ p^{(n)}[k,i] g(n,k)$ can be computed;
\cite{Blath2016} apply \eqref{eq:master} to $\Xi$-coalescents
admitting simultaneous multiple mergers.  On the implementation side,
$p^{(n)}[k,i] $ involves nested loops. 



Figure~\ref{fig:graph} records an example 


\begin{SCfigure}[0.8][htb]
    \begin{tikzpicture}
      \begin{axis}[
        xlabel = $\log(i/n) - \log(1 - i/n)$,
        axis line style = {draw = none},
        tick style = {draw=none},
        xticklabels={draw=none},
        yticklabels={draw=none}]
        \addplot table[col sep=comma] {forplottingfile1};
        %\addplot table[col sep=comma] {forplottingfile2};
       \end{axis}
       \end{tikzpicture}
       \caption{\textcolor{violet}{\it  An example approximation of $\EE{R_{i}(n)}$
graphed as logits against $\log(i/n) - log(1-i/n)$ where $n$ is sample size} }
      \label{fig:graph}
       \end{SCfigure}






\bibliographystyle{alpha}
\begin{thebibliography}{Tan11}


\bibitem[Blath2016]{Blath2016}
	{Jochen Blath and Mathias Christensen Cronj\"{a}ger and Bjarki Eldon and Matthias Hammer}.
	\newblock {Theoretical Population Biology}
	\newblock {36--50},
	\newblock {The site-frequency spectrum associated with {$\Xi$}-coalescents},
	\url{https://dx.doi.org/10.1016/j.tpb.2016.04.002},
	\newblock {110}, {2016}.


\bibitem [BBE2013a]{BBE2013a}
	{M Birkner and J Blath and B Eldon}.
	\newblock {Genetics}, {1037--1053}.
	\newblock {Statistical properties of the site-frequency spectrum associated with {$\Lambda$}-coalescents}.
	\newblock  195, 2013.



\bibitem[Eld24]{eldon24:_coales}
Bjarki Eldon.
\newblock Gene genealogies in haploid populations evolving according to sweepstakes reproduction.
\newblock In preparation, 2024+.

\bibitem[KL94]{knuth1994cweb}
Donald~Ervin Knuth and Silvio Levy.
\newblock {\em The CWEB system of structured documentation: version 3.0}.
\newblock Addison-Wesley Longman Publishing Co., Inc., Reading, Massachusetts,
  1994.

\bibitem[KR88]{kernighan1988c}
Brian~W Kernighan and Dennis~M Ritchie.
\newblock The {C} programming language, 1988.

\bibitem[Sch03]{schweinsberg03}
J~Schweinsberg.
\newblock Coalescent processes obtained from supercritical {G}alton-{W}atson
  processes.
\newblock {\em Stoch Proc Appl}, 106:107--139, 2003.

\bibitem[Tan11]{tange11:_gnu_paral}
O~Tange.
\newblock {GNU} parallel -- the command-line power tool.
\newblock The USENIX Magazine, 2011.

\end{thebibliography}


@
\end{document}
