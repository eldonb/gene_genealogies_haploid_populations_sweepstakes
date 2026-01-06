%$ NAFN=quenched_mixed_bounded_sfswplotting
%$ echo 'const bool MODEL = 0;' > $NAFN.hpp
%$ echo 'const double ALPHA = 1;' >> $NAFN.hpp
%$ echo 'const double KAPPA = 2;' >> $NAFN.hpp
%$ echo 'const std::size_t POP_SIZE = 1e2;' >> $NAFN.hpp
%$ echo 'const std::size_t CUTOFF = POP_SIZE;' >> $NAFN.hpp
%$ echo 'const double EPSILON = log((double)POP_SIZE)/((double)POP_SIZE);' >> $NAFN.hpp
%$ echo 'const std::size_t SAMPLE_SIZE = 10;' >> $NAFN.hpp
%$ echo 'const int EXPERIMENTS = 2500;' >> $NAFN.hpp
%$ ctangle $NAFN.w
%$ g++ -std=c++26 -m64 -march=native -O3 -x c++ $NAFN.c -lm -lgsl -lgslcblas
%$ rm -f tmp_runs
%$ for i in $(seq 4); do echo "./a.out 10 $(shuf -i 432293-31911194833 -n1) > resout"$i >> tmp_runs ; done
%$ parallel --gnu -j4 :::: ./tmp_runs
%$ paste resout* -d, | sed '1d' | awk -F, '{s=0; for (i=1;i<=NF;i++) {s+=$i} print log(s/1e4) - log(1-(s/1e4))}' > logitresout
%$ seq 9 | awk '{S=10;print log($1/S) - log(1 - ($1/S))}' > nlogits
%$ paste -d',' nlogits logitresout > forplottingfile1
%$ sed -i 's/SAMPLE_SIZE = 10/SAMPLE_SIZE = 100/g' $NAFN.hpp
%$ ctangle $NAFN.w
%$ g++ -std=c++26 -m64 -march=native -O3 -x c++ $NAFN.c -lm -lgsl -lgslcblas
%$ rm -f tmp_runs
%$ for i in $(seq 4); do echo "./a.out 100 $(shuf -i 432293-31911194833 -n1) > resout"$i >> tmp_runs ; done
%$ parallel --gnu -j4 :::: ./tmp_runs
%$ paste resout* -d, | sed '1d' | awk -F, '{s=0; for (i=1;i<=NF;i++) {s+=$i} print log(s/1e4) - log(1-(s/1e4))}' > logitresout
%$ seq 99 | awk '{S=100;print log($1/S) - log(1 - ($1/S))}' > nlogits
%$ paste -d',' nlogits logitresout > forplottingfile2
%$ NAFN=quenched_mixed_bounded_sfswplotting
%$ cweave $NAFN.w
%$ tail -n4 $NAFN.tex > endi
%$ for i in $(seq 5); do $(sed -i '$d' $NAFN.tex) ; done
%$ cat endi >> $NAFN.tex
%$ emacs --version | head -n1 > innleggemacs
%$ g++ --version | head -n1 > innleggcpp
%$ xelatex --version | head -n1  > innleggxelatex
%$ cweave  --version | head -n1 > innleggcweave
%$ ctangle  --version | head -n1  > innleggctangle
%$ uname  --kernel-release -o  > innleggop
%$ bash --version | head -n1 > innleggbash
%$ sed -i 's/x86/x86\\/g' innleggbash
%$ gsl-config --version > innlegggsl
%$ parallel --version | head -n1 > innleggparallel
%$ awk --version | head -n1 > innleggawk
%$ spix --version > innleggspix
%$ xelatex $NAFN.tex
%$ xelatex $NAFN.tex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Copyright (C) 2025 Bjarki Eldon
\documentclass[a4paper,10pt,final]{cweb}
\usepackage{fancyvrb,  graphicx, url}
\usepackage[svgnames]{xcolor}
\usepackage{minted}
%% a4 paper size 210 x 297 millimeters
%% for xelatex
\usepackage{xunicode}
\usepackage{fontspec}
\usepackage{xltxtra}
\usepackage{lineno}
\usepackage[all]{xy}
%%\usepackage[bigdelims,vvarbb]{newtxmath}
\usepackage{amsfonts, amsmath, amssymb}
\usepackage{fullpage}
\usepackage{marvosym}
\usepackage{bm}
\usepackage{natbib}
%\usepackage[backend=biber,style=authoryear-icomp,
%    sortlocale=de_DE,
%    natbib=true,
%    url=false, 
%    doi=true,
%    eprint=false]{biblatex}
%    \addbibresource{refs.bib}
    \usepackage[]{hyperref}
%    \hypersetup{colorlinks=true}
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
\newcommand{\x}{{\circ}}
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
%%\pagecolor{DarkGray}
\title{Gene genealogies in haploid populations evolving according to
sweepstakes reproduction \\ --- approximating $\EE{\widetilde R_{i}^N(n)}$}
\author{Bjarki Eldon\footnote{\href{beldon11@@gmail.com}{beldon11@@gmail.com}} \footnote{
%%  Deutsche Forschungsgemeinschaft (DFG) - Projektnummer 273887127 
%% „funded by the Deutsche Forschungsgemein-schaft (DFG, German Research Foundation) –Projektnummer(n)“.
%%  
 compiled @@
{\DTMcurrenttime} on  {\today}  \\
\input{innleggctangle} \\
 \input{innleggcpp} \\  \input{innleggop} \\  \input{innleggbash} \\
GSL \input{innlegggsl} \\ \input{innleggcweave} \\  {\XeLaTeX}
\input{innleggxelatex}  \\ \input{innleggparallel} \\ \input{innleggawk} \\ \input{innleggspix}\\  \input{innleggemacs}}
\orcidlink{0000-0001-9354-2391} }

\begin{document}
\maketitle
\renewcommand{\abstractname}{\vspace{-\baselineskip}}

\begin{abstract}
Consider a haploid panmictic population of constant size $N$ evolving
according to some given model.  Suppose we {\it (i)} evolve the
population forward in time for an infinite amount of time and record
the ancestral relations of everyone in the population at any time;
{\it (ii)} at some given time we record a random sample from the
population. Since we know the entire population ancestry the tree
relating the sampled individuals is fixed; {\it (iii)} trace the
sample tree and record the branch lengths of the fixed sample tree.
After erasing the population ancestry and the information about the
sampled individuals we repeat steps {\it (i)}--{\it (iii)} a given
number of times, each time starting from scratch with a new population
ancestry independent of any previous (and  future) ancestries.  We
denote the mean relative branch lengths obtained in this way with
$\EE{\widetilde R_{i}^{N}(n)}$. With this C++ simulation  code we
 approximate $\EE{\widetilde R_{i}^{N}(n)}$
\end{abstract}



\tableofcontents


@* {\bf intro}.
\label{sec:intro}



Gene genealogies are generated through inheritance as the population
evolves  forward in time. This means that the individuals (gene
copies) in a given sample are related through  one  fixed
tree. Coalescents  (Markov processes on  partitions of positive
integers  describing the random ancestral relations of  sampled gene
copies)  are derived by describing the ancestral relations of a sample
one generation at a time going backwards in time from the sample, thus
ignoring the fact that the sampled  gene copies are related through a
fixed complete  tree.    


Here we provide the means to investigate functionals of fixed complete
trees.  An effort to illustrate the idea may be found in
Figure~\ref{fig:illsampling}.  Suppose a haploid panmictic population
of constant size has evolved according to a given model for an
infinite number of generations
(Figure~\ref{fig:illsampling}$(i)$). The ancestry of all the
individuals in the population at any time is known and is stored in a
sigma field. At time 0 we record a random sample (the $\s$ in
Figure~\ref{fig:illsampling}$(ii)$).  Since the ancestry of the
population is known arbitrarily far back in time, the sampled gene
copies will have a fixed complete tree given the population ancestry
(marked out by the $\s$ in Figure~\ref{fig:illsampling}$(iii)$).  



Write $[n] \equiv \set{1,2,\ldots,n}$ for any
$n\in \IN \equiv \set{1,2,\ldots}$.  Let $ \widetilde L_{i}^{N}(n)$
denote the random length of branches supporting $i$ leaves (sampled
gene copies) for all $i \in [n-1]$ read off a fixed tree given a
random population ancestry.  Write $\widetilde R_{i}^{N}(n) \equiv
\widetilde L_{i}^{N}(n)/\sum_{j}\widetilde L_{j}^{N}(n)$. We are
interested in  approximating  $\EE{ \widetilde R_{i}^{N}(n) }$ by
averaging realisations of   $\widetilde R_{i}^{N}(n)$ over
independent ancestries (Figure~\ref{fig:illancestries}).   In
Figure~\ref{fig:illancestries} we repeat the process in
Figure~\ref{fig:illsampling} some given number  $(M)$ of times, each
time reading branch lengths off fixed complete trees.     For example,
for the gene genealogy in Figure~\ref{fig:forw} we have $\widetilde L_{1}^{N}(3)
= 4$, and $\widetilde L_{2}^{N}(3) = 1$. 



\begin{figure}[htp]
\centering
\[
\xymatrix @@R=0.05pc @@C=0.4pc @@M=1pt {%
& \text{$(i)$  population ancestry} & & \text{$(ii)$ random sample} &&
\text{$(iii)$ complete sample tree} \\ 
\text{past} \ar[ddd] &  \x\x\x\x\x\x\x\x\x\x\x\x\x\x &&  \x\x\x\x\x\x\x\x\x\x\x\x\x\x &&  \x\x\x\x\x\x\x\x\x\x\widehat\s\x\x\x  \\
 &  \x\x\x\x\x\x\x\x\x\x\x\x\x\x && \vdots  &&   \vdots  \\ 
 &  \vdots &&  \x\x\x\x\x\x\x\x\x\x\x\x\x\x &&    \x\s\x\x\x\x\x\s\x\x\x\x\x\x \\
\text{present}  &  \x\x\x\x\x\x\x\x\x\x\x\x\x\x && \to  \x\x\x\s\x\x\s\x\x\x\s\x\x\x && \vdots  \\
&&&&&  \to  \x\x\x\s\x\x\s\x\x\x\s\x\x\x  \\
}
\]
\caption{A fixed complete sample tree is obtained in the  three steps
shown; the most recent common ancestor of the sampled gene copies
($\s$) is marked  $\widehat \s$ }
\label{fig:illsampling}
\end{figure}


\begin{figure}[htp]
\centering
\[
\xymatrix @@R=0.1pc @@C=0.4pc @@M=1pt  {%
&  \text{ancestry 1} && \text{ancestry 2}  & \ldots & \text{ancestry $M$} \\
\text{past}&\x\x\x\x\x\widehat\s\x\x\x\x\x\x\x\x  & &    \x\x\x\x\x\widehat\s\x\x\x\x\x\x\x\x &  &   \x\x\x\x\x\widehat\s\x\x\x\x\x\x\x\x \\
& \x\s\x\x\x\s\x\x\x\x\x\x\x\x  &&    \x\x\s\x\x\s\x\x\x\s\x\x\x\x &  &    \x\x\x\x\x\s\x\x\x\x\x\s\x\x \\
&\vdots                        && \vdots & &  \vdots  \\
& \x\s\x\x\x\x\x\x\x\x\s\x\x\x &&  \x\s\x\x\x\x\s\x\x\x\x\s\x\x &&   \x\x\x\x\s\x\x\x\x\x\x\x\s\x   \\
& \x\x\s\x\x\x\s\x\x\x\s\x\x\x  &&  \to \x\x\s\x\x\x\s\x\s\x\x\x\x\x  &&   \to \x\x\x\x\x\x\x\x\s\x\s\x\s\x   \\
& \vdots                       & & &   \\
\text{present} \ar[uuuuuu]^{\text{\rotatebox{90}{time in generations}}} & \to  \x\s\x\x\x\x\s\x\x\x\x\x\s\x  & -\ar@@{-}[uuuuuuu]  &   & -\ar@@{-}[uuuuuuu] \\
}
\]
\caption{Complete sample trees within independent  population ancestries}
\label{fig:illancestries}
\end{figure}




At any given time the current individuals independently produce potential offspring according to a given law. 
For all $k \in [\zeta(N)]$ suppose 
\begin{equation}
\label{eq:1}
\prb{X = k} =  C\left( k^{-\alpha} - (1+k)^{-\alpha}\right)
\end{equation}
where $C > 0$ is chosen such that $\prb{X\in [\zeta(N)]} = 1$.  The
quantity $\zeta(N) \in \IN$ is an upper bound on the number of
potential offspring an individual may produce.  From the pool of
potential offspring we sample (uniformly without replacement) a fixed
number to survive.  Let $\svigi{\varepsilon_{N}}_{N}$ be a sequence of
positive constants where $0 < \varepsilon_{N} < 1$ and it may hold
that   $\varepsilon_{N}\to 0$ as $N\to \infty$. Write
\begin{displaymath}
X \vartriangleright L(a,\zeta(N))
\end{displaymath}
when $X$ is distributed according to \eqref{eq:1} with $a$ and
$\zeta(N)$ as given.   


We will consider two models for the law on the number of potential
offspring, `Model 0` and `Model 1`.   We now define  Model 0.  
Let $E$ be the event
\begin{equation}
\label{eq:2}
\set{ X_{i} \vartriangleright L(\alpha,\zeta(N) \text{ for all }  i \in [N] }
\end{equation}
and $E^{\sf c}$ the event where $\kappa$ replaces $\alpha$ in $E$ in \eqref{eq:2}. 
Suppose $\prb{E} = \varepsilon_{N}$ and $\prb{E^{\sf c}} = 1 -
\varepsilon_N$.  We call this `Model 0'  ($|MODEL| = 0$).


We now define   Model 1.  
Let $E_{1}$ be the event
\begin{multline}
\label{eq:3}
\left \{ \text{there exists exactly one $i \in [N]$ with } X_{i}
\vartriangleright L(\alpha,\zeta(N)); \right. \\  \left. \text{moreover,
it holds that  }  X_{j} \vartriangleright
L(\kappa,\zeta(N)) \text{ for all }  j \in [N]\setminus\set{i}  \right\}
\end{multline}
Let $E_{1}^{\sf c}$ be the event  where $\kappa$ replaces  $\alpha$ in
$E_{1}$; when $E_{1}^{\sf c}$ occurs it holds that
$X_{i}\vartriangleright L(\kappa,\zeta(N))$ for all $i \in [N]$.  
Suppose $\prb{E_{1}} = \overline \varepsilon_{N}$ and $\prb{E_{1}^{\sf c}}
= 1 - \overline\varepsilon_N$  and it may hold that $\overline
\varepsilon_N \to 0$ as $N\to \infty$.   This is `Model 1' ($|MODEL| =
1$).   



We record the ancestry, i.e.\ the immediate
ancestor of every individual in the population at any time. Thus, as
the population evolves forward in time we generate gene
genealogies. Given a random sample, the sampled gene copies share a
fixed complete tree.








In \S~\ref{sec:code} we summarise the algorithm, the code follows in
\S~\ref{sec:includes} -- \S~\ref{sec:main}, we conclude in
\S~\ref{sec:concl}. Comments within the code are in \aths{this font
and colour}


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



To produce this  document one can use the shell tool  {\tt spix} on the source  file

{\tt spix <prefix>.w}

thus running the shell script appearing  before 
the preamble in the source file  ({\tt .w} file),  the lines beginning with \%\$


One may also copy the shell script (minus the leading \%\$) into a
file (say, scriptfile) and use {\tt parallel} \cite{tange11:_gnu_paral}


{\tt parallel --gnu -j1 :::: ./scriptfile}



Use {\tt valgrind} to check for memory leaks:


{\tt valgrind -v --leak-check=full --show-leak-kinds=all <program call>}


Use {\tt cppcheck} to check the code


{\tt cppcheck  --enable=all --language=c++ <prefix>.c}



Figure~\ref{fig:res} holds  an example output 



@* {\bf code}.
\label{sec:code}

we  briefly summarise the main steps of the algorithm 


\begin{enumerate}
\item $(r_{1}^{N}(n), \ldots, r_{n-1}^{N}(n)) \leftarrow (0,\ldots,
0)$
\item for each of $M$ experiments \S~\ref{sec:approximate}
\begin{enumerate}
\item $\svigi{\ell_{1}^{N}(n), \ldots, \ell_{n-1}^{N}(n)} \leftarrow
(0,\ldots, 0)$
\item initialise population ancestry $a_{i}(0) = i$ for $i=0,1,\ldots,N-1$
\item sample (uniformly without replacement)  random number $n$ of
levels \S~\ref{sec:randomsample}
\item {\bf while} sample tree is not complete: \S~\ref{sec:checkifcomplete}
\begin{enumerate}
\item add one generation to the population ancestry
\S~\ref{sec:addgeneration}
\item take a new random sample  \S~\ref{sec:randomsample}
\end{enumerate}
\item read lengths $\ell_{i}^{N}(n)$ off the fixed complete sample
tree \S~\ref{sec:onexperiment}
\item given the lengths $\ell_{i}^{N}(n)$ update the  approximations $r_{i}^{N}(n)$ \S~\ref{sec:updateri} 
\end{enumerate}
\item return $r_{i}^{N}(n)/M$ as an approximation of $\EE{\widetilde
R_{i}^{N}(n)}$ for $i = 1,2,\ldots, n-1$
\end{enumerate}


The population evolves as a supercritical branching process
conditioned to be of a fixed size.  We record the ancestry  of the
individuals in the population at any time.  The individuals at any
time are seen as occupying levels, one individual on each of $N$
levels.  The  individual on level $i$ at time $g$ has immediate
ancestor  $a_{i}(g)$, where $a_{i}(g)$  is the level of the parent of
the given  individual.  The offspring  of the individual on level $i$
will all have value $i$; each  offspring `points' to its' immediate
ancestor (Figure~\ref{fig:anc})
\begin{figure}[htp]
\centering
\[
\xymatrix @@R=0.1pc @@C=0.4pc  {%
\text{time} & \text{levels} \\
     & 0  & 1 & \ldots & \ell-1 &  \ell & \ell+1 &  \ldots & N-1 \\
g-1  & \x & \ldots & &\x & \x & \x &\ldots & \x \\ 
g & \x & \ldots  & \x & \x \ar[ur] & \x \ar[ur] & \x \ar[ul] & \ldots & \x  \\
}
\]
\caption{The immediate ancestor of the individuals on levels  $\ell-1$
and $\ell + 1$ at time $g$ occupies level $\ell$ at time $g-1$
($a_{\ell-1}(g) = a_{\ell+1}(g) = \ell$); and
$a_{\ell}(g) = \ell + 1$ }
\label{fig:anc}
\end{figure}
Gene genealogies evolve forward in time (Figure~\ref{fig:forw}).     The individuals (gene
copies) of any given sample will share one complete gene
genealogy. Given that we know the population ancestry we can use this
information to  read gene genealogical statistics off fixed complete
sample trees.
\begin{figure}[htp]
\centering
\[
\xymatrix @@R=0.05pc @@C=0.2pc  {%
\text{time} & \text{levels} \\
&   0  & 1 & \ldots & \ell-1 &  \ell & \ell+1 &  \ldots & N-1 \\
g-2& \x & \ldots &&  \x  & \widehat \s \ar[dl] \ar[dr] & \x &\ldots & \x \\
g-1& \x & \ldots & \x & \s\ar[dl]  & \x & \s \ar[dl]  \ar[dr]  & \ldots & \x \\
g& \x & \ldots &  \s & \x &\s & \x &\s &\x &  \\
}
\]
\caption{The gene genealogy of the 3 individuals sampled at time $g$}
\label{fig:forw}
\end{figure}



@*1 {\bf the included libraries}.
\label{sec:includes}

the included libraries

@<included libraries@>=@#
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
#include <gsl/gsl_math.h>
#include "quenched_mixed_bounded_sfswplotting.hpp"


@*1 {\bf the random number generators}.
\label{sec:rngs}

define the random number generators

@<random number generators@>=@#
/* \newline \aths{ obtain a seed out of thin air for the random number engine} */
  std::random_device randomseed;
  /* \newline \aths{ Standard mersenne twister  random number engine seeded with |rng()|} */
  std::mt19937_64 rng(randomseed());

gsl_rng * rngtype ;

@#


void setup_rng( unsigned long int s )
{

@#

        const gsl_rng_type *T ; 
        gsl_rng_env_setup(); 
        T = gsl_rng_default ;
        rngtype = gsl_rng_alloc(T);
        gsl_rng_set( rngtype,  s) ;
}




@*1 {\bf probability mass function for number of offspring}.
\label{sec:mass}


The mass function \eqref{eq:1} 

@<mass function@>=@#
double  kernel(const  std::size_t k,  const double a)
{
 /* \newline \aths{|a| is either $\alpha$  |ALPHA| or $\kappa$  |KAPPA|  } */
  return ( pow(1./static_cast<double>(k), a) - pow( 1./static_cast<double>(k+1), a) ) ;
}


@*1 {\bf generate cumulative mass function}. 
\label{sec:cmf}


generate the cumulative mass function for sampling numbers of
potential offspring 


@<cmf for number of offspring@>=@#
void generate_cmfs(  std::vector<double>& cmfa,  std::vector<double>&  cmfA  )
{
@#

 double s {} ;
  double d {} ;
  for( std::size_t i = 1; i <= CUTOFF; ++i){
    cmfa[i] = cmfa[i-1] + kernel( i, ALPHA);
    cmfA[i] = cmfA[i-1] + kernel( i, KAPPA);
    s += kernel(i, ALPHA);
    d += kernel(i, KAPPA);}

@#

  assert(s > 0.);

@#

  assert(d > 0);


@#

  std::transform( cmfa.begin(), cmfa.end(), cmfa.begin(), [&s](const double &x){return x/s ;} );
  std::transform( cmfA.begin(), cmfA.end(), cmfA.begin(), [&d](const double &x){return x/d ;} );

@#

  assert( abs( cmfa.back() - 1.) < 0.999999 );
  assert( abs( cmfA.back() - 1.) < 0.999999 );

  assert( abs(cmfa[0]) < 0.000001 );
  assert( abs(cmfA[0]) < 0.000001 );
}


@*1 {\bf sample one bounded random number of potential offspring}.
\label{sec:boundedX}


sample one bounded random number of potential offspring using
\eqref{eq:1}


@<one bounded number@>=@# 
static std::size_t randomX( const  std::vector<double>& f )
{
@#
 /* \newline \aths{|rngtype| \S~\ref{sec:rngs}} */
  const double u = gsl_rng_uniform_pos( rngtype);

@#

  std::size_t j {1} ;
 
  while( u > f[j]){ ++j ; assert(j <= CUTOFF); }

@#
  return j ;
}



@*1 {\bf one unbounded random number of potential offspring}.
\label{sec:unboundedX}


sample one unbounded random number of potential offspring as  $
\lfloor U^{-1/\alpha}\rfloor $ where $U$ is a standard random uniform
and  $\alpha > 0$ is fixed 

@<an unbounded $X$@>=@#
static std::size_t unboundedrandomX(  )
{
  return  static_cast<std::size_t>( floor( pow( 1./gsl_rng_uniform_pos( rngtype), 1./ALPHA) ) ) ;
}



@*1 {\bf add one generation to population ancestry}.
\label{sec:addgeneration}


add one generation to population ancestry; the offspring (if any) of
the individual   on level $i$ will have value $i$ (point to
$i$). Recall ${X}_{1}, \ldots, X_{N}$ are the random numbers  of
potential offspring produced by the current $N$ individuals, by
\eqref{eq:1} $X_{1} + \cdots + X_{N} \ge N$ almost surely. We add to
the ancestry in the following 3 steps:
\begin{enumerate}
\item \label{item:1}  sample $X_{1},\ldots, X_{N}$ and append indexes  $( \ldots, \underset{X_{1}}{\underbrace{1, \ldots, 1}}, \ldots,   \underset{X_{N}}{\underbrace{N, \ldots, N}})$
\item \label{item:2}  shuffle the new $X_{1}+\cdots + X_{N}$ indexes  $\left( \ldots,
\sigma_{1}, \ldots, \sigma_{X_{1} + \cdots + X_{N}} \right)$
\item \label{item:3}  keep the first $N$ of the new  (shuffled) indexes $ \left ( \ldots, \sigma_{1}, \ldots, \sigma_{N} \right) $
\end{enumerate}


@<one generation@>=@#
void  addgeneration_modelnull( std::vector<std::size_t>& tree, const std::vector<double>& vcmf )
{


@#

  const std::size_t e = tree.size() ; 

  /* \newline \aths{step \ref{item:1}} */
  std::size_t x {} ;

@#

  for( std::size_t i = 0 ; i < POP_SIZE; ++i){
  /* \newline \aths{ |randomX| \S~\ref{sec:boundedX} } */
    x = randomX(vcmf) ;
    tree.reserve(tree.size() + x);
    /* \newline  \aths{ the |x| offspring of individual on
    level $i$ point to $i$; recall Figure~\ref{fig:forw}} */
    tree.insert( tree.end(), x, i); }

  assert(tree.size() >= e + POP_SIZE); 

  /* \newline \aths{step \ref{item:2}}  */
  std::shuffle(tree.begin() + e, tree.end(), rng); 

  /* \newline \aths{step \ref{item:3} } */ 
  tree.resize(e + POP_SIZE) ;

  /* \newline \aths{population size $N$   |POP_SIZE|
   is fixed so  |POP_SIZE|  must divide
  the size of the ancestry} */
   assert( (tree.size() % POP_SIZE) == 0);
}



@*1 {\bf add one generation to population ancestry when Model 1}.
\label{sec:onegenerationM1}


add one generation to the population ancestry when  Model 1 holds   

@<add a generation for model 1@>=@#
void addgeneration_modelone( std::vector<std::size_t>& tree, const std::vector<double>& vcmf_a, const std::vector<double>& vcmf_A )
{
    const std::size_t e = tree.size() ;

    const double u = gsl_rng_uniform_pos(rngtype);
    std::size_t x  =  randomX( u < EPSILON ? vcmf_a : vcmf_A ) ;
  
  for( std::size_t i = 1 ; i < POP_SIZE; ++i){
    x = randomX( vcmf_A ) ;
    tree.reserve(tree.size() +  x);
    tree.insert( tree.end(), x, i);
  }

  assert( tree.size() >= e + POP_SIZE);
  
  std::shuffle( tree.begin() + e, tree.end(), rng);

  tree.resize( e + POP_SIZE);


  assert( (tree.size() % POP_SIZE) == 0 );
}





@*1 {\bf get a random sample}.
\label{sec:randomsample}


 get a random sample; return the levels of the sampled individuals
 (gene copies) 


@<random sample@>=@#
static void randomsample( std::vector< std::pair< std::size_t, std::size_t>>& sample )
{

@#

  /* \newline  \aths{pair is (size of block,  level of block)} */
  std::vector< std::size_t> V (POP_SIZE, 0);
  std::iota( V.begin(), V.end(), 0);
  std::shuffle( V.begin(), V.end(), rng);
  sample.clear();
  assert( sample.size() < 1);

  sample.reserve(SAMPLE_SIZE); 
    for( std::size_t i = 0 ; i < SAMPLE_SIZE ; ++i){
    /* \newline \aths{ pair is (size of block, level of block)} */
    sample.push_back( std::make_pair( 1, V[i]   ) ) ; }
   assert(sample.size() == SAMPLE_SIZE);
}


@*1 {\bf look up the immediate ancestor}.
\label{sec:getagi}

look up the immediate ancestor; $a_{i}(g) \in \set{0,1,\ldots, N-1} $ is the level of the
immediate ancestor of  the individual on level $i$ at time $g$. The
population  ancestry is the vector ($N$ is population size)
\begin{displaymath}
\left(\ldots,  a_{0}(g-1), a_{0}(g-1),\ldots, a_{N-1}(g-1), \ldots,  a_{0}(g),\ldots, a_{N-1}(g),  a_{0}(g+1),\ldots, a_{N-1}(g+1), \ldots  \right)
\end{displaymath}
where  $a_{i}(0)=i $ for $i = 0,1,\ldots, N-1$; 
\begin{displaymath}
\begin{split}
\text{time } & \quad  \text{level}\\
0 & \quad  0, 1, \ldots , N-1 \\
1 & \quad  a_{0}(1), a_{1}(1), \ldots, a_{N-1}(1) \\
\vdots \\
g-1 &  \quad  a_{0}(g-1), a_{1}(g-1), \ldots, a_{N-1}(g-1) \\
g &  \quad  a_{0}(g), a_{1}(g), \ldots, a_{N-1}(g) \\
g+1 &  \quad  a_{0}(g+1), a_{1}(g+1), \ldots, a_{N-1}(g+1) \\
\vdots \\
\end{split}
\end{displaymath}



@<get immediate ancestor@>=@#
static std::size_t getagi( const std::size_t &g, const std::size_t &i,  const std::vector<std::size_t>& tree )
{

@#

/* \newline \aths{return $a_{i}(g) = $ |tree[(gN) + i]|  , the level of the immediate
ancestor of the individual on level |i| at time |g| } */

  assert( i < POP_SIZE);
  assert( g < (tree.size() * POP_SIZE) );
  assert( ((g*POP_SIZE) + i) < tree.size() );
  return ( tree[ (g*POP_SIZE) + i]) ;  
}



@*1 {\bf check if sample tree is complete}.
\label{sec:checkifcomplete}


check if a sample tree is complete, if we find a common ancestor of
all the sampled leaves in the given population ancestry; when in
generation  $g$  the goal is
to carry out the steps
\begin{displaymath}
\begin{split}
&  (i,\ldots,i, j, \ldots, j, k, \ldots, k) \\
\to &  (a_{i}(g), \ldots, a_{i}(g), a_{j}(g), \ldots, a_{j}(g), a_{k}(g), \ldots, a_{k}(g)) \\
\to & (a_{i}(g), a_{j}(g), a_{k}(g)) \\
g & \leftarrow g - 1 
\end{split} 
\end{displaymath}
where  $i = a_{i}(g+1)$,  $j = a_{j}(g+1)$,  $k = a_{k}(g+1)$




@<is sample tree complete?@>=@#
static bool allcoalesced(  const std::vector<std::pair< std::size_t, std::size_t>>& sample, const std::vector< std::size_t>& tre )
{

 @#

  assert(POP_SIZE <=  tre.size() );
  assert( (tre.size() % POP_SIZE) < 1 );
  int g =  static_cast<int>(tre.size() / POP_SIZE)  - 1 ;
  assert( sample.size() == SAMPLE_SIZE);

  
  std::vector< std::size_t> a(SAMPLE_SIZE, 0);
    
  /* \newline \aths{ pair is (size of block, level of block)} */
  /* \newline \aths{ copy ancestor indexes into |a|; then work with
  |a| to check if tree is complete} */
  
  std::transform( sample.begin(), sample.end(), a.begin(), [](const auto &x){return x.second;});
 
  while( (a.size() > 1) && (g > -1)){
  /* \aths{|getagi| \S~\ref{sec:getagi}; $(i,j,k) \to (a_{i}(g),
  a_{j}(g), a_{k}(g))$ } */
  std::transform(a.begin(), a.end(), a.begin(), [&g, &tre](const auto &i){return getagi(g, i, tre);});

  /* \newline \aths{remove duplicates from |a|; for example assuming
  $i < j < k$ then  $(k,i,k,j, \ldots , i) \to (i,j,k)$ } */
    std::sort( a.begin(), a.end());
    a.erase( std::unique( a.begin(), a.end()), a.end() );
    
    --g ;
  }

@#

 /* \newline  \aths{return True if the tree is complete} */
  return  (a.size() < 2) ;
}



@*1 {\bf sort in descending order}.
\label{sec:bigger}


sort the sample in descending order according to size of block; each
element of sample is (size of block, level  of block)


@<bigger@>=@#
static bool comp( const std::pair< std::size_t, std::size_t> &a, const std::pair< std::size_t, std::size_t> &b )
{
  return a.first > b.first ;
}



@*1 {\bf update sample}.
\label{sec:updatesample}



here we update the sample over one generation and merge  blocks with
same immediate ancestor;
suppose $a_{1}^{\prime}(g),\ldots, a_{k}^{\prime}(g)$ are the unique
immediate ancestors 

\begin{displaymath}
\begin{split}
\text{look up immediate ancestors} & \set{ \ldots,  \svigi{s,i }, \ldots } \to \set{ \ldots, \svigi{s, a_{i}(g)}, \ldots } \\
\text{sum block sizes of blocks with same ancestor} & \set{\ldots, \svigi{ \sum\nolimits_{j}\one{ a_{j}(g) = a_{\ell}^{\prime}(g) }s_{j}, a_{\ell}^{\prime}(g)}, \ldots } \\
\end{split}
\end{displaymath}
The goal here is to end up with the sample configuration
\begin{displaymath}
 \set{ \svigi{ \sum\nolimits_{j}\one{ a_{j}(g) = a_{1}^{\prime}(g) }s_{j}, a_{1}^{\prime}(g)},    \ldots, \svigi{ \sum\nolimits_{j}\one{ a_{j}(g) = a_{k}^{\prime}(g) }s_{j}, a_{k}^{\prime}(g)} }
\end{displaymath}


@<merge blocks and update sample@>=@#
static void  updatesample( const std::size_t& g, std::vector<std::pair< std::size_t, std::size_t>>& sample, const std::vector<std::size_t>& tree )
{

  /* \newline \aths{ recall each sample `block'  is (size of block, level of block)} */
   std::size_t s {} ;
   const std::size_t e = sample.size();
   std::size_t k {} ;
  std::size_t t {} ;
   for( std::size_t i = 0 ; i < e; ++i){
    s = 0 ;
    for( std::size_t j = i ; j < e ; ++j){
     /* \newline \aths{ sum size of blocks with same immediate
     ancestor;     |getagi| \S~\ref{sec:getagi}}  */
      s += ( getagi(g, sample[i].second ,tree) == getagi(g, sample[j].second, tree) ? sample[j].first : 0) ;
      /* \newline \aths{ if block has already been merged set size of block to zero} */
      sample[j].first = ( getagi(g,  sample[i].second, tree) == getagi(g, sample[j].second, tree) ? 0 : sample[j].first); } 
    if( s > 0 ){
      /* \newline \aths{ record only a new merged block or a current  block not merging} */
      ++k ;
      /* \newline \aths{record the continuing block with the block size} */
      /* \newline \aths{ and the level of the immediate ancestor of the block} */
      sample.push_back( std::make_pair( s, getagi( g, sample[i].second, tree) ) ) ; }
      t += s ;
      }

 /* \newline \aths{push blocks with size 0 to the back;  |comp| \S~\ref{sec:bigger}} */
   std::sort( sample.begin(), sample.end(), comp);
   /* \newline \aths{ |k| counts the new/continuing blocks} */
   sample.resize(k);
   assert( sample.size() == k); 
   assert( t == SAMPLE_SIZE);

   assert( sample.size() > 1 ? (sample.back().first < SAMPLE_SIZE) :   (sample.back().first == SAMPLE_SIZE) ); 
}



@*1 {\bf update $\ell_{i}^{N}(n)$}.
\label{sec:updatelengths}


the goal here is to update the lengths $\ell_{i}^{N}(n)$; if $s_{1},
\ldots, s_{m}$ are the current block sizes then
\begin{displaymath}
\ell_{s} \leftarrow \ell_{s} + 1
 \end{displaymath}
for $s = s_{1},\ldots, s_{m}$


@<update $\ell_{i}^{N}(n)$@>=@#
static void updatevbi( const std::vector< std::pair< std::size_t, std::size_t>>& sample, std::vector<std::size_t>& vbi )
{
  /* \newline \aths{ pair is (size of block, level of block)} */
  for( std::size_t i = 0; i < sample.size(); ++i){
      assert( sample[i].first > 0);
    assert( sample[i].first < SAMPLE_SIZE ); 
    ++vbi[0] ; 
    ++vbi[ sample[i].first] ; }
}



@*1 {\bf update $r_{i}^{N}(n)$}.
\label{sec:updateri}

given the lengths $\ell_{1}^{N}(n), \ldots, \ell_{n-1}^{N}(n)$  of one
complete tree  update the
approximation of $\EE{\widetilde R_{i}^{N}(n)}$ 

\begin{displaymath}
r_{i}^{N}(n) \leftarrow r_{i}^{N}(n) +   \frac{\ell_{i}^{N}(n)}{ \ell_{1}^{N}(n) + \cdots + \ell_{n-1}^{N}(n) }
\end{displaymath}
and it must hold that  $\sum_{j} \ell_{j}^{N}(n) \ge n $ 

@<update $r_{i}^{N}(n)$@>=@#
static void updaterestimate( const std::vector< std::size_t>& b, std::vector<double>& r,  const std::size_t SAMPLE_SIZE )
{

@#

  assert( b[0] >= SAMPLE_SIZE);
  for( std::size_t i = 1; i < SAMPLE_SIZE; ++i){
    r[i] += static_cast<double>( b[i])/static_cast<double>(b[0]) ; }
}



@*1 {\bf one experiment}.
\label{sec:onexperiment}


the goal here is to generate one complete sample tree and read the
branch lengths off the tree and update the approximation of
$\EE{R_{i}^{N}(n)}$; recall Figure~\ref{fig:illsampling}


@<generate one experiment@>=@#
static void onexperiment( std::vector<double>& vr, const
std::vector<double>& vcmf_alpha,  const std::vector<double>& vcmf_kappa )
{

@#

   std::vector<std::pair< std::size_t, std::size_t>>  sample {} ;
  std::vector< std::size_t>  tre (POP_SIZE,0);
  std::iota( tre.begin(), tre.end(), 0);
  std::vector<std::size_t> currentbi (SAMPLE_SIZE, 0);
  
  std::size_t g {} ;

  /* \newline \aths{sample levels for the sample  \S~\ref{sec:randomsample}} */
  randomsample(sample );
  /* \newline \aths{\S~\ref{sec:checkifcomplete}} */
  while( !allcoalesced( sample, tre) ){
  /* \newline \aths{add to the population ancestry \S \ref{sec:addgeneration}} */
    MODEL ? addgeneration_modelone(tre, vcmf_alpha, vcmf_kappa) : addgeneration_modelnull(tre, (gsl_rng_uniform(rngtype) < EPSILON ? vcmf_alpha : vcmf_kappa)) ;
    randomsample(sample);
    ++g ; }


/* \newline \aths{we have a complete sample tree so read the branch
lengths $\ell_{i}^{N}(n)$ } */
  while( sample.back().first < SAMPLE_SIZE){
  /* \newline \aths{ update $\ell_{i}^{N}(n)$   \S~\ref{sec:updatelengths} } */
    updatevbi( sample, currentbi);
    /* \newline \aths{ merge blocks and update sample configuration  \S~\ref{sec:updatesample}} */
    updatesample( g, sample, tre);
    --g ; }


 /* \newline \aths{given lengths $\ell_{i}^{N}(n)$ update
 approximation of 
  $\EE{\widetilde R_{i}^{N}(n)}  $ \S~\ref{sec:updateri}} */
  updaterestimate( currentbi, vr, SAMPLE_SIZE );
}



@*1 {\bf approximate $\EE{\widetilde  R_{i}^{N}(n)}$}. 
\label{sec:approximate}


 approximate $\EE{\widetilde  R_{i}^{N}(n)}$  by sampling complete
 sample trees and reading branch lengths and averaging;   recall Figure~\ref{fig:illancestries}


@<go ahead -- approximate $\EE{\widetilde R_i^N(n)}$@>=@#
static void approximate_eri( )
{

@#


  std::vector<double> vcmf_smalla (CUTOFF + 1 ) ;
  std::vector<double> vcmf_biga (CUTOFF + 1 ) ;
  std::vector< double> vR (SAMPLE_SIZE);

  /* \newline \aths{generate a cumulative mass function for sampling
  bounded numbers of potential offspring using \eqref{eq:1};    \S~\ref{sec:cmf} } */
  generate_cmfs( vcmf_smalla, vcmf_biga); 
  int r = EXPERIMENTS + 1 ;
  while( --r > 0){
  /* \newline \aths{\S~\ref{sec:onexperiment}} */
    onexperiment(vR, vcmf_smalla, vcmf_biga); }


  /* \newline \aths{return $r_{i}^{N}(n)$ summed over the experiments} */ 
  for( const auto & z: vR){ std::cout << z << '\n';}
}





@*1 {\bf the main module}.
\label{sec:main}

@C

/* \newline \aths{\S~\ref{sec:includes}} */
@<included libraries@>@#
/* \newline \aths{\S~\ref{sec:rngs}} */
@<random number generators@>@#
/* \newline \aths{\S~\ref{sec:mass}} */
@<mass function@>@#
/* \newline \aths{\S~\ref{sec:cmf}} */
@<cmf for number of offspring@>@#
/* \newline \aths{\S~\ref{sec:boundedX}} */
@<one bounded number@>@# 
/* \newline \aths{\S~\ref{sec:unboundedX}} */
@<an unbounded $X$@>@#
/* \newline \aths{\S~\ref{sec:addgeneration}} */
@<one generation@>@#
/* \newline \aths{\S~\ref{sec:onegenerationM1} } */
@<add a generation for model 1@>@#
/* \newline \aths{\S~\ref{sec:randomsample}} */
@<random sample@>@#
/* \newline \aths{\S~\ref{sec:getagi}} */
@<get immediate ancestor@>@#
/* \newline \aths{\S~\ref{sec:checkifcomplete}} */
@<is sample tree complete?@>@#
/* \newline \aths{\S~\ref{sec:bigger}} */
@<bigger@>@#
/* \newline \aths{\S~\ref{sec:updatesample}} */
@<merge blocks and update sample@>@#
/* \newline \aths{\S~\ref{sec:updatelengths} } */
@<update $\ell_{i}^{N}(n)$@>@#
/* \newline \aths{\S~\ref{sec:updateri}} */
@<update $r_{i}^{N}(n)$@>@#
/* \newline \aths{\S~\ref{sec:onexperiment}} */
@<generate one experiment@>@# 
/* \newline \aths{\S~\ref{sec:approximate}} */
@<go ahead -- approximate $\EE{\widetilde R_i^N(n)}$@>@#


int main(int argc, const char * argv[])
{

/* \newline \aths{input random seed  \S~\ref{sec:rngs}} */
setup_rng( static_cast<std::size_t>(atoi(argv[1])) );

/* \newline \aths{\S~\ref{sec:approximate}} */
  approximate_eri( );


 gsl_rng_free(rngtype);

return 0 ;
} 





@* {\bf conclusion and bibliography}.
\label{sec:concl}


Gene genealogies are generated as the population evolves forward in
time (Figure~\ref{fig:forw}).    We 
evolve a haploid panmictic population of constant size forward in time
and record the ancestry (Figure~\ref{fig:anc}).  We read the branch
lengths off independent fixed complete sample  trees to approximate
$\EE{\widetilde R_{i}^{N}(n)}$.  In Figure~\ref{fig:res} are  examples of
approximations of $\EE{\widetilde R_{i}^{N}(n)}$.


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
       %%
       \end{axis}
       \end{tikzpicture}
       \caption{\textcolor{violet}{ \it  An example approximation of
       $\EE{\widetilde  R_{i}(n)}$
 for the  given  parameter values   and graphed as logits against
 $\log(i/n) - \log(1 - i/n)$ for $i = 1,2,\ldots, n-1$ where $n$ is
 sample size }}
       \label{fig:res}
       \end{SCfigure}




Here we consider the simplest case of a single haploid panmictic
population of constant size. The population may evolve according to
sweepstakes reproduction (heavy-tailed offspring number
distribution). We model sweepstakes reproduction through the law
\eqref{eq:1} on the number of potential offspring produced
independently by the current individuals;  \eqref{eq:1} is
an extension of the one in \cite{schweinsberg03}.  There is evidence
that  predictions of genetic diversity in populations evolving
according to sweepstakes reproduction may  depend on how one views
gene genealogies \cite{D2024}



\bibliographystyle{plain}
%%\bibliography{refs}
%% use with biber
%% xelatex, biber %, xelatex
%% \textcite{}  \parencite{}
%% \printbibliography
%% kaffi

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
