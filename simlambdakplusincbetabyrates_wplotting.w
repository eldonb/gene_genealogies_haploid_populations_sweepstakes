%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Copyright (C) 2025 Bjarki Eldon
%$ NAFN=simlambdakplusincbetabyrates_wplotting
%$ ctangle $NAFN
%$ g++ -std=c++26 -m64 -march=native -O3 -x c++ $NAFN.c -lm -lgsl -lgslcblas
%$ rm -f tmp_runs
%$ for i in $(seq 5); do echo "./a.out "  $(shuf -i 39393-282929191 -n1) "> resout"$i >> tmp_runs; done
%$ parallel --gnu -j5 :::: ./tmp_runs
%$ paste resout* -d, | sed '1d' | awk -F, '{M=1e6;s=0; for (i=1;i<=NF;i++) {s+=$i} print log(s/M) - log(1-(s/M))}' > logitresout
%$ seq 99 | awk '{S=100;print log($1/S) - log(1 - ($1/S))}' > nlogits
%$ paste -d',' nlogits logitresout > forplottingfile1
%$ g++ --version | head -n1 > innleggcpp
%$ gsl-config --version > innleggsl
%$ xelatex --version | head -n1  > innleggxelatex
%$ lualatex --version | head -n2 > innlegglualatex
%$ cweave  --version | head -n1 > innleggcweave
%$ uname  --kernel-release -o  > innleggop
%$ ctangle  --version | head -n1  > innleggctangle
%$ emacs --version | head -n1 > innleggemacs
%$ parallel --version | head -n1 > innleggparallel
%$ awk --version | head -n1 > innleggawk
%$ spix --version > innleggspix
%$ NAFN=simlambdakplusincbetabyrates_wplotting
%$ cweave $NAFN
%$ tail -n4 $NAFN.tex > endi
%$ for i in $(seq 5); do $(sed -i '$d' $NAFN.tex) ; done
%$ cat endi >> $NAFN.tex
%$ lualatex $NAFN
%$ lualatex $NAFN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\documentclass[a4paper,10pt]{cweb}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
%% for xelatex
%\usepackage{xunicode}
%\usepackage{fontspec}
%\usepackage{xltxtra}
%%\usepackage[lf]{Baskervaldx}
%%\usepackage[bigdelims,vvarbb]{newtxmath}
\usepackage{amsfonts, amsmath}
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
\newcommand{\EE}[1]{\ensuremath{\mathds{E}\left[ #1 \right]}}%
\newcommand{\im}{\ensuremath{\imath} }%
\newcommand{\jm}{\ensuremath{\jmath} }%
\newcommand{\prob}[1]{\ensuremath{\mathds{P}\left( #1 \right) } }%
\newcommand{\h}[1]{\ensuremath{\uptheta_{ #1 } } }%
\newcommand{\VV}[1]{\ensuremath{ \mathbb{V}\left( #1 \right)}}%
\newcommand{\hp}{\ensuremath{\theta_1}}
\newcommand{\hs}{\ensuremath{\theta_2}}
\newcommand{\D}{\ensuremath{\mathbb{D}}}
\newcommand{\F}{\ensuremath{\mathbb{F}} }
\newcommand{\bt}[1]{\textcolor{blue}{\tt #1}}
\newcommand{\set}[1]{\ensuremath{\left\{ #1 \right\}}}
\newcommand{\aths}[1]{\textcolor{violet}{\textbf{\small #1}}}
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
\title{Gene genealogies in  haploid populations evolving according to sweepstakes reproduction\\  --- estimating $\EE{R_{i}(n)}$ for the $\delta_{0}$-Beta($\gamma, 2-\alpha, \alpha$)-coalescent}
\author{Bjarki Eldon\footnote{\href{mailto:beldon11@@gmail.com}{beldon11@@gmail.com}} 
%%  Deutsche Forschungsgemeinschaft (DFG) - Projektnummer 273887127 
%% „funded by the Deutsche Forschungsgemein-schaft (DFG, German Research Foundation) –Projektnummer(n)“.
%% 
%through DFG SPP 1819: Rapid Evolutionary Adaptation grant STE 325/17
%to Wolfgang Stephan; acknowledge funding by the Icelandic Centre of
%Research (Rann\'is) through an Icelandic Research Fund
%(Ranns\'oknasj\'o{\dh}ur) Grant of Excellence no.\ 185151-051 to Einar
%\'Arnason, Katr\'in Halld\'orsd\'ottir, Alison M.\ Etheridge, WS, and
%BE; Start-up module grants through SPP 1819 with
%Jere Koskela and Maite Wilke-Berenguer, and with Iulia Dahmer. 
\footnote{compiled @@ {\DTMcurrenttime} on  {\today}  \\ \input{innleggop}\\   \input{innleggctangle}  \\  \input{innleggcpp} \\ GSL \input{innleggsl}\\  \input{innleggcweave} \\  \input{innlegglualatex}  \\   \input{innleggparallel}  \\ \input{innleggawk}  \\ \input{innleggspix} \\   \input{innleggemacs}}
\orcidlink{0000-0001-9354-2391}}
\begin{document}
\maketitle
\renewcommand{\abstractname}{\vspace{-\baselineskip}}

%%\rule{\textwidth}{.8pt}

\begin{abstract}
  Let $\{ \xi^{n}(t) : t \ge 0 \}$ be the $\delta_{0}$-Beta($\gamma,
2-\alpha, \alpha$)-coalescent, a Markov chain tracking the random
relations of sampled gene copies, $\# A$ is the cardinality of a given
set $A$, $L_{i}(n) \equiv \int_{0}^{\tau(n) } \# \left\{ \xi \in
\xi^{n}(t) : \#\xi = i \right\}dt $ and $L(n) \equiv
\int_{0}^{\tau(n)} \# \xi^{n}(t)dt$ and $\tau(n) \equiv \inf \left\{ t
\ge 0 : \# \xi^{n}(t) = 1 \right\} $ for $i \in \{1, 2, \ldots,
n-1\}$; $R_{i}(n) \equiv L_{i}(n)/L(n)$  for
$i=1,2,\ldots, n-1$.  Then $L_{i}(n)$ is interpreted as the random
total length of branches supporting $i \in \{1, 2, \ldots, n-1\}$
leaves, with the length measured in coalescent time units,  and $n$ sample size.
We then have $L(n) = L_{1}(n) + \cdots + L_{n-1}(n)$.  With this
C++ code one estimates $\EE{R_{i}(n)}$ predicted by the
$\delta_{0}$-Beta($\gamma, 2-\alpha, \alpha$)-coalescent for $0 <
\gamma \le 1$ and $0 < \alpha < 2$. This family of multiple-merger
coalescents can be shown to describe the gene genealogies of a sample
of $n$ gene copies from a haploid panmictic population of constant
size evolving in a random environment (Definitions \ref{def:haplrandomalpha} and \ref{def:alpha-random-one})  according to sweepstakes
reproduction \cite{eldon24:_coales}.
\end{abstract}

\tableofcontents



@* {\bf Copyright}. 


Copyright {\textcopyright} {\the\year}  Bjarki Eldon \newline



This document and any source code it contains  is distributed under the terms of the GNU General Public Licence (version $\ge 3$).  You
should have received a copy of the licence along with this file (see file COPYING).  


    The source codes  described in this document  are  free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This document and the code it contains   is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this file (see COPYING).  If not, see \url{http://www.gnu.org/licenses/}.


@* {\bf Compilation,  output and execution}. 
\label{compile}

 This CWEB
      \citep{knuth1994cweb} document (the {\tt .w} file) can be
      compiled with {\tt cweave} to generate a {\tt .tex} file, and
      with {\tt ctangle} to generate a {\tt .c} \citep{kernighan1988c}
      file.


Compiles on Linux debian 6.12.6-amd64 with  CTANGLE 4.11 (TeX Live 2025/dev/Debian) and  g++ 14.2 and  GSL 2.8

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
\tab\quad\quad\quad\quad        g++ -Wall -Wextra -pedantic -std=c++26 -O3 -march=native -m64 -x c++ iguana.c -lm -lgsl -lgslcblas \\
        
       
clean :  \\
\tab\quad\quad\quad\quad        rm -vf iguana.c iguana.tex \\
}



Use {\tt valgrind} to check for memory leaks:

{\tt valgrind -v --leak-check=full --show-leak-kinds=all --vgdb=full --leak-resolution=high --num-callers=50 <program call>}




Use {\tt cppcheck} to check the code

{\tt cppchek  ---enable=all ----language=c++ <prefix>.c}

To generate estimates on a computer with several CPUs it may be
convenient to put  in a text file ({\tt simfile}):


{\tt ./a.out \$(shuf -i 484433-83230401 -n1) > resout<i>}


for  $i = 1,\ldots, y$ and use  {\tt
parallel}\cite{tange11:_gnu_paral}

{\tt parallel ---gnu -jy :::: ./simfile}





@* {\bf intro}.  
\label{sec:intro}

Suppose the population evolves as in Definition~\ref{hschwpop}.  

\begin{defn}[Evolution of a haploid panmictic population]
\label{hschwpop}%
Consider a haploid, i.e.\ each individual carries exactly one gene
copy, panmictic population of constant size $N$ evolving in discrete
(non-overlapping) generations.  In any given generation each
individual independently produces a random number of potential
offspring according to some given law.  If the total number of
potential offspring produced in this way is at least $N$, then $N$ of
them sampled uniformly at random without replacement survive to
maturity and form a new set of reproducing individuals, and the
remaining potential offspring perish (before reaching maturity).
Otherwise we will assume an unchanged population over the generation
(all the potential offspring perish before reaching maturity).
\end{defn}%


Let $\left\{\xi^{N,n}(t); t \ge 0 \right\}$ denote the Markov sequence describing
the random ancestral relations of $n$ gene copies (leaves) sampled
from a population of size $N$ evolving according to  Def~\ref{hschwpop}.   
If  the random number of potential offspring ($X$)   of an arbitrary
individual is  distributed according to ($C, \alpha> 0$ fixed)
\begin{equation}
\label{eq:4}
\lim_{x \to \infty}Cx^{\alpha}\prob{X \ge x}  = 1
\end{equation}
\cite[Eq~11]{schweinsberg03} then $\left\{\xi^{N,n}(\lfloor t/c_{N} \rfloor
); t \ge 0 \right\}$ converges as $N\to \infty$ in the Skorokhod topology to
the Beta$(2-\alpha,\alpha)$-coalescent provided $1 \le \alpha < 2$
\cite[Thm 4c]{schweinsberg03}.


Suppose $X$ is distributed according to ($[n] \equiv \set{1,2,\ldots, n}$ for $n\in \mathds N \equiv  \set{1,2,\ldots}$)
\begin{equation}
\label{eq:PXiJ}
g(k) \left( \frac{1}{k^a} - \frac{1}{(1+k)^a} \right)  \leq
\prob{X = k} \leq f(k) \left( \frac{1}{k^a} -
\frac{1}{(1+k)^a} \right), \quad k  \in
[\uN],
\end{equation}
where the functions $g,f$ are independent of $N$ and so that
$\EE{X} > 1$ and $\prob{X \le \uN} = 1$.  Write
$X\vartriangleright \mathds L(a,\zeta(N))$ when $X$ is distributed
according to  \eqref{eq:PXiJ} with $a$ and $\zeta(N)$ as given each time.     Furthermore, the population evolves according to
Definition~\ref{def:haplrandomalpha} or
Definition~\ref{def:alpha-random-one}.
\begin{defn}[The random environment of type $A$]
\label{def:haplrandomalpha}
%%  by Eq~\eqref{eq:PXiJ}.
Suppose a population evolves according to Def.\ \ref{hschwpop} with
the number of potential offspring produced by each individual
distributed according to \eqref{eq:4} or \eqref{eq:PXiJ}.  Fix
$0 < \alpha < 2 \le \at $.  Write $E $ for the event
$\left\{ X_i \vartriangleright \mathds L(\alpha, \uN) \text{ for all } i \in [N] \right\}$ and $E^{\sf c}$
when $\at$ replaces $\alpha$ in $ E$ for some given $\uN$.  Suppose,
with $(\varepsilon_{N})_{N}$ a positive sequence with
$0 < \varepsilon_{N} < 1$ and $\varepsilon_{N}\to 0$,
\begin{displaymath}
\prob{ E}  = \varepsilon_N, \quad \prob{ E^{\sf c}}  = 1 - \varepsilon_N
\end{displaymath}
%% Further, it will be assumed that Assumption~\ref{a:mN1plus} holds so that $\EE{X_{1}} > 1$ for all $N$. 
\end{defn}%



\begin{defn}[The random  environment of type $B$]
\label{def:alpha-random-one}
Suppose a population evolves according to Definition~\ref{hschwpop} with
the law for the number of potential offspring given by
\eqref{eq:PXiJ}.  Fix $0 < \alpha < 2 \le \at$.  Write $ E_{1}$ for
the event
\begin{displaymath}
\set{\text{there exists exactly one   $i \in [N]$ where   $X_i \vartriangleright \mathds L(\alpha, \uN)$, and $X_{j} \vartriangleright \mathds L(\at,
\uN) $  all  $j \in [N]\setminus \set{i}$}}
\end{displaymath}
with $i$ picked
uniformly at random, and $ E _{1}^{\sf c}$ when $\at$ replaces
$\alpha$ in $ E_{1}$ so that $ E_{1}^{\sf c} = E^{\sf c}$ from
Definition~\ref{def:haplrandomalpha}. Suppose
%%write $E_{2}$ for the event  $\{ X_i^{N} \sim L(\at, \uN)\, \forall\, i \in [N] \}$.
\begin{displaymath}
\begin{split}
\prob{ E_{1} } & = \overline \varepsilon_N, \quad \prob{ E_{1}^{\sf c}}  = 1 - \overline \varepsilon_N
\end{split}
\end{displaymath}
where $(\overline \varepsilon_{N})_{N}$ is a positive sequence with
$0 < \overline \varepsilon_{N} < 1$ for all $N$ and it may hold that
$\overline \varepsilon_{N}\to 0$ as $N\to \infty$.
\end{defn}
Then $\{\xi^{N,n}(\lfloor t/c_{N} \rfloor ); t \ge 0\}$ converges as
$N\to \infty$ in the Skorohod topology to the
$\delta_{0}$-Beta$(\gamma,2-\alpha, \alpha)$-coalescent provided
$\liminf_{N\to \infty} \uN/N > 0$.  One may identify conditions on
$\varepsilon_{N}$ and $\overline \varepsilon_{N}$ such that the models
described in Definition~\ref{hschwpop} and
Definition~\ref{def:alpha-random-one} and
Definition~\ref{def:haplrandomalpha} are in the domain of attraction of the
$\delta_{0}$-Beta$(\gamma,2-\alpha,\alpha)$-coalescent
\cite{eldon24:_coales}.


 Let  $\{ \xi^{n}(t) : t \ge 0 \}$ be  a given
coalescent, a Markov chain tracking the random relations of sampled
gene copies, $\# A$  the cardinality of a given set $A$,   
$L_{i}^{N}(n) \equiv \int_{0}^{\tau(n) } \# \left\{ \xi \in \xi^{n}(t)
: \#\xi = i \right\}dt $ and
$\tau(n) \equiv \inf \left\{ t \ge 0 : \# \xi^{n}(t) = 1 \right\} $
for $i \in \{1, 2, \ldots, n-1\}$; 
$R_{i}(n) \equiv L_{i}(n)/\sum_{j}L_{j}(n) $ for
$i=1,2,\ldots, n-1$.  Then $L_{i}^{N}(n)$ is interpreted as the random
total length of branches supporting $i \in \{1, 2, \ldots, n\}$
leaves, with the length measured in generations, and $n$ sample size.
We then have $L(n) = L_{1}(n) + \cdots + L_{n-1}(n)$. 



   We estimate 
\begin{equation}
\label{eq:varphi}
\varrho_{i}(n) :=  \EE{ \frac{L_{i}(n) }{ L(n) } }, \quad i = 1,2,\ldots, n-1
\end{equation}
for the $\delta_{0}$-Beta($\gamma, 2-\alpha, \alpha$)-coalescent with
transition rates
\begin{subequations}
\begin{align}
\label{eq:ratesarandall}
\lambda_{n,k} & =  \one{k=2}\binom{n}{k}\frac{ C_{\kappa} }{C_{\alpha,\gamma}} +   \binom{n}{k} \frac{ \alpha c_{\alpha}^{\prime}  }{ C_{\alpha,\gamma} m_{\infty}^{\alpha} } B(\gamma, k-\alpha,n - k+ \alpha)  \\
\label{eq:27}
%% \mathfrak{w}
\gamma & =  \one{\frac{\uN}{N} \to K } \frac{K}{m_{\infty} + K} + \one{ \frac{\uN}{N} \to \infty }  \\
\label{eq:23}
C_{\alpha,\gamma} &  =   C_{\kappa}   +  \frac{ \alpha c }{m_{\infty}^{\alpha}} \int_{\one{\uN/N \to K}\frac{m_{\infty} }{K + m_{\infty} }  }^{1} (1-y)^{1-\alpha}y^{\alpha - 1}dy  =  C_{\kappa} +   \frac{ \alpha c }{m_{\infty}^{\alpha}} B(\gamma, 2-\alpha,\alpha)   \\
%%  C_{\alpha,\gamma} & =    C_{\kappa}   +  \frac{ \alpha c }{m_{\infty}^{\alpha}} B(\gamma, 2-\alpha,\alpha) \\
\label{eq:1}
C_{\kappa} & =   \one{\kappa = 2} m_{\infty}^{-2}    +  \one{\kappa > 2 } m_{\infty}^{-2} c_{\kappa}
\end{align}
\end{subequations}
The range of $\alpha$ in Def~\ref{def:haplrandomalpha} and
Def~\ref{def:alpha-random-one} intersects at 1 so one must decide
which of the two definitions holds when $\alpha = 1$.  Let where
$0 < c^{\prime} < 1$ and $c > 0$. Write $c_{\alpha}^{\prime} = c$ if
Def~\ref{def:haplrandomalpha} holds when $\alpha = 1$, and
$c_{\alpha}^{\prime} = \one{\alpha = 1}c^{\prime} + \one{\alpha \in
(0,2)\setminus \{1\} }c $ when Def~\ref{def:alpha-random-one} holds
when $\alpha = 1$ (and $\kappa > 2$).  We approximate $m_{\infty}$
with $m_{\infty} \approx (2 + (1 + 2^{1 - \kappa })/(\kappa -
1))/2$. For $\kappa > 2$ we can approximate $c_{\kappa} \approx  (c_{\kappa}^{\prime} + c_{\kappa}^{\prime\prime})/2$ where
\begin{displaymath}
c_{\kappa}^{\prime} =   \frac{\kappa 2^{2-\kappa} }{\kappa - 2} - \frac{3\kappa 2^{1-\kappa} }{\kappa - 1 } + 2^{-\kappa}, \quad  c_{\kappa}^{\prime\prime} =  \frac{\kappa 2^{2-\kappa} }{\kappa - 2} - \frac{\kappa 2^{1-\kappa} }{\kappa - 1 }
\end{displaymath}
When $\zeta(N)/N\to 0$ we have, as one would expect, $\gamma = 0$ and
so $B(\gamma,2-\alpha,\alpha) = 0$ and
$C_{\alpha,\gamma} = C_{\kappa}$ and thus
$\lambda_{n,k} = \one{\kappa = 2}\binom{n}{k}$, the transition rate
of  the Kingman-coalescent.


In \S~\ref{sec:code} we summarize the algorithm.  Our goal is simply
to arrive at a working code that is correct. The code is partitioned
into the functions described in
\S~\ref{sec:includes}--\S~\ref{sec:main}; we conclude in
\S~\ref{sec:conclusions}.  Comments within the code are shown in
\aths{this font and colour.}


@* {\bf code}. 
\label{sec:code}


We use  the GSL library. It suffices to keep track of
the current block sizes $(b_{1}, \ldots, b_{m})$ where $b_{j} \ge 1$
and $\sum_{j}b_{j} = n$. Let $n$ denote the sample size and $m$ the
current number of blocks. Let $\ell_{i}(n)$ for $i = 1,2,\ldots, n-1$
denote realised branch lengths.

\begin{enumerate}
\item $\left( r_{1}(n),\ldots, r_{n-1}(n) \right) \leftarrow (0,\ldots, 0)$
\item for each of $M$ experiments \S~\ref{sec:estimate} :
\begin{enumerate}
\item initialize block sizes  $(b_{1},\ldots, b_{n}) \leftarrow (1,\ldots, 1)$
\item initialise branch lengths  $(\ell_{i}(n), \ldots, \ell_{n-1}(n)) \leftarrow (0, \ldots, 0)$
\item set the current number of blocks to sample size  $m \leftarrow n$
\item {\bf while} $m > 1$:
\begin{enumerate}
\item sample exponential time $t$ with rate \eqref{eq:ratesarandall}  $\lambda_{m} = \lambda_{m,2} + \cdots + \lambda_{m,m}$  \S~\ref{sec:genealogy}
\item update branch lengths $\ell_{b_{j}}(n) \leftarrow t +  \ell_{b_{j}}(n) $ for $j = 1,2,\ldots, m$  \S~\ref{sec:updatecurrentL}
\item sample merger size $k = \min\{j : u \le   \sum_{i=2}^{j}\lambda_{i,m}/\lambda_{m} \} $ \S~\ref{sec:samplemergersize}
\item merge blocks and record size $b_{\sigma(1)}  + \cdots  + b_{\sigma(k)}$ of continuing block  \S~\ref{sec:genealogy}
\item given merger size update current number of blocks  $m\leftarrow m - k + 1$
\end{enumerate}
\item given a realisation of branch lengths  update estimate of $\varrho_{i}(n)$ \S~\ref{sec:updateri}  :
\begin{displaymath}
r_{i}(n) \leftarrow r_{i}(n) + \frac{ \ell_{i}(n)}{ \sum_{j=1}^{n-1}\ell_{j}(n) }
\end{displaymath}
\end{enumerate}
\item return an estimate  $(1/M)r_{i}(n)$   of  $\varrho_{i}(n) = \EE{R_{i}(n)}$   for $i = 1,2,\ldots, n-1$
\end{enumerate}


@*1 {\bf the included libraries}.
\label{sec:includes}


the included C++ libraries 

@<includes@>=@#
#include <iostream>
#include <vector>
#include <random>
#include <functional>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <list>
#include <forward_list>
#include <cassert>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_cdf.h>


@*1 {\bf GSL random number generator}.
\label{sec:gslrng}


@<rngs@>=@#

std::random_device randomseed;
  /* \newline \aths{  Standard mersenne twister  random number engine seeded with |rng()|} */
  std::mt19937_64 rng(randomseed());
/* \newline \aths{  define a GSL random number generator} */
gsl_rng * rngtype;

static  void setup_rng(  const unsigned long s)
  {
    const gsl_rng_type *T ; 
    gsl_rng_env_setup(); 
    T = gsl_rng_default ;
    rngtype = gsl_rng_alloc(T);
    gsl_rng_set( rngtype,  s) ;
  }



@*1 {\bf the beta function}. 
\label{sec:betafunc}

return the logarithm of the (incomplete) beta function
$B(x,a,b) = \int_{0}^{x} t^{a-1}(1-t)^{b-1}dt$ for $0 < x \le 1$ and
$a,b > 0$ using the Gauss hypergeometric function
\begin{displaymath}
B(x,a,b) = \frac 1a x^{a}(1-x)^{b}F(a+b, 1, a+1, x)
\end{displaymath}

@<beta function@>=@#
static long double betafunc( const double& x, const double& a, const  double&  b )
{
  /* \newline  \aths{ the GSL incomplete beta function is normalised by the complete beta function}  */
  assert( x > 0);
  assert( a > 0);
  assert( b > 0);
  /* \newline \aths{ the standard way would be | gsl_sf_beta_inc( a, b, x) * gsl_sf_beta( a, b) |} */
  /* \newline \aths{ return  the logarithm  of the beta function  as $\log\Gamma(a) + \log\Gamma(b) - \log\Gamma(a+b)$ } */
  const long  double f = static_cast<long double>( (x < 1 ? gsl_sf_hyperg_2F1(a + b, 1, a+1,x) : 1) );
  
  assert( f > 0);
  /* \newline  \aths{ return $\one{x < 1}( \log f  +a\log x + b\log(1-x) - \log a) + \one{x = 1} (\log\Gamma(a) + \log\Gamma(b) - \log\Gamma(a+b))$} */
  return( x < 1 ?  (logl( f ) +  ( static_cast<long double>( (a*log(x)) + (b * log(1-x)) - log(a) ) ) ) : (lgammal( static_cast<long double>(a)) + lgammal( static_cast<long double>(b) ) - lgammal( static_cast<long double>(a + b) ) ) ); 
}



@*1 {\bf binomial constant}.
\label{SEC:binomialconstant}

compute the log of the  binomial constant $\binom{n}{k}$  

@<binomial@>=@#
static long double binom( const int& n, const int& k)
{

@#

assert( k <= n);
assert( k >= 0);
assert( n >= 0);
/* \newline  \aths{ $\log \binom{n}{k} =  \log \Gamma(n+1) - \log\Gamma(k+1) - \log\Gamma(n-k + 1)$} */

   return( lgammal( static_cast<long double>(n+1) ) - lgammal( static_cast<long double>(k+1) ) - lgammal( static_cast<long double>(n-k + 1) ) ) ;
}



@*1 {\bf total beta rate}. 
\label{sec:totalbetarate}

return the total beta rate  $ \binom{n}{k}B(\gamma, k-\alpha, n+\alpha - k)$. 

@<betarate@>=@#
static double lbetank( const int& n, const int& k, const double& a, const double& g)
{
/* \newline \aths{ return $\exp( \log \binom{n}{k} +  \log B(\gamma, k-\alpha, n-k + \alpha) )$} */
/* \newline \aths{ |binom|  \S~\ref{SEC:binomialconstant} and |betafunc|  \S~\ref{sec:betafunc}} */

  return( static_cast<double>( expl( binom(n, k) + betafunc( g, static_cast<double>(k) - a, static_cast<double>(n-k) + a) ) ) ) ;
}



@*1 {\bf the constant $C_{\alpha,\gamma}$ \eqref{eq:23}}. 
\label{sec:C}

compute the constant $C_{\alpha,\gamma}$  \eqref{eq:23}  where we
approximate
$m_{\infty} \approx (2 + (1 + 2^{1-\kappa})/(\kappa - 1))/2 $

@<constantC@>=@#
static double constantC( const double& a, const double& kappa, const double& c, const double& cprime, const double&  g)
{
  /* \newline \aths{ $a = \alpha$, $k = \kappa$, $c =c$, $cprime = c^{\prime}$} */
  /* \newline \aths{ $g = \gamma =  \one{\uN/N \to K} \tfrac{K}{m_{\infty} + K} + \one{\uN/N \to \infty}$  and $ \tfrac{m_{\infty}}{K + m_{\infty}} = 1 - \gamma$} */

  const double m = (2. + (1. + pow(2, 1-kappa))/(kappa-1))/2. ;

  /* \newline \aths{ \eqref{eq:1}; first taking $\kappa = 2$} */
  const double Ck = 1/pow(m, 2.) ;

/* \newline \aths{ $c_{\alpha}^{\prime} =  \one{\alpha = 1}c^{\prime} + \one{ \alpha \in (0,2)\setminus\{1\} }c $} */
  const double caprime =   (a < 1 ? c : (a > 1 ? c : cprime) );

  /* \newline \aths{ the integral in \eqref{eq:23}  $B(\gamma,2-\alpha,\alpha) =  \int_{\one{\uN/N \to K}\frac{m_{\infty} }{K + m_{\infty} }  }^{1} (1-y)^{1-\alpha}y^{\alpha - 1}dy $} */
 /* \newline \aths{ the incomplete beta function  is normalised by the complete beta function in GSL} */
const double tegur = (g < 1 ? (gsl_sf_beta_inc( 2 - a, a, g)*gsl_sf_beta(2-a,a)) : gsl_sf_beta(2-a,a) ) ;

  return (Ck +  (a*caprime*tegur/pow(m, a))) ;
}



@*1 {\bf $\lambda_{n,k}$  \eqref{eq:ratesarandall}}.
\label{sec:lambdank}

compute the transition rate $\lambda_{n,k}$ in   \eqref{eq:ratesarandall}

@<lambdank@>=@#
static  double lambdank( const int& n, const int& k, const double& a, const double& kappa, const double& c, const double& cprime, const double& gm)
{

@#

   assert( k > 1);
   assert(n > 1) ;
   assert( k <= n);
   /* \newline \aths{ we approximate $m_{\infty} \approx  \tfrac 12 \left(2 +  \tfrac{1 + 2^{1-\kappa}}{\kappa - 1} \right) $} */
    const double m = (2. + (1. + pow(2., 1.- kappa) )/(kappa-1.) ) / 2. ;
    /* \newline \aths{ |Ck| is $C_{\kappa}$   \eqref{eq:1}; first taking $\kappa = 2$} */
     const double Ck = 2/pow(m, 2.) ;
      const double caprime =   (a > 1 ?  c : (a < 1 ? c :  cprime) );

/*  \newline  \aths{  return  \eqref{eq:ratesarandall} using |lbetank|  \S~\ref{sec:totalbetarate} and  |constantC| \S~\ref{sec:C} }  */
  return( ((k < 3 ?  (Ck * ( static_cast<double>(k*(k-1)))/2.) : 0) + (a*caprime*lbetank(n, k, a, gm)/pow(m, a))) / constantC(a, kappa,c, cprime, gm) ); 
}


@*1 {\bf the total rate $\lambda_{n}$}. 
\label{sec:lambdan}

compute the total jump rate $\lambda_{n} = \sum_{k}\lambda_{n,k}$ with
$\lambda_{n,k}$ as in \eqref{eq:ratesarandall}


@<jump rate $\lambda_{n}$@>=@#
static void lambdan( const int& n, const double &ia, const double &ikappa, const double &ic, const double &icpr, const double &ig, std::vector<double>& l )
{

@#

   assert( n > 1); 
   for( int i = 2; i <= n; ++i){
   /* \newline \aths{  compute $\lambda_{i}$ for $i = 2,3, \ldots, n$} */
   l[i] = 0; 
   for( int j = 2; j <= i; ++j){
   /* \newline \aths{ compute $\lambda_{i}$} */
     /* \newline \aths{ |lambdank| \S~\ref{sec:lambdank}} */
     
   l[i] +=  lambdank( i, j, ia, ikappa, ic, icpr, ig); }
   assert( l[i] > 0);
   }
}





@*1 {\bf sample merger size}.
\label{sec:samplemergersize}

sample merger size $\min\{ j \in \{2,3, \ldots, m\}:  \sum_{k=2}^j\lambda_{n,k}/\sum_{i=2}^m \lambda_{i} \ge U \}$ for a given random uniform $U$

@<merger size@>=@#
static int samplemerger( const int& m, const double &ia, const double &ikappa, const double &ic, const double &icpr, const double &ig, const std::vector<double>& l)
{

 /* \newline \aths{ |m| is current number of blocks} */
 assert( m > 1);
  int j = 2 ;
  const double u = gsl_rng_uniform( rngtype);
  /* \newline \aths{ |lambdank| \S~\ref{sec:lambdank}} */ 
  double  s = lambdank(m, 2, ia, ikappa, ic, icpr, ig) ;
  assert( l[m] > 0);
  while(  u > (s/l[m]) ){
    ++j ;
    s += lambdank(m, j,  ia, ikappa, ic, icpr, ig) ; }

  /* \newline \aths{  return the size of the merger} */
  return j ;
}


@*1 {\bf update estimate of $\varrho_{i}(n)$ \eqref{eq:varphi}}.
\label{sec:updateri}

update estimate of $\varrho_{i}(n)$; i.e.\ given a realisation $r(n)$
of $(R_{1}(n), \ldots, R_{n-1}(n))$ update the estimate
$e_{i}(n) \leftarrow e_{i}(n) + r_{i}(n)$ where
$r_{i}(n) = \ell_{i}(n)/\sum_{j}\ell_{j}(n)$ 


@<update $\varrho_{i}(n)$@>=@#
static void updateri( const int& leaves, const std::vector<double>& tLi, std::vector<double>& eri )
{
/* \newline \aths{  |tLi| is the the current realisation $\ell_{i}(n)$   of  $L_{i}(n)$ } */
/* \newline \aths{ |tLi[0]| is the current total tree length $\sum_{i}\ell_{i}(n)$} */
  assert( tLi[0] > 0);

   for( int i = 1; i < leaves ; ++i){
   /* \newline \aths{ for checking return estimates of $\EE{L_i(n)} / \EE{L(n)}$; |tLi[0]| records estimates of $\EE{L(n)}$} */
   /* \newline \aths{ update with  | tLi[i]/tLi[0]| when comparing with $\EE{R_i^N(n)}$}  */
    eri[i] += (tLi[i] / tLi[0]) ; }
}


@*1 {\bf update current branch length}. 
\label{sec:updatecurrentL}


update current branch lengths $\ell_{i}(n)$

@<update current lengths@>=@#
static void updatells( const double &qi,   const std::vector<int>& tree, std::vector<double>& elli)
{

/* \newline \aths{ |qi| is the total jump rate $\lambda_{m} = \sum_{j=2}^{m}\lambda_{j,m}$ \eqref{eq:ratesarandall}}  */
  assert(qi > 0);
  /* \newline \aths{ sample sojourn time in current state} */
  const double  timi = gsl_ran_exponential( rngtype,  1./qi );
  assert(timi > 0);

  for( const auto &b : tree){
  /* \newline \aths{ each block size |b| must be at least 1} */ 
  assert( b > 0);
      elli[0] +=  timi;
      elli[b] +=  timi; }

  assert( elli[0] > 0);
}


@*1 {\bf gene genealogy and branch lengths}.
\label{sec:genealogy}


generate one gene genealogy or tree and record the branch lengths of
the tree


@<genealogy@>=@#
static void genealogy( const int &leaves,  const double &ia, const double &ikappa, const double &ic, const double &icpr, const double &ig,  const std::vector<double>& lambdab, std::vector<double>& eri )
{
  /* \newline \aths{ |lambdab| is  $\lambda_{n}$ the total jump rate} */

  std::vector<int> tree (leaves, 1) ;
  std::vector<double> tli  (leaves + 1, 0);
  double timi {} ;
  int mergersize {} ;
  int sizenewblock {} ;
  int b {} ;
  assert( tree.size() == leaves); 
  while( tree.size() > 1)
  {
   /* \newline \aths{ update the current branch lengths $\ell_{i}(n)$ \S~\ref{sec:updatecurrentL}}  */
    updatells(lambdab[ tree.size()], tree, tli);

   /* \newline \aths{ sample merger size number of blocks merging \S~\ref{sec:samplemergersize}} */
   mergersize =  samplemerger( static_cast<int>(tree.size()), ia, ikappa, ic, icpr, ig, lambdab);
   assert( mergersize > 1);

   std::shuffle( std::begin(tree), std::end(tree), rng);

   sizenewblock = static_cast<int>( std::accumulate( tree.rbegin(), std::next( tree.rbegin(), mergersize), 0) );
   assert( sizenewblock > 1);
   b = static_cast<int>( tree.size() );
   /* \newline \aths{ remove the merged blocks from the tree} */
   assert( mergersize <= b);
   tree.resize(b - mergersize);
   /* \newline \aths{ add the new block to the tree} */
   tree.push_back( sizenewblock); }

  assert( tree.size() < 2);
  /* \newline \aths{ generated one realisation of branch lengths; update estimate of  $\varrho_{i}(n)$   \S~\ref{sec:updateri}} */
  updateri( leaves, tli, eri );
}



@*1 {\bf estimate $\varrho_{i}(n)$ }.
\label{sec:estimate}


estimate  $\varrho_{i}(n)$

@<estimate@>=@#
static void estimate(  const int& n, const int &nexperiments, const double &ia, const double &ikappa, const double &ic, const double &icpr, const double &ig )
{

@#

  assert( n > 1);

 /* \newline  \aths{ initialise the vector for the total jump rate $\lambda_{i}$} */
  std::vector<double> il (n + 1, 0);
  /* \newline \aths{ record the jump rates |lambdan| \S~\ref{sec:lambdan}} */
   lambdan( n, ia, ikappa, ic, icpr, ig,  il);

  /* \newline \aths{ |veri| records the estimates of $\varrho_{i}(n)$} */
   std::vector< double> veri (n); 
   int r = nexperiments + 1;
   while( --r > 0)
   {
   /* \newline \aths{ \S~\ref{sec:genealogy}} */
      genealogy( n, ia, ikappa, ic, icpr, ig, il, veri);}

  /* \newline \aths{ print the estimate of $\varrho_{i}(n)$} */
  for( const auto &x : veri){
    std::cout << x << '\n' ; }
}


@*1 {\bf the main function}. 
\label{sec:main}

@C

/* \newline \S~\ref{sec:includes} */
@<includes@>@#
/* \newline \S~\ref{sec:gslrng} */
@<rngs@>@#
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
/* \newline \S~\ref{sec:lambdan} */
@<jump rate $\lambda_{n}$@>@#
/* \newline \S~\ref{sec:samplemergersize} */
@<merger size@>@#
/* \newline \S~\ref{sec:updateri} */
@<update  $\varrho_{i}(n)$@>@#
/* \newline \S~\ref{sec:updatecurrentL} */
@<update current lengths@>@#
/* \newline  \S~\ref{sec:genealogy} */
@<genealogy@>@#
/* \newline \S~\ref{sec:estimate} */
@<estimate@>@#


int main( int argc, char *argv[])
{

/* \newline  \aths{ |setup_rng|  \S~\ref{sec:gslrng}} */
   setup_rng( static_cast<unsigned long>(atoi(argv[1])) );

/* \newline \aths{set the parameter values} */
/* \newline \aths{set $\alpha$} */
 const double setalpha  = 0.01;
/* \newline \aths{set $\kappa$} */
 const double setkappa  = 2;
 const double setc      = 1;
 const double setcprime = 1;
 /* \newline \aths{set $\gamma$} */
 const double setgamma  = 1.0 ;

/* \newline \aths{  |estimate| from   \S~\ref{sec:estimate}} */
 estimate( 100, 200000, setalpha, setkappa, setc, setcprime, setgamma );

 gsl_rng_free( rngtype);



 return GSL_SUCCESS ;
}


@* {\bf conclusions and bibliography}. 
\label{sec:conclusions}


We estimate the functionals $\EE{R_{i}(n)}$ when the coalescent
$\set{\xi^{n}; t \ge 0}$ is the
$\delta_{0}$-Beta$(\gamma,2-\alpha,\alpha)$-coalescent.  The
$\delta_{0}$-Beta$(\gamma,2-\alpha,\alpha)$-coalescent can be obtained
from a population model of a haploid panmictic population of constant
size evolving according to sweepstakes reproduction (heavy
right-tailed offspring number distribution) in a random environment
(Definitions~\ref{def:haplrandomalpha} and
\ref{def:alpha-random-one}).  The estimates of $\EE{R_{i}(n)}$ can be
compared to estimates of $\EE{R_{i}^{N}(n)}$, funcionals of an
ancestral process in the domain-of-attraction of $\{\xi^{n}\}$.

Figure~\ref{fig:graph} records an example of an approximation of $\EE{R_i(n)}$


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
       \caption{\it  An example approximation of $\EE{R_{i}(n)}$
graphed as logits against $\log(i/n) - log(1-i/n)$ where $n$ is sample size }
      \label{fig:graph}
       \end{SCfigure}




%%\bibliographystyle{alpha}
%%\bibliography{refs}

\bibliographystyle{plain}
\begin{thebibliography}{Tan11}

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
