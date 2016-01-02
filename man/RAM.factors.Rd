\name{RAM.factors}
\alias{RAM.factors}
\title{Factor Formatting for RAM}
\description{This help page details how to pass metadata 
arguments in RAM.}
\details{Many functions will expect arguments such as 
\code{meta} and \code{factors} (possibly \code{row.factor} or 
\code{col.factor}). These functions are expecting the full 
metadata table for \code{meta} (which you probably read into R 
using \code{\link{read.meta}}).
The other argument, \code{factor}, should be a \emph{named} 
character vector. The values of this vector should be the 
columns to be used from \code{meta}, while the names of the 
vector should be the labels you wish to have displayed in the 
plots. There are several ways to name a character vector:
% NOTE: in an extremely frustrating bug (which was supposedly
% fixed in R 2.13) Sexpr will crash on multi-line input, so we
% break it in to many, many lines.
\Sexpr[echo=TRUE, results=verbatim]{
my.vec <- c(This="is", a="named", character="vector")}
\Sexpr[echo=TRUE, results=verbatim]{names(my.vec)}
\Sexpr[echo=TRUE, results=verbatim]{cat(my.vec)}
Notice that \code{myvec} has \emph{names} "\code{This}", 
"\code{a}", "\code{character}", but has \emph{values} 
"\code{is}", "\code{named}", "\code{vector}". It is the names 
that will be used to label graphs in \code{RAM}, but the values 
that will be used to extract the actual data. This is useful if 
you have more complicated names; say we want the data from the 
column named "Precip_14d_before_harvest", but we want a nicer 
label for the plot. We can do the following:
\Sexpr[echo=TRUE, results=verbatim]{
my.vec <- "Precip_14d_before_harvest"}
\Sexpr[echo=TRUE, results=verbatim]{
names(my.vec) <- "Precipitation (14 d. prior to Harvest, mm)"}
Now we will be able to plot the value of the 
"\code{Precip_14d_before_harvest}" column, but we will have the 
(much nicer!) label "Precipitation (14 d. prior to
Harvest, mm)" appear in our plots.
}
