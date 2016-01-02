\name{RAM.plotting}
\alias{RAM.plotting}
\title{
Creating Plots with RAM
}
\description{
This help page details the standards for RAM plotting functions.
}
\section{Overview}{
The RAM package contains many functions to produce plots and
visualizations for metagenomic data. Currently, the plotting
functions are grouped into 'casual' and 'publication' categories.
The 'casual' plotting functions only accept a \code{file} 
argument and produce a \code{.tiff} file automatically. They 
are meant to quickly highlight certain aspects of the data, but 
are not meant to be published. The 'publication' quality plots 
accept many more graphing parameters, and should be of suitable 
quality for future publication. All 'publication' plots should 
accept the following parameters, in addition to those required 
to produce the plot:
\itemize{
\item \code{"file"} the file name for the plot.
\item \code{"ext"} the file type for the plot (see below).
\item \code{"height"} the height of the plot, in inches.
\item \code{"width"} the width of the plot, in inches.
}
Additionally, the following parameters are accepted by some 
functions:
\itemize{
\item \code{"bw"} should the plot be in black and white?
}
For 'casual' plots, if \code{file} is not provided, the plot is
displayed to the default graphics device (usually a new window),
otherwise a \code{.tiff} file is created.
For 'publication' plots, if neither \code{file} nor \code{ext} 
are provided, the plot is displayed to the default graphics 
device (usually a new window). If both \code{file} and 
\code{ext} are provided, a file of type \code{ext} is created 
at \code{file}.
If only one of \code{file} or \code{ext} is given, an error is
raised.
In either case, the file extension will automatically be 
appended to the end of \code{file}, if \code{file} does not 
already end in the appropriate extension.
For example, if \code{file = ~/my/path.tiff} and 
\code{ext = png}, the file will be called 
\code{~/my/path.tiff.png}; but if \code{file = ~/my/path.png}, 
the file will just be called \code{~/my/path.png}. 
Finally, if \code{file = ~/my/path}, the file will be called 
\code{~/my/path.png}.
}
\section{ggplot2}{
Furthermore, some of the 'publication' quality plots allow for 
a \code{ggplot2} argument. If \code{ggplot2} is \code{TRUE}, 
then the plot will be produced using the \code{ggplot2} package,
and a \code{ggplot} object will be returned. This allows for 
additional, personal customization of the plot for those who 
have used the package before. If \code{ggplot2} is \code{FALSE},
then the plot will be created using the base plotting functions.
}
\section{File Types}{
For 'publication' quality plots, the following file types are
supported (use any of the following values for the \code{ext}
argument): \code{"pdf"}, \code{"png"}, \code{"tiff"}, 
\code{"svg"}, \code{"bmp"}, \code{"jpeg"}.
}
\section{Note}{
If \code{file} is given without a directrory
(e.g. \code{file} = my_fancy_file), then that file will be 
created in the current working directory (see \code{?getwd} and 
\code{?setwd} for more information).
The current default resolution is 1000 dpi for all plots.
}
\section{See Also}{
\code{\link{ggplot}}
}
\section{Author(s)}{
Wen Chen and Joshua Simpson.
}
