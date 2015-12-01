\name{seq_var}
\alias{seq_var}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Intra And Inter Specific Variation For An Alignment
}
\description{
This function calculates and plots inter/intra specific variation
for an alignment.
}
\usage{
seq_var(taxon=NULL, region="ITS", align=NULL, file.align=NULL, 
        file.align.format="fasta", outgroup.name=NULL, 
        taxa.sep=":", col.total=5, col.genus=1, col.genus.syn=3,
        col.species=4, col.strain=5,
        box.cex.axis.text.x=6, box.cex.axis.text.y=6, 
        box.cex.xlab=8, box.cex.ylab=8, 
        intra.fill="purple", inter.fill="orange",
        den.cex.axis.text.x=6, den.cex.axis.text.y=6,
        den.cex.xlab=6, den.cex.ylab=6, frame.col="blue",
        table.cex=8, main=FALSE, main.cex=14,
        file=NULL, width=8, height=8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{taxon}{
the name of the genus to be analyzed. 
}
  \item{region}{
the DNA region of the alignment
}
  \item{align}{
optional, an alignment object in R environment
}
  \item{file.align}{
optional, a file of an alignment
}        
  \item{file.align.format}{ 
alignment file format, default is in FASTA format
}
  \item{outgroup.name}{
the genus/genera names
}
  \item{taxa.sep}{
the separator used to split alignment ID, default is ":". For 
example, as in "genus_name:seq_identifier:genus_synonym:species:strain"
}
  \item{col.total}{
using the \code{taxa.sep}, how many columns the alignment ID can 
be split in to.
}
  \item{col.genus}{
genus name location
}
  \item{col.genus.syn}{
generic synonyms location
}        
  \item{col.species}{
species name location
} 
  \item{col.strain}{
strain information location
}        
  \item{box.cex.axis.text.x}{
size of labels on x-axis of the boxplots
}
  \item{box.cex.axis.text.y}{
size of labels on y-axis of the boxplot
}        
  \item{box.cex.xlab}{
size of the xlab of boxplots
} 
  \item{box.cex.ylab}{
size of ylab of boxplots, 
}        
  \item{intra.fill}{
color to fill the intra specific variation boxplot and density plot.
Default is "purple".
}
  \item{inter.fill}{
color to fill the intra specific variation boxplot and density plot.
Default is "orange",
}        
  \item{den.cex.axis.text.x}{
size of labels on x-axis of the density plots 
}
  \item{den.cex.axis.text.y}{
size of labels on y-axis of the density plots
}        
  \item{den.cex.xlab}{
size of the xlab of density plots
} 
  \item{den.cex.ylab}{
size of the ylab of boxplots
}
  \item{frame.col}{
color of the frames separating inter and intra -specific variation
plots and tables
}        
  \item{table.cex}{
the font size of the summary table
} 
  \item{main}{
whether or not to give a title to the plot, can be "FALSE", "TRUE",
or a string
}
  \item{main.cex}{
font size of the title of the plot
}        
  \item{file}{
whether or not save the plot to a PDF file. 
} 
  \item{width}{
width of the plot
} 
  \item{height}{
height of the plot
}
}
\details{
The sequence ID of the alignment to be analyzed should be in a 
similar format as follows:
"Alternaria:gb_AF229477:Bipolaris:tetramera:strain_BMP_51_31_01",
the separator, ":", and the location of genus / species / strain
names can be changes or rearranged, as long as you provide such 
information to the function properly.  
}
\value{
This function calculate and plot inter/intra specific variation 
for an alignment.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen & C. Andre Levesque
}
% \note{
%% ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
% ~~further links~~
%}
\examples{
data(alignment)
\dontrun{

seq_var(taxon="Claviceps", region="ITS", align=alignment, 
        file.align=NULL, outgroup.name=c("Talaromyces"), 
        box.cex.axis.text.x=6, box.cex.axis.text.y=6, box.cex.xlab=8, 
        box.cex.ylab=8, intra.fill="purple", inter.fill="orange",
        den.cex.axis.text.x=6, den.cex.axis.text.y=6, den.cex.xlab=6, 
        den.cex.ylab=6, table.cex=8, frame.col="blue", main=TRUE, 
        main.cex=14)

seq_var(taxon="Claviceps", region="ITS", align=NULL, 
        file.align="/path/to/alignment_file", outgroup.name=c("Talaromyces"), 
        box.cex.axis.text.x=6, box.cex.axis.text.y=6, box.cex.xlab=8, 
        box.cex.ylab=8, intra.fill="purple", inter.fill="orange",
        den.cex.axis.text.x=6, den.cex.axis.text.y=6, den.cex.xlab=6, 
        den.cex.ylab=6, table.cex=8, frame.col="blue", main=TRUE, 
        main.cex=14)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
