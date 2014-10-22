\name{RAM.dates}
\alias{RAM.dates}

\title{Date Formatting for RAM}

\description{
This help page details the expected format for dates in RAM.
}

\details{For all functions expecting some type of date data, you 
will need to specify which column of the metadata table contains 
that information. The date information will likely be encoded as 
a character vector from read.meta, so these functions will try to 
coerce it to a Date object (see \code{\link{Date}} and 
\code{\link{as.Date}}), with a warning. These functions are 
expecting the date information to be in YYYY-MM-DD format.}
