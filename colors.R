#colors
rainbow(n, s = 1, v = 1, start = 0, end = max(1, n - 1)/n, alpha = 1)
heat.colors(n, alpha = 1)
terrain.colors(n, alpha = 1)
topo.colors(n, alpha = 1)
cm.colors(n, alpha = 1)
###vidrisi
img <- function(raster, col) {
  image(raster, col=col, asp=1, axes=FALSE, xaxs="i", xaxt='n', yaxt='n', ann=FALSE)
}
img(raster, rev(viridis(64)))

library(dichromat)
img(AS, dichromat(rev(viridis(64)), "tritan"))

# replacing NA's by zero
r[is.na(r[])] <- 0 
####
# training and test data sets
set.seed(22)
sampled_rows <- sample(seq_len(nrow(Gary)), 1529, replace=FALSE)
leftover_rows <- setdiff(seq_len(nrow(Gary)), selected_rows)
train <- Gary[sampled_rows, ]
leftover <- Gary[leftover_rows, ]
str(train)
str(leftover)
#Bootstrapping is often used for: Calculation of confidence intervals (and estimation of the standard errors)
#Estimation of the bias of the point estimates
#https://ropensci.org/tutorials/getlandsat_tutorial.html
#https://aws.amazon.com/public-datasets/landsat/
#landsat data frm ra package
library("getlandsat")
#list scences
(res <- lsat_scenes(n_max = 50))

#list scence files 
lsat_scene_files(x = res$download_url[1])
#Get an image
#Returns path to the image
lsat_image(x = "LC80101172015002LGN00_B4.TIF")
#Visualize
x <- lsat_cache_details()
img <- raster(x[[1]]$file)
plot(img)

setwd("C:/Users/uqrdissa/Downloads")
wnk <- read.csv("wildnetkoalalocations-2017-06-02.csv", header = TRUE) # centroids for full study area as some sros are required.

#When I get that error I issue plot erros  # r remove add= TRUE.
plot.new() 
#and if I don't see a plot window, I issue  quartz() as well.
# write a raster stack in layer names.
 library(sos)
de <- findFn('ROC random forest')

### how to address NA in linear regression
#This is how you can get the raster data

library(raster)
names = c('...','...','...','...','...')
s <- stack(names)
y <- values(s)
#You could now do something like this.

x <- log(c(10,20,30,40,50))
# need to exclude the rows that are all NA
i <- rowSums(is.na(y)) < ncol(y)
coef <- apply(y[i, ], 1, function(y) lm(y ~ x)$coefficients)
aa <- matrix(NA, ncol=2, nrow=length(i))
aa[i, ] <- coef
b <- brick(s, nl=2)
values(b) <- aa
#But you do not need to do that. To do regression like this, I would do

fun <- function(y) { lm(y ~ x)$coefficients }
r <- calc(s, fun)
#But because you have cells with only NA values (across the layers) this will fail (like in the apply above). You need to write a function to catch these cases:
  
  funa <- function(y) { 
    if(all(is.na(y))) {
      c(NA, NA)
    } else {
      lm(y ~ x)$coefficients 
    }
  }
r <- calc(s, funa)
#Or for a much faster approach

X <- cbind(1, y)
invXtX <- solve(t(X) %*% X) %*% t(X)
quickfun <- function(i) (invXtX %*% i)
m <- calc(s, quickfun) 
names(m) <- c('intercept', 'slope')

# set margins labels 
par(mgp=c(3,1,0),mar(5,4,4,2)+0.1)


##aggregate from 40x40 resolution to 120x120 (factor = 3)
meuse.raster.aggregate <- aggregate(meuse.raster, fact=3)
res(meuse.raster.aggregate)
#[1] 120 120

#disaggregate from 40x40 resolution to 10x10 (factor = 4)
meuse.raster.disaggregate <- disaggregate(meuse.raster, fact=4)
res(meuse.raster.disaggregate)


#save plots as objects Use recordPlot and replayPlot:
  
  plot(BOD)
plt <- recordPlot()
plot(0)
replayPlot(plt)


#load all packages and download inf any needed.
lapply(.packages(all.available = TRUE), function(xx) library(xx,     character.only = TRUE))

####R help
R Help: help() and ?
The help() function and ? operator are useful only if you already know the name of the function that you wish to use.
The help() function and ? help 
To access documentation for the standard lm (linear model) function,
help(lm) 	help("lm")	 ?lm  	?"lm" 
To access help for a function in a package that's not currently loaded, specify in addition the name of the package:
  help(rlm, package="MASS").
Use the help() function to access information about a package in your library 
help(package="MASS")  displays an index of available help pages for the package 
Execute examples illustrating how the functions work. 
example() command: e.g., example(lm).
Vignettes and Code demons:
  browseVignettes(), vignette(), demo()
Vignettes illustrate and explain facilities in the package. 
browseVignettes() function
browseVignettes() opens a list of vignettes from all of your installed packages in your browser
browseVignettes(package=package-name)  browseVignettes(package="survival")) 
To shows the vignettes, if any, for a particular package. 
vignette() displays a list of vignettes in text form.
vignette("vignette-name") 
If the vignette name is not unique): vignette("timedep") or vignette("timedep", package="survival") 
Packages code demonstrations ("demos").
demo()
demo() lists all demos for all packages in your library
demo(package="package-name")  demo(package="stats")) lists demos in a particular package. 
demo("nlm")
demo("nlm", package="stats"), 
Searching for Help Within R
Discovering functions and other objects. Use the help system to obtain complete documentation for these functions: for example, ?apropos.
apropos()
The apropos() function searches for objects, including functions, directly accessible in the current R session This may be a literal string or a regular expression to be used for pattern-matching (see ?"regular expression"). By default, string matching by apropos() is case-insensitive. 
apropos("^glm") returns the names of all accessible objects that start with the (case-insensitive) characters "glm".
help.search() and ??
The help.search() function scans the documentation for packages installed in your library. For help.search("^glm") searches for help pages, vignettes, and code demos that have help "aliases," "concepts," or titles that begin (case-insensitively) with the characters "glm". The ?? operator is a synonym for help.search(): for example, ??"^glm".
RSiteSearch()
RSiteSearch() uses an internet search engine to search for information in function help pages and vignettes for all CRAN packages
RSiteSearch() requires an active internet connection and doesn't employ regular expressions. Braces may be used to specify multi-word terms; otherwise matches for individual words are included. For example, RSiteSearch("{generalized linear model}") 
findfn() and ??? in the sos package, which is not part of the standard R distribution but is available on CRAN, provide an alternative interface to RSiteSearch().
help.start()
help.start() starts and displays a hypertext based version of R's online documentation in your default browser that provides links to locally installed versions of the R manuals, a listing of your currently installed packages and other documentation resources.
R Help on the Internet
search.r-project.org (which is the site used by RSiteSearch) and Rseek.org.
For recent packages and information go to 
rviews.rstudio.com
#######
This answers the question in the title "how to count the number of points from a shapefile within each cell"

library(rgdal)           # this package to read and manipulate shapefiles
library(raster)          # this package for rasters
shp <- readOGR("my_shape_file.shp")        # read in your points from the shape file
ras <- raster("my_raster_file")            # read in your raster
shp <- spTransform(shp, projection(ras))   # make sure that the two spatial objects share the same coordinate system
cells <- cellFromXY(ras,shp)     # find which cells the points are in
table(cells)                     # and make a table of number of occurences in each cell

#I developed an R package to clear console that will do this, borrowing from the suggestions above. The package is called called mise, as in "mise en place." You can install and run it using
install.packages("mise")
library(mise)
mise()
Note that mise() also deletes all variables and functions and closes all figures by default. 
#To just clear the console, use mise(vars = FALSE, figs = FALSE)

# looping
for (i in 10:0) {
   cat(i, " \r")
   flush.console()
   Sys.sleep(1)
}
For vectors there is seq_along, for DataFrames you may use seq_len

for(i in seq_len(nrow(the.table)){
do.stuff()
}
#
#another way to do it.
rasters <- stack(myfullstack)
values <- getValues(rasters) # all raster valuves
is.na(r.val) <- sapply(values, is.infinite) # convert infinite to na.
r.val = na.omit(r.val) # remvoe na
pca <- prcomp(r.val, scale = TRUE)
summary(myPCA)