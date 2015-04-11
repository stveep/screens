###
# Annotation and analysis of Plate 11/12 screening data
# Steve Pettitt spettitt@icr.ac.uk
# 23/8/2013
###
# First, need to get the data into R and make sure it is annotated and structured correctly.
# Change the following lines as appropriate for your data:
# Cell line names in the order they appear in the data
cellnames = c("BR13","BR8","JM8Cas9","TdgL7")
# Number of cell lines (set automatically based on the number of names given - therefore this script should work for any number of cell lines)
nlines = length(cellnames)

# This command reads the data from the plate list tab of the Victor excel output (need to save just this sheet as a CSV file) into R. Change file name as needed. ~ represents your home directory.
df <- read.csv("~/work/screens/tdg-parp1/victor.csv",colClasses=c("integer","integer","character","factor","factor","numeric"))

#	Example of the CSV file from the Victor - the first line is automatically read as the column headings:
#	Plate,Repeat,Well,Type,Time,0.1sec (CPS)
#	1,1,A01,M,00:00:06.83,15300
#	1,1,A02,M,00:00:07.09,1967850
#	...

# Now set up vectors that will form new columns in the table, annotating the data
# This gives each replicate a number - rep is an R function to "rep"eat a number a given no. of times, reps is the variable name. If the replicates weren't read together, this needs to be changed:
reps = rep(c(rep(1,384),rep(2,384),rep(3,384)),2*nlines)	# Length = total wells in screen.
# What this looks like:
#> rep(c(rep(1,384),rep(2,384),rep(3,384)),2*nlines)
#   [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#...
# [371] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
#...
# [741] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3
#...
#[1148] 3 3 3 3 3 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#...

# This column shows whether the plate is a plate 11 or a plate 12 (here assumes all 11s are read for one cell line, then all 12s, then similarly for next cell line). If the order of reading was different, this needs changing.
plates = rep(c(11,11,11,12,12,12),nlines)		# Length of this vector should equal total plates in screen.
#> plates
# [1] 11 11 11 12 12 12 11 11 11 12 12 12 11 11 11 12 12 12

# no. of replicates (not readily changeable in this script)
nreps = 3
####
## Now add our new columns to the data frame (df) that was formed by reading the Victor data:
###
# This makes a column with the cell line names directly in the data frame. i gives the position in the vector of cell line names, above. Again, change if the plates were read in a different order.
for (i in 1:length(cellnames)) {
	df$cell[(1+(384*((i-1)*2*nreps))):(384*i*2*nreps)] <- cellnames[i]
}	# This subsetting command takes blocks of the data frame and puts the cell line name in.

# Plate (with a capital P) is the plate number from the Victor data. The "Plate" column heading comes from the csv file that was read in. plate (lower case) is the name of the new column that tells us whether it is a plate 11 or 12 by looking up the Victor plate number in the plates vector. 
df$plate <- plates[df$Plate]
# The replicates column was already formed above and can be added with no further processing
df$rep <- reps

# Now give the columns some new names to make things a bit more R friendly and identify what we are interested in:
names(df) <- c("oldplate","oldrep","well","junk","junk2","ctg","cell","plate","rep")
# This extracts just the columns that we are interested in: i.e. the plate, well no., replicate, cell line name and the raw CTG value, throwing away the old Victor designations etc. The inital comma means it includes all rows.
df <- df[,c("plate","well","rep","cell","ctg")]
# This is what the data frame now looks like:
#> head(df)
#   plate well rep   cell     ctg
#      11  A01   1 Arid5b   15300
#   ...


###
# Reading drug information
###
# These commands read drug names, wells and concentrations from a text file (modified from Rich). In this case it is important to specify the data type in each column - in particular we need to make sure that the concentration is numeric or a lot of things become difficult later. The drug is also read as a character type rather than factors (the default) - this makes things easier when combining data frames. The three factors columns at the end aren't going to be used for anything.
# Start of one of the files to be read (tab-separated):
#  Plate	Well	GeneID	CAS.	Collection number	Concentration [These headings are not meaningful] 
#  1	A01				
#  1	A02				
#  1	A03	5-FU(anti-metabolite)	1000		ND	ND	5mM
#  1	A04	5-FU(anti-metabolite)	100		ND	ND	5mM
drugs11 <- read.table("~/work/hapdip/CompoundLibrary_plate11.txt",head=T,colClasses=c("integer","character","character","numeric","factor","factor","factor"),fill=T,flush=T) # tab separated works by default, in some cases might need to use sep="\t" to specify explicitly.
drugs12 <- read.table("~/work/hapdip/CompoundLibrary_plate12.txt",head=T,colClasses=c("integer","character","character","numeric","factor","factor","factor"),fill=T,flush=T)
# Add a column in these data frames specifying which plate they represent so we don't lose this information when we combine them
drugs11$Plate <- 11
drugs12$Plate <- 12
# Now combine the two data frames into one
drugs <- rbind(drugs11,drugs12)
# Sort out the column names and select just the ones we're interested in, as before
names(drugs) <- c("plate","well","drug","conc","junk","jonk","jank")
drugs <- drugs[,1:4]

# The other information required is which wells are the negative controls (This could probably also be included in the drugs list above but wasn't in the original file at the time so I did it in R). Read in a list of which wells are negative controls:
negwells <- read.csv("~/work/hapdip/negwells.csv",head=F)
# Now set the appropriate wells in the drugs data frame to dmso or puromycin. The "%in%" function is a set intersect operation, used to choose only the wells that are in the list of negative wells.
drugs$drug[drugs$well %in% negwells$V1] = "dmso"
# Just a few positive controls so can specify "on the fly"
poswells <- c("I02","I23","J02","J23")
drugs$drug[drugs$well %in% poswells] = "puromycin"

# > head(drugs)
#   plate well                  drug conc
# 1    11  A01                         NA
# 2    11  A02                  dmso   NA
# 3    11  A03 5-FU(anti-metabolite) 1000
# 4    11  A04 5-FU(anti-metabolite)  100
# 5    11  A05 5-FU(anti-metabolite)  500
# 6    11  A06 5-FU(anti-metabolite)   50

# Now, add the drug information to the main data frame
df <- merge(df,drugs,by=c("plate","well"))  # Both data frames need to have columns called plate and well for this to work.
# Note that merging the data frames can change the order of the rows for dull computer reasons, so lets sort the data:
df <- df[order(df$cell,df$plate,df$rep,df$well),]
# Now we have the following:
#> head(df)
#   plate well rep   cell     ctg                  drug conc
#1     11  A01   1 Arid5b   15300                         NA
#18    11  A02   1 Arid5b 1967850                  dmso   NA
#25    11  A03   1 Arid5b   84810 5-FU(anti-metabolite) 1000
#34    11  A04   1 Arid5b 1365860 5-FU(anti-metabolite)  100
# ...

# This gets a vector containg the names of all the drugs, which is useful later. 
dlist <- levels(as.factor(drugs$drug))  # Could also use unique(drugs$drug).

###
# Now we have all our data in a table where each row is one well, and also has information about the cell line, drug and concentration. This annotation of the data makes it easy (in R terms!) to analyse by breaking it down by any of these. Below, I've used the "by" function to do the normalisation and calculate medians of the replicates.
###

# Make a key to allow platewise operations. This is a unique identifier for each plate in the screen, formed by combining the cell line, plate number (11 or 12) and replicate (1, 2 or 3).
df$key <- paste(as.character(df$cell),as.character(df$plate),as.character(df$rep))
# Log-transform the ctg values
df$log <- log(df$ctg,2)

# Calculate the plate medians. This "by" says: Take the "log" column from df, split it up by the "key" column - i.e. by plate - and apply the "median" function to each of these sets of numbers. The results are stored in a list, which is accessed by the double bracket notation e.g. list[[x]]. Rather than using the key, a more complicated 2nd argument to "by" could also be used (see medians below).
# platemed is the list of plate medians for negative control wells
platemed <- by(df[df$drug=="dmso","log"],		# Values to use (log CTG in this case)
		df[df$drug=="dmso","key"],		# How to split them up (by unique plate ID key)
		median)
			# What function to apply
# What the platemed list looks like:
# > head(platemed)
# df[df$drug == "dmso", "key"]
# Arid5b 11 1 Arid5b 11 2 Arid5b 11 3 Arid5b 12 1 Arid5b 12 2 Arid5b 12 3 
#    21.35274    21.44151    21.35135    21.30562    21.30115    21.18216 
# To get an individual value:
# > platemed[["HAP3 12 1"]]
# [1] 21.36462


###### Simple way of doing it - this loops through the whole data frame and sets values one-by-one, but is very slow for a big data frame; generally loops are not a smart way of doing things in R.
##### for (i in 1:nrow(df)) {
##### df$pcen[i] <- df$log[i] - platemed[[df$key[i]]]
##### }

# A better way, I think, is to make a simple function and use one of the "apply" series of functions to do this for each row of the data frame: 
pcentre <- function(dflog,key,...) {
		dflog - platemed[[key]]
}

df$pcen <- mapply(pcentre,df$log,df$key)

# Calculate medians - in this case the "by" function takes the plate centered values (pcen) and splits them up by cell line, plate AND well (specified in a list) - this means that we have just the three replicates in each subset. Then the median function is applied to these three values: 

reps <- by(df[,"pcen"],					# Take plate centered values
	list(df[,"cell"],df[,"plate"],df[,"well"]),	# Split them by cell line, plate type (11 or 12) and well no.
	median)						# Calculate median of the resulting sets of 3 reps.

# Now put them in the data frame, again using mapply:

dfrep <- function(cell,plate,well,...) {
		reps[[cell,as.character(plate),well]]
}

df$median <- mapply(dfrep,df$cell,df$plate,df$well)

# Now we have the data ready to analyse. Z scores etc., can be calculated by similar methods.
# NB at this point, each replicate has its own line in the table still, but all 3 lines have the same value in the median column.

## A few useful analyses...

# We can use the "reps" "by" function output above to calculate the difference in survival between two given cell lines.
# e.g. Arid5b and Chd6:

aridch <- function(plate,well) {
	reps[["Arid5b",as.character(plate),well]] - reps[["Chd6",as.character(plate),well]]
}

df$aridch <- mapply(aridch,df$plate,df$well)


## A function for a simple kill curve, using the lattice package. To specify which drug to use, set d to the number of the drug in the dlist vector.
library(lattice)
# Default colours are ugly, make some nicer ones:
library(RColorBrewer)
kcol <- brewer.pal(8,"Set1")
# Start our function definition:
killc <- function(d) {
xyplot(								#Make a scatter plot
	df[df$drug==dlist[d],]$pcen ~ log(df[df$drug==dlist[d],]$conc,10), 
	# [Above] plots the plate centered log CTG for the given drug, against log concentration.
	group = df[df$drug==dlist[d],]$cell,			# Grouped by cell line
	xlab=dlist[d],ylab="Plate-centered CTG value (log2)",		# Label axes
	ylim=c(-6,1),
	auto.key=T,						# Add legend
	type=c("p","smooth"),					# Plot "p"oints and a smoothed curve
	# A potential improvement is to add a proper sigmod fit...
	par.settings=simpleTheme(pch=19,							# Plot points as small filled circles	
	col=kcol))						# Colour points and lines with our nice colours
}	# End of function definition.				


# Make a PDF of all kill curves in the data.
kcpdf <- function(f="all.pdf") {
	pdf(f)
	for (i in 1:length(dlist)) {
		print(killc(i))
	}
	dev.off()
}

# What about just getting kill curves for drugs that show a difference between the two cell lines?
intpdf <- function(f="interest.pdf",thresh=1.5) {				# Default minimum diff. is +/- 1.5
	interest <- unique(df[abs(df$aridch) > thresh,]$drug)	# Get the drugs of interest
	pdf(f)							# Open a new PDF
	for (i in interest) {					# Go through list of drugs of interest
		d = which(dlist == i)				# Which number in the list is this drug?
		print(killc(d))					# Print appropriate kill curve
	}							# End for loop
	dev.off()						# Finish writing the PDF
}								# End of function	



#### Below is my attempt at sigmoid fit - needs work! There is probably a package at CRAN that does this much more easily...
#### Argument d is the number of the drug in dlist vector.
###kc <- function(d,prms=list(a=-4,b=-5,c=1)) {
###	require(lattice)
###	sigmoid <- function(params,x) {
###		params[1]/(1+exp(params[2]*(x-params[3]))) 
###	}
###	pdf(d)
###	lns <- list()
###	for (i in levels(as.factor(df$cell))) {
###		x=log(df[df$drug==dlist[d],]$conc,10)
###		y=df[df$drug==dlist[d],]$pcen
###		curve <- nls(y~a/(1+exp(b*(x-c))),start=prms)
###		params <- coef(curve)
###		lns[[i]] <- sigmoid(params,x)
###	}
###	print(xyplot(df[df$drug==dlist[d],]$pcen ~ log(df[df$drug==dlist[d],]$conc,10), group = df[df$drug==dlist[d],]$cell,pch=22))	
###	for (j in lns) {
###		points(j[order(x)],x[order(x)],type="l")
###	}
###	dev.off()
###} 
