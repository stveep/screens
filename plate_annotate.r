wellrange <- function(range,first="row"){
	#Check format
	if (!grepl("[A-Za-z]\\d+-[A-Za-z]\\d+",range)) {
		warning("Range must be of format A2-C4\n")
	}
	#Split out the start and end of the range
	range <- unlist(strsplit(range,"-"))
	if (length(range) > 2) { warning("More than one dash in range given") }
	range_start = range[1]
	range_end = range[2]
	range_start_row = toupper(sub("\\d+","",range_start))
	range_start_col = as.integer(sub("[A-Za-z]+","",range_start))
	range_end_row = toupper(sub("\\d+","",range_end))
	range_end_col = as.integer(sub("[A-Za-z]+","",range_end))
	row_range = LETTERS[which(LETTERS==range_start_row):which(LETTERS==range_end_row)]
	col_range = seq(from=range_start_col, to=range_end_col)

	if (first=="row") {
		first = row_range
		second = col_range
		label1="Row"
		label2="Column"
	} else {
		first = col_range
		second = row_range
		label2="Row"
		label1="Column"
	}
	output = data.frame(matrix(nrow=1,ncol=2))
	for (i in first) {
		for (j in second) {
			output <- rbind(output,c(i,j))
		}
	}
	names(output) <- c(label1,label2)
	#First row has NA values
	output <- output[-1,]
	#Sort by first then second (i.e. row then column default)
	output[order(output[,1],output[,2]),]
}

condition_wells <- function(wells,cname,value) {
	df <- wellrange(as.character(wells))
	df[,as.character(cname)]=value
	df
	
}

annotate_plate <- function(df,type=96,plate=1,area=2,condition=3,value=4) {
# For stuff replicated across plates, can use this with by functions	
# DF provided needs Plate, Area, Condition, Value columns with indexes as in default args

	# set up the plate types 
	nplates = length(unique(df[,plate]))
	if (type == 96) {
		wells = wellrange("A1-H12")
	} else {
		if (type == 384) {
			wells = wellrange("A1-P24")
		}
	}	
	
	# Apply the annotation
	# If Plates field is empty, apply to all plates
	# To apply to all wells, need to enter the full range e.g. A1-H12) - could add something to do this automatically
	plate_empty <- df[is.na(df[,plate]),]

	if (nrow(plate_empty) != 0) {
	# Remove old empty plate rows
		df <- df[!is.na(df[,plate]),]

	# Makes a list of duplicate data frames, one for each plate
		new_rows <- lapply(unique(df[,plate]), function(x,pe) { temp = pe; temp$Plate = x; temp },pe=plate_empty)

	# Add the new ones
		new_rows$new = df
		df <- do.call(rbind,new_rows)

	}
	# Remove any lines where the condition or value is NA	

	# Set up the plates data frame	
	lst = list()
	for (i in unique(df[,plate])) {
		lst[[i]] = wells
		lst[[i]][,"Plate"] <- i	
	}
	plates <- do.call(rbind, lst)

	#anno is a list of data frames	
	#Now merge all of these with the plates data frame
	annolist = list()
	for (i in unique(df[,condition])) {
		anno <- df[df[,condition] == i,]
		# Expand well ranges
		anno <- apply(anno,1,function(x) {dft <- condition_wells(x[area],x[condition],x[value]); dft$Plate = x[plate]; dft})
	#	anno$plates <- plates
	#	annolist[[i]] <- Reduce(function(x,y) merge(x,y,by=c("Plate","Row","Column"),all=T), anno)
		tempplates <- plates
		for (j in anno) {
			tempplates <- merge(tempplates,j,all=T)
		}
		# Generates NAs for some reason, remove
		tempplates <- tempplates[!is.na(tempplates[,4]),]
		annolist[[i]] <- tempplates
	}
	red <- Reduce(function(x,y) merge(x,y,all=T),annolist)	
	red[!apply(red,1,function(x) any(is.na(x))),]
	
}
library(lattice)
library(RColorBrewer)
kcol <- c(brewer.pal(9,"Set1"),brewer.pal(9,"Set2"))
killc <- function(df,...) {
xyplot(                                                         #Make a scatter plot
        df$Survival ~ log(df$Concentration,10),
        # [Above] plots the plate centered log CTG for the given drug, against log concentration.
        group = df$siRNA,                    # Grouped by cell line
 #       ylim=c(-6,1),
        auto.key=T,                                             # Add legend
        type=c("p","smooth"),                                   # Plot "p"oints and a smoothed curve
        # A potential improvement is to add a proper sigmod fit...
        par.settings=simpleTheme(                                                        # Plot points as small filled circles   
        col=kcol,
	pch=seq(from=15,to=25),
	),...)                                              # Colour points and lines with our nice colours
}  
