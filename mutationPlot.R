
library(biomaRt)
library(ggplot2)
library(gtable)
library(grid) 
require(ggthemes)
library(FRACTION)
library(IRanges)
library(plyr)
library(plotly)
require(stats)
options(warn = -1)



mutationPlot <- function(typeOfGeneName = NULL, 
						NCBI_Build = "37", varForAxisYup = NULL,
						varForAxisYdown = NULL, varForSize = NULL, 
						varForShape = NULL, drawOnlyTranscript = FALSE,
						outputFormat = NULL, onlyExons = FALSE, 
						databaseUp = NULL, databaseDown = NULL, 
						useDavoliData = FALSE, drawInteractive=FALSE,
						mutationTypeUp = NULL, mutationTypeDown = NULL,
						windowStart = NULL, DataFile = NULL,
						windowEnd = NULL, gene = NULL){
						
	
	if(!is.null(DataFile)){
		tsvData <- read.csv(DataFile, header=T, sep="\t", stringsAsFactors=F)
	} else {
		useDavoliData <- TRUE
		tsvData <- NULL
	}
	if(useDavoliData){
		tumorDaten <- read.csv("Mutation_Dataset.txt", header=T, sep="\t", stringsAsFactors=F)
	}
	
	if(!is.null(tsvData)){
		genesInFile <- unique(tsvData$gene)
		transcriptsInFile <- unique(tsvData$transcript)
	} else {
		genesInFile <- NULL
		transcriptsInFile <- NULL
	}
	
	if (drawOnlyTranscript & !is.null(transcriptsInFile)){
		data <- transcriptsInFile
		allElementData <- tsvData$transcript
		xAxisAnno <- 'Transcript: '
	} else if(!is.null(gene)) {
		if(typeOfGeneName=="hgnc"){
			data <- c()
			for(j in 1:length(gene)){
				name <- getEnsemblGeneName(gene[j], NCBI_Build)
				if(nrow(name)==0){
					message(paste('Gene with name', gene,'not found in biomaRt.'))
				} else {
					data <- c(data, name$ensembl_gene_id[1])
				}
			}
		} else if(typeOfGeneName=="ensembl"){
			data <- gene
		} else {
			stop('Please specify the typeOfGeneName (choose "hgnc" or "ensembl")')
		}
		if(!is.null(tsvData)){
			allElementData <- tsvData$gene
		}
		xAxisAnno <- 'Gene: '
	} else {
		if(!is.null(tsvData)){
			data <- genesInFile
			allElementData <- tsvData$gene
			xAxisAnno <- 'Position in genome, Gene: '
		}
	}

	for (i in 1:length(data)){
		transcript <- drawOnlyTranscript
		name <- getHGNC_symbol(data[i], transcript, NCBI_Build)
		if(!is.null(tsvData)){
			elData <-tsvData[allElementData == data[i],]
			chromosome <- unique(elData$chromosome)
		}
		if(useDavoliData){
			if(!is.null(tsvData)){
				elementDataUp <- prepareMutations(elData, up=T, mutationType=mutationTypeUp,
					varForAxisYup, varForAxisYdown, varForShape, varForSize)
				elementDataDown <- processTumorFileDavoli(name$hgnc_symbol,tumorDaten)
				elementDataDown <- prepareMutations(elementDataDown, up=F, mutationType=mutationTypeDown,
					varForAxisYup, varForAxisYdown, varForShape, varForSize)
			} else {
				elementData <- processTumorFileDavoli(name$hgnc_symbol,tumorDaten)
				chromosome <- unique(elementData$chromosome)
				elementDataUp <- prepareMutations(elementData, up=T, mutationType=mutationTypeUp,
					varForAxisYup, varForAxisYdown, varForShape, varForSize)
				elementDataDown <- prepareMutations(elementData, up=F, mutationType=mutationTypeDown,
					varForAxisYup, varForAxisYdown, varForShape, varForSize)
				
			}
		} else {
			elementDataUp <- prepareMutations(elData, up=T, mutationType=mutationTypeUp,
					varForAxisYup, varForAxisYdown, varForShape, varForSize)
			elementDataDown <- prepareMutations(elData, up=F, mutationType=mutationTypeDown,
					varForAxisYup, varForAxisYdown, varForShape, varForSize)
		}
		position <- getPosition(data[i], transcript, chrNr=chromosome,
								NCBI_Build)
		exonPositions <- getPosition(data[i], transcript, chrNr=chromosome, 
									 NCBI_Build, drawAllExons=T)
		userCoord <- NULL
		if(!is.null(windowStart) & !is.null(windowEnd)){
			if(findInterval(windowStart[i], c(position$start_position,position$end_position) &
				windowEnd[i], c(position$start_position,position$end_position))){
					userCoord <- data.frame(position=seq(windowStart[i],windowEnd[i],by=1))
			} else {
				message(paste("Interval for gene", data[i],"must be within interval [",
					position$start_position,"-",position$end_position,"]"))
			}
		}
		leerUp <- isEmpty(elementDataUp)
		leerDown <- isEmpty(elementDataDown)
		
		if(is.null(varForShape)){
			varForShape <- "Not selected"
		}
		
		if(is.null(varForSize)){
			varForSize <- "Not selected"
		}
		
		if(drawInteractive){
			void <- drawPlot(data[i], elementDataUp, elementDataDown, position, 
							chrNr=chromosome, NCBI_Build, xAxisAnno=paste(
							xAxisAnno, data[i],"\nShape: ",varForShape,
							" Size: ", varForSize),transcript, exonPositions, leerUp,
							leerDown, onlyExons, databaseUp, databaseDown,
							drawPlotly=T, userCoord=NULL)
		}

		if(onlyExons){
			exonRanges <- getExonRanges(exonPositions)
			if(!leerUp){
				elementDataUp <- removeIntronMut(exonRanges, elementDataUp)
			}
			if(!leerDown){
				elementDataDown <- removeIntronMut(exonRanges, elementDataDown)
			}
			if(!is.null(userCoord)){
				userCoord <- removeIntronMut(exonRanges, userCoord)
			}
			relExonRanges <- getRelativeExons(exonRanges)
			if(nrow(elementDataUp)!=0){
				elementDataUp <- getRelativeDataPos(elementDataUp, exonRanges, relExonRanges)
			}
			if(nrow(elementDataDown)!=0){
				elementDataDown <- getRelativeDataPos(elementDataDown, exonRanges, relExonRanges)
			}
			if(!is.null(userCoord)){
			userCoord <- getRelativeDataPos(userCoord, exonRanges, relExonRanges)
			}
			exonPositions <- data.frame(exon_chrom_start=start(relExonRanges),
				exon_chrom_end=end(relExonRanges))
			position <- range(relExonRanges)
			name <-getHGNC_symbol(data[i], transcript, NCBI_Build)
			position <- data.frame(hgnc_symbol=name$hgnc_symbol, start_position=start(position),
				end_position=end(position))
		} else {
			exonRanges <- NULL
			relExonRanges <-NULL
		}
		leerUp <- isEmpty(elementDataUp)
		leerDown <- isEmpty(elementDataDown)

		drawing <- drawPlot(data[i], elementDataUp, elementDataDown, position, 
							chrNr=chromosome, NCBI_Build, xAxisAnno=paste(
							xAxisAnno, data[i],"\nShape: ",varForShape,
							" Size: ", varForSize),transcript, exonPositions, leerUp,
							leerDown, onlyExons, databaseUp, databaseDown,
							drawPlotly=F, userCoord=userCoord, oldExons=exonRanges, 
							exRanges=relExonRanges)
		postscript(paste("chr",chromosome,"plot",i,position$hgnc_symbol,data[i],"up",
			varForAxisYup,"down",varForAxisYdown,outputFormat,sep="."),
			width = 14, height = 8)
		grid.newpage()
		grid.draw(drawing)
		dev.off()
	}
}

prepareMutations <- function(elementData, up, mutationType, varForAxisYup,
							varForAxisYdown, varForShape, varForSize){
	if(!isEmpty(elementData)){
		elementDataNamed <- renameColumns(elementData, up, varForAxisYup,
									varForAxisYdown, varForShape, varForSize)
	}
	if(!is.null(mutationType)){
			elementData <- elementDataNamed[grep(paste(mutationType,collapse="|"),
											elementData$annotation),]
	} else {
		elementData <- elementDataNamed
	}
	return(elementData)
}

isEmpty <- function(dataset){
	if(nrow(dataset)==0){
		leer <- TRUE
	} else {
		leer <- FALSE
	}
	return(leer)
}

renameColumns <-function(dataset, up, varForAxisYup = NULL, varForAxisYdown = NULL,
						varForShape = NULL, varForSize = NULL){
	if(up) {
		colnames(dataset)[which(colnames(dataset) == varForAxisYup)] = "varForAxisYup"
	} else {
		dataset["varForAxisYdown"] <- dataset[varForAxisYdown]*-1
	}
	if(!is.null(varForShape)){
		colnames(dataset)[which(colnames(dataset) == varForShape)] = "varForShape"
	} else {
		dataset["varForShape"] <- sample(1, nrow(dataset), replace=T)
	}
	if(!is.null(varForShape)){
		colnames(dataset)[which(colnames(dataset) == varForSize)] = "varForSize"
	} else {
		dataset["varForSize"] <- sample(1, nrow(dataset), replace=T)
	}
	return(dataset)
}

processTumorFileDavoli <-function(gene, tumorDaten){
	testset <- tumorDaten[tumorDaten$Gene == gene,]
	if(isEmpty(testset)){
		message(paste('Gene with name', gene,'not found in tumor data.'))
	}
	colnames(testset)[1] <- "gene" 
	colnames(testset)[6] <- "annotation"

	m <- gregexpr("^[0-9]{1,2}", testset$Genome.position.hg19, perl=TRUE)
	match <-regmatches(testset$Genome.position.hg19, m)
	testset["chromosome"] <-unlist(match)

	m <- gregexpr(":([0-9]+)-", testset$Genome.position.hg19, perl=TRUE)
	match <-regmatches(testset$Genome.position.hg19, m)
	match <- substring(match, 2)
	match <- substring(match,1,nchar(match)-1)
	testset["position"] <- as.integer(match)

	withCounts <-ddply(testset,.(gene,Reference, Mutation, Protein_Change, annotation, chromosome, position),nrow)
	colnames(withCounts)[ncol(withCounts)] <- "counts"

	return(withCounts)
}

getPosition <- function(elementName, transcript = FALSE, chrNr, NCBI_Build, drawAllExons = FALSE){
	ensembl <- getEnsemblObject(NCBI_Build)
	if (transcript & drawAllExons==F){
		position <- getBM(attributes=c("hgnc_symbol","ensembl_transcript_id", "transcript_start", "transcript_end"), 
				filters=c("chromosome_name","ensembl_transcript_id"), values=list(chrNr, elementName), 
				mart=ensembl)
	} else if (transcript & drawAllExons){
		position <- getBM(attributes=c("exon_chrom_start", "exon_chrom_end"),
				filters=c("chromosome_name","ensembl_transcript_id"), values=list(chrNr, elementName),
				mart=ensembl)
	} else if (transcript==F & drawAllExons){
		position <- getBM(attributes=c("exon_chrom_start", "exon_chrom_end"),
				filters=c("chromosome_name","ensembl_gene_id"), values=list(chrNr, elementName),
				mart=ensembl)
	} else {
		position <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id", "start_position", "end_position"), 
				filters=c("chromosome_name","ensembl_gene_id"), values=list(chrNr, elementName), 
				mart=ensembl)
	}
	return(position)
}

getEnsemblObject <- function(NCBI_Build){
	if(length(grep("7|19", NCBI_Build))){
		ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "grch37.ensembl.org",
			path = "/biomart/martservice", dataset = "hsapiens_gene_ensembl")
	} else if( length(grep("8|20", NCBI_Build)) ){
		ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "www.ensembl.org")
	}  else {
		stop('NCBI reference assembly number is not recognized')
	}
	return(ensembl) 
}

getHGNC_symbol <- function(ensembl_gene_id, transcript, NCBI_Build){
	ensembl <- getEnsemblObject(NCBI_Build)
	if(transcript==TRUE){
	name <- getBM(attributes=c("hgnc_symbol"), filters=c("ensembl_transcript_id"),
				values=list(ensembl_gene_id), mart=ensembl)
	} else {
	name <- getBM(attributes=c("hgnc_symbol"), 
				filters=c("ensembl_gene_id"), values=list(ensembl_gene_id), 
				mart=ensembl)
	}
	return(name)
}

getEnsemblGeneName <- function(hgnc_symbol, NCBI_Build){
	ensembl <- getEnsemblObject(NCBI_Build)
	name <- getBM(attributes=c("ensembl_gene_id"), 
				filters=c("hgnc_symbol"), values=list(hgnc_symbol), 
				mart=ensembl)
	return(name)
}

drawPlot <- function(element, datasetUp, datasetDown, position, chrNr, NCBI_Build,
					xAxisAnno = NULL, transcript = FALSE, exonPositions, leerUp, leerDown,
					onlyExons = FALSE, databaseUp = NULL, databaseDown = NULL,
					drawPlotly = FALSE, userCoord = NULL, oldExons=NULL, 
					exRanges=NULL){
					
	prot_length <-(as.numeric(position$end_position)-as.numeric(position$start_position))
	if(!is.null(userCoord)){
		prot_length <- userCoord$position[length(userCoord$position)] - userCoord$position[1]
	}
	ratio_prot <-prot_length*0.1
	df <- data.frame(x=seq(from=position$start_position,to=position$end_position, length.out=6), y=c(0:5),v1=c(0:5))

	if(varForAxisYup=="MAF"){
		y_up <- log(datasetUp$varForAxisYup)
	} else {
		y_up <- datasetUp$varForAxisYup
	}
	
	if(leerUp){
		p1 <-ggplot() + 
		geom_point(aes(x = df$x, y = df$y, size = df$v1), fill=NA, color=NA, shape=21) + 
		theme(legend.position="none", legend.text =  element_blank()) 
	} else {
		p1 <-ggplot() + 
		geom_point(data = datasetUp, mapping = aes(x = datasetUp$position,
			y = y_up, color=datasetUp$annotation , 
			shape=factor(datasetUp$varForShape)
			)) +
		geom_segment(aes(x = datasetUp$position, xend = datasetUp$position, 
			y = min(y_up)-1,
			yend = y_up), linetype = "dotted", colour="gray") + 
		theme(legend.text = element_text(size=8, face="bold"),
			legend.position = 'right') +
		scale_x_continuous(limits=c(position$start_position,position$end_position))
	}
	
	if(leerDown){
		p3 <- ggplot() + 
		geom_point(aes(x = df$x, y = df$y, size = df$v1), fill=NA, color=NA, shape=21) + 
		theme(legend.position="none", legend.text =  element_blank()) + 
		scale_y_continuous(breaks=c(0:5), labels=c(as.character(5:0)))
	} else {
		p3 <- ggplot() + 
		geom_point(data = datasetDown, mapping = aes(x = datasetDown$position,
			y = datasetDown$varForAxisYdown, 
			color=datasetDown$annotation
			)) +
		theme(legend.text = element_text(size=8, face="bold"),
			legend.position = 'right') +
		geom_segment(aes(x = datasetDown$position,
			xend = datasetDown$position, y = datasetDown$varForAxisYdown,
			yend =0), linetype = "dotted", colour="gray")+
		scale_x_continuous(limits=c(position$start_position,position$end_position))
	}

	p1 <-p1 + cowplot::theme_cowplot() +
		scale_shape_manual(values =c(16,17), labels = c("No information about LOF","LOF")) +         
		theme(axis.line = axis_custom_up(), axis.line.x = element_blank(), 
			plot.margin=unit(c(1,1,-0.2,1), "cm"), legend.title = element_blank(),
			axis.ticks.x = element_blank(), axis.text.x=element_blank(), 
			axis.title.x = element_blank()) +
		labs(y=paste(databaseDown,"\nVariable: ",varForAxisYdown), x=paste(xAxisAnno))
			
	if(!is.null(exRanges) & !is.null(oldExons)){
		pos_for_ranges <-removeIntronMut(exRanges, data.frame(position=seq(position$start_position,position$end_position)))
		native_pos_whole_gene <-getNativeDataPos(pos_for_ranges, oldExons, exRanges)
	} else {
		native_pos_whole_gene <- data.frame(position=c(position$start_position,position$end_position)) 
	}
			
	plot_title <- paste("Gene:",position$hgnc_symbol,
			"Chromosome:",chrNr,"NCBI_Build:", NCBI_Build,"\n")
	plot_subtitle <- paste("All exons from",native_pos_whole_gene$position[1],"to",
			native_pos_whole_gene$position[length(native_pos_whole_gene$position)],"bp in genome\n")
			
	if(varForAxisYup=="MAF"){
		p1 <- p1 + labs(y=paste(databaseUp,"\nVariable: ln(",varForAxisYup,")"))
	} else {
		p1 <- p1 + labs(y=paste(databaseUp,"\nVariable:",varForAxisYup))
	}
	
	if(!drawPlotly & !leerUp){
		p1 <- p1 + ggtitle(bquote(atop(.(plot_title), atop(.(plot_subtitle), ""))))
		p1 <- scaleAxisUp(datasetUp$varForAxisYup,p1)
	}
	
	p2 <-ggplot() + coord_fixed(ratio = ratio_prot) + 
		geom_rect( mapping=aes(
			xmin = position$start_position, xmax = position$end_position,
			ymin = 0, ymax = 1, fill = 'a')) + 
		theme(axis.line.x = element_blank(), plot.margin=unit(c(-0.2,1,-0.2,1), "cm"), 
			axis.ticks.x = element_blank(), axis.text.x=element_blank(),
			axis.title.x = element_blank(), axis.text.y=element_blank(),
			axis.title.y = element_blank(), panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(), panel.border = element_blank(),
			panel.background = element_blank(), legend.title=element_blank()) + 
		geom_rect(aes(xmin = exonPositions$exon_chrom_start,
			xmax = exonPositions$exon_chrom_end, ymin = 0, ymax = 1, fill = 'b')) +
		scale_fill_manual(values = c("orange", "magenta4"), labels=c("Intron","Exon"))

	p3 <- p3  + cowplot::theme_cowplot() +
		theme(axis.line = axis_custom_down(), axis.line.x = element_blank(), 
			plot.margin=unit(c(-0.2,1,1,1), "cm"), 
			legend.text = element_text(size=8, face="bold"), legend.position = 'right',
			axis.ticks.x = element_blank(), axis.text.x=element_blank(), 
            axis.title.x = element_text(), legend.title = element_blank(),
			legend.position="none") +
		labs(y=paste(databaseDown,"\nVariable: ",varForAxisYdown), x=paste(xAxisAnno))
		
	if(!drawPlotly & !leerDown){
		p3 <- scaleAxisDown(datasetDown$varForAxisYdown,p3)
	}
	
	if(drawPlotly) {
		xstart <-position$start_position
		xend <-position$end_position
		if(!leerUp){
			p4 <-p1 + geom_rect( mapping=aes(
					xmin = xstart, xmax = xend,
					ymin =min(y_up)-1, ymax = min(y_up)-0.7,
					fill = 'a')) + geom_rect(aes(xmin = exonPositions$exon_chrom_start,
					xmax = exonPositions$exon_chrom_end, ymin =min(y_up)-1,
					ymax = min(y_up)-0.7, fill = 'b')) +
				scale_fill_manual(values = c("orange", "magenta4"), labels=c("Intron","Exon")) +
				theme(plot.margin=unit(c(1,1,1,1), "cm"), axis.text.x = element_text(face="bold"),
				 axis.ticks.x=element_line())+
				ggtitle(label = paste("Gene:",position$hgnc_symbol, "Chromosome:",chrNr,"NCBI_Build:", NCBI_Build))
		} else {
			p4 <- p1 + geom_rect( mapping=aes(
					xmin = xstart, xmax = xend,
					ymin = 0, ymax = 0.3, fill = 'a')) + geom_rect(aes(xmin = exonPositions$exon_chrom_start,
					xmax = exonPositions$exon_chrom_end, ymin = 0, ymax = 0.3, fill = 'b')) +
				scale_fill_manual(values = c("orange", "magenta4"), labels=c("Intron","Exon")) +
				ggtitle(label = paste("Gene:",position$hgnc_symbol, "Chromosome:",chrNr,"NCBI_Build:", NCBI_Build)) + 
				theme(axis.text.x = element_text(face="bold"), axis.ticks.x=element_line())
		}
		
		p5 <-p3 + geom_rect( mapping=aes(
				xmin = xstart, xmax = xend,
				ymin = 0, ymax = 0.3, fill = 'a')) + geom_rect(aes(xmin = exonPositions$exon_chrom_start,
				xmax = exonPositions$exon_chrom_end, ymin = 0, ymax = 0.3, fill = 'b')) +
			scale_fill_manual(values = c("orange", "magenta4"), labels=c("Intron","Exon")) +
			ggtitle(label = paste("Gene:",position$hgnc_symbol, "Chromosome:",chrNr,"NCBI_Build:", NCBI_Build)) +
			theme(axis.text.x = element_text(face="bold"), axis.ticks.x=element_line(),axis.ticks.y=element_blank())
		
		gg1 <- ggplotly(p4)
		filenameUp <- paste(position$hgnc_symbol,"chr",chrNr,element,"up",
				varForAxisYup,"html",sep=".")
		htmlwidgets::saveWidget(as.widget(gg1), filenameUp)
	
		gg2 <- ggplotly(p5)
		filenameDown <- paste(position$hgnc_symbol,"chr",chrNr,element,"down",
				varForAxisYdown,"html",sep=".")
		htmlwidgets::saveWidget(as.widget(gg2), filenameDown)
	}		
	
	if(leerUp){
		p1 <- p1+theme(legend.text =  element_blank())
	}
	if(leerDown){
		p3 <- p3+theme(legend.text =  element_blank())
	}
	
	if(!is.null(userCoord)){
		datasetUp <-datasetUp[order(datasetUp$position),]
		idx <- which(datasetUp$position >= userCoord$position[1] &
				datasetUp$position <= userCoord$position[length(userCoord$position)])
		if(length(idx)>0){
			varY <- datasetUp$varForAxisYup[idx[1]:idx[length(idx)]]
			if(varForAxisYup!="MAF"){
				p1 <- scaleAxisUp(varY,p1)
			}
		}
		if(!is.null(oldExons)){
			nativeCoord <- getNativeDataPos(userCoord, oldExons, exRanges)
		} else {
			nativeCoord <- userCoord
		}
		plot_subtitle <- paste("Position from",nativeCoord$position[1],"to",
			nativeCoord$position[length(nativeCoord$position)],"bp\n")
		
		p1 <- p1 + coord_cartesian(xlim = c(userCoord$position[1],
			userCoord$position[length(userCoord$position)]), expand = FALSE) +
		ggtitle(bquote(atop(.(plot_title), atop(.(plot_subtitle), ""))))

		p2 <- p2 + coord_fixed(xlim = c(userCoord$position[1],
			userCoord$position[length(userCoord$position)]), ratio=ratio_prot, expand = FALSE)
			
		datasetDown <-datasetDown[order(datasetDown$position),]
		idx <- which(datasetDown$position >= userCoord$position[1] &
				datasetDown$position <= userCoord$position[length(userCoord$position)])
		if(length(idx)>0){
			varY <- datasetDown$varForAxisYdown[idx[1]:idx[length(idx)]]
			p3 <- scaleAxisDown(varY,p3)
		}
		p3 <- p3 + coord_cartesian(xlim = c(userCoord$position[1],
			userCoord$position[length(userCoord$position)]), expand = FALSE) +
		theme(axis.ticks.x = element_line(), axis.text.x=element_text())
			
		p3 <- scaleX(position$start_position,position$end_position,p3,prot_length, oldExons, exRanges)
	}
	g1 <- ggplotGrob(p1)
	g2 <- ggplotGrob(p2)
	g3 <- ggplotGrob(p3)
	g <-rbind.gtable(g1, g2, g3)
	return(g)
}

scaleAxisDown <-function(vect,p) {
	if(abs(min(vect))+1 < 16){
		p <- p + scale_y_continuous(breaks=c((min(vect)-1):0),
			labels=c(as.character((abs(min(vect))+1):0)), limits=c(min(vect)-1,0))
	} else if(abs(min(vect))+1 < 30){
		p <- p + scale_y_continuous(breaks=seq((min(vect)-1),0,5),
			labels=c(as.character(seq((abs(min(vect))+1),0,-5))), limits=c(min(vect)-1,0))
	} else if (abs(min(vect))+1 < 100){
		p <- p + scale_y_continuous(breaks=seq((min(vect)-1),0,10),
			labels=c(as.character(seq((abs(min(vect))+1),0,-10))), limits=c(min(vect)-1,0))
	} else {
		p <- p + scale_y_continuous(breaks=seq((min(vect)-1),0,50),
			labels=c(as.character(seq((abs(min(vect))+1),0,-50))), limits=c(min(vect)-1,0))
	}
	return(p)
}

scaleAxisUp <-function(vect,p) {
	if(max(vect)+1 < 16){
		p <- p + scale_y_continuous(breaks=c(0:(max(vect)+1)),
			labels=c(as.character(0:(max(vect)+1))), limits=c(0, max(vect)+1))
	} else if(max(vect)+1 < 30){
		p <- p + scale_y_continuous(breaks=seq(0,(max(vect)+1),5),
			labels=c(as.character(seq(0,(max(vect)+1),5))), limits=c(0, max(vect)+1))
	} else if (max(vect)+1 < 100){
		p <- p + scale_y_continuous(breaks=seq(0,(max(vect)+1),10),
			labels=c(as.character(seq(0,(max(vect)+1),10))), limits=c(0, max(vect)+1))
	} else if(max(vect)+1 < 1000){
		p <- p + scale_y_continuous(breaks=seq(0,(max(vect)+1),50),
			labels=c(as.character(seq(0,(max(vect)+1),50))), limits=c(0, max(vect)+1))
	} else {
		p <- p
	}
	return(p)
}

scaleX <-function(st_pos,end_pos,p,len, oldExons, exRanges){
	pos <- seq(st_pos,end_pos,round(len*0.4))

	if(!is.null(oldExons)&!is.null(exRanges)){
		data <-removeIntronMut(exRanges,data.frame(position=pos))
		native_pos <- getNativeDataPos(data, oldExons, exRanges)
		p <- p + scale_x_continuous(breaks=data,
			labels=native_pos,
			limits=c(st_pos,end_pos))
	} else {
		p <- p + scale_x_continuous(breaks=pos,
			labels=pos,
			limits=c(st_pos,end_pos))
	}
	
	return(p)
}
	
getExonRanges <-function(exonPositions) {
	ex <- IRanges(start=exonPositions$exon_chrom_start,
		 end=exonPositions$exon_chrom_end)
	ex <- reduce(ex)
	return(ex)
}

getRelativeExons <-function(exons){
	st_pos <- start(exons)
	width <- width(exons)
	pos <- c(0)
	if(length(st_pos)>1){
		for (i in 2:length(st_pos)){
			pos[i] <-pos[i-1]+width[i-1]+5
		}
	} 
	
	new_exons_pos <- IRanges(start=pos, width=width)
	return(new_exons_pos)
}

removeIntronMut <- function(exons, dataset){
	mutInExons <-findOverlaps(as.integer(dataset$position),
			 exons, select="first")
	dataset["inExon"] <- mutInExons
	dataset <-na.omit(dataset)
	return(dataset)
}

getRelativeDataPos <- function(dataset, oldExonPos, newExonPos){
	for(i in 1:nrow(dataset)){
		dataset$position[i] <- start(newExonPos)[dataset$inExon[i]] + 
		dataset$position[i] - start(oldExonPos)[dataset$inExon[i]]
	}
	return(dataset)
}

getNativeDataPos <- function(dataset, oldExonPos, newExonPos){
	for(i in 1:nrow(dataset)){
		dataset$position[i] <-dataset$position[i] - start(newExonPos)[dataset$inExon[i]] + 
		start(oldExonPos)[dataset$inExon[i]]
	}
	return(dataset)
}

rbind.gtable <- function(..., size = "max", z = NULL){
	gtables <- list(...)
	if (!is.null(z)) {
		gtables <- z_arrange_gtables(gtables, z)
	}
	Reduce(function(x, y) rbind_gtable(x, y, size = size), gtables)
}

rbind_gtable <- function(x, y, size = "max"){
	stopifnot(ncol(x) == ncol(y))
	if (nrow(x) == 0) return(y)
	if (nrow(y) == 0) return(x)
	
	y$layout$t <- y$layout$t + nrow(x)
	y$layout$b <- y$layout$b + nrow(x)
	x$layout <- rbind(x$layout, y$layout)
	
	x$heights <- insert.unit(x$heights, y$heights)
	x$rownames <- c(x$rownames, y$rownames)
	
	size <- match.arg(size, c("first", "last", "max", "min"))
	x$widths <- switch(size,
					first = x$widths,
					last = y$widths,
					min = unit.pmin(x$widths, y$widths),
					max = unit.pmax(x$widths, y$widths)
	)
	
	x$grobs <- append(x$grobs, y$grobs)
	
	return(x)
}

insert.unit <- function(x, values, after = length(x)){
	lengx <- length(x)
	if (lengx == 0) return(values)
	if (length(values) == 0) return(x)
	
	if (after <= 0) {
		unit.c(values, x)
	} else if (after >= lengx) {
		unit.c(x, values)
	} else {
		unit.c(x[1L:after], values, x[(after + 1L):lengx])
	}
}

element_grob.element_custom_up <- function(element, ...){
	grid::segmentsGrob(1,0,1,1, arrow = arrow())
}

axis_custom_up <- function(...){
	structure(
		list(...),
		class = c("element_custom_up","element_blank", "element") 
	) 

}

element_grob.element_custom_down <- function(element, ...){
	grid::segmentsGrob(1,1,1,0, arrow = arrow())
}

axis_custom_down <- function(...){
	structure(
		list(...), 
		class = c("element_custom_down","element_blank", "element") 
	) 

}












