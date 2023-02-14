## edit
plot_comutation <- function(dat.integrate, out.dir="results/", cache=T, frequent=T, Reordering=NULL, cutoff.frequent=1)
{
#	cex_heatmap_axis = 0.66
#	cex_legend_axis = 0.66
	cex_legend_axis = 1.5
	cex_heatmap_axis = 1.5
	image_width = 9
	image_height = 9
        system(paste("mkdir -p", paste0(out.dir,"/cache")))
        #dat.integrate <- dat.integrate.impute

        # cache <- F
        ## TODO: make a real function
        # - especially automatizing the color scale variation according to the range of the plot

        if (frequent)
        {
                dat.integrate <- dat.integrate[, colSums(dat.integrate, na.rm=T)>=cutoff.frequent*nrow(dat.integrate)/100]
        }

        ## Reorder dat.integrate
        if (!is.null(Reordering))
        {
                dat.integrate <- dat.integrate[, Reordering$order]
        }

        # TODO: GROUPING OF THE FEATURES
        genomicGroups <- factor(colnames(dat.integrate),levels=colnames(dat.integrate))

        ##
        if (!cache)
        {
                logPInt <- sapply(1:ncol(dat.integrate), function(i) sapply(1:ncol(dat.integrate), function(j) {f<- try(fisher.test(dat.integrate[,i], dat.integrate[,j]), silent=TRUE); if(class(f)=="try-error") 0 else ifelse(f$estimate>1, -log10(f$p.val),log10(f$p.val))} ))
                odds <- sapply(1:ncol(dat.integrate), function(i) sapply(1:ncol(dat.integrate), function(j) {f<- try(fisher.test(table(dat.integrate[,i], dat.integrate[,j])), silent=TRUE); if(class(f)=="try-error") f=NA else f$estimate} ))

                # TODO: Verify
                pairs <- sapply(1:ncol(dat.integrate), function(i) colMeans(dat.integrate * dat.integrate[,i], na.rm=TRUE))

                #
                save(logPInt, file=paste0(out.dir,"/cache/logPInt.RData"))
                save(odds, file=paste0(out.dir,"/cache/odds.RData"))
                save(pairs, file=paste0(out.dir,"/cache/pairs.RData"))

        } else
        {
                logPInt <- get(load(paste0(out.dir,"/cache/logPInt.RData")))
                odds <- get(load(paste0(out.dir,"/cache/odds.RData")))
                pairs <- get(load(paste0(out.dir,"/cache/pairs.RData")))

        }

        diag(logPInt) <- 0
        diag(odds) <- 1
        colnames(odds) <- rownames(odds) <- colnames(logPInt) <- rownames(logPInt) <- colnames(dat.integrate)
        # Thresholding the odd ratios to stand between -4 and 4
        odds[odds<1e-3] = 1e-4
        odds[odds>1e3] = 1e4
        odds[10^-abs(logPInt) > 0.05] = 1
        logOdds=log10(odds)
        diag(logPInt) <- NA

        ########### HEATMAP ALL 
        library(RColorBrewer)
        library(gplots)

        if (is.null(Reordering))
        {
                pdf(paste0(out.dir,"/comutation_pattern.pdf"), width=image_width, height=image_height, pointsize=8) 
        } else
        {
                pdf(paste0(out.dir,"/comutation_pattern_reordering.pdf"), width=9, height=9, pointsize=8) 
        }

        par(bty="n", mgp = c(2,.5,0), mar=c(6,6,1,8)+.1, las=2, tcl=-.33)
        ix = TRUE#colnames(interactions) %in% colnames(all_genotypes)
        d <- dist(t(dat.integrate[,ix]) + 10*as.numeric(genomicGroups))
        h = hclust(d, method="average")
        o = order(genomicGroups,-colSums(dat.integrate, na.rm=TRUE))#order(cmdscale(d, k=1))#h$order #c(h$order,(length(h$order) +1):ncol(interactions))
        M <-  matrix( NA, ncol=ncol(odds), nrow=nrow(odds))

        ####
        M[lower.tri(M)] <- cut(logOdds[o,o][lower.tri(M)], breaks = c(-4:0-.Machine$double.eps,0:4), include.lowest=TRUE)  # TODO: Look at how many points
        # Breaks is size 10

        # tmp <- logOdds[o,o][lower.tri(M)]

        # toto <-M[lower.tri(M)]
        # summary(toto)

        # M[upper.tri(M, diag=TRUE)] <- as.numeric(cut(pairs[o,o][upper.tri(M, diag=TRUE)]*nrow(dat.integrate), breaks=c(-1,0,5,10,20,50,100,200,600))) + 9 # TODO: Look at how many points
        M[upper.tri(M, diag=TRUE)] <- as.numeric(cut(pairs[o,o][upper.tri(M, diag=TRUE)]*nrow(dat.integrate), breaks=c(-1,0,5,10,20,50,100,200,max(pairs[o,o][upper.tri(M, diag=TRUE)]*nrow(dat.integrate))+1))) + 9 # TODO: Look at how many points

        # toto <- M[upper.tri(M, diag=TRUE)]
        # tmp <- pairs[o,o][upper.tri(M, diag=TRUE)]*nrow(dat.integrate)

        # Size is 7

        ## Need to print max(M) then adjust colours accordingly

        colorset <- c(
                      brewer.pal(9, "RdBu"), #BrBG
                      "white",
                      brewer.pal(max(na.omit(M[upper.tri(M,diag=T)]))-10, "Purples")) #Blues

        # print(max(M))
        # image(x=1:ncol(logPInt), y=1:nrow(logPInt), M, col=c(brewer.pal(9,"BrBG"), c("white",brewer.pal(7,"Blues"))), breaks=0:max(M,na.rm=TRUE), xaxt="n", yaxt="n", xlab="",ylab="", xlim=c(0, ncol(logPInt)+3), ylim=c(0, ncol(logPInt)+3))
        image(x=1:ncol(logPInt), y=1:nrow(logPInt), M, col=colorset, breaks=0:max(M,na.rm=TRUE), xaxt="n", yaxt="n", xlab="",ylab="", xlim=c(0, ncol(logPInt)+3), ylim=c(0, ncol(logPInt)+3))
        #image(x=1:ncol(logPInt), y=1:nrow(logPInt), M, col=c(brewer.pal(8,"BrBG"), c("white",brewer.pal(7,"Blues"))), breaks=0:max(M,na.rm=TRUE), xaxt="n", yaxt="n", xlab="",ylab="", xlim=c(0, ncol(logPInt)+3), ylim=c(0, ncol(logPInt)+3))
        l <- colnames(logPInt)[o]

        if (is.null(Reordering))
        {
                mtext(side=1, at=1:ncol(logPInt), l, cex=cex_heatmap_axis, font=ifelse(grepl("^[A-Z]",l),3,1))
                mtext(side=2, at=1:ncol(logPInt), l, cex=cex_heatmap_axis, font=ifelse(grepl("^[A-Z]",l),3,1))
        } else
        {
                mtext(side=1, at=1:ncol(logPInt), l, cex=cex_heatmap_axis, font=ifelse(grepl("^[A-Z]",l),3,1), col=Reordering$color.features)
                mtext(side=2, at=1:ncol(logPInt), l, cex=cex_heatmap_axis, font=ifelse(grepl("^[A-Z]",l),3,1), col=Reordering$color.features)
        }



        abline(h=0:ncol(logPInt)+.5, col="white", lwd=.5)
        abline(v=0:ncol(logPInt)+.5, col="white", lwd=.5)
        P <- 10^-abs(logPInt[o,o])
        P[upper.tri(P)] <- NA
        w = arrayInd(which(p.adjust(P, method="BH") < .1), rep(nrow(logPInt),2))
        points(w, pch=".", col="salmon3", cex=3)
        w = arrayInd(which(p.adjust(P) < .05), rep(nrow(logPInt),2))
        points(w, pch="*", col="salmon3", cex=2)

        # Frequency position
        # image(y = 1:9 +18, x=rep(ncol(logPInt),2)+c(2,3), z=matrix(c(1:8), nrow=1), col=c("white",brewer.pal(7,"Blues")), add=TRUE)
        # axis(side = 4, at = seq(1,7) + 19, cex.axis=cex_legend_axis, tcl=-.15, label=c(1,5,10,20,50,100,200), las=1, lwd=.5)
        # mtext(side=4, at=28, "Frequency", las=2, line=-1,cex=cex_legend_axis)

        image(y = 1:9 +12, x=rep(ncol(logPInt),2)+c(2,3), z=matrix(c(1:8), nrow=1), col=c("white",brewer.pal(max(na.omit(M[upper.tri(M,diag=T)]))-10, "Purples")), add=TRUE)
        axis(side = 4, at = seq(1,7) + 13, cex.axis=cex_legend_axis, tcl=-.15, label=c(1,5,10,20,50,100,200), las=1, lwd=.5)
        mtext(side=4, at=22, "Frequency", las=2, line=-1,cex=cex_legend_axis)

        # Odd-ratio position
        # image(y = 1:8 +5, x=rep(ncol(logPInt),2)+c(2,3), z=matrix(c(1:8), nrow=1), col=brewer.pal(8,"BrBG"), add=TRUE)
        # axis(side = 4, at = seq(1,7) + 5.5, cex.axis=cex_legend_axis, tcl=-.15, label=10^seq(-3,3), las=1, lwd=.5)
        # mtext(side=4, at=14, "Odds ratio", las=2, line=-1,cex=cex_legend_axis)

        image(y = 1:8 +4, x=rep(ncol(logPInt),2)+c(2,3), z=matrix(c(1:8), nrow=1), col=brewer.pal(9,"RdBu"), add=TRUE)
        axis(side = 4, at = seq(1,7) + 4.5, cex.axis=cex_legend_axis, tcl=-.15, label=10^seq(-3,3), las=1, lwd=.5)
        mtext(side=4, at=13, "Odds ratio", las=2, line=-1,cex=cex_legend_axis)


        # Significance position
        # mtext(side=4, at=4, "Significance", las=2, line=-1,cex=cex_legend_axis)
        # points(x=rep(ncol(logPInt),2)+2.5, y=1:2, pch=c("*","."))
        # image(x=rep(ncol(logPInt),2)+c(2,3), y=(2:3) +0.5, z=matrix(1), col=brewer.pal(3,"BrBG"), add=TRUE)
        # mtext(side=4, at=3:1, c("P > 0.05", "FDR < 0.1", "FWER < 0.05"), cex=cex_legend_axis, line=0.2)

        mtext(side=4, at=4, "Significance", las=2, line=-1, cex=cex_legend_axis)
        points(x=rep(ncol(logPInt),2)+2.5, y=1:2, pch=c("*","."), cex=c(2,3))
        image(x=rep(ncol(logPInt),2)+c(2,3), y=(2:3)+0.5, z=matrix(1), col=brewer.pal(3,"BrBG"), add=TRUE)
        mtext(side=4, at=3:1, c("P > 0.05", "FDR < 0.1", "FWER < 0.05"), cex=cex_legend_axis, line=0.2)

        ## Add Rectangles
        if (!is.null(Reordering))
        {
                rect(Reordering$groups.begin,Reordering$groups.begin,Reordering$groups.end,Reordering$groups.end, border=Reordering$cols, lwd=5)
        }

        dev.off()

}
