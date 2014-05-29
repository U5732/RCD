# svc (milsav@gmail.com)

library("igraph")

# cp - community partition
print_info <- function(cp) {
	cat("Algorithm: ", algorithm(cp), 
          "\nNumber of communities: ", length(cp), 
	    "\nModularity: ", modularity(cp), 
          "\nCommunity sizes: \n", sort(sizes(cp)), "\n\n")
}

# p    - community partitions
# attr - normalized mutual information, rand, adjusted-rand, etc.
compare_communities <- function(p, attr) {
	cat("Comparing communities using ", attr, "\n")
	for (i in 2:length(p)) {
		for (j in 1:(i-1)) {	
			sim <- compare(p[[i]], p[[j]], attr);
			cat(algorithm(p[[i]]), " -- ",
                      algorithm(p[[j]]), ", sim = ", sim, "\n")
		}
	}
	cat("\n\n")
}

# g - graph
# p - community partitions
community_quality <- function(g, p, debug=FALSE) {
 	comms <- list();
	length(comms) <- length(p)
	for (i in 1:length(V(g))) {
		c <- membership(p)[i]
		comms[[c]] <- append(comms[[c]], list(i))
	}

	cat("Algorithm -- ", algorithm(p), "\n")	

	for (i in 1:length(p)) {
		if (debug) {
			cat("Community ", i, ", size: ", length(comms[[i]]), "\n")
		}
		size <- length(comms[[i]])

		rStrong <- TRUE
		totalIntraWeight <- 0
		totalInterWeight <- 0

		for (j in 1:length(comms[[i]])) {
			node_id <- comms[[i]][[j]]
			intraW <- 0
			interW <- 0			

			current_node <- V(g)[node_id]
			neis <- neighbors(g, current_node)

			# cat("Current node -- ", current_node$id, "\n")
			for (k in 1:length(neis)) {
				n <- neis[k]
				e <- get.edge.ids(g, c(node_id, n))
				# cat(V(g)[n]$id, " --, w = ", E(g)[e]$weight, "\n")				
				w <- E(g)[e]$weight
				if (membership(p)[current_node] == membership(p)[n]) {
					intraW <- intraW + w	
				} 
				else {
					interW <- interW + w
				}
			}

			if (debug) {
				cat(current_node$id, ", InterW: ", interW, ", IntraW: ", intraW, "\n")
			}

			if (interW >= intraW) {
				rStrong <- FALSE
			}

			totalIntraWeight <- totalIntraWeight + intraW
			totalInterWeight <- totalInterWeight + interW
		}

		totalIntraWeight <- totalIntraWeight / 2
		conductance <- totalInterWeight / (totalIntraWeight + totalInterWeight)

		rWeak = totalIntraWeight > totalInterWeight		

		cat("Community ", i, ", size = ", size, ", CON: ", conductance, ", RStrong = ", rStrong, 
                "RWeak = ", rWeak, "\n") 
	}
}

# p - community partitions
# g - graph
export <- function(p, g, out_file) {
	sink(file=out_file, append=FALSE)

	cat("NUM,ID,")
	for (j in 1:length(p)) {
		cat(algorithm(p[[j]]))
		if (j < length(p)) {
			cat(",")
		}
	}
	cat("\n")

	for (i in 1:length(V(g))) {
      	cat(V(g)[i]) 
		cat(",") 
		cat(V(g)[i]$id) 
		cat(",")
		for (j in 1:length(p)) {
			cat(membership(p[[j]])[i])
			if (j < length(p)) {
				cat(",")
			}
		}
		cat("\n")
	}

	sink()
}

perform_community_detection <- function(input_file, file_format, out_file) {	
	cat("Community detection for ", input_file, "\n")
	g <- read.graph(input_file, file_format)
	
	c_fg <- fastgreedy.community(g)
	c_wt <- walktrap.community(g)
	c_gn <- edge.betweenness.community(g)
	c_im <- infomap.community(g)
	c_lp <- label.propagation.community(g)
	c_ei <- leading.eigenvector.community(g)
	c_ml <- multilevel.community(g)

	p <- list(c_fg, c_wt, c_gn, c_im, c_lp, c_ei, c_ml)

	for (i in 1:length(p)) {
		print_info(p[[i]])
	}	

	cat("Performing comparison... \n\n")
	compare_communities(p, "rand")	
	compare_communities(p, "adjusted.rand")
	compare_communities(p, "nmi")
	compare_communities(p, "vi")

	cat("Community quality... \n\n")
	for (i in 1:length(p)) {
		community_quality(g, p[[i]], debug=FALSE)
	}

	cat("Community detection finished...\n")
	export(p, g, out_file)
}

perform_community_detection("Comp_0.net", "pajek", "Comp_0_nodes.dat")