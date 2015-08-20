chooseMsDataParser <- function(wanted_version) {
	if (wanted_version == '1') {
		return(parser_v1);
	}
}

getMsDataVersionId <- function(msdata_version,data) {
	if (msdata_version == '1') {
		metadata = data$metadata
		if ( is.list(metadata) ) {
			metadata = data$metadata[[1]]
		}
		return(paste(sapply(metadata$software,function(soft) { paste(soft$name,soft$version,sep=':',collapse=':') }),collapse=':'))
	}
}

join_composition <- function(composition) {
	paste(composition,collapse=';')
}


parser_v1 <- function(rows,rKeys) {
	row <- rows[[1]]
	do.call("rbind", lapply(rows,function(row) {
		peptide_id <- substring(tempfile(pattern="peptide",tmpdir=''),2)

		base <- data.frame(	peptide=row$sequence,
					peptide.key=peptide_id,
					peptide.start=row$peptide_start,
					peptide.end=row$peptide_start + nchar(row$sequence),
					site=NA,
					site.composition=NA,
					ambiguous=NA,
					multiple.protein=row$multi_protein,
					quantification=NA,
					quantification_mad=NA,
					composition=join_composition(row$composition),
					source=row$source
				)
		if ('sites' %in% names(row)) {
			num_sites <- length(row$sites)
			base <- base[rep(seq_len(nrow(base)), num_sites), ]
			base$site <- sapply(row$sites,function(site) { site[[1]] })
			base$site.composition <- sapply(row$sites,function(site) { site[[2]] })
		}

		if ('ambiguous_mods' %in% names(row)) {
			base$ambiguous <- paste(row$ambiguous_mods,collapse=';')
		}

		if ('quant' %in% names(row)) {
			base$quantification <- row$quant$quant
			if ('mad' %in% names(row$quant)) {
				base$quantification_mad <- row$quant$mad
			}
		}

		base
	}));
}