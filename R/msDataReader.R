chooseMsDataParser <- function(wanted_version) {
	if (wanted_version == "1") {
		wanted_version = '1.0'
	}
	version_info = package_version(wanted_version)
	if (version_info <= '1.1.999') {
		return(parser_v1);
	}
}

getMsDataVersionId <- function(msdata_version,data) {
	if (msdata_version == "1") {
		msdata_version = '1.0'
	}
	version_info = package_version(msdata_version)
	if (version_info <= '1.1.999') {
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

spectrum_key <- function(spectrum) {
	if (is.null(spectrum$ppm)) {
		spectrum$ppm = Inf
	}
	paste(c(round(spectrum$rt,3),spectrum$scan,spectrum$charge,round(spectrum$ppm,2),round(spectrum$score,2)),collapse='|')
}


parser_v1 <- function(rows,rKeys) {
	row <- rows[[1]]
	as.data.frame(do.call("rbind", lapply(rows,function(row) {
		peptide_id <- substring(tempfile(pattern="peptide",tmpdir=''),2)
		base <- data.frame(	peptide=row$sequence,
					peptide.key=peptide_id,
					peptide.start=row$peptide_start,
					peptide.end=row$peptide_start + nchar(row$sequence),
					site=NA,
					site.composition=NA,
					ambiguous=NA,
					multiple.protein=ifelse(is.null(row$multi_protein),FALSE,row$multi_protein),
					quantification=NA,
					quantification_mad=NA,
					quant_confidence='high',
					composition=join_composition(row$composition),
					site_ambiguity=NA,
					source=row$source,
					spectra=NA
				)
		if ('spectra' %in% names(row)) {
			base$spectra <- digest::digest( join_composition(sapply( row$spectra[ order(sapply(row$spectra,function(x) { x$rt })) ], spectrum_key )),"crc32")
		}
		if ('sites' %in% names(row)) {
			num_sites <- length(row$sites)
			base <- base[rep(seq_len(nrow(base)), num_sites), ]
			base$site <- sapply(row$sites,function(site) { site[[1]] })
			base$site.composition <- sapply(row$sites,function(site) { site[[2]] })
		}

		if ('ambiguous_mods' %in% names(row)) {
			base$ambiguous <- paste(row$ambiguous_mods,collapse=';')
		}

		if ('made_ambiguous' %in% names(row)) {
			base$site_ambiguity = row[['made_ambiguous']]
		}
		if ('quant' %in% names(row)) {
			base$quantification <- row$quant$quant
			if ('mad' %in% names(row$quant)) {
				base$quantification_mad <- row$quant$mad
			}
			if ('singlet_confidence' %in% names(row$quant)) {
				base$quant_confidence <- row$quant$singlet_confidence
			}
		}

		base
	})));
}