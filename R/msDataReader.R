
hexnac.ratios = function(msdata) {
	all_ratios = attributes(msdata)$hexnac.ratios
	results = do.call(rbind,lapply(names(all_ratios),function(specid) {
		ratios = suppressWarnings(as.numeric(unlist(strsplit(all_ratios[[specid]],',',fixed=T))))
		data.frame(spectra=specid,ratio=ratios)
	}))
	results[results$spectra %in% msdata$spectra,]
}

quant.areas = function(msdata) {
	all_areas = attributes(msdata)$quant.areas
	areas_list = unlist(unlist(all_areas,recursive = F),recursive=F)
	results = do.call(rbind,lapply(names(areas_list),function(spec_channel) {
		areas = areas_list[[spec_channel]]
		areas[areas == 0] = NA
		id_parts = unlist(strsplit(spec_channel,'.',fixed=T))
		data.frame(spectra=id_parts[1],channel=id_parts[2],area=sum(areas))
	}))
	results[results$spectra %in% msdata$spectra,]
}

peptide.identifiers = function(msdata) {
	ids = with(msdata,
		lapply(split(
			ifelse(is.na(site),composition,paste(site,site.composition,sep='-')),
			as.factor(spectra)
			),
			function(items) { paste(sort(unique(items)),collapse=';')  }
		)
	)
	data.frame(spectra=msdata$spectra,identifier=paste(msdata$peptide,ids[msdata$spectra]))
}


chooseMsDataParser <- function(wanted_version) {
	if (wanted_version == "1") {
		wanted_version = '1.0'
	}
	version_info = package_version(wanted_version)
	if (version_info <= '1.2.999') {
		return(parser_v1);
	}
	if (version_info <= '1.3.999') {
		return(parser_v1_3);
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


parser_v1 <- function(rows,rKeys,attribs) {
	quants = list()
	hexnac_ratios = list()
	result = as.data.frame(do.call("rbind", lapply(rows,function(row) {
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
			if (is.null(quants[[base$spectra]])) {
				quants[[base$spectra]] <<- list()
			}
			quants[[base$spectra]][[length(quants[[base$spectra]])+1]] <<- row$quant_areas
			hexnac_ratios[[base$spectra]] <<- row$hexnac_ratio
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
	attribs$quant.areas = c(attribs$quant.areas, quants)
	attribs$hexnac.ratios = c(attribs$hexnac.ratios, hexnac_ratios)
	result
}

parser_v1_3 <- function(rows,rKeys,attribs) {
	quants = list()
	hexnac_ratios = list()
	result = as.data.frame(do.call("rbind", lapply(rows,function(row) {
		peptide_id <- substring(tempfile(pattern="peptide",tmpdir=''),2)
		base = NULL
		base_template <- data.frame(	peptide=row$sequence,
					peptide.key=peptide_id,
					peptide.start=row$peptide_start,
					peptide.end=row$peptide_start + nchar(row$sequence),
					site=NA,
					site.composition=NA,
					ambiguous.site.start=NA,
					ambiguous.site.end=NA,
					ambiguous.site.composition=NA,
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
			base_template$spectra <- digest::digest( join_composition(sapply( row$spectra[ order(sapply(row$spectra,function(x) { x$rt })) ], spectrum_key )),"crc32")
			if (is.null(quants[[base_template$spectra]])) {
				quants[[base_template$spectra]] <<- list()
			}
			quants[[base_template$spectra]][[length(quants[[base_template$spectra]])+1]] <<- row$quant_areas
			hexnac_ratios[[base_template$spectra]] <<- row$hexnac_ratio
		}
		if ('sites' %in% names(row) && length(row$sites) > 0) {
			num_sites <- length(row$sites)
			base <- base_template[rep(seq_len(nrow(base_template)), num_sites), ]
			base$site <- sapply(row$sites,function(site) { site[[1]] })
			base$site.composition <- sapply(row$sites,function(site) { site[[2]] })
		}
		if ('ambiguous_sites' %in% names(row) && length(row$ambiguous_sites) > 0) {
			num_sites <- length(row$ambiguous_sites)
			base_ambig <- base_template[rep(seq_len(nrow(base_template)), num_sites), ]
			base_ambig$ambiguous.site.start <- sapply(row$ambiguous_sites,function(site) { site[[1]][1] })
			base_ambig$ambiguous.site.end <- sapply(row$ambiguous_sites,function(site) { site[[1]][2] })
			base_ambig$ambiguous.site.composition <- sapply(row$ambiguous_sites,function(site) { site[[2]] })
			if (is.null(base) || nrow(base) > 0) {
				base = rbind(base,base_ambig)
			} else {
				base = base_ambig
			}
		}
		if (is.null(base)) {
			base = base_template
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
	attribs$quant.areas = c(attribs$quant.areas, quants)
	attribs$hexnac.ratios = c(attribs$hexnac.ratios, hexnac_ratios)
	result
}
