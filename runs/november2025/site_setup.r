settings_path <- "/projectnb/dietzelab/dongchen/anchorSites/NA_runs/SDA_8k_site/pecan.xml"
multisettings <- PEcAn.settings::read.settings(settings_path)

# ------------------------------------------------------------------------------
# Specify sites:
#   At present, sites are chosen by name and the names are used to fetch the 
#   IDs. This can be generalized when we want to include unnamed sites, in
#   which case sites can be chosen based on ID or geographic coordinates.
# ------------------------------------------------------------------------------

site_info_path <- file.path("/projectnb", "dietzelab", "dongchen", "anchorSites",
                            "NA_runs", "SDA_8k_site", "site_info.Rdata")

site_names <- c("Harvard Forest EMS Tower/HFR1 (US-Ha1)",
                "Bartlett Experimental Forest (US-Bar)")

met_time_horizon <- c(met.start = "2012/01/01",
                      met.end = "2024/12/31")

# ------------------------------------------------------------------------------
# Validate sites exist and fetch site IDs
# ------------------------------------------------------------------------------

site_idx <- match(site_names, site_info$site_name)
if(anyNA(site_idx)) {
  stop("Specified site(s) not found in ", site_info_path)
}

site_ids <- site_info$site_ids[site_idx]


# ------------------------------------------------------------------------------
# Helper function for constructing site block of settings
# ------------------------------------------------------------------------------

construct_site_settings_block <- function(site_id, site_info) {
  site_idx <- match(site_id, site_info$site_ids)
  if(length(site_idx) > 1) {
    stop("construct_site_settings_block(site_id) matches multiple sites.",
         "site_id = ", site_id)
  }
  if(is.na(site_idx)) {
    stop("construct_site_settings_block(site_id) failed to find site info.",
         "site_id = ", site_id)
  }
  
  list(
    id = site_id,
    name = site_info$site_name[site_idx],
    lat = site_info$lat[site_idx],
    lon = site_info$lon[site_idx],
    met.start = met_time_horizon["met.start"],
    met.end =  met_time_horizon["met.end"]
  )
}

