#' Creates a set of diagnostic plots for Bayesian models run on lists of imputed data sets
#' 
#' @param model a \code{brms.fit} object
#' 
#' @param data_list a \code{list} containing the imputed data sets used to fit the model
#' 
#' @param n_samples a \code{integer} indicating the number of posterior samples to used for the ppc_density plot
#' 
#' @importFrom brms pp_check
#' @importFrom cowplot plot_grid
#' @importFrom ggplot2 labs
#' @importFrom glue glue

create_diagnostic_plots <- function(model, impute_list, n_samples, dir_path, height = 6, width = 9, device = "png", 
                                    dpi = 600) {
    model_name <- deparse(substitute(model))
    ppc_dens_plotlist <- list() 
    ppc_hist_plotlist <- list()
    for(i in seq_along(impute_list)) {
        # Posterior density plots
        ppc_dens_plotlist[[i]] <- brms::pp_check(
            model, 
            newdata = impute_list[[i]], 
            nsamples = n_samples
        )
        
        # Poster residual plots
        ppc_hist_plotlist[[i]] <- brms::pp_check(
            model, 
            newdata = impute_list[[i]], 
            nsamples = 1, 
            type = "error_hist"
        )
    }
    
    dens_plot <- cowplot::plot_grid(
        plotlist = ppc_dens_plotlist, 
        align = "hv", 
        axis = "tblr",
        labels = paste("m =", seq_along(impute_list))
    ) 
    
    hist_plot <- cowplot::plot_grid(
        plotlist = ppc_hist_plotlist, 
        align = "hv", 
        axis = "tblr",
        labels = paste("m =", seq_along(impute_list)) 
    ) 
    
    ggplot2::ggsave(
        filename = glue::glue("{dir_path}/{model_name}_dens.{device}"), 
        plot = dens_plot, 
        device = device, 
        height = height, 
        width = width, 
        units = "in", 
        dpi = dpi
    )
}