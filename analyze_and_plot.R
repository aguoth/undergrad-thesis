library(magrittr)
options(stringsAsFactors = FALSE);
WGCNA::enableWGCNAThreads()

net_type <- "signed"
cor <- WGCNA::cor

data <-
    readr::read_tsv("DGEx_P1e-3_matrix_log2centered.txt") %>%
    tidyr::pivot_longer(!Probe, names_to="Sample", values_to="Expression") %>%
    tidyr::pivot_wider(names_from=Probe, values_from=Expression) %>%
    dplyr::select(-Sample)

traits <-
    readr::read_csv("traits.csv") %>%
    dplyr::select(-sample)

sft = WGCNA::pickSoftThreshold(data, RsquaredCut = 0.85, networkType = net_type)

net <- 
    data %>%
    WGCNA::adjacency(
        power=sft$powerEstimate,
        type=net_type
    )

inet <-
    net %>%
    igraph::graph_from_adjacency_matrix(weighted=TRUE,mode="directed",diag=FALSE)

mods <-
    data %>%
    WGCNA::blockwiseModules(
        power=sft$powerEstimate,
        networkType=net_type
    )

net_params <- 
    data %>%
    WGCNA::networkConcepts(
        power=sft$powerEstimate,
        networkType=net_type
    )

tom <-
    net %>%
    WGCNA::TOMdist()

dendro <-
    tom %>%
    as.dist() %>%
    flashClust::flashClust()

eg <-
    data %>%
    WGCNA::moduleEigengenes(mods$colors)

cols <-
    mods$colors %>%
    unique()

sizes <-
    cols %>% 
    purrr::map_int(function(col) sum(col == unname(mods$colors))) %>%
    setNames(cols) %>%
    tibble::enframe(name="Module", value="Size")

hubs <-
    data %>%
    WGCNA::chooseTopHubInEachModule(
        mods$colors,
        power=sft$powerEstimate,
        type=net_type
    ) %>%
    tibble::enframe(name="Module", value="Hub Node")

edge_colour_lookup <-
    cols %>%
    purrr::map(
        function(colour) {
            names(mods$colors[mods$colors == colour]) %>%
            igraph::induced_subgraph(inet,.) %>% 
            igraph::E() %>% 
            igraph::as_ids() %>%
            setNames(rep(colour,length(.)),.)
        }
    ) %>%
    unlist()

mod_trait_cor_raw <-
    cor(eg$eigengenes %>% WGCNA::orderMEs(), traits, use = "p")

mod_trait_cor <-
    mod_trait_cor_raw %>%
    tibble::as_tibble(rownames="eigengene")

mod_trait_p <-
    WGCNA::corPvalueStudent(mod_trait_cor_raw, traits %>% nrow()) %>%
    tibble::as_tibble(rownames="eigengene")

# Analysis
dplyr::inner_join(sizes,hubs)%>%
    readr::write_csv("module_properties.csv")

net_params$Summary %>%
    tibble::as_tibble(rownames="Parameter") %>%
    readr::write_csv("network_properties.csv")

eg$eigengenes %>%
    tibble::as_tibble() %>%
    readr::write_csv("module_eigengenes.csv")

net %>% 
    tibble::as_tibble(rownames="Probe") %>% 
    readr::write_csv("correlation_adjacency_net.csv")

mod_trait_cor %>%
    readr::write_csv("eigengene_trait_correlation.csv")

mod_trait_p %>%
    readr::write_csv("eigengene_trait_pvalues.csv")

dplyr::inner_join(
    mod_trait_cor %>%
    tidyr::pivot_longer(!eigengene,  names_to="trait", values_to="cor"),
    mod_trait_p %>%
    tidyr::pivot_longer(!eigengene,  names_to="trait", values_to="p_val") 
) %>%
dplyr::mutate(abs_cor = abs(cor)) %>%
dplyr::group_by(eigengene) %>%
dplyr::filter(abs_cor == max(abs_cor),eigengene != "MEgrey") %>%
dplyr::select(!abs_cor) %>%
dplyr::ungroup() %>%
readr::write_csv("eigengene_best_trait_by_cor.csv")

# Figures
png("tom_plot.png")
tom %>% WGCNA::TOMplot(dendro)
dev.off()

png("dendro_plot.png")
dendro %>% WGCNA::plotDendroAndColors(mods$colors)
dev.off()

mds <-
    tom %>%
    as.dist() %>%
    cmdscale()

mds_plot <-
    tibble::tibble(SD1=mds[,1],SD2=mds[,2]) %>%
    dplyr::mutate(
        module = mods$colors %>% unname()
    ) %>%
    ggplot2::ggplot() +
        ggplot2::geom_point(ggplot2::aes(x=SD1, y=SD2, colour=module)) +
        ggplot2::scale_color_manual(breaks=mods$colors %>% unique() %>% sort (), values=mods$colors %>% unique() %>% sort ())

"mds_plot.png" %>%
    ggplot2::ggsave(plot=mds_plot)

cor_raster_plot <-
    mod_trait_cor %>%
    dplyr::filter(eigengene != "MEgrey") %>%
    tidyr::pivot_longer(!eigengene,names_to="class", values_to="correlation") %>%
    dplyr::mutate(
        class = forcats::fct_relevel(class, "nymph", "soldier", "worker", "boston", "raleigh", "toronto"),
        eigengene = forcats::fct_relevel(eigengene, "MEyellow", "MEbrown", "MEturquoise", "MEblue")
    ) %>%
    ggplot2::ggplot(ggplot2::aes(class, eigengene, fill = correlation)) + 
        ggplot2::geom_raster() +
        ggplot2::scale_fill_gradient2(
            low = scales::muted("red"),
            mid = "white",
            high = scales::muted("green"),
            midpoint = 0
        )

"cor_raster_plot.png" %>%
    ggplot2::ggsave(plot=cor_raster_plot)

degree_dist_plot <-
    igraph::strength(inet, mode="out") %>%
    unname() %>%
    tibble::enframe(name = NULL, value = "degree") %>%
    ggplot2::ggplot(ggplot2::aes(x=degree)) + 
        ggplot2::geom_histogram(ggplot2::aes(y=..density..), colour="black", fill="white")+
        ggplot2::geom_density()

"degree_dist_plot.png" %>%
    ggplot2::ggsave(plot=degree_dist_plot)

prop_var_plot <-
    eg$varExplained[1,] %>%
    unlist() %>%
    setNames(eg$eigengenes %>% colnames()) %>%
    tibble::enframe(name="eigengene", value="variance") %>%
    dplyr::filter(eigengene != "MEgrey") %>%
    ggplot2::ggplot() +
        ggplot2::geom_point(ggplot2::aes(x=eigengene, y=variance, color=eigengene, size=2)) +
        ggplot2::scale_color_manual(values=eg$eigengenes %>% colnames() %>% stringr::str_remove("ME") %>% purrr::keep(function(s) s != "grey")) +
        ggplot2::theme(legend.position="none", panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
        ggplot2::coord_fixed(ratio = 15)

"prop_var_plot.png" %>%
    ggplot2::ggsave(plot=prop_var_plot)

inet_pruned <- 
    inet %>%
    igraph::delete_edges(which(igraph::E(.)$weight < 0.3)) %>%
    igraph::as.undirected(mode="collapse")

inet_pruned_with_attr <-
    inet_pruned %>%
	igraph::set_vertex_attr("label", value="") %>%
	igraph::set_vertex_attr("size", value=
		inet_pruned %>%
		igraph::degree(normalized=TRUE) %>%
		magrittr::multiply_by(5)
	) %>%
	igraph::set_vertex_attr("color", value=
		inet_pruned %>% 
        igraph::V() %>% 
        igraph::as_ids() %>% 
        {mods$colors[.]}
	) %>%
	igraph::set_vertex_attr("frame.color", value=NA) %>%
	igraph::set_edge_attr("color", value=
		inet_pruned %>%
		igraph::E() %>%
		igraph::as_ids() %>%
		{edge_colour_lookup[.]} %>%
        purrr::map_chr(function(col) scales::alpha(ifelse(is.na(col),"gray70",col),0.1))
	) %>%
	igraph::set_edge_attr("width", value=
        inet_pruned %>%
		{igraph::E(.)$weight} %>%
		magrittr::multiply_by(0.5)
    ) %>%
	igraph::set_graph_attr("layout", igraph::layout_with_kk)

png(file="net_plot.png")
inet_pruned_with_attr %>% plot() %>% print()
dev.off()