psicomA <- function(model, bd, ordered = FALSE, stdlv = TRUE, vd = TRUE, sw = NULL, ...) {

  # ============================================================================
  # DOCUMENTAÇÃO RÁPIDA
  # ============================================================================
  # model   : modelo lavaan (caracter ou objeto)
  # bd      : base de dados (data.frame)
  # ordered : lógico — variáveis ordinais?
  # stdlv   : lógico — standardizar variáveis latentes?
  # vd      : lógico — calcular validade discriminante?
  # sw      : pesos amostrais (opcional)
  # ============================================================================
  
  # ---- 0. INSTALAÇÃO DO PACMAN (se necessário) ---------------------------------
  if (!requireNamespace("pacman", quietly = TRUE)) {
    install.packages(
      "pacman", 
      repos = "https://cloud.r-project.org/",
      dependencies = TRUE,
      quiet = TRUE
    )
  }
  
    
  pacman::p_load(
    lavaan, semTools, dplyr, skimr, psych,
    knitr, kableExtra, DiagrammeR
  )
  
  #--------------------------------------------------
  # DESCRIPTIVOS
  #--------------------------------------------------
  myskim <- skim_with(
    numeric = sfl(
      n  = ~sum(!is.na(.)),
      sk = ~psych::skew(., na.rm = TRUE),
      ku = ~psych::kurtosi(., na.rm = TRUE)
    ),
    append = TRUE
  )
  
  descr <- bd |> myskim() |> as.data.frame()
  names(descr) <- gsub("numeric.","",names(descr))
  names(descr) <- gsub("skim_","",names(descr))
  descr <- descr |> dplyr::select(variable,type,n,everything(),-hist,hist)
  
  descT <- knitr::kable(
    descr,
    caption="Item Descriptives",
    digits=3
  ) |> kableExtra::kable_classic() |> kableExtra::kable_styling(full_width=TRUE)
  
  #--------------------------------------------------
  # CFA
  #--------------------------------------------------
  fit <- if(ordered){
    lavaan::sem(
      model=model,
      data=bd,
      std.lv=stdlv,
      ordered=TRUE,
      parameterization="theta",
      sampling.weights=sw
    )
  } else {
    lavaan::sem(
      model=model,
      data=bd,
      std.lv=stdlv,
      estimator="MLR",
      sampling.weights=sw
    )
  }
  
  #--------------------------------------------------
  # LOADINGS
  #--------------------------------------------------
  tablePesos <- parameterEstimates(fit, standardized=TRUE) |>
    filter(op=="=~") |>
    mutate(
      stars=dplyr::case_when(
        pvalue<.001~"***",
        pvalue<.01~"**",
        pvalue<.05~"*",
        TRUE~""
      )
    ) |>
    select(
      "Latent Factor"=lhs,
      Indicator=rhs,
      Std.Loadings=std.all,
      SE=se,
      Z=z,
      sig=stars
    )
  
  TabLoad <- knitr::kable(
    tablePesos,
    caption="Factor Loadings for measurement model",
    digits=3
  ) |> kableExtra::kable_classic() |> kableExtra::kable_styling(full_width=TRUE)
  
 
  #--------------------------------------------------
  # DETECTA Fatores de 1ª, 2ª e 3ª ordem
  #--------------------------------------------------
  ptable <- parTable(fit)
  meas <- ptable[ptable$op=="=~",c("lhs","rhs")]
  latent <- unique(meas$lhs)
  latent_indicators <- meas[meas$rhs %in% latent,]
  higher <- unique(latent_indicators$lhs)
  first_order <- setdiff(latent,higher)
  second_order <- setdiff(higher, latent_indicators$lhs[latent_indicators$rhs %in% higher])
  third_order <- setdiff(higher,second_order)
  
  #--------------------------------------------------
  # DIAGRAMA SEM (ajuste visual avançado)
  #--------------------------------------------------
  std_sol <- standardizedSolution(fit)
  measurement_edges <- std_sol |> filter(op=="=~")
  regression_edges  <- std_sol |> filter(op=="~")
  covariance_edges  <- std_sol |> filter(op=="~~", lhs!=rhs)
  
  observed_vars <- setdiff(unique(measurement_edges$rhs), latent)
  
  #--------------------------------------------------
  # Calcula erros: e = 1 - est.std^2
  #--------------------------------------------------
  errors <- sapply(observed_vars, function(var) {
    meas_row <- measurement_edges[measurement_edges$rhs == var, ]
    1 - (meas_row$est.std)^2
  })
  names(errors) <- observed_vars
  error_nodes <- paste0("e_", observed_vars)
  
  #--------------------------------------------------
  # Disturbances = fatores latentes que são indicadores de outro latente
  #--------------------------------------------------
  latent_as_indicators <- unique(as.character(measurement_edges[measurement_edges$rhs %in% latent, "rhs"]))
  disturbance_factors <- latent_as_indicators
  
  disturbances <- as.numeric(sapply(disturbance_factors, function(factor){
    rows <- measurement_edges[measurement_edges$rhs == factor, ]
    var_explained <- sum(as.numeric(rows$est.std)^2)
    1 - var_explained
  }))
  names(disturbances) <- disturbance_factors
  disturbance_nodes <- paste0("d_", disturbance_factors)
  
  #--------------------------------------------------
  # Código DOT
  #--------------------------------------------------
  
  dot_code <- "
digraph SEM {
  rankdir=LR
  bgcolor=transparent
  nodesep=0.25
  ranksep=0.6
  splines=false

  node [fontname=Helvetica]
  edge [fontname=Helvetica]

  # Nós latentes
  node [shape=circle style=filled color=lightblue width=1 fixedsize=true fontname=Helvetica]
"
  
  # Distribuição por ordem dos fatores
  if(length(third_order) > 0){
    dot_code <- paste0(dot_code, "{rank=same; ", paste(third_order, collapse='; '), ";}\n")
  }
  if(length(second_order) > 0){
    dot_code <- paste0(dot_code, "{rank=same; ", paste(second_order, collapse='; '), ";}\n")
  }
  if(length(first_order) > 0){
    dot_code <- paste0(dot_code, "{rank=same; ", paste(first_order, collapse='; '), ";}\n")
  }
  
  # Nós observáveis (quadrados)
  dot_code <- paste0(dot_code, "
node [shape=box style=filled color=lightgray width=1 fixedsize=true fontname=Helvetica]
", paste(observed_vars, collapse='; '), ";
")
  
  # Nós de erro (círculos pequenos cor-de-rosa)
  error_labels <- paste0(error_nodes,
                         " [shape=circle style=filled color=pink width=0.5 fixedsize=true fontname=Helvetica label=\"",
                         round(errors,2), "\"]")
  dot_code <- paste0(dot_code, paste(error_labels, collapse=';\n'), ";\n")
  
  # Nós de disturbance (só se existirem)
  if(length(disturbance_factors) > 0){
    for(i in seq_along(disturbance_factors)){
      dot_code <- paste0(dot_code,
                         disturbance_nodes[i], 
                         " [shape=circle style=filled color=orange width=0.6 fixedsize=true fontname=Helvetica label=\"",
                         round(disturbances[i], 2), "\"];\n")
    }
    
    # Subgraph para forçar cada disturbance no mesmo rank que o seu fator
    for(i in seq_along(disturbance_factors)){
      dot_code <- paste0(dot_code,
                         "subgraph { rank=same; ", disturbance_nodes[i], "; ", disturbance_factors[i], "; }\n")
    }
    
    # Força nós de disturbance no mesmo rank que os fatores endógenos
    dot_code <- paste0(dot_code,
                       "{rank=same; ", paste(c(disturbance_nodes, disturbance_factors), collapse="; "), ";}\n")
  }
  
  # Força círculos de erro à direita dos quadrados (invisível)
  for(i in seq_along(observed_vars)){
    dot_code <- paste0(dot_code,
                       observed_vars[i], " -> ", error_nodes[i], " [style=invis];\n")
  }
  
  # Arestas de mensuração (fatores -> observáveis)
  meas_code <- apply(measurement_edges, 1, function(x){
    sprintf("%s:e -> %s:w [label=<<TABLE BORDER='0' CELLBORDER='0' CELLPADDING='1' BGCOLOR='white'><TR><TD>%.2f</TD></TR></TABLE>> arrowsize=.6 labelfloat=true fontname=Helvetica]",
            x["lhs"], x["rhs"], as.numeric(x["est.std"]))
  })
  # Arestas de regressão
  reg_code <- apply(regression_edges, 1, function(x){
    sprintf("%s -> %s [label=<<TABLE BORDER='0' CELLBORDER='0' CELLPADDING='1' BGCOLOR='white'><TR><TD>%.2f</TD></TR></TABLE>> arrowsize=.6 fontname=Helvetica]",
            x["lhs"], x["rhs"], as.numeric(x["est.std"]))
  })
  
  # Arestas de covariância
  cov_code <- apply(covariance_edges, 1, function(x){
    sprintf("%s -> %s [dir=both arrowhead=normal arrowtail=normal label=<<TABLE BORDER='0' CELLBORDER='0' CELLPADDING='2' BGCOLOR='white'><TR><TD><FONT FACE='Helvetica'>%.2f</FONT></TD></TR></TABLE>> arrowsize=0.6 fontname=Helvetica fontcolor=black constraint=false]",
            x["lhs"], x["rhs"], as.numeric(x["est.std"]))
  })
  
  # Arestas de erro (dos nós de erro para os observáveis)
  error_edges <- sapply(seq_along(observed_vars), function(i){
    sprintf("%s -> %s [arrowsize=.6 tailport=w headport=e constraint=false]", 
            error_nodes[i], observed_vars[i])
  })
  
  # Arestas de disturbance (só se existirem)
  disturbance_edges <- if(length(disturbance_factors) > 0){
    sapply(seq_along(disturbance_factors), function(i){
      sprintf("%s -> %s [arrowsize=.6 headport=n]",
              disturbance_nodes[i], disturbance_factors[i])
    })
  } else { character(0) }
  
  # Alinha nós de erro à direita dos observáveis
  for(i in seq_along(observed_vars)){
    dot_code <- paste0(dot_code,
                       observed_vars[i], " -> ", error_nodes[i], 
                       " [style=invis weight=10];\n")
  }
  
  # Junta tudo
  dot_code <- paste(dot_code,
                    paste(meas_code, collapse="\n"),
                    paste(reg_code, collapse="\n"),
                    paste(cov_code, collapse="\n"),
                    paste(error_edges, collapse="\n"),
                    if(length(disturbance_edges) > 0) paste(disturbance_edges, collapse="\n"),
                    "}", sep="\n")
  
  # Plot
  #cat(dot_code)  # para debug
  plotMod <- DiagrammeR::grViz(dot_code)
  
  #--------------------------------------------------
  # FIT INDICES
  #--------------------------------------------------
  GOF <- c("df","chisq","pvalue","cfi","tli","nfi","srmr","rmsea","rmsea.ci.lower","rmsea.ci.upper","rmsea.pvalue")
  listGOF <- round(fitmeasures(fit, GOF),3)
  dfGOF <- data.frame(
    "χ²" = listGOF[2],          
    "df" = listGOF[1],
    "CFI" = listGOF[4],
    "TLI" = listGOF[5],
    "NFI" = listGOF[6],
    "SRMR" = listGOF[7],
    "RMSEA" = listGOF[8],
    "IC90" = paste0("[", listGOF[9], ";", listGOF[10], "]"),
    "P(RMSEA≤0.05)" = listGOF[11],
    check.names = FALSE  # prevent R from renaming columns
  )
  
  # Determine overall fit
  poor_indices <- c()
  if(dfGOF$CFI <= 0.9) poor_indices <- c(poor_indices, "CFI")
  if(dfGOF$TLI <= 0.9) poor_indices <- c(poor_indices, "TLI")
  if(dfGOF$SRMR >= 0.08) poor_indices <- c(poor_indices, "SRMR")
  if(dfGOF$RMSEA >= 0.06) poor_indices <- c(poor_indices, "RMSEA")
  
  if(length(poor_indices) == 0){
    footnote_text <- "Overall fit: Good (CFI & TLI ≥ 0.9, SRMR ≤ 0.08, RMSEA ≤ 0.06)."
  } else {
    footnote_text <- paste0(
      "Overall fit: Poor fit. Indices indicating poor fit: ",
      paste(poor_indices, collapse=", "), "."
    )
  }
  
  # Create table with footnote
  TabGOF <- knitr::kable(
    dfGOF,
    caption = "Fit indices for measurement model",
    row.names = FALSE,
    format = "html"
  ) |> 
    kableExtra::kable_classic() |>
    kableExtra::kable_styling(full_width = TRUE) |>
    kableExtra::footnote(general = footnote_text)
  
  
  #--------------------------------------------------
  # FIABILIDADE
  #--------------------------------------------------
  fiab <- semTools::reliability(fit)
  latent_vars <- unique(ptable$lhs[ptable$op=="=~"])
  indicators <- ptable[ptable$op=="=~",]
  second_order_factors <- unique(indicators$lhs[indicators$rhs %in% latent_vars])
  
  summary_table <- as.data.frame(fiab)
  
  if(length(second_order_factors)>0){
    Omega2ndOrd <- semTools::reliabilityL2(fit, second_order_factors)
    second_table <- as.data.frame(Omega2ndOrd)
    final_table <- dplyr::bind_rows(summary_table,second_table)
  } else {
    final_table <- summary_table
  }
  
  Tabfiab <- knitr::kable(final_table, caption="Reliability for model factors", digits=3) |> kableExtra::kable_classic() |> kableExtra::kable_styling(full_width=TRUE)
  
  #--------------------------------------------------
  # INICIALIZA RESULTADOS
  #--------------------------------------------------
  results <- list(
    fit = fit,
    Sensibilidade = descT,
    PesosFatoriais = TabLoad,
    IndicesGOF = TabGOF,
    PathPlot = plotMod,
    Fiabilidade = Tabfiab
  )
  
  #--------------------------------------------------
  # VALIDADES DISCRIMINANTE E CONVERGENTE
  #--------------------------------------------------
  # ---- 9. VALIDADE DISCRIMINANTE E CONVERGENTE (opcional) ------------------------
  if (vd) {
    
    # AVE e matriz de correlações latentes
    ave_values <- fiab["avevar", ]
    cor_matrix <- lavInspect(fit, "cor.lv")
    
    first_order_factors <- setdiff(latent_vars, second_order)
    cor_matrix <- cor_matrix[first_order_factors, first_order_factors, drop = FALSE]
    
    # Matriz Fornell-Larcker
    sqrt_ave <- sqrt(ave_values[first_order_factors])
    fl_matrix <- cor_matrix
    diag(fl_matrix) <- sqrt_ave
    
    # Formatar: manter apenas triângulo inferior
    fl_display <- round(fl_matrix, 3)
    fl_display[upper.tri(fl_display)] <- ""
    
    # Verificar validade convergente (AVE ≥ 0.5)
    conv_ok <- ave_values[first_order_factors] >= 0.5
    conv_msg <- if (all(conv_ok)) {
      "All constructs meet convergent validity (AVE ≥ 0.5)"
    } else {
      paste("Convergent validity issues in…:", 
            paste(names(conv_ok)[!conv_ok], collapse = ", "))
    }
    
    tab_convergent <- data.frame(
      Construct = names(ave_values[first_order_factors]),
      AVE = round(ave_values[first_order_factors], 3),
      `Validade Convergente` = ifelse(conv_ok, "Sim", "Não"),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    
    TabCV <- kable(tab_convergent, caption = "Convergent Validity (AVE ≥ 0.5)") |>
      kable_classic() |>
      footnote(general = enc2utf8(conv_msg))
    
    # Verificar validade discriminante (√AVE > correlações fora da diagonal)
    disc_issues <- character()
    for (i in seq_len(nrow(fl_matrix))) {
      for (j in seq_len(ncol(fl_matrix))) {
        if (i != j && diag(fl_matrix)[i] < abs(fl_matrix[i, j])) {
          disc_issues <- c(disc_issues, 
                           paste0("Problem between ", rownames(fl_matrix)[i], " and ", colnames(fl_matrix)[j]))
        }
      }
    }
    disc_msg <- if (length(disc_issues) == 0) {
      "All constructs meet discriminant validity."
    } else {
      paste(" ", paste(disc_issues, collapse = "; "))
    }
    
    # Nota de rodapé unificada
    footnote_txt <- paste(
      "• Discriminant Validity: √AVE should be greater than the off-diagonal correlations.",
      "• Convergent Validity: AVE ≥ 0.5 (or √AVE ≥ 0.707).",
      paste("• Discriminant Validity?", ifelse(length(disc_issues)==0, "Yes", "No")),
      disc_msg,
      paste("• Convergent Validity?", ifelse(all(conv_ok), "Yes", "No")),
      conv_msg,
      sep = "\n"
    )
    
    TabVD <- kable(fl_display, digits = 3, format = "html",
                   caption = "Matriz Fornell-Larcker (Diagonal = $\\sqrt{AVE}$, Off-diagonal = correlations)") |>
      kable_classic() |>
      footnote(general = enc2utf8(footnote_txt))
    
    # Adiciona ao results
    results$ValidadeDiscriminante <- TabVD
    results$ValidadeConvergente <- TabCV
  }
  
  return(results)
}