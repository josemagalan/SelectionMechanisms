# app.R  ---  Selection mechanisms demo (Shiny)
# UI: English
# Methods included:
# 1) Tournament (k + checkbox "without replacement within tournament")
# 2) Roulette wheel (scaling: none / linear)
# 3) SUS (Stochastic Universal Sampling) (scaling: none / linear)
# 4) Rank selection (Baker linear rank; parameter s in [1,2])
# 5) Truncation (top_p)
# 6) Boltzmann selection (temperature T)

library(shiny)
library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)
library(rlang)

`%||%` <- rlang::`%||%`

# ============================
# Helpers
# ============================

shannon_index <- function(x) {
  p <- table(x) / length(x)
  -sum(p * log(p))
}

geom_line_compat <- function(width = 1.1, ...) {
  if ("linewidth" %in% names(formals(ggplot2::geom_line))) {
    ggplot2::geom_line(linewidth = width, ...)
  } else {
    ggplot2::geom_line(size = width, ...)
  }
}

summarise_bands <- function(df, value_col) {
  value_col <- rlang::ensym(value_col)
  df %>%
    group_by(gen) %>%
    summarise(
      mean = mean(!!value_col, na.rm = TRUE),
      p05  = quantile(!!value_col, 0.05, na.rm = TRUE),
      p95  = quantile(!!value_col, 0.95, na.rm = TRUE),
      .groups = "drop"
    )
}

apply_linear_scaling <- function(w, scaling = "none", a = 1, b = 0) {
  if (scaling == "linear") return(a * w + b)
  w
}

safe_probs <- function(w, eps = 1e-9) {
  if (any(!is.finite(w))) stop("Non-finite weights (Inf/NaN).")
  w <- pmax(w, eps)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) stop("Invalid sum of weights.")
  w / sw
}

sus_sample <- function(ids, prob, n_select) {
  if (length(ids) != length(prob)) stop("SUS: ids/prob length mismatch.")
  if (n_select <= 0) return(integer(0))
  
  prob <- prob / sum(prob)
  cum <- cumsum(prob)
  
  step <- 1 / n_select
  u0 <- runif(1, 0, step)
  ptrs <- u0 + step * (0:(n_select - 1))
  
  idx <- vapply(ptrs, function(u) which(cum >= u)[1], integer(1))
  ids[idx]
}

# ============================
# Selection core (without elitism)
# Returns vector of length n_select
# ============================

select_core <- function(pop_ids, fit, method, params, n_select) {
  
  if (n_select <= 0) return(integer(0))
  
  if (method == "tournament") {
    k <- params$k
    no_replace <- isTRUE(params$no_replace)
    
    replicate(n_select, {
      cand <- sample(pop_ids, size = k, replace = !no_replace)
      cand[which.max(fit[cand])]
    })
  }
  
  else if (method == "roulette") {
    w <- fit[pop_ids]
    scaling <- params$scaling %||% "none"
    w <- apply_linear_scaling(w, scaling, a = params$a %||% 1, b = params$b %||% 0)
    prob <- safe_probs(w, eps = 1e-9)
    sample(pop_ids, size = n_select, replace = TRUE, prob = prob)
  }
  
  else if (method == "sus") {
    w <- fit[pop_ids]
    scaling <- params$scaling %||% "none"
    w <- apply_linear_scaling(w, scaling, a = params$a %||% 1, b = params$b %||% 0)
    prob <- safe_probs(w, eps = 1e-9)
    sus_sample(pop_ids, prob, n_select = n_select)
  }
  
  else if (method == "rank") {
    # Baker linear rank selection:
    # s in [1,2], expected copies of best = s, worst = 2 - s
    s <- params$s
    N <- length(pop_ids)
    if (N < 2) return(rep(pop_ids, length.out = n_select))
    if (!is.finite(s) || s < 1 || s > 2) stop("Rank parameter s must be in [1,2].")
    
    # ranks: 1 = worst, N = best
    ord <- order(fit[pop_ids], decreasing = FALSE)
    r <- seq_len(N)
    
    p_r <- (2 - s) / N + (2 * (r - 1) * (s - 1)) / (N * (N - 1))
    
    prob <- numeric(N)
    prob[ord] <- p_r
    prob <- safe_probs(prob, eps = 1e-12)
    
    sample(pop_ids, size = n_select, replace = TRUE, prob = prob)
  }
  
  else if (method == "truncation") {
    top_p <- params$top_p
    N <- length(pop_ids)
    n_keep <- max(1L, ceiling(top_p * N))
    keep <- pop_ids[order(fit[pop_ids], decreasing = TRUE)][1:n_keep]
    sample(keep, size = n_select, replace = TRUE)
  }
  
  else if (method == "boltzmann") {
    T <- params$T
    if (!is.finite(T) || T <= 0) stop("Temperature T must be > 0.")
    f <- fit[pop_ids]
    f0 <- max(f)
    w <- exp((f - f0) / T)  # stable
    prob <- safe_probs(w, eps = 1e-12)
    sample(pop_ids, size = n_select, replace = TRUE, prob = prob)
  }
  
  else {
    stop("Method not implemented: ", method)
  }
}

# ============================
# Elitism wrapper
# ============================

select_next_population <- function(pop_ids, fit, method, params, elite_n) {
  N <- length(pop_ids)
  elite_n <- max(0L, min(as.integer(elite_n), N))
  
  elites <- integer(0)
  if (elite_n > 0) {
    elites <- pop_ids[order(fit[pop_ids], decreasing = TRUE)][1:elite_n]
  }
  
  rest <- select_core(pop_ids, fit, method, params, n_select = N - elite_n)
  c(elites, rest)
}

# ============================
# Simulation
# ============================

simulate_one <- function(N, G, method, params, elite_n,
                         fitness_mode = c("permutation", "random"),
                         seed = NULL) {
  fitness_mode <- match.arg(fitness_mode)
  if (!is.null(seed)) set.seed(seed)
  
  fit <- switch(
    fitness_mode,
    permutation = sample(1:N, N, replace = FALSE),
    random      = sample(1:100, N, replace = TRUE)
  )
  
  pop <- 1:N  # sequential init
  
  out <- vector("list", G + 1)
  for (g in 0:G) {
    out[[g + 1]] <- tibble(
      gen = g,
      unique_ids   = dplyr::n_distinct(pop),
      shannon_ids  = shannon_index(pop),
      mean_fitness = mean(fit[pop]),
      max_fitness  = max(fit[pop])
    )
    if (g < G) pop <- select_next_population(pop, fit, method, params, elite_n)
  }
  bind_rows(out)
}

simulate_many <- function(n_sims, N, G, method, params, elite_n, fitness_mode, seed = NULL) {
  seeds <- if (is.null(seed)) rep(NA_integer_, n_sims) else seed + seq_len(n_sims)
  map_dfr(seq_len(n_sims), function(i) {
    simulate_one(
      N = N, G = G,
      method = method,
      params = params,
      elite_n = elite_n,
      fitness_mode = fitness_mode,
      seed = if (is.na(seeds[i])) NULL else seeds[i]
    ) %>% mutate(sim = i)
  })
}

# ============================
# UI
# ============================

ui <- fluidPage(
  titlePanel("Selection mechanisms: convergence and loss of genetic diversity"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Common settings"),
      numericInput("N", "Population size (N)", value = 100, min = 20, step = 10),
      sliderInput("G", "Generations (G)", min = 10, max = 200, value = 60, step = 10),
      sliderInput("nsims", "Simulations", min = 10, max = 500, value = 30, step = 10),
      
      selectInput(
        "fitness_mode", "Fitness distribution per ID",
        choices = c(
          "Permutation 1..N (no repetition)" = "permutation",
          "Random 1..100 (with repetition)" = "random"
        ),
        selected = "permutation"
      ),
      numericInput("seed", "Seed (optional)", value = NA, min = 1),
      
      hr(),
      h4("Selection method"),
      selectInput(
        "method", "Method",
        choices = c(
          "Tournament" = "tournament",
          "Roulette wheel" = "roulette",
          "SUS (Stochastic Universal Sampling)" = "sus",
          "Rank selection (Baker linear rank)" = "rank",
          "Truncation selection" = "truncation",
          "Boltzmann selection" = "boltzmann"
        ),
        selected = "tournament"
      ),
      
      uiOutput("method_params"),
      
      hr(),
      actionButton("run", "▶ Run", class = "btn-primary", width = "100%"),
      fluidRow(
        column(6, actionButton("add_comp", "➕ Add to comparison", width = "100%")),
        column(6, actionButton("clear_comp", "🧹 Clear comparison", width = "100%"))
      )
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Diversity",
                 plotOutput("plot_unique", height = 320),
                 plotOutput("plot_shannon", height = 320)),
        tabPanel("Fitness",
                 plotOutput("plot_fit_mean", height = 320),
                 plotOutput("plot_fit_max", height = 320)),
        tabPanel("Summary",
                 tableOutput("final_table")),
        tabPanel("Comparison",
                 helpText("Overlays saved experiments (Add to comparison)."),
                 plotOutput("plot_comp_unique", height = 320),
                 plotOutput("plot_comp_shannon", height = 320),
                 plotOutput("plot_comp_fit", height = 320),
                 tableOutput("comp_table"))
      )
    )
  )
)

# ============================
# Server
# ============================

server <- function(input, output, session) {
  
  output$method_params <- renderUI({
    m <- input$method
    
    if (m == "tournament") {
      tagList(
        sliderInput("k", "Tournament size (k)", min = 2, max = 10, value = 2),
        checkboxInput("tournament_nr", "Without replacement within tournament", value = FALSE),
        sliderInput("elite", "Elitism (n elites)", min = 0, max = 20, value = 1)
      )
      
    } else if (m %in% c("roulette", "sus")) {
      tagList(
        sliderInput("elite", "Elitism (n elites)", min = 0, max = 20, value = 1),
        selectInput(paste0(m, "_scaling"), "Scaling",
                    choices = c("None" = "none", "Linear (f' = a f + b)" = "linear"),
                    selected = "none"),
        conditionalPanel(
          condition = sprintf("input.%s_scaling == 'linear'", m),
          sliderInput("a", "a (slope)", min = 0.1, max = 5, value = 1, step = 0.1),
          sliderInput("b", "b (intercept)", min = -200, max = 200, value = 0, step = 1),
          helpText("Weights are w = max(a·fitness + b, ε).")
        ),
        helpText(if (m == "roulette")
          "Roulette: N independent draws with probability proportional to weight."
          else
            "SUS: equidistant pointers (lower sampling variance than roulette).")
      )
      
    } else if (m == "rank") {
      tagList(
        sliderInput("elite", "Elitism (n elites)", min = 0, max = 20, value = 1),
        sliderInput("rank_s", "Selection pressure s (Baker) [1..2]", min = 1, max = 2, value = 1.5, step = 0.05),
        helpText("s=1 ~ near-uniform by rank; s=2 ~ highest pressure by rank.")
      )
      
    } else if (m == "truncation") {
      tagList(
        sliderInput("top_p", "Top proportion (p)", min = 0.05, max = 0.95, value = 0.5, step = 0.05),
        sliderInput("elite", "Elitism (n elites)", min = 0, max = 20, value = 1),
        helpText("Only the best p·N compete; then sampling with replacement.")
      )
      
    } else if (m == "boltzmann") {
      tagList(
        sliderInput("elite", "Elitism (n elites)", min = 0, max = 20, value = 1),
        sliderInput("T", "Temperature T", min = 0.1, max = 50, value = 5, step = 0.1),
        helpText("Weights: exp((f - max(f))/T). High T ~ soft; low T ~ aggressive.")
      )
      
    } else {
      helpText("Configure a method.")
    }
  })
  
  get_engine_method_and_params <- reactive({
    m <- input$method
    
    if (m == "tournament") {
      req(input$k)
      list(
        method = "tournament",
        params = list(
          k = as.integer(input$k),
          no_replace = isTRUE(input$tournament_nr)
        )
      )
      
    } else if (m == "roulette") {
      sc <- input$roulette_scaling %||% "none"
      if (sc == "linear") {
        req(input$a, input$b)
        list(method = "roulette", params = list(scaling = "linear", a = as.numeric(input$a), b = as.numeric(input$b)))
      } else {
        list(method = "roulette", params = list(scaling = "none"))
      }
      
    } else if (m == "sus") {
      sc <- input$sus_scaling %||% "none"
      if (sc == "linear") {
        req(input$a, input$b)
        list(method = "sus", params = list(scaling = "linear", a = as.numeric(input$a), b = as.numeric(input$b)))
      } else {
        list(method = "sus", params = list(scaling = "none"))
      }
      
    } else if (m == "rank") {
      req(input$rank_s)
      list(method = "rank", params = list(s = as.numeric(input$rank_s)))
      
    } else if (m == "truncation") {
      req(input$top_p)
      list(method = "truncation", params = list(top_p = as.numeric(input$top_p)))
      
    } else if (m == "boltzmann") {
      req(input$T)
      list(method = "boltzmann", params = list(T = as.numeric(input$T)))
      
    } else {
      list(method = "tournament", params = list(k = 2L, no_replace = FALSE))
    }
  })
  
  results <- eventReactive(input$run, {
    N <- as.integer(input$N)
    G <- as.integer(input$G)
    ns <- as.integer(input$nsims)
    elite_n <- if (is.null(input$elite) || is.na(input$elite)) 0L else as.integer(input$elite)
    
    validate(
      need(!is.na(N) && N >= 2, "Invalid N."),
      need(!is.na(G) && G >= 1, "Invalid G."),
      need(!is.na(ns) && ns >= 1, "Invalid number of simulations."),
      need(elite_n >= 0, "Invalid elitism."),
      need(elite_n <= N, "Elitism cannot exceed N.")
    )
    
    mp <- get_engine_method_and_params()
    method <- mp$method
    params <- mp$params
    
    if (method == "tournament") {
      validate(need(params$k >= 2, "k must be >= 2."))
      validate(need(params$k <= N, "k cannot exceed N."))
    }
    if (method == "truncation") {
      validate(need(params$top_p > 0 && params$top_p < 1, "Top p must be in (0,1)."))
    }
    if (method == "boltzmann") {
      validate(need(is.finite(params$T) && params$T > 0, "Temperature must be > 0."))
    }
    if (method %in% c("roulette", "sus") && (params$scaling %||% "none") == "linear") {
      validate(
        need(is.finite(params$a) && params$a > 0, "a must be > 0."),
        need(is.finite(params$b), "Invalid b.")
      )
    }
    if (method == "rank") {
      validate(need(params$s >= 1 && params$s <= 2, "s must be in [1,2]."))
    }
    
    seed_val <- if (is.na(input$seed)) NULL else as.integer(input$seed)
    
    simulate_many(
      n_sims = ns,
      N = N, G = G,
      method = method,
      params = params,
      elite_n = elite_n,
      fitness_mode = input$fitness_mode,
      seed = seed_val
    )
  })
  
  summary_gen <- reactive({
    req(results())
    list(
      unique   = summarise_bands(results(), unique_ids),
      shannon  = summarise_bands(results(), shannon_ids),
      fit_mean = summarise_bands(results(), mean_fitness),
      fit_max  = summarise_bands(results(), max_fitness)
    )
  })
  
  output$plot_unique <- renderPlot({
    df <- summary_gen()$unique
    ggplot(df, aes(gen, mean)) +
      geom_line_compat(1.1) +
      geom_ribbon(aes(ymin = p05, ymax = p95), alpha = 0.25) +
      labs(title = "Diversity: unique IDs (mean and p05–p95)", x = "Generation", y = "Unique IDs") +
      theme_minimal()
  })
  
  output$plot_shannon <- renderPlot({
    df <- summary_gen()$shannon
    ggplot(df, aes(gen, mean)) +
      geom_line_compat(1.1) +
      geom_ribbon(aes(ymin = p05, ymax = p95), alpha = 0.25) +
      labs(title = "Diversity: Shannon index on IDs (mean and p05–p95)", x = "Generation", y = "Shannon") +
      theme_minimal()
  })
  
  output$plot_fit_mean <- renderPlot({
    df <- summary_gen()$fit_mean
    ggplot(df, aes(gen, mean)) +
      geom_line_compat(1.1) +
      geom_ribbon(aes(ymin = p05, ymax = p95), alpha = 0.25) +
      labs(title = "Convergence: mean fitness (mean and p05–p95)", x = "Generation", y = "Mean fitness") +
      theme_minimal()
  })
  
  output$plot_fit_max <- renderPlot({
    df <- summary_gen()$fit_max
    ggplot(df, aes(gen, mean)) +
      geom_line_compat(1.1) +
      geom_ribbon(aes(ymin = p05, ymax = p95), alpha = 0.25) +
      labs(title = "Max fitness (mean and p05–p95)", x = "Generation", y = "Max fitness") +
      theme_minimal()
  })
  
  output$final_table <- renderTable({
    req(results())
    last_gen <- max(results()$gen)
    
    results() %>%
      filter(gen == last_gen) %>%
      summarise(
        sims = n(),
        unique_ids_mean = mean(unique_ids),
        unique_ids_p05 = quantile(unique_ids, 0.05),
        unique_ids_p95 = quantile(unique_ids, 0.95),
        shannon_mean = mean(shannon_ids),
        mean_fitness_mean = mean(mean_fitness),
        max_fitness_mean = mean(max_fitness)
      )
  }, digits = 3)
  
  # ============================
  # Comparison store
  # ============================
  
  rv <- reactiveValues(comps = list(), meta = tibble())
  
  observeEvent(input$add_comp, {
    req(results())
    
    mp <- get_engine_method_and_params()
    method <- mp$method
    params <- mp$params
    
    method_name <- switch(
      input$method,
      tournament = "Tournament",
      roulette = "Roulette",
      sus = "SUS",
      rank = "Rank",
      truncation = "Truncation",
      boltzmann = "Boltzmann",
      input$method
    )
    
    param_str <- switch(
      input$method,
      tournament = paste0("k=", input$k, ifelse(isTRUE(input$tournament_nr), " (NR)", "")),
      roulette = {
        if ((params$scaling %||% "none") == "linear") paste0("sc=linear a=", round(params$a, 2), " b=", round(params$b, 2))
        else "sc=none"
      },
      sus = {
        if ((params$scaling %||% "none") == "linear") paste0("sc=linear a=", round(params$a, 2), " b=", round(params$b, 2))
        else "sc=none"
      },
      rank = paste0("s=", round(params$s, 2)),
      truncation = paste0("top_p=", input$top_p),
      boltzmann = paste0("T=", round(params$T, 2)),
      ""
    )
    
    label <- paste0(
      method_name, " | ",
      "N=", input$N,
      " ", param_str,
      " e=", input$elite,
      " G=", input$G,
      " ns=", input$nsims,
      " fit=", input$fitness_mode,
      ifelse(is.na(input$seed), "", paste0(" seed=", input$seed))
    )
    
    comp_df <- results() %>%
      group_by(gen) %>%
      summarise(
        unique_mean  = mean(unique_ids, na.rm = TRUE),
        shannon_mean = mean(shannon_ids, na.rm = TRUE),
        fit_mean     = mean(mean_fitness, na.rm = TRUE),
        fit_max_mean = mean(max_fitness, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(model = label, method = method_name)
    
    rv$comps <- append(rv$comps, list(comp_df))
    rv$meta  <- bind_rows(rv$meta, tibble(model = label, method = method_name, added_at = Sys.time()))
  })
  
  observeEvent(input$clear_comp, {
    rv$comps <- list()
    rv$meta  <- tibble()
  })
  
  comp_all <- reactive({
    if (length(rv$comps) == 0) return(NULL)
    bind_rows(rv$comps)
  })
  
  output$plot_comp_unique <- renderPlot({
    df <- comp_all(); req(df)
    ggplot(df, aes(gen, unique_mean, color = model, group = model)) +
      geom_line_compat(1.1) +
      labs(title = "Comparison: unique IDs (mean)", x = "Generation", y = "Unique IDs (mean)") +
      theme_minimal()
  })
  
  output$plot_comp_shannon <- renderPlot({
    df <- comp_all(); req(df)
    ggplot(df, aes(gen, shannon_mean, color = model, group = model)) +
      geom_line_compat(1.1) +
      labs(title = "Comparison: Shannon index (mean)", x = "Generation", y = "Shannon (mean)") +
      theme_minimal()
  })
  
  output$plot_comp_fit <- renderPlot({
    df <- comp_all(); req(df)
    ggplot(df, aes(gen, fit_mean, color = model, group = model)) +
      geom_line_compat(1.1) +
      labs(title = "Comparison: mean fitness (mean)", x = "Generation", y = "Mean fitness (mean)") +
      theme_minimal()
  })
  
  output$comp_table <- renderTable({
    req(rv$meta)
    rv$meta
  }, digits = 0)
}

shinyApp(ui, server)