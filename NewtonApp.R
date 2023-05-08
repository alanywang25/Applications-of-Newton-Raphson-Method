# Load required libraries
library(shiny)
library(shinydashboard)
library(ggplot2)
library(tidyverse)
library(Deriv)

# Define the UI
ui <- dashboardPage(
  dashboardHeader(title = "Combined App"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Newton's Method Calculator", tabName = "newton_method", 
               icon = icon("calculator")),
      menuItem("Damped Driven Oscillator", tabName = "damped_oscillator", 
               icon = icon("wave-square")),
      menuItem("IRR Calculator", tabName = "irr_calculator", 
               icon = icon("chart-line"))
    )
  ),
  dashboardBody(
    tabItems(
      # Newton's Method Calculator
      tabItem(
        tabName = "newton_method",
        fluidPage(
          titlePanel("Newton's Method Calculator"),
          sidebarLayout(
            sidebarPanel(
              helpText(
              "Enter the function f(x), its derivative f'(x), initial guess x,
               tolerance, and maximum number of iterations (Note: when entering
               a term w/ coefficients, include the multiplication symbol. For
               example, enter 2x as 2*x):"),
              textInput("func", "Function f(x)", "x^2"),
              numericInput("x", "Initial Guess", 0),
              numericInput("tol", "Tolerance", 1e-8),
              numericInput("n", "Maximum number of iterations", 100),
              numericInput("xlim_min", "x-axis min bound", -10),
              numericInput("xlim_max", "x-axis max bound", 10),
              h4("Optional Inputs: "),
              checkboxInput("options", label = "Show options", value = FALSE),
              conditionalPanel(
                "input.options == true",
                numericInput("sigfig", "# of Significant Figures", 3),
                textInput("color", "Color of Estimated Root Point", "black"),
                numericInput("size", "Size of Estimated Root Point", 2)
              ),
              actionButton("calculate", "Calculate Root")
            ),
            mainPanel(
              h3("Result:"),
              verbatimTextOutput("result"),
              plotOutput("plot"),
              plotOutput("plot2")
            )
          )
        )
      ),
      
      # Damped Driven Oscillator
      tabItem(
        tabName = "damped_oscillator",
        fluidPage(
          titlePanel("Damped Driven Harmonic Oscillator Simulation"),
          sidebarLayout(
            sidebarPanel(
              numericInput("mass", "Mass (m)", value = 1),
              numericInput("damping", "Damping Coefficient (c)", value = 0.1),
              numericInput("spring", "Spring Constant (k)", value = 10),
              numericInput("force_amplitude", "Force Amplitude (F0)", 
                           value = 1),
              numericInput("angular_frequency", "Angular Frequency (w)", 
                           value = 1),
              actionButton("simulate", "Simulate")
            ),
            mainPanel(
              plotOutput("oscillator_plot")
            )
          )
        )
      ),
      
      # IRR Calculator
      tabItem(
        tabName = "irr_calculator",
        fluidPage(
          titlePanel("Internal Rate of Return (IRR) Calculator"),
          sidebarLayout(
            sidebarPanel(
              numericInput("initial_investment", "Initial Investment:", 
                           value = 3000),
              numericInput("guess", "Initial IRR Guess:", value = 0.1, min = -1, 
                           max = 1, step = 0.01),
              helpText("Enter cash flows separated by commas."),
              textInput("cash_flows", "Cash Flows:", 
                        value = "-1000, 200, 400, 600, 800"),
              actionButton("calculate_irr", "Calculate")
            ),
            mainPanel(
              textOutput("irr_result"),
              plotOutput("npv_plot")
            )
          )
        )
      )
    )
  )
)

# Newton Method/Bisection method fallback for IRR
newton_method <- function(f, f_prime, x0, tol = 1e-8, max_iter = 10000) {
  x <- x0
  
  for (i in 1:max_iter) {
    if (abs(f_prime(x)) < tol) {
      break
    }
    
    x_new <- x - f(x) / f_prime(x)
    
    if (abs(x_new - x) < tol) {
      return(x_new)
    }
    
    x <- x_new
  }
  
  return(x)
}

# Newton Method for Newton's Method Root Finder
newton <- function(f, fp, x, tol, n) {
  x_values <- c()
  for (i in 1:n) {
    if (fp(x) == 0) {
      x_values <- c(x_values, x)
      break
    }
    x <- x - f(x) / fp(x)
    x_values <- c(x_values, x)
    if (abs(f(x)) < tol) {
      break
    }
  }
  data.frame(Iteration = 1:length(x_values), X = x_values)
}

# Define server
server <- function(input, output) {
  # Newton's Method Calculator
  f <- reactive({
    func <- input$func
    eval(parse(text = paste0("f <- function(x){ ", func,"}")))
    environment(f) <- .GlobalEnv
    f
  })
  
  fp <- reactive({
    func_prime <- deriv(parse(text = input$func), 'x')
    eval(parse(text = paste0("fp <- function(x){ ", sub(".*<- ", "", 
                                                        deparse(func_prime)[4]),
                             "}")))
    environment(fp) <- .GlobalEnv
    fp
  })
  
  newton_output <- reactiveVal()
  
  observeEvent(input$calculate, {
    if (!is.null(f()) && !is.null(fp())) {
      newton_output(newton(f(), fp(), input$x, input$tol, input$n))
      output$result <- renderText({
        paste0("The root found is:", signif(tail(newton_output()$X, 1), 
                                            input$sigfig))
      })
    }
  })
  
  output$plot <- renderPlot({
    if (!is.null(newton_output()) && nrow(newton_output()) > 0) {
      ggplot(newton_output(), aes(x = Iteration, y = X)) +
        geom_point() +
        geom_line() +
        labs(title = "X Movement with Each Iteration", x = "Iteration", 
             y = "X") +
        theme_minimal()
    }
  })
  output$plot2 <- renderPlot({
    if (!is.null(f()) && !is.null(newton_output()) 
        && nrow(newton_output()) > 0) {
      ggplot() +
        stat_function(fun = f()) +
        xlim(input$xlim_min, input$xlim_max) +
        geom_point(aes(x = tail(newton_output()$X, 1), y = 0), 
                   color = input$color, size = input$size) +
        labs(x = "x",y = "y", title = "Function f(x) with Estimated Root") +
        theme_minimal()
    }
  })


  # Damped Driven Oscillator
  observeEvent(input$simulate, {
    m <- input$mass
    c <- input$damping
    k <- input$spring
    F0 <- input$force_amplitude
    w <- input$angular_frequency
    
    newtons_method <- function(func, func_prime, x0, max_iter = 100, 
                               tol = 1e-6) {
      x_current <- x0
      for (i in 1:max_iter) {
        x_next <- x_current - func(x_current) / func_prime(x_current)
        if (abs(x_next - x_current) < tol) {
          break
        }
        x_current <- x_next
      }
      return(x_current)
    }
    
    dt <- 0.01
    num_steps <- 2000
    
    x <- 0
    v <- 0
    
    solution <- matrix(0, nrow = num_steps + 1, ncol = 2)
    solution[1, ] <- c(x, v)
    colnames(solution) <- c("x", "v")
    
    # Iterate through the time steps
    for (i in 2:(num_steps + 1)) {
      t <- (i - 1) * dt
      
      # Define the function and its derivative for displacement
      func_x <- function(x_next) {
        x_next - x - dt * v - 0.5 * dt^2 * (F0 * sin(w * t) - 
                                              c * v - k * x_next) / m
      }
      func_x_prime <- function(x_next) {
        1 - 0.5 * dt^2 * (-k) / m
      }
      
      # Update the displacement using Newton's method
      x <- newtons_method(func_x, func_x_prime, x)
      
      # Update the velocity using the new displacement
      v <- v + 0.5 * dt * ((F0 * sin(w * t) - c * v - k * x) / m + 
                             (F0 * sin(w * (t + dt)) - c * v - k * x) / m)
      
      # Store the new state
      solution[i, ] <- c(x, v)
    }
    
    # Create a data frame from the solution
    solution_df <- data.frame(time = seq(0, num_steps * dt, by = dt),
                              x = solution[, "x"])
    output$oscillator_plot <- renderPlot({
      ggplot(data.frame(time = (0:num_steps) * dt, x = solution[, "x"])) +
        geom_line(aes(x = time, y = x)) +
        labs(title = "Damped Driven Harmonic Oscillator Simulation",
             x = "Time",
             y = "Displacement") +
        theme_minimal()
    })
  })
  
  # IRR Calculator
  irr_calculation <- eventReactive(input$calculate_irr, {
    validate(
      need(input$guess >= -1 && input$guess <= 1, 
           "Please enter a valid IRR guess between -1 and 1.")
    )
    cash_flows <- c(-input$initial_investment, 
                    as.numeric(unlist(strsplit(input$cash_flows, ","))))
    n <- length(cash_flows)
    
    f <- function(r) {
      sum(cash_flows / ((1 + r) ^ (0:(n - 1))))
    }
    
    f_prime <- function(r) {
      sum(-((0:(n - 1)) * cash_flows) / ((1 + r) ^ (1:(n))))
    }
    
    irr <- tryCatch({
      newton_method(f, f_prime, input$guess)
    }, error = function(e) {
      NULL
    })
    
    irr
  })
  
  output$irr_result <- renderText({
    result <- irr_calculation()
    if (!is.null(result)) {
      paste("IRR:", round(result * 100, 2), "%")
    } else {
      "Error: Unable to calculate IRR. 
      Please try a different initial guess or check the cash flows."
    }
  })
  
  output$npv_plot <- renderPlot({
    cash_flows <- as.numeric(unlist(strsplit(input$cash_flows, ",")))
    validate(
      need(length(cash_flows) > 0, 
           "Please enter valid cash flows separated by commas.")
    )
    n <- length(cash_flows)
    
    f <- function(r) {
      sum(cash_flows / ((1 + r) ^ (0:(n - 1))))
    }
    
    r_range <- seq(-0.5, 0.5, length.out = 100)
    npv_values <- sapply(r_range, f)
    
    plot_data <- data.frame(r = r_range, NPV = npv_values)
    irr_value <- irr_calculation()
    
    if (!is.null(irr_value)) {
      plot_data$highlight <- ifelse(abs(plot_data$r - irr_value) < 1e-4, "IRR", 
                                    "Other")
      irr_label <- paste("IRR:", round(irr_value * 100, 2), "%")
      
      ggplot(plot_data, aes(x = r, y = NPV, color = highlight)) +
        geom_line(size = 1) +
        scale_color_manual(values = c("Other" = "gray","IRR" = "red")) +
        labs(title = "Net Present Value (NPV) vs. Discount Rate",
             x = "Discount Rate",
             y = "Net Present Value (NPV)") +
        theme_minimal() +
        geom_vline(xintercept = irr_value, linetype = "dashed", color = "red") +
        annotate("text", x = irr_value, y = min(plot_data$NPV), 
                 label = irr_label, color = "red", hjust = -0.5, vjust = -1, 
                 size = 4, fontface = "bold")
    } else {
      ggplot() +
        labs(title = "Error: Unable to plot NPV vs. Discount Rate",
             x = "Discount Rate",
             y = "Net Present Value (NPV)") +
        theme_minimal() +
        theme_void()
    }
  })
}

# Run app
shinyApp(ui, server)

                      
