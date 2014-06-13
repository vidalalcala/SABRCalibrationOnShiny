library(shiny)
# Define UI for application that plots random distributions 
shinyUI(pageWithSidebar(
  headerPanel("SABR calibration on Shiny"),
  #a('https://s3.amazonaws.com/github/ribbons/forkme_right_darkblue_121621.png', href='https://github.com/timelyportfolio/rCharts_nvd3_perf'),
  sidebarPanel(
    numericInput("forward", "Forward:", 93.5, min=0.0, max=1000, step=0.01),
    numericInput("maturity", "Maturity(year):", 1/12, min=0.0, max=4, step=0.0001),
    numericInput("r", "Risk free continuous rate:", 0.0015, min=0.0, max=1.00, step=0.0001),
    numericInput("minVol", "Minimum volume:", 400, min=0, max=1000000, step=1),
    numericInput("maxVol", "Maximum volume:", 10000, min=0, max=10000000, step=1),
    numericInput("nu0", "Initial nu:", 2.7 , min=0.0, max=1000.0, step=0.001),
    numericInput("alpha0", "Initial alpha:", 0.22 , min=0.0, max=100.0, step=0.001),
    numericInput("rho0", "Initial rho:", -0.12 , min=-1.0, max=0.0, step=0.001),
    numericInput("beta0", "Initial beta:", 1.0 , min=0.0, max=2.0, step=0.001),
    numericInput("minSigma", "Minimum volatility:", 0.001 , min=0.0, max=100.0, step=0.001),
    numericInput("maxSigma", "Maximum volatility:", 1.0 , min=0.0, max=100.0, step=0.0001)
  ),
  mainPanel(
    tabsetPanel(
      tabPanel("Hagan", 
        h3("Implied volatility (market and SABR model)"),
        plotOutput("distPlotHagan"),
        h3("Calibrated parameters via market prices"),
        verbatimTextOutput("summaryHagan"),
        h3("Delta hedge (SABR model)"),
        plotOutput("deltaPlotHagan"),
        tableOutput("deltaTableHagan")
      ),
      tabPanel("Nu expansion", 
               h3("Implied volatility (market and nuSABR model)"),
               plotOutput("distPlotSABR"),
               h3("Option prices (market and nuSABR model)"),
               plotOutput("pricesPlotSABR"),
               h3("Calibrated parameters via market prices"),
               verbatimTextOutput("summarySABR"),
               h3("Delta hedge (nuSABR model)"),
               plotOutput("deltaPlotSABR"),
               tableOutput("deltaTableSABR")
      ),
      tabPanel("About",
        p('This application demonstrates to what extent',
         a("SABRnu model", href="http://en.wikipedia.org/wiki/SABR_volatility_model", target="_blank"),
         'can fit the market prices.'),
        br(),
        strong('Code'),
        p('Source code for this application at',
          a('GitHub', href='https://github.com/vidalalcala/SABRCalibrationOnShiny', target="_blank")),
        p('If you want to run this code on your computer, run the code below:',
          br(),
          code('library(shiny)'),br(),
          code('runGitHub("SABRCalibrationOnShiny","vidalalcala")')
        ),br(),
        strong('References'),
        p(HTML('<ul>'),
        HTML('<li>'),a('Managing Smile Risk, P. Hagan et al(pdf)', href="http://www.math.columbia.edu/~lrb/sabrAll.pdf", target="_blank"),HTML('</li>'),
        HTML('</ul>'))
      )
    )
  )
))
