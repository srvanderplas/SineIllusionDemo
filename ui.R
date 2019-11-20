library(shiny)

# Define UI for application that plots random distributions
shinyUI(pageWithSidebar(

  # Application title
  headerPanel("Sine Illusion"),

  sidebarPanel(
    sliderInput("obs",
                "Number of lines:",
                min = 20,
                max = 100,
                value = 40, step=5),
    sliderInput("ell", "Line Length:", min=.5, max=5, value=1, step=.1),
    sliderInput("amp", "Amplitude:", min=.5, max=5, value=1, step=.1),
    conditionalPanel(condition="input.whichtab=='Sine Illusion'",
                     checkboxInput("hidelines", "Hide lines", FALSE)
    ),
    conditionalPanel(condition="input.whichtab=='Correct Y Axis'",
                     radioButtons("correct", "Correction Type:", c("No Correction" = "none",
                                                                   "Trigonometry Correction" = "geom",
                                                                   # "Linear Correction" = "linear",
                                                                   "Quadratic Correction" = "quad"))
    ),
    conditionalPanel(condition="input.whichtab!= 1",
                     sliderInput("weight", "Shrink by 1-w", min=0, max=1, value=1, step=.05)
    )
  ),
  mainPanel(
    tabsetPanel(id="whichtab",
                tabPanel("Sine Illusion", plotOutput("illusion", width="auto")),
                tabPanel("Correct X Axis", plotOutput("xcorrect", width="auto"), height="600px"),
                tabPanel("Correct Y Axis", plotOutput("ycorrect", width="auto", height="600px"))
    )
  )
))
