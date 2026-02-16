
### function to round the results to show with the second decimal place
scaleFUN <- function(x) {
  a <- sprintf("%.1f", x)
  a2 <- sub('^-', '\U2212', format(a))
  a2 <- trimws(a2)
  return(a2)
}

scaleFUN2 <- function(x) {
  a <- sprintf("%.2f", x)
  a2 <- sub('^-', '\U2212', format(a))
  a2 <- trimws(a2)
  return(a2)
}

scaleFUN_int <- function(x) {
  a <- sprintf("%.0f", x)
  a2 <- sub('^-', '\U2212', format(a))
  a2 <- trimws(a2)
  return(a2)
}

### function to convert p-values to asterisks　
p_to_sig <- function(p, type = "asterisk") {
  if (type == "asterisk") {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "n.s.")
         )
  )
  } else if (type == "threshold") {
    ifelse(p < 0.001, 'italic(p)*" < 0.001"',
           ifelse(p < 0.01, 'italic(p)*" < 0.01"',
                  ifelse(p < 0.05, 'italic(p)*" < 0.05"', '"n.s."')
                  )
           )
  } else { 
    stop('Please specify type = "asterisk" or "threshold".')
    }
}

