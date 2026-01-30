# check for tasks to do
findR::findRscript(pattern = 'TODO', path = '.', case.sensitive = FALSE,
                   comments = TRUE)

# check for place markers
findR::findRscript(pattern = 'HERE', path = '.', case.sensitive = TRUE,
                   comments = TRUE)

findR::findRscript(pattern = 'DELETE', path = '.', case.sensitive = TRUE,
                   comments = TRUE)

