# .onLoad = function(libname, pkgname) {
#     reticulate::source_python(system.file('python', 'brca-exchange-query.py',
#                                           package = 'annotateMaf'),
#                               envir = globalenv())
# }

# Add SystemRequirements: Python (>= 2.7.0) back to DESCRIPTION if re-activating