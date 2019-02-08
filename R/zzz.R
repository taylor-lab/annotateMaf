.onLoad = function(libname, pkgname) {
    reticulate::source_python(system.file('python', 'brca-exchange-query.py',
                                          package = 'annotateMaf'),
                              envir = globalenv())
}
