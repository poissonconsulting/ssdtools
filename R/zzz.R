#    Copyright 2015 Province of British Columbia
#
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
#
#       https://www.apache.org/licenses/LICENSE-2.0
#
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.

.onLoad <- function(...) {
  register_s3_method("stats", "nobs", "fitdist")
  register_s3_method("stats", "nobs", "fitdists")
  register_s3_method("stats", "nobs", "fitdistcens")
  register_s3_method("stats", "coef", "fitdists")
  register_s3_method("stats", "predict", "fitdists")
  register_s3_method("graphics", "plot", "fitdists")
  register_s3_method("ggplot2", "autoplot", "fitdist")
  register_s3_method("ggplot2", "autoplot", "fitdists")
  register_s3_method("ggplot2", "autoplot", "fitdistcens")
  invisible()
}

register_s3_method <- function(pkg, generic, class) {
  stopifnot(is.character(pkg), length(pkg) == 1)
  stopifnot(is.character(generic), length(generic) == 1)
  stopifnot(is.character(class), length(class) == 1)

  fun <- get(paste0(generic, ".", class), envir = parent.frame())

  if (pkg %in% loadedNamespaces()) {
    registerS3method(generic, class, fun, envir = asNamespace(pkg))
  }

  setHook(
    packageEvent(pkg, "onLoad"),
    function(...) {
      registerS3method(generic, class, fun, envir = asNamespace(pkg))
    }
  )
}

.onUnload <- function (libpath) {
  library.dynam.unload("ssdtools", libpath)
}

