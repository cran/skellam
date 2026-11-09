## --------------------------------------------------------------------------
## utils.R -- internal helpers shared by the distribution functions
##
## Internal functions (not exported):
##   .skellam_recycle()   recycle arguments to a common length, warning
##                        as base R does when lengths are not multiples
##   .skellam_bad_rate()  flag rates that are not finite and non-negative
##   .skellam_nonint()    flag values that are not (near-)integers
##
## Original author: Montasser Ghachem and Oguz Ersan, 2026
## License:         GPL (>= 2)
## Package:         skellam, https://github.com/monty-se/skellam
## --------------------------------------------------------------------------

# Recycle all arguments to the longest length. Warn, like `1:3 + 1:2`, when
# a non-empty argument's length does not divide the longest length.
# A zero-length argument makes the result zero-length (as in dpois).
.skellam_recycle <- function(...) {
  args <- list(...)
  lens <- lengths(args)
  if (any(lens == 0L)) return(lapply(args, function(a) a[0]))
  len <- max(lens)
  if (any(len %% lens != 0L))
    warning("longer object length is not a multiple of shorter object length", call. = FALSE)
  lapply(args, rep_len, length.out = len)
}

# TRUE where a rate is present but not a finite non-negative number.
# Missing rates (NA) are not "bad": they propagate as NA without a warning.
.skellam_bad_rate <- function(lambda) {
  !is.na(lambda) & !(is.finite(lambda) & lambda >= 0)
}

# TRUE where x is finite and further than 1e-7 (relative) from an integer,
# the tolerance base R uses for dpois() and friends.
.skellam_nonint <- function(x) {
  is.finite(x) & abs(x - round(x)) > 1e-7 * pmax(1, abs(x))
}
