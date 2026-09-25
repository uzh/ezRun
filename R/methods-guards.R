## Guards on the LLM-written Methods text, run by write_methods() on the Description:
## numbers the run's configuration does not contain, steps whose parameter was off,
## and a template Description for when the writer fails or keeps failing a guard.

## A standalone number: not glued to letters (GRCh38, CD45, 10x, log2FC, p28409), except
## a unit suffix or %. Thousands commas, scientific forms (1e5, 10^-5) and ranges
## (1-20, 1:20, "1 to 20": each endpoint is its own token) are read as numbers.
METHODS_NUMBER_CORE <- paste0("(?:\\d{1,3}(?:,\\d{3})+(?:\\.\\d+)?|\\d+(?:\\.\\d+)*)",
                              "(?:[eE][-+]?\\d+|\\^[-+]?\\d+)?")
METHODS_NUMBER_TOKEN <- paste0("(?<![A-Za-z0-9_.])(?:[vV](?=\\d))?", METHODS_NUMBER_CORE,
                               "(?:%|\\s?(?:kb|Kb|bp|Mb|Gb)(?![A-Za-z0-9_])|(?:k|K|M)(?![A-Za-z0-9_]))?",
                               "(?![A-Za-z0-9_])")
METHODS_UNIT_SCALE <- c(kb = 1e3, Kb = 1e3, bp = 1, Mb = 1e6, Gb = 1e9, k = 1e3, K = 1e3, M = 1e6)

## The values a number token can stand for: "5%" is 5 or 0.05, "2 kb" is 2000.
methodsNumberValues <- function(tok) {
  x <- sub("^[vV]", "", tok)
  if (endsWith(x, "%")) {
    v <- methodsNumberValues(sub("%$", "", x))
    return(c(v, v / 100))
  }
  unit <- regmatches(x, regexpr("(kb|Kb|bp|Mb|Gb|k|K|M)$", x))
  scale <- if (length(unit)) METHODS_UNIT_SCALE[[unit]] else 1
  x <- gsub(",", "", trimws(sub("(kb|Kb|bp|Mb|Gb|k|K|M)$", "", x)))
  v <- if (grepl("^", x, fixed = TRUE)) {
    p <- as.numeric(strsplit(x, "^", fixed = TRUE)[[1]])
    p[1]^p[2]
  } else suppressWarnings(as.numeric(x))
  v * scale
}

## Every numeric value in the configuration text; loose on purpose (Release_48 gives 48,
## "2kb" gives 2 and 2000), so a configured value is never flagged for its spelling.
methodsConfigValues <- function(config_text) {
  toks <- regmatches(config_text, gregexpr(paste0("(?<![0-9.])", METHODS_NUMBER_CORE,
                                                  "(?:%|\\s?(?:kb|Kb|bp|Mb|Gb|k|K|M))?"),
                                           config_text, perl = TRUE))[[1]]
  vals <- unlist(lapply(toks, function(t) c(methodsNumberValues(t),
                                            methodsNumberValues(sub("[^0-9]+$", "", t)))))
  unique(vals[is.finite(vals)])
}

## Numbers in the Description that the run's configuration does not contain. config_text is
## what configured the run (job script, parameters, app defaults, config.csv, app facts,
## citation candidates); all_text adds the logs and is only used for versions, because logs
## hold results: nearly every number an LLM writes appears somewhere in them.
methods_check_numbers <- function(description, config_text, all_text, sample_count, allow_small = 0:10) {
  text <- paste(description, collapse = "\n")
  m <- gregexpr(METHODS_NUMBER_TOKEN, text, perl = TRUE)[[1]]
  if (m[1] == -1) return(character(0))
  toks <- regmatches(text, list(m))[[1]]
  config <- methodsConfigValues(paste(config_text, collapse = "\n"))
  all <- paste(all_text, collapse = "\n")
  flagged <- vapply(seq_along(toks), function(i) {
    tok <- toks[i]
    core <- sub("^[vV]", "", tok)
    before <- substr(text, max(1, m[i] - 8), m[i] - 1)
    ## x.y.z, v-prefixed, or after "version" / "R": a version, looked up as a string
    if (grepl("^\\d+(\\.\\d+){2,}$", core) || grepl("^[vV]\\d", tok) ||
        grepl("(version|\\bR)[ -]$", before)) {
      return(!grepl(paste0("(?<![0-9])", gsub(".", "\\.", core, fixed = TRUE), "(?![0-9])"),
                    all, perl = TRUE))
    }
    v <- methodsNumberValues(tok)
    v <- v[is.finite(v)]
    if (!length(v)) return(FALSE)
    if (any(v == round(v) & v %in% c(allow_small, sample_count))) return(FALSE)
    !any(vapply(v, function(x) any(abs(config - x) <= 1e-9 * pmax(abs(config), abs(x))), logical(1)))
  }, logical(1))
  unique(toks[flagged])
}
