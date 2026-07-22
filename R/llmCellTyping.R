###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Annotate clusters with an internal (self-hosted) LLM
##' @description
##' Cluster-level cell type annotation that talks DIRECTLY to an FGCZ-internal,
##' OpenAI-compatible LLM endpoint (e.g. the Qwen3.6 vLLM server on
##' \code{fgcz-c-056}). This is a local-only alternative to
##' \code{CyteTypeR::CyteTypeR()}, which uploads marker data to Nygen's hosted
##' cloud server. No data leaves FGCZ, no API key, no per-call cost.
##'
##' The call pattern mirrors the fgcz-monitoring \code{collect_llm_analysis.py}
##' collector: a single \code{POST} to \code{\{base_url\}/chat/completions} with
##' thinking disabled and \code{response_format=json_object}, so the whole token
##' budget goes to a structured answer.
##'
##' The returned data.frame uses the SAME column names as
##' \code{scData@misc$cytetype_results} so existing downstream mapping and the
##' ScSeurat.Rmd reporting work unchanged.
##'
##' @param markers data.frame of cluster markers (e.g. from
##'   \code{Seurat::FindAllMarkers}). Must contain columns \code{cluster},
##'   \code{gene} and \code{avg_log2FC}. Expected to be pre-filtered to positive
##'   markers; the top \code{n_top_genes} per cluster (by \code{avg_log2FC}) are
##'   sent to the model.
##' @param study_context character string describing the biological context
##'   (species, tissue, disease, experimental conditions). Strongly improves
##'   annotation quality.
##' @param cluster_metadata optional named list (one entry per cluster) of
##'   metadata composition to give the model extra context. Each entry is a
##'   named list of \code{column -> (value -> percent)}. Passed through as-is
##'   from \code{PrepareCyteTypeR}'s \code{clusterMetadata} if available.
##' @param n_top_genes integer, number of top marker genes per cluster to send.
##' @param base_url OpenAI-compatible base URL (must end in \code{/v1}).
##'   Defaults to env var \code{FGCZ_LLM_BASE_URL} or the Qwen3.6 server.
##' @param model model id served by the endpoint. Defaults to env var
##'   \code{FGCZ_LLM_MODEL} or \code{"Qwen3.6-35B-A3B-FP8"}.
##' @param api_key bearer token. vLLM does not enforce auth; defaults to env var
##'   \code{FGCZ_LLM_API_KEY} or \code{"dummy"}.
##' @param timeout request timeout in seconds.
##' @param temperature sampling temperature (low = deterministic).
##' @param max_tokens maximum completion tokens.
##' @param enable_thinking logical; keep \code{FALSE} for reasoning models so the
##'   token budget is spent on the JSON answer rather than reasoning traces.
##' @return data.frame with columns \code{clusterId}, \code{annotation},
##'   \code{granularAnnotation}, \code{ontologyTerm}, \code{cellState},
##'   \code{justification}, \code{supportingMarkers}, \code{conflictingMarkers}.
##'   One row per cluster present in \code{markers}.
##' @author Functional Genomics Center Zurich
##' @export
##' @examples
##' \dontrun{
##' markers <- Seurat::FindAllMarkers(scData, only.pos = TRUE)
##' res <- ezAnnotateClustersLLM(
##'   markers,
##'   study_context = "Human PBMCs from healthy donors, 10x 3' scRNA-seq"
##' )
##' }
ezAnnotateClustersLLM <- function(markers,
                                  study_context = "",
                                  cluster_metadata = NULL,
                                  n_top_genes = 15,
                                  base_url = Sys.getenv("FGCZ_LLM_BASE_URL", "http://fgcz-c-056:8081/v1"),
                                  model = Sys.getenv("FGCZ_LLM_MODEL", "Qwen3.6-35B-A3B-FP8"),
                                  api_key = Sys.getenv("FGCZ_LLM_API_KEY", "dummy"),
                                  timeout = 300,
                                  temperature = 0.1,
                                  max_tokens = 12000,
                                  enable_thinking = FALSE) {
  stopifnot(is.data.frame(markers))
  required_cols <- c("cluster", "gene", "avg_log2FC")
  missing_cols <- setdiff(required_cols, colnames(markers))
  if (length(missing_cols) > 0) {
    stop("markers is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  markers$cluster <- as.character(markers$cluster)
  cluster_ids <- unique(markers$cluster)
  if (length(cluster_ids) == 0) {
    stop("No clusters found in markers.")
  }

  ## Top-N marker genes per cluster (by descending avg_log2FC) --------------
  top_genes <- lapply(cluster_ids, function(cl) {
    sub <- markers[markers$cluster == cl, , drop = FALSE]
    sub <- sub[order(-sub$avg_log2FC), , drop = FALSE]
    head(sub$gene, n_top_genes)
  })
  names(top_genes) <- cluster_ids

  ## Build the per-cluster digest the model reasons over --------------------
  digest_lines <- vapply(cluster_ids, function(cl) {
    line <- sprintf('- clusterId "%s": %s', cl, paste(top_genes[[cl]], collapse = ", "))
    meta <- cluster_metadata[[cl]]
    if (!is.null(meta) && length(meta) > 0) {
      meta_str <- paste(
        vapply(names(meta), function(k) {
          vals <- meta[[k]]
          sprintf("%s=%s", k, paste(sprintf("%s:%s%%", names(vals), unlist(vals)), collapse = "/"))
        }, character(1)),
        collapse = "; "
      )
      line <- paste0(line, "  [composition: ", meta_str, "]")
    }
    line
  }, character(1))

  system_prompt <- paste(
    "You are an expert computational biologist specialising in single-cell",
    "RNA-seq cell type annotation. You are given, for each cluster, its top",
    "marker genes (ranked by log2 fold-change) and the biological study",
    "context. Assign the most specific defensible cell type to each cluster.",
    "Ground every call strictly in the supplied marker genes and study context;",
    "never invent genes that were not provided. Prefer canonical lineage markers",
    "over ambiguous ones. If the markers do not support a confident call, say so",
    "in the annotation (e.g. 'Unresolved / low-confidence') rather than guessing."
  )

  schema_instruction <- paste0(
    "Return ONLY a JSON object with this exact shape, no prose outside it:\n",
    "{\n",
    '  "annotations": [\n',
    "    {\n",
    '      "clusterId": "<the cluster id exactly as given>",\n',
    '      "annotation": "<concise canonical cell type>",\n',
    '      "granularAnnotation": "<finer subtype / cell state, or same as annotation>",\n',
    '      "ontologyTerm": "<Cell Ontology id like CL:0000624 if confident, else empty string>",\n',
    '      "cellState": "<activation/cycling/exhausted/etc. or empty string>",\n',
    '      "justification": "<1-2 sentences citing the supporting markers>",\n',
    '      "supportingMarkers": "<comma-separated genes from the input that support the call>",\n',
    '      "conflictingMarkers": "<comma-separated input genes that conflict, or empty string>"\n',
    "    }\n",
    "  ]\n",
    "}\n",
    "Produce exactly one entry per cluster listed below. Use the cluster ids verbatim."
  )

  user_prompt <- paste(
    schema_instruction,
    paste0("Study context: ", if (nzchar(study_context)) study_context else "not specified"),
    "Clusters and their top marker genes:",
    paste(digest_lines, collapse = "\n"),
    sep = "\n\n"
  )

  ## --- LLM call (OpenAI-compatible /chat/completions) ---------------------
  raw <- .ezLLMChat(
    base_url = base_url, model = model, api_key = api_key,
    system_prompt = system_prompt, user_prompt = user_prompt,
    temperature = temperature, max_tokens = max_tokens,
    enable_thinking = enable_thinking, timeout = timeout
  )

  ## --- Parse + normalise into a cytetype_results-shaped data.frame --------
  annotations <- .ezParseLLMAnnotations(raw)
  if (length(annotations) == 0) {
    stop("LLM returned no parseable annotations. Raw response head:\n",
         substr(raw, 1, 500))
  }

  out_cols <- c("clusterId", "annotation", "granularAnnotation", "ontologyTerm",
                "cellState", "justification", "supportingMarkers", "conflictingMarkers")
  rows <- lapply(annotations, function(a) {
    vals <- lapply(out_cols, function(col) .ezFlattenField(a[[col]]))
    names(vals) <- out_cols
    as.data.frame(vals, stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, rows)
  res$clusterId <- as.character(res$clusterId)

  ## Remap returned ids back to the canonical input ids. Models often echo the
  ## label ("Cluster 0") instead of the bare id ("0"); normalise both sides
  ## (drop a leading "cluster" prefix, lowercase, trim) and match.
  .norm <- function(x) tolower(trimws(sub("^cluster[ _]*", "", x, ignore.case = TRUE)))
  canon <- setNames(cluster_ids, .norm(cluster_ids))
  remapped <- canon[.norm(res$clusterId)]
  res$clusterId[!is.na(remapped)] <- remapped[!is.na(remapped)]
  res <- res[!duplicated(res$clusterId), , drop = FALSE]

  ## Ensure every input cluster has a row; backfill any the model dropped ---
  missing_clusters <- setdiff(cluster_ids, res$clusterId)
  if (length(missing_clusters) > 0) {
    futile.logger::flog.warn(
      "LLM did not annotate cluster(s): %s -- marked 'Unknown Annotation'",
      paste(missing_clusters, collapse = ", ")
    )
    filler <- data.frame(
      clusterId = missing_clusters,
      annotation = "Unknown Annotation",
      granularAnnotation = "Unknown Annotation",
      ontologyTerm = "", cellState = "",
      justification = "No annotation returned by the model for this cluster.",
      supportingMarkers = "", conflictingMarkers = "",
      stringsAsFactors = FALSE
    )
    res <- rbind(res, filler)
  }

  ## Keep input cluster order
  res <- res[match(cluster_ids, res$clusterId), , drop = FALSE]
  rownames(res) <- NULL
  res
}


##' @title Low-level OpenAI-compatible chat completion (internal)
##' @description POSTs a single chat completion to an OpenAI-compatible endpoint
##'   and returns the assistant message content. Disables vLLM "thinking" and
##'   requests a JSON object; on HTTP 400 (endpoint rejects the vLLM extensions)
##'   it retries with a plain body.
##' @keywords internal
.ezLLMChat <- function(base_url, model, api_key, system_prompt, user_prompt,
                       temperature = 0.1, max_tokens = 12000,
                       enable_thinking = FALSE, timeout = 300) {
  url <- paste0(sub("/$", "", base_url), "/chat/completions")
  base_body <- list(
    model = model,
    messages = list(
      list(role = "system", content = system_prompt),
      list(role = "user", content = user_prompt)
    ),
    temperature = temperature,
    max_tokens = max_tokens
  )
  enriched_body <- c(
    base_body,
    list(
      response_format = list(type = "json_object"),
      chat_template_kwargs = list(enable_thinking = enable_thinking)
    )
  )

  do_post <- function(body) {
    req <- httr2::request(url) |>
      httr2::req_headers(
        "Content-Type" = "application/json",
        Authorization = paste("Bearer", api_key)
      ) |>
      httr2::req_body_json(body, auto_unbox = TRUE) |>
      httr2::req_timeout(timeout) |>
      httr2::req_error(is_error = function(resp) FALSE) # handle status manually
    resp <- httr2::req_perform(req)
    status <- httr2::resp_status(resp)
    if (status >= 400) {
      return(list(ok = FALSE, status = status,
                  body = tryCatch(httr2::resp_body_string(resp), error = function(e) "")))
    }
    payload <- httr2::resp_body_json(resp)
    list(ok = TRUE, content = payload$choices[[1]]$message$content)
  }

  res <- do_post(enriched_body)
  if (!res$ok && res$status == 400) {
    futile.logger::flog.warn(
      "LLM endpoint rejected vLLM extensions (400); retrying with a plain body."
    )
    res <- do_post(base_body)
  }
  if (!res$ok) {
    stop(sprintf("LLM endpoint %s returned HTTP %s: %s",
                 url, res$status, substr(res$body, 1, 300)))
  }
  if (is.null(res$content) || !nzchar(res$content)) {
    stop("LLM endpoint returned an empty message content.")
  }
  res$content
}


##' @title Parse LLM JSON annotation response (internal)
##' @description Tolerates <think> blocks and markdown code fences, then extracts
##'   the \code{annotations} array.
##' @keywords internal
.ezParseLLMAnnotations <- function(content) {
  text <- gsub("(?s)<think>.*?</think>", "", content, perl = TRUE)
  text <- trimws(text)
  text <- sub("^```(json)?", "", text)
  text <- sub("```$", "", trimws(text))
  text <- trimws(text)

  obj <- tryCatch(jsonlite::fromJSON(text, simplifyVector = FALSE),
                  error = function(e) NULL)
  if (is.null(obj)) {
    ## Fallback: grab the outermost {...} block
    m <- regmatches(text, regexpr("(?s)\\{.*\\}", text, perl = TRUE))
    if (length(m) == 1) {
      obj <- tryCatch(jsonlite::fromJSON(m, simplifyVector = FALSE),
                      error = function(e) NULL)
    }
  }
  if (is.null(obj)) return(list())
  ann <- obj$annotations
  if (is.null(ann)) {
    ## Some models return a bare array
    if (is.list(obj) && is.null(names(obj))) ann <- obj else return(list())
  }
  ann
}


##' @title Flatten a possibly-list LLM field to a single string (internal)
##' @keywords internal
.ezFlattenField <- function(x) {
  if (is.null(x) || length(x) == 0) return("")
  if (is.list(x)) x <- unlist(x, use.names = FALSE)
  x <- x[!is.na(x)]
  if (length(x) == 0) return("")
  paste(as.character(x), collapse = ", ")
}
