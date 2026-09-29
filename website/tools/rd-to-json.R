# Converts the package's help pages (man/*.Rd) into src/data/reference.json, which the
# reference pages render with the site's own components:
#
#   Rscript tools/rd-to-json.R          (from website/; needs R and jsonlite)
#
# It runs in GitHub Actions whenever man/ changes (.github/workflows/site-data.yml), and
# the JSON is committed, so Netlify's build never needs R. Do not edit
# reference.json by hand: edit the roxygen comments in R/, run roxygen, then this.
#
# Each topic becomes an object with its name, aliases, title, the R file its docs come
# from, and its sections. Text sections (description, details, value, arguments, ...)
# become a small subset of HTML: <p>, <code>, <a>, <em>, <strong>, <ul>/<ol>/<li>,
# <dl>/<dt>/<dd> and <pre>. Usage and examples stay plain R code, for Shiki to highlight
# at build time. Links to other topics in the package point at their reference pages.
#
# Rd markup it does not know stops the script with the file and the macro, rather than
# being dropped or passed through: a new macro in man/ should fail the pull request that
# introduces it, not quietly vanish from the site.

main <- function(man_dir = file.path("..", "man"), out = file.path("src", "data", "reference.json")) {
  files <- sort(list.files(man_dir, pattern = "[.]Rd$", full.names = TRUE))
  if (!length(files)) stop("no .Rd files in ", man_dir)
  rds <- lapply(files, tools::parse_Rd)
  names(rds) <- basename(files)

  # What the package exports, so the reference can mark the help pages of internal
  # functions (qdecr, for one) as reachable only through QDECR:::.
  pkg <- normalizePath(file.path(man_dir, ".."))
  ns <- parseNamespaceFile(basename(pkg), dirname(pkg))
  exported <- c(ns$exports, apply(ns$S3methods[, 1:2, drop = FALSE], 1, paste, collapse = "."))

  # Every alias, to the topic it lives in, so \link{qdecr_read_p} finds qdecr_read.
  topic_of <- list()
  for (rd in rds) {
    name <- plain(field(rd, "name"))
    for (alias in fields(rd, "alias")) topic_of[[plain(alias)]] <- name
  }

  topics <- lapply(names(rds), function(file) {
    tryCatch(topic(rds[[file]], topic_of, exported), error = function(e) stop(file, ": ", conditionMessage(e), call. = FALSE))
  })
  topics <- topics[order(vapply(topics, `[[`, "", "name"))]

  json <- as.character(jsonlite::toJSON(list(topics = topics), auto_unbox = TRUE, pretty = 2, null = "null"))
  # Two quirks of jsonlite that only make the diffs noisier: empty arrays spread over three
  # lines, and every / escaped as \/ (valid JSON, but every closing tag reads <\/p>).
  json <- gsub("\\[\\s*\n\\s*\n\\s*\\]", "[]", json, perl = TRUE)
  json <- gsub("(?<!\\\\)\\\\/", "/", json, perl = TRUE)
  # Binary mode, so Windows writes \n rather than \r\n and the file does not churn.
  con <- file(out, "wb")
  on.exit(close(con))
  writeLines(enc2utf8(json), con, sep = "\n", useBytes = TRUE)
  message("wrote ", length(topics), " topics to ", out)
}

# ---------- the parse tree ----------

tag <- function(x) {
  t <- attr(x, "Rd_tag")
  if (is.null(t)) "" else t
}

# The top-level sections of a help page named `name` (\name, \alias, ...).
fields <- function(rd, name) Filter(function(x) tag(x) == paste0("\\", name), rd)
field <- function(rd, name) {
  found <- fields(rd, name)
  if (length(found)) found[[1]] else NULL
}

# The text of a node with all markup dropped: for names, titles and link targets.
plain <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.character(x)) return(paste(x, collapse = ""))
  squish(paste(vapply(x, plain, ""), collapse = ""))
}

# The topic's URL segment: its name with dots as hyphens. /reference/print.vw would end in
# what looks like a file extension, and a static host may not add .html to find the page.
# R names cannot contain hyphens, so no two topics can collide.
slug <- function(name) gsub(".", "-", name, fixed = TRUE)

squish <- function(s) trimws(gsub("[[:space:]]+", " ", s))

escape <- function(s) {
  s <- gsub("&", "&amp;", s, fixed = TRUE)
  s <- gsub("<", "&lt;", s, fixed = TRUE)
  gsub(">", "&gt;", s, fixed = TRUE)
}

# ---------- R code: usage and examples ----------

# Usage and examples as R code. \method{print}{vw}(x) is how Rd writes an S3 method; R's
# own help prints it as print(x) under a comment naming the class, and so does this.
code <- function(x) {
  if (is.character(x)) return(paste(x, collapse = ""))
  switch(tag(x),
    "COMMENT" = "",
    "\\method" = ,
    "\\S3method" = paste0("## S3 method for class '", plain(x[[2]]), "'\n", plain(x[[1]])),
    "\\dots" = ,
    "\\ldots" = "...",
    # Its contents open and close with line breaks of their own.
    "\\dontrun" = paste0("## Not run:", code_list(x), "## End(Not run)"),
    "\\donttest" = code_list(x),
    "\\dontshow" = ,
    "\\testonly" = "",
    "\\R" = "R",
    "\\code" = ,
    "\\link" = ,
    "\\var" = code_list(x),
    "RCODE" = ,
    "TEXT" = ,
    "VERB" = paste(x, collapse = ""),
    stop("unsupported markup in code: ", tag(x))
  )
}

code_list <- function(x) paste(vapply(x, code, ""), collapse = "")

# Code with the blank lines at either end dropped and indentation kept.
code_block <- function(x) {
  if (is.null(x)) return(NULL)
  text <- code_list(x)
  text <- sub("^[[:space:]]*\n", "", text)
  sub("[[:space:]]+$", "", text)
}

# ---------- text: everything else ----------

# Block-level elements are set aside as placeholders while the text around them is cut
# into paragraphs, then put back, so a blank line inside a <pre> does not split it.
new_blocks <- function() {
  blocks <- character()
  list(
    add = function(html) {
      blocks[[length(blocks) + 1]] <<- html
      paste0("\n\n\001", length(blocks), "\001\n\n")
    },
    restore = function(s) {
      for (i in seq_along(blocks)) s <- sub(paste0("\001", i, "\001"), blocks[[i]], s, fixed = TRUE)
      s
    }
  )
}

# A text section as HTML paragraphs and blocks.
html <- function(x, topic_of) {
  if (is.null(x)) return(NULL)
  blocks <- new_blocks()
  inline <- render_list(x, topic_of, blocks)
  paragraphs <- trimws(strsplit(inline, "\n[[:space:]]*\n")[[1]])
  paragraphs <- paragraphs[nzchar(paragraphs)]
  out <- vapply(paragraphs, function(p) {
    if (grepl("^\001[0-9]+\001$", p)) p else paste0("<p>", squish(p), "</p>")
  }, "", USE.NAMES = FALSE)
  out <- blocks$restore(paste(out, collapse = ""))
  if (nzchar(out)) out else NULL
}

render <- function(x, topic_of, blocks) {
  if (is.character(x)) return(text(x))
  t <- tag(x)
  inner <- function(node = x) render_list(node, topic_of, blocks)
  if (t == "") return(inner())
  switch(t,
    "TEXT" = text(x),
    "RCODE" = ,
    "VERB" = escape(paste(x, collapse = "")),
    "COMMENT" = "",
    "\\code" = ,
    "\\samp" = ,
    "\\kbd" = ,
    "\\option" = ,
    "\\env" = ,
    "\\command" = ,
    "\\file" = ,
    "\\var" = paste0("<code>", inner(), "</code>"),
    "\\pkg" = ,
    "\\strong" = ,
    "\\bold" = paste0("<strong>", inner(), "</strong>"),
    "\\emph" = ,
    "\\dfn" = ,
    "\\cite" = paste0("<em>", inner(), "</em>"),
    "\\acronym" = inner(),
    "\\sQuote" = paste0("‘", inner(), "’"),
    "\\dQuote" = paste0("“", inner(), "”"),
    "\\dots" = ,
    "\\ldots" = "…",
    "\\R" = "R",
    "\\cr" = "<br>",
    "\\tab" = " ",
    "\\link" = link(x, topic_of, inner()),
    "\\href" = paste0('<a href="', escape(plain(x[[1]])), '">', render_list(x[[2]], topic_of, blocks), "</a>"),
    "\\url" = ,
    "\\email" = {
      target <- escape(plain(x))
      href <- if (t == "\\email") paste0("mailto:", target) else target
      paste0('<a href="', href, '">', target, "</a>")
    },
    "\\eqn" = ,
    "\\deqn" = paste0("<code>", escape(plain(if (length(x) > 1) x[[2]] else x[[1]])), "</code>"),
    "\\preformatted" = blocks$add(paste0("<pre>", escape(code_block(x)), "</pre>")),
    "\\itemize" = blocks$add(list_items(x, "ul", topic_of)),
    "\\enumerate" = blocks$add(list_items(x, "ol", topic_of)),
    "\\describe" = blocks$add(description_list(Filter(function(n) tag(n) == "\\item", x), topic_of)),
    stop("unsupported markup: ", t)
  )
}

# Children in order. A run of two-part \item{term}{description} outside a list (roxygen
# writes a \value's components this way) becomes a description list.
render_list <- function(x, topic_of, blocks) {
  out <- character()
  run <- list()
  flush <- function() {
    if (length(run)) out[[length(out) + 1]] <<- blocks$add(description_list(run, topic_of))
    run <<- list()
  }
  for (node in x) {
    if (tag(node) == "\\item" && length(node) == 2) {
      run[[length(run) + 1]] <- node
    } else if (length(run) && tag(node) == "TEXT" && !nzchar(trimws(node))) {
      next # the whitespace between items
    } else {
      flush()
      out[[length(out) + 1]] <- render(node, topic_of, blocks)
    }
  }
  flush()
  paste(out, collapse = "")
}

# Rd text, escaped. roxygen here is not in Markdown mode, so `code` is written with
# literal backticks and reaches the Rd as text; it is set as code here.
text <- function(x) {
  s <- escape(paste(x, collapse = ""))
  gsub("`([^`]+)`", "<code>\\1</code>", s)
}

# \link{topic} goes to that topic's page when it is in this package. \link[pkg]{topic} and
# anything unknown stay as text: base R's help is not on this site.
link <- function(x, topic_of, label) {
  option <- attr(x, "Rd_option")
  target <- if (!is.null(option) && startsWith(plain(option), "=")) sub("^=", "", plain(option)) else plain(x)
  if (!is.null(option) && !startsWith(plain(option), "=")) return(label)
  name <- topic_of[[target]]
  if (is.null(name)) return(label)
  paste0('<a href="/reference/', name, '">', label, "</a>")
}

# \itemize and \enumerate: an \item marker, then that item's content, up to the next one.
list_items <- function(x, element, topic_of) {
  items <- list()
  for (node in x) {
    if (tag(node) == "\\item") {
      items[[length(items) + 1]] <- list()
    } else if (length(items)) {
      items[[length(items)]][[length(items[[length(items)]]) + 1]] <- node
    }
  }
  body <- vapply(items, function(item) {
    blocks <- new_blocks()
    paste0("<li>", blocks$restore(squish(render_list(item, topic_of, blocks))), "</li>")
  }, "")
  paste0("<", element, ">", paste(body, collapse = ""), "</", element, ">")
}

description_list <- function(items, topic_of) {
  body <- vapply(items, function(item) {
    blocks <- new_blocks()
    term <- blocks$restore(squish(render_list(item[[1]], topic_of, blocks)))
    desc <- html(item[[2]], topic_of)
    paste0("<dt>", term, "</dt><dd>", if (is.null(desc)) "" else desc, "</dd>")
  }, "")
  paste0("<dl>", paste(body, collapse = ""), "</dl>")
}

# ---------- one topic ----------

topic <- function(rd, topic_of, exported) {
  known <- c(
    "\\name", "\\alias", "\\title", "\\description", "\\usage", "\\arguments", "\\value",
    "\\details", "\\note", "\\author", "\\references", "\\seealso", "\\examples",
    "\\section", "\\keyword", "\\concept", "\\format", "\\source", "\\docType", "\\encoding",
    "COMMENT", "TEXT"
  )
  unknown <- setdiff(unique(vapply(rd, tag, "")), known)
  if (length(unknown)) stop("unsupported section: ", paste(unknown, collapse = ", "))

  # roxygen notes the source file in a comment at the top: "Please edit documentation in
  # R/plotting.R". The reference links to it, since that is where a fix goes.
  comments <- vapply(Filter(function(x) tag(x) == "COMMENT", rd), plain, "")
  source <- regmatches(comments, regexpr("R/[^ ,]+[.][Rr]", comments))

  arguments <- field(rd, "arguments")
  arguments <- if (is.null(arguments)) list() else lapply(
    Filter(function(x) tag(x) == "\\item", arguments),
    function(item) list(name = plain(item[[1]]), html = html(item[[2]], topic_of))
  )

  text_section <- function(name) html(field(rd, name), topic_of)
  sections <- lapply(fields(rd, "section"), function(s) list(title = plain(s[[1]]), html = html(s[[2]], topic_of)))

  aliases <- vapply(fields(rd, "alias"), plain, "")

  name <- plain(field(rd, "name"))

  list(
    name = name,
    slug = slug(name),
    aliases = I(aliases),
    # Whether any function on the page is exported: an internal one is not on the search
    # path after library(QDECR).
    exported = any(aliases %in% exported),
    title = plain(field(rd, "title")),
    source = if (length(source)) source[[1]] else NULL,
    description = text_section("description"),
    usage = code_block(field(rd, "usage")),
    arguments = arguments,
    value = text_section("value"),
    details = text_section("details"),
    sections = sections,
    note = text_section("note"),
    author = text_section("author"),
    references = text_section("references"),
    seealso = text_section("seealso"),
    examples = code_block(field(rd, "examples")),
    keywords = I(vapply(fields(rd, "keyword"), plain, ""))
  )
}

main()
