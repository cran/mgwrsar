## ----setup, include=FALSE-----------------------------------------------------
target_url <- "https://ggeniaux.github.io/mgwrsar_vignettes/GWR-and-Mixed-GWR-with-spatial-autocorrelation.html"

## ----redirect, echo=FALSE, results='asis'-------------------------------------
cat(sprintf(
  '<p><strong>Online version:</strong> <a href="%s">%s</a></p>',
  target_url, target_url
))

# Meta refresh (works without JS); short delay to keep the page readable
cat(sprintf('<meta http-equiv="refresh" content="3; url=%s">', target_url))

# JS redirect as a fallback
cat(sprintf("<script>setTimeout(function(){ window.location.href = '%s'; }, 3000);</script>", target_url))

