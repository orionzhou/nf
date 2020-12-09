#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

ps <- ArgumentParser(description = 'Make predictions using trained model on given dataset')
ps$add_argument("fm", help = "(ML) model file")
ps$add_argument("fi", help = "input dataset")
ps$add_argument("fo", help = "output file to save predictions")
args <- ps$parse_args()

require(tidyverse)
require(tidymodels)
##require(progress)

ti = read_tsv(args$fi)
mdl = readRDS(args$fm)

model_predict <- function(mdl, ti) {
  #{{{
  #pb$tick()
  features = names(mdl$fit$variable.importance)
  new.cols = features[!features %in% colnames(ti)]
  if (length(new.cols) > 0) {
    ti2 = matrix(0, nrow=nrow(ti), ncol=length(new.cols), dimnames=list(NULL, new.cols))
    ti = ti %>% bind_cols(as.data.frame(ti2))
  }
  x1 = predict(mdl, ti, type='prob') %>% select(prob = .pred_1)
  x2 = predict(mdl, ti, type='class') %>% select(pred=1)
  tibble(gid=ti$gid, pred=x2$pred, prob=x1$prob)
  #}}}
}

#pb = progress_bar$new(total=nrow(ti))
#pb$tick()
to = model_predict(mdl, ti)

write_tsv(to, args$fo)



