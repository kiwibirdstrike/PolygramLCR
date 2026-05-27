classification_metrics = function(y, prob, threshold = 0.5) {
  y = as.integer(y)
  pred = as.integer(prob >= threshold)
  tp = sum(y == 1 & pred == 1)
  tn = sum(y == 0 & pred == 0)
  fp = sum(y == 0 & pred == 1)
  fn = sum(y == 1 & pred == 0)

  sensitivity = tp / (tp + fn)
  specificity = tn / (tn + fp)
  accuracy = mean(pred == y)
  balanced_accuracy = mean(c(sensitivity, specificity))
  log_loss = safe_log_loss(y, prob)
  auc = binary_auc(y, prob)

  data.frame(
    accuracy = accuracy,
    balanced_accuracy = balanced_accuracy,
    sensitivity = sensitivity,
    specificity = specificity,
    log_loss = log_loss,
    auc = auc,
    stringsAsFactors = FALSE
  )
}

binary_auc = function(y, prob) {
  y = as.integer(y)
  ranks = rank(prob, ties.method = "average")
  n_pos = as.numeric(sum(y == 1))
  n_neg = as.numeric(sum(y == 0))
  (sum(ranks[y == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}
