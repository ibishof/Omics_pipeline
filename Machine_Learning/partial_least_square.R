
###############################################################################

proteinGroups.numeric <- proteinGroups[,sapply(proteinGroups, is.numeric)]


library(pls)
pls.model = plsr(Diagnosis ~ ., data = proteinGroups, validation = "CV")

# Find the number of dimensions with lowest cross validation error
cv = RMSEP(pls.model)
best.dims = which.min(cv$val[estimate = "adjCV", , ]) - 1

# Rerun the model
pls.model = plsr(Diagnosis ~ ., data = proteinGroups, ncomp = best.dims)

coefficients = coef(pls.model)
sum.coef = sum(sapply(coefficients, abs))
coefficients = coefficients * 100 / sum.coef
coefficients = sort(coefficients[, 1 , 1])
barplot(tail(coefficients, 5))


saveRDS(pls.model, "partial_least_square_modle.rds")
