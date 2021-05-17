seasonal_matrix <- function(frequency)
{
    if (frequency < 2) stop("seasonal_matrix: frequency must be > 1")
    Fmat <- matrix(0, frequency, frequency)
    Fmat[1, frequency] <- 1
    Fmat[2:frequency, 1:(frequency - 1)] <- diag(frequency - 1)
    return(Fmat)
}
pars2list <- function(x, etsenv)
{
    names(x) <- etsenv$parameter_names
    parmatrix <- etsenv$parmatrix
    rnames <- rownames(parmatrix)
    parmatrix[names(x),"init"] <- x
    list(alpha = parmatrix["alpha","init"], beta = parmatrix["beta","init"],
         phi = parmatrix["phi","init"], gamma = parmatrix["gamma","init"],
         l0 = parmatrix["l0","init"], b0 = parmatrix["b0","init"],
         s0 = c(parmatrix[grepl("s[1-9]",rnames),"init"]),
         b = c(parmatrix[grepl("rho[1-9]",rnames),"init"]))
}