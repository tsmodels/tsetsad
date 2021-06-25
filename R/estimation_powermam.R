powermam_constraint_violations <- function(x, etsenv)
{
    names(x) <-  etsenv$estimation_names
    parmatrix <- etsenv$parmatrix
    parmatrix[names(x),"init"] <- x
    v <- FALSE
    if (parmatrix["beta","estimate"] == 1) {
        if (parmatrix["beta","init"] > parmatrix["alpha","init"]) {
            v <- TRUE
            return(v)
        }
    }
    return(v)
}

prepare_inputs_powermam_ad <- function(spec, solver = "nloptr")
{
    y <- c(0, as.numeric(spec$target$y))
    x <- rbind(matrix(0, ncol = ncol(spec$xreg$xreg), nrow = 1), spec$xreg$xreg)
    good <- c(0, rep(1, NROW(y) - 1))
    if (any(is.na(y))) {
        good[which(is.na(y))] <- 0
        y[which(is.na(y))] <- 0
    }
    pmatrix <- spec$model$parmatrix
    seasonal_frequency <- spec$seasonal$frequency
    seasonal_normalized <- as.integer(spec$model$normalized_seasonality)
    seasonal_include <- as.integer(spec$model$include_seasonal)
    fmat <- seasonal_matrix(seasonal_frequency)
    alpha <- pmatrix["alpha", 1]
    beta <- pmatrix["beta", 1]
    gamma <- pmatrix["gamma", 1]
    phi <- pmatrix["phi", 1]
    theta <- pmatrix["theta", 1]
    delta <- pmatrix["delta", 1]
    l0 <- as.numeric(pmatrix["l0",1])
    b0 <- as.numeric(pmatrix["b0",1])
    s0 <- unname(pmatrix[grepl("s[0-9]",rownames(pmatrix)), 1])
    b <- as.numeric(pmatrix[grepl("rho[0-9]",rownames(pmatrix)), 1])
    fixed <- list()
    if (pmatrix["alpha","estimate"] == 0) {
        fixed$alpha <- factor(NA)
    }
    if (pmatrix["beta","estimate"] == 0) {
        fixed$beta <- factor(NA)
    }
    if (pmatrix["gamma","estimate"] == 0) {
        fixed$gamma <- factor(NA)
    }
    if (pmatrix["phi","estimate"] == 0) {
        fixed$phi <- factor(NA)
    }
    if (pmatrix["theta","estimate"] == 0) {
        fixed$theta <- factor(NA)
    }
    if (pmatrix["delta","estimate"] == 0) {
        fixed$delta <- factor(NA)
    }
    if (pmatrix["l0","estimate"] == 0) {
        fixed$level_init <- factor(NA)
    }
    if (pmatrix["b0","estimate"] == 0) {
        fixed$b0 <- factor(NA)
    }
    if (spec$model$seasonal_init == "fixed" | spec$model$include_seasonal == 0) {
        s_fixed <- rep(factor(NA), length(s0))
        fixed$s0 <- s_fixed
    }
    rho_coef <- pmatrix[grepl("rho[0-9]",rownames(pmatrix)),,drop = FALSE]
    if (any(rho_coef[,"estimate"] == 0)) {
        ix <- which(rho_coef[,"estimate"] == 0)
        b_fixed <- rep(1:NROW(rho_coef))
        b_fixed[ix] <- factor(NA)
        b_fixed <- as.factor(b_fixed)
        fixed$b <- b_fixed
    }
    # create the functions for the inequalities
    ineq_fun <- NULL
    ineq_jac <- NULL
    # only certain solvers accept inequalities and inequality jacobians
    if (solver %in% c("nloptr")) {
        # check if included and not fixed
        if (spec$model$parmatrix["beta","estimate"] == 1) {
            ineq_fun <- function(x, fun, etsenv)
            {
                pmatrix <- etsenv$parmatrix
                names(x) <- etsenv$estimation_names
                pmatrix[etsenv$estimation_names,"init"] <- x
                ineq <- rbind(pmatrix["beta","init"] +  1e-4 - pmatrix["alpha","init"])
                return(ineq)
            }
            ineq_jac <- function(x, fun, etsenv)
            {
                jac <- matrix(0, ncol = length(etsenv$estimation_names), nrow = 1)
                colnames(jac) <- etsenv$estimation_names
                jac[1,c("alpha","beta")]  <- c(-1, 1)
                return(jac)
            }
        }
    }
    if (solver %in% c("nloptr")) {
        llh_fun <- function(x, fun, etsenv)
        {
            names(x) <- etsenv$parameter_names
            llh <- fun$fn(x)
            etsenv$llh <- llh
            return(llh)
        }
    } else {
        llh_fun <- function(x, fun, etsenv)
        {
            llh <- etsenv$llh
            if (powermam_constraint_violations(x, etsenv)) {
                llh <- llh + 0.2 * abs(llh)
            } else {
                names(x) <- etsenv$parameter_names
                llh <- fun$fn(x)
            }
            etsenv$llh <- llh
            return(llh)
        }
    }
    grad_fun <- function(x, fun, etsenv)
    {
        names(x) <- etsenv$parameter_names
        grad <- fun$gr(x)
        return(grad)
    }
    hess_fun <- function(x, fun, etsenv)
    {
        names(x) <- etsenv$parameter_names
        hess <- fun$he(x)
        return(hess)
    }
    inc <- which(spec$model$parmatrix[,"estimate"] == 1)
    start <- spec$model$parmatrix[inc,"init"]
    names(start) <- rownames(spec$model$parmatrix[inc,])
    rnames <- rownames(pmatrix[which(pmatrix[,"estimate"] == 1), ])
    lower <- pmatrix[which(pmatrix[,"estimate"] == 1), "lower"]
    upper <- pmatrix[which(pmatrix[,"estimate"] == 1), "upper"]
    names(lower) <- rnames
    names(upper) <- rnames
    out <- list(data = list(model = "powermam", y = y, x = x, good = good, fmat = fmat, 
                            seasonal_include = seasonal_include,
                            seasonal_normalized = seasonal_normalized,
                            trend_include = spec$model$include_trend,
                            seasonal_frequency = seasonal_frequency),
                par_list = list(alpha = alpha, 
                                beta = beta,
                                gamma = gamma, 
                                phi = phi,
                                theta = theta,
                                delta = delta,
                                l0 = l0, 
                                b0 = b0,
                                s0 = s0, b = b),
                map = fixed, llh_fun = llh_fun, grad_fun = grad_fun, hess_fun = hess_fun, 
                ineq_fun = ineq_fun, ineq_jac = ineq_jac, start = start, lower = lower, upper = upper)
    return(out)
}


estimate_powermam_ad <- function(object, solver = "nloptr", control = list(trace = 1), ...)
{
    spec_list <- prepare_inputs_powermam_ad(object, solver = solver)
    other_opts <- list(...)
    if (!is.null(other_opts$silent)) {
        silent <- other_opts$silent
    } else {
        silent <- TRUE
    }
    fun <- try(MakeADFun(data = spec_list$data, parameters = spec_list$par_list, DLL = "tsetsad_TMBExports", 
                         map = spec_list$map, trace = FALSE, silent = silent), silent = FALSE)
    if (inherits(fun, 'try-error')) {
        stop("\nestimate_ad found an error. Please use non ad version of estimator and contact developer with reproducible example.")
    }
    etsenv <- new.env()
    etsenv$llh <- 1
    etsenv$grad <- NULL
    etsenv$parameter_names <- names(fun$par)
    etsenv$estimation_names <- rownames(object$model$parmatrix[which(object$model$parmatrix[,"estimate"] == 1 ),])
    etsenv$parmatrix <- object$model$parmatrix
    
    if (solver == "nlminb") {
        sol <- nlminb(start = spec_list$start, objective = spec_list$llh_fun, 
                      gradient = spec_list$grad_fun, hessian = spec_list$hess_fun, 
                      lower = spec_list$lower,
                      upper = spec_list$upper,  control = control, 
                      fun = fun, etsenv = etsenv)
        pars <- sol$par
        llh <- spec_list$llh_fun(pars, fun, etsenv)
        gradient <- spec_list$grad_fun(pars, fun, etsenv)
        hessian <- spec_list$hess_fun(pars, fun, etsenv)
        
        names(pars) <- etsenv$estimation_names
        colnames(gradient) <- etsenv$estimation_names
        colnames(hessian) <- rownames(hessian) <- etsenv$estimation_names
        
        out <- list(pars = pars, llh = llh, gradient = gradient, hessian = hessian, solver_out = sol)
    } else if (solver == "nloptr") {
        if (control$trace == 1) print_level = 1 else print_level = 0
        sol <- nloptr(x0 = spec_list$start, eval_f = spec_list$llh_fun, 
                      eval_grad_f = spec_list$grad_fun, eval_g_ineq = spec_list$ineq_fun,
                      eval_jac_g_ineq = spec_list$ineq_jac, 
                      lb = spec_list$lower,
                      ub = spec_list$upper,  
                      opts = list(algorithm = "NLOPT_LD_SLSQP", maxeval = 5000, print_level = print_level), 
                      fun = fun, etsenv = etsenv)
        pars <- sol$solution
        llh <- spec_list$llh_fun(pars, fun, etsenv)
        gradient <- spec_list$grad_fun(pars, fun, etsenv)
        hessian <- spec_list$hess_fun(pars, fun, etsenv)
        
        names(pars) <- etsenv$estimation_names
        colnames(gradient) <- etsenv$estimation_names
        colnames(hessian) <- rownames(hessian) <- etsenv$estimation_names
        out <- list(pars = pars, llh = llh, gradient = gradient, hessian = hessian, solver_out = sol)
    }
    return(out)
}