BaumWelchBin <- function (x, series.length, m = 2, Q, p,
     initial.prob, maxiter = 500, tol = 1e-05, dig = 3) {

###############################################
#              OPTION PROCESSING              #
###############################################

    if (is.matrix(x) && ncol(x) == 1) x = as.vector(x)

    xnona = x[!is.na(x)]
    stopifnot(all(xnona %in% c(0L,1L)) ||
                  all(xnona %in% c(FALSE,TRUE)))

    n <- length(x)

    if (missing(series.length)) {
        series.length <- n;
    }

    stopifnot (sum(series.length) == n);

    if (!missing(Q)) {
        if (!all.equal(dim(Q), c(m,m)))
            stop("Q must be an m by m matrix");
    }
    else {
        Q <- matrix(0.1/(m-1), ncol = m, nrow = m);
        diag(Q) <- 0.9;
    }

    if (missing(p)) {
       p <- (1:m)/(m+1)
    }

    stopifnot(length(p) == m)

    adjust.initial.prob <- missing(initial.prob);




###############################################
#            VARIABLE DEFINITIONS             #
###############################################


    old.loglik <- -Inf;
    iter <- 0;

    # Updated in E.step (Q as well).
    emission.prob <- matrix(NA_real_, nrow = n, ncol = m)
    phi <- matrix(NA_real_, nrow = n, ncol = m)
    loglik <- NA_real_


###############################################
#           FUNCTION DEFINITIONS              #
###############################################


    initial.steady.state.probabilities <- function () {
    # Compute the steady-state initial probabilities.

        spectre <- eigen(t(Q))

        if (is.complex(spectre$values[1]))
            return(rep(1/m,m))
        if (spectre$values[1] > 0.99)
            return(spectre$vectors[,1] / sum(spectre$vectors[,1]))

        return(rep(1/m,m))

    }

    E.step <- function () {

        # Emission probabilities.
        for (i in 1:m) {
            emission.prob[,i] <<- p[i]*x + (1-p[i])*(1-x)
        }

        # Emission probabilities of NAs are set to 1 for every states.
        emission.prob[is.na(emission.prob)] <<- 1;

        if (adjust.initial.prob) {
            initial.prob <- initial.steady.state.probabilities();
        }

        loglik <<- 0;
        cumulative.n <- 0;
        transitions <- matrix(double(m * m), nrow = m);

        counter = 0;

        for (n.i in series.length) {

            counter = counter+1;

            forwardback <- .Fortran("fwdb",
                as.integer(m),
                as.integer(n.i),
                initial.prob,
                emission.prob[(cumulative.n + 1):(cumulative.n + n.i),], Q,
                double(n.i),
                matrix(double(n.i * m), nrow = n.i),
                matrix(double(m^2), nrow = m), double(1),
                PACKAGE = "HMMt");

            loglik <<- loglik + forwardback[[9]];
            phi[(cumulative.n + 1):(cumulative.n + n.i),] <<-
                forwardback[[7]];
            transitions <- transitions + forwardback[[8]];

            cumulative.n <- cumulative.n + n.i;

        }

        Q <<- transitions / rowSums(transitions);

    }

    # Update the probabilities.
    update.p <- function() {
       colSums(phi*x + (1-phi)*(1-x), na.rm=TRUE) / n
    }



###############################################
#                 MAIN LOOP                   #
###############################################

    for (iter in 1:maxiter) {

        cat(paste("\riteration:", iter))

        # Update Q and emission.prob.
        E.step()
        p <- update.p();

        if (abs(loglik - old.loglik) < tol)
            break

        old.loglik <- loglik

    } # for (iter in 1:maxiter)

    cat("\n")

    if (adjust.initial.prob) {
        initial.prob <- initial.steady.state.probabilities()
    }
        
    vPath <- Viterbi(Q, initial.prob, emission.prob, series.length)

    # Returns an object of thethe class BaumWelchTfit.
    return(list(Q = Q, p = p, ViterbiPath = vPath, phi = phi,
        logL = loglik, iterations = iter))

}
