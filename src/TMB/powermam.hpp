/// @file powermam.hpp
#ifndef powermam_hpp
#define powermam_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type powermam(objective_function<Type>* obj) {
    DATA_VECTOR(y); // data vector
    DATA_MATRIX(x); // data vector
    DATA_MATRIX(fmat);
    DATA_VECTOR(good);
    DATA_INTEGER(seasonal_frequency);
    DATA_INTEGER(seasonal_include);
    DATA_INTEGER(trend_include);
    PARAMETER(alpha);
    PARAMETER(beta);
    PARAMETER(gamma);
    PARAMETER(phi);
    PARAMETER(theta);
    PARAMETER(delta);
    PARAMETER(l0);
    PARAMETER(b0);
    PARAMETER_VECTOR(s0);
    PARAMETER_VECTOR(b);
    const int timesteps = y.rows();
    vector<Type> level(timesteps);
    vector<Type> slope(timesteps);
    matrix<Type> seasonal(timesteps, seasonal_frequency);
    vector<Type> fit(timesteps);
    vector<Type> fitpower(timesteps);
    vector<Type> error(timesteps);
    vector<Type> mu(timesteps);
    vector<Type> mu_pow_theta(timesteps);
    vector<Type> s_pow_delta(timesteps);
    vector<Type> s_pow_delta_m1(timesteps);
    error.setZero();
    level.setZero();
    slope.setZero();
    seasonal.setOnes();
    fit.setZero();
    level(0) = l0;
    slope(0) = b0;
    for(int j=0;j<(seasonal_frequency - 1);j++){
        seasonal(0,j) = s0(j);
    }
    //seasonal.block(0,0,0,seasonal.cols() - 1) = seasonal_init;
    seasonal(0, seasonal_frequency - 1) = Type(seasonal_frequency) - Type(s0.sum());
    Type llh = 0.0;
    Type fsum = 0.0;
    Type seasonal_update = 0.0;
    vector<Type> xreg = x * b;
    vector<Type> tmp(seasonal_frequency);
    vector<Type> stmp(seasonal_frequency);
    for(int i = 1;i<timesteps;i++){
        mu(i) = level(i-1) + phi * slope(i-1);
        mu_pow_theta(i) = pow(mu(i), theta);
        s_pow_delta(i) = pow(seasonal(i-1, seasonal_frequency-1), delta);
        s_pow_delta_m1(i) = pow(seasonal(i-1, seasonal_frequency-1), delta - Type(1.0));
        fit(i) =  (mu(i) + xreg(i)) * seasonal(i-1, seasonal_frequency-1);
        fitpower(i) = pow(mu(i) + xreg(i), theta) * s_pow_delta(i);
        if (good(i) > 0.5) {
            error(i) = (y(i) - fit(i))/fitpower(i);
        } else {
            error(i) = 0.0;
        }
        level(i) = mu(i) + alpha * mu_pow_theta(i) * s_pow_delta_m1(i) * error(i);
        if (trend_include == 1){
            slope(i) = phi * slope(i-1) + beta * mu_pow_theta(i) * s_pow_delta_m1(i) * error(i);
        }
        if (seasonal_include == 1) {
            seasonal_update = seasonal(i - 1, seasonal_frequency - 1) + gamma * s_pow_delta(i) * pow(level(i - 1) + phi * slope(i-1), theta - Type(1.0)) * error(i);
            stmp = seasonal.row(i - 1);
            tmp = (fmat * stmp);
            seasonal.row(i) = tmp;
            seasonal(i,0) = seasonal_update;
        }
        llh += CppAD::CondExpGe(good(i), Type(0.5), error(i) * error(i), Type(0.0));
        fsum += CppAD::CondExpGe(good(i), Type(0.5), Type(log(fabs(fitpower(i)))), Type(0.0));
    }
    Type good_timesteps = good.sum();
    llh = (good_timesteps) * log(llh) + Type(2.0) * fsum;
    return(llh);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif