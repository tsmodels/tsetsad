/// @file aaa.hpp
#ifndef aaa_hpp
#define aaa_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type aaa(objective_function<Type>* obj) {
    DATA_VECTOR(y); // data vector
    DATA_MATRIX(x); // data vector
    DATA_MATRIX(fmat);
    DATA_VECTOR(good);
    DATA_INTEGER(seasonal_frequency);
    DATA_INTEGER(seasonal_normalized);
    DATA_INTEGER(seasonal_include);
    PARAMETER(alpha);
    PARAMETER(beta);
    PARAMETER(gamma);
    PARAMETER(phi);
    PARAMETER(l0);
    PARAMETER(b0);
    PARAMETER_VECTOR(s0);
    PARAMETER_VECTOR(b);
    const int timesteps = y.rows();
    vector<Type> level(timesteps);
    vector<Type> slope(timesteps);
    matrix<Type> seasonal(timesteps, seasonal_frequency);
    vector<Type> fit(timesteps);
    vector<Type> error(timesteps);
    vector<Type> mu(timesteps);
    error.setZero();
    level.setZero();
    slope.setZero();
    seasonal.setZero();
    fit.setZero();
    level(0) = l0;
    slope(0) = b0;
    for(int j = 0;j< (seasonal_frequency - 1);j++){
        seasonal(0,j) = s0(j);
    }
    //seasonal.block(0,0,0,seasonal.cols() - 1) = seasonal_init;
    seasonal(0, seasonal_frequency - 1) = (Type(-1.0) * Type(s0.sum()));
    Type llh = 0;
    Type a = 0.0;
    vector<Type> onestmp(seasonal_frequency);
    onestmp.setZero();
    Type seasonal_update = 0.0;
    vector<Type> xreg = x * b;
    vector<Type> tmp(seasonal_frequency);
    vector<Type> stmp(seasonal_frequency);
    for(int i = 1;i<timesteps;i++){
        mu(i) = level(i-1) + phi * slope(i-1);
        fit(i) =  mu(i) + seasonal(i-1, seasonal_frequency-1) + xreg(i);
        if (good(i) > 0.5) {
            error(i) = y(i) - fit(i);
        } else {
            error(i) = 0.0;
        }
        if(seasonal_normalized == 1) {
            a = (gamma/seasonal_frequency) * error(i);
        } else {
            a = 0.0;
        }
        level(i) = (mu(i) + a) + alpha * error(i);
        slope(i) = phi * slope(i-1) + beta * error(i);
        onestmp.setConstant(a);
        if (seasonal_include == 1) {
            seasonal_update = (seasonal(i - 1, seasonal_frequency - 1)  - a) + gamma * error(i) ;
            stmp = seasonal.row(i - 1);
            tmp = (fmat * stmp);
            tmp = tmp - onestmp;
            seasonal.row(i) = tmp;
            seasonal(i,0) = seasonal_update;
        }
        llh += CppAD::CondExpGe(good(i), Type(0.5), error(i) * error(i), Type(0.0));
    }
    Type good_timesteps = good.sum();
    llh = (good_timesteps) * log(llh);
    return(llh);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
