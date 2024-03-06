#include <math.h>
#include <stdio.h>

double truncated_hernquist_energy(double t, double *pars, double *q, int n_dim) {
    /*  pars:
            G : Gravitational constant
            m : Hernquist M constant
            c : Hernquist a scale length
         rmax : Cutoff radius
    */
    double R;
    R = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);
    if (R > pars[3]) {
        double menc = pars[1] * pars[3] * pars[3] / pow((pars[3] + pars[2]),2);
        return -pars[0] * menc / R;
    } else {
        return -pars[0] * pars[1] / (R + pars[2]);
    }
}

void truncated_hernquist_gradient(double t, double *pars, double *q, int n_dim, double *grad) {
    /*  pars:
            G : Gravitational constant
            m : Hernquist M constant
            c : Hernquist a scale length
         rmax : Cutoff radius
    */
    double R, fac;
    R = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);

    if (R > pars[3]) {
        double menc = pars[1] * pars[3] * pars[3] / pow((pars[3] + pars[2]),2);
        fac = pars[0] * menc / (R*R*R);

        grad[0] = grad[0] + fac*q[0];
        grad[1] = grad[1] + fac*q[1];
        grad[2] = grad[2] + fac*q[2];
    } else {
        fac = pars[0] * pars[1] / ((R + pars[2]) * (R + pars[2]) * R);

        grad[0] = grad[0] + fac*q[0];
        grad[1] = grad[1] + fac*q[1];
        grad[2] = grad[2] + fac*q[2];
    }
}

double truncated_hernquist_density(double t, double *pars, double *q, int n_dim) {
    /*  pars:
            G : Gravitational constant
            m : Hernquist M constant
            c : Hernquist a scale length
         rmax : Cutoff radius
    */
    double R, rho0;
    R = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);

    if (R > pars[3]) {
        return 0.0;
    } else {
        rho0 = pars[1]/(2*M_PI*pars[2]*pars[2]*pars[2]);
        return rho0 / ((R/pars[2]) * pow(1+R/pars[2],3));
    }
}
