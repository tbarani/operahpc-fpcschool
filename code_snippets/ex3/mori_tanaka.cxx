const mfem_mgis::real hydrostatic_pressure = XXXXXXX;
const mfem_mgis::real elastic_modulus = 150e-3;
const mfem_mgis::real poisson_modulus = 0.3;
const mfem_mgis::real matrix_bulk_modulus = elastic_modulus / (3. * (1. - 2. * poisson_ratio));
const mfem_mgis::real matrix_shear_modulus = elastic_modulus / (2. * (1. + poisson_ratio));

const mfem_mgis::real effective_bulk_modulus = 
    matrix_bulk_modulus * 
    (1. - porosity) /
    (1. + porosity *
          matrix_bulk_modulus / 
         (matrix_bulk_modulus + 4./3. * matrix_shear_modulus));

// Formula: P = -K * tr(eps) => eps_xx = eps_yy = eps_zz = -P / 3K
const mfem_mgis::real eps_hydro = -hydrostatic_pressure /  (3. * effective_bulk_modulus);

// Fill the flattened tensor
std::vector<mfem_mgis::real> e = {eps_hydro,eps_hydro,eps_hydro,0,0,0};
// We set this “evolution” method by having it returning always the same quantity (no t input)
problem.setMacroscopicGradientsEvolution(
                            [e](const mfem_mgis::real) { return e; });
