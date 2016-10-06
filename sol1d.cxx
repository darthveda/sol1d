/*
 1D simulation of plasma-neutral interactions
 
  Normalisations
  --------------
 
  Ne   (density) normalised to Nnorm [m^-3]
  T    (temperature) normalised to Tnorm [eV]
  B    (magnetic field) normalised to Bnorm [eV]
  
  t    (time) normalised using ion cyclotron frequency Omega_ci [1/s]
  Vi   (velocity) normalised to sound speed Cs [m/s]
  L    (lengths) normalised to hybrid Larmor radius rho_s = Cs/Omega_ci [m]
  
 */

#include <mpi.h>

#include <bout/physicsmodel.hxx>
#include <bout/constants.hxx>
#include <derivs.hxx>
#include <field_factory.hxx>
#include <invert_parderiv.hxx>

#include "div_ops.hxx"
#include "loadmetric.hxx"
#include "radiation.hxx"

class SD1D : public PhysicsModel {
protected:
  int init(bool restarting) {
    Options *opt = Options::getRoot()->getSection("sol1d");
    // Normalisation
    OPTION(opt, Tnorm, 100);  // Reference temperature [eV]
    OPTION(opt, Nnorm, 1e19); // Reference density [m^-3]
    OPTION(opt, Bnorm, 1.0);  // Reference magnetic field [T]
    OPTION(opt, AA, 2.0);     // Ion mass
    SAVE_ONCE4(Tnorm, Nnorm, Bnorm, AA); // Save normalisations
    
    // Model parameters
    OPTION(opt, vwall, 1./3);   // 1/3rd Franck-Condon energy at wall
    OPTION(opt, frecycle, 1.0); // Recycling fraction 100%
    OPTION(opt, fredistribute, 0.0); // Fraction of neutrals redistributed evenly along leg
    OPTION(opt, sheath_density_drop, true); // Does density drop at the sheath?
    OPTION(opt, gaspuff,  0.0); // Additional gas flux at target
    OPTION(opt, dneut, 1.0);    // Scale neutral gas diffusion
    OPTION(opt, nloss, 0.0);    // Neutral gas loss rate
    OPTION(opt, fimp,  0.01);   // 1% impurity
    OPTION(opt, Eionize,   30);     // Energy loss per ionisation (30eV)
    OPTION(opt, sheath_gamma, 6.5); // Sheath heat transmission
    OPTION(opt, neutral_gamma, 0.25); // Neutral heat transmission
    
    // Plasma anomalous transport
    OPTION(opt, anomalous_D, -1);
    OPTION(opt, anomalous_chi, -1);
    
    if(sheath_gamma < 6) 
      throw BoutException("sheath_gamma < 6 not consistent");

    OPTION(opt, tn_floor, 3.5);  // Minimum neutral gas temperature [eV]

    OPTION(opt, atomic, true);

    OPTION(opt, neutral_f_pn, true);

    OPTION(opt, hyper, -1);            // Numerical hyper-diffusion
    OPTION(opt, ADpar, -1);            // Added Dissipation scheme
    OPTION(opt, viscos, -1);           // Parallel viscosity
    OPTION(opt, ion_viscosity, false); // Braginskii parallel viscosity
    
    // Read the flux-tube area from input file
    // This goes into the Jacobian.
    string area_string;
    FieldFactory ffact(mesh); 
    
    // Calculate normalisation factors
    
    Cs0      = sqrt(SI::qe*Tnorm / (AA*SI::Mp)); // Reference sound speed [m/s]
    Omega_ci = SI::qe*Bnorm / (AA*SI::Mp);       // Ion cyclotron frequency [1/s]
    rho_s0   = Cs0 / Omega_ci;                   // Length scale [m]
  
    mi_me  = AA*SI::Mp/SI::Me;

    BoutReal Coulomb = 6.6 - 0.5*log(Nnorm * 1e-20) + 1.5*log(Tnorm);
    tau_e0 = 1. / (2.91e-6 * (Nnorm / 1e6) * Coulomb * pow(Tnorm, -3./2));
    
    // Save normalisation factors
    SAVE_ONCE5(Cs0, Omega_ci, rho_s0, tau_e0, mi_me);
    
    OPTION(opt, volume_source, true);
    if(volume_source) {
      // Volume sources of particles and energy
      
      string source_string;
    
      Options* optne = Options::getRoot()->getSection("Ne");
      optne->get("source", source_string, "0.0");
      NeSource = ffact.create2D(source_string, optne);
      //SAVE_ONCE(NeSource);
      
      Options* optpe = Options::getRoot()->getSection("P");
      optpe->get("source", source_string, "0.0");
      PeSource = ffact.create2D(source_string, optpe);
      SAVE_ONCE(PeSource);
      
      // Normalise sources
      NeSource /= Nnorm * Omega_ci;
      PeSource /= SI::qe * Nnorm * Tnorm * Omega_ci;
    }else {
      // Point sources, fixing density and specifying energy flux
      
      
      Options* optpe = Options::getRoot()->getSection("P");
      OPTION(optpe, powerflux, 2e7); // Power flux in W/m^2
      powerflux /= rho_s0 * SI::qe * Tnorm * Nnorm * Omega_ci; // Normalised energy flux
    }

    /////////////////////////
    // Density controller
    OPTION(opt, density_upstream, -1); // Fix upstream density? [m^-3]
    if(density_upstream > 0.0) {
      // Fixing density
      density_upstream /= Nnorm;
      
      // Controller
      OPTION(opt, density_controller_p, 1e-2);
      OPTION(opt, density_controller_i, 1e-3);

      density_error_lasttime = -1.0;  // Signal no value

      // Save and load error integral from file, since
      // this determines the source function
      solver->addToRestart(density_error_integral, "density_error_integral");
      
      if(!restarting) {
        density_error_integral = 0.0;
        
        if(volume_source) {
          // Set density_error_integral so that
          // the input source is used
          density_error_integral = 1./density_controller_i;
        }
      }
    }

    if(volume_source) {
      if(density_upstream > 0.0) {
        // Evolving NeSource
        SAVE_REPEAT(NeSource);

        NeSource0 = NeSource; // Save initial value
      }else {
        // Fixed NeSource
        SAVE_ONCE(NeSource);
      }
    }
    
    Options::getRoot()->getSection("NVn")->get("evolve", evolve_nvn, true);
    Options::getRoot()->getSection("Pn")->get("evolve", evolve_pn, true);
    
    nloss /= Omega_ci;

    // Specify variables to evolve
    solver->add(Ne, "Ne");
    solver->add(NVi, "NVi");
    solver->add(P, "P");
    solver->add(Nn, "Nn");
    if(evolve_nvn) {
      solver->add(NVn, "NVn");
    }
    if(evolve_pn) {
      solver->add(Pn, "Pn");
    }
    
    // Load the metric tensor
    LoadMetric(rho_s0, Bnorm);
    
    opt->get("area", area_string, "1.0");
    mesh->J = ffact.create2D(area_string, Options::getRoot());
    
    dy4 = SQ(SQ(mesh->dy));
    
    // Use carbon radiation for the impurity
    rad = new HutchinsonCarbonRadiation();

    // Add extra quantities to be saved
    
    SAVE_REPEAT4(S, R, E, F); // Save net plasma particle source, radiated power, energy transfer, friction
    
    SAVE_REPEAT(kappa_epar); // Save coefficient of thermal conduction
    SAVE_REPEAT2(Dn, kappa_n);  // Neutral diffusion coefficients
    SAVE_REPEAT(flux_ion);   // Flux of ions to target
    
    bool diagnose;
    OPTION(opt, diagnose, true);
    if(diagnose) {
      // Output extra variables
      SAVE_REPEAT2(Srec,Siz);       // Save particle sources
      SAVE_REPEAT3(Frec,Fiz,Fcx);   // Save momentum sources
      SAVE_REPEAT3(Rrec,Riz,Rzrad); // Save radiation sources
      SAVE_REPEAT3(Erec,Eiz,Ecx);   // Save energy transfer
      
      SAVE_REPEAT(Vi);
      if(evolve_nvn) {
        SAVE_REPEAT(Vn);
      }
    }
    
    if(ion_viscosity) 
      SAVE_REPEAT(eta_i);
    
    kappa_epar = 0.0;

    Srec = 0.0; Siz = 0.0; S = 0.0;
    Frec = 0.0; Fiz = 0.0; Fcx = 0.0;   F = 0.0;
    Rrec = 0.0; Riz = 0.0; Rzrad = 0.0; R = 0.0;
    Erec = 0.0; Eiz = 0.0; Ecx = 0.0;   E = 0.0;
    
    flux_ion = 0.0;
    
    // Neutral gas diffusion and heat conduction
    Dn = 0.0;
    kappa_n = 0.0;

    // Anomalous transport
    if(anomalous_D > 0.0) {
      // Normalise
      anomalous_D /= rho_s0*rho_s0*Omega_ci; // m^2/s
      output.write("\tnormalised anomalous D_perp = %e\n", anomalous_D);
    }
    if(anomalous_chi > 0.0) {
      // Normalise
      anomalous_chi /= rho_s0*rho_s0*Omega_ci; // m^2/s
      output.write("\tnormalised anomalous chi_perp = %e\n", anomalous_chi);
    }
    
    // Calculate neutral gas redistribution weights over the domain
    string redist_string;
    opt->get("redist_weight", redist_string, "1.0");
    redist_weight = ffact.create2D(redist_string, opt);
    BoutReal localweight = 0.0;
    for(int j=mesh->ystart;j<=mesh->yend;j++) {
      localweight += redist_weight(mesh->xstart, j);
    }

    MPI_Comm ycomm = mesh->getYcomm(mesh->xstart); // MPI communicator

    // Calculate total weight by summing over all processors
    BoutReal totalweight;
    MPI_Allreduce(&localweight, &totalweight, 1, MPI_DOUBLE, MPI_SUM, ycomm);
    // Normalise redist_weight so sum over domain is 1
    redist_weight /= totalweight;
    
    setPrecon( (preconfunc) &SD1D::precon );

    //////////////////////////////////////////
    // Split operator (IMEX) schemes
    // Use combination of explicit and implicit methods
    //
    // Boolean flags rhs_explicit and rhs_implicit
    // turn on explicit and implicit terms 

    bool split_operator;
    OPTION(opt, split_operator, false);
    if(!split_operator) {
      // Turn on all terms in rhs
      rhs_explicit = rhs_implicit = true;
      update_coefficients = true;
    }
    setSplitOperator(split_operator);
   
    return 0;
  }

  /*!
   * This function calculates the time derivatives
   * of all evolving quantities
   *
   */
  int rhs(BoutReal time) {
    fprintf(stderr, "\rTime: %e", time);

    mesh->communicate(Ne, NVi, P, Nn);
    if(evolve_nvn) {
      mesh->communicate(NVn);
    }
    if(evolve_pn) {
      mesh->communicate(Pn);
    }
    
    // Floor small values
    P = floor(P, 1e-10);
    Ne = floor(Ne, 1e-10);
    Nn = floor(Nn, 1e-10);
    
    Field3D Nelim = floor(Ne, 1e-5);
    Field3D Nnlim = floor(Nn, 1e-5);
    
    Vi = NVi / Ne;
    
    if(evolve_nvn) {
      Vn = NVn / Nnlim;
    }else {
      Vn = - vwall * sqrt(3.5/Tnorm);
      NVn = Nn * Vn;
    }
    
    Field3D Te = 0.5*P / Ne; // Assuming Te = Ti
    
    Field3D Tn;
    if(evolve_pn) {
      Tn = Pn / Nnlim;
      //Tn = floor(Tn, 0.025/Tnorm); // Minimum tn_floor
      Tn = floor(Tn, 1e-12); 
    }else {
      Tn = Te; // Strong CX coupling
      Pn = Tn * floor(Nn, 0.0);
      Tn = floor(Tn, tn_floor/Tnorm); // Minimum of tn_floor
    }
    
    if(update_coefficients) {
      // Update diffusion coefficients
      
      tau_e = Omega_ci * tau_e0 * (Te^1.5)/Ne;

      kappa_epar = 3.2 * mi_me * 0.5*P * tau_e;
      kappa_epar.applyBoundary("neumann");
      
      // Neutral diffusion rate
      
      for(int i=0;i<mesh->ngx;i++)
        for(int j=0;j<mesh->ngy;j++)
          for(int k=0;k<mesh->ngz-1;k++) {
            // Charge exchange frequency, normalised to ion cyclotron frequency
            BoutReal sigma_cx = Nelim(i,j,k)*Nnorm*hydrogen.chargeExchange(Te(i,j,k)*Tnorm)/Omega_ci;

            // Ionisation frequency
            BoutReal sigma_iz = Nelim(i,j,k)*Nnorm*hydrogen.ionisation(Te(i,j,k)*Tnorm)/Omega_ci;
            
            // Neutral thermal velocity
            BoutReal tn = Tn(i,j,k);
            if(tn < tn_floor/Tnorm)
              tn = tn_floor/Tnorm;
            BoutReal vth_n = sqrt(tn); // Normalised to Cs0
            
            // Neutral-neutral mean free path
            BoutReal Lmax = 1.0; // meters
            BoutReal a0 = PI*SQ(5.29e-11);
            BoutReal lambda_nn = 1. / (Nnorm*Nnlim(i,j,k)*a0); // meters
            if(lambda_nn > Lmax) {
              // Limit maximum mean free path
              lambda_nn = Lmax;
            }
            
            lambda_nn /= rho_s0; // Normalised length to rho_s0
            // Neutral-Neutral collision rate
            BoutReal sigma_nn = vth_n / lambda_nn;
	  
            // Total neutral collision frequency
            BoutReal sigma = sigma_cx + sigma_iz + sigma_nn;
            
            
            // Neutral gas diffusion
            if(dneut > 0.0) {
              Dn(i,j,k) = dneut * SQ(vth_n) / sigma;
            }
              
            // Neutral gas heat conduction
            kappa_n(i,j,k) = Nnlim(i,j,k) * SQ(vth_n) / sigma;
          }
      
      kappa_n.applyBoundary("Neumann");
      Dn.applyBoundary("dirichlet_o2");
      mesh->communicate(kappa_n, Dn);
    }
    
    // Set sheath boundary condition on flow
    
    ddt(P) = 0.0; // Need to set heat flux

    if(evolve_pn) {
      ddt(Pn) = 0.0;
    }
    
    for(RangeIterator r=mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      int jz = 0;
      
      // Outward flow velocity to >= Cs
      BoutReal Vout = sqrt(2.0*Te(r.ind, mesh->yend, jz)); // Sound speed outwards
      if(Vi(r.ind, mesh->yend, jz) > Vout)
        Vout = Vi(r.ind, mesh->yend, jz); // If plasma is faster, go to plasma velocity
      
      BoutReal Nout;
      if( sheath_density_drop ) {
        // Zero gradient particle flux N*Vi* J*dx*dz
        // Since Vi increases into the sheath, density should drop
        Nout = Ne(r.ind, mesh->yend, jz) * mesh->J(r.ind, mesh->yend) * Vi(r.ind, mesh->yend, jz) / (0.5*(mesh->J(r.ind, mesh->yend) + mesh->J(r.ind, mesh->yend+1)) * Vout);
      }else {
        // Free boundary on density (constant gradient)
        Nout = 0.5*( 3.*Ne(r.ind, mesh->yend, jz) - Ne(r.ind, mesh->yend-1, jz) );
        // Zero gradient
        //Nout = Ne(r.ind, mesh->yend, jz); 
      }
      
      if(Nout < 0.0)
        Nout = 0.0; // Prevent Nout from going negative -> Flux is always to the wall
      
      // Flux of particles is Ne*Vout
      BoutReal flux = Nout * Vout;

      if(rhs_explicit) {
        // Additional cooling
        BoutReal q = (sheath_gamma - 6) * Te(r.ind, mesh->yend, jz) * flux;
        
        // Multiply by cell area to get power
        BoutReal heatflux = q * (mesh->J(r.ind, mesh->yend)+mesh->J(r.ind, mesh->yend+1))/(sqrt(mesh->g_22(r.ind, mesh->yend)) + sqrt(mesh->g_22(r.ind, mesh->yend+1)));
        
        // Divide by volume of cell, and 2/3 to get pressure
        ddt(P)(r.ind, mesh->yend, jz) -= (2./3)*heatflux / (mesh->dy(r.ind, mesh->yend)*mesh->J(r.ind, mesh->yend));
      }
      
      // Set boundary half-way between cells
      for(int jy=mesh->yend+1; jy<mesh->ngy; jy++) {
        
        ///// Plasma model
        
        // Vi fixed value (Dirichlet)
        Vi(r.ind, jy, jz)  = 2.*Vout - Vi(r.ind, mesh->yend, jz);
        
        // Ne set from flux (Dirichlet)
        Ne(r.ind, jy, jz)  = 2*Nout - Ne(r.ind, mesh->yend, jz);
        
        // NVi. This can be negative, so set this to the flux
        // going out of the domain (zero gradient)
        NVi(r.ind, jy, jz) = Nout * Vout;
        
        // Te zero gradient (Neumann)
        Te(r.ind, jy, jz) = Te(r.ind, mesh->yend, jz);
        
        P(r.ind, jy, jz)   = 2. * Ne(r.ind, jy, jz) * Te(r.ind, jy, jz);
        
        ///// Neutral model
        // Flux of neutral particles, momentum, and energy are set later
        // Here the neutral velocity is set to no-flow conditions
	
        // Vn fixed value (Dirichlet)
        Vn(r.ind, jy, jz)  = - Vn(r.ind, mesh->yend, jz);
        
        // Nn free boundary (constant gradient)
        Nn(r.ind, jy, jz) = 2.*Nn(r.ind, mesh->yend, jz) - Nn(r.ind, mesh->yend-1, jz);
        
        if(evolve_pn) {
          // Tn fixed value (Dirichlet)
          //Tn(r.ind, jy, jz) = 3.5/Tnorm - Tn(r.ind, mesh->yend, jz);
          
          // Tn zero gradient. Heat flux set by gamma
          Tn(r.ind, jy, jz) = Tn(r.ind, mesh->yend, jz);

          if(rhs_explicit && (neutral_gamma > 0.0)) {
            // Density at the target
            BoutReal Nnout = 0.5*(Nn(r.ind, mesh->yend, jz) + Nn(r.ind, mesh->yend+1, jz));
            // gamma * n * T * cs
            BoutReal q = neutral_gamma * Nnout * Tn(r.ind, jy, jz) * sqrt(Tn(r.ind, jy, jz));
            
            // Multiply by cell area to get power
            BoutReal heatflux = q * (mesh->J(r.ind, mesh->yend)+mesh->J(r.ind, mesh->yend+1))/(sqrt(mesh->g_22(r.ind, mesh->yend)) + sqrt(mesh->g_22(r.ind, mesh->yend+1)));
            
            // Divide by volume of cell, and 2/3 to get pressure
            ddt(Pn)(r.ind, mesh->yend, jz) -= (2./3)*heatflux / (mesh->dy(r.ind, mesh->yend)*mesh->J(r.ind, mesh->yend));
          }
        }else {
          Tn(r.ind, jy, jz) = Te(r.ind, jy, jz);
        }
        Pn(r.ind, jy, jz) = Nn(r.ind, jy, jz) * Tn(r.ind, jy, jz); 
        NVn(r.ind, jy, jz) = -NVn(r.ind, mesh->yend, jz);
      }
    }
    
    for(RangeIterator r=mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      // No-flow boundary condition on left boundary
      
      for(int jz=0; jz<mesh->ngz-1; jz++) {
        for(int jy=0; jy<mesh->ystart; jy++) {
          Te(r.ind, jy, jz) = Te(r.ind, mesh->ystart, jz);
          Ne(r.ind, jy, jz) = Ne(r.ind, mesh->ystart, jz);
          P(r.ind, jy, jz) = P(r.ind, mesh->ystart, jz);
          Vi(r.ind, jy, jz) = -Vi(r.ind, mesh->ystart, jz);
          
          Vn(r.ind, jy, jz) = -Vn(r.ind, mesh->ystart, jz);
        }
      }
    }

    if((density_upstream > 0.0) && rhs_explicit) {
      ///////////////////////////////////////////////
      // Set velocity on left boundary to set density
      //
      // This calculates a source needed in the first grid cell, to relax towards
      // the desired density value. 
      //

      BoutReal source;
      for(RangeIterator r=mesh->iterateBndryLowerY(); !r.isDone(); r++) {
        int jz = 0;
          
        // Density source, so dn/dt = source
        BoutReal error = density_upstream - Ne(r.ind, mesh->ystart, jz);
        
        // PI controller, using crude integral of the error
        if(density_error_lasttime < 0.0) {
          // First time
          density_error_lasttime = time;
          density_error_last = error;
        }
        
        // Integrate using Trapezium rule
        if(time > density_error_lasttime) { // Since time can decrease
          density_error_integral += (time - density_error_lasttime)*
            0.5*(error + density_error_last);
        }
        
        
        // Calculate source from combination of error and integral
        source = density_controller_p * error + density_controller_i * density_error_integral;
        
        //output.write("\n Source: %e, %e : %e, %e -> %e\n", time, (time - density_error_lasttime), error, density_error_integral, source);

        density_error_last = error;
        density_error_lasttime = time;

        if(!volume_source) {
          // Convert source into a flow velocity
          // through the boundary, based on a zero-gradient boundary on the density.
          // This ensures that the mass and momentum inputs are consistent,
          // but also carries energy through the boundary. This flux
          // of energy is calculated, and subtracted from the pressure equation,
          // so that the density boundary does not contribute to energy balance.
          
          // Calculate needed input velocity
          BoutReal Vin = source * sqrt(mesh->g_22(r.ind, mesh->ystart))*mesh->dy(r.ind, mesh->ystart) / Ne(r.ind, mesh->ystart, jz);
          
          // Limit at sound speed
          BoutReal cs = sqrt(Te(r.ind, mesh->ystart, jz));
          if( fabs(Vin) > cs ) {
            Vin *= cs / fabs(Vin); // + or - cs
          }
          Vi(r.ind, mesh->ystart-1, jz) = 2.*Vin - Vi(r.ind, mesh->ystart, jz);
          
          // Power flux is v * (5/2 P + 1/2 m n v^2 )
          BoutReal inputflux = Vin * ( 2.5 * P(r.ind, mesh->ystart, jz) + 0.5*Ne(r.ind, mesh->ystart, jz) * Vin*Vin );  // W/m^2 (normalised)
          
          // Subtract input energy flux from P equation
          // so no net power input
          ddt(P)(r.ind, mesh->ystart, jz) -= (2./3)*inputflux / ( mesh->dy(r.ind, mesh->ystart) * sqrt(mesh->g_22(r.ind, mesh->ystart)));
        }
      }

      if(volume_source) {
        if(source < 0.0)
          source = 0.0; // Don't remove particles
        
        // Broadcast the value of source from processor 0
        MPI_Bcast(&source, 1, MPI_DOUBLE, 0, BoutComm::get());
        
        // Scale NeSource
        NeSource = source * NeSource0;
      }
    }
    
    if(atomic && rhs_explicit) {
      // Atomic physics
      
      // Impurity radiation 
      Rzrad = rad->power(Te*Tnorm, Ne*Nnorm, Ne*(Nnorm*fimp)); // J / m^3 / s
      Rzrad /= SI::qe*Tnorm*Nnorm * Omega_ci; // Normalise
      
      E = 0.0; // Energy transfer to neutrals

      // Lower floor on Nn for atomic rates
      Field3D Nnlim2 = floor(Nn, 0.0);
      
      for(int i=0;i<mesh->ngx;i++)
        for(int j=mesh->ystart;j<=mesh->yend;j++)
          for(int k=0;k<mesh->ngz-1;k++) {
            
            // Integrate rates over each cell using Simpson's rule
            // Calculate cell centre (C), left (L) and right (R) values
            
            BoutReal Te_C = Te(i,j,k), Te_L = 0.5*(Te(i,j-1,k) + Te(i,j,k)), Te_R = 0.5*(Te(i,j,k) + Te(i,j+1,k));
            BoutReal Ne_C = Ne(i,j,k), Ne_L = 0.5*(Ne(i,j-1,k) + Ne(i,j,k)), Ne_R = 0.5*(Ne(i,j,k) + Ne(i,j+1,k));
            BoutReal Vi_C = Vi(i,j,k), Vi_L = 0.5*(Vi(i,j-1,k) + Vi(i,j,k)), Vi_R = 0.5*(Vi(i,j,k) + Vi(i,j+1,k));
            BoutReal Tn_C = Tn(i,j,k), Tn_L = 0.5*(Tn(i,j-1,k) + Tn(i,j,k)), Tn_R = 0.5*(Tn(i,j,k) + Tn(i,j+1,k));
            BoutReal Nn_C = Nnlim2(i,j,k), Nn_L = 0.5*(Nnlim2(i,j-1,k) + Nnlim2(i,j,k)), Nn_R = 0.5*(Nnlim2(i,j,k) + Nnlim2(i,j+1,k));
            BoutReal Vn_C = Vn(i,j,k), Vn_L = 0.5*(Vn(i,j-1,k) + Vn(i,j,k)), Vn_R = 0.5*(Vn(i,j,k) + Vn(i,j+1,k));

            // Jacobian (Cross-sectional area)
            BoutReal J_C = mesh->J(i,j), J_L = 0.5*(mesh->J(i,j-1) + mesh->J(i,j)), J_R = 0.5*(mesh->J(i,j) + mesh->J(i,j+1));
            
            ///////////////////////////////////////
            // Charge exchange
            
            BoutReal R_cx_L = Ne_L*Nn_L*hydrogen.chargeExchange(Te_L*Tnorm) * (Nnorm / Omega_ci);
            BoutReal R_cx_C = Ne_C*Nn_C*hydrogen.chargeExchange(Te_C*Tnorm) * (Nnorm / Omega_ci);
            BoutReal R_cx_R = Ne_R*Nn_R*hydrogen.chargeExchange(Te_R*Tnorm) * (Nnorm / Omega_ci);
            
            // Ecx is energy transferred to neutrals
            Ecx(i,j,k) = (3./2)* (
                                       J_L * (Te_L - Tn_L)*R_cx_L
                                  + 4.*J_C * (Te_C - Tn_C)*R_cx_C
                                  +    J_R * (Te_R - Tn_R)*R_cx_R
                                  ) / (6. * J_C);
            
            // Fcx is friction between plasma and neutrals 
            Fcx(i,j,k) = (
                               J_L * (Vi_L - Vn_L)*R_cx_L
                          + 4.*J_C * (Vi_C - Vn_C)*R_cx_C
                          +    J_R * (Vi_R - Vn_R)*R_cx_R
                          ) / (6. * J_C);
            
            ///////////////////////////////////////
            // Recombination
            
            BoutReal R_rc_L  = hydrogen.recombination(Ne_L*Nnorm, Te_L*Tnorm)*SQ(Ne_L) * Nnorm / Omega_ci;
            BoutReal R_rc_C  = hydrogen.recombination(Ne_C*Nnorm, Te_C*Tnorm)*SQ(Ne_C) * Nnorm / Omega_ci;
            BoutReal R_rc_R  = hydrogen.recombination(Ne_R*Nnorm, Te_R*Tnorm)*SQ(Ne_R) * Nnorm / Omega_ci;
            
            // Rrec is radiated energy, Erec is energy transferred to neutrals
            // Factor of 1.09 so that recombination becomes an energy source at 5.25eV
            Rrec(i,j,k) = (
                                J_L * (1.09*Te_L - 13.6/Tnorm)*R_rc_L
                           + 4.*J_C * (1.09*Te_C - 13.6/Tnorm)*R_rc_C
                           +    J_R * (1.09*Te_R - 13.6/Tnorm)*R_rc_R
                           ) / (6. * J_C);
            
            Erec(i,j,k) = (3./2) * (
                                         J_L * Te_L * R_rc_L
                                    + 4.*J_C * Te_C * R_rc_C
                                    +    J_R * Te_R * R_rc_R
                                    ) / (6. * J_C);
	    
            Frec(i,j,k) = (
                                 J_L * Vi_L * R_rc_L
                           + 4.* J_C * Vi_C * R_rc_C
                           +     J_R * Vi_R * R_rc_R
                           ) / (6. * J_C);

            Srec(i,j,k) = (
                                 J_L * R_rc_L
                           + 4.* J_C * R_rc_C
                           +     J_R * R_rc_R
                           ) / (6. * J_C);
            
            ///////////////////////////////////////      
            // Ionisation
            BoutReal R_iz_L = Ne_L*Nn_L*hydrogen.ionisation(Te_L*Tnorm) * Nnorm / Omega_ci;
            BoutReal R_iz_C = Ne_C*Nn_C*hydrogen.ionisation(Te_C*Tnorm) * Nnorm / Omega_ci;
            BoutReal R_iz_R = Ne_R*Nn_R*hydrogen.ionisation(Te_R*Tnorm) * Nnorm / Omega_ci;
            
            Riz(i,j,k) = (Eionize/Tnorm) * (    // Energy loss per ionisation
                                                 J_L * R_iz_L
                                            + 4.*J_C * R_iz_C
                                            +    J_R * R_iz_R
                                             ) / (6. * J_C);   
            Eiz(i,j,k) = -(3./2)* (   // Energy from neutral atom temperature
                                          J_L * Tn_L * R_iz_L
                                   + 4. * J_C * Tn_C * R_iz_C
                                   +      J_R * Tn_R * R_iz_R
                                  ) / (6. * J_C);

            // Friction due to ionisation
            Fiz(i,j,k) = - (
                                   J_L * Vn_L * R_iz_L
                            + 4. * J_C * Vn_C * R_iz_C
                            +      J_R * Vn_R * R_iz_R
                            ) / (6. * J_C);
            
            // Plasma sink due to ionisation (negative)
            Siz(i,j,k) = - (
                                 J_L * R_iz_L
                           + 4.* J_C * R_iz_C
                           +     J_R * R_iz_R
                            ) / (6. * J_C);

            // Total energy lost from system
            R(i,j,k) = Rzrad(i,j,k)     // Radiated power from impurities
                     + Rrec(i,j,k)      // Recombination
                     + Riz(i,j,k);      // Ionisation
            
            // Total energy transferred to neutrals
            E(i,j,k) = Ecx(i,j,k)       // Charge exchange
                     + Erec(i,j,k)      // Recombination
                     + Eiz(i,j,k);      // ionisation

            // Total friction
            F(i,j,k) = Frec(i,j,k) + Fiz(i,j,k) + Fcx(i,j,k);

            // Total sink of plasma, source of neutrals
            S(i,j,k) = Srec(i,j,k) + Siz(i,j,k);
          }
      
      if(!evolve_nvn && neutral_f_pn) {
        // Not evolving neutral momentum
        F = Grad_par(Pn);
      }
    }
    

    ///////////////////////////////////////////////////
    // Plasma model

    /// Density
    
    if(rhs_explicit) {
      
      ddt(Ne) = 
        - Div_par_FV(Ne, Vi) // Mass flow
        - S                  // Sink to recombination
        ;
      
      if(volume_source) {
        ddt(Ne) += NeSource; // External volume source
      }
      
    }else {
      ddt(Ne) = 0.0;
    }

    if(rhs_implicit && (anomalous_D > 0.0)) {
      ddt(Ne) += Div_Par_Diffusion(anomalous_D, Ne);
    }
    
    if((hyper > 0.0) && (rhs_implicit)) {
      ddt(Ne) += D(Ne, hyper);
    }
    
    if((ADpar > 0.0) && (rhs_implicit)) {
      ddt(Ne) += ADpar * AddedDissipation(1.0, P, Ne, true);
    }
    
    /// Momentum

    if(rhs_explicit) {
      ddt(NVi) = 
        - Div_par_FV(NVi, Vi) // Momentum flow
        - Grad_par(P)         // Pressure gradient
        - F                   // Friction with neutrals
        ;
    }else {
      ddt(NVi) = 0.0; 
    }
    
    if((viscos > 0.) && (rhs_implicit)) {
      //ddt(NVi) += Div_Par_Diffusion(viscos*SQ(mesh->dy)*Ne, Vi);
      ddt(NVi) += Div_Par_Diffusion(viscos*SQ(mesh->dy), Vi);
    }
    
    if(rhs_implicit && (anomalous_D > 0.0)) {
      ddt(NVi) += Div_Par_Diffusion(anomalous_D*Vi, Ne);
    }
    
    if(ion_viscosity) {
       // Braginskii ion viscosity
      if(rhs_explicit) {
        // Update viscosity
        
        Field3D tau_i = sqrt(2 * mi_me) * tau_e;
        eta_i = (4./3) * 0.96 * Ne * tau_i * Te;  // Ti = Te
        eta_i.applyBoundary("neumann");
      }
      if(rhs_implicit) {
        ddt(NVi) += Div_Par_Diffusion(eta_i, Vi);
      }
    }
    
    if((hyper > 0.0) && (rhs_implicit)) {
      ddt(NVi) += D(NVi, hyper);
    }
    
    if((ADpar > 0.0) && (rhs_implicit)) {
      ddt(NVi) += ADpar * AddedDissipation(1.0, P, NVi, true);
    }

    if(rhs_explicit) {
      ddt(P) += // ddt(P) set earlier for sheath
        - Div_par_FV(P, Vi)         // Advection
        - (2./3)*P*Div_par(Vi)      // Compression
        - (2./3) * R          // Radiated power
        - (2./3) * E          // Energy transferred to neutrals
        ;
      
      if(volume_source) {
        // Volumetric source
      
        ddt(P) += PeSource;            // External source of energy
      }else {
        // Insert power into the first grid point
        for(RangeIterator r=mesh->iterateBndryLowerY(); !r.isDone(); r++)
          for(int jz=0; jz<mesh->ngz-1; jz++) {
            ddt(P)(r.ind, mesh->ystart, jz) += (2./3)*powerflux / ( mesh->dy(r.ind, mesh->ystart) * sqrt(mesh->g_22(r.ind, mesh->ystart)));
          }
      }
    }
    
    if(rhs_implicit) {
      ddt(P) += (2./3)*Div_Par_Diffusion(kappa_epar, Te);  // Parallel heat conduction

      if(anomalous_D > 0.0) {
        ddt(P) += Div_Par_Diffusion(anomalous_D*2.*Te, Ne);
      }
      if(anomalous_chi > 0.0) {
        ddt(P) += Div_Par_Diffusion(anomalous_chi, Te);
      }
    }
    
    if((hyper > 0.0) && (rhs_implicit)) {
      ddt(P) += D(P, hyper);
    }
    
    if((ADpar > 0.0) && (rhs_implicit)) {
      ddt(P) += ADpar * AddedDissipation(1.0, P, P, true);
    }

    // Switch off evolution at very low densities
    for(int i=0;i<mesh->ngx;i++)
      for(int j=mesh->ystart;j<=mesh->yend;j++)
        for(int k=0;k<mesh->ngz-1;k++) {
          if((Ne(i,j,k) < 1e-5) && (ddt(Ne)(i,j,k) < 0.0)) {
            ddt(Ne)(i,j,k) = 0.0;
            ddt(NVi)(i,j,k) = 0.0;
            ddt(P)(i,j,k) = 0.0;
          }
        }

    ///////////////////////////////////////////////////
    // Neutrals model
    //
    // 

    Field3D logPn = log(Pn);
    logPn.applyBoundary("neumann");
    
    if(rhs_explicit) {
      ddt(Nn) = 
        - Div_par_FV(Nn, Vn)        // Advection
        + S                         // Source from recombining plasma
        - nloss*Nn                  // Loss of neutrals from the system
        ;
    }else {
      ddt(Nn) = 0.0;
    }

    if(rhs_implicit) {
      if(dneut > 0.0) {
        ddt(Nn) += Div_Par_Diffusion(Dn * Nn, logPn); // Diffusion
      }
    }
    
    if((hyper > 0.0) && (rhs_implicit)) {
      ddt(Nn) += D(Nn, hyper);
    }

    if(evolve_nvn) {
      // Evolving momentum of the neutral gas

      if(rhs_explicit) {
        ddt(NVn) = 
          - Div_par_FV(NVn, Vn)        // Momentum flow
          + F                          // Friction with plasma
          - nloss*NVn                  // Loss of neutrals from the system
          - Grad_par(Pn)               // Pressure gradient
          ;
      }else {
        ddt(NVn) = 0.0;
      }
      
      if(rhs_implicit) {
        if(viscos > 0.) {
          // Note no factor of Nn
          ddt(NVn) += Div_Par_Diffusion(viscos*SQ(mesh->dy), Vn);
        }
        
        if(hyper > 0.) {
          // Numerical dissipation
          ddt(NVn) += D(NVn, hyper);
        }

        if(ion_viscosity) {
          // Relationship between heat conduction and viscosity for neutral gas
          // Chapman, Cowling "The Mathematical Theory of Non-Uniform Gases", CUP 1952
          // Ferziger, Kaper "Mathematical Theory of Transport Processes in Gases", 1972
          //
          Field3D eta_n = (2./5) * kappa_n;
          
          ddt(NVn) += Div_Par_Diffusion(eta_n, Vn);
        }
        
        if(dneut > 0.0)
          ddt(NVn) += Div_Par_Diffusion(NVn*Dn, logPn); // Diffusion
      }
    }

    if(evolve_pn) {
      // Evolving temperature of neutral gas
      // Essentially the same as the plasma equation

      if(rhs_explicit) {
        ddt(Pn) +=
          - Div_par_FV(Pn, Vn)                     // Advection
          - (2./3)*Pn*Div_par(Vn)                  // Compression
          + (2./3) * E                             // Energy transferred to neutrals
          - nloss * Pn                             // Loss of neutrals from the system
	  ;
      }

      if(rhs_implicit) {
        ddt(Pn) += (2./3)*Div_Par_Diffusion(kappa_n, Tn);  // Parallel heat conduction
        
        if(dneut > 0.0) {
          ddt(Pn) += Div_Par_Diffusion(Dn*Pn, logPn); // Perpendicular diffusion
        }
      }

      if((hyper > 0.0) && (rhs_implicit)) {
        ddt(Pn) += D(Pn, hyper);
      }
      
      // Switch off evolution at very low densities
      // This seems to be necessary to get through initial transients
      
      for(int i=0;i<mesh->ngx;i++)
        for(int j=mesh->ystart;j<=mesh->yend;j++)
          for(int k=0;k<mesh->ngz-1;k++) {
            if(Nn(i,j,k) < 1e-5) {
              // Relax to the plasma temperature
              ddt(Pn)(i,j,k) = -1e-2*(Pn(i,j,k) - Te(i,j,k)*Nn(i,j,k));
            }
          }
    }

    if(rhs_explicit) {
      // Boundary condition on fluxes
      BoutReal nredist;
      for(RangeIterator r=mesh->iterateBndryUpperY(); !r.isDone(); r++) {
        int jz = 0; // Z index
	
        flux_ion = 0.0;
        BoutReal flux_neut = 0.0;
        
        for(int j = mesh->yend+1; j<mesh->ngy; j++) {
          flux_ion += ddt(Ne)(r.ind, j, jz) * mesh->J(r.ind,j) * mesh->dy(r.ind,j);        
          flux_neut += ddt(Nn)(r.ind, j, jz) * mesh->J(r.ind,j) * mesh->dy(r.ind,j);
          
          ddt(Ne)(r.ind, j, jz) = 0.0;
          ddt(Nn)(r.ind, j, jz) = 0.0;
        }
        
        // Make sure that mass is conserved
        
        // Total amount of neutral gas to be added
        BoutReal nadd = flux_ion*frecycle + flux_neut + gaspuff;
        
        // Neutral gas arriving at the target
        BoutReal ntarget = (1 - fredistribute) * nadd / ( mesh->J(r.ind,mesh->yend) * mesh->dy(r.ind,mesh->yend) );
        
        ddt(Nn)(r.ind, mesh->yend, jz) += ntarget;
        
        if(evolve_nvn) {
          // Set velocity of neutrals coming from the wall to a fraction of the 
          // Franck-Condon energy
          BoutReal Vneut = - vwall * sqrt(3.5/Tnorm);
          ddt(NVn)(r.ind, mesh->yend, jz) += ntarget* Vneut;
        }
        
        if(evolve_pn) {
          // Set temperature of the incoming neutrals to F-C
          ddt(Pn)(r.ind, mesh->yend, jz) += ntarget * (3.5/Tnorm);
        }
        
        // Re-distribute neutrals
        nredist = fredistribute * nadd;
        
        // Divide flux_ion by J so that the result in the output file has units of flux per m^2
        flux_ion /= mesh->J(mesh->xstart, mesh->yend+1);
      }
      // Now broadcast redistributed neutrals to other processors
      MPI_Comm ycomm = mesh->getYcomm(mesh->xstart); // MPI communicator
      int np; MPI_Comm_size(ycomm, &np); // Number of processors
      
      // Broadcast from final processor (presumably with target)
      // to all other processors
      MPI_Bcast(&nredist, 1, MPI_DOUBLE, np-1, ycomm);
      
      // Distribute along length
      for(int j=mesh->ystart;j<=mesh->yend;j++) {
        // Neutrals into this cell
        BoutReal ncell = nredist * redist_weight(mesh->xstart,j) / ( mesh->J(mesh->xstart,j) * mesh->dy(mesh->xstart,j) );
        
        ddt(Nn)(mesh->xstart, j, 0) += ncell;
        
        // No momentum
        
        if(evolve_pn) {
          // Set temperature of the incoming neutrals to F-C
          ddt(Pn)(mesh->xstart, j, 0) += ncell * (3.5/Tnorm);
        }
      }
    }
    
    return 0;
  }
  
  int precon(BoutReal t, BoutReal gamma, BoutReal delta) {
    // Preconditioner
    /*
    if(ion_viscosity) {
      output << "PRECON\n";
      if(evolve_nvn) {
        // Viscosity in neutrals
        
      }
    }
    */
    
    static InvertPar *inv = NULL;
    if(!inv) {
      // Initialise parallel inversion class
      inv = InvertPar::Create();
      inv->setCoefA(1.0);
    }
    // Set the coefficient in front of Grad2_par2
    inv->setCoefB(-(2./3)*gamma*kappa_epar);
    Field3D dT = ddt(P);
    dT.applyBoundary("neumann");
    ddt(P) = inv->solve(dT);

    if(evolve_pn) {
      // Neutral pressure
      inv->setCoefB(-(2./3)*gamma*kappa_n);
      Field3D dT = ddt(Pn);
      dT.applyBoundary("neumann");
      ddt(Pn) = inv->solve(dT);
    }

    return 0;
  }

  int convective(BoutReal t) {
    // Turn on explicit terms only
    rhs_explicit = true;
    rhs_implicit = false;
    update_coefficients = true;
    return rhs(t);
  }

  int diffusive(BoutReal t, bool linear) {
    // turn on implicit terms only
    rhs_explicit = false;
    rhs_implicit = true;
    update_coefficients = !linear; // Don't update coefficients in linear solve
    return rhs(t);
  }
  
private:
  // Normalisation parameters
  BoutReal Tnorm, Nnorm, Bnorm, AA;
  BoutReal Cs0, Omega_ci, rho_s0, tau_e0, mi_me;

  /////////////////////////////////////////////////////////////////
  // Evolving quantities
  Field3D Ne, NVi, P; // Plasma (electron) density, momentum, and pressure
  Field3D Nn, NVn, Pn; // Neutral density, momentum, pressure
  
  Field3D Vi, Vn;  // Ion and neutral velocities
  
  bool evolve_nvn; // Evolve neutral momentum?
  bool evolve_pn;  // Evolve neutral pressure?
  
  /////////////////////////////////////////////////////////////////
  // Diffusion and viscosity coefficients
  
  Field3D Dn;           // Neutral gas diffusion
  BoutReal dneut;       // Neutral gas diffusion multiplier
  
  Field3D kappa_n;      // Neutral gas thermal conduction
  Field3D kappa_epar;   // Plasma thermal conduction

  Field3D tau_e;
  Field3D eta_i;        // Braginskii ion viscosity
  bool ion_viscosity;   // Braginskii ion viscosity on/off

  BoutReal nloss;       // Neutral loss rate (1/timescale)

  BoutReal anomalous_D, anomalous_chi; // Anomalous transport
  
  /////////////////////////////////////////////////////////////////
  // Atomic physics transfer channels
  
  bool atomic;     // Include atomic physics? 
  
  Field3D Srec, Siz;        // Plasma particle sinks due to recombination and ionisation
  Field3D Frec, Fiz, Fcx;   // Plasma momentum sinks due to recombination, ionisation, charge exchange and total
  Field3D Rrec, Riz, Rzrad; // Plasma power sinks due to recombination, ionisation, impurity radiation, charge exchange and total
  Field3D Erec, Eiz, Ecx;   // Transfer of power from plasma to neutrals
  
  Field3D S, F, E; // Exchange of particles, momentum and energy from plasma to neutrals
  Field3D R;       // Radiated power
  
  RadiatedPower *rad;            // Impurity atomic rates
  UpdatedRadiatedPower hydrogen; // Atomic rates
  
  BoutReal fimp;     // Impurity fraction (of Ne)
  BoutReal Eionize;  // Ionisation energy loss
  
  bool neutral_f_pn; // When not evolving NVn, use F = Grad_par(Pn)

  ///////////////////////////////////////////////////////////////
  // Sheath boundary
  
  BoutReal sheath_gamma;   // Sheath heat transmission factor
  BoutReal neutral_gamma;  // Neutral heat transmission
  
  bool sheath_density_drop; // Does density fall at the sheath?
                            // True->Flux constant. False->Free boundary

  BoutReal frecycle; // Recycling fraction
  BoutReal gaspuff;  // Additional source of neutral gas at the target plate
  BoutReal vwall;    // Velocity of neutrals coming from the wall
                     // as fraction of Franck-Condon energy

  BoutReal flux_ion; // Flux of ions to target (output)

  // Re-distribution of recycled neutrals
  Field2D redist_weight; // Weighting used to decide redistribution
  BoutReal fredistribute; // Fraction of recycled neutrals re-distributed along length
  
  ///////////////////////////////////////////////////////////////
  // Sources
  
  bool volume_source;         // Include volume sources?
  Field2D NeSource, PeSource; // Volume sources
  Field2D NeSource0; // Used in feedback control
  BoutReal powerflux; // Used if no volume sources
  
  // Upstream density controller
  BoutReal density_upstream;
  BoutReal density_controller_p, density_controller_i; // Controller settings
  
  BoutReal density_error_lasttime, density_error_last; // Value and time of last error
  BoutReal density_error_integral; // Integral of error
  
  ///////////////////////////////////////////////////////////////
  // Numerical dissipation
  
  BoutReal tn_floor; // Minimum neutral gas temperature [eV]
  
  BoutReal hyper, viscos;     // Numerical dissipation terms
  BoutReal ADpar;             // Added Dissipation numerical term
  
  Field2D dy4;
  
  // Numerical diffusion
  const Field3D D(const Field3D &f, BoutReal d) {
    if(d < 0.0)
      return 0.0;
    return Div_Par_Diffusion(d*SQ(mesh->dy), f);
    //return -D4DY4_FV(d*dy4,f);
  }
  
  ///////////////////////////////////////////////////////////////
  // Splitting into implicit and explicit
  bool rhs_implicit, rhs_explicit;
  bool update_coefficients;
};

BOUTMAIN(SD1D);
