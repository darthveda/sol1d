#include <mpi.h>

#include "div_ops.hxx"

#include <bout/mesh.hxx>
#include <globals.hxx>
#include <derivs.hxx>
#include <output.hxx>
#include <utils.hxx>
#include <bout/assert.hxx>

#include <cmath>

// Div ( n * a * b x Grad(f) / B ), where n is the quantity being advected
const Field3D Div_n_a_bxGrad_f_B(const Field3D &n, const Field3D &a, const Field3D &f, bool xflux, bool yderiv) {
  Field3D result = 0.0;

  // X flux
 
  Field3D Vx = a * mesh->J*DDZ(f,DIFF_C2);
  if(yderiv)
    Vx -= a * mesh->J*(mesh->g_23/mesh->g_22)*DDY(f,DIFF_C2);
  
  mesh->communicate(Vx);
  if(xflux) {
    Vx.applyBoundary("neumann");
  }else 
    Vx.applyBoundary("dirichlet");
  
  Field3D ns = n;
  Field3D fs = f;
  if(mesh->ShiftXderivs) {
    Vx = Vx.shiftZ(true); // Shift to orthogonal coordinates for X derivatives
    ns = n.shiftZ(true);
    fs = f.shiftZ(true);
  }
  
  int xs = mesh->xstart;
  int xe = mesh->xend;
  if(!xflux) {
    // Need to use 1st order on boundaries
    if(mesh->firstX()) {
      int jx = xs;
      for(int jy=mesh->ystart;jy<=mesh->yend;jy++)
        for(int jz=0;jz<mesh->ngz-1;jz++) {
          result(jx,jy,jz) += ( ns(jx+1,jy,jz)*Vx(jx+1,jy,jz) + ns(jx,jy,jz)*Vx(jx,jy,jz) ) / (2.*mesh->dx(jx,jy)*mesh->J(jx,jy));
          /*
          output.write("%d, %d, %d: += %e*%e + %e*%e -> %e\n", 
                       jx, jy, jz,
                       ns(jx+1,jy,jz), Vx(jx+1,jy,jz), 
                       ns(jx,jy,jz), Vx(jx,jy,jz),
                       result(jx,jy,jz));
          */
      }
      xs++;
    }
    if(mesh->lastX()) {
      int jx = xe;
      for(int jy=mesh->ystart;jy<=mesh->yend;jy++) 
        for(int jz=0;jz<mesh->ngz-1;jz++) {
          result(jx,jy,jz) -= ( ns(jx,jy,jz)*Vx(jx,jy,jz) + ns(jx-1,jy,jz)*Vx(jx-1,jy,jz) ) / (2.*mesh->dx(jx,jy)*mesh->J(jx,jy));

 //       if(jy == mesh->ystart)
 //         output.write("%d: + %e  - %e\n", jx, n(jx,jy,jz)*Vx(jx,jy,jz), n(jx-1,jy,jz)*Vx(jx-1,jy,jz));
        }
       xe--;
     }
  }

  //output << "RESULT X1: " << result(2,2,0) << endl;

  for(int jx=xs; jx<=xe;jx++) {
    for(int jy=mesh->ystart;jy<=mesh->yend;jy++) {
      for(int jz=0;jz<mesh->ngz-1;jz++) {
        int jxm = jx - 1;
        if(mesh->firstX() && (jx == mesh->xstart))
          jxm = jx;
        int jxp = jx + 1;
        if(mesh->lastX() && (jx == mesh->xend))
          jxp = jx;
        
        result(jx, jy, jz) += ( ns(jxp,jy,jz)*Vx(jxp,jy,jz) - ns(jxm,jy,jz)*Vx(jxm,jy,jz) ) / (2.*mesh->dx(jx, jy) * mesh->J(jx,jy));
  
     //   if(jy == mesh->ystart)
     //     output.write("%d: + %e  - %e\n", jx, n(jx+1,jy,jz)*Vx(jx+1,jy,jz), n(jx-1,jy,jz)*Vx(jx-1,jy,jz) );     
        //result(jx, jy, jz) += ( Vx(jx,jy,jz)*(n(jx+1,jy,jz) - n(jx-1,jy,jz)) + n(jx,jy,jz)*(Vx(jx+1,jy,jz) - Vx(jx-1,jy,jz)) ) / (2.*mesh->dx(jx, jy) * mesh->J(jx,jy));
        
      }
    }
  }

  //output << "RESULT X2: " << result(2,2,0) << endl;

  if(mesh->ShiftXderivs)
    result = result.shiftZ(false); // Back to field-aligned

  //output << "RESULT X3: " << result(2,2,0) << endl;

  Field3D dxf;
  dxf.allocate();
  for(int jx=mesh->xstart; jx<=mesh->xend;jx++) {
    for(int jy=mesh->ystart;jy<=mesh->yend;jy++) {
      for(int jz=0;jz<mesh->ngz-1;jz++) {
        int jxm = jx - 1;
        BoutReal fac = 2.0;
        if(mesh->firstX() && (jx == mesh->xstart)) {
          jxm = jx;
          //fac = 1.0;
        }
        int jxp = jx + 1;
        if(mesh->lastX() && (jx == mesh->xend)) {
          jxp = jx;
          //fac = 1.0;
        }
        
        dxf(jx,jy,jz) = (fs(jxp,jy,jz) - fs(jxm,jy,jz))/(fac*mesh->dx(jx,jy));
      }
    }
  }
  if(mesh->ShiftXderivs)
    dxf = dxf.shiftZ(false);
  
  if(yderiv) {
    // Y flux
    Field3D Vy = a * (mesh->J*mesh->g_23/mesh->g_22)*dxf;
    //Vy.applyBoundary("neumann");
    mesh->communicate(Vy);
    
    for(int jx=mesh->xstart; jx<=mesh->xend;jx++) {
      for(int jy=mesh->ystart;jy<=mesh->yend;jy++) {
        for(int jz=0;jz<mesh->ngz-1;jz++) {
          
          //if(jx == mesh->xstart) {
          //  output.write("Y %d: + %e  - %e\n", jy, n(jx,jy+1,jz)*Vy(jx,jy+1,jz), n(jx,jy-1,jz)*Vy(jx,jy-1,jz));
          //}
          //result(jx, jy, jz) += ( Vy(jx,jy,jz)*(n(jx,jy+1,jz) - n(jx,jy-1,jz)) + n(jx,jy,jz)*(Vy(jx,jy+1,jz) - Vy(jx,jy-1,jz)) ) / (2.*mesh->dy(jx, jy) * mesh->J(jx,jy));
          
          result(jx, jy, jz) += ( n(jx,jy+1,jz)*Vy(jx,jy+1,jz) - n(jx,jy-1,jz)*Vy(jx,jy-1,jz) ) / (2.*mesh->dy(jx, jy) * mesh->J(jx,jy));
          
        }
      }
    } 
  }
  
  //output << "RESULT Y: " << result(2,2,0) << endl;

  // Z flux
  Field3D Vz = -a*mesh->J*dxf;
  
  int ncz = mesh->ngz-1;
  for(int jx=mesh->xstart; jx<=mesh->xend;jx++) {
    for(int jy=mesh->ystart;jy<=mesh->yend;jy++) {
      for(int jz=0;jz<ncz;jz++) {
        
        //output.write("OLD: (%d,%d,%d) : (%e, %e)\n", jx,jy,jz,Vx(jx,jy,jz), Vz(jx,jy,jz));

        int jzp = (jz + 1) % ncz;
        int jzm = (jz - 1 + ncz) % ncz;
        
        result(jx, jy, jz) += ( n(jx,jy,jzp)*Vz(jx,jy,jzp) - n(jx,jy,jzm)*Vz(jx,jy,jzm) ) / (2.*mesh->dz * mesh->J(jx,jy));
        
      }
    }
  }
  
  //output << "RESULT Z: " << result(2,2,0) << endl;

  return result;
}


// Div ( n * b x Grad(f) / B ), where n is the quantity being advected
// This version uses central differences
const Field3D Div_n_bxGrad_f_B(const Field3D &n, const Field3D &f) {
  Field3D result = 0.0;
  
  // X flux
  
  Field3D g = n * mesh->J*(mesh->g_23/mesh->g_22);
  Field3D Fx = 0.5* (
                     //n * a * mesh->J*(DDZ(f) - (mesh->g_23/mesh->g_22)*DDY(f,DIFF_C2)) +
                     - g * DDY(f,DIFF_C2)
                     + f * DDY(g,DIFF_C2)
                     );

  Fx.applyBoundary("neumann");

  if(mesh->ShiftXderivs)
    Fx = Fx.shiftZ(true); // Shift to orthogonal coordinates for X derivatives
  
  for(int jx=mesh->xstart; jx<=mesh->xend;jx++) {
    for(int jy=mesh->ystart;jy<=mesh->yend;jy++) {
      for(int jz=0;jz<mesh->ngz-1;jz++) {
        
        result(jx, jy, jz) += ( Fx(jx+1,jy,jz) - Fx(jx-1,jy,jz) )
          / (2.*mesh->dx(jx, jy) * mesh->J(jx,jy));
      }
    }
  }
  if(mesh->ShiftXderivs)
    result = result.shiftZ(false); // Back to field-aligned
  
  // Y flux
  Field3D Fy;
  Fy.allocate();
  for(int jx = mesh->xstart; jx<=mesh->xend;jx++) {
    for(int jy=0;jy<mesh->ngy;jy++) {
      for(int jz=0;jz<mesh->ngz-1;jz++) {
        Fy(jx,jy,jz) = 0.5*(   g(jx,jy,jz)*(f(jx+1,jy,jz) - f(jx-1,jy,jz))
                             - f(jx,jy,jz)*(g(jx+1,jy,jz) - g(jx-1,jy,jz)) )
          /(2.*mesh->dx(jx,jy));
      }
    }
  }
  //mesh->communicate(Vy);
  //Fy.applyBoundary("neumann");
  
  for(int jx=mesh->xstart; jx<=mesh->xend;jx++) {
    for(int jy=mesh->ystart;jy<=mesh->yend;jy++) {
      for(int jz=0;jz<mesh->ngz-1;jz++) {
        result(jx, jy, jz) += ( Fy(jx,jy+1,jz) - Fy(jx,jy-1,jz) )
          / (2.*mesh->dy(jx, jy) * mesh->J(jx,jy));
      }
    }
  }
  /*
  // Z flux
  Field3D Vz = -DDX(f);
  
  int ncz = mesh->ngz-1;
  for(int jx=mesh->xstart; jx<=mesh->xend;jx++) {
    for(int jy=mesh->ystart;jy<=mesh->yend;jy++) {
      for(int jz=0;jz<ncz;jz++) {
        // Flux at lower boundary
        BoutReal fl = 0.5*(Fz(jx,jy,(jz-1+ncz)%ncz) + Fz(jx,jy,jz));
        // Flux at upper boundary
        BoutReal fu = 0.5*(Fz(jx,jy,jz) + Fz(jx,jy,(jz+1)%ncz));
        
        result(jx, jy, jz) += (fu - fl) / mesh->dz;
        
        //output.write("Val: %d,%d,%d : %e\n", jx,jy,jz, result(jx, jy, jz));
      }
    }
  }
  */
  return result;
}

// Div ( a * Grad_perp(f) )
const Field3D Div_Perp_Lap(const Field3D &a, const Field3D &fin) {
  Field3D result;
  result.allocate();

  Field3D f = fin;
  if(mesh->ShiftXderivs)
    f = fin.shiftZ(true); // Shift so X-Z is orthogonal
  
  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=mesh->ystart;j<=mesh->yend;j++)
      for(int k=0;k<mesh->ngz-1;k++) {
        // Z indices
        int kp = (k + 1) % (mesh->ngz-1);
        int km = (k - 1 + mesh->ngz-1) % (mesh->ngz-1);
        int kpp = (kp + 1) % (mesh->ngz-1);
        int kmm = (km - 1 + mesh->ngz-1) % (mesh->ngz-1);
        
        BoutReal fm, fp;
        
        // X flux
        
        // Calculate derivatives at (i-1,j,k)
        fm = mesh->g11(i-1,j)*(f(i,j,k) - f(i-2,j,k))/(2.*mesh->dx(i-1,j))
          //+ mesh->g12(i-1,j)*(f(i-1,j+1,k) - f(i-1,j-1,k))/(2.*mesh->dy(i-1,j))
          + mesh->g13(i-1,j)*(f(i-1,j,kp) - f(i-1,j,km))/(2.*mesh->dz);

        // Derivatives at (i,j,k)
        
        
        // Derivatives at (i+1,j,k)
        fp = mesh->g11(i+1,j)*(f(i+2,j,k) - f(i,j,k))/(2.*mesh->dx(i+1,j))
          //+ mesh->g12(i+1,j)*(f(i+1,j+1,k) - f(i+1,j-1,k))/(2.*mesh->dy(i-1,j))
          + mesh->g13(i+1,j)*(f(i+1,j,kp) - f(i+1,j,km))/(2.*mesh->dz);
          
        result(i,j,k) = (mesh->J(i+1,j)*fp*a(i+1,j,k) - mesh->J(i-1,j)*fm*a(i-1,j,k))/(2.*mesh->J(i,j)*mesh->dx(i,j));
        
        /*
        if((i == 90) && (j == 8) && (k == 0)) {
          output.write("Div: %e, %e, %e, %e, %e, %e\n",
                       f(i,j,k),
                       f(i,j,k) - f(i-2,j,k),
                       f(i+2,j,k) - f(i,j,k),
                       fm,
                       fp,
                       result(i,j,k));
        }
        */

        // Y flux
        
        // Derivatives at (i,j-1,k)
        fm = //mesh->g12(i,j-1)*(f(i+1,j-1,k) - f(i-1,j-1,k))/(2.*mesh->dx(i,j-1))
          + (mesh->g22(i,j-1) - 1./mesh->g_22(i,j-1))*(f(i,j+1,k) - f(i,j-1,k))/(2.*mesh->dy(i,j-1))
          + mesh->g23(i,j-1)*(f(i,j-1,kp) - f(i,j-1,km))/(2.*mesh->dz);
        
        
        // Derivatives at (i,j+1,k)
        fp = //mesh->g12(i,j+1)*(f(i+1,j-1,k) - f(i-1,j-1,k))/(2.*mesh->dx(i,j-1))
          + (mesh->g22(i,j+1) - 1./mesh->g_22(i,j+1))*(f(i,j+2,k) - f(i,j,k))/(2.*mesh->dy(i,j+1))
          + mesh->g23(i,j+1)*(f(i,j+1,kp) - f(i,j+1,km))/(2.*mesh->dz);
        
        //result(i,j,k) += (mesh->J(i,j+1)*fp*a(i,j+1,k) - mesh->J(i,j-1)*fm*a(i,j-1,k))/(mesh->J(i,j)*mesh->dy(i,j));
        
        // Z flux
        
        // Derivatives at (i,j,k-1)
        fm = mesh->g13(i,j)*(f(i+1,j,km) - f(i-1,j,km))/(2.*mesh->dx(i,j))
          + mesh->g23(i,j)*(f(i,j+1,km) - f(i,j-1,km))/(2.*mesh->dy(i,j))
          + mesh->g33(i,j)*(f(i,j,k) - f(i,j,kmm))/(2.*mesh->dz);
        
        fp = mesh->g13(i,j)*(f(i+1,j,kp) - f(i-1,j,kp))/(2.*mesh->dx(i,j))
          + mesh->g23(i,j)*(f(i,j+1,kp) - f(i,j-1,kp))/(2.*mesh->dy(i,j))
          + mesh->g33(i,j)*(f(i,j,kpp) - f(i,j,k))/(2.*mesh->dz);
        
        //result(i,j,k) += (fp*a(i,j,k+1) - fm*a(i,j,k-1))/(2.*mesh->dz);
      }
  
  return result;
}

// Div ( a * Grad_perp(f) )
const Field3D Div_Perp_Lap_x3(const Field3D &a, const Field3D &f, bool xflux) {
  Field3D result;
  result.allocate();
  
  Field3D fs = f;
  Field3D as = a;
  if(mesh->ShiftXderivs) {
    fs = f.shiftZ(true);
    as = a.shiftZ(true);
  }
  
  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=mesh->ystart;j<=mesh->yend;j++)
      for(int k=0;k<mesh->ngz-1;k++) {
        // Z indices
        int kp = (k + 1) % (mesh->ngz-1);
        int km = (k - 1 + mesh->ngz-1) % (mesh->ngz-1);
        
        BoutReal fm, fp;
        
        // X flux
        
        // Calculate derivatives at (i-1/2,j,k)
        fm = 0.5*(mesh->g11(i-1,j) + mesh->g11(i,j))*(fs(i,j,k) - fs(i-1,j,k))/(0.5*(mesh->dx(i-1,j) + mesh->dx(i,j)))
          //+ mesh->g12(i-1,j)*(f(i-1,j+1,k) - f(i-1,j-1,k))/(2.*mesh->dy(i-1,j))
          + 0.5*(mesh->g13(i-1,j) + mesh->g13(i,j))*(fs(i-1,j,kp) - fs(i-1,j,km) + f(i,j,kp) - f(i,j,km))/(4.*mesh->dz);
        
        // Derivatives at (i+1/2,j,k)
        fp = 0.5*(mesh->g11(i,j) + mesh->g11(i+1,j))*(fs(i+1,j,k) - fs(i,j,k))/(0.5*(mesh->dx(i+1,j) + mesh->dx(i,j)))
          //+ mesh->g12(i+1,j)*(f(i+1,j+1,k) - f(i+1,j-1,k))/(2.*mesh->dy(i-1,j))
          + 0.5*(mesh->g13(i,j) + mesh->g13(i+1,j))*(fs(i+1,j,kp) - fs(i+1,j,km) + fs(i,j,kp) - fs(i,j,km))/(4.*mesh->dz);
        
        if(mesh->firstX() && (i == mesh->xstart) && (!xflux))
          fm = 0.0;
        if(mesh->lastX() && (i == mesh->xend) && (!xflux))
          fp = 0.0;

        result(i,j,k) = ( + fp*0.5*(mesh->J(i+1,j)*as(i+1,j,k) + mesh->J(i,j)*as(i,j,k)) 
                          - fm*0.5*(mesh->J(i-1,j)*as(i-1,j,k) + mesh->J(i,j)*as(i,j,k)) )
          /(mesh->J(i,j)*mesh->dx(i,j));
        /*

        // Y flux
        
        // Derivatives at (i,j-1,k)
        fm = //mesh->g12(i,j-1)*(f(i+1,j-1,k) - f(i-1,j-1,k))/(2.*mesh->dx(i,j-1))
          + (mesh->g22(i,j-1) - 1./mesh->g_22(i,j-1))*(f(i,j+1,k) - f(i,j-1,k))/(2.*mesh->dy(i,j-1))
          + mesh->g23(i,j-1)*(f(i,j-1,kp) - f(i,j-1,km))/(2.*mesh->dz);
        
        
        // Derivatives at (i,j+1,k)
        fp = //mesh->g12(i,j+1)*(f(i+1,j-1,k) - f(i-1,j-1,k))/(2.*mesh->dx(i,j-1))
          + (mesh->g22(i,j+1) - 1./mesh->g_22(i,j+1))*(f(i,j+2,k) - f(i,j,k))/(2.*mesh->dy(i,j+1))
          + mesh->g23(i,j+1)*(f(i,j+1,kp) - f(i,j+1,km))/(2.*mesh->dz);
        
        //result(i,j,k) += (mesh->J(i,j+1)*fp*a(i,j+1,k) - mesh->J(i,j-1)*fm*a(i,j-1,k))/(mesh->J(i,j)*mesh->dy(i,j));
        */
        // Z flux
        
        // Derivatives at (i,j,k-1)
        fm = mesh->g13(i,j)*(fs(i+1,j,km) - fs(i-1,j,km))/(2.*mesh->dx(i,j))
          //+ mesh->g23(i,j)*(f(i,j+1,km) - f(i,j-1,km))/(2.*mesh->dy(i,j))
          + mesh->g33(i,j)*(fs(i,j,k) - fs(i,j,km))/(2.*mesh->dz);
        
        fp = mesh->g13(i,j)*(fs(i+1,j,kp) - fs(i-1,j,kp))/(2.*mesh->dx(i,j))
          //+ mesh->g23(i,j)*(f(i,j+1,kp) - f(i,j-1,kp))/(2.*mesh->dy(i,j))
          + mesh->g33(i,j)*(fs(i,j,kp) - fs(i,j,k))/(2.*mesh->dz);
        
        result(i,j,k) += (fp*a(i,j,kp) - fm*a(i,j,km))/(2.*mesh->dz);
      }

  if(mesh->ShiftXderivs)
    result = result.shiftZ(false);
  
  return result;
}

// Div ( a * Grad_perp(f) )
const Field3D Div_Perp_Lap_XYZ(const Field3D &a, const Field3D &f, bool bndryflux) {
  Field3D result = 0.0;
  
  Field3D fs = f;
  Field3D as = a;
  if(mesh->ShiftXderivs) {
    fs = f.shiftZ(true);
    as = a.shiftZ(true);
  }

  // X flux
  
  int xfirst = mesh->xstart;
  int xlast = mesh->xend+1;
  if(!bndryflux && !mesh->periodicX) {
    if(mesh->firstX())
      ++xfirst;
    if(mesh->lastX())
      --xlast;
  }
  
  for(int i=xfirst;i<=xlast;i++)
    for(int j=mesh->ystart;j<=mesh->yend;j++)
      for(int k=0;k<mesh->ngz-1;k++) {
        // Z indices
        int kp = (k + 1) % (mesh->ngz-1);
        int km = (k - 1 + mesh->ngz-1) % (mesh->ngz-1);
        
        // X flux
        
        // Calculate derivatives at (i-1/2,j,k)
        BoutReal fm = 0.5*(mesh->g11(i-1,j) + mesh->g11(i,j))*(fs(i,j,k) - fs(i-1,j,k))/(0.5*(mesh->dx(i-1,j) + mesh->dx(i,j)))
          //+ mesh->g12(i-1,j)*(f(i-1,j+1,k) - f(i-1,j-1,k))/(2.*mesh->dy(i-1,j))
          //+ 0.5*(mesh->g13(i-1,j) + mesh->g13(i,j))*(fs(i-1,j,kp) - fs(i-1,j,km) + f(i,j,kp) - f(i,j,km))/(4.*mesh->dz)
          ;

        BoutReal flux = -fm * 0.5*(mesh->J(i-1,j) + mesh->J(i,j)) * 0.5*(as(i-1,j,k) + as(i,j,k));
        
        result(i,j,k) += flux / (mesh->J(i,j)*mesh->dx(i,j));
        result(i-1,j,k) -= flux  / (mesh->J(i-1,j)*mesh->dx(i-1,j));
      }

  if(mesh->ShiftXderivs)
    result = result.shiftZ(false);
  
  // Y flux

  for(int i=mesh->xstart;i<=mesh->xend;i++) {
    int yfirst = mesh->ystart;
    int ylast = mesh->yend+1;
    if(!bndryflux && !mesh->periodicY(i)) {
      if(mesh->firstY(i))
        ++yfirst;
      if(mesh->lastY(i))
        --ylast;
    }
    for(int j=yfirst;j<=ylast;j++)
      for(int k=0;k<mesh->ngz-1;k++) {
        // Z indices
        int kp = (k + 1) % (mesh->ngz-1);
        int km = (k - 1 + mesh->ngz-1) % (mesh->ngz-1);
                                                         
        // Derivatives at (i,j-1/2,k)
        BoutReal fm = //mesh->g12(i,j-1)*(f(i+1,j-1,k) - f(i-1,j-1,k))/(2.*mesh->dx(i,j-1))
          + (mesh->g22(i,j-1) + mesh->g22(i,j) - 1./mesh->g_22(i,j-1) - 1./mesh->g_22(i,j))*(f(i,j,k) - f(i,j-1,k))/(mesh->dy(i,j-1) + mesh->dy(i,j))
          + 0.5*(mesh->g23(i,j-1) + mesh->g23(i,j)) * (f(i,j-1,kp) +  f(i,j,kp) - f(i,j-1,km) - f(i,j,km))/(4.*mesh->dz);
        
        BoutReal flux = -fm * 0.25*(mesh->J(i,j-1) + mesh->J(i,j))*(a(i,j-1,k) + a(i,j,k));
        
        result(i,j,k) += flux / (mesh->J(i,j)*mesh->dy(i,j));
        result(i,j-1,k) -= flux / (mesh->J(i,j-1)*mesh->dy(i,j-1));
        //if(i == 30) {
        //  output.write("%d: %e, %e, %e, %e\n", j, fm, flux, flux / (mesh->J(i,j)*mesh->dy(i,j)), flux / (mesh->J(i,j-1)*mesh->dy(i,j-1)));
        //}
      }
  }
  // Z flux
  
  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=mesh->ystart;j<=mesh->yend;j++)
      for(int k=0;k<mesh->ngz-1;k++) {
        // Z indices
        int kp = (k + 1) % (mesh->ngz-1);
        int km = (k - 1 + mesh->ngz-1) % (mesh->ngz-1);
        
        // Derivatives at (i,j,k-1)
        BoutReal fm = //mesh->g13(i,j)*(fs(i+1,j,km) - fs(i-1,j,km))/(2.*mesh->dx(i,j))
          + mesh->g23(i,j)* ( (f(i,j+1,km) + f(i,j+1,k)) - (f(i,j-1,km) + f(i,j-1,k)) ) / (4.*mesh->dy(i,j))
          + mesh->g33(i,j)*(f(i,j,k) - f(i,j,km))/mesh->dz;

        BoutReal flux = -fm * 0.5*(a(i,j,km) + a(i,j,k)) / mesh->dz;
        
        result(i,j,k) += flux;
        result(i,j,km) -= flux;
      }
  
  return result;
}

const Field3D Div_Par_Diffusion(const Field3D &K, const Field3D &f, bool bndry_flux) {
  Field3D result;
  result = 0.0;
  
  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=mesh->ystart-1;j<=mesh->yend;j++)
      for(int k=0;k<mesh->ngz-1;k++) {
        // Calculate flux at upper surface
        
        if(!bndry_flux && !mesh->periodicY(i)) {
          if((j == mesh->yend) && mesh->lastY(i))
            continue;
        
          if((j == mesh->ystart-1) && mesh->firstY(i))
            continue;
        }
        BoutReal c = 0.5*(K(i,j,k) + K(i,j+1,k)); // K at the upper boundary
        BoutReal J = 0.5*(mesh->J(i,j) + mesh->J(i,j+1)); // Jacobian at boundary
        
        //if((i == mesh->xstart) && (j == mesh->ystart) && (k==0))
        //  output << "c = " << c << endl;
        
        BoutReal g_22 = 0.5*(mesh->g_22(i,j) + mesh->g_22(i,j+1));
        
        BoutReal gradient = 2.*(f(i,j+1,k) - f(i,j,k)) / (mesh->dy(i,j) + mesh->dy(i,j+1));
        
        BoutReal flux = c * J * gradient / g_22;
        
        result(i,j,k) += flux / (mesh->dy(i,j) * mesh->J(i,j));
        result(i,j+1,k) -= flux / (mesh->dy(i,j+1) * mesh->J(i,j+1));
      }
  return result;
}

const Field3D Div_Par_Diffusion_Index(const Field3D &f, bool bndry_flux) {
  Field3D result;
  result = 0.0;
  
  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=mesh->ystart-1;j<=mesh->yend;j++)
      for(int k=0;k<mesh->ngz-1;k++) {
        // Calculate flux at upper surface
        
        if(!bndry_flux && !mesh->periodicY(i)) {
          if((j == mesh->yend) && mesh->lastY(i))
            continue;
        
          if((j == mesh->ystart-1) && mesh->firstY(i))
            continue;
        }
        BoutReal J = 0.5*(mesh->J(i,j) + mesh->J(i,j+1)); // Jacobian at boundary
        
        BoutReal gradient = f(i,j+1,k) - f(i,j,k);
        
        BoutReal flux = J * gradient;
        
        result(i,j,k) += flux / mesh->J(i,j);
        result(i,j+1,k) -= flux / mesh->J(i,j+1);
      }
  return result;
}

/*
 * Added Dissipation scheme (related to Momentum Interpolation)
 *
 * This uses a 3rd-order derivative of the pressure as
 * a correction to the velocity. 
 *
 * This should appear in the form
 * 
 * df/dt = ... + AddedDissipation(N, P, f);
 */
const Field3D AddedDissipation(const Field3D &N, const Field3D &P, const Field3D f, bool bndry_flux) {
  Field3D result = 0.0;
  
  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=mesh->ystart-1;j<=mesh->yend;j++)
      for(int k=0;k<mesh->ngz-1;k++) {
        // Calculate flux at upper surface
        
        if(!bndry_flux && !mesh->periodicY(i)) {
          if((j >= mesh->yend-1) && mesh->lastY(i))
            continue;
          
          if((j <= mesh->ystart) && mesh->firstY(i))
            continue;
        }
        
        // At upper boundary
        
        BoutReal d = 0.5*(1./N(i,j,k) + 1./N(i,j+1,k));
        
        // Velocity 
        BoutReal v = - 0.25*d*( (P(i,j-1,k) + P(i,j+1,k) - 2.*P(i,j,k)) - (P(i,j,k) + P(i,j+2,k)-2.*P(i,j+1,k)) );
        
        // Variable being advected. Could use different interpolation?
        BoutReal var = 0.5*(f(i,j,k) + f(i,j+1,k));

        BoutReal flux = var * v * (mesh->J(i,j) + mesh->J(i,j+1)) / (sqrt(mesh->g_22(i,j))+ sqrt(mesh->g_22(i,j+1)));
        
        result(i,j,k) -= flux / (mesh->dy(i,j) * mesh->J(i,j));
        result(i,j+1,k) += flux / (mesh->dy(i,j+1) * mesh->J(i,j+1));
      }
  return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// XPPM methods

BoutReal BOUTMIN(const BoutReal &a, const BoutReal &b, const BoutReal &c, const BoutReal &d) {
  BoutReal r1 = (a < b) ? a : b;
  BoutReal r2 = (c < d) ? c : d;
  return (r1 < r2) ? r1 : r2;
}

struct Stencil1D {
  // Cell centre values
  BoutReal c, m, p, mm, pp;
  
  // Left and right cell face values
  BoutReal L, R;
};

// First order upwind for testing
void Upwind(Stencil1D &n, const BoutReal h) {
  n.L = n.R = n.c;
}

// Fromm method
void Fromm(Stencil1D &n, const BoutReal h) {
  n.L = n.c - 0.25*(n.p - n.m);
  n.R = n.c + 0.25*(n.p - n.m);
}

void XPPM(Stencil1D &n, const BoutReal h) {
  // 4th-order PPM interpolation in X
  
  const BoutReal C = 1.25; // Limiter parameter
  
  BoutReal h2 = h*h;
  
  n.R = (7./12)*(n.c + n.p) - (1./12)*(n.m + n.pp);
  n.L = (7./12)*(n.c + n.m) - (1./12)*(n.mm + n.p);
  
  // Apply limiters
  if( (n.c - n.R)*(n.p - n.R) > 0.0 ) {
    // Calculate approximations to second derivative
    
    BoutReal D2 = (3./h2)*(n.c - 2*n.R + n.p);
    BoutReal D2L = (1./h2)*(n.m - 2*n.c + n.p);
    BoutReal D2R = (1./h2)*(n.c - 2.*n.p + n.pp);
    
    BoutReal D2lim; // Value to be used in limiter
    
    // Check if they all have the same sign
    if( (D2*D2L > 0.0) && (D2*D2R > 0.0) ) {
      // Same sign
	    
      D2lim = SIGN(D2) * BOUTMIN( C*fabs(D2L), C*fabs(D2R), fabs(D2) );
    }else {
      // Different sign
      D2lim = 0.0;
    }
    
    n.R = 0.5*(n.c + n.p) - (h2/6)*D2lim;
  }
  
  if( (n.m - n.L)*(n.c - n.L) > 0.0 ) {
    // Calculate approximations to second derivative
    
    BoutReal D2 = (3./h2)*(n.m - 2*n.L + n.c);
    BoutReal D2L = (1./h2)*(n.mm - 2*n.m + n.c);
    BoutReal D2R = (1./h2)*(n.m - 2.*n.c + n.p);
    
    BoutReal D2lim; // Value to be used in limiter
    
    // Check if they all have the same sign
    if( (D2*D2L > 0.0) && (D2*D2R > 0.0) ) {
      // Same sign
      
      D2lim = SIGN(D2) * BOUTMIN( C*fabs(D2L), C*fabs(D2R), fabs(D2) );
    }else {
      // Different sign
      D2lim = 0.0;
    }
    
    n.L = 0.5*(n.m + n.c) - (h2/6)*D2lim;
  }
  
  if( ( (n.R - n.c)*(n.c - n.L) <= 0.0 ) || ( (n.m - n.c)*(n.c - n.p) <= 0.0 ) ) {
    // At a local maximum or minimum
    
    BoutReal D2 = (6./h2)*(n.L - 2.*n.c + n.R);
    
    if(fabs(D2) < 1e-10) {
      n.R = n.L = n.c;
    }else {
      BoutReal D2C = (1./h2)*(n.m - 2.*n.c + n.p);
      BoutReal D2L = (1./h2)*(n.mm - 2*n.m + n.c);
      BoutReal D2R = (1./h2)*(n.c - 2.*n.p + n.pp);
    
      BoutReal D2lim;
      // Check if they all have the same sign
      if( (D2*D2C > 0.0) && (D2*D2L > 0.0) && (D2*D2R > 0.0) ) {
	// Same sign
	
	D2lim = SIGN(D2) * BOUTMIN( C*fabs(D2L), C*fabs(D2R), C*fabs(D2C), fabs(D2) );
	n.R = n.c + (n.R - n.c)*D2lim / D2;
	n.L = n.c + (n.L - n.c)*D2lim / D2;
      }else {
	// Different signs
	n.R = n.L = n.c;
      }
    }
  }
}

// Div (n * b x Grad(f)/B)
const Field3D Div_n_bxGrad_f_B_XPPM(const Field3D &n_in, const Field3D &f_in, bool bndry_flux) {
  Field3D result = 0;
  
  //////////////////////////////////////////
  // X-Z advection.
  // 
  //             Z
  //             |
  // 
  //    fmp --- vU --- fpp
  //     |      nU      |
  //     |               |
  //    vL nL        nR vR    -> X
  //     |               |
  //     |      nD       |
  //    fmm --- vD --- fpm
  //

  Field3D n = n_in;  // Done in orthogonal X-Z coordinates
  Field3D f = f_in;

  if(mesh->ShiftXderivs) {
    n = n_in.shiftZ(true);  // Done in orthogonal X-Z coordinates
    f = f_in.shiftZ(true);
  }
  
  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=mesh->ystart;j<=mesh->yend;j++)
      for(int k=0;k<mesh->ngz-1;k++) {
	int kp = (k+1) % (mesh->ngz-1);
	int kpp = (kp+1) % (mesh->ngz-1);
	int km = (k-1+mesh->ngz-1) % (mesh->ngz-1);
	int kmm = (km-1+mesh->ngz-1) % (mesh->ngz-1);
	
	// 1) Interpolate stream function f onto corners fmp, fpp, fpm
	
	BoutReal fmm = 0.25*(f(i,j,k) + f(i-1,j,k) + f(i,j,km) + f(i-1,j,km));
	BoutReal fmp = 0.25*(f(i,j,k) + f(i,j,kp) + f(i-1,j,k) + f(i-1,j,kp)); // 2nd order accurate
	BoutReal fpp = 0.25*(f(i,j,k) + f(i,j,kp) + f(i+1,j,k) + f(i+1,j,kp));
	BoutReal fpm = 0.25*(f(i,j,k) + f(i+1,j,k) + f(i,j,km) + f(i+1,j,km));
	
	// 2) Calculate velocities on cell faces
	
	BoutReal vU = mesh->J(i,j)*(fmp - fpp)/mesh->dx(i,j); // -J*df/dx
	BoutReal vD = mesh->J(i,j)*(fmm - fpm)/mesh->dx(i,j); // -J*df/dx
	
	BoutReal vR = 0.5*(mesh->J(i,j)+mesh->J(i+1,j))*(fpp - fpm)/mesh->dz; // J*df/dz 
	BoutReal vL = 0.5*(mesh->J(i,j)+mesh->J(i-1,j))*(fmp - fmm)/mesh->dz; // J*df/dz 
	
        //output.write("NEW: (%d,%d,%d) : (%e/%e, %e/%e)\n", i,j,k,vL,vR, vU,vD);

	// 3) Calculate n on the cell faces. The sign of the
	//    velocity determines which side is used.
	
	// X direction
	Stencil1D s;
	s.c  = n(i,  j,k);
	s.m  = n(i-1,j,k);
	s.mm = n(i-2,j,k);
	s.p  = n(i+1,j,k);
	s.pp = n(i+2,j,k);
	
	//Upwind(s, mesh->dx(i,j));
	//XPPM(s, mesh->dx(i,j)); 
	Fromm(s, mesh->dx(i,j));
	
        // Right side
        if((i==mesh->xend) && (mesh->lastX())) {
          // At right boundary in X
          
          if(bndry_flux) {
            BoutReal flux;
            if(vR > 0.0) {
              // Flux to boundary
              flux = vR * s.R;
            }else {
              // Flux in from boundary
              flux = vR * 0.5*(n(i+1,j,k) + n(i,j,k));
            }
            result(i,j,k)   += flux / (mesh->dx(i,j) * mesh->J(i,j));
            result(i+1,j,k) -= flux / (mesh->dx(i+1,j) * mesh->J(i+1,j));
          }
        }else {
          // Not at a boundary
          if(vR > 0.0) {
            // Flux out into next cell
            BoutReal flux = vR * s.R;
            result(i,j,k)   += flux / (mesh->dx(i,j) * mesh->J(i,j));
            result(i+1,j,k) -= flux / (mesh->dx(i+1,j) * mesh->J(i+1,j));

            //if(i==mesh->xend)
            //  output.write("Setting flux (%d,%d) : %e\n", j,k,result(i+1,j,k));
          }
	}
        
        // Left side
          
	if((i==mesh->xstart) && (mesh->firstX())) {
          // At left boundary in X
          
          if(bndry_flux) {
            BoutReal flux;
            
            if(vL < 0.0) {
              // Flux to boundary
              flux = vL * s.L;
              
            }else {
              // Flux in from boundary
              flux = vL * 0.5*(n(i-1,j,k) + n(i,j,k));
            }
            result(i,j,k)   -= flux / (mesh->dx(i,j) * mesh->J(i,j));
            result(i-1,j,k) += flux / (mesh->dx(i-1,j) * mesh->J(i-1,j));
          }
        }else {
          // Not at a boundary
          
          if(vL < 0.0) {
            BoutReal flux = vL * s.L;
            result(i,j,k)   -= flux / (mesh->dx(i,j) * mesh->J(i,j));
            result(i-1,j,k) += flux / (mesh->dx(i-1,j) * mesh->J(i-1,j));
          }
        }

	/// NOTE: Need to communicate fluxes

	// Z direction
	s.m  = n(i,j,km);
	s.mm = n(i,j,kmm);
	s.p  = n(i,j,kp);
	s.pp = n(i,j,kpp);
	
	//Upwind(s, mesh->dz);
	//XPPM(s, mesh->dz);
	Fromm(s, mesh->dz);

	if(vU > 0.0) {
	  BoutReal flux = vU * s.R / (mesh->J(i,j)*mesh->dz);
	  result(i,j,k)   += flux;
	  result(i,j,kp)  -= flux;
	}
	if(vD < 0.0) {
	  BoutReal flux = vD * s.L / (mesh->J(i,j)*mesh->dz);
	  result(i,j,k)   -= flux;
	  result(i,j,km)  += flux;
	}
      }
  communicateFluxes(result);

  if(mesh->ShiftXderivs)
    result = result.shiftZ(false); // Shift back to field-aligned

#if 0
  
  //////////////////////////////////////////
  // X-Y advection.
  // 
  //             Y
  //             |
  // 
  //    fmp --- vU --- fpp
  //     |      nU      |
  //     |              |
  //    vL nL        nR vR    -> X
  //     |              |
  //     |      nD      |
  //    fmm --- vD --- fpm
  //

  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=mesh->ystart;j<=mesh->yend;j++)
      for(int k=0;k<mesh->ngz-1;k++) {
        
        // 1) Interpolate stream function f_in onto corners fmm, fmp, fpp, fpm
        // This is complicated because the corner point in X-Y is not communicated
        // and at an X-point it is shared with 8 cells, rather than 4
        // (being at the X-point itself)
        
        
      }
  
#endif
  
  return result;
}

/// Div ( n * v )  -- Magnetic drifts
const Field3D Div_f_v_XPPM(const Field3D &n, const Vector3D &v, bool bndry_flux) {

  Field3D result = 0;

  Field3D vx = v.x;
  Field3D vz = v.z;
  
  if(mesh->ShiftXderivs) {
    vx = v.x.shiftZ(true);
    vz = v.z.shiftZ(true);
  }
    
  // X-Z advection
  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=mesh->ystart;j<=mesh->yend;j++)
      for(int k=0;k<mesh->ngz-1;k++) {
        int kp = (k+1) % (mesh->ngz-1);
	int kpp = (kp+1) % (mesh->ngz-1);
	int km = (k-1+mesh->ngz-1) % (mesh->ngz-1);
	int kmm = (km-1+mesh->ngz-1) % (mesh->ngz-1);
        
        // Calculate velocities
        BoutReal vU = 0.5*(vz(i,j,kp) + vz(i,j,k))*mesh->J(i,j);
        BoutReal vD = 0.5*(vz(i,j,km) + vz(i,j,k))*mesh->J(i,j);
        BoutReal vL = 0.25*(vx(i-1,j,k) + vx(i,j,k))*(mesh->J(i-1,j) + mesh->J(i,j));
        BoutReal vR = 0.25*(vx(i+1,j,k) + vx(i,j,k))*(mesh->J(i+1,j) + mesh->J(i,j));

        // X direction
	Stencil1D s;
	s.c  = n(i,  j,k);
	s.m  = n(i-1,j,k);
	s.mm = n(i-2,j,k);
	s.p  = n(i+1,j,k);
	s.pp = n(i+2,j,k);
        
        //Upwind(s, mesh->dx(i,j));
	//XPPM(s, mesh->dx(i,j)); 
	Fromm(s, mesh->dx(i,j));
        
        if((i==mesh->xend) && (mesh->lastX())) {
          // At right boundary in X
          
          if(bndry_flux) {
            BoutReal flux;
            if(vR > 0.0) {
              // Flux to boundary
              flux = vR * s.R;
            }else {
              // Flux in from boundary
              flux = vR * 0.5*(n(i+1,j,k) + n(i,j,k));
            }
            result(i,j,k)   += flux / (mesh->dx(i,j) * mesh->J(i,j));
            result(i+1,j,k) -= flux / (mesh->dx(i+1,j) * mesh->J(i+1,j));
          }
        }else {
          // Not at a boundary
          if(vR > 0.0) {
            // Flux out into next cell
            BoutReal flux = vR * s.R;
            result(i,j,k)   += flux / (mesh->dx(i,j) * mesh->J(i,j));
            result(i+1,j,k) -= flux / (mesh->dx(i+1,j) * mesh->J(i+1,j));
            //if(i==mesh->xend)
            //  output.write("Setting flux (%d,%d) : %e\n", j,k,result(i+1,j,k));
          }
	}
        
        // Left side
          
	if((i==mesh->xstart) && (mesh->firstX())) {
          // At left boundary in X
          
          if(bndry_flux) {
            BoutReal flux;
            
            if(vL < 0.0) {
              // Flux to boundary
              flux = vL * s.L;
              
            }else {
              // Flux in from boundary
              flux = vL * 0.5*(n(i-1,j,k) + n(i,j,k));
            }
            result(i,j,k)   -= flux / (mesh->dx(i,j) * mesh->J(i,j));
            result(i-1,j,k) += flux / (mesh->dx(i-1,j) * mesh->J(i-1,j));
          }
        }else {
          // Not at a boundary
          
          if(vL < 0.0) {
            BoutReal flux = vL * s.L;
            result(i,j,k)   -= flux / (mesh->dx(i,j) * mesh->J(i,j));
            result(i-1,j,k) += flux / (mesh->dx(i-1,j) * mesh->J(i-1,j));
          }
        }
        
        /// NOTE: Need to communicate fluxes

	// Z direction
	s.m  = n(i,j,km);
	s.mm = n(i,j,kmm);
	s.p  = n(i,j,kp);
	s.pp = n(i,j,kpp);
	
	//Upwind(s, mesh->dz);
	//XPPM(s, mesh->dz);
	Fromm(s, mesh->dz);

	if(vU > 0.0) {
	  BoutReal flux = vU * s.R / (mesh->J(i,j)*mesh->dz);
	  result(i,j,k)   += flux;
	  result(i,j,kp)  -= flux;
	}
	if(vD < 0.0) {
	  BoutReal flux = vD * s.L / (mesh->J(i,j)*mesh->dz);
	  result(i,j,k)   -= flux;
	  result(i,j,km)  += flux;
	}
        
      }
  communicateFluxes(result);
  if(mesh->ShiftXderivs) {
    result = result.shiftZ(false);  // Shift to field-aligned
  }
  
  return result;
}

void communicateFluxes(Field3D &f) {
  // Communicate fluxes between processors
  // Takes values in guard cells, and adds them to cells
  
  // Use X=0 as temporary buffer 
  if(mesh->xstart != 2)
    throw BoutException("communicateFluxes: Sorry!");
  
  int size = mesh->ngy * mesh->ngz;
  comm_handle xin, xout;
  if(!mesh->firstX())
    xin = mesh->irecvXIn(f[0][0], size, 0);
  
  if(!mesh->lastX())
    xout = mesh->irecvXOut(f[mesh->ngx-1][0], size, 1);
  
  // Send X=1 values
  if(!mesh->firstX())
    mesh->sendXIn(f[1][0], size, 1);
  
  if(!mesh->lastX())
    mesh->sendXOut(f[mesh->ngx-2][0], size, 0);
  
  // Wait
  if(!mesh->firstX()) {
    mesh->wait(xin);
    // Add to cells
    for(int y=mesh->ystart;y<=mesh->yend;y++) 
      for(int z=0;z<mesh->ngz-1;z++) {
        //output.write("Adding flux (%d,%d) : %e\n", y,z,f(0,y,z));
        f(2,y,z) += f(0,y,z);
      }
  }
  if(!mesh->lastX()) {
    mesh->wait(xout);
    // Add to cells
    for(int y=mesh->ystart;y<=mesh->yend;y++) 
      for(int z=0;z<mesh->ngz-1;z++) {
        f(mesh->ngx-3,y,z) += f(mesh->ngx-1,y,z);
      }
  }

  
}

const Field3D Div_Perp_Lap_FV_Index(const Field3D &a, const Field3D &f, bool xflux) {
  
  Field3D result = 0.0;

  //////////////////////////////////////////
  // X-Z diffusion in index space
  // 
  //            Z
  //            |
  // 
  //     o --- gU --- o
  //     |     nU     |
  //     |            |
  //    gL nL      nR gR    -> X
  //     |            |
  //     |     nD     |
  //     o --- gD --- o
  //

  
  Field3D fs = f;
  Field3D as = a;
  
  if(mesh->ShiftXderivs) {
    fs = f.shiftZ(true);
    as = a.shiftZ(true);
  }

  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=mesh->ystart;j<=mesh->yend;j++)
      for(int k=0;k<mesh->ngz-1;k++) {
        int kp = (k+1) % (mesh->ngz-1);
	int km = (k-1+mesh->ngz-1) % (mesh->ngz-1);
        
        // Calculate gradients on cell faces
        
        BoutReal gR =  fs(i+1,j,k) - fs(i,j,k);
        
        BoutReal gL = fs(i,j,k) - fs(i-1,j,k);
        
        BoutReal gD = fs(i,j,k) - fs(i,j,km);
        
        BoutReal gU = fs(i,j,kp) - fs(i,j,k);
        
        // Flow right
        BoutReal flux = gR * 0.25*(mesh->J(i+1,j) + mesh->J(i,j))*(mesh->dx(i+1,j) + mesh->dx(i,j))*(as(i+1,j,k) + as(i,j,k));
        
        result(i,j,k)   += flux / (mesh->dx(i,j)*mesh->J(i,j));
        
        // Flow left
        flux = gL * 0.25*(mesh->J(i-1,j) + mesh->J(i,j))*(mesh->dx(i-1,j) + mesh->dx(i,j)) *(as(i-1,j,k) + as(i,j,k));
        
        result(i,j,k)   -= flux / (mesh->dx(i,j)*mesh->J(i,j));
        
        // Flow up
          
        flux = gU * 0.5*(as(i,j,k) + as(i,j,kp));
        result(i,j,k) += flux;
        
        flux = gD * 0.5*(as(i,j,k) + as(i,j,km));
        result(i,j,k) -= flux;
      }
  
  if(mesh->ShiftXderivs) {
    result = result.shiftZ(false);
  }
  
  return result;
}

const Field3D Div_Perp_Lap_FV(const Field3D &a, const Field3D &f, bool xflux) {
  
  Field3D result = 0.0;

  //////////////////////////////////////////
  // X-Z diffusion
  // 
  //            Z
  //            |
  // 
  //     o --- gU --- o
  //     |     nU     |
  //     |            |
  //    gL nL      nR gR    -> X
  //     |            |
  //     |     nD     |
  //     o --- gD --- o
  //

  
  Field3D fs = f;
  Field3D as = a;
  
  if(mesh->ShiftXderivs) {
    fs = f.shiftZ(true);
    as = a.shiftZ(true);
  }

  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=mesh->ystart;j<=mesh->yend;j++)
      for(int k=0;k<mesh->ngz-1;k++) {
        int kp = (k+1) % (mesh->ngz-1);
	int km = (k-1+mesh->ngz-1) % (mesh->ngz-1);
        
        // Calculate gradients on cell faces
        
        BoutReal gR = (mesh->g11(i,j) + mesh->g11(i+1,j)) * (fs(i+1,j,k) - fs(i,j,k))/(mesh->dx(i+1,j) + mesh->dx(i,j))
          + 0.5*(mesh->g13(i,j) + mesh->g13(i+1,j))*(fs(i+1,j,kp) - fs(i+1,j,km) + fs(i,j,kp) - fs(i,j,km))/(4.*mesh->dz);
        
        BoutReal gL = (mesh->g11(i-1,j) + mesh->g11(i,j))*(fs(i,j,k) - fs(i-1,j,k))/(mesh->dx(i-1,j) + mesh->dx(i,j))
          + 0.5*(mesh->g13(i-1,j) + mesh->g13(i,j))*(fs(i-1,j,kp) - fs(i-1,j,km) + f(i,j,kp) - f(i,j,km))/(4.*mesh->dz);
        
        BoutReal gD = mesh->g13(i,j)*(fs(i+1,j,km) - fs(i-1,j,km) + fs(i+1,j,k) - fs(i-1,j,k))/(4.*mesh->dx(i,j))
          + mesh->g33(i,j)*(fs(i,j,k) - fs(i,j,km))/mesh->dz;
        
        BoutReal gU = mesh->g13(i,j)*(fs(i+1,j,kp) - fs(i-1,j,kp) + fs(i+1,j,k) - fs(i-1,j,k))/(4.*mesh->dx(i,j))
          + mesh->g33(i,j)*(fs(i,j,kp) - fs(i,j,k))/mesh->dz;
          
        
        // Flow right
        BoutReal flux = gR * 0.25*(mesh->J(i+1,j) + mesh->J(i,j)) *(as(i+1,j,k) + as(i,j,k));
        
        result(i,j,k)   += flux / (mesh->dx(i,j)*mesh->J(i,j));
        //result(i+1,j,k) -= flux / (mesh->dx(i+1,j)*mesh->J(i+1,j));
        
        // Flow left
        flux = gL * 0.25*(mesh->J(i-1,j) + mesh->J(i,j)) *(as(i-1,j,k) + as(i,j,k));
        
        result(i,j,k)   -= flux / (mesh->dx(i,j)*mesh->J(i,j));
        //result(i-1,j,k) += flux / (mesh->dx(i+1,j)*mesh->J(i+1,j));
        
        
        // Flow up
          
        flux = gU * 0.5*(as(i,j,k) + as(i,j,kp)) / mesh->dz;
        result(i,j,k) += flux;
        //result(i,j,kp) -= flux;
        
        flux = gD * 0.5*(as(i,j,k) + as(i,j,km)) / mesh->dz;
        result(i,j,k) -= flux;
        //result(i,j,km) += flux;
      }
  //communicateFluxes(result);

  if(mesh->ShiftXderivs) {
    result = result.shiftZ(false);
  }
  
  return result;
}


const Field3D Div_par_FV(const Field3D &f, const Field3D &v) {
  // Finite volume parallel divergence
  Field3D result = 0.0;
  
  bool highorder = false; // using a high order method

  int ys = mesh->ystart;
  int ye = mesh->yend;
  if(!highorder) {
    // Only need one guard cell, so no need to communicate fluxes
    // Instead calculate in guard cells
    
    //if(!mesh->firstY())
    ys--;
    //if(!mesh->lastY())
    ye++;
  }
  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=ys;j<=ye;j++) {
      for(int k=0;k<mesh->ngz-1;k++) {
        
        // Reconstruct f at the cell faces
        Stencil1D s;
        s.c  = f(i,j,  k);
        s.m  = f(i,j-1,k);
        s.p  = f(i,j+1,k);
        
        if(highorder) {
          s.mm = f(i,j-2,k);
          s.pp = f(i,j+2,k);
        }
        
        //Upwind(s, mesh->dy(i,j));  // 1st order accurate
        Fromm(s, mesh->dy(i,j));     // 2nd order, some upwinding
	//XPPM(s, mesh->dy(i,j));    // High order
        
        // Calculate velocity at right boundary (y+1/2)
        BoutReal vpar = 0.5*(v(i,j,k) + v(i,j+1,k));
        /*
        if((vpar < 0.0) && (j == mesh->yend)) {
          // In from right boundary
          BoutReal flux = 0.5*(f(i,j,k)+f(i,j+1,k)) * vpar * (mesh->J(i,j) + mesh->J(i,j+1)) / (sqrt(mesh->g_22(i,j)) + sqrt(mesh->g_22(i,j+1)));

          result(i,j,k)   += flux / (mesh->dy(i,j)*mesh->J(i,j));
          result(i,j+1,k) -= flux / (mesh->dy(i,j+1)*mesh->J(i,j+1));
        }else 
        */
        if(vpar > 0.0) {
          // Out of this cell; use s.R
	  
	  if(mesh->lastY(i) && (j == mesh->yend) && !mesh->periodicY(i)) {
	    // Last point in domain: Use mid-point to be consistent with boundary conditions
	    s.R = 0.5*(s.c + s.p);
	  }
	  BoutReal flux = s.R * vpar * (mesh->J(i,j) + mesh->J(i,j+1)) / (sqrt(mesh->g_22(i,j))+ sqrt(mesh->g_22(i,j+1)));
          
          result(i,j,k)   += flux / (mesh->dy(i,j)*mesh->J(i,j));
          result(i,j+1,k) -= flux / (mesh->dy(i,j+1)*mesh->J(i,j+1));
          
        }
        
        // Calculate at left boundary
        vpar = 0.5*(v(i,j,k) + v(i,j-1,k));
        
        if(vpar < 0.0) {
          // Out of this cell; use s.L
	  
          
	  if(mesh->firstY(i) && (j == mesh->ystart) && !mesh->periodicY(i)) {
	    // First point in domain
	    s.L = 0.5*(s.c + s.m);
	  }
          BoutReal flux = s.L * vpar * (mesh->J(i,j) + mesh->J(i,j-1)) / (sqrt(mesh->g_22(i,j)) + sqrt(mesh->g_22(i,j-1)));
          
          result(i,j,k)   -= flux / (mesh->dy(i,j)*mesh->J(i,j));
          result(i,j-1,k) += flux / (mesh->dy(i,j-1)*mesh->J(i,j-1));
        }
        
      }
      
    }

  if(highorder) {
    // Communicate fluxes.
    // Not implemented yet, so fail if NYPE > 1
    if(!(mesh->firstY() && mesh->lastY()))
      throw BoutException("Not implemented for NYPE > 1");
  }
  return result;
}

const Field3D D4DY4_FV(const Field3D &d, const Field3D &f) {
  Field3D result = 0.0;

  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=mesh->ystart;j<=mesh->yend;j++) {
      BoutReal dy3 = SQ(mesh->dy(i,j))*mesh->dy(i,j);
      for(int k=0;k<mesh->ngz-1;k++) {
        
        // 3rd derivative at right boundary
        
        BoutReal d3fdx3 = (
                                f(i,j+2,k)
                           - 3.*f(i,j+1,k)
                           + 3.*f(i,j,  k)
                           -    f(i,j-1,k)
                                ) / dy3;
                           
        BoutReal flux = 0.5*(d(i,j,k) + d(i,j+1,k))*(mesh->J(i,j) + mesh->J(i,j+1)) * d3fdx3;
        
        /*
        if(mesh->lastY() && (j == mesh->yend)) {
          // Boundary? Only if not periodic
          flux = 0.0;
        }
        */

        result(i,j,  k) += flux / (mesh->J(i,j) * mesh->dy(i,j));
        result(i,j+1,k) -= flux / (mesh->J(i,j+1) * mesh->dy(i,j+1));
        
        if(j == mesh->ystart && (!mesh->firstY())) {
          // Left cell boundary, no flux through boundaries
          d3fdx3 = (
                         f(i,j+1,k)
                    - 3.*f(i,j,  k)
                    + 3.*f(i,j-1,k)
                    -    f(i,j-2,k)
                    ) / dy3;
          
          flux = 0.5*(d(i,j,k) + d(i,j-1,k))*(mesh->J(i,j) + mesh->J(i,j-1)) * d3fdx3;
          
          result(i,j,  k) -= flux / (mesh->J(i,j) * mesh->dy(i,j));
          result(i,j-1,k) += flux / (mesh->J(i,j-1) * mesh->dy(i,j-1));
        }
      }
    }
  return result;
}

const Field3D Div_parV_FV(const Field3D &v) {
  // Finite volume parallel divergence of a flow velocity
  Field3D result = 0.0;
  
  bool highorder = false; // using a high order method

  int ys = mesh->ystart;
  int ye = mesh->yend;
  if(!highorder) {
    // Only need one guard cell, so no need to communicate fluxes
    // Instead calculate in guard cells
    
    //if(!mesh->firstY())
    ys--;
    //if(!mesh->lastY())
    ye++;
  }
  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=ys;j<=ye;j++) {
      for(int k=0;k<mesh->ngz-1;k++) {
        
        // Reconstruct v at the cell faces
        Stencil1D s;
        s.c  = v(i,j,  k);
        s.m  = v(i,j-1,k);
        s.p  = v(i,j+1,k);
        
        if(highorder) {
          s.mm = v(i,j-2,k);
          s.pp = v(i,j+2,k);
        }
        
        //Upwind(s, mesh->dy(i,j));  // 1st order accurate
        Fromm(s, mesh->dy(i,j));     // 2nd order, some upwinding
	//XPPM(s, mesh->dy(i,j));    // High order
        
        if(s.R > 0.0) {
          // Out of this cell; use s.R
	  
	  if(mesh->lastY(i) && (j == mesh->yend) && !mesh->periodicY(i)) {
	    // Last point in domain: Use 1st order to avoid overshoots at sheath
	    s.R = s.c;
	  }
	  BoutReal flux = s.R * (mesh->J(i,j) + mesh->J(i,j+1)) / (sqrt(mesh->g_22(i,j))+ sqrt(mesh->g_22(i,j+1)));
          
          result(i,j,k)   += flux / (mesh->dy(i,j)*mesh->J(i,j));
          result(i,j+1,k) -= flux / (mesh->dy(i,j+1)*mesh->J(i,j+1));
        }
        
        if(s.L < 0.0) {
          // Out of this cell; use s.L
	  
	  if(mesh->firstY() && (j == mesh->ystart)) {
	    // First point in domain
	    s.L = s.c;
	  }
          BoutReal flux = s.L * (mesh->J(i,j) + mesh->J(i,j-1)) / (sqrt(mesh->g_22(i,j)) + sqrt(mesh->g_22(i,j-1)));
          
          result(i,j,k)   -= flux / (mesh->dy(i,j)*mesh->J(i,j));
          result(i,j-1,k) += flux / (mesh->dy(i,j-1)*mesh->J(i,j-1));
        }
      }
      
    }

  if(highorder) {
    // Communicate fluxes.
    // Not implemented yet, so fail if NYPE > 1
    if(!(mesh->firstY() && mesh->lastY()))
      throw BoutException("Not implemented for NYPE > 1");
  }
  return result;
}

/*!
 * X-Y diffusion
 *
 * NOTE: Assumes g^12 = 0, so X and Y are orthogonal. Otherwise
 * we would need the corner cell values to take Y derivatives along X edges
 * 
 */
const Field2D Laplace_FV(const Field2D &k, const Field2D &f) {
  Field2D result;
  result.allocate();
  
  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=mesh->ystart;j<=mesh->yend;j++) {
      
      // Calculate gradients on cell faces
      
      BoutReal gR = (mesh->g11(i,j) + mesh->g11(i+1,j)) * (f(i+1,j) - f(i,j))/(mesh->dx(i+1,j) + mesh->dx(i,j)); 
      
      BoutReal gL = (mesh->g11(i-1,j) + mesh->g11(i,j))*(f(i,j) - f(i-1,j))/(mesh->dx(i-1,j) + mesh->dx(i,j));
      
      BoutReal gU = (mesh->g22(i,j) + mesh->g22(i,j+1)) * (f(i,j+1) - f(i,j))/(mesh->dx(i,j+1) + mesh->dx(i,j));
      
      BoutReal gD = (mesh->g22(i,j-1) + mesh->g22(i,j)) * (f(i,j) - f(i,j-1))/(mesh->dx(i,j) + mesh->dx(i,j-1)); 
      
      // Flow right
      
      BoutReal flux = gR * 0.25*(mesh->J(i+1,j) + mesh->J(i,j))*(k(i+1,j) + k(i,j));
      
      result(i,j)  = flux / (mesh->dx(i,j)*mesh->J(i,j));
      
      // Flow left
      
      flux = gL * 0.25*(mesh->J(i-1,j) + mesh->J(i,j))*(k(i-1,j) + k(i,j));
      result(i,j) -= flux / (mesh->dx(i,j)*mesh->J(i,j));
      
      // Flow up
      
      flux = gU * 0.25*(mesh->J(i,j+1) + mesh->J(i,j))*(k(i,j+1) + k(i,j));
      result(i,j) += flux / (mesh->dx(i,j)*mesh->J(i,j));
      
      // Flow down
      
      flux = gD * 0.25*(mesh->J(i,j-1) + mesh->J(i,j))*(k(i,j-1) + k(i,j));
      result(i,j) -= flux / (mesh->dx(i,j)*mesh->J(i,j));
    }
  return result;
}
