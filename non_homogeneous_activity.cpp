//
// Created by Birte Geerds in 2023-2026, based on Code from Abhinav Singh on 15.03.20.

// include all necessary libraries
#include "config.h"
#define BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
#define BOOST_MPL_LIMIT_VECTOR_SIZE 50
#include <iostream>
#include <fstream>
#include "DCPSE/DCPSE_op/DCPSE_op.hpp"
#include "DCPSE/DCPSE_op/DCPSE_Solver.hpp"
#include "Operators/Vector/vector_dist_operators.hpp"
#include "Vector/vector_dist_subset.hpp"
#include "DCPSE/DCPSE_op/EqnsStruct.hpp"
#include "OdeIntegrators/OdeIntegrators.hpp"
#include <vector>

//define namespace
using namespace std;

//define constexpr
constexpr int x = 0;
constexpr int y = 1;
constexpr int POLARIZATION= 0,VELOCITY = 1, VORTICITY = 2, NORMAL = 3,PRESSURE = 4, STRAIN_RATE = 5, STRESS = 6, MOLFIELD = 7, DPOL = 8, DV = 9, VRHS = 10, F1 = 11, F2 = 12, F3 = 13, F4 = 14, F5 = 15, F3_OLD = 16, V_T = 17, DIV = 18, DELMU = 19, HPB = 20, FE = 21, R = 22, PID = 23, POLD = 24, FESSPLAY = 25, FESBEND = 26, ZETA = 27;

//define parameters
double eta = 1;
double nu = 2.;
double gama = 0.5;
double lambda;
double Ks = 5;
double chi2 = -0.5;
double chi4 = -chi2;
double k = .0;
double dkK;
double Kb;
double sum_fe;
double sumsplay;
double sumbend;
double anglesum;
double anglesum_old;
double angle_weighted_sum;
double pol_at_boundary;
double boxsizehalf;
double delmu=-0.05714285714285714;
double zeta;

int wr_f;
int wr_at;
double V_err_eps;
//and timers
timer gt;//timer solving PDEs in one iteration before testing steady state condition
timer tt2;//
size_t Gd,vsz; //grid size, size of bulk

//define parameters
int flag = 1;
double steady_tol=7e-8;
double timetol;
double spacing;

//define vector names and respective types
void *vectorGlobal=nullptr,*vectorGlobal_bulk=nullptr,*vectorGlobal_boundary=nullptr,*SolverPointer=nullptr,*SolverPointerpetsc=nullptr;
const openfpm::vector<std::string>
PropNAMES={"00-Polarization","01-Velocity","02-Vorticity","03-Normal","04-Pressure","05-StrainRate","06-Stress","07-MolecularField","08-DPOL","09-DV","10-VRHS","11-f1","12-f2","13-f3","14-f4","15-f5","16-f3_old","17-V_T","18-DIV","19-DELMU","20-HPB","21-FrankEnergy","22-R","23-particleID","24-P_old", "25-FesSplay", "26-FesBend", "27-Zeta"};
typedef aggregate<VectorS<2, double>,VectorS<2, double>,double[2][2],VectorS<2, double>,double,double[2][2],double[2][2],VectorS<2, double>,VectorS<2, double>,VectorS<2, double>,VectorS<2, double>,double,double,double,double,double,double,VectorS<2, double>,double,double,double,double,double,int,VectorS<2, double>,double, double, double> Activegels;
typedef vector_dist_ws<2, double,Activegels> vector_type;
typedef vector_dist_subset<2, double, Activegels> vector_type2;


openfpm::vector<aggregate<vect_dist_key_dx[2]>> CorrVec;        // vector to store the Ids for the Dirichlet BC, 0: boundary 1: bulk

std::vector<double> avangle;
std::vector<double> fes;
std::vector<double> fesBend;
std::vector<double> fesSplay;
std::vector<double> times;
std::vector<double> polar;
std::ofstream AvAngle("avangle.txt", std::ofstream::app);
std::ofstream FesBend("fesbend.txt", std::ofstream::app);
std::ofstream FesSplay("fessplay.txt", std::ofstream::app);
std::ofstream Fes("fes.txt", std::ofstream::app);
std::ofstream Time("time.txt", std::ofstream::app);
std::ofstream Polar("polar.txt", std::ofstream::app);
std::ofstream Tracking("tracking.txt", std::ofstream::app);

//Functor to Compute RHS of the time derivative of the polarity
//solve PDEs and calculate velocity
template<typename DX,typename DY,typename DXX,typename DXY,typename DYY>
struct PolarEv
{
    DX &Dx;
    DY &Dy;
    DXX &Dxx;
    DXY &Dxy;
    DYY &Dyy;
    //Constructor
    PolarEv(DX &Dx,DY &Dy,DXX &Dxx,DXY &Dxy,DYY &Dyy):Dx(Dx),Dy(Dy),Dxx(Dxx),Dxy(Dxy),Dyy(Dyy)
    {}

    void operator()( const state_type_2d_ofp &X , state_type_2d_ofp &dxdt , const double t ) const
    {

        vector_type &Particles= *(vector_type *) vectorGlobal; //for all particles
        vector_type2 &Particles_bulk= *(vector_type2 *) vectorGlobal_bulk; //for bulk particles
        vector_type2 &Particles_boundary= *(vector_type2 *) vectorGlobal_boundary; //for boundary particles
        DCPSE_scheme<equations2d2, vector_type> &Solver= *(DCPSE_scheme<equations2d2, vector_type> *) SolverPointer; //define DCPS_scheme
        petsc_solver<double> &solverPetsc= *(petsc_solver<double> *) SolverPointerpetsc; //define petsc_solver

        //create cluster, bulk, boundary pointer
        auto & v_cl=create_vcluster();
        auto & bulk = Particles_bulk.getIds();
        auto & boundary = Particles_boundary.getIds();


        auto Dyx=Dxy; 
        //define parameters from fields with getV
        auto Pol=getV<POLARIZATION>(Particles);
        auto Pol_bulk=getV<POLARIZATION>(Particles_bulk);
        auto V = getV<VELOCITY>(Particles);
        auto Vbdry = getV<VELOCITY>(Particles_boundary);
        auto dVbdry = getV<DV>(Particles_boundary);
        auto h = getV<MOLFIELD>(Particles);
        auto u = getV<STRAIN_RATE>(Particles);
        auto dPol = getV<DPOL>(Particles);
        auto W = getV<VORTICITY>(Particles);
        auto r = getV<R>(Particles);
        auto dPol_bulk = getV<DPOL>(Particles_bulk);
        auto dPol_boundary = getV<DPOL>(Particles_boundary);

        auto sigma = getV<STRESS>(Particles);
        auto FranckEnergyDensity = getV<FE>(Particles);
        auto FranckEnergyDensitySplay = getV<FESSPLAY>(Particles);
        auto FranckEnergyDensityBend = getV<FESBEND>(Particles);
        auto f1 = getV<F1>(Particles);
        auto f2 = getV<F2>(Particles);
        auto f3 = getV<F3>(Particles);
        auto f4 = getV<F4>(Particles);
        auto f5 = getV<F5>(Particles);
        auto dV = getV<DV>(Particles);
        auto P = getV<PRESSURE>(Particles);
        auto P_bulk = getV<PRESSURE>(Particles_bulk);
        auto RHS = getV<VRHS>(Particles);
        auto RHS_bulk = getV<VRHS>(Particles_bulk);
        auto div = getV<DIV>(Particles);
        auto V_t = getV<V_T>(Particles);
        auto delmu = getV<DELMU>(Particles);
        auto zeta = getV<ZETA>(Particles);

        //impose Dirichlet BC for polarisation given the boundary and bulk pairs
        for (int j = 0; j < CorrVec.size(); j++) {
                auto p_out = CorrVec.get<0>(j)[0];
                auto p = CorrVec.get<0>(j)[1];
                Particles.getProp<0>(p)[x] = Particles.getProp<0>(p_out)[x];
                Particles.getProp<0>(p)[y] = Particles.getProp<0>(p_out)[y];
                Particles.getProp<F3>(p) = Particles.getProp<F3>(p_out);
            }

        //X data for the solver is polarization
        Pol[x]=X.data.get<0>();
        Pol[y]=X.data.get<1>();
        //ids for equations
        eq_id x_comp, y_comp;
        x_comp.setId(0);
        y_comp.setId(1);

        int n = 0,nmax = 500,errctr = 0, Vreset = 0; //namx...maximal loops aload; errctr...maximal loop for Verr>Verr_old; Vreset...1 if next loop allowed, 0 if not
        double V_err = 1, V_err_old,sum_velocity_diff, sum_velocity_square; //set parameters needed for pressure correction
        timer tt; //timer for initialization for each time step
        tt.start(); //start timer
        Particles.ghost_get<POLARIZATION>(SKIP_LABELLING); //update ghost particles for polarization
        //define reoccuring derivative operators
        texp_v<double> dxpx=Dx(Pol[x]),dxpy=Dx(Pol[y]),dypx=Dy(Pol[x]),dypy=Dy(Pol[y]),
                       dxxpx=Dxx(Pol[x]),dyypx=Dyy(Pol[x]),dxxpy=Dxx(Pol[y]),dyypy=Dyy(Pol[y]),
                       dxypx=Dxy(Pol[x]),dxypy=Dxy(Pol[y]);
        auto dyxpx=dxypx; //symmetry
        auto dyxpy=dxypy;

        //erickson stress from Ramaswamy paper
        sigma[x][x] = - Ks * dxpx * dxpx - Ks * dypy * dxpx
                      - Kb * dxpy * dxpy + Kb * dypx * dxpy
                      - k * dxpx;
        sigma[x][y] = 0.5*((- Ks * dypy * dxpy - Ks * dxpx * dxpy
                              - Kb * dypx * dxpx + Kb * dxpy * dxpx
                              - k * dxpy)
                        +(- Ks * dxpx * dypx - Ks * dypy * dypx
                              - Kb * dxpy * dypy + Kb * dypx * dypy
                              - k * dypx));
        sigma[y][x] = 0.5*((- Ks * dypy * dxpy - Ks * dxpx * dxpy
                              - Kb * dypx * dxpx + Kb * dxpy * dxpx
                              - k * dxpy)
                        +(- Ks * dxpx * dypx - Ks * dypy * dypx
                              - Kb * dxpy * dypy + Kb * dypx * dypy
                              - k * dypx));
        sigma[y][y] = - Ks * dypy * dypy - Ks * dxpx * dypy 
                      - Kb * dypx * dypx + Kb * dxpy * dypx
                      - k * dypy;


        Particles.ghost_get<STRESS>(SKIP_LABELLING); //update ghost particles for stress

        // calulate FranckEnergyDensity: bend, splay and total=bend+splay+softconstraint
        FranckEnergyDensityBend = (Kb/2.) * ((dypx * dypx)
                                          + (dxpy * dxpy)
                                          - 2. * dypx * dxpy);
        FranckEnergyDensitySplay = (Ks/2.) * ((dxpx * dxpx)
                                          + (dypy * dypy)
                                          + 2 * (dxpx) * (dypy));
        FranckEnergyDensity = FranckEnergyDensityBend
                              + FranckEnergyDensitySplay
                              + chi2/2. * (Pol[x] * Pol[x] + Pol[y] * Pol[y])
                              + chi4/4. * (Pol[x] * Pol[x] * Pol[x] * Pol[x] + Pol[y] * Pol[y] * Pol[y] * Pol[y] + 2 * Pol[x] * Pol[x] * Pol[y] * Pol[y])
                              + k * (dxpx + dypy)
                              ;
        Particles.ghost_get<FE, FESSPLAY, FESBEND>(SKIP_LABELLING); //update ghost particles for free energy density

        //equation for h_x and h_y (components of molecular field)
        f1 = -(chi2 * Pol[x] + chi4 * Pol[x] * Pol[x] * Pol[x] + chi4 * Pol[y] * Pol[y] * Pol[x])
            + Ks * (dxxpx + dyxpy)
            + Kb * (dyypx - dyxpy);
        f2 = -(chi2 * Pol[y] + chi4 * Pol[y] * Pol[y] * Pol[y] + chi4 * Pol[x] * Pol[x] * Pol[y])
            + Ks * (dyypy + dyxpx)
            + Kb * (dxxpy - dxypx);

        Particles.ghost_get<F1, F2>(SKIP_LABELLING); //update ghost particles for molecular field
        texp_v<double> Dxf1 = Dx(f1),Dxf2 = Dx(f2), Dyf1 = Dy(f1), Dyf2 = Dy(f2);

        //H_perpendicular
        f4 = f2 * Pol[x] - f1 * Pol[y];
        //H_parallel
        f5 = f1 * Pol[x] + f2 * Pol[y];
        Particles.ghost_get<F4, F5>(SKIP_LABELLING); //update ghost particles for molecular field

        // calulate RHS of Stokes equ (without pressure (because pressure correction will be done later)
        dV[x] = (- Dx(sigma[x][x]) - Dy(sigma[x][y]) //erickson stress
                - .5 * (f2 * dypx + Pol[x] * Dyf2 
                        - f1 * dypy - Pol[y] * Dyf1) // antisymmetric
                + zeta * Dx(delmu * Pol[x] * Pol[x]) + zeta * Dy(delmu * Pol[x] * Pol[y]) - zeta * Dx(delmu * ((Pol[x] * Pol[x] + Pol[y] * Pol[y])/2))       //active stress
                - nu/2. * (f1 * dxpx + Pol[x] * Dxf1
                            - f2 * dxpy - Pol[y] * Dxf2
                            + f2 * dypx + Pol[x] * Dyf2
                            + f1 * dypy + Pol[y] * Dyf1)); // flow alignment

        dV[y] = (- Dy(sigma[y][y]) - Dx(sigma[y][x]) //erickson stress
                - .5 * (f1 * dxpy + Pol[y] * Dxf1
                        - f2 * dxpx - Pol[x] * Dxf2)
                + zeta * Dy(delmu * Pol[y] * Pol[y]) + zeta * Dx(delmu * Pol[x] * Pol[y]) - zeta * Dy(delmu * ((Pol[x] * Pol[x] + Pol[y] * Pol[y])/2))      //active stress
                - nu/2. * (f2 * dypy + Pol[y] * Dyf2
                            - f1 * dypx - Pol[x] * Dyf1
                            + f1 * dxpy + Pol[y] * Dxf1
                            + f2 * dxpx + Pol[x] * Dxf2));

        tt.stop(); //stop timer for initialization
        if (v_cl.rank() == 0) {
        std::cout << "Init of Velocity took " << tt.getwct() << " seconds." << std::endl;
        }
        tt.start(); //start timer again
        V_err = 1;
        n = 0;
        errctr = 0;
        if (Vreset == 1) {
            P = 0; //consider no pressure
            Vreset = 0;
        }

        double divmax;
        // approximate velocity
        while ((V_err >= V_err_eps) && n <= nmax) {
            Particles.ghost_get<PRESSURE>(SKIP_LABELLING);
            //pressure correction: use divergence of velocity -> if not zero gets corrected by pressure to fulfill incompressibility
            RHS[x] = dV[x] + Dx(P);
            RHS[y] = dV[y] + Dy(P);
            Particles.ghost_get<VRHS>(SKIP_LABELLING);
            Solver.reset_b();
            Solver.impose_b(bulk, RHS[0], x_comp); //update b Seite von Gleichung
            Solver.impose_b(bulk, RHS[1], y_comp);
            Solver.impose_b(boundary, 0, x_comp);
            Solver.impose_b(boundary, 0, y_comp);
            Solver.reset_x_ig();
            Solver.impose_x_ig(bulk, V[0], x_comp); //update x_ig Seite von Gleichung
            Solver.impose_x_ig(bulk, V[1], y_comp);
            Solver.impose_x_ig(boundary, 0, x_comp);
            Solver.impose_x_ig(boundary, 0, y_comp);
            Solver.solve_with_solver_ig_successive(solverPetsc, V[x], V[y]);
            Vbdry=0; //no slip BC
            Particles.ghost_get<VELOCITY>(SKIP_LABELLING); //update ghost particles for velocity

            div = -(Dx(V[x]) + Dy(V[y])); //calculate divergence of velocity
            P_bulk = P + div;//pressure correction
            divmax = 0;

            for(int l=0; l<=bulk.size(); l++){
              auto p = bulk.get<0>(l);
              if(divmax<Particles.getProp<DIV>(p)){
                divmax = (Particles.getProp<DIV>(p));
              }

            }


            sum_velocity_diff = 0; //
            sum_velocity_square = 0;
            for (int j = 0; j < bulk.size(); j++) {
                auto p = bulk.get<0>(j);
                sum_velocity_diff += (Particles.getProp<V_T>(p)[0] - Particles.getProp<VELOCITY>(p)[0]) *
                       (Particles.getProp<V_T>(p)[0] - Particles.getProp<VELOCITY>(p)[0]) +
                       (Particles.getProp<V_T>(p)[1] - Particles.getProp<VELOCITY>(p)[1]) *
                       (Particles.getProp<V_T>(p)[1] - Particles.getProp<VELOCITY>(p)[1]);
                sum_velocity_square += Particles.getProp<VELOCITY>(p)[0] * Particles.getProp<VELOCITY>(p)[0] +
                        Particles.getProp<VELOCITY>(p)[1] * Particles.getProp<VELOCITY>(p)[1];
            }
            V_t = V;
            v_cl.sum(sum_velocity_diff);
            v_cl.sum(sum_velocity_square);
            v_cl.max(divmax);
            v_cl.execute();
            sum_velocity_diff = sqrt(sum_velocity_diff);
            sum_velocity_square = sqrt(sum_velocity_square);
            V_err = sum_velocity_diff / sum_velocity_square;
            //if error of velocity is higher then V_err_old set errctr +1, else =0
            if (V_err > V_err_old){
                errctr++;
            } else {
                errctr = 0;
            }

            V_err_old = V_err;//update V_err_old
            //f errctr >6: convergence loop broken and break
            if (n > 6) {
                if (errctr > 6) {
                    std::cout << "CONVERGENCE LOOP BROKEN DUE TO INCREASE/VERY SLOW DECREASE IN DIVERGENCE" << std::endl;
                    Vreset = 1;
                    break;
                  } else {
                    Vreset = 0;
                }
            }
            n++;
        }
        tt.stop();//stop timer
        if (v_cl.rank() == 0) {
            std::cout <<"Relative cgs error:"<<V_err<< ". Relative Divergence = " << divmax/sum_velocity_square/double(vsz) << " and took " << tt.getwct() << " seconds with " << n
                      << " iterations."
                      << std::endl;
        }

        Particles.ghost_get<VELOCITY>(SKIP_LABELLING);//update ghost for velocity
        //define reoccuring derivatives
        texp_v<double> dxvx=Dx(V[x]),dyvx=Dy(V[x]),dxvy=Dx(V[y]),dyvy=Dy(V[y]);
        //calculate strain rate
        u[x][x] = dxvx;
        u[x][y] = 0.5 * (dxvy + dyvx);
        u[y][x] = 0.5 * (dyvx + dxvy);
        u[y][y] = dyvy;

        // calculate vorticity
        W[x][x] = 0;
        W[x][y] = 0.5 * (dxvy - dyvx); 
        W[y][x] = 0.5 * (dyvx - dxvy); 
        W[y][y] = 0;

        //evolution equation for poilarization
        dPol_bulk[x] = f1/gama
                  - W[x][y] * Pol[y] 
                   - nu * u[x][x] * Pol[x] - nu * u[x][y] * Pol[y]
                   + 0.5*nu * u[x][x] * Pol[x] + 0.5*nu * u[y][y] * Pol[x] 
                   - (V[x] * dxpx + V[y] * dypx)
                   + lambda * delmu * Pol[x]
                   ;

        dPol_bulk[y] = f2/gama
                  - W[y][x] * Pol[x]
                   - nu * u[y][x] * Pol[x] - nu * u[y][y] * Pol[y]
                   + 0.5*nu * u[x][x] * Pol[y] + 0.5*nu * u[y][y] * Pol[y]
                   - (V[x] * dxpy + V[y] * dypy)
                   + lambda * delmu * Pol[y]
                   ;


        //impose Dirichlet BC with CorrVec
       for (int i = 0; i < CorrVec.size(); ++i)
               {
                   auto p_boundary = CorrVec.get<0>(i)[0];
                   auto p_bulk = CorrVec.get<0>(i)[1];

                   Particles.getProp<0>(p_bulk)[x] = Particles.getProp<0>(p_boundary)[x];
                   Particles.getProp<0>(p_bulk)[y] = Particles.getProp<0>(p_boundary)[y];
                   Particles.getProp<F3>(p_bulk) = Particles.getProp<F3>(p_boundary);
               }
        Particles.ghost_get<POLARIZATION, F3>(SKIP_LABELLING); //update ghost particles for polarization and angle
        //set dxdt for solver
        dxdt.data.get<0>()=dPol[x];
        dxdt.data.get<1>()=dPol[y];
    }
};

// Functor to find steady state
template<typename DX,typename DY,typename DXX,typename DXY,typename DYY>
struct ObserverFunctor
{

    DX &Dx;
    DY &Dy;
    DXX &Dxx;
    DXY &Dxy;
    DYY &Dyy;

    double t_old;//for calculating time for one step
    int ctr; //counter for wiritng (into txt or pvtp) files

    //Constructor
    ObserverFunctor(DX &Dx,DY &Dy,DXX &Dxx,DXY &Dxy,DYY &Dyy):Dx(Dx),Dy(Dy),Dxx(Dxx),Dxy(Dxy),Dyy(Dyy)
    {
        t_old = 0.0;
        ctr = 0;
    }

    struct MinInfo{
        int pb;
        double rank;
        double lowest_pol;
    };

    void operator() (state_type_2d_ofp &state, double t)
    {
        vector_type &Particles= *(vector_type *) vectorGlobal;
        vector_type2 &Particles_bulk= *(vector_type2 *) vectorGlobal_bulk;
        vector_type2 &Particles_boundary= *(vector_type2 *) vectorGlobal_boundary;

        //create cluster, and get ids
        auto & v_cl=create_vcluster();
        auto & bulk = Particles_bulk.getIds();
        auto & boundary = Particles_boundary.getIds();

        //define parameters from fields with getV
        auto Pol=getV<POLARIZATION>(Particles);
        auto Pol_bulk=getV<POLARIZATION>(Particles_bulk);
        auto dPol = getV<DPOL>(Particles);
        auto dPol_bulk = getV<DPOL>(Particles_bulk);
        auto dPol_boundary = getV<DPOL>(Particles_boundary);
        auto FranckEnergyDensity = getV<FE>(Particles);
        auto Pol_old = getV<POLD>(Particles);
        auto f5 = getV<F5>(Particles);
        auto f3 = getV<F3>(Particles);
        auto f3_old = getV<F3_OLD>(Particles);
        auto P = getV<PRESSURE>(Particles);

        P=0;//set pressure to zero to meet incompressibility condition

        double angle_weighted;//introduce a weighted anglesum to accoutn for numerical error
        gt.stop();//stop wall clock time for one step

        Point<2,double> Xpn,X2pn,xpd1={boxsizehalf,boxsizehalf};//define points and middle points for weighting 
        
        double x_def_local = boxsizehalf; // Initialize local x-coordinate
        double y_def_local = boxsizehalf; // Initialize local y-coordinate
        double x_def_global = boxsizehalf; // Initialize global x-coordinate
        double y_def_global = boxsizehalf; // Initialize global y-coordinate
        Point<2, double> xpdef;

        // get point for defect as the lowest value for polarity in bulk
        double lowest_pol_local = 1.0; // Initialize the local minimum with a high value
        double lowest_pol_global = 1.0; // Initialize the global minimum with the same high value

        MinInfo n_mins_local{-1, -1, std::numeric_limits<double>::max()};
        MinInfo n_mins_global{-1, -1, std::numeric_limits<double>::max()};

        for (int n = 0; n < bulk.size(); n++) {
            auto pb = bulk.get<0>(n);
            double p1 = Particles.getProp<POLARIZATION>(n)[x];
            double p2 = Particles.getProp<POLARIZATION>(n)[y];
            double mag_pol = sqrt(p1 * p1 + p2 * p2);

            if (mag_pol <= lowest_pol_local) {
                lowest_pol_local = mag_pol;
                n_mins_local = {static_cast<int>(pb), static_cast<double>(v_cl.rank()), static_cast<double>(mag_pol)};// Store local minimum, rank, and the third value
                x_def_local = Particles.getPos(n)[0];
                y_def_local = Particles.getPos(n)[1];
            }
        }

        //another try

        openfpm::vector<double> minpols;
        openfpm::vector<double> minranks;
        openfpm::vector<double> minxdefs;
        openfpm::vector<double> minydefs;

        MinInfo localMinInfo;
        localMinInfo.lowest_pol = lowest_pol_local;
        localMinInfo.rank = v_cl.rank();

        minpols.add(localMinInfo.lowest_pol);
        minranks.add(localMinInfo.rank);
        minxdefs.add(x_def_local);
        minydefs.add(y_def_local);

        v_cl.execute();
        v_cl.allGather(localMinInfo.lowest_pol, minpols);
        v_cl.allGather(localMinInfo.rank, minranks);
        v_cl.allGather(x_def_local, minxdefs);
        v_cl.allGather(y_def_local, minydefs);
        v_cl.execute();


        for (size_t i = 0; i < minpols.size(); ++i)
        {
            double gatheredMin = minpols.get(i);
            size_t processorID = minranks.get(i);
            double x_def = minxdefs.get(i);
            double y_def = minydefs.get(i);
            size_t min_rank_global;

            if(gatheredMin<=lowest_pol_global){
                lowest_pol_global= gatheredMin;
                x_def_global = x_def;
                y_def_global = y_def;
            }

            std::cout << "Processor " << v_cl.rank() << ": Local minimum from Processor " << processorID << ": " << gatheredMin << " at " << x_def << " , " << y_def << std::endl;
        }
        xpdef = {x_def_global, y_def_global};

        // Print the information on each processor
        std::cout << "Processor " << v_cl.rank() << ": Found local minimum " << n_mins_local.lowest_pol << " at particle index " << n_mins_local.pb << " on rank " << n_mins_local.rank << std::endl;
        
	    cout << "Global minimum polarization value: " << lowest_pol_global << '\n';
	    cout << "Global minimum coordinates: (" << x_def_global << ", " << y_def_global << ")\n";
		
        std::cout <<"L. 585" << std::endl;
        if (t != 0) {
            //update values for solver
            Pol[x]=state.data.get<0>();
            Pol[y]=state.data.get<1>();
            
            Particles.ghost_get<POLARIZATION>(SKIP_LABELLING); //update ghost particles for polarization

            

            //normalize polarization for outer region
            for (int s=0; s<bulk.size(); s++){
              auto p = bulk.get<0>(s);
              Xpn={Particles.getPos(p)[0],Particles.getPos(p)[1]};
              double dist=Xpn.distance(xpd1);
              dist = (dist == 0) ? 1 : dist;
              double dist_def =  Xpn.distance(xpdef);
              dist_def = (dist_def == 0) ? 1 : dist_def;

              double p1 = Particles.getProp<POLARIZATION>(s)[x];
              double p2 = Particles.getProp<POLARIZATION>(s)[y];
              //compute normal
              double x1 = (xpdef[0] - Xpn[0])/dist_def;
              double x2 = (xpdef[1] - Xpn[1])/dist_def;
              double scalar = x1 * p1 + x2 * p2;
              double pmag=sqrt(p1*p1 + p2*p2);
              double nmag = sqrt(x1*x1 + x2*x2);
              double norm = nmag * pmag;
              //calculate weighted angle as angle times fermi function
              double angle = acos(scalar/norm);
              Particles.getProp<F3>(s) = angle;
              angle_weighted = angle * (1- 1/(1 + exp(2 * (dist - 2))) -(1 - 1/(1 + exp(2 * (- 2)))));
            }
            Particles.ghost_get<POLARIZATION, F3>(SKIP_LABELLING);
            //impose Dirichlet BC
            for (int j = 0; j < CorrVec.size(); j++) {
                auto p_out = CorrVec.get<0>(j)[0];
                auto p = CorrVec.get<0>(j)[1];
                //redirecting the boundary to trap the defect in the plane
                Xpn={Particles.getPos(p)[0],Particles.getPos(p)[1]};
                double dist=Xpn.distance(xpd1);
                dist = (dist == 0) ? 1 : dist;

                double p1 = Particles.getProp<POLARIZATION>(p_out)[x];
                double p2 = Particles.getProp<POLARIZATION>(p_out)[y];

                //enables to fix boundary independendtly from bulk
                double x1 = (boxsizehalf - Xpn[0])/dist;
                double x2 = (boxsizehalf - Xpn[1])/dist;
                double scalar = x1 * p1 + x2 * p2;
                double pmag=sqrt(p1*p1 + p2*p2);

                //update boundary independently of bulk
                Particles.getProp<POLARIZATION>(p_out)[x]=Particles.getProp<POLARIZATION>(p_out)[x]/pmag;
                Particles.getProp<POLARIZATION>(p_out)[y]=Particles.getProp<POLARIZATION>(p_out)[y]/pmag;

                //impose Dirichlet BC
                Particles.getProp<0>(p)[x]=Particles.getProp<0>(p_out)[x];
                Particles.getProp<0>(p)[y]=Particles.getProp<0>(p_out)[y];
                Particles.getProp<F3>(p)=Particles.getProp<F3>(p_out);
            }
        }
        Particles.ghost_get<POLARIZATION, F3>(SKIP_LABELLING);
        
        //print information over simulation step
        if (v_cl.rank() == 0) {
            std::cout << "point of the defect is: x="<< xpdef[0] << " , y=" <<xpdef[1]<< "\n";
            std::cout << "dt for the stepper is " << t-t_old << " Time Taken: "<<gt.getwct()
                      << "\n"
                      << "Time t = " << t << "\n";
        }
        //write pvtp file
        if (ctr%wr_at==0 || ctr==wr_f){
            Particles.deleteGhost();
            Particles.write_frame("Polar",ctr,t);
            Particles.ghost_get<POLARIZATION>();
        }
        ctr++;//increment counter

        int ctr_normalization = 1;//introduce counter for normalization of sums

        //difference for polarization
        dPol[x]=Pol[x]-Pol_old[x];
        dPol[y]=Pol[y]-Pol_old[y];

        //introduce maxrate of change
        double MaxRateOfAngleChange =0;
        for (int j = 0; j < bulk.size(); j++) {
            //build grid
            auto p = bulk.get<0>(j);
            auto xp=Particles.getPos(p);

            //calculate magnitude
            double magnitude = sqrt(Particles.getProp<POLARIZATION>(j)[0] * Particles.getProp<POLARIZATION>(j)[0]
                                  + Particles.getProp<POLARIZATION>(j)[1] * Particles.getProp<POLARIZATION>(j)[1]);
            //only consider, if magnitude of polarization is >0.99, which is outer region, defined before
            if(magnitude>1e-11){
              sum_fe += Particles.getProp<FE>(j);
              sumsplay += Particles.getProp<FESSPLAY>(j);
              sumbend += Particles.getProp<FESBEND>(j);
              anglesum += Particles.getProp<F3>(j);
              ctr_normalization++;
                //get maximum change of angle in respective region
                if(fabs(Particles.getProp<F3>(p) - Particles.getProp<F3_OLD>(p))>MaxRateOfAngleChange){
                    MaxRateOfAngleChange = fabs(Particles.getProp<F3>(p) - Particles.getProp<F3_OLD>(p));
                }
            }
        }
        for (int k=0; k<boundary.size();k++){
            pol_at_boundary +=sqrt(Particles.getProp<POLARIZATION>(k)[0] * Particles.getProp<POLARIZATION>(k)[0]
                                  + Particles.getProp<POLARIZATION>(k)[1] * Particles.getProp<POLARIZATION>(k)[1]);
        }

        pol_at_boundary=pol_at_boundary/boundary.size();

        double normalise = v_cl.size(); //get total number of points included
        v_cl.sum(anglesum); //sum over all angles
        v_cl.execute();
        //normalize sums
        if(ctr_normalization!=0){
            sum_fe = sum_fe/ctr_normalization;
            sumsplay = sumsplay/ctr_normalization;
            sumbend = sumbend/ctr_normalization;
            anglesum = anglesum/(normalise*ctr_normalization);
        }
        v_cl.max(MaxRateOfAngleChange);//get maximum of angle change on cluster
        v_cl.execute();
        //print out values

        std::cout << "ctr = "<< ctr << '\n';
        std::cout << "v_cl_rank = " << v_cl.rank() << '\n';
        std::cout << "FEsum_new = " << sum_fe << '\n';
        std::cout << "angle = " << anglesum << '\n';
        std::cout << "pol_at_boundary= " << pol_at_boundary << '\n';
        std::string anglesumSt = std::to_string(anglesum);
        AvAngle << anglesumSt << std::endl;
        if(v_cl.rank()==0)
        std::cout << "MaxRateOfAngleChange: "<< MaxRateOfAngleChange << std::endl;
        
        //check for minimum of free energy
        bool minimum=false;
        for (int i=(fes.size()-40); i<fes.size();i++){
          if(sum_fe<=fes[i]){
            minimum=true;}
          else if (fabs(fes[i-1])<steady_tol){
            minimum=true;}
          else
            minimum = false;
        }

        //write to files
        if(v_cl.rank()==0){
        std::string sumbendSt = std::to_string(sumbend);
        FesBend << sumbendSt << std::endl;
        std::string sumsplaySt = std::to_string(sumsplay);
        FesSplay << sumsplaySt << std::endl;
        std::string sum2St = std::to_string(sum_fe);
        Fes << sum2St << std::endl;
        std::string timeSt = std::to_string(t);
        Time << timeSt << std::endl;
        std::string pol_at_boundarySt = std::to_string(pol_at_boundary);
        Polar << pol_at_boundarySt << std::endl;
        Tracking << xpdef[0] << "; " << xpdef[1] << std::endl;
        std::cout << "textfiles written" << '\n';
        }
        t_old=t; //set time

        if (minimum){std::cout << "minimum true" << '\n';} //print out whether minimum is reached
        if (ctr>800){std::cout << "ctr true" << '\n';} //print, whether respective countr is reached
        double anglesum_diff = fabs(anglesum_old - anglesum); //calculate difference in anglesum
        std::cout << "anglesum_diff = " << anglesum_diff << '\n'; //...and print this out

        //if not just start next step
        if (v_cl.rank()==0)
        {
            std::cout << "-----------------------------------------------------------------\n";
        }

        //reset values before next step
        anglesum_old = anglesum;
        anglesum =0.0;
        sum_fe = 0.0;
        sumsplay = 0.0;
        sumbend =0.0;
        Pol_old = Pol;
        f3_old = f3;
        t_old=t;

        gt.start(); //start wall clock timer again
    }
};

int main(int argc, char* argv[])
{
    {   openfpm_init(&argc,&argv); //initialize openfpm

        auto &v_cl = create_vcluster(); //create cluster

        tt2.start(); //start timer for whole simulation
        Gd = int(std::atof(argv[1])); //grid size
        double tf = std::atof(argv[2]); //maximum time
        double dt = tf/std::atof(argv[3]); //time step size
        wr_f=int(std::atof(argv[3])); 
        wr_at=int(std::atof(argv[4]));
        //give dimensionless value for slpay/bend as 5th value to program
        dkK = double(std::atof(argv[5]));
        //give f3 as 6th value for program
        double phi = M_PI*double(std::atof(argv[6]));
        V_err_eps=std::atof(argv[7]); //change dependent of Gd
        timetol=std::atof(argv[8]);
        double timereltol=std::atof(argv[9]);
        //give parameter fro underlying grid
        double boxsize = std::atof(argv[10]);
        //calculate Kb from given dimensionles parameter dkK
        Kb = (dkK + 1.) * Ks;
        zeta = std::atof(argv[11]);
        double ri = std::atof(argv[12])*Gd/2;
        lambda = std::atof(argv[13]);
        double zetadelmu = (zeta + nu*gama*lambda)*(delmu*(Gd/2)*(Gd/2))/Ks;
        

        //create grid
        const size_t sz[2] = {Gd, Gd};
        Box<2, double> box({0.0, 0.0}, {boxsize, boxsize});
        double Lx = box.getHigh(0),Ly = box.getHigh(1);
        size_t bc[2] = {NON_PERIODIC,NON_PERIODIC}; //set boundary condition
        spacing = box.getHigh(0) / (sz[0]); //get spacing
        double rCut = 3.6 * spacing; //define cut off radius
        // std::cout << "spacing in main = " << spacing << std::endl;
        int ord = 2; //order used in discretization for differential operators
        Ghost<2, double> ghost(rCut); //define ghost layer
        vector_dist_ws<2, double,Activegels> Particles(0, box, bc, ghost); //define vector for particles
        Particles.setPropNames(PropNAMES); //set names to names created above

        if(v_cl.rank()==0){
            std::cout<< "Gd="<<Gd<<'\n'<<"tf="<<tf<<'\n'<<"dt="<<dt<<'\n';
            std::cout<<"boxsize="<<boxsize<<'\n';
            std::cout<< "wr_f="<<wr_f<<'\n'<<"wr_at="<<wr_at<<'\n';
            std::cout << "Kb = " << Kb <<" ,dkK="<<dkK<<'\n';
            std::cout<< "phi(or psi0)="<<phi<<'\n';
            std::cout<<"V_err_eps="<<V_err_eps<<'\n'<<"timetol="<<timetol<<'\n'<< "timereltol="<<timereltol<<'\n';
            std::cout<< "zetadelmu="<<zetadelmu<<'\n'<<" delmu="<<delmu<<'\n'<<" zeta="<<zeta<<'\n'<<" lambda="<<lambda<<'\n';
            std::cout<< "ri/ro="<< std::atof(argv[12]) << "ri="<<ri<<'\n';
            std::cout << "spacing in main = " << spacing << std::endl;
        }

        double x0=15.0, y0 = 15.0, x1=boxsize, y1=boxsize; //basic box
        boxsizehalf = boxsize/2;
        double x_defect=boxsizehalf+2*spacing;
        double y_defect=boxsizehalf;
        Point<2,double> Xpn,xpd1={boxsizehalf,boxsizehalf}, xpd2={x_defect, y_defect}; //middle point

        //random number generator
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> gridrandomnes(-0.2, 0.2);

        //kernel
        auto it = Particles.getGridIterator(sz);
        while (it.isNext()) {
            auto key = it.get();
            double x = key.get(0) * spacing;
            double y = key.get(1) * spacing;
            //define cut off radius
            if(sqrt((boxsizehalf-x)*(boxsizehalf-x)+(boxsizehalf-y)*(boxsizehalf-y))<boxsizehalf-7*spacing/6.0)
            {
                Particles.add();
                //add randomness: othwise there divergence in x-/y-direction ar weighted more
                Particles.getLastPos()[0] = x + gridrandomnes(gen); 
                Particles.getLastPos()[1] = y + gridrandomnes(gen);
                Particles.getLastProp<NORMAL>()[0] = 0;
                Particles.getLastProp<NORMAL>()[1] = 0;
                Particles.getLastProp<PID>() = 0;
                Particles.getLastSubset(0);
            }
            ++it; //increment iterator
        }

        //define circle as grid for outer particles
        if (v_cl.rank()==0){
            int n_b=int(Gd)*boxsizehalf;
            double radius = boxsizehalf - 3*spacing/4.0;
            //double Golden_f3=M_PI * (3.0 - sqrt(5.0));
            double Golden_f3=2.0*M_PI/double(n_b);
                for(int i=1;i<=n_b;i++)
                {
                    double Golden_theta = Golden_f3 * i;
                    double x = boxsizehalf+cos(Golden_theta) * radius;
                    double y = boxsizehalf+sin(Golden_theta) * radius;
                    Particles.add();
                    Particles.getLastPos()[0] = x;
                    Particles.getLastPos()[1] = y;
                    Particles.getLastSubset(0);
                    Particles.getLastProp<NORMAL>()[0] = (x-boxsizehalf)/sqrt((x-boxsizehalf)*(x-boxsizehalf)+(y-boxsizehalf)*(y-boxsizehalf));
                    Particles.getLastProp<NORMAL>()[1] = (y-boxsizehalf)/sqrt((x-boxsizehalf)*(x-boxsizehalf)+(y-boxsizehalf)*(y-boxsizehalf));
                    Particles.getLastProp<PID>() = -1;
                }
        }
        Particles.map(); //map particles onto grid

        size_t pctr=1;
        vector_dist_subset<2, double, Activegels> Particles_bulk(Particles,0);

        //Dirichlet Boundary Conditions: set angle for beginning: Dirichlet BC imposed by CorrVec.
        auto & Bulk = Particles_bulk.getIds();
        double max = .1;
        std::uniform_real_distribution<> dis(0.0, max);
        for (int j = 0; j < Bulk.size(); j++) {
            auto p = Bulk.get<0>(j);
            Xpn={Particles.getPos(p)[0],Particles.getPos(p)[1]};
            double dist=Xpn.distance(xpd1);
            dist = (dist == 0) ? 1 : dist;
            double dist2 = Xpn.distance(xpd2);
            dist2 = (dist2==0) ? 1 : dist2;

            //calculate distance of defect to boundary along normal to defect
            double dist3;
            double theta_d = atan2(Xpn[1]-xpd2[1], Xpn[0]-xpd2[0]);
            double a= boxsizehalf-xpd2[0];
            dist3= a*cos(theta_d) + 0.5*sqrt(2*a*a*(cos(2*theta_d)-1)+4*pow(boxsizehalf, 2));



            //use commented if you want to have activity gradient following a certain region
            //change x*spacing for different radius of passive region
            if(dist>(boxsizehalf)||dist<ri*spacing){
                Particles.getProp<DELMU>(p)=0.0; //change to delmu, if you want to have continious activity
                Particles.getProp<ZETA>(p)=zeta;
            }else{Particles.getProp<DELMU>(p)=delmu;
            Particles.getProp<ZETA>(p)=zeta;
            }
            if(j==0){
                if(v_cl.rank()==0){
                std::cout << "zetadelmu= " << Particles.getProp<DELMU>(p) * (pow((Gd/2.), 2.) * (Particles.getProp<ZETA>(p) + nu * gama * lambda))/Ks << std::endl;
                }
            }
            
            //compute normal
            if(dist2>1){
                double snormal = tanh(0.5*(dist2))/(tanh(0.5*(dist3)));
                double x1 = (xpd2[0] - Xpn[0])/dist2 * snormal;
                double x2 = (xpd2[1] - Xpn[1])/dist2 * snormal;
                //set angle of rotation
                //rotate normal by angle
                Particles.getProp<POLARIZATION>(p)[x] = cos(phi)*x1 - sin(phi)*x2;
                Particles.getProp<POLARIZATION>(p)[y] = sin(phi)*x1 + cos(phi)*x2;
                Particles.getProp<F3>(p) = phi;

                //noise for angle in bulk
                if (Particles.getProp<PID>(p)==0){
                  double smooth = 1/(exp((dist - boxsize/4)*3) +1);
                  double v = (dis(gen)-max/2) * M_PI * smooth;//use for positive delta K
                  //double v = (dis(gen)) * M_PI * smooth;//use for negative delta K
                  Particles_bulk.getProp<POLARIZATION>(p)[x] = cos(v)*Particles_bulk.getProp<POLARIZATION>(p)[x] - sin(v)*Particles_bulk.getProp<POLARIZATION>(p)[y];
                  Particles_bulk.getProp<POLARIZATION>(p)[y] = sin(v)*Particles_bulk.getProp<POLARIZATION>(p)[x] + cos(v)*Particles_bulk.getProp<POLARIZATION>(p)[y];
                  Particles.getProp<F3>(p) = phi+v;
                }
            }
            //put defect
            if (dist2<=5){
                double snormal_defect = tanh(0.5*(dist2))/(tanh(0.5*(dist3)));
                double xd = (xpd2[0] - Xpn[0])/dist2 * snormal_defect;
                double yd = (xpd2[1] - Xpn[1])/dist2 * snormal_defect;
                double r = sqrt(x*x + y*y);
                Particles.getProp<POLARIZATION>(p)[x] = cos(phi)*xd - sin(phi)*yd;
                Particles.getProp<POLARIZATION>(p)[y] = sin(phi)*xd + cos(phi)*yd;
                Particles.getProp<F3>(p) = phi;
            }
            

            //connect particles for Dirichlet BC
            if (Particles.getProp<PID>(p)!=0){
                double xp = Xpn[0];
                double yp = Xpn[1];
                double theta=atan2(yp-boxsizehalf,xp-boxsizehalf);
                Particles.add();
                Particles.getLastPos()[0] = boxsizehalf*(1.0+cos(theta));
                Particles.getLastPos()[1] = boxsizehalf*(1.0+sin(theta));
                Particles.getLastSubset(1);
                double xb = (boxsizehalf - Xpn[0])/dist;
                double yb = (boxsizehalf - Xpn[1])/dist;
                double phi_bound = phi - 0.*M_PI;
                Particles.getLastProp<0>()[0] = cos(phi_bound)*xb - sin(phi_bound)*yb;
                Particles.getLastProp<0>()[1] = sin(phi_bound)*xb + cos(phi_bound)*yb;
                Particles.getLastProp<F3>() = phi_bound;
                Particles.getLastProp<DELMU>() = Particles.getProp<DELMU>(p);
                Particles.getLastProp<ZETA>() = Particles.getProp<ZETA>(p);
                Particles.getLastProp<PID>() = pctr;
                Particles.getProp<PID>(p) = pctr;
                pctr++;
            }
        }


        vector_dist_subset<2, double, Activegels> Particles_boundary(Particles, 1);
        Particles.ghost_get<0,12,DELMU,NORMAL, ZETA>();
        Particles.ghost_get<F3>();
        //insert randomness to angle -> not already start with steady state
        auto & bulk = Particles_bulk.getIds();
        auto & boundary = Particles_boundary.getIds();

        Particles.ghost_get<POLARIZATION, F3>(SKIP_LABELLING); //update ghost particles for polarization and angle

        //define parameter from fields with getV
        auto Pol = getV<POLARIZATION>(Particles);
        auto V = getV<VELOCITY>(Particles);
        auto P = getV<PRESSURE>(Particles);
        auto dPol = getV<DPOL>(Particles);
        auto Pol_old = getV<POLD>(Particles);
        auto f3 = getV<F3>(Particles);

        //set parameters to zero in the beginning
        P = 0;V = 0,dPol=0,Pol_old=0;
        sum_fe =0.;
        sumsplay = 0.0;
        sumbend = 0.0;
        anglesum = 0.0;
        anglesum_old =0.0;

        Particles.ghost_get<POLARIZATION,NORMAL,DELMU,DPOL, F3>(SKIP_LABELLING); //update ghost particles for polarization, nomral, delmu, dpol and angle

        // now link the partners of boundary and bulk
        //before PID for both particles at same time, bt still need to iterate, bc jusr to find partners
        //find or access single particles -> CorrVec
        CorrVec.clear(); //clear to make sure
        for(int i = 0; i < bulk.size(); i++)
        {
                auto p = bulk.get<0>(i); //0 for id by convention
                if(Particles.getProp<NORMAL>(p)[0]!=0 || Particles.getProp<NORMAL>(p)[1]!=0){
                    for(int j = 0; j < boundary.size(); j++) {
                        auto p_out = boundary.get<0>(j);
                        if(Particles.getProp<PID>(p_out)==Particles.getProp<PID>(p)){
                            CorrVec.add();
                            CorrVec.get<0>(CorrVec.size()-1)[0]=p_out;
                            CorrVec.get<0>(CorrVec.size()-1)[1]=p;
                        }
                    }
                }
        }

        std::cout << "CORRVEC SIZE " << CorrVec.size() << "\n";

        //define parameter from fields with getV
        auto P_bulk = getV<PRESSURE>(Particles_bulk);//Pressure only on inside
        auto Pol_bulk = getV<POLARIZATION>(Particles_bulk);
        auto dPol_bulk = getV<DPOL>(Particles_bulk);
        auto dPol_boundary = getV<DPOL>(Particles_boundary);
        auto dV_bulk = getV<DV>(Particles_bulk);
        auto RHS_all = getV<VRHS>(Particles);
        auto RHS_bulk = getV<VRHS>(Particles_bulk);
        auto div_bulk = getV<DIV>(Particles_bulk);


        Particles.write("Init"); //write initialization pvtp file

        //define differential operators
        Derivative_x Dx(Particles,ord,rCut,3.1,support_options::RADIUS); 
        Derivative_y Dy(Particles, ord, rCut,3.1,support_options::RADIUS); 
        Derivative_xy Dxy(Particles, ord, rCut,3.1,support_options::RADIUS);
        auto Dyx = Dxy; //symmetry
        Derivative_xx Dxx(Particles, ord, rCut,1.9,support_options::RADIUS);
        Derivative_yy Dyy(Particles, ord, rCut,1.9,support_options::RADIUS);

        //call numerical methods from odeint library
        boost::numeric::odeint::adams_bashforth_moulton<2,state_type_2d_ofp,double,state_type_2d_ofp,double,boost::numeric::odeint::vector_space_algebra_ofp> abm;
        boost::numeric::odeint::adaptive_adams_bashforth_moulton<2, state_type_2d_ofp,double,state_type_2d_ofp,double,boost::numeric::odeint::vector_space_algebra_ofp > abmA;
        
        //ids for components of equations
        eq_id x_comp, y_comp; 
        x_comp.setId(0);
        y_comp.setId(1);

        petsc_solver<double> solverPetsc; //define solver
        DCPSE_scheme<equations2d2, vector_type> Solver(Particles); //define solver woth DCPSE scheme

        //sink and source free fluid
        auto Stokes1 = eta * Dxx(V[x]) + eta * Dyy(V[x]);
        auto Stokes2 = eta * Dxx(V[y]) + eta * Dyy(V[y]);
        //impose equations
        Solver.impose(Stokes1, bulk, 0, x_comp);
        Solver.impose(Stokes2, bulk, 0, y_comp);
        Solver.impose(V[x], boundary, 0, x_comp);
        Solver.impose(V[y], boundary, 0, y_comp);
        //...and solve them with petsc
        Solver.solve_with_solver(solverPetsc, V[x], V[y]);

        //get size of bulk
        vsz=bulk.size(); 
        v_cl.sum(vsz);
        v_cl.execute();

        //define vectors to call particles/solver
        vectorGlobal=(void *) &Particles;
        vectorGlobal_bulk=(void *) &Particles_bulk;
        vectorGlobal_boundary=(void *) &Particles_boundary;
        vectorGlobal_boundary=(void *) &Particles_boundary;
        SolverPointer=(void *) &Solver;
        SolverPointerpetsc=(void *) &solverPetsc;


        //create instances of functors
        PolarEv<Derivative_x,Derivative_y,Derivative_xx,Derivative_xy,Derivative_yy> System(Dx,Dy,Dxx,Dxy,Dyy);
        ObserverFunctor<Derivative_x,Derivative_y,Derivative_xx,Derivative_xy,Derivative_yy> Observer(Dx,Dy,Dxx,Dxy,Dyy);

        //define p and derivative of p for solver
        state_type_2d_ofp tPol;
        tPol.data.get<0>()=Pol[x];
        tPol.data.get<1>()=Pol[y];

        //equation ids
        eq_id vx, vy;
        vx.setId(0);
        vy.setId(1);

        timer tt; //timer for whole simulation

        //set start values for parameter
        dPol = Pol;
        double V_err = 1, V_err_old;
        double tim=0; //for giving to odeint library

        // for(int k=0; k<100; k++){

        // }
        //call odeint library for solving the equations
        size_t steps; //declare type
        steps=boost::numeric::odeint::integrate_adaptive(boost::numeric::odeint::make_controlled(timetol,timereltol,abmA),System,tPol,tim,tf,dt,Observer);

        if(v_cl.rank()==0){std::cout << "Time steps: " << steps << std::endl;} //print time steps

        //define polarization at bulk
        Pol_bulk[x]=tPol.data.get<0>();
        Pol_bulk[y]=tPol.data.get<1>();

        //write pvtp file
        Particles.deleteGhost();
        Particles.write("Polar_Last");
        Dx.deallocate(Particles);
        Dy.deallocate(Particles);
        Dxy.deallocate(Particles);
        Dxx.deallocate(Particles);
        Dyy.deallocate(Particles);
        std::cout.precision(17); //specify precision for writing
        tt2.stop(); //stop timer for whole simulation
        //...and print info
        if (v_cl.rank() == 0) {
            std::cout << "The simulation took " << tt2.getcputime() << "(CPU) ------ " << tt2.getwct()
                      << "(Wall) Seconds.";
        }
    }
    //close txt files
    AvAngle.close();
    FesSplay.close();
    FesBend.close();
    Fes.close();
    Time.close();

    openfpm_finalize(); //finalize OpenFPM
}
