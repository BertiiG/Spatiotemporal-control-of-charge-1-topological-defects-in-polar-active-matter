//
// Created by Birte Geerds in 2023, based on Code from Abhinav Singh on 15.03.20.

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
#include <thread>
#include <mutex>

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
double delmu;
double zeta=-1.;
double passive_old;
double passive_setpoint;
int ctr=0;

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

class PIController {
public:
    PIController(double kp, double ki, double outputMin, double outputMax)
        : kp(kp), ki(ki), integral(0.0), outputMin(outputMin), outputMax(outputMax) {}

    // Explicitly delete copy constructor and copy assignment operator
    // PIController(const PIController&) = delete;
    // PIController& operator=(const PIController&) = delete;

    double compute(double setpoint, double processVariable) {
        //std::lock_guard<std::mutex> lock(mutex);

        double error = setpoint - processVariable;
        integral += error;

        //anti-windup machanism
        // Adjust integral term based on saturation limits
        if (integral > outputMax) {
            integral = outputMax;  // Limit the integral term
        } else if (integral < outputMin) {
            integral = outputMin;  // Limit the integral term
        }

        double output = kp * error + ki * integral;
        output = std::max(outputMin, std::min(outputMax, output));
        return output;
    }

private:
    double kp, ki;
    double integral;
    double outputMin, outputMax;
};

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
        //auto zeta = getV<ZETA>(Particles);

        //impose Dirichlet BC for polarisation given the boundary and bulk pairs
        for (int j = 0; j < CorrVec.size(); j++) {
                auto p_out = CorrVec.get<0>(j)[0];
                auto p = CorrVec.get<0>(j)[1];
                Particles.getProp<0>(p)[x] = Particles.getProp<0>(p_out)[x];
                Particles.getProp<0>(p)[y] = Particles.getProp<0>(p_out)[y];
                Particles.getProp<F3>(p) = Particles.getProp<F3>(p_out);
                Particles.getProp<DELMU>(p_out) = Particles.getProp<DELMU>(p);
            }

        //X data for the solver is polarization
        Pol[x]=X.data.get<0>();
        Pol[y]=X.data.get<1>();
        //ids for equations
        eq_id x_comp, y_comp;
        x_comp.setId(0);
        y_comp.setId(1);

        int n = 0,nmax = 10000,errctr = 0, Vreset = 0; //namx...maximal loops aload; errctr...maximal loop for Verr>Verr_old; Vreset...1 if next loop allowed, 0 if not
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
                   + 0.5*nu * u[x][x] * Pol[x] + 0.5*nu * u[y][y] * Pol[x] //positive should be right?
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
    // int ctr; //counter for wiritng (into txt or pvtp) files
    double x_old;
    double y_old;

    // PI Controller
    //PIController piController_delmu;
    PIController piController_dist;

    //Constructor
    //ObserverFunctor(DX &Dx,DY &Dy,DXX &Dxx,DXY &Dxy,DYY &Dyy, double kp, double ki, double outputMin, double outputMax, double ki_dist, double kp_dist, double outputMin_dist, double outputMax_dist):Dx(Dx),Dy(Dy),Dxx(Dxx),Dxy(Dxy),Dyy(Dyy), t_old(0.0), ctr(0), piController_delmu(kp, ki, outputMin, outputMax), piController_dist(kp_dist, ki_dist, outputMin_dist, outputMax_dist)
    ObserverFunctor(DX &Dx,DY &Dy,DXX &Dxx,DXY &Dxy,DYY &Dyy, double ki_dist, double kp_dist, double outputMin_dist, double outputMax_dist):Dx(Dx),Dy(Dy),Dxx(Dxx),Dxy(Dxy),Dyy(Dyy), t_old(0.0), piController_dist(kp_dist, ki_dist, outputMin_dist, outputMax_dist)
    {
        t_old = 0.0;
        // ctr = 0;
        x_old=15*spacing;
        y_old=boxsizehalf;
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

        //find global minimum and broadcast to all cores
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
        }
        xpdef = {x_def_global, y_def_global};
		
        // std::cout <<"L. 585" << std::endl;
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
        //std::cout << "point of the defect is: x="<< xpdef[0] << " , y=" <<xpdef[1]<< "\n";
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

        //calculate velocity of defect
        //calculate distance between last time step and actual time step
        double travelled_distance = sqrt(pow((x_def_global-x_old), 2) + pow((y_def_global-y_old),2));
        double velocity = travelled_distance/(t-t_old);

        size_t v_ctr=0;
        double v_x = 0, v_y=0;
        auto it=Particles.getDomainIterator();
        std::cout << "spacing= " << spacing << std::endl;

        while(it.isNext()){
            auto p=it.get();
            double xp = Particles.getPos(p)[x];
            double yp = Particles.getPos(p)[y];
            double dist_to_def = sqrt(pow(xp-x_def_global, 2) + pow(yp-y_def_global, 2));
            //std::cout<< "dist_to_def= " << dist_to_def << std::endl;
            if(dist_to_def<spacing*1.5){
                //std::cout << "region around defect defined" << std::endl;
                v_x+=Particles.getProp<VELOCITY>(p)[x];
                v_y+=Particles.getProp<VELOCITY>(p)[y];
                v_ctr++;
            }

            ++it;
        }
        v_cl.sum(v_x);
        v_cl.sum(v_y);
        v_cl.sum(v_ctr);
        v_cl.execute();

        std::cout<<"v_ctr = " << v_ctr << std::endl;
        double v_norm=0.0;
        v_norm =sqrt(pow((v_x/v_ctr), 2) + pow((v_y/v_ctr), 2));
        if(v_cl.rank()==0){std::cout << "v_norm = " << v_norm << std::endl;}

        double zetadelmu_old= delmu * (pow((Gd/2.), 2.) * (zeta + nu * gama * lambda))/Ks;
        if(v_cl.rank()==0){std::cout << "zetadelmu_old= " << zetadelmu_old << std::endl;}
        double zetadelmu_new=zetadelmu_old;

        //example on how to influence angular velocity of the defect by setting the activity if the active region based on no-crossterms ansatz
        //and hwo to change the setpoint for r_d for different times
        double setpoint1 = 7;
        double setpoint2 = 14;
        if(t<60){
            passive_setpoint = setpoint2;
        }

        if(t>=60){
            passive_setpoint = setpoint1;
        }
        if(t>=180){
            passive_setpoint = setpoint2;
        }
        if(t>=240){
            passive_setpoint = setpoint1;
        }
        if(t>=360){
            passive_setpoint = setpoint2;
        }
        if(t>=420){
            passive_setpoint = setpoint1;
        }
        if(t>=540){
            passive_setpoint = setpoint2;
        }
        if(t>=600){
            passive_setpoint = setpoint1;
        }
        if(t>=720){
            passive_setpoint = setpoint2;
        }
        if(t>=780){
            passive_setpoint = setpoint1;
        }
        if(t>=900){
            passive_setpoint = setpoint2;
        }
        if(t>=960){
            passive_setpoint = setpoint1;
        }
        if(t>=1080){
            passive_setpoint = setpoint2;
        }
        if(t>=1140){
            passive_setpoint = setpoint1;
        }
        if(t>=1260){
            passive_setpoint = setpoint2;
        }
        if(t>=1320){
            passive_setpoint = setpoint1;
        }
        if(t>=1440){
            passive_setpoint = setpoint2;
        }
        if(t>=1500){
            passive_setpoint = setpoint1;
        }
        if(t>=1620){
            passive_setpoint = setpoint2;
        }
        if(t>=1680){
            passive_setpoint = setpoint1;
        }
        if(t>=1800){
            passive_setpoint = setpoint2;
        }
        if(t>=1860){
            passive_setpoint = setpoint1;
        }
        if(t>=1980){
            passive_setpoint = setpoint2;
        }
        if(t>=2040){
            passive_setpoint = setpoint1;
        }
        if(t>=2160){
            passive_setpoint = setpoint2;
        }

        //PI controller for distance
        double dist_to_center = sqrt(pow((x_def_global-boxsizehalf), 2) + pow((y_def_global-boxsizehalf), 2));
        if(v_cl.rank()==0){std::cout << "passive_old= " << passive_old << std::endl;}

        if(t>=1. && std::fmod(ctr, 10.)<=1e-4 && v_norm!=0){
            double output_dist = piController_dist.compute(passive_setpoint, dist_to_center);
            if(v_cl.rank()==0){std::cout<< "output_dist= " <<output_dist << std::endl;}
            double passive_new = passive_old + output_dist;
            delmu = (zetadelmu_new * Ks)/(pow((Gd/2.), 2.) * (-1. + nu * gama * lambda));

            //define center point
            Point<2,double> x_particle,x_center={boxsizehalf,boxsizehalf};

            for (int j = 0; j < bulk.size(); j++) {
                auto p = bulk.get<0>(j);
                //calculate distance of point
                x_particle={Particles.getPos(p)[0],Particles.getPos(p)[1]};
                double dist=x_particle.distance(x_center);
                dist = (dist == 0) ? 1 : dist;

                if(dist>(boxsizehalf)||dist<passive_new){
                    Particles.getProp<DELMU>(p)=0.0;
                    Particles.getProp<ZETA>(p)=-1;
                }else{Particles.getProp<DELMU>(p)=delmu;
                Particles.getProp<ZETA>(p)=-1;
                }
                if(v_cl.rank()==0){
                    if(j==0){
                        std::cout << "zetadelmu_new= " << Particles.getProp<DELMU>(p) * (pow((Gd/2.), 2.) * (Particles.getProp<ZETA>(p) + nu * gama * lambda))/Ks << std::endl;
                        std::cout << "passive_new= " << passive_new << std::endl;
                    }
                }
            }

            passive_old = passive_new;

        }
        

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
        x_old=x_def_global;
        y_old=y_def_global;

        

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
        //calculate delmu from given dimensionles parameter zetadelmu
        zeta = std::atof(argv[11]);
        double ri = std::atof(argv[12])*Gd/2;
        lambda = std::atof(argv[13]);
        delmu = std::atof(argv[14]);
        double zetadelmu = (zeta + nu*gama*lambda)*(delmu*(Gd/2)*(Gd/2))/Ks;
        
        //give set point for the distance r_d (defect center to plane center) as 12th value
        passive_setpoint = std::atof(argv[15]);
        

        //create grid
        const size_t sz[2] = {Gd, Gd};
        Box<2, double> box({0.0, 0.0}, {boxsize, boxsize});
        double Lx = box.getHigh(0),Ly = box.getHigh(1);
        size_t bc[2] = {NON_PERIODIC,NON_PERIODIC}; //set boundary condition
        spacing = box.getHigh(0) / (sz[0]); //get spacing
        double rCut = 3.6 * spacing; //define cut off radius
        int ord = 2; //order used in discretization for differential operators
        Ghost<2, double> ghost(rCut); //define ghost layer
        vector_dist_ws<2, double,Activegels> Particles(0, box, bc, ghost); //define vector for particles
        Particles.setPropNames(PropNAMES); //set names to names created above

        passive_old = 3.9*spacing;
        //PIController for distance
        double kp_dist=0.0005;
        double ki_dist=0.005;
        double outputMin_dist = -2.*spacing;
        double outputMax_dist = 2.*spacing;
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
            std::cout << "passive_old= " << passive_old << std::endl;
            std::cout << "kp_dist= " << kp_dist << "ki_dist" << ki_dist << std::endl;
        }
   

        double x0=0.0, y0 = 0.0, x1=boxsize, y1=boxsize; //basic box
        boxsizehalf = boxsize/2;
        double x_defect=1.28*boxsizehalf; //already start with defect at rd=7
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



            //use commented when you want to have activity gradient
            if(dist>(boxsizehalf)||dist<passive_old){
                Particles.getProp<DELMU>(p)=0.0;
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
        ObserverFunctor<Derivative_x,Derivative_y,Derivative_xx,Derivative_xy,Derivative_yy> Observer(Dx,Dy,Dxx,Dxy,Dyy, kp_dist, ki_dist, outputMin_dist, outputMax_dist);

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
        double t_change1=120;
        double phi_bound2 = 0.*M_PI;
        double t_change2=180;
        double phi_bound3 = phi;//0.025*M_PI;
        double t_change3=300;
        double phi_bound4 = 0.0*M_PI;
        double t_change4=360;
        double phi_bound5 = phi;//0.025*M_PI;
        double t_change5=480;
        double phi_bound6 = 0.0*M_PI;
        double t_change6=540;
        double phi_bound7 = phi;//0.025*M_PI;
        double t_change7=660;
        double phi_bound8 = 0.0*M_PI;
        double t_change8=720;
        double phi_bound9 = phi;//0.025*M_PI;
        double t_change9=840;
        double phi_bound10 = 0.0*M_PI;
        double t_change10=900;
        double phi_bound11 = phi;//0.025*M_PI;
        double t_change11=1020;
        double phi_bound12 = 0.0*M_PI;
        double t_change12=1080;
        double phi_bound13 = phi;//0.025*M_PI;
        double t_change13=1200;
        double phi_bound14 = 0.0*M_PI;
        double t_change14=1260;
        double phi_bound15 = phi;//0.025*M_PI;
        double t_change15=1380;
        double phi_bound16 = 0.0*M_PI;
        double t_change16=1440;
        double phi_bound17 = phi;//0.025*M_PI;
        double t_change17=1560;
        double phi_bound18 = 0.0*M_PI;
        double t_change18=1620;
        double phi_bound19 = phi;//0.025*M_PI;
        double t_change19=1740;
        double phi_bound20 = 0.0*M_PI;
        double t_change20=1800;
        double phi_bound21 = phi;//0.025*M_PI;
        double t_change21=1920;
        double phi_bound22 = 0.0*M_PI;
        double t_change22=1980;
        double phi_bound23 = phi;//0.025*M_PI;
        double t_change23=2100;
        double phi_bound24 = 0.0*M_PI;
        double t_change24=2160;
        double phi_bound25 = phi;//0.025*M_PI;
        double t_change25=tf;

        //call odeint library for solving the equations
        size_t steps; //declare type
        steps=boost::numeric::odeint::integrate_adaptive(boost::numeric::odeint::make_controlled(timetol,timereltol,abmA),System,tPol,tim,t_change1,dt,Observer);
        
        std::cout << "Time steps: " << steps << std::endl; //print time steps
        //define polarization at bulk
        Pol_bulk[x]=tPol.data.get<0>();
        Pol_bulk[y]=tPol.data.get<1>();

        std::vector<std::pair<double,double>> changes = {
            {t_change2, phi_bound2},
            {t_change3, phi_bound3},
            {t_change4, phi_bound4},
            {t_change5, phi_bound5},
            {t_change6, phi_bound6},
            {t_change7, phi_bound7},
            {t_change8, phi_bound8},
            {t_change9, phi_bound9},
            {t_change10, phi_bound10},
            {t_change11, phi_bound11},
            {t_change12, phi_bound12},
            {t_change13, phi_bound13},
            {t_change14, phi_bound14},
            {t_change15, phi_bound15},
            {t_change16, phi_bound16},
            {t_change17, phi_bound17},
            {t_change18, phi_bound18},
            {t_change19, phi_bound19},
            {t_change20, phi_bound20},
            {t_change21, phi_bound21},
            {t_change22, phi_bound22},
            {t_change23, phi_bound23},
            {t_change24, phi_bound24},
            {t_change25, phi_bound25},
        };

        double t_current = t_change1;

        for (auto &[t_next, phi_bound] : changes)
        {
            // --- Update boundary at t_current with phi_bound ---
            auto &boundary = Particles_boundary.getIds();

            if (t_next > tf) {
                // last integration up to tf, then stop
                t_next = tf;
            }
            for (int j = 0; j < boundary.size(); j++) {
                auto p_out = boundary.get<0>(j);

                double Xpos = Particles.getPos(p_out)[0];
                double Ypos = Particles.getPos(p_out)[1];
                double dist = sqrt((Xpos-boxsizehalf)*(Xpos-boxsizehalf) + (Ypos-boxsizehalf)*(Ypos-boxsizehalf));
                dist = (dist == 0) ? 1 : dist;

                double x1 = (boxsizehalf - Xpos)/dist;
                double x2 = (boxsizehalf - Ypos)/dist;

                double new_px = cos(phi_bound)*x1 - sin(phi_bound)*x2;
                double new_py = sin(phi_bound)*x1 + cos(phi_bound)*x2;

                Particles.getProp<POLARIZATION>(p_out)[0] = new_px;
                Particles.getProp<POLARIZATION>(p_out)[1] = new_py;
                Particles.getProp<F3>(p_out)              = phi_bound;
            }

            // enforce boundary → bulk
            for (int i = 0; i < CorrVec.size(); ++i) {
                auto p_boundary = CorrVec.get<0>(i)[0];
                auto p_bulk     = CorrVec.get<0>(i)[1];
                Particles.getProp<POLARIZATION>(p_bulk) = Particles.getProp<POLARIZATION>(p_boundary);
                Particles.getProp<F3>(p_bulk)           = Particles.getProp<F3>(p_boundary);
            }
            Particles.ghost_get<POLARIZATION, F3>(SKIP_LABELLING);

            // sync back into ODE state
            tPol.data.get<0>() = Pol[x];
            tPol.data.get<1>() = Pol[y];

            // --- Integrate until next change ---
            size_t steps = boost::numeric::odeint::integrate_adaptive(
                boost::numeric::odeint::make_controlled(timetol, timereltol, abmA),
                System, tPol, t_current, t_next, dt, Observer);

            std::cout << "Integrated " << steps << " steps from t=" << t_current
                      << " to t=" << t_next << " with phi=" << phi_bound << std::endl;

            // update bulk polarization storage
            Pol_bulk[x] = tPol.data.get<0>();
            Pol_bulk[y] = tPol.data.get<1>();

            // advance current time
            t_current = t_next;
        }


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