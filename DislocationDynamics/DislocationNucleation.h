#include <iostream>
#include <Eigen/Dense>
#include <Polycrystal.h>


namespace model
{
    template <int dim>
    class DislocationNucleation
    {
    private:
        /* data */
    public:

        typedef Eigen::Matrix<double,1,6> VogitStress;
        typedef Eigen::Matrix<double,dim,1> VectorDim;

        const Polycrystal<dim>& poly;




        const double kB_eV=8.617e-5;       // Boltzmann constant in [eV]
        const double kB_SI=1.38064852e-23; // Boltzmann constant in SI units [J/K]
        double Tm_K=poly.Tm;                //1141;
        double sigma_th=5.2/86.0;
        double v0=1.0e13;
        double A_eV=4.8;
        double alpha=4.1;
        double vvSt=0.0;


    
        DislocationNucleation(  const double & dt_in,
                                const Polycrystal<dim>& _poly,

        ):
        /* init */
        {

        }

        void equivlantStress(VogitStress& stressV)
        {
            Eigen::Matrix<double,3,3> stress;
            stress<<stressV[1],stressV[4],stressV[6],
                    stressV[4],stressV[2],stressV[5],
                    stressV[6],stressV[5],stressV[3];
            // std::mt19937 generator(std::chrono::system_clock::now().time_since_epoch().count());
            // const int rSS=distribution(generator); // a random SlipSystem ID
            const int rSS=19; // a specific SlipSystem ID
            const int grainID=0;
            const SlipSystem& slipSystem(*poly.grain(grainID).slipSystems()[rSS]);
            const VectorDim b=slipSystem.s.cartesian();    // Burgers vector

            
        }


        void nucleationProbability(const double& dt)
        {
            double T_K=poly.T;
            double    sigma=externalLoadController->externalStress().trace();
            //double    sigma= 1.0/86.0;
            std::cout<<"sigma= "<<sigma<<std::endl;
            double Q=(1.0-T_K/Tm_K)*A_eV*pow(1-sigma/sigma_th, alpha);
            double vv=v0*exp(-Q/kB_eV/T_K);
            //double prob=vv*dt*7.506e-14;
            double vvS=vv*1.005e9;//Fe
            // double vvSt=vvSt+vvS*dt;//1e-8stressRate//1e-10-2
            //double vvSt=vvSt+vvS*dt*0.01;//1e-10stressRate//not good
            //std::cout<<"externalLoadController->externalStressRate().trace() = "<<externalLoadController->externalStressRate().trace()<<std::endl;
            //double vvSt=vvSt+vvS*dt*externalLoadController->externalStressRate().trace()/1.0e-8;//1e-10-3
            //double vvSt=vvSt+vvS*dt*0.1;//1e-10-4
            double vvSt=vvSt+vvS*dt*0.01;//1e-10-5
            std::cout<<"Q = "<<Q<<std::endl;
            std::cout<<"vv = "<<vv<<std::endl;
            std::cout<<"vvS = "<<vvS<<std::endl;
            std::cout<<"vvSt = "<<vvSt<<std::endl;

            while(vvSt>1.0)
            {
                addShearLoopWhileLoad(1.0,150.0);
                vvSt=vvSt-1.0;
                std::cout<<"vvSt = "<<vvSt<<std::endl;   
            }
            /* //for (const auto& eIter :  bvpSolver->finiteElement()->elements())
               std::cout<<"000"<< bvpSolver->stress().onBoundary()<<std::endl;
             for (const auto&  eIter : bvpSolver->finiteElement().elements())
            {
                std::cout<<"AAA"<<std::endl;
                if(eIter.second.isBoundaryElement())
                {
                    std::cout<<"BBB"<<std::endl;
                }
           }*/ 
        }










        ~DislocationNucleation();
    };

}