/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDtimeIntegrator_H_
#define model_DDtimeIntegrator_H_

#include <chrono>

#include <Eigen/Dense>


namespace model
{
	
	/**********************************************************************/
	/**********************************************************************/
	template <int N>
	struct DDtimeIntegrator
    {
		
		
	};
    
    /**********************************************************************/
    /**********************************************************************/
    template <>
    struct DDtimeIntegrator<0>
    {
        static constexpr auto tag="vMax integrator";
        static double dxMax;
        static double shearWaveSpeedFraction;
        
        
        /******************************************************************/
        static void initFromFile(const std::string& fileName)
        {
            
            dxMax=TextFileParser(fileName).readScalar<double>("dxMax",true);
            assert(dxMax>0.0);
            
        }
        
        /**********************************************************************/
        template <typename DislocationNetworkType>
        static double getGlideTimeIncrement(const DislocationNetworkType& DN)
        {
            //! Compute and store DislocaitonNode velocities
//            DN.assembleAndSolve(runID,straightSegmentsDeq);
            
            /*! Computes the time step size \f$dt\f$ for the current simulation step,
             *  based on maximum nodal velocity \f$v_{max}\f$.
             *
             *  The time step is calculated according to:
             *	\f[
             *  dt=
             *  \begin{cases}
             *		\frac{dx}{v_{max}} & v_{max} > fc_s\\
             *      \frac{dx}{fc_s} & v_{max} \le fc_s\\
             *  \end{cases}
             *	\f]
             *  where \f$c_s\f$ is the shear velocity and \f$f=0.1\f$ is a constant.
             */

            double vmax=0.0;

            for (const auto& nodeIter : DN.nodes())
            {
                if(

                   nodeIter.second->glidePlanes().size()<3

                   )
                {

                    const double vNorm(nodeIter.second->get_V().norm());

                    if (vNorm>vmax)
                    {
                        vmax=vNorm;

                    }
                }
            }
            
            return vmax > DN.poly.cs*shearWaveSpeedFraction? dxMax/vmax : dxMax/(DN.poly.cs*shearWaveSpeedFraction);
            
        }
        
        /**********************************************************************/
        template <typename DislocationNetworkType>
        static double getClimbTimeIncrement(const DislocationNetworkType& DN)
        {
            //! Compute and store DislocaitonNode velocities
            //            DN.assembleAndSolve(runID,straightSegmentsDeq);
            
            /*! Computes the time step size \f$dt\f$ for the current simulation step,
             *  based on maximum nodal velocity \f$v_{max}\f$.
             *
             *  The time step is calculated according to:
             *    \f[
             *  dt=
             *  \begin{cases}
             *        \frac{dx}{v_{max}} & v_{max} > fc_s\\
             *      \frac{dx}{fc_s} & v_{max} \le fc_s\\
             *  \end{cases}
             *    \f]
             *  where \f$c_s\f$ is the shear velocity and \f$f=0.1\f$ is a constant.
             */

            double vmax=0.0;

            for (const auto& nodeIter : DN.nodes())
            {
                if(

                   nodeIter.second->glidePlanes().size()<3

                   )
                {

                    const double vNorm(nodeIter.second->get_V().norm());

                    if (vNorm>vmax)
                    {
                        vmax=vNorm;

                    }
                }
            }
            
            return dxMax/vmax;

        }
        
    };

#ifndef _MODEL_GREATWHITE_
    double DDtimeIntegrator<0>::dxMax=10.0;
    double DDtimeIntegrator<0>::shearWaveSpeedFraction=1.0e-3;
#endif

} // end namespace
#endif

