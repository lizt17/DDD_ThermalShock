/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2017 by Yinan Cui <cuiyinan@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DefectiveCrystal_H_
#define model_DefectiveCrystal_H_

#include <iostream>
#include <vector>
#include <memory>
#include <DefectiveCrystalParameters.h>                 
#include <DislocationNetwork.h>

namespace model
{
    
    template <int _dim, short unsigned int corder, typename InterpolationType>
    class DefectiveCrystal
    {
        
    public:
        static constexpr int dim=_dim; // make dim available outside class
        typedef DefectiveCrystal<dim,corder,InterpolationType> DefectiveCrystalType;
        typedef DislocationNetwork<dim,corder,InterpolationType> DislocationNetworkType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;         
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef BVPabaqus<dim,1> BVPabaqusType;
        typedef typename BVPabaqusType::ElementType ElementType;  // use ABAQUS to solve FEM
        
        DefectiveCrystalParameters simulationParameters;    // DefectiveCrystalParameters为结构体，视为特殊的类   
                
        const SimplicialMesh<dim> mesh;
        const std::vector<VectorDim> periodicShifts;
        const Polycrystal<dim> poly;
        const std::unique_ptr<DislocationNetworkType> DN;
        const std::unique_ptr<BVPabaqusType> bvpAbaqus;
        //const std::unique_ptr<ExternalLoadControllerBase<dim>> externalLoadController;

        /**********************************************************************/
        static std::vector<VectorDim> getPeriodicShifts(const SimplicialMesh<dim>& m,
                                                        const DefectiveCrystalParameters& params)
        {
            // Set up periodic shifts
            std::vector<VectorDim> temp;
            if(params.simulationType==DefectiveCrystalParameters::PERIODIC_IMAGES)
            {
                const VectorDim meshDimensions(m.xMax()-m.xMin());
                model::cout<<"meshDimensions="<<meshDimensions.transpose()<<std::endl;
                for(int i=-params.periodicImages_x;i<=params.periodicImages_x;++i)
                {
                    for(int j=-params.periodicImages_y;j<=params.periodicImages_y;++j)
                    {
                        for(int k=-params.periodicImages_z;k<=params.periodicImages_z;++k)
                        {
                            const Eigen::Array<int,dim,1> cellID((Eigen::Array<int,dim,1>()<<i,j,k).finished());
                            temp.push_back((meshDimensions.array()*cellID.template cast<double>()).matrix());
                        }
                    }
                }
            }
            else
            {
                temp.push_back(VectorDim::Zero());
            }
            
            model::cout<<"periodic shift vectors:"<<std::endl;
            for(const auto& shift : temp)
            {
                model::cout<<shift.transpose()<<std::endl;
                
            }
            
            return temp;
            
        }

        /**********************************************************************/
        void updateAbaqusInput(const long int& runID)
        {/*! Updates bvpAbaqus using the stress and displacement fields of the
          *  current DD configuration.
          *  updates inc plastic strain
          */
            if(bvpAbaqus)
            {
                if (!(runID%bvpAbaqus->para.stepsBetweenBVPupdates))
                {// enter the if statement if use_bvp!=0 and runID is a multiple of use_bvp
                    model::cout<<"		Updating bvp ABAQUS stress... "<<std::endl;
                    bvpAbaqus->updateStress(runID);
                    //bvpAbaqus->clearPlasticStrain();
                    //bvpAbaqus->updateDisplacement();
                }
                // if (runID==0)
                // {
                //     bvpAbaqus->updateIntPoint();
                // }
                
            }

        }

        /**********************************************************************/
        DefectiveCrystal(int& argc, char* argv[]) :
        /* init */ simulationParameters(argc,argv)
        /* init */,mesh(TextFileParser("./inputFiles/polycrystal.txt").readString("meshFile",true))
        /* init */,periodicShifts(getPeriodicShifts(mesh,simulationParameters))
        /* init */,poly("./inputFiles/polycrystal.txt",mesh)
        /* init */,DN(simulationParameters.useDislocations? new DislocationNetworkType(argc,argv,simulationParameters,mesh,poly,bvpAbaqus,periodicShifts,simulationParameters.runID) : nullptr)
        /* init */,bvpAbaqus(simulationParameters.simulationType==DefectiveCrystalParameters::ABAQUS? new BVPabaqusType(mesh,*DN) : nullptr)

        {
            assert(mesh.simplices().size() && "MESH IS EMPTY.");
            
            
            if(  simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_IMAGES
               || simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_FEM  )
               
            {
                assert(poly.grains().size()==1 && "ONLY SINGLE-CRYSTAL PERIODIC SIMULATIONS SUPPORTED.");
                
                for(const auto& rIter : mesh.regions())
                {
                    for(const auto& pair : rIter.second->parallelFaces())
                    {
                        model::cout<<"Checking if parallel faces "<<pair.first<<"<->"<<pair.second<<" are commensurate"<<std::endl;
                        const PlanarMeshFace<dim>& face1(*rIter.second->faces().at(pair.first));
                        const PlanarMeshFace<dim>& face2(*rIter.second->faces().at(pair.second));
                        const VectorDim cc(face1.center()-face2.center());
                        const VectorDim ccc(cc.dot(face1.outNormal())*face1.outNormal());
                        
                        const LatticeDirection<dim> ld(poly.grains().begin()->second.latticeDirection(face1.outNormal()));
                        const double normRatio(ccc.norm()/ld.cartesian().norm());
                        if(std::fabs(std::round(normRatio)-normRatio)>FLT_EPSILON)
                        {
//                            std::cout<<"Face outNormal="<<std::setprecision(15)<<std::scientific<<face1.outNormal().transpose()<<std::endl;
                            std::cout<<"Mesh in direction "<< std::setprecision(15)<<std::scientific<<ld.cartesian().normalized().transpose()<<" is not commensurate for periodicity"<<std::endl;
                            std::cout<<"Mesh size in that direction must be a multiple of "<< std::setprecision(15)<<std::scientific<<ld.cartesian().norm()<<std::endl;
                            std::cout<<"Size detected="<< std::setprecision(15)<<std::scientific<<ccc.norm()<<std::endl;
                            std::cout<<"Closest commensurate size="<< std::setprecision(15)<<std::scientific<<std::round(normRatio)*ld.cartesian().norm()<<std::endl;
                            assert(false && "MESH NOT COMMENSURATE");
                        }
                    }
                }
            }
        }

        /**********************************************************************/
        void singleGlideStep()
        {
            model::cout<<blueBoldColor<< "runID="<<simulationParameters.runID<<" (of "<<simulationParameters.Nsteps<<")"
            /*                    */<< ", time="<<simulationParameters.totalTime;
            if(DN)
            {
                model::cout<< ": nodes="<<DN->nodes().size()
                /*                    */<< ", segments="<<DN->links().size()
                /*                    */<< ", loopSegments="<<DN->loopLinks().size()
                /*                    */<< ", loops="<<DN->loops().size()
                /*                    */<< ", components="<<DN->components().size();
            }
            model::cout<< defaultColor<<std::endl;
            
            if(DN)
            {
                // calculate plastic strain in volume
                DN->updateGeometry(simulationParameters.dt);
                
                updateAbaqusInput(simulationParameters.runID);

                // DN->dislocationNucleation(simulationParameters.dt);
                
                DN->assembleAndSolveGlide(simulationParameters.runID);
                // control time increment
                simulationParameters.dt=DDtimeIntegrator<0>::getGlideTimeIncrement(*DN); // TO DO: MAKE THIS std::min between DN and CrackSystem
                if ( simulationParameters.userdt > 0.0 )
                {
                    simulationParameters.dt=simulationParameters.userdt;
                }
                // DN->localPlasticStrain(simulationParameters.dt);
                // output
                DN->io().output(simulationParameters.runID);
                // move dislocation
                DN->moveGlide(simulationParameters.dt);               
                // consider velocity reduce
                DN->localPlasticStrain_2(simulationParameters.dt);
                // // output
                // DN->io().output(simulationParameters.runID);
                // menage discrete topological events，Junction
                DN->singleGlideStepDiscreteEvents(simulationParameters.runID);
                if (!(simulationParameters.runID%simulationParameters.stepsBetweenNucleation))
                {
                    DN->dislocationNucleation(simulationParameters.dt);

                    //DN->NucSite.clear();
                }
            }
            simulationParameters.totalTime+=simulationParameters.dt;
            ++simulationParameters.runID;
        }
        
        /**********************************************************************/
        void runGlideSteps()
        {/*! Runs a number of simulation time steps defined by simulationParameters.Nsteps
          */
            const auto t0= std::chrono::system_clock::now();
            // double incTime(0.0);
            while (simulationParameters.runID<simulationParameters.Nsteps)
            // for (int runk=0;runk<simulationParameters.Nsteps;runk++)
            {
                model::cout<<std::endl; // leave a blank line
                singleGlideStep();

                // incTime+=simulationParameters.dt;
            }
            model::cout<<"total slipped aera = "<<DN->totalS<<std::endl;
            model::cout<<greenBoldColor<<std::setprecision(3)<<std::scientific<<simulationParameters.Nsteps<< " simulation steps completed in "<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" [sec]"<<defaultColor<<std::endl;
        }

        /**********************************************************************/
        VectorDim displacement(const VectorDim& x) const
        {/*!\param[in] P position vector
          * \returns The displacement field in the DefectiveCrystal at P
          */
            VectorDim temp(VectorDim::Zero());
            if(DN)
            {
                temp+=DN->displacement(x);
            }
            return temp;
        }
        
        /**********************************************************************/
        void displacement(std::vector<FEMnodeEvaluation<ElementType,dim,1>>& fieldPoints) const
        {
            if(DN)
            {
                DN->displacement(fieldPoints);
            }
        }
        
        /**********************************************************************/
        MatrixDim stress(const VectorDim& x) const
        {/*!\param[in] P position vector
          * \returns The stress field in the DefectiveCrystal at P
          * Note:
          */
            MatrixDim temp(MatrixDim::Zero());
            if(DN)
            {
                temp+=DN->stress(x);
            }
            return temp;
        }
        
        /**********************************************************************/
        MatrixDim plasticDistortion() const
        {/*!\param[in] P position vector
          * \returns The stress field in the DefectiveCrystal at P
          * Note:
          */
            MatrixDim temp(MatrixDim::Zero());
            if(DN)
            {
                temp+=DN->plasticDistortion();
            }
            return temp;
        }
        
        /**********************************************************************/
        MatrixDim plasticDistortionRate() const
        {/*!\param[in] P position vector
          * \returns The stress field in the DefectiveCrystal at P
          * Note:
          */
            MatrixDim temp(MatrixDim::Zero());
            if(DN)
            {
                temp+=DN->plasticDistortionRate();
            }
            return temp;
        }
        
        /**********************************************************************/
        MatrixDim plasticStrainRate() const
        {/*!\param[in] P position vector
          * \returns The stress field in the DefectiveCrystal at P
          * Note:
          */
            MatrixDim temp(plasticDistortionRate());
            return 0.5*(temp+temp.transpose());
        }
        
    };
}
#endif       