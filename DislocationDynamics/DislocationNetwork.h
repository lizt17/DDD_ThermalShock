/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po             <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez       <ramirezbrf@gmail.com>
 * Copyright (C) 2011 by Mamdouh Mohamed        <msm07d@fsu.edu>
 * Copyright (C) 2011 by Tamer Crsoby           <tamercrosby@gmail.com>
 * Copyright (C) 2011 by Can Erel               <canerel55@gmail.com>
 * Copyright (C) 2011 by Yinan Cui              <cuiyinan@ucla.edu>
 * Copyright (C) 2017 by Sabyasachi Chatterjee  <sabyasac@andrew.cmu.edu>
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

// valgrind --leak-check=full --show-leak-kinds=all ./DDomp

#ifndef model_DISLOCATIONNETWORK_H_
#define model_DISLOCATIONNETWORK_H_

#ifdef _MODEL_MPI_
#define _MODEL_DD_MPI_
#endif


#ifdef _OPENMP
#include <omp.h>
#define EIGEN_DONT_PARALLELIZE // disable Eigen Internal openmp Parallelization
#endif

#include <vector>
#include <chrono>
#include <map>
#include <memory>
#include <Eigen/Dense>
#include <LoopNetwork.h>
#include <TerminalColors.h>
#include <DislocationNetworkTraits.h>
#include <DislocationNetworkComponent.h>
#include <DislocationNode.h>
#include <DislocationSegment.h>
#include <DislocationLoop.h>
#include <DislocationNetworkRemesh.h>
#include <DislocationJunctionFormation.h>
#include <DislocationCrossSlip.h>
#include <DislocationNetworkIO.h>
#include <DislocationParticle.h>
#include <DislocationStress.h>
#include <MPIcout.h>
#include <DDtimeIntegrator.h>
#include <EqualIteratorRange.h>
#include <GrainBoundaryTransmission.h>
#include <BVPabaqus.h>

#include <Polycrystal.h>
#include <DislocationNodeContraction.h>
#include <EshelbyInclusion.h>
#include <TextFileParser.h>
#include <DefectiveCrystalParameters.h>

#include <DislocationInjector.h>
#include <PeriodicDislocationLoop.h>


#ifdef _MODEL_GREATWHITE_
#include <MooseSolution.h>
#endif


namespace model
{
    
    
    
    template <int _dim, short unsigned int _corder, typename InterpolationType>
    class DislocationNetwork :
    /*                      */public LoopNetwork<DislocationNetwork<_dim,_corder,InterpolationType> >
    /*                      */,public std::map<size_t,EshelbyInclusion<_dim>>
    {
        
        
    public:
        
        
        static constexpr int dim=_dim; // make dim available outside class
        static constexpr int corder=_corder; // make dim available outside class
        typedef DislocationNetwork<dim,corder,InterpolationType> DislocationNetworkType;
        typedef LoopNetwork<DislocationNetworkType> LoopNetworkType;
        typedef DislocationNode<dim,corder,InterpolationType> NodeType;
        typedef DislocationSegment<dim,corder,InterpolationType> LinkType;
        typedef DislocationNetworkComponent<NodeType,LinkType> DislocationNetworkComponentType;
        typedef NetworkComponent<NodeType,LinkType> NetworkComponentType;
        typedef Eigen::Matrix<double,dim,dim>	MatrixDimD;
        typedef Eigen::Matrix<double,dim,1>		VectorDim;
        typedef Eigen::Matrix<double,1,6> VoigtStress;
        typedef typename TypeTraits<DislocationNetworkType>::LoopType LoopType;

        typedef BVPabaqus<dim,1> BvpAbaqusType;
        typedef typename BvpAbaqusType::FiniteElementType FiniteElementType;
        typedef typename FiniteElementType::ElementType ElementType;
        typedef typename LoopNetworkType::IsNodeType IsNodeType;
        typedef DislocationNetworkIO<DislocationNetworkType> DislocationNetworkIOType;
        typedef Polycrystal<dim> PolycrystalType;

        typedef std::map<size_t,EshelbyInclusion<_dim>> EshelbyInclusionContainerType;
        typedef NetworkLinkObserver<LinkType> NetworkLinkObserverType;
        typedef typename NetworkLinkObserverType::LinkContainerType NetworkLinkContainerType;
        typedef PeriodicDislocationLoop<DislocationNetworkType> PeriodicDislocationLoopType;
        typedef Eigen::Matrix<double,dim,1> VectorDimD;


        
#ifdef DislocationNucleationFile
#include DislocationNucleationFile
#endif
        
#ifdef _MODEL_GREATWHITE_
#include <DislocationNetworkGreatWhite.h>
#endif
        
    private:
        
        
        /**********************************************************************/
        void updateVirtualBoundaryLoops()
        {
            
            
            // if(   simulationParameters.simulationType==DefectiveCrystalParameters::FINITE_FEM
            //    || simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_FEM
            //    || simulationParameters.simulationType==DefectiveCrystalParameters::ABAQUS       )
            if(   simulationParameters.simulationType==DefectiveCrystalParameters::FINITE_FEM
               || simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_FEM     )
               
            {               
                                    
                const auto t0= std::chrono::system_clock::now();
                model::cout<<"        Updating virtual boundary loops "<<std::flush;
                                
                // First clean up outdated boundary loops
                std::set<size_t> removeLoops;
                for(const auto& loop : this->loops())
                {
                    if((loop.second->isVirtualBoundaryLoop() && loop.second->links().size()!=4) || loop.second->isPureVirtualBoundaryLoop())
                    {// clean up left over loops from topological operations
                        removeLoops.insert(loop.second->sID);
                    }
                }
                
                for(const size_t& loopID : removeLoops)
                {// Remove the virtual loops with ID in removeLoops
                    this->deleteLoop(loopID);
                }
                 
                                    
                // Now reconstruct virtual boundary loops
                std::vector<std::tuple<std::vector<std::shared_ptr<NodeType>>,VectorDim,size_t,int>> virtualLoopVector;
                for(const auto& link : this->links())
                {
                    if(link.second->isBoundarySegment() && !link.second->hasZeroBurgers())
                    {
                        virtualLoopVector.emplace_back(std::vector<std::shared_ptr<NodeType>>{link.second->sink,link.second->source,link.second->source->virtualBoundaryNode(),link.second->sink->virtualBoundaryNode()},
                                                    link.second->burgers(),
                                                    (*link.second->grains().begin())->grainID,
                                                    DislocationLoopIO<dim>::VIRTUALLOOP);
                    }
                }

                
                for(const auto& tup : virtualLoopVector)
                {// Insert the new virtual loops
                    this->insertLoop(std::get<0>(tup),std::get<1>(tup),std::get<2>(tup),std::get<3>(tup));
                }
                model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;

            }
        }

    public:
        
        const DefectiveCrystalParameters& simulationParameters;
        const SimplicialMesh<dim>& mesh;
        const Polycrystal<dim>& poly;
        GlidePlaneFactory<dim> glidePlaneFactory;
        const std::unique_ptr<PeriodicDislocationLoopFactory<DislocationNetworkType>> periodicDislocationLoopFactory;
        const std::unique_ptr<BVPabaqus<dim,1>>& bvpAbaqus;

        const std::vector<VectorDim>& periodicShifts;
        DislocationNetworkRemesh<DislocationNetworkType> networkRemesher;
        DislocationJunctionFormation<DislocationNetworkType> junctionsMaker;
        DislocationNodeContraction<DislocationNetworkType> nodeContractor;
        GrainBoundaryTransmission<DislocationNetworkType> gbTransmission;

        MatrixDimD _plasticDistortionFromAreas;
        MatrixDimD _plasticDistortionRateFromVelocities;
        MatrixDimD _plasticDistortionRateFromAreas;
        int ddSolverType;
        bool computeDDinteractions;
        int crossSlipModel;
        int  outputFrequency;
        bool outputBinary;
        bool outputGlidePlanes;

        bool outputElasticEnergy;
        bool outputMeshDisplacement;
        bool outputFEMsolution;
        bool outputDislocationLength;
        bool outputPlasticDistortion;
        bool outputPlasticDistortionRate;
        bool outputQuadraturePoints;
        bool outputLinkingNumbers;
        bool outputLoopLength;
        bool outputSegmentPairDistances;
        bool outputPeriodicConfiguration;

        bool use_stochasticForce;
        double surfaceAttractionDistance;

        std::string folderSuffix;
        std::vector<VectorDim> NucSite;

        //luosc
        const double epsilon; //used to regluarized the plastic strain accprding to YP JMPS paper. 
        const double expandheight;
        const double ele_h;
        const double nc;
        const double stress_nuc;

        // 2021-5-14, lizt for test
        double totalS;
        // Frank-Read sources
        //const double targetFrankReadDislocationDensityLoad;
        const double FrankReadSizeMeanLoad;
        const double FrankReadSizeStdLoad;
        const double FrankReadAspectRatioMeanLoad;
        const double FrankReadAspectRatioStdLoad;
        
        /**********************************************************************/
        DislocationNetwork(int& argc, char* argv[],
                           const DefectiveCrystalParameters& _simulationParameters,
                           const SimplicialMesh<dim>& _mesh,
                           const Polycrystal<dim>& _poly,
                           const std::unique_ptr<BVPabaqus<dim,1>>& _bvpAbaqus,
                           const std::vector<VectorDim>& _periodicShifts,
                           long int& runID) :
        /* init */ simulationParameters(_simulationParameters)
        /* init */,mesh(_mesh)
        /* init */,poly(_poly)
        /* init */,glidePlaneFactory(poly)
        /* init */,periodicDislocationLoopFactory(simulationParameters.isPeriodicSimulation()? new PeriodicDislocationLoopFactory<DislocationNetworkType>(poly,glidePlaneFactory) : nullptr)
        /* init */,bvpAbaqus(_bvpAbaqus)
        /* init */,periodicShifts(_periodicShifts)
        /* init */,networkRemesher(*this)
        /* init */,junctionsMaker(*this)
        /* init */,nodeContractor(*this)
        /* init */,gbTransmission(*this)
        /* init */,_plasticDistortionFromAreas(MatrixDimD::Zero())
        /* init */,_plasticDistortionRateFromVelocities(MatrixDimD::Zero())
        /* init */,_plasticDistortionRateFromAreas(MatrixDimD::Zero())
        /* init */,ddSolverType(TextFileParser("./inputFiles/DD.txt").readScalar<int>("ddSolverType",true))
        /* init */,computeDDinteractions(TextFileParser("./inputFiles/DD.txt").readScalar<int>("computeDDinteractions",true))
        /* init */,crossSlipModel(TextFileParser("./inputFiles/DD.txt").readScalar<int>("crossSlipModel",true))
        /* init */,outputFrequency(TextFileParser("./inputFiles/DD.txt").readScalar<int>("outputFrequency",true))
        /* init */,outputBinary(TextFileParser("./inputFiles/DD.txt").readScalar<int>("outputBinary",true))
        /* init */,outputGlidePlanes(TextFileParser("./inputFiles/DD.txt").readScalar<int>("outputGlidePlanes",true))
        /* init */,outputElasticEnergy(TextFileParser("./inputFiles/DD.txt").readScalar<int>("outputElasticEnergy",true))
        /* init */,outputMeshDisplacement(TextFileParser("./inputFiles/DD.txt").readScalar<int>("outputMeshDisplacement",true))
        /* init */,outputFEMsolution(TextFileParser("./inputFiles/DD.txt").readScalar<int>("outputFEMsolution",true))
        /* init */,outputDislocationLength(TextFileParser("./inputFiles/DD.txt").readScalar<int>("outputDislocationLength",true))
        /* init */,outputPlasticDistortion(TextFileParser("./inputFiles/DD.txt").readScalar<int>("outputPlasticDistortion",true))
        /* init */,outputPlasticDistortionRate(TextFileParser("./inputFiles/DD.txt").readScalar<int>("outputPlasticDistortionRate",true))
        /* init */,outputQuadraturePoints(TextFileParser("./inputFiles/DD.txt").readScalar<int>("outputQuadraturePoints",true))
        /* init */,outputLinkingNumbers(TextFileParser("./inputFiles/DD.txt").readScalar<int>("outputLinkingNumbers",true))
        /* init */,outputLoopLength(TextFileParser("./inputFiles/DD.txt").readScalar<int>("outputLoopLength",true))
        /* init */,outputSegmentPairDistances(TextFileParser("./inputFiles/DD.txt").readScalar<int>("outputSegmentPairDistances",true))
        /* init */,outputPeriodicConfiguration(simulationParameters.isPeriodicSimulation()? TextFileParser("./inputFiles/DD.txt").readScalar<int>("outputPeriodicConfiguration",true) : false)
        /* init */,use_stochasticForce(TextFileParser("./inputFiles/DD.txt").readScalar<int>("use_stochasticForce",true))
        /* init */,surfaceAttractionDistance(TextFileParser("./inputFiles/DD.txt").readScalar<double>("surfaceAttractionDistance",true))
        /* init */,folderSuffix("")
        //luosc//
         /* init  */,epsilon(TextFileParser("./inputFiles/DD.txt").readScalar<double>("epsilon",true))  
        /* init  */,expandheight(TextFileParser("./inputFiles/DD.txt").readScalar<double>("expandheight",true)) 
        /* init  */,ele_h(TextFileParser("./inputFiles/DD.txt").readScalar<double>("ele_h",true)) 
       /* init  */,nc(TextFileParser("./inputFiles/DD.txt").readScalar<double>("g_nc",true))
        /* init*/,FrankReadSizeMeanLoad(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("FrankReadSizeMean",true))
        /* init*/,FrankReadSizeStdLoad(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("FrankReadSizeStd",true) )
        /* init*/,FrankReadAspectRatioMeanLoad( TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("FrankReadAspectRatioMean",true) )
        /* init*/,FrankReadAspectRatioStdLoad( TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("FrankReadAspectRatioStd",true))
        /* init*/,stress_nuc( TextFileParser("./inputFiles/DD.txt").readScalar<double>("stress_nuc",true))
        {
            
            // Some sanity checks
            
            // Initialize static variables
            LinkType::initFromFile("./inputFiles/DD.txt");
            NodeType::initFromFile("./inputFiles/DD.txt");
            LoopType::initFromFile("./inputFiles/DD.txt");
            PeriodicDislocationBase::initFromFile("./inputFiles/DD.txt");
            DislocationNetworkComponentType::initFromFile("./inputFiles/DD.txt");
            DislocationStressBase<dim>::initFromFile("./inputFiles/DD.txt");
            DDtimeIntegrator<0>::initFromFile("./inputFiles/DD.txt");
            DislocationCrossSlip<DislocationNetworkType>::initFromFile("./inputFiles/DD.txt");
            int stochasticForceSeed=TextFileParser("./inputFiles/DD.txt").readScalar<int>("stochasticForceSeed",true);
            if(stochasticForceSeed<0)
            {
                StochasticForceGenerator::init(std::chrono::system_clock::now().time_since_epoch().count());
            }
            else
            {
                StochasticForceGenerator::init(stochasticForceSeed);
            }
            
            if(argc>1)
            {
                folderSuffix=argv[1];

            }

            // Read Vertex and Edge information
            DDconfigIO<dim> evl(folderSuffix);
            evl.read(runID);
            setConfiguration(evl);
            createEshelbyInclusions();
        }
        
        /**********************************************************************/
        void setConfiguration(const DDconfigIO<dim>& evl)
        {
            this->loopLinks().clear(); // erase base network
            std::map<size_t,std::shared_ptr<NodeType>> tempNodes;
            createVertices(evl,tempNodes);
            createEdges(evl,tempNodes);
            updatePlasticDistortionFromAreas(simulationParameters.dt);
#ifdef _MODEL_MPI_
            // Avoid that a processor starts writing before other are done reading
            MPI_Barrier(MPI_COMM_WORLD);
#endif
            // Initializing configuration
//            this->io().output(simulationParameters.runID);
            moveGlide(0.0);    // initial configuration
        }


        /**********************************************************************/
        void createVertices(const DDconfigIO<dim>& evl,std::map<size_t,std::shared_ptr<NodeType>>& tempNodes)
        {/*!Creates DislocationNode(s) based on the data read by the DDconfigIO<dim>
          * object.
          */
            size_t kk(1);
            for (const auto& node : evl.nodes())
            {
                const size_t nodeIDinFile(node.sID);
                NodeType::set_count(nodeIDinFile);
                if(node.sID==node.masterID)
                {// a regular node is created
                    model::cout<<"Creating DislocationNode "<<nodeIDinFile<<" ("<<kk<<" of "<<evl.nodes().size()<<")"<<std::endl;
                    
                    const size_t nodeID(StaticID<NodeType>::nextID());
                    const auto inserted(tempNodes.emplace(std::piecewise_construct,
                                                          std::make_tuple(nodeID),
                                                          std::make_tuple(new NodeType(this,node.P,node.V,node.velocityReduction))));
                    assert(inserted.second && "COULD NOT INSERT NETWORK VERTEX IN VERTEX CONTAINER.");
                    assert(inserted.first->first == nodeID && "KEY != nodeID");
                    assert(inserted.first->second->sID == nodeID && "sID != nodeID");
                    assert(nodeID==nodeIDinFile);
                }
                else
                {
                    model::cout<<"Creating DislocationNode "<<nodeIDinFile<<" ("<<kk<<" of "<<evl.nodes().size()<<"), virtual of "<<node.masterID<<std::endl;
                    const auto isNode(this->node(node.masterID));
                    assert(isNode.first);
                    isNode.second->resetVirtualBoundaryNode();
                }
                kk++;
            }
        }
        
        std::shared_ptr<NodeType> getSharedNode(const size_t& nodeID,
                                                const std::map<size_t,std::shared_ptr<NodeType>>& tempNodes)
        {
            const auto isNode(this->node(nodeID));
            assert(isNode.first);
            if(isNode.second->masterNode)
            {// a virtual node
                return isNode.second->masterNode->virtualBoundaryNode();
                //                loopNodes.push_back(isNode.second->masterNode->virtualBoundaryNode());
            }
            else
            {
                const auto tempNodesFound(tempNodes.find(nodeID));
                if(tempNodesFound==tempNodes.end())
                {
                    model::cout<<"node "<<nodeID<<" not found"<<std::endl;
                    assert(false && "node shared pointer not found");
                    return nullptr;
                }
                else
                {
                    return tempNodesFound->second;
                }
                //                loopNodes.push_back(isSharedNode.second);
            }
        }
        
        /**********************************************************************/
        void createEdges(const DDconfigIO<dim>& evl,const std::map<size_t,std::shared_ptr<NodeType>>& tempNodes)
        {/*!
          */
            
            
//            std::map<size_t,std::shared_ptr<PeriodicDislocationLoopType> > periodicDislocationLoopMap;
//            if(   simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_IMAGES
//               || simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_FEM)
//            {
//                for(const auto& pLoop : evl.periodicLoops())
//                {// Collect LoopLinks by loop IDs
//                    StaticID<PeriodicDislocationLoopType>::set_count(pLoop.sID);
//                    periodicDislocationLoopMap.emplace(pLoop.sID,new PeriodicDislocationLoopType(this));
//                }
//            }
            
            std::map<size_t,std::map<size_t,size_t>> loopMap;
            for(const auto& looplink : evl.links())
            {// Collect LoopLinks by loop IDs
                loopMap[looplink.loopID].emplace(looplink.sourceID,looplink.sinkID);
            }
            
            assert(loopMap.size()==evl.loops().size());
            
            size_t loopLumber=1;
            for(const auto& loop : evl.loops())
            {// for each loop in the DDconfigIO<dim> object
                
                const auto loopFound=loopMap.find(loop.sID); // there must be an entry with key loopID in loopMap
                assert(loopFound!=loopMap.end());
                std::vector<std::shared_ptr<NodeType>> loopNodes;
                loopNodes.push_back(getSharedNode(loopFound->second.begin()->first,tempNodes));
                for(size_t k=0;k<loopFound->second.size();++k)
                {
                    const auto nodeFound=loopFound->second.find(loopNodes.back()->sID);
                    if(k<loopFound->second.size()-1)
                    {
                        loopNodes.push_back(getSharedNode(nodeFound->second,tempNodes));
                    }
                    else
                    {
                        assert(nodeFound->second==loopNodes[0]->sID);
                    }
                }
                
                LoopType::set_count(loop.sID);
                
                const bool faulted= poly.grain(loop.grainID).rationalLatticeDirection(loop.B).rat.asDouble()!=1.0? true : false;
                
                model::cout<<"Creating Dislocation Loop "<<loop.sID<<" ("<<loopLumber<<" of "<<evl.loops().size()<<"), type="<<loop.loopType<<", faulted="<<faulted<<", |b|="<<loop.B.norm()<<std::endl;
                switch (loop.loopType)
                {
                    case DislocationLoopIO<dim>::GLISSILELOOP:
                    {
                        LatticePlane loopPlane(loop.P,poly.grain(loop.grainID).reciprocalLatticeDirection(loop.N));
                        GlidePlaneKey<dim> loopPlaneKey(loop.grainID,loopPlane);
                        if(simulationParameters.isPeriodicSimulation())
                        {
                            const size_t newLoopID=this->insertLoop(loopNodes,loop.B,glidePlaneFactory.get(loopPlaneKey),loop.periodicShift)->sID;
                            assert(loop.sID==newLoopID);
                        }
                        else
                        {
                            const size_t newLoopID=this->insertLoop(loopNodes,loop.B,glidePlaneFactory.get(loopPlaneKey))->sID;
                            assert(loop.sID==newLoopID);
                        }
//                        if(loop.periodicLoopID>=0)
//                        {// a loop belonging to a periodic loop
//                            const auto pLoopIter(periodicDislocationLoopMap.find(loop.periodicLoopID));
//                            if(pLoopIter!=periodicDislocationLoopMap.end())
//                            {
//                                const size_t newLoopID=this->insertLoop(loopNodes,loop.B,glidePlaneFactory.get(loopPlaneKey),pLoopIter->second,loop.periodicShift)->sID;
//                                assert(loop.sID==newLoopID);
//                            }
//                            else
//                            {
//                                assert(false && "PeriodicLoop not found in map");
//                            }
//                        }
//                        else
//                        {
//                            const size_t newLoopID=this->insertLoop(loopNodes,loop.B,glidePlaneFactory.get(loopPlaneKey))->sID;
//                            assert(loop.sID==newLoopID);
//                        }
                        break;
                    }
                        
                    case DislocationLoopIO<dim>::SESSILELOOP:
                    {
                        const size_t newLoopID=this->insertLoop(loopNodes,loop.B,loop.grainID,loop.loopType)->sID;
                        assert(loop.sID==newLoopID);
                        break;
                    }
                        
                    case DislocationLoopIO<dim>::VIRTUALLOOP:
                    {
                        const size_t newLoopID=this->insertLoop(loopNodes,loop.B,loop.grainID,loop.loopType)->sID;
                        assert(loop.sID==newLoopID);
                        break;
                    }
                        
                    default:
                        assert(false && "Unknown DislocationLoop type");
                        break;
                }
                
                loopLumber++;
            }
            
            model::cout<<std::endl;
        }
        
        /**********************************************************************/
        void createEshelbyInclusions()
        {
            IDreader<'E',1,14,double> inclusionsReader;
            inclusionsReader.read(0,true);
            
            const std::vector<double> inclusionsMobilityReduction(TextFileParser("./inputFiles/initialMicrostructure.txt").readArray<double>("inclusionsMobilityReduction",true));
            for(const auto& pair : inclusionsReader)
            {
                
                const size_t& inclusionID(pair.first);
                Eigen::Map<const Eigen::Matrix<double,1,14>> row(pair.second.data());
                
                const VectorDim C(row.template segment<dim>(0));
                const double a(row(dim+0));
                MatrixDimD eT(MatrixDimD::Zero());
                const int typeID(row(13));
                int k=dim+1;
                for(int i=0;i<dim;++i)
                {
                    for(int j=0;j<dim;++j)
                    {
                        eT(i,j)=row(k);
                        k++;
                    }
                }
                
                
                
                EshelbyInclusion<dim>::set_count(inclusionID);
                eshelbyInclusions().emplace(std::piecewise_construct,
                                            std::make_tuple(inclusionID),
                                            std::make_tuple(C,a,eT,poly.nu,poly.mu,inclusionsMobilityReduction[typeID],typeID) );
            }
        }
        
        /**********************************************************************/
        bool remove(const size_t& nodeID)
        {
            const auto isNode=this->node(nodeID);
            if(isNode.first)
            {// remove virtual node together with current node
                if(isNode.second->virtualBoundaryNode())
                {
                    LoopNetworkType::remove(isNode.second->virtualBoundaryNode()->sID);
                }
            }
            return LoopNetworkType::remove(nodeID);
        }
        
        /**********************************************************************/
        void updateGeometry(const double& dt)
        {
            for(auto& loop : this->loops())
            {// copmute slipped areas and right-handed normal // TODO: PARALLELIZE THIS LOOP
                loop.second->updateGeometry();
            }
            updatePlasticDistortionFromAreas(dt);
        }
        
        //        /**********************************************************************/
        //        double get_dt() const
        //        {
        //            switch (timeIntegrationMethod)
        //            {
        //                case 0:
        //                    return DDtimeIntegrator<0>::get_dt(*this);
        //                    break;
        //
        //                default:
        //                    assert(0 && "time integration method not implemented");
        //                    return 0;
        //                    break;
        //            }
        //        }
        
        /**********************************************************************/
        void removeZeroAreaLoops()
        {
            const auto t0= std::chrono::system_clock::now();
            model::cout<<"        Removing zero-area loops "<<std::flush;
            std::deque<size_t> loopIDs;
            for(const auto& loop : this->loops())
            {
                if(loop.second->slippedArea()<FLT_EPSILON)
                {
                    loopIDs.push_back(loop.second->sID);
                }
            }
            
            for(const auto& loopID : loopIDs)
            {
                this->deleteLoop(loopID);
            }
            model::cout<<"("<<loopIDs.size()<<" removed)"<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
        }
        
        /**********************************************************************/
        void singleGlideStepDiscreteEvents(const long int& runID)
        {
            //! 10- Cross Slip (needs upated PK force)
            DislocationCrossSlip<DislocationNetworkType>(*this);
            
            //! 12- Form Junctions
           junctionsMaker.formJunctions(DDtimeIntegrator<0>::dxMax);                    
            
            //! 13- Node redistribution
            networkRemesher.remesh(runID);

            updateVirtualBoundaryLoops();

        }
        
        /**********************************************************************/
        const EshelbyInclusionContainerType& eshelbyInclusions() const
        {
            return *this;
        }
        
        EshelbyInclusionContainerType& eshelbyInclusions()
        {
            return *this;
        }

        /**********************************************************************/
        bool contract(std::shared_ptr<NodeType> nA,
                      std::shared_ptr<NodeType> nB)
        {
            return nodeContractor.contract(nA,nB);
        }      
        
        
        /**********************************************************************/
        void assembleAndSolveGlide(const long int& runID)
        {/*! Performs the following operatons:
          */
#ifdef _OPENMP
            const size_t nThreads = omp_get_max_threads();
#else
            const size_t nThreads = 1;
#endif
            
            //! -1 Compute the interaction StressField between dislocation particles
            if(corder==0)
            {// For straight segments use analytical expression of stress field

                const auto t1= std::chrono::system_clock::now();
                model::cout<<"		Computing analytical stress field at quadrature points ("<<nThreads<<" threads) "<<std::flush;
#ifdef _OPENMP
                //                const size_t nThreads = omp_get_max_threads();
                EqualIteratorRange<typename NetworkLinkContainerType::iterator> eir(this->links().begin(),this->links().end(),nThreads);
                //                SOMETHING WRONG HERE? CPU USE SAYS 100% ONLY?
                
#pragma omp parallel for
                for(size_t thread=0;thread<eir.size();thread++)
                {
                    for (auto& linkIter=eir[thread].first;linkIter!=eir[thread].second;linkIter++)
                    {
                        linkIter->second->assembleGlide();
                    }
                }
#else
                for (auto& linkIter : this->links())
                {
                    linkIter.second->assembleGlide();
                }
#endif
                model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]."<<defaultColor<<std::endl;
                
            }
            else
            {// For curved segments use quandrature integration of stress field
                
                assert(0 && "ALL THIS MUST BE RE-IMPLEMENTED FOR CURVED SEGMENTS");
                
            }
            
            //! -3 Loop over DislocationSubNetworks, assemble subnetwork stiffness matrix and force vector, and solve
            const auto t2= std::chrono::system_clock::now();
            model::cout<<"		Assembling NetworkComponents and solving "<<std::flush;       
            
            switch (ddSolverType)
            {
                case 1: // iterative solver
                {
                    model::cout<<"(MINRES "<<nThreads<<" threads)..."<<std::flush;
#ifdef _OPENMP // MINRES is not multi-threaded. So parallelize loop over NetworkComponents.
#pragma omp parallel for
                    for (unsigned int k=0;k<this->components().size();++k)
                    {
                        auto snIter(this->components().begin());
                        std::advance(snIter,k);
                        DislocationNetworkComponentType(*snIter->second).iterativeSolve();
                    }
#else
                    for (const auto& networkComponent : this->components())
                    {
                        DislocationNetworkComponentType(*networkComponent.second).iterativeSolve();
                    }
#endif
                    break;
                }
                    
                case 2: // direct solver
                {
#ifdef _MODEL_PARDISO_SOLVER_
                    model::cout<<"(PardisoLDLT "<<nThreads<<" threads)..."<<std::flush;
                    for (const auto& networkComponent : this->components())
                    {
                        DislocationNetworkComponentType(*networkComponent.second).directSolve();
                    }
#else
                    model::cout<<"(SimplicialLDLT "<<nThreads<<" threads)..."<<std::flush;
#ifdef _OPENMP // SimplicialLDLT is not multi-threaded. So parallelize loop over NetworkComponents.
#pragma omp parallel for
                    for (unsigned int k=0;k<this->components().size();++k)
                    {
                        auto snIter(this->components().begin());
                        std::advance(snIter,k);
                        DislocationNetworkComponentType(*snIter->second).directSolve();
                    }
#else
                    for (const auto& networkComponent : this->components())
                    {
                        DislocationNetworkComponentType(*networkComponent.second).directSolve();
                    }
#endif
#endif
                    break;
                }
                    
                default: // lumped solver
                {
                    model::cout<<"(lumpedSolver "<<nThreads<<" threads)..."<<std::flush;
#ifdef _OPENMP
#pragma omp parallel for
                    for (unsigned int k=0;k<this->components().size();++k)
                    {
                        auto snIter(this->components().begin());
                        std::advance(snIter,k);
                        DislocationNetworkComponentType(*snIter->second).lumpedSolve(runID);
                    }
#else
                    for (const auto& networkComponent : this->components())
                    {
                        DislocationNetworkComponentType(*networkComponent.second).lumpedSolve(runID);
                    }
#endif
                    break;
                }
            }
            
            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t2)).count()<<" sec]."<<defaultColor<<std::endl;
        }
                
        
        /**********************************************************************/
        void moveGlide(const double & dt_in)
        {/*! Moves all nodes in the DislocationNetwork using the stored velocity and current dt
          */
            
            model::cout<<"        Moving DislocationNodes by glide (dt="<<dt_in<< ")... "<<std::flush;
            const auto t0= std::chrono::system_clock::now();
            
            if(simulationParameters.isPeriodicSimulation())
            {
                for (auto& nodeIter : this->nodes())
                {
                    static_cast<typename NodeType::NodeBaseType* const>(nodeIter.second)->set_P(nodeIter.second->get_P()+nodeIter.second->get_V()*dt_in);
                }
                
                for(const auto& pair : *periodicDislocationLoopFactory)
                {// output periodic glide planes too
                    
                    if(!pair.second.expired())
                    {
                        const auto periodicLoop(pair.second.lock());
                        periodicLoop->updateRVEloops(*this);
                    }
                }
                       
                // remove mesh faces from boundary nodes
                
                for (auto& nodeIter : this->nodes())
                {// trigger regular calls in DislocationNode::moveGlide, but with zero motion
                    nodeIter.second->moveGlide(0.0);
                }
            }
            else
            {
                for (auto& nodeIter : this->nodes())
                {
                    nodeIter.second->moveGlide(dt_in);
                }
            }
            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
        }

        // 2021-4-2, lizt
        void  localPlasticStrain(const double & dt_in )
        {
            model::cout<<"        Calculating Plastic Strain Increment by Glide (dt="<<dt_in<<")..."<<std::flush;
            const auto t0= std::chrono::system_clock::now();
            if (simulationParameters.isPeriodicSimulation())
            {
                assert(0 && "periodical case is uncompleted ...");
            }
            else
            {
                // get quadrature position and velocity from nodal velocity; 
                for (const auto& link : this->links())
                {
                    if (  !link.second->hasZeroBurgers() 
                    && !link.second->isBoundarySegment()
                    && !link.second->isVirtualBoundarySegment())
                    {
                        std::vector<VectorDim> quadPosition;
                        std::vector<VectorDim> quadVelocity;
                        VectorDim sourcePosition(link.second->source->get_P());
                        VectorDim sourceVelocity(link.second->source->get_V());
                        VectorDim sinkPosition(link.second->sink->get_P());
                        VectorDim sinkVelocity(link.second->sink->get_V());
                        // const int quadNum(bvpAbaqus->para.quadratureNumber==-1? link.second->quadraturePoints().size() : bvpAbaqus->para.quadratureNumber);
                        const int quadNum(bvpAbaqus->para.quadratureNumber*link.second->quadraturePoints().size());
                        //quadPosition.push_back(link.second->source->get_P());
                        //quadVelocity.push_back(link.second->source->get_V());            
                        //for (const auto& qPoint : link.second->quadraturePoints())
                        //{
                            //qPositionNew.push_back(qPoint.r);
                            //quadVelocity.push_back(qPoint.glideVelocity);
                        //}
                        VectorDim unitPosition=(sinkPosition-sourcePosition)/(quadNum+1);
                        VectorDim unitVelocity=(sinkVelocity-sourceVelocity)/(quadNum+1);
                        const VectorDim faceNormal((sinkPosition-sourcePosition).cross(sourceVelocity+sinkVelocity).normalized());
                        const VectorDim linkb(link.second->burgers().normalized());
                        const VectorDim linkn(faceNormal);
                        //faceNormal.normalize();
                        // model::cout<<"faceNormal = "<<faceNormal.transpose()<<std::endl;
                        for (int i=0;i<quadNum+2;i++)
                        {
                            // linear distribution
                            quadPosition.push_back( sourcePosition+unitPosition*i );
                            quadVelocity.push_back( sourceVelocity+unitVelocity*i );
                        }
                        //model::cout<<"quadrature position size="<<quadPosition.size()<<std::endl;
                        //model::cout<<"quadrature point  size="<<quadNum<<std::endl;
                        //qPositionNew.push_back(link.second->sink->get_P());
                        //quadVelocity.push_back(link.second->sink->get_V());
                        for(int j=0 ; j<quadNum+1 ; j++)
                        {
                            double incSlipArea=0;
                            incSlipArea+=0.5*unitPosition.cross(quadVelocity[j]*dt_in).norm();
                            incSlipArea+=0.5*(unitPosition+dt_in*unitVelocity).cross(quadVelocity[j+1]*dt_in).norm();
                            std::vector<VectorDim> Point;
                            int searchOrder=10;
                            for (int n=0; n<searchOrder; n++)
                            {
                                Point.push_back(quadPosition[j]      + quadVelocity[j]     *dt_in/(searchOrder-1)*n);
                                Point.push_back(quadPosition[j+1] + quadVelocity[j+1]*dt_in/(searchOrder-1)*n);
                                // Point.push_back(quadPosition[j]     +quadVelocity[j]     *dt_in);
                                // Point.push_back(quadPosition[j+1]+quadVelocity[j+1]*dt_in);
                            }
                            // model::cout<<"Point size = "<<Point.size()<<std::endl;
                            // out test
                            // for (auto& temp : Point)
                            // {
                            //     model::cout<<"Point = "<<temp.transpose()<<std::endl;
                            // }
                            // find simplex and element
                            //std::vector<int> elementID;
                            std::map<int,double> ele_plas; // elementID and distance between integral point and slipface
                            double volume(0);
                            int id; // element ID in DD
                            bvpAbaqus->findElementSet(Point,faceNormal);
                            ele_plas.clear();
                            for( const auto& iter : bvpAbaqus->getElementSet().second)
                            {
                                //model::cout<<"integral point position="<<iter->simplex.integralPoint<<std::endl;
                                //elementID.clear();                            

                                id = bvpAbaqus->finiteElement().s2eID.find(iter->simplex.xID)->second;
                                // VectorDim intP2Quad( Point.back() - bvpAbaqus->integralPoints.row(id).transpose() );//vector from integral point to dislocation quadrature point, P->Q
                                // VectorDim subP2Quad(intP2Quad - intP2Quad.dot(faceNormal)*faceNormal);// vector from integral's subpoint to dislocation quadrature point, A->Q
                                // if ( quadVelocity[j+1].dot(subP2Quad)>FLT_EPSILON)
                                // {
                                //     ele_plas.emplace(id,abs((bvpAbaqus->integralPoints.row(id).transpose()-Point[0]).dot(faceNormal)));
                                // }
                                // ele_plas.emplace(id,abs((bvpAbaqus->integralPoints.row(id).transpose()-Point[0]).dot(faceNormal)));
                                ele_plas.emplace(id,(bvpAbaqus->integralPoints.row(id).transpose()-Point[0]).dot(faceNormal));
                                model::cout<<"point0 = "<<Point[0].transpose()<<std::endl;
                                //elementID.push_back();                            
                                //iter->incPlasticStrain=Eigen::Matrix<double,6,1>::Ones();
                            }
                            // model::cout<<"elementID size = "<<ele_plas.size()<<std::endl;
                            volume=bvpAbaqus->getElementSet().first;
                            // assert(abs(volume)>FLT_EPSILON && "characteristic volume = 0 !");
                            assert(volume>FLT_EPSILON && "characteristic volume is Zero or negative !");

                            // total increment of plastic strain, inc_gamma=b*S/V
                            double incGamma;
                            // model::cout<<"SlipArea = "<<incSlipArea<<std::endl;
                            incGamma=incSlipArea/volume;
                            // model::cout<<"gamma increment = "<<incGamma<<std::endl;
                            // 2021-4-30, luosc
                            // const VectorDim linkb(link.second->burgers().normalized());
                            // const VectorDim linkn(faceNormal);
                            // model::cout<<"linkb = "<<linkb.transpose()<<std::endl;
                            // model::cout<<"linkn = "<<linkn.transpose()<<std::endl;
                            Eigen::Matrix<double,1,6> gamma(Eigen::Matrix<double,1,6>::Zero());
                            gamma(0,0) = incGamma*(linkb(0)*linkn(0)+linkb(0)*linkn(0))/2;
                            gamma(0,1) = incGamma*(linkb(1)*linkn(1)+linkb(1)*linkn(1))/2;
                            gamma(0,2) = incGamma*(linkb(2)*linkn(2)+linkb(2)*linkn(2))/2;
                            gamma(0,3) = incGamma*(linkb(0)*linkn(1)+linkb(1)*linkn(0))/2;
                            gamma(0,4) = incGamma*(linkb(1)*linkn(2)+linkb(2)*linkn(1))/2;
                            gamma(0,5) = incGamma*(linkb(0)*linkn(2)+linkb(2)*linkn(0))/2;
                            // model::cout<<"gamma = "<<gamma<<std::endl;
                            // for (const auto& ele : ele_plas)
                            // {
                            //     bvpAbaqus->PlasticStrain().row(ele.first)+=gamma;
                            // }

                            // 2021-4-30, luosc  distribution function
                            if (incGamma > FLT_EPSILON)
                            {
                                double  g_ep=0;
                                for (const auto& ele : ele_plas)
                                {
                                    if(abs(ele.second)<expandheight )
                                    {
                                        //bvpAbaqus->PlasticStrain().row(ele.first) +=d_Gamma(ele.second)*gamma;
                                    
                                        g_ep = g_ep + d_Gamma(ele.second)*incGamma;
                                    }
                                }
                                // to sure the ep is equal to the incGamma                            
                                for (const auto& ele : ele_plas)
                                {
                                    if(abs(ele.second)<expandheight )
                                    // if (g_ep<=incGamma)
                                    {
                                        // model::cout<<"g_ep = "<<g_ep<<std::endl;
                                        double no_enough= incGamma/g_ep;
                                        //assert(abs(g_ep)>FLT_EPSILON && "total plastic strain is Zero ");
                                        // model::cout<<"no_enough = "<<no_enough<<std::endl;
                                        // assert(no_enough<=100.0 && "total plastic strain is Zero ");
                                        // bvpAbaqus->PlasticStrain().row(ele.first) = bvpAbaqus->PlasticStrain().row(ele.first)*no_enough;
                                        bvpAbaqus->PlasticStrain().row(ele.first) += d_Gamma(ele.second)*gamma*no_enough;
                                        // model::cout<<"h = "<<ele.second<<std::endl;

                                        // for data transfer test
                                        // Eigen::Matrix<double,1,6> tempTest;
                                        // tempTest<<0.001,0.001,0.001,0.002,0.003,ele.second;
                                        // bvpAbaqus->PlasticStrain().row(ele.first) = tempTest;
                                        // bvpAbaqus->PlasticStrain().row(ele.first) = ele.second;

                                        // model::cout<<"1"<<d_Gamma(ele.second)*gamma*no_enough<<std::endl;
                                        // model::cout<<"2"<<bvpAbaqus->PlasticStrain().row(ele.first)<<std::endl;
                                        
                                    }
                                }
                            }
                            


                        }
                        // output for test 
                        //model::cout<<"quadrature position size="<<qPositionNew.size()<<std::endl;
                        //model::cout<<"quadrature velocity  size="<<quadVelocity.size()<<std::endl;
                    // calculate total inc_gamma for each slip aera
                    }
                }

            } 
            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
        }



        /*************************************************************************/
        // 2021-5-14, lizt constider velocity reduction
        void  localPlasticStrain_1(const double & dt_in )
        {
            model::cout<<"        Calculating Plastic Strain Increment by Glide (dt="<<dt_in<<")..."<<std::flush;
            const auto t0= std::chrono::system_clock::now();
            if (simulationParameters.isPeriodicSimulation())
            {
                assert(0 && "periodical case is uncompleted ...");
            }
            else
            {
                // get quadrature position and velocity from nodal velocity; 
                for (const auto& link : this->links())
                {
                    if (  !link.second->hasZeroBurgers() 
                    && !link.second->isBoundarySegment()
                    && !link.second->isVirtualBoundarySegment())
                    {
                        std::vector<VectorDim> quadPosition;
                        std::vector<VectorDim> quadVelocity;
                        VectorDim sourcePosition(link.second->source->get_P());
                        VectorDim sourceVelocity(link.second->source->get_V());
                        VectorDim sinkPosition(link.second->sink->get_P());
                        VectorDim sinkVelocity(link.second->sink->get_V());
                        const int quadNum(bvpAbaqus->para.quadratureNumber==-1? link.second->quadraturePoints().size() : bvpAbaqus->para.quadratureNumber);
                        //quadPosition.push_back(link.second->source->get_P());
                        //quadVelocity.push_back(link.second->source->get_V());            
                        //for (const auto& qPoint : link.second->quadraturePoints())
                        //{
                            //qPositionNew.push_back(qPoint.r);
                            //quadVelocity.push_back(qPoint.glideVelocity);
                        //}
                        VectorDim unitPosition=(sinkPosition-sourcePosition)/(quadNum+1);
                        VectorDim unitVelocity=(sinkVelocity-sourceVelocity)/(quadNum+1);
                        const VectorDim faceNormal((sinkPosition-sourcePosition).cross(sourceVelocity+sinkVelocity).normalized());
                        const VectorDim linkb(link.second->burgers().normalized());
                        const VectorDim linkn(faceNormal);
                        //faceNormal.normalize();
                        // model::cout<<"faceNormal = "<<faceNormal.transpose()<<std::endl;
                        for (int i=0;i<quadNum+2;i++)
                        {
                            // linear distribution
                            quadPosition.push_back( sourcePosition+unitPosition*i );
                            quadVelocity.push_back( sourceVelocity+unitVelocity*i );
                        }
                        //model::cout<<"quadrature position size="<<quadPosition.size()<<std::endl;
                        //model::cout<<"quadrature point  size="<<quadNum<<std::endl;
                        //qPositionNew.push_back(link.second->sink->get_P());
                        //quadVelocity.push_back(link.second->sink->get_V());
                        for(int j=0 ; j<quadNum+1 ; j++)
                        {
                            double incSlipArea=0;
                            incSlipArea+=0.5*unitPosition.cross(quadVelocity[j]*dt_in).norm();
                            incSlipArea+=0.5*(unitPosition-dt_in*unitVelocity).cross(quadVelocity[j+1]*dt_in).norm();
                            std::vector<VectorDim> Point;
                            const int searchOrder(bvpAbaqus->para.searchOrder);
                            for (int n=0; n<searchOrder; n++)
                            {
                                Point.push_back(quadPosition[j]      - quadVelocity[j]     *dt_in/(searchOrder-1)*n);
                                Point.push_back(quadPosition[j+1] - quadVelocity[j+1]*dt_in/(searchOrder-1)*n);
                                // Point.push_back(quadPosition[j]     +quadVelocity[j]     *dt_in);
                                // Point.push_back(quadPosition[j+1]+quadVelocity[j+1]*dt_in);
                            }
                            // model::cout<<"Point size = "<<Point.size()<<std::endl;
                            // out test
                            // for (auto& temp : Point)
                            // {
                            //     model::cout<<"Point = "<<temp.transpose()<<std::endl;
                            // }
                            // find simplex and element
                            //std::vector<int> elementID;
                            std::map<int,double> ele_plas; // elementID and distance between integral point and slipface
                            double volume(0);
                            int id; // element ID in DD
                            bvpAbaqus->findElementSet(Point,faceNormal);
                            ele_plas.clear();
                            for( const auto& iter : bvpAbaqus->getElementSet().second)
                            {
                                //model::cout<<"integral point position="<<iter->simplex.integralPoint<<std::endl;
                                //elementID.clear();                            

                                id = bvpAbaqus->finiteElement().s2eID.find(iter->simplex.xID)->second;
                                // VectorDim intP2Quad( Point.back() - bvpAbaqus->integralPoints.row(id).transpose() );//vector from integral point to dislocation quadrature point, P->Q
                                // VectorDim subP2Quad(intP2Quad - intP2Quad.dot(faceNormal)*faceNormal);// vector from integral's subpoint to dislocation quadrature point, A->Q
                                // if ( quadVelocity[j+1].dot(subP2Quad)>FLT_EPSILON)
                                // {
                                //     ele_plas.emplace(id,abs((bvpAbaqus->integralPoints.row(id).transpose()-Point[0]).dot(faceNormal)));
                                // }
                                // ele_plas.emplace(id,abs((bvpAbaqus->integralPoints.row(id).transpose()-Point[0]).dot(faceNormal)));
                                ele_plas.emplace(id,(bvpAbaqus->integralPoints.row(id).transpose()-Point[0]).dot(faceNormal));
                                // model::cout<<"point0 = "<<Point[0].transpose()<<std::endl;
                                //elementID.push_back();                            
                                //iter->incPlasticStrain=Eigen::Matrix<double,6,1>::Ones();
                            }
                            // model::cout<<"elementID size = "<<ele_plas.size()<<std::endl;
                            volume=bvpAbaqus->getElementSet().first;
                            // assert(abs(volume)>FLT_EPSILON && "characteristic volume = 0 !");
                            assert(volume>FLT_EPSILON && "characteristic volume is Zero or negative !");

                            // total increment of plastic strain, inc_gamma=b*S/V
                            double incGamma;
                            // model::cout<<"SlipArea = "<<incSlipArea<<std::endl;
                            incGamma=incSlipArea/volume;
                            // model::cout<<"gamma increment = "<<incGamma<<std::endl;
                            // 2021-4-30, luosc
                            // const VectorDim linkb(link.second->burgers().normalized());
                            // const VectorDim linkn(faceNormal);
                            // model::cout<<"linkb = "<<linkb.transpose()<<std::endl;
                            // model::cout<<"linkn = "<<linkn.transpose()<<std::endl;
                            Eigen::Matrix<double,1,6> gamma(Eigen::Matrix<double,1,6>::Zero());
                            gamma(0,0) = incGamma*(linkb(0)*linkn(0)+linkb(0)*linkn(0))/2;
                            gamma(0,1) = incGamma*(linkb(1)*linkn(1)+linkb(1)*linkn(1))/2;
                            gamma(0,2) = incGamma*(linkb(2)*linkn(2)+linkb(2)*linkn(2))/2;
                            gamma(0,3) = incGamma*(linkb(0)*linkn(1)+linkb(1)*linkn(0))/2;
                            gamma(0,4) = incGamma*(linkb(1)*linkn(2)+linkb(2)*linkn(1))/2;
                            gamma(0,5) = incGamma*(linkb(0)*linkn(2)+linkb(2)*linkn(0))/2;
                            // model::cout<<"gamma = "<<gamma<<std::endl;
                            // for (const auto& ele : ele_plas)
                            // {
                            //     bvpAbaqus->PlasticStrain().row(ele.first)+=gamma;
                            // }

                            // 2021-4-30, luosc  distribution function
                            if (incGamma > FLT_EPSILON)
                            {
                                double  g_ep1=0;
                                double  g_ep2=0;
                                for (const auto& ele : ele_plas)
                                {
                                    if(ele.second<= expandheight && ele.second>=0.0)
                                    {
                                        //bvpAbaqus->PlasticStrain().row(ele.first) +=d_Gamma(ele.second)*gamma;                                    
                                        g_ep1 = g_ep1 + d_Gamma(ele.second);
                                    }
                                    if(ele.second>=-expandheight && ele.second<=0.0)
                                    {
                                        //bvpAbaqus->PlasticStrain().row(ele.first) +=d_Gamma(ele.second)*gamma;                                    
                                        g_ep2 = g_ep2 + d_Gamma(ele.second);
                                    }
                                }
                                double no_enough1= 0.5/g_ep1;
                                double no_enough2= 0.5/g_ep2;
                                // to sure the ep is equal to the incGamma                            
                                for (const auto& ele : ele_plas)
                                {
                                    if(ele.second<=expandheight && ele.second>=0.0 )
                                    // if (g_ep<=incGamma)
                                    {
                                        // model::cout<<"g_ep = "<<g_ep<<std::endl;
                                        
                                        //assert(abs(g_ep)>FLT_EPSILON && "total plastic strain is Zero ");
                                        // model::cout<<"no_enough = "<<no_enough<<std::endl;
                                        // assert(no_enough<=100.0 && "total plastic strain is Zero ");
                                        // bvpAbaqus->PlasticStrain().row(ele.first) = bvpAbaqus->PlasticStrain().row(ele.first)*no_enough;
                                        bvpAbaqus->PlasticStrain().row(ele.first) += d_Gamma(ele.second)*gamma*no_enough1;
                                        // model::cout<<"h = "<<ele.second<<std::endl;

                                        // for data transfer test
                                        // Eigen::Matrix<double,1,6> tempTest;
                                        // tempTest<<0.001,0.001,0.001,0.002,0.003,ele.second;
                                        // bvpAbaqus->PlasticStrain().row(ele.first) = tempTest;
                                        // bvpAbaqus->PlasticStrain().row(ele.first) = ele.second;

                                        // model::cout<<"1"<<d_Gamma(ele.second)*gamma*no_enough<<std::endl;
                                        // model::cout<<"2"<<bvpAbaqus->PlasticStrain().row(ele.first)<<std::endl;
                                        
                                    }
                                    else if(ele.second>-expandheight && ele.second<=0.0 )
                                    {
                                        
                                        bvpAbaqus->PlasticStrain().row(ele.first) += d_Gamma(ele.second)*gamma*no_enough2;
                                    }
                                }
                            }

                        }
                        // output for test 
                        //model::cout<<"quadrature position size="<<qPositionNew.size()<<std::endl;
                        //model::cout<<"quadrature velocity  size="<<quadVelocity.size()<<std::endl;
                    // calculate total inc_gamma for each slip aera
                    }
                }

            } 
            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
        }

        /**********************************************************************/
        // 2021-5-14, lizt  another method to localize plastic strain
        // 2022-11-16, update, for C3D8 element
       void  localPlasticStrain_2(const double & dt_in )
        {
            model::cout<<"        Calculating Plastic Strain Increment by Glide (dt="<<dt_in<<")..."<<std::flush;
            const auto t0= std::chrono::system_clock::now();
            if (simulationParameters.isPeriodicSimulation())
            {
                assert(0 && "periodical case is uncompleted ...");
            }
            else
            {
                auto plasticStrain = bvpAbaqus->PlasticStrain();
                // get quadrature position and velocity from nodal velocity; 
#ifdef _OPENMP  // parallelize loop over DislocationLinks
            const size_t nThreads = omp_get_max_threads();

#pragma omp parallel for reduction(+:totalS) shared(plasticStrain)
                for (unsigned int k=0;k<this->links().size();++k)
                {
                    auto snIter(this->links().begin());
                    std::advance(snIter,k);
                    auto link(*snIter);
                    if (  !link.second->hasZeroBurgers() 
                    // && !link.second->isBoundarySegment()
                    && !link.second->isVirtualBoundarySegment())
                    {
                        std::vector<VectorDim> quadPosition;
                        std::vector<VectorDim> quadVelocity;
                        VectorDim sourcePosition(link.second->source->get_P());
                        VectorDim sourceVelocity(link.second->source->get_V());
                        VectorDim sinkPosition(link.second->sink->get_P());
                        VectorDim sinkVelocity(link.second->sink->get_V());
                        const int quadNum(bvpAbaqus->para.quadratureNumber* link.second->quadraturePoints().size());

                        VectorDim unitPosition=(sinkPosition-sourcePosition)/(quadNum+1);
                        VectorDim unitVelocity=(sinkVelocity-sourceVelocity)/(quadNum+1);
                        const VectorDim faceNormal((sinkPosition-sourcePosition).cross(sourceVelocity+sinkVelocity).normalized());
                        const VectorDim linkb(link.second->burgers().normalized());
                        const VectorDim linkn(faceNormal);

                        for (int i=0;i<quadNum+2;i++)
                        {
                            // linear distribution
                            quadPosition.push_back( sourcePosition+unitPosition*i );
                            quadVelocity.push_back( sourceVelocity+unitVelocity*i );
                        }

                        for(int j=0 ; j<quadNum+1 ; j++)
                        {
                            double incSlipArea=0;
                            incSlipArea+=0.5*unitPosition.cross(quadVelocity[j]*dt_in).norm();
                            incSlipArea+=0.5*(unitPosition-dt_in*unitVelocity).cross(quadVelocity[j+1]*dt_in).norm();
                            totalS+=incSlipArea;
                            // model::cout<<"slip aera increament = "<<incSlipArea<<std::endl;
                            std::vector<VectorDim> Point;
                            const int searchOrder(bvpAbaqus->para.searchOrder);
                            for (int n=0; n<searchOrder; n++)
                            {
                                Point.push_back(quadPosition[j]      - quadVelocity[j]     *dt_in/(searchOrder-1)*n);
                                Point.push_back(quadPosition[j+1] - quadVelocity[j+1]*dt_in/(searchOrder-1)*n);

                            }

                            // find simplex and element
                            //std::vector<int> elementID;
                            std::map<int,double> ele_plas(bvpAbaqus->findElementSetAndWeight(Point,faceNormal)); // elementID and plastic strain weight
                            double volume(0);
                            // bvpAbaqus->findElementSet(Point,faceNormal);
                            volume=bvpAbaqus->volumeC*6.0;
                            assert(volume>FLT_EPSILON && "characteristic volume is Zero or negative !");

                            // total increment of plastic strain, inc_gamma=b*S/V
                            double incGamma;
                            incGamma=incSlipArea/volume;

                            Eigen::Matrix<double,1,6> gamma(Eigen::Matrix<double,1,6>::Zero());
                            gamma(0,0) = incGamma*(linkb(0)*linkn(0)+linkb(0)*linkn(0))/2;
                            gamma(0,1) = incGamma*(linkb(1)*linkn(1)+linkb(1)*linkn(1))/2;
                            gamma(0,2) = incGamma*(linkb(2)*linkn(2)+linkb(2)*linkn(2))/2;
                            gamma(0,3) = incGamma*(linkb(0)*linkn(1)+linkb(1)*linkn(0))/2;
                            gamma(0,4) = incGamma*(linkb(1)*linkn(2)+linkb(2)*linkn(1))/2;
                            gamma(0,5) = incGamma*(linkb(0)*linkn(2)+linkb(2)*linkn(0))/2;

                            for (const auto& ele : ele_plas)
                            {
                                plasticStrain.row(ele.first-1) += gamma*ele.second;
                                // model::cout<<"weight = "<<ele.second<<std::endl;
                            }
                        }
                    }
                }
#else
                for (const auto& link : this->links())
                {
                    if (  !link.second->hasZeroBurgers() 
                    // && !link.second->isBoundarySegment()
                    && !link.second->isVirtualBoundarySegment())
                    {
                        std::vector<VectorDim> quadPosition;
                        std::vector<VectorDim> quadVelocity;
                        VectorDim sourcePosition(link.second->source->get_P());
                        VectorDim sourceVelocity(link.second->source->get_V());
                        VectorDim sinkPosition(link.second->sink->get_P());
                        VectorDim sinkVelocity(link.second->sink->get_V());
                        const int quadNum(bvpAbaqus->para.quadratureNumber* link.second->quadraturePoints().size());

                        VectorDim unitPosition=(sinkPosition-sourcePosition)/(quadNum+1);
                        VectorDim unitVelocity=(sinkVelocity-sourceVelocity)/(quadNum+1);
                        const VectorDim faceNormal((sinkPosition-sourcePosition).cross(sourceVelocity+sinkVelocity).normalized());
                        const VectorDim linkb(link.second->burgers().normalized());
                        const VectorDim linkn(faceNormal);

                        for (int i=0;i<quadNum+2;i++)
                        {
                            // linear distribution
                            quadPosition.push_back( sourcePosition+unitPosition*i );
                            quadVelocity.push_back( sourceVelocity+unitVelocity*i );
                        }

                        for(int j=0 ; j<quadNum+1 ; j++)
                        {
                            double incSlipArea=0;
                            incSlipArea+=0.5*unitPosition.cross(quadVelocity[j]*dt_in).norm();
                            incSlipArea+=0.5*(unitPosition-dt_in*unitVelocity).cross(quadVelocity[j+1]*dt_in).norm();
                            totalS+=incSlipArea;
                            // model::cout<<"slip aera increament = "<<incSlipArea<<std::endl;
                            std::vector<VectorDim> Point;
                            const int searchOrder(bvpAbaqus->para.searchOrder);
                            for (int n=0; n<searchOrder; n++)
                            {
                                Point.push_back(quadPosition[j]      - quadVelocity[j]     *dt_in/(searchOrder-1)*n);
                                Point.push_back(quadPosition[j+1] - quadVelocity[j+1]*dt_in/(searchOrder-1)*n);

                            }

                            // find simplex and element
                            //std::vector<int> elementID;
                            std::map<int,double> ele_plas(bvpAbaqus->findElementSetAndWeight(Point,faceNormal)); // elementID and plastic strain weight
                            // double total;
                            // for (const auto& ele : ele_plas)
                            // {
                            //     total+=ele.second;
                            //     // model::cout<<"weight = "<<ele.second<<std::endl;
                            // }
                            // model::cout<<"total  = "<<total<<std::endl;
                            // for(const auto& iter : ele_plas)
                            // {//  test
                            //     std::cout<<iter.first<<"  "<<iter.second<<std::endl;
                            // }  
                            double volume(0);
                            // bvpAbaqus->findElementSet(Point,faceNormal);
                            volume=bvpAbaqus->volumeC*6.0;
                            assert(volume>FLT_EPSILON && "characteristic volume is Zero or negative !");

                            // total increment of plastic strain, inc_gamma=b*S/V
                            double incGamma;
                            incGamma=incSlipArea/volume;
                            // model::cout<<"inc gamma = "<<incGamma<<std::endl;

                            Eigen::Matrix<double,1,6> gamma(Eigen::Matrix<double,1,6>::Zero());
                            gamma(0,0) = incGamma*(linkb(0)*linkn(0)+linkb(0)*linkn(0))/2;
                            gamma(0,1) = incGamma*(linkb(1)*linkn(1)+linkb(1)*linkn(1))/2;
                            gamma(0,2) = incGamma*(linkb(2)*linkn(2)+linkb(2)*linkn(2))/2;
                            gamma(0,3) = incGamma*(linkb(0)*linkn(1)+linkb(1)*linkn(0))/2;
                            gamma(0,4) = incGamma*(linkb(1)*linkn(2)+linkb(2)*linkn(1))/2;
                            gamma(0,5) = incGamma*(linkb(0)*linkn(2)+linkb(2)*linkn(0))/2;

                            for (const auto& ele : ele_plas)
                            {
                                plasticStrain.row(ele.first-1) += gamma*ele.second;
                                // bvpAbaqus->PlasticStrain().row(ele.first) += gamma*ele.second;
                                // model::cout<<"weight = "<<ele.second<<std::endl;
                            }
                            // for (size_t i = 0; i <bvpAbaqus-> PlasticStrain().rows(); i++)
                            // {
                            //     model::cout<<"plastic strain = "<<bvpAbaqus->PlasticStrain().row(i)<<std::endl;
                            // }
                            
                        }
                    }
                }
#endif
                bvpAbaqus->PlasticStrain() = plasticStrain;
            } 

            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
        }


        /**********************************************************************/
        //luosc g(x) =(x/epsilon)/(1+x/epsilon)+arctan(x/epsilon)
        double g_function(const double& ele_x)
       {
           double elex(ele_x/epsilon);
           double  g = elex/(1+elex*elex)+atan(elex);

           return g;

       }
        /**********************************************************************/
        //luosc delta_Gamma
         double d_Gamma(const double& ele_R)
       {
           double  dG = (g_function(ele_R+ele_h/2)-g_function(ele_R-ele_h/2))/(2*g_function(nc*ele_h));

           return dG;

       }

        // 2023-1-3
        /**********************************************************************/
        MatrixDimD equivlantStress(VoigtStress stressV)
        {
            MatrixDimD stress;
            stress<<stressV(0),stressV(3),stressV(5),
                    stressV(3),stressV(1),stressV(4),
                    stressV(5),stressV(4),stressV(2);

            return stress;
        }


        /**********************************************************************/
        void dislocationNucleation(const double& dt)
        {
            double kB_SI=1.38064852e-23;     // Boltzmann constant in SI units [J/K]
            double kB_eV=8.617e-5;           // Boltzmann constant in [eV]
            double Tm_K=poly.Tm;             // 1141;
            double sigma_th=poly.sigma_th/poly.mu0_SI;      // code units
            double v0=1.0e13;
            double A_eV=4.8;
            double alpha=poly.alpha;
            double vvSt=0.0;
            double T_K=poly.T;
            double dt_SI=dt*poly.b_SI/poly.cs_SI;
            double shearLoopSize=simulationParameters.shearLoopSize;
            double nucleationSiteDistance_s=simulationParameters.nucleationSiteDistance_s;
            double nucleationSiteDistance_n=simulationParameters.nucleationSiteDistance_n;
            double targetShearLoopDensity=simulationParameters.targetShearLoopDensity;
            double targetShearLoopDensity_1=simulationParameters.targetShearLoopDensity_1;
            double targetFRDensity=simulationParameters.targetFRDensity;
            double shearLoopDensity = 0.0;
            double shearLoopDensity_1 = 0.0;
            double FRDensity = 0.0;
            double sigma_c = stress_nuc/poly.mu0_SI;

           

            // find all possible nucleation sites & slipSystems
            // stochastic sites to dislocation nucleation

            // std::mt19937 generator(std::chrono::system_clock::now().time_since_epoch().count());
            // const int rSS=distribution(generator); // a random SlipSystem ID
            // const int rSS=19; // a specific SlipSystem ID
            const int grainID=1;
            // const SlipSystem& slipSystem(*poly.grain(grainID).slipSystems()[rSS]);
            // const VectorDim b=slipSystem.s.cartesian();     // Burgers vector
            // const VectorDim n=slipSystem.unitNormal;        // slip plane normal

            // std::uniform_real_distribution<> dis(0.0, 0.5*M_PI);
            // const double theta=dis(generator); // random angle of the dislocation line in the plane from screw orientation.

            // bvpAbaqus->IDmap.initNodes();
            std::map<int,std::pair<int,VectorDim>> nucleationSite;    // map< orderID, < slipSystemID, sitePosition >>
           //  if (!((simulationParameters.runID)%outputFrequency))
           //  {
           //  NucSite.clear();
          //   }
            int ID = 0;
// To do: parallel
            for (size_t i = 0; i < bvpAbaqus->stressAbaqus.rows(); i++)
            {
                double tau = 0.0;
                int slipSystemID;
                for (size_t rSS = 0; rSS < poly.grain(grainID).slipSystems().size(); rSS++)
                {
                    const SlipSystem& slipSystem(*poly.grain(grainID).slipSystems()[rSS]);
                    const VectorDim b=slipSystem.s.cartesian();     // Burgers vector
                    const VectorDim n=slipSystem.unitNormal;        // slip plane normal                    
                    double sigma=(equivlantStress(bvpAbaqus->stressAbaqus.row(i))*n).transpose()*b;
                    if (sigma > tau)
                    {
                        tau = sigma;
                        slipSystemID = rSS;
                    }
                }                
                // double sigma=0.01;
                double Q=(1.0-T_K/Tm_K)*A_eV*pow(1.0-tau/sigma_th, alpha);
                double vv=v0*exp(-Q/kB_eV/T_K);
                double vvS=vv;          // single nucleation site
                // double vvS=vv*1.005e9;//Fe
                double vvSt=0.0;
                vvSt=vvSt+vvS*dt_SI;
                // vvSt=vvSt+vvS*dt*0.01;//1e-10-5
                // std::cout<<"i = "<<i<<std::endl;
                // std::cout<<"Q = "<<Q<<std::endl;
                // std::cout<<"vv = "<<vv<<std::endl;
                // std::cout<<"vvS = "<<vvS<<std::endl;
                // std::cout<<"vvSt = "<<vvSt<<std::endl;
                // if(vvSt>1.0)
                if(tau>sigma_c)
                {
                    std::cout<<"sigma = "<<tau<<std::endl;
                    // VectorDim P_n(50.0*i,50.0*i,0.0);
                    VectorDim P_n(bvpAbaqus->IDmap.getNodes()[bvpAbaqus->IDmap.NodeMap()[i+1]]);
                    nucleationSite.emplace(ID,std::make_pair(slipSystemID,P_n));
                    ++ID;
                    vvSt=vvSt-1.0;
                    std::cout<<"vvSt = "<<vvSt<<std::endl;
                }
            }
            if(ID>=1)
            {
                std::mt19937 generator(std::chrono::system_clock::now().time_since_epoch().count());
                std::uniform_int_distribution<> distribution(0,ID-1);
                while(shearLoopDensity<targetShearLoopDensity && shearLoopDensity<ID)
                {
                    int siteID = distribution(generator);       // first nucleation site, random
                    auto firstSite = nucleationSite[siteID];
                    auto ssID0 = firstSite.first;
                    auto P_0 = firstSite.second;
                    NucSite.push_back(P_0);
                    addShearLoopWhileLoad(shearLoopSize,ssID0,P_0);
                    shearLoopDensity=shearLoopDensity+1.0;
                }
                while(shearLoopDensity_1<targetShearLoopDensity_1 && shearLoopDensity_1<ID)
                {
                    int siteID = distribution(generator);       // first nucleation site, random
                    auto firstSite = nucleationSite[siteID];
                    auto ssID0 = firstSite.first;
                    auto P_0 = firstSite.second;
                    NucSite.push_back(P_0);
                    addShearLoopWhileLoad_1(shearLoopSize,ssID0,P_0);
                    shearLoopDensity_1=shearLoopDensity_1+1.0;
                }
                while(FRDensity<targetFRDensity && FRDensity<ID)
                {
                    int siteID = distribution(generator);       // first nucleation site, random
                    auto firstSite = nucleationSite[siteID];
                    auto ssID0 = firstSite.first;
                    auto P_0 = firstSite.second;
                    NucSite.push_back(P_0);
                    addFRLoad(shearLoopSize,ssID0,P_0);
                    FRDensity=FRDensity+1.0;
                }
                // std::mt19937 generator(std::chrono::system_clock::now().time_since_epoch().count());
                // std::uniform_int_distribution<> distribution(0,ID-1);
                // int siteID = distribution(generator);       // first nucleation site, random
                // auto firstSite = nucleationSite[siteID];
                // auto ssID0 = firstSite.first;
                // auto P_0 = firstSite.second;
                // addShearLoopWhileLoad(shearLoopSize,ssID0,P_0);

                // for (const auto& site : nucleationSite)
                // {
                //     auto ssID(site.second.first);           // slipSystem ID
                //     auto P_s(site.second.second);           // site position
                //     double distance_n;                      // distance on n-direction
                //     double distance_s;                      // distance in slip plane
                //     auto n = poly.grain(grainID).slipSystems()[ssID]->unitNormal;
                //     distance_n = ((P_s - P_0).transpose()*n).norm();
                //     distance_s = ((P_s - P_0) -distance_n*n).norm();
                //     if (distance_n >= nucleationSiteDistance_n)
                //     {
                //         addShearLoopWhileLoad(shearLoopSize,ssID,P_s);
                //     }
                //     else if( distance_s > nucleationSiteDistance_s )
                //     {
                //         addShearLoopWhileLoad(shearLoopSize,ssID,P_s);
                //     }
                // }
            }            

            //double    sigma= 1.0/86.0;
            //double prob=vv*dt*7.506e-14;
            //double vvSt=vvSt+vvS*dt;//1e-8stressRate//1e-10-2
            //double vvSt=vvSt+vvS*dt*0.01;//1e-10stressRate//not good
            //std::cout<<"externalLoadController->externalStressRate().trace() = "<<externalLoadController->externalStressRate().trace()<<std::endl;
            //double vvSt=vvSt+vvS*dt*externalLoadController->externalStressRate().trace()/1.0e-8;//1e-10-3
            //double vvSt=vvSt+vvS*dt*0.1;//1e-10-4
            // std::cout<<"Q = "<<Q<<std::endl;
            // std::cout<<"vv = "<<vv<<std::endl;
            // std::cout<<"vvS = "<<vvS<<std::endl;
            // std::cout<<"vvSt = "<<vvSt<<std::endl;

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
//            if (!((simulationParameters.runID)%outputFrequency))
//          {
//            NucSite.clear();
//          }
        }

        /**********************************************************************/
        void addShearLoopWhileLoad( double shearLoopSize,
                                    const int rSS,
                                    VectorDim P_n)
        {
            std::cout<<"        addding shear loop~~~"<<std::endl;
            
            size_t nodeID(StaticID<NodeType>::nextID());
            size_t snID(0);
            size_t loopID(StaticID<LoopType>::nextID()); 
            
            for(const auto & component : this->components())
            {
                if(component.first>snID)
                {
                    snID=component.first;
                }
            }
            snID++;          

            size_t nodeID0(nodeID);
            size_t snID0(snID);
            size_t loopID0(loopID);
            DDconfigIO<3> configIO;            

            std::cout<<magentaBoldColor<<"Generating shear loops"<<defaultColor<<std::endl;
            // VectorDim Pn(0.0,0.0,0.0);
            const std::pair<LatticeVector<dim>,int> rp=poly.specificLatticePointInMesh(P_n);
            const int& grainID=rp.second;           // random grain ID
            const LatticeVector<dim>& L0=rp.first;  // random lattice position in the grain
            const VectorDim P0(L0.cartesian());     // cartesian position of L0
            std::mt19937 generator(std::chrono::system_clock::now().time_since_epoch().count());
            // std::uniform_int_distribution<> distribution(0,poly.grain(grainID).slipSystems().size()-1);
            // const int rSS=distribution(generator); // a random SlipSystem ID
            const SlipSystem& slipSystem(*poly.grain(grainID).slipSystems()[rSS]);
            const VectorDim b=slipSystem.s.cartesian()*(-1);    // Burgers vector
            const VectorDim n=slipSystem.unitNormal; // slip plane normal
            std::uniform_real_distribution<> dis(0.0, 0.5*M_PI);
            const double theta=dis(generator); // random angle of the dislocation line in the plane from screw orientation.
            VectorDim d1=Eigen::AngleAxisd(theta, n)*b.normalized();
            //std::cout<<"d1.norm() = "<<d1.norm()<<std::endl;
            VectorDim d2=Eigen::AngleAxisd(theta+0.5*M_PI, n)*b.normalized();
            VectorDim P1;
            VectorDim P2;
            VectorDim P3;
            VectorDim P4;
            P1[0]=P0[0]+shearLoopSize*d1[0];
            P1[1]=P0[1]+shearLoopSize*d1[1];
            P1[2]=P0[2]+shearLoopSize*d1[2];
            P2[0]=P0[0]+shearLoopSize*d2[0];
            P2[1]=P0[1]+shearLoopSize*d2[1];
            P2[2]=P0[2]+shearLoopSize*d2[2];
            P3[0]=P0[0]-shearLoopSize*d1[0];
            P3[1]=P0[1]-shearLoopSize*d1[1];
            P3[2]=P0[2]-shearLoopSize*d1[2];
            P4[0]=P0[0]-shearLoopSize*d2[0];
            P4[1]=P0[1]-shearLoopSize*d2[1];
            P4[2]=P0[2]-shearLoopSize*d2[2];
            const auto search1(mesh.search(P1));
            const auto search2(mesh.search(P2));
            const auto search3(mesh.search(P3));
            const auto search4(mesh.search(P4));

            if(    search1.first && search1.second->region->regionID==grainID
                && search2.first && search2.second->region->regionID==grainID
                && search3.first && search3.second->region->regionID==grainID
                && search4.first && search4.second->region->regionID==grainID)
            {
                std::vector<VectorDim> nodePosShear;
                nodePosShear.push_back(P1);
                nodePosShear.push_back(P2);
                nodePosShear.push_back(P3);
                nodePosShear.push_back(P4);
                for(size_t k=0;k<nodePosShear.size();++k)
                {
                    configIO.nodes().emplace_back(nodeID+k,nodePosShear[k],Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                    const int nextNodeID=(k+1)<nodePosShear.size()? nodeID+k+1 : nodeID;
                    configIO.links().emplace_back(loopID,nodeID+k,nextNodeID,0);
                }
                nodeID+=nodePosShear.size();

                // write loop file
                configIO.loops().emplace_back(loopID+0, b,n,P1,grainID,DislocationLoopIO<dim>::GLISSILELOOP,-1,VectorDim::Zero());
                loopID+=1;
                snID+=1;
                std::cout<<" shear loop +1, rSS =  "<<rSS<<std::endl;

                std::map<size_t,std::shared_ptr<NodeType>> tempNodes;
                //std::cout<<"in A"<<std::endl;
                const size_t& initNodeSize(nodeID0);
                const size_t& initLoopSize(loopID0);            
                createVertices2(configIO,tempNodes,initNodeSize);
                //std::cout<<"in B"<<std::endl;            
                createEdges2(configIO,tempNodes,initNodeSize,initLoopSize);
            }
        }
        /**********************************************************************/
        void addShearLoopWhileLoad_1( double shearLoopSize,
                                    const int rSS,
                                    VectorDim P_n)
        {
            std::cout<<"        addding shear loop~~~"<<std::endl;
            
            size_t nodeID(StaticID<NodeType>::nextID());
            size_t snID(0);
            size_t loopID(StaticID<LoopType>::nextID()); 
            
            for(const auto & component : this->components())
            {
                if(component.first>snID)
                {
                    snID=component.first;
                }
            }
            snID++;          

            size_t nodeID0(nodeID);
            size_t snID0(snID);
            size_t loopID0(loopID);
            DDconfigIO<3> configIO;            

            std::cout<<magentaBoldColor<<"Generating shear loops"<<defaultColor<<std::endl;
            // VectorDim Pn(0.0,0.0,0.0);
            const std::pair<LatticeVector<dim>,int> rp=poly.specificLatticePointInMesh(P_n);
            const int& grainID=rp.second;           // random grain ID
            const LatticeVector<dim>& L0=rp.first;  // random lattice position in the grain
            const VectorDim P0(L0.cartesian());     // cartesian position of L0
            std::mt19937 generator(std::chrono::system_clock::now().time_since_epoch().count());
            // std::uniform_int_distribution<> distribution(0,poly.grain(grainID).slipSystems().size()-1);
            // const int rSS=distribution(generator); // a random SlipSystem ID
            const SlipSystem& slipSystem(*poly.grain(grainID).slipSystems()[rSS]);
            const VectorDim b=slipSystem.s.cartesian()*(1);    // Burgers vector
            const VectorDim n=slipSystem.unitNormal; // slip plane normal
            std::uniform_real_distribution<> dis(0.0, 0.5*M_PI);
            const double theta=dis(generator); // random angle of the dislocation line in the plane from screw orientation.
            VectorDim d1=Eigen::AngleAxisd(theta, n)*b.normalized();
            //std::cout<<"d1.norm() = "<<d1.norm()<<std::endl;
            VectorDim d2=Eigen::AngleAxisd(theta+0.5*M_PI, n)*b.normalized();
            VectorDim P1;
            VectorDim P2;
            VectorDim P3;
            VectorDim P4;
            P1[0]=P0[0]+shearLoopSize*d1[0];
            P1[1]=P0[1]+shearLoopSize*d1[1];
            P1[2]=P0[2]+shearLoopSize*d1[2];
            P2[0]=P0[0]+shearLoopSize*d2[0];
            P2[1]=P0[1]+shearLoopSize*d2[1];
            P2[2]=P0[2]+shearLoopSize*d2[2];
            P3[0]=P0[0]-shearLoopSize*d1[0];
            P3[1]=P0[1]-shearLoopSize*d1[1];
            P3[2]=P0[2]-shearLoopSize*d1[2];
            P4[0]=P0[0]-shearLoopSize*d2[0];
            P4[1]=P0[1]-shearLoopSize*d2[1];
            P4[2]=P0[2]-shearLoopSize*d2[2];
            const auto search1(mesh.search(P1));
            const auto search2(mesh.search(P2));
            const auto search3(mesh.search(P3));
            const auto search4(mesh.search(P4));

            if(    search1.first && search1.second->region->regionID==grainID
                && search2.first && search2.second->region->regionID==grainID
                && search3.first && search3.second->region->regionID==grainID
                && search4.first && search4.second->region->regionID==grainID)
            {
                std::vector<VectorDim> nodePosShear;
                nodePosShear.push_back(P1);
                nodePosShear.push_back(P2);
                nodePosShear.push_back(P3);
                nodePosShear.push_back(P4);
                for(size_t k=0;k<nodePosShear.size();++k)
                {
                    configIO.nodes().emplace_back(nodeID+k,nodePosShear[k],Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                    const int nextNodeID=(k+1)<nodePosShear.size()? nodeID+k+1 : nodeID;
                    configIO.links().emplace_back(loopID,nodeID+k,nextNodeID,0);
                }
                nodeID+=nodePosShear.size();

                // write loop file
                configIO.loops().emplace_back(loopID+0, b,n,P1,grainID,DislocationLoopIO<dim>::GLISSILELOOP,-1,VectorDim::Zero());
                loopID+=1;
                snID+=1;
                std::cout<<" shear loop +1, rSS =  "<<rSS<<std::endl;

                std::map<size_t,std::shared_ptr<NodeType>> tempNodes;
                //std::cout<<"in A"<<std::endl;
                const size_t& initNodeSize(nodeID0);
                const size_t& initLoopSize(loopID0);            
                createVertices2(configIO,tempNodes,initNodeSize);
                //std::cout<<"in B"<<std::endl;            
                createEdges2(configIO,tempNodes,initNodeSize,initLoopSize);
            }
        }
        /**********************************************************************/
        void addFRLoad( double shearLoopSize,
                                    const int rSS,
                                    VectorDim P_n)
        {
           
            
            size_t nodeID(StaticID<NodeType>::nextID());
            size_t snID(0);
            size_t loopID(StaticID<LoopType>::nextID()); 
            
            for(const auto & component : this->components())
            {
                if(component.first>snID)
                {
                    snID=component.first;
                }
            }
            snID++;          

            size_t nodeID0(nodeID);
            size_t snID0(snID);
            size_t loopID0(loopID);
            DDconfigIO<3> configIO;            
            std::mt19937 generator(std::chrono::system_clock::now().time_since_epoch().count());
                double edgeDensity=0.0;

                const double fractionEdge=1.0; // TEMPORARY
                std::normal_distribution<double> sizeDistribution(FrankReadSizeMeanLoad/poly.b_SI,FrankReadSizeStdLoad/poly.b_SI);
                std::normal_distribution<double> aspectRatioDistribution(FrankReadAspectRatioMeanLoad,FrankReadAspectRatioStdLoad);

                  std::cout<<magentaBoldColor<<"Generating FR source"<<defaultColor<<std::endl;
                    const std::pair<LatticeVector<dim>,int> rp=poly.specificLatticePointInMesh(P_n);
                   // if (std::abs(P_n(1))<200&&std::abs(P_n(2))<200)
                    
                    const LatticeVector<dim> L0=rp.first;
                    const int grainID=rp.second;
                    
                    std::uniform_int_distribution<> distribution(16,poly.grain(grainID).slipSystems().size()-1);
                    // std::uniform_int_distribution<> distribution(0,poly.grain(grainID).slipSystems().size()-1);

                   // const int rSS=distribution(generator); // a random SlipSystem
                    //  2021-6-15, lizt
                    // const int rSS=19;

                    const auto& slipSystem(*poly.grain(grainID).slipSystems()[rSS]);
                    const VectorDimD b(slipSystem.s.cartesian());

                    // Compute the ReciprocalLatticeDirection corresponding to s
                    ReciprocalLatticeDirection<3> sr(poly.grain(grainID).reciprocalLatticeDirection(b));

                    bool isEdge=true;



                    LatticeDirection<3> d1(LatticeVector<dim>(sr.cross(slipSystem.n)*randomSign()));
                    double d1cNorm(d1.cartesian().norm());
                    //                    const double size = distribution(generator)*inclusionsDiameterLognormalDistribution_A[f]/poly.b_SI;

                    int a1=sizeDistribution(generator)/d1cNorm;
                    // 2021-6-15, lizt
                    // int a1=50.0/d1cNorm;
                    if(a1>0)
                    {
                        LatticeVector<dim> L1=L0+d1*a1;



                        // Compute the LatticeDireciton corresponding to -n
                        LatticeDirection<3> d2(poly.grain(grainID).latticeDirection(-slipSystem.n.cartesian()*randomSign()));
                        double d2cNorm(d2.cartesian().norm());

                        const int a2=aspectRatioDistribution(generator)*a1; // aspect ratio of double FR source
                        if(a2>0)
                        {
                            LatticeVector<dim> L2=L1+d2*a2;
                            LatticeVector<dim> L3=L0+d2*a2;

                            const VectorDimD P0=L0.cartesian();
                            const VectorDimD P1=L1.cartesian();
                            const VectorDimD P2=L2.cartesian();
                            const VectorDimD P3=L3.cartesian();
                            const auto search1(mesh.search(P1));
                            const auto search2(mesh.search(P2));
                            const auto search3(mesh.search(P3));

                            double dh=0.0;
                            std::vector<VectorDimD> nodePos;

                            if(   search1.first && search1.second->region->regionID==grainID
                               && search2.first && search2.second->region->regionID==grainID
                               && search3.first && search3.second->region->regionID==grainID
                               )
                            {
                                 std::cout<<"        addding FR source~~~"<<std::endl;
                                //density += 2.0*(d1cNorm*a1 + d2cNorm*a2)/mesh.volume()/std::pow(poly.b_SI,2);
                                if(isEdge)
                                {
                                    edgeDensity+=2.0*(d1cNorm*a1 + d2cNorm*a2)/mesh.volume()/std::pow(poly.b_SI,2);
                                }
                                //std::cout<<"density="<<density<<" (edge density="<<edgeDensity<<")"<<std::endl;

                                const VectorDimD P4=0.5*(P0+P1);
                                const VectorDimD P5=0.5*(P2+P3);

                                const VectorDimD n1=slipSystem.unitNormal;
                                const VectorDimD n2=d2.cross(d1).cartesian().normalized();

                                configIO.nodes().emplace_back(nodeID+0,P0,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                                configIO.nodes().emplace_back(nodeID+1,P1,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                                configIO.nodes().emplace_back(nodeID+2,P2,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                                configIO.nodes().emplace_back(nodeID+3,P3,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                                configIO.nodes().emplace_back(nodeID+4,P4,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                                configIO.nodes().emplace_back(nodeID+5,P5,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);

                                configIO.loops().emplace_back(loopID+0,b,n1,P0,grainID,DislocationLoopIO<dim>::GLISSILELOOP,-1,VectorDimD::Zero());
                                configIO.loops().emplace_back(loopID+1,b,n2,P0,grainID,DislocationLoopIO<dim>::SESSILELOOP,-1,VectorDimD::Zero());
                                configIO.loops().emplace_back(loopID+2,b,n1,P3,grainID,DislocationLoopIO<dim>::SESSILELOOP,-1,VectorDimD::Zero());

                                configIO.links().emplace_back(loopID+0,nodeID+0,nodeID+1,0);
                                configIO.links().emplace_back(loopID+0,nodeID+1,nodeID+4,0);
                                configIO.links().emplace_back(loopID+0,nodeID+4,nodeID+0,0);

                                configIO.links().emplace_back(loopID+1,nodeID+0,nodeID+3,0);
                                configIO.links().emplace_back(loopID+1,nodeID+3,nodeID+2,0);
                                configIO.links().emplace_back(loopID+1,nodeID+2,nodeID+1,0);
                                configIO.links().emplace_back(loopID+1,nodeID+1,nodeID+0,0);

                                configIO.links().emplace_back(loopID+2,nodeID+3,nodeID+5,0);
                                configIO.links().emplace_back(loopID+2,nodeID+5,nodeID+2,0);
                                configIO.links().emplace_back(loopID+2,nodeID+2,nodeID+3,0);

                                nodeID+=6;
                                loopID+=3;
                                snID++;
                            }
                        }
                    
                    }
        }
        /**********************************************************************/
        void createVertices2(   const DDconfigIO<dim>& evl,
                                std::map<size_t,std::shared_ptr<NodeType>>& tempNodes,
                                const size_t& initialNodesize)
        {/*!Creates DislocationNode(s) based on the data read by the DDconfigIO<dim>
          * object.
          */
            size_t kk(initialNodesize);
            for (const auto& node : evl.nodes())
            {
                const size_t nodeIDinFile(node.sID);
                NodeType::set_count(nodeIDinFile);
                if(node.sID==node.masterID)
                {// a regular node is created
                    model::cout<<"Creating DislocationNode "<<nodeIDinFile<<" ("<<kk<<" of "<<evl.nodes().size()+initialNodesize<<")"<<std::endl;
                    
                    const size_t nodeID(StaticID<NodeType>::nextID());
                    const auto inserted(tempNodes.emplace(std::piecewise_construct,
                                                          std::make_tuple(nodeID),
                                                          std::make_tuple(new NodeType(this,node.P,node.V,node.velocityReduction))));
                    assert(inserted.second && "COULD NOT INSERT NETWORK VERTEX IN VERTEX CONTAINER.");
                    assert(inserted.first->first == nodeID && "KEY != nodeID");
                    assert(inserted.first->second->sID == nodeID && "sID != nodeID");
                    assert(nodeID==nodeIDinFile);
                }
                else
                {
                    model::cout<<"Creating DislocationNode "<<nodeIDinFile<<" ("<<kk<<" of "<<evl.nodes().size()+initialNodesize<<"), virtual of "<<node.masterID<<std::endl;
                    const auto isNode(this->node(node.masterID));
                    assert(isNode.first);
                    isNode.second->resetVirtualBoundaryNode();
                }
                kk++;
            }
        }


        /**********************************************************************/
        void createEdges2(  const DDconfigIO<dim>& evl,
                            const std::map<size_t,
                            std::shared_ptr<NodeType>>& tempNodes,
                            const size_t& initialNodeSize,
                            const size_t& initialLoopSize)
        {
            std::map<size_t,std::map<size_t,size_t>> loopMap;
            for(const auto& looplink : evl.links())
            {// Collect LoopLinks by loop IDs
                loopMap[looplink.loopID].emplace(looplink.sourceID,looplink.sinkID);
            }
            
            assert(loopMap.size()==evl.loops().size());
            
            size_t loopLumber=initialLoopSize;
            for(const auto& loop : evl.loops())
            {// for each loop in the DDconfigIO<dim> object
                
                const auto loopFound=loopMap.find(loop.sID); // there must be an entry with key loopID in loopMap
                assert(loopFound!=loopMap.end());
                std::vector<std::shared_ptr<NodeType>> loopNodes;
                loopNodes.push_back(getSharedNode(loopFound->second.begin()->first,tempNodes));
                for(size_t k=0;k<loopFound->second.size();++k)
                {
                    const auto nodeFound=loopFound->second.find(loopNodes.back()->sID);
                    if(k<loopFound->second.size()-1)
                    {
                        loopNodes.push_back(getSharedNode(nodeFound->second,tempNodes));
                    }
                    else
                    {
                        assert(nodeFound->second==loopNodes[0]->sID);
                    }
                }
                
                LoopType::set_count(loop.sID);
                
                const bool faulted= poly.grain(loop.grainID).rationalLatticeDirection(loop.B).rat.asDouble()!=1.0? true : false;
                
                model::cout<<"Creating Dislocation Loop "<<loop.sID<<" ("<<loopLumber<<" of "<<evl.loops().size()<<"), type="<<loop.loopType<<", faulted="<<faulted<<", |b|="<<loop.B.norm()<<std::endl;
                switch (loop.loopType)
                {
                    case DislocationLoopIO<dim>::GLISSILELOOP:
                    {
                        LatticePlane loopPlane(loop.P,poly.grain(loop.grainID).reciprocalLatticeDirection(loop.N));
                        GlidePlaneKey<dim> loopPlaneKey(loop.grainID,loopPlane);
                        if(simulationParameters.isPeriodicSimulation())
                        {
                            const size_t newLoopID=this->insertLoop(loopNodes,loop.B,glidePlaneFactory.get(loopPlaneKey),loop.periodicShift)->sID;
                            assert(loop.sID==newLoopID);
                        }
                        else
                        {
                            const size_t newLoopID=this->insertLoop(loopNodes,loop.B,glidePlaneFactory.get(loopPlaneKey))->sID;
                            assert(loop.sID==newLoopID);
                        }

                        break;
                    }
                        
                    case DislocationLoopIO<dim>::SESSILELOOP:
                    {
                        const size_t newLoopID=this->insertLoop(loopNodes,loop.B,loop.grainID,loop.loopType)->sID;
                        assert(loop.sID==newLoopID);
                        break;
                    }
                        
                    case DislocationLoopIO<dim>::VIRTUALLOOP:
                    {
                        const size_t newLoopID=this->insertLoop(loopNodes,loop.B,loop.grainID,loop.loopType)->sID;
                        assert(loop.sID==newLoopID);
                        break;
                    }
                        
                    default:
                        assert(false && "Unknown DislocationLoop type");
                        break;
                }
                
                loopLumber++;
            }
            
            model::cout<<std::endl;
        }
        
        /**********************************************************************/
        void updatePlasticDistortionFromAreas(const double& dt)
        {

            const MatrixDimD old(_plasticDistortionFromAreas);
            _plasticDistortionFromAreas.setZero();
            for(const auto& loop : this->loops())
            {
                _plasticDistortionFromAreas+= loop.second->plasticDistortion();
            }
            if(dt>0.0)
            {
                _plasticDistortionRateFromAreas=(_plasticDistortionFromAreas-old)/dt;
            }

        }    
        
        /**********************************************************************/
        const MatrixDimD& plasticDistortionRate() const
        {
            return  _plasticDistortionRateFromAreas;
        }
        
        /**********************************************************************/
        const MatrixDimD& plasticDistortion() const
        {
            return _plasticDistortionFromAreas;
        }
        
        /**********************************************************************/
        MatrixDimD plasticStrainRate() const
        {/*!\returns the plastic strain rate tensor generated during the last time step.
          */
            return (plasticDistortionRate()+plasticDistortionRate().transpose())*0.5;
        }
        
        /**********************************************************************/
        std::tuple<double,double,double,double> networkLength() const
        {/*!\returns the total line length of the DislocationNetwork. The return
          * value is a tuple, where the first value is the length of bulk glissile
          * dislocations, the second value is the length of bulk sessile
          * dislocations, and the third value is the length accumulated on
          * the mesh boundary.
          */
            double bulkGlissileLength(0.0);
            double bulkSessileLength(0.0);
            double boundaryLength(0.0);
            double grainBoundaryLength(0.0);
            
            for(auto& loop : this->loops())
            {
                for(const auto& loopLink : loop.second->links())
                {
                    if(!loopLink.second->pLink->hasZeroBurgers())
                    {
                        if(loopLink.second->pLink->isBoundarySegment())
                        {
                            boundaryLength+=loopLink.second->pLink->chord().norm();
                        }
                        else if(loopLink.second->pLink->isGrainBoundarySegment())
                        {
                            grainBoundaryLength+=loopLink.second->pLink->chord().norm();
                        }
                        else
                        {
                            if(loopLink.second->pLink->isSessile())
                            {
                                bulkSessileLength+=loopLink.second->pLink->chord().norm()/loopLink.second->pLink->loopLinks().size();
                            }
                            else
                            {
                                bulkGlissileLength+=loopLink.second->pLink->chord().norm()/loopLink.second->pLink->loopLinks().size();
                            }
                        }
                    }
                }
            }
            return std::make_tuple(bulkGlissileLength,bulkSessileLength,boundaryLength,grainBoundaryLength);
        }
        
        /**********************************************************************/
        MatrixDimD stress(const VectorDim& x) const
        {/*!\param[in] P position vector
          * \returns The stress field generated by the DislocationNetwork at P
          *
          * Note:
          */
            MatrixDimD temp(MatrixDimD::Zero());
            for(const auto& link : this->links())
            {// sum stress field per segment
                if(   !link.second->hasZeroBurgers()
                   && !(link.second->isBoundarySegment() && simulationParameters.simulationType==DefectiveCrystalParameters::FINITE_NO_FEM) // exclude boundary segments even if they are non-zero Burgers
                   //                   && !(link.second->isVirtualBoundarySegment() && simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_IMAGES)
                   )
                {
                    for(const auto& shift : periodicShifts)
                    {
                        temp+=link.second->straight.stress(x+shift);
                    }
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        void stress(std::deque<FEMfaceEvaluation<ElementType,dim,dim>>& fieldPoints) const
        {
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for(size_t k=0;k<fieldPoints.size();++k)
            {
                fieldPoints[k]=stress(fieldPoints[k].P);
            }
        }
        
        /**********************************************************************/
        VectorDim displacement(const VectorDim& x) const
        {/*!\param[in] P position vector
          * \returns The stress field generated by the DislocationNetwork at P
          *
          * Note:
          */
            
            VectorDim temp(VectorDim::Zero());
            
            for(const auto& loop : this->loops())
            {// sum solid angle of each loop
                for(const auto& shift : periodicShifts)
                {
                    temp-=loop.second->solidAngle(x+shift)/4.0/M_PI*loop.second->burgers();
                }

            }
            
            for(const auto& link : this->links())
            {// sum line-integral part of displacement field per segment
                if(   !link.second->hasZeroBurgers() )
                {
                    for(const auto& shift : periodicShifts)
                    {
                        temp+=link.second->straight.displacement(x+shift);
                    }
                }
            }
            
            return temp;
        }

        /**********************************************************************/
        int randomSign()
        {
            std::uniform_int_distribution<> dis(0,1);
            std::mt19937 generator(std::chrono::system_clock::now().time_since_epoch().count());
            return  dis(generator)*2-1;
        }
        /**********************************************************************/
        void displacement(std::vector<FEMnodeEvaluation<ElementType,dim,1>>& fieldPoints) const
        {
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for(size_t k=0;k<fieldPoints.size();++k)
            {
                fieldPoints[k]=displacement(fieldPoints[k].P);
            }
        }
        
        
        /**********************************************************************/
        DislocationNetworkIOType io()
        {
            return DislocationNetworkIOType(*this,folderSuffix);
        }
        
        /**********************************************************************/
        DislocationNetworkIOType io() const
        {
            return DislocationNetworkIOType(*this,folderSuffix);
        }
        
    };
    
}
#endif