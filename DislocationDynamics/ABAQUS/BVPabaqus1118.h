/* 2021-4-26, lizt
 *
 * store Abaqus solution
 * 
 * 2022-11-14, update
 * Recombine triangular mesh to hexahedron
 * 
 */

#ifndef model_BVPabaqus_H_
#define model_BVPabaqus_H_

#include <array>
#include <deque>
#include <chrono>
#include <memory> // unique_ptr
#ifdef _MODEL_PARDISO_SOLVER_
#include <Eigen/PardisoSupport>
#endif
#include<Eigen/SparseCholesky>
#include <FiniteElement.h>

#include <SimplicialMesh.h>
#include <JGNselector.h>
#include <SingleFieldPoint.h>
#include <FEMnodeEvaluation.h>
#include <FEMfaceEvaluation.h>

#include <LinearWeakList.h>

#include <NullSpaceSolver.h>

#include <TextFileParser.h>

#include <IDmatch.h>
#include <CoupleParameters.h>

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <unistd.h>
#include <string>



namespace model
{
    
    template <int dim, int sfOrder>
    class BVPabaqus
    {
        
        
        
    public:
        
        
        typedef LagrangeElement<dim,sfOrder> ElementType;
        typedef FiniteElement<ElementType> FiniteElementType;

        typedef IDmatch<dim,FiniteElementType> IDmatchType;

        typedef TrialFunction<'u',dim,FiniteElementType> TrialFunctionType;
        typedef TrialFunction<'s',6,FiniteElementType> StressTrialType;
        typedef TrialFunction<'Tem',6,FiniteElementType> TemTrialType;

        constexpr static int dofPerNode=TrialFunctionType::dofPerNode;
        typedef TrialGrad<TrialFunctionType> TrialGradType;
        typedef TrialDef<TrialFunctionType> TrialDefType;
        typedef Eigen::Matrix<double,6,6> CmatrixType;
        typedef Constant<CmatrixType,6,6> CconstantType;
        typedef TrialProd<CconstantType,TrialDefType> TrialStressType;

        typedef BilinearForm<TrialDefType,TrialStressType> BilinearFormType;
        typedef IntegrationDomain<FiniteElementType,0,4,GaussLegendre> IntegrationDomainType;
        typedef BilinearWeakForm<BilinearFormType,IntegrationDomainType> BilinearWeakFormType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
        typedef Eigen::SparseMatrix<double> SparseMatrixType;       
        

        const SimplicialMesh<dim>& mesh;

        CoupleParameters<dim> para;
        size_t nodeSize;
        size_t elementSize;
        
        //        static bool apply_DD_displacement;
        
    private:
        Eigen::Matrix<double,6,6> C; // matrix of elastic moduli
        FiniteElementType fe;
        std::pair<double,std::set< const ElementType*>> _elementSet; // key=characteristic volume, value=elementSet

        TrialFunctionType  u;  // displacement field u=[u1 u2 u3]'
        // TrialGradType  b;      // displacement gradient b=[u11 u12 u13 u21 u22 u23 u31 u32 u33]'
        // TrialDefType  e;       // strain e=[e11 e22 e33 e12 e23 e13]'
        StressTrialType s;     // stress s=[s11 s22 s33 s12 s23 s13]'
        TemTrialType Tem;
        
        // SparseMatrixType A;
        // SparseMatrixType T;
        // SparseMatrixType A1;
        
        /**********************************************************************/
        Eigen::Matrix<double,dim+1,1> face2domainBary(const Eigen::Matrix<double,dim,1>& b1,
        /*                                         */ const int& boundaryFace) const
        {
            // Transform to barycentric coordinate on the volume, adding a zero on the boundaryFace-face
            Eigen::Matrix<double,dim+1,1> bary;
            for (int k=0;k<dim;++k)
            {
                bary((k<boundaryFace)? k : k+1)=b1(k);
            }
            bary(boundaryFace)=0.0;
            return bary;
        }

        /**********************************************************************/
        Eigen::Matrix<double,6,6> get_C(const double& mu, const double& nu) const
        {
//            const double  mu=1.0;  // dimensionless
//            const double  nu=Material<dim,Isotropic>::nu;
            const double lam=2.0*mu*nu/(1.0-2.0*nu);
            const double C11(lam+2.0*mu);
            const double C12(lam);
            const double C44(mu); // C multiplies engineering strain
            
            Eigen::Matrix<double,6,6> temp;
            temp<<C11, C12, C12, 0.0, 0.0, 0.0,
            /***/ C12, C11, C12, 0.0, 0.0, 0.0,
            /***/ C12, C12, C11, 0.0, 0.0, 0.0,
            /***/ 0.0, 0.0, 0.0, C44, 0.0, 0.0,
            /***/ 0.0, 0.0, 0.0, 0.0, C44, 0.0,
            /***/ 0.0, 0.0, 0.0, 0.0, 0.0, C44;
            return temp;
        }

    public:
        Eigen::Matrix<double,Eigen::Dynamic,6> stressAbaqus; 
        Eigen::Matrix<double,Eigen::Dynamic,6> TemAbaqus; 
        Eigen::Matrix<double,Eigen::Dynamic,6> stressIntAba;    // intrgral point stress in ABAQUS
        Eigen::Matrix<double,Eigen::Dynamic,6> TemIntAba;   
        Eigen::VectorXd stressV;                                // voigt stress in all node
        Eigen::VectorXd TemV;
        Eigen::Matrix<double,Eigen::Dynamic,dim> integralPoints;
        Eigen::Matrix<double,Eigen::Dynamic,6> plasticStrain;

        IDmatchType IDmap;


        //luosc
        const double epsilon; //used to regluarized the plastic strain accprding to YP JMPS paper. 
        const double expandheight;
        // const double ele_h;
        // const double nc;
        double volumeC;


        /**********************************************************************/
        template <typename DislocationNetworkType>
        BVPabaqus(const SimplicialMesh<dim>& mesh_in,const DislocationNetworkType& DN) :
        /* init  */ mesh(mesh_in)
        /* init  */,para(mesh)
        /* init  */,nodeSize(0)
        /* init  */,elementSize(0)
        /* init  */,C(get_C(DN.poly.mu,DN.poly.nu))
        /* init  */,fe(mesh)        
        /* init  */,u(fe.template trial<'u',dim>())
        /* init  */,s(fe.template trial<'s', 6 >())
        /* init  */,IDmap(fe)
                //luosc//
        /* init  */,epsilon(TextFileParser("./inputFiles/DD.txt").readScalar<double>("epsilon",true))  
        /* init  */,expandheight(TextFileParser("./inputFiles/DD.txt").readScalar<double>("expandheight",true)) 
        // /* init  */,ele_h(TextFileParser("./inputFiles/DD.txt").readScalar<double>("ele_h",true)) 
    //    /* init  */,nc(TextFileParser("./inputFiles/DD.txt").readScalar<double>("g_nc",true))


        {
            
            
            std::cout<<"Initializing BVPabaqus"<<std::endl;
            nodeSize = TrialBase<TrialFunctionType>::nodeSize();
            elementSize = TrialBase<TrialFunctionType>::elementSize()/6;    // element size for C3D8
            
            stressV = Eigen::VectorXd::Zero(nodeSize*6);
            TemV = Eigen::VectorXd::Zero(nodeSize*6);
            stressAbaqus = Eigen::Matrix<double,Eigen::Dynamic,6>::Zero(nodeSize,6);
            TemAbaqus = Eigen::Matrix<double,Eigen::Dynamic,6>::Zero(nodeSize,6);        
            stressIntAba = Eigen::Matrix<double,Eigen::Dynamic,6>::Zero(elementSize,6);    
            TemIntAba = Eigen::Matrix<double,Eigen::Dynamic,6>::Zero(elementSize,6);    
            integralPoints = Eigen::Matrix<double,Eigen::Dynamic,dim>::Zero(elementSize,dim);    
            

            //model::cout<<"          reading integral points from ABAQUS ..."<<std::endl;
            
            clearPlasticStrain();


        }


        /**********************************************************************/
        void calNodalStress(const long int& runID)
        {
            // calculate Nodal Stress from Integral point stress
            //std::string fileNameFlag("flag/flagS_0.txt");
             std::string fileNameFlag("flag/flagS_"+std::to_string(runID)+".txt");
           // std::string fileNameS("stress/stress_0.txt");
            std::string fileNameS("stress/stress_"+std::to_string(runID)+".txt");
            std::string fileDir("./");
            int num=0;
            std::ifstream openFile(fileDir+fileNameFlag);
            if(!openFile.is_open())
            {
                double waitTime(clock());
                model::cout<<fileNameS<<" is not ready yet, sleeping... "<<std::flush;
                while (!openFile.is_open())
                {
                    openFile.open(fileDir+fileNameFlag,std::ios::in);
                    sleep(2);
                    num = num+1;
                    assert(num<1000 && "sleep, good night~ ");
                }
                model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-waitTime)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
            }
            openFile.close();
            // const auto t0= std::chrono::system_clock::now(); 
            double t0(clock());
            std::ifstream openFileS(fileDir+fileNameS);
            model::cout<<"				Start reading stress from "<<fileNameS<<std::flush;
            // stressIntAba = TextFileParser(fileDir).readMatrix<double>("integral stress",elementSize,6,false); 
            std::string line;
            for (size_t id=0;id<elementSize;id++)
            {
                std::getline(openFileS,line);
                std::stringstream ss(line);
                for (size_t i = 0; i < 6; i++)
                {
                    ss>>stressIntAba(id,i);
                }
            }
            // model::cout<<stressIntAba<<std::endl;
            // model::cout<<"				read integral point stress from "<<fileNameS<<" is done... "<<std::flush;
            model::cout<<" , done ! "<<std::flush;
            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
            // const std::map<int,std::vector<int>> list(IDmap.NodeEleList);
            stressAbaqus = Eigen::Matrix<double,Eigen::Dynamic,6>::Zero(nodeSize,6);    
            model::cout<<"				clear nodal stress ABAQUS ... "<<std::flush;
            const auto t1= std::chrono::system_clock::now(); 
            auto list(IDmap.NodeEleList());
            for(const auto& iter : list)    // ABAQUS ID start from 1, <nodeID, <eleID1, eleID2, ....>>
            {
                // stressAbaqus.row(iter.first-1) = Eigen::Matrix<double,1,6>::Zero(1,6);
                // model::cout<<"1"<<std::endl;
                for(const auto& eleID : iter.second)
                {
                    stressAbaqus.row(iter.first-1) =  stressAbaqus.row(iter.first-1) + stressIntAba.row(eleID-1);
                }
                stressAbaqus.row(iter.first-1) = stressAbaqus.row(iter.first-1)/iter.second.size();
            }
            // for(int i=0;i<nodeSize;i++)   // nodeID
            // {
            //     for(int j=0;j<6;j++)    // stress tensor 
            //     {
            //         stressAbaqus(i,j) = 0.0;
            //         for (const auto& it : list[i])  // elementID
            //        {
            //            stressAbaqus(i,j) =+ stressIntAba(it,j);
            //         }
            //     }
            // }
            // Test: output nodal stress
            // std::ofstream outfile("./AbaqusInput/nodalS.txt");
            // outfile<<stressAbaqus;
            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]."<<defaultColor<<std::endl;

        }

        /**********************************************************************/
        void calNodalTem(const long int& runID)
        {
            // calculate Nodal Stress from Integral point stress
            //std::string fileNameFlag("flag/flagS_0.txt");
             std::string fileNameFlag("flag/flagT_"+std::to_string(runID)+".txt");
           // std::string fileNameS("stress/stress_0.txt");
            std::string fileNameT("Tem/Tem_"+std::to_string(runID)+".txt");
            std::string fileDir("./");
            int num=0;
            std::ifstream openFile(fileDir+fileNameFlag);
            if(!openFile.is_open())
            {
                double waitTime(clock());
                model::cout<<fileNameT<<" is not ready yet, sleeping... "<<std::flush;
                while (!openFile.is_open())
                {
                    openFile.open(fileDir+fileNameFlag,std::ios::in);
                    sleep(2);
                    num = num+1;
                    assert(num<1000 && "sleep, good night~ ");
                }
                model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-waitTime)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
            }
            openFile.close();
            // const auto t0= std::chrono::system_clock::now(); 
            double t0(clock());
            std::ifstream openFileT(fileDir+fileNameT);
            model::cout<<"				Start reading stress from "<<fileNameT<<std::flush;
            // stressIntAba = TextFileParser(fileDir).readMatrix<double>("integral stress",elementSize,6,false); 
            std::string line;
            for (size_t id=0;id<elementSize;id++)
            {
                std::getline(openFileT,line);
                std::stringstream tt(line);
                for (size_t i = 0; i < 6; i++)
                {
                    tt>>TemIntAba(id,i);
                }
            }
            // model::cout<<stressIntAba<<std::endl;
            // model::cout<<"				read integral point stress from "<<fileNameS<<" is done... "<<std::flush;
            model::cout<<" , done ! "<<std::flush;
            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
            // const std::map<int,std::vector<int>> list(IDmap.NodeEleList);
            TemAbaqus = Eigen::Matrix<double,Eigen::Dynamic,6>::Zero(nodeSize,6);    
            model::cout<<"				clear nodal Tem ABAQUS ... "<<std::flush;
            const auto t1= std::chrono::system_clock::now(); 
            auto list(IDmap.NodeEleList());
            for(const auto& iter : list)    // ABAQUS ID start from 1, <nodeID, <eleID1, eleID2, ....>>
            {
                // stressAbaqus.row(iter.first-1) = Eigen::Matrix<double,1,6>::Zero(1,6);
                // model::cout<<"1"<<std::endl;
                for(const auto& eleID : iter.second)
                {
                    TemAbaqus.row(iter.first-1) =  TemAbaqus.row(iter.first-1) + TemIntAba.row(eleID-1);
                }
                TemAbaqus.row(iter.first-1) = TemAbaqus.row(iter.first-1)/iter.second.size();
            }
            // for(int i=0;i<nodeSize;i++)   // nodeID
            // {
            //     for(int j=0;j<6;j++)    // stress tensor 
            //     {
            //         stressAbaqus(i,j) = 0.0;
            //         for (const auto& it : list[i])  // elementID
            //        {
            //            stressAbaqus(i,j) =+ stressIntAba(it,j);
            //         }
            //     }
            // }
            // Test: output nodal stress
            // std::ofstream outfile("./AbaqusInput/nodalS.txt");
            // outfile<<stressAbaqus;
            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]."<<defaultColor<<std::endl;

        }
        /**********************************************************************/
        void updateStress(const long int& runID)
        {
            // stressAbaqus = TextFileParser("./AbaqusInput/stress.txt").readMatrix<double>("nodal stress",nodeSize,6,false);
            calNodalStress(runID);
            for (size_t i=0;i<nodeSize;i++)
            {
                for (int j=0;j<6;j++)
                {
                    stressV(IDmap.NodeMap()[i+1]*6+j) = stressAbaqus(i,j); 
                }
            }
            stress()=stressV;
        }
        /**********************************************************************/
       //void updateTem(const long int& runID)
       // {
            // stressAbaqus = TextFileParser("./AbaqusInput/stress.txt").readMatrix<double>("nodal stress",nodeSize,6,false);
         //   calNodalTem(runID);
         //   for (size_t i=0;i<nodeSize;i++)
         //   {
         //       for (int j=0;j<6;j++)
          //      {
           //         TemV(IDmap.NodeMap()[i+1]*6+j) = TemAbaqus(i,j); 
           //     }
           // }
           // Tem()=TemV; 
    //    }
        /**********************************************************************/
        void updateIntPoint()
        {
            // Eigen::Matrix<double,Eigen::Dynamic,dim> intPAbaqus;
            // model::cout<<"                reading integral points position from file ..."<<std::endl;
            // intPAbaqus = TextFileParser("./AbaqusInput/integralPoints.txt").readMatrix<double>("integral points",elementSize,dim,false);
            //model::cout<<"get!"<<std::endl;
            model::cout<<"                get integral points position by calculating ..."<<std::endl;
            // intPAbaqus = IDmap.intPointsAba;
            for (size_t i=0;i<elementSize;i++)
            {
                // model::cout<<"i = "<<std::endl;
                // integralPoints.row(IDmap.ElementMap()[i+1]) = IDmap.intPointsAba.row(i);
                integralPoints = IDmap.intPointsAba;
            }
        }

        
        /**********************************************************************/
        const FiniteElementType& finiteElement() const
        {
            return fe;
        }
        
        /**********************************************************************/
        FiniteElementType& finiteElement()
        {
            return fe;
        }
        
        /**********************************************************************/
        const TrialFunctionType& displacement() const
        {
            return u;
        }
        
        /**********************************************************************/
        TrialFunctionType& displacement()
        {
            return u;
        }
        
        /**********************************************************************/
        const StressTrialType& stress() const
        {
            return s;
        }

        /**********************************************************************************/
        StressTrialType& stress()
        {
            return s;
        }

        /**********************************************************************************/
        /**********************************************************************/
 //       const TemTrialType& Tem() const
 //       {
  //          return Tem;
 //       }

        /**********************************************************************************/
//        TemTrialType& Tem()
//        {
 //           return Tem;
//        }

        /**********************************************************************************/
        Eigen::Matrix<double,dim,1> displacement(const Eigen::Matrix<double,dim,1> P,
                                                const Simplex<dim,dim>* guess) const
        {
            return eval(u)(P,guess);
        }

        /**********************************************************************/
        Eigen::Matrix<double,dim,dim> stress(const Eigen::Matrix<double,dim,1> P,
                                             const Simplex<dim,dim>* guess) const
        {
            Eigen::Matrix<double,6,1> tempV(eval(s)(P,guess));
            Eigen::Matrix<double,dim,dim> tempM;
            tempM(0,0)=tempV(0); // s11
            tempM(1,1)=tempV(1); // s22
            tempM(2,2)=tempV(2); // s33
            tempM(1,0)=tempV(3); // s21
            tempM(2,1)=tempV(4); // s32
            tempM(2,0)=tempV(5); // s31
            tempM(0,1)=tempM(1,0); //symm
            tempM(1,2)=tempM(2,1); //symm
            tempM(0,2)=tempM(2,0); //symm
            
            return tempM;
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,dim,dim> stress(const ElementType& ele,
                                             const Eigen::Matrix<double,dim+1,1>& bary) const
        {
            Eigen::Matrix<double,6,1> tempV(eval(s)(ele,bary));
            Eigen::Matrix<double,dim,dim> tempM;
            tempM(0,0)=tempV(0); // s11
            tempM(1,1)=tempV(1); // s22
            tempM(2,2)=tempV(2); // s33
            tempM(1,0)=tempV(3); // s21
            tempM(2,1)=tempV(4); // s32
            tempM(2,0)=tempV(5); // s31
            tempM(0,1)=tempM(1,0); //symm
            tempM(1,2)=tempM(2,1); //symm
            tempM(0,2)=tempM(2,0); //symm
            
            return tempM;
        }

        // 2021-4-2, lizt
        void  findElementSet(const std::vector<VectorDim>& Point, const VectorDim& faceNormal)
        {
            std::set<const ElementType*> elementSet;
            const ElementType* firstElement;
            double findShift(para.findShift);  // distance in a step for finding simplex, relative with element size
            double volume(0);
            int layers(para.layers);

            firstElement=finiteElement().search(Point[0]).second;
            elementSet.insert(firstElement);
            for (size_t si=1;si<Point.size();si++)
            {
                elementSet.insert(finiteElement().searchWithGuess(Point[si],&(firstElement->simplex)).second);
            }
            for (const auto& iter : elementSet )
            {
                volume+=iter->simplex.vol0;
            }
            volume=volume/elementSet.size(); // there is not just 1 element
            // model::cout<<"element set size = "<<elementSet.size()<<std::endl;

            for(int layer=-layers;layer<0;layer++)
            {
                for(size_t si=0;si<Point.size();si++)
                {
                    elementSet.insert(finiteElement().searchWithGuess(Point[si]+layer*findShift*faceNormal,&(firstElement->simplex)).second);
                }
            }
            // model::cout<<"element set size = "<<elementSet.size()<<std::endl;

            for(int layer=1;layer<=layers;layer++)
            {
                for(size_t si=0;si<Point.size();si++)
                {
                    elementSet.insert(finiteElement().searchWithGuess(Point[si]+layer*findShift*faceNormal,&(firstElement->simplex)).second);
                }
            }
            assert(elementSet.size() && "Element Set is Empty ...");
            // model::cout<<"element set size = "<<elementSet.size()<<std::endl;
            getElementSet()=std::make_pair(volume,elementSet);

        }


        // 2021-5-14, lizt
/*         std::map<int,double>  findElementSetAndWeight(const std::vector<VectorDim>& Point, const VectorDim& faceNormal )
        {
            std::map<int,double> weightMap; // key is id of dd element, value is element plastic weight
            const Simplex<dim,dim>* firstElement;
            double findShift(para.findShift);  // distance in a step for finding simplex, relative with element size
            double corFactor; // truncation error correction factor
            double g_ep(0.0); 
            int layers(para.layers);

            for (int n=-layers; n<=layers; n++)
            {
                g_ep+=d_Gamma(fabs(n*findShift));
            }
            corFactor = 1.0/Point.size()/g_ep;

            firstElement=mesh.search(Point[0]).second;
            volumeC=firstElement->vol0;

            for(int layer=-layers;layer<=layers;layer++)
            {
                double weight(d_Gamma(fabs(layer*findShift)));
                           
                for(size_t si=0;si<Point.size();si++)
                {
                    int id(-1);
                    id = finiteElement().s2eID[mesh.searchWithGuess(Point[si]+layer*findShift*faceNormal,firstElement).second->xID];
                    assert(id!=-1 && "find simplex failed .");
                    weightMap[id]+=weight*corFactor;                   
                }
            }
            return weightMap;

        } */

        // 2022-11-14, lizt
        // std::map<int,double>  findElementSetAndWeight(const std::vector<VectorDim>& Point, const VectorDim& faceNormal )
        // {
        //     std::map<int,double> weightMap; // key is id of C3D8 element, value is element plastic weight
        //     const Simplex<dim,dim>* firstElement;
        //     double findShift(para.findShift);  // distance in a step for finding simplex, relative with element size
        //     double corFactor; // truncation error correction factor
        //     double g_ep(0.0); 
        //     int layers(para.layers);
        //     double deviation=0.01;

        //     for (int n=-layers; n<=layers; n++)
        //     {
        //         if(n==0)
        //         {
        //             g_ep+=d_Gamma(fabs(deviation))*2;
        //         }
        //         else
        //         {
        //             g_ep+=d_Gamma(fabs(n*findShift));
        //         }
        //         // g_ep+=d_Gamma(fabs(n*findShift));
        //     }
        //     corFactor = 1.0/Point.size()/g_ep;

        //     firstElement=mesh.search(Point[0]).second;
        //     volumeC=firstElement->vol0;
        //     // test
        //     // model::cout<<"volume = "<<volumeC<<std::endl;

        //     for(int layer=-layers;layer<=layers;layer++)
        //     {
        //         if (layer==0)
        //         {
        //             for(size_t si=0;si<Point.size();si++)
        //             {
        //                 int id1(-1);     // C3D4 element ID in DD
        //                 int id2(-1);
        //                 double weight(d_Gamma(fabs(deviation)));
        //                 id1 = finiteElement().s2eID[mesh.searchWithGuess(Point[si]+deviation*faceNormal,firstElement).second->xID];
        //                 id2 = finiteElement().s2eID[mesh.searchWithGuess(Point[si]-deviation*faceNormal,firstElement).second->xID];
        //                 assert(id1!=-1 && id2!=-1 && "find simplex failed .");
        //                 weightMap[IDmap.ElementMap()[id1]]+=weight*corFactor;   
        //                 weightMap[IDmap.ElementMap()[id2]]+=weight*corFactor;   
        //             }
        //         }
        //         else
        //         {
        //             for(size_t si=0;si<Point.size();si++)
        //             {
        //                 int id(-1);
        //                 double weight(d_Gamma(fabs(layer*findShift)));
        //                 id = finiteElement().s2eID[mesh.searchWithGuess(Point[si]+layer*findShift*faceNormal,firstElement).second->xID];
        //                 assert(id!=-1 && "find simplex failed .");
        //                 weightMap[IDmap.ElementMap()[id]]+=weight*corFactor;   
        //             }
        //         }
                
                           
        //         // for(size_t si=0;si<Point.size();si++)
        //         // {
        //         //     id = finiteElement().s2eID[mesh.searchWithGuess(Point[si]+layer*findShift*faceNormal,firstElement).second->xID];
        //         //     assert(id!=-1 && "find simplex failed .");
        //         //     weightMap[IDmap.ElementMap()[id]]+=weight*corFactor;   
        //         // }
                
        //         }
        //         return weightMap;
        //     }

        // 2022-11-22, lizt
        std::map<int,double>  findElementSetAndWeight(const std::vector<VectorDim>& Point, const VectorDim& faceNormal )
        {
            std::map<int,double> weightMap; // key is id of C3D8 element, value is element plastic weight
            const Simplex<dim,dim>* firstElement;
            double findShift(para.findShift);  // distance in a step for finding simplex, relative with element size
            double corFactor; // truncation error correction factor
            double g_ep(0.0); 
            int layers(para.layers);
            double deviation=0.01;

            for (int n=-layers; n<=layers; n++)
            {
                if(n==0)
                {
                    g_ep+=d_Gamma(fabs(deviation))*2;
                }
                else
                {
                    g_ep+=d_Gamma(fabs(n*findShift));
                }
                // g_ep+=d_Gamma(fabs(n*findShift));
            }
            corFactor = 1.0/Point.size()/g_ep;

            firstElement=mesh.search(Point[0]).second;
            volumeC=firstElement->vol0;
            // test
            // model::cout<<"volume = "<<volumeC<<std::endl;

            for(int layer=-layers;layer<=layers;layer++)
            {
                if (layer==0)
                {
                    for(size_t si=0;si<Point.size();si++)
                    {
                        int id1(-1);     // C3D4 element ID in DD
                        int id2(-1);
                        double weight(d_Gamma(fabs(deviation)));
                        id1 = finiteElement().s2eID[mesh.searchWithGuess(Point[si]+deviation*faceNormal,firstElement).second->xID];
                        id2 = finiteElement().s2eID[mesh.searchWithGuess(Point[si]-deviation*faceNormal,firstElement).second->xID];
                        // assert(id1!=-1 && id2!=-1 && "find simplex failed .");
                        weightMap[IDmap.ElementMap()[id1]]+=weight*corFactor;   
                        weightMap[IDmap.ElementMap()[id2]]+=weight*corFactor;   
                    }
                }
                else
                {
                    for(size_t si=0;si<Point.size();si++)
                    {
                        int id(-1);
                        double weight(d_Gamma(fabs(layer*findShift)));
                        // auto SearchElement = mesh.searchWithGuess(Point[si]+layer*findShift*faceNormal,firstElement);
                        // assert(id!=-1 && "find simplex failed .");
                        // if (SearchElement.first)
                        // {
                        id = finiteElement().s2eID[mesh.searchWithGuess(Point[si]+layer*findShift*faceNormal,firstElement).second->xID];
                        weightMap[IDmap.ElementMap()[id]]+=weight*corFactor;   
                            // model::cout<<"id = "<<id<<" ID map = "<<IDmap.ElementMap()[id]<<std::endl;
                        // }
                        
                        // model::cout<<"done "<<std::endl;
                    }
                }
                
                           
                // for(size_t si=0;si<Point.size();si++)
                // {
                //     id = finiteElement().s2eID[mesh.searchWithGuess(Point[si]+layer*findShift*faceNormal,firstElement).second->xID];
                //     assert(id!=-1 && "find simplex failed .");
                //     weightMap[IDmap.ElementMap()[id]]+=weight*corFactor;   
                // }
                
                }
                return weightMap;
            }


        /**********************************************************************/
        std::pair<double,std::set< const ElementType*>>& getElementSet()
        {
            return _elementSet;
        }

        /**********************************************************************/
        void clearPlasticStrain()
        {
            model::cout<<"                clear plastic strian ..."<<std::endl;
            plasticStrain = Eigen::Matrix<double,Eigen::Dynamic,6>::Zero(elementSize,6);
        }

        /**********************************************************************/
        Eigen::Matrix<double,Eigen::Dynamic,6>& PlasticStrain()
        {
            return plasticStrain;
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
           double ele_h = para.findShift;
           double  dG = (g_function(ele_R+ele_h/2)-g_function(ele_R-ele_h/2))/(2*g_function(para.findShift*para.layers));

           return dG;

       }


    };

}       
        
#endif
