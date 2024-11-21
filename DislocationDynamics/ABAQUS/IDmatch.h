/* 2021-4-26, lizt
 *
 * store Abaqus solution
 * 
 * 
 * 
 */

#ifndef model_IDmatch_H_
#define model_IDmatch_H_

#include <TextFileParser.h>
#include <Eigen/Dense>
#include <map>
#include <vector>
#include <array>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>

#include <FiniteElement.h>

#ifdef _OPENMP
#include <omp.h>
#define EIGEN_DONT_PARALLELIZE // disable Eigen Internal openmp Parallelization
#endif

namespace model
{
    template <int dim, typename FiniteElementType>
    class IDmatch
    {
    public:
        const size_t nodeSize;
        const size_t elementSize;
        // typedef std::array<double,dim> Position;
        typedef std::map<int,std::vector<int>> listContainer; // 
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef std::map<int,VectorDim> mapContainer;

    private:
        mapContainer abaqusNodes; // nodeID in ABAQUS
        // mapContainer abaqusElements; // Position is barycenter of element
        mapContainer Nodes; // nodeID in DD
        listContainer Elements; // elementID in DD
        std::map<int,int> nodeMap;
        std::map<int,int> elementMap;   // <DD elementID, ABAQUS elementID>
        listContainer eleNodeList; // map<elementID,elementNodeID>,ID in ABAQUS
        listContainer nodeEleList;  //map<nodeID,elementID>, ID in ABAQUS


    public:

        std::vector<VectorDim> nodePosition;
        std::unique_ptr<FiniteElementType> fe_ptr;
        Eigen::Matrix<double,Eigen::Dynamic,dim> intPointsAba;


        IDmatch( FiniteElementType& fe) : 
        /* init  */ nodeSize(fe.nodeSize())
        /* init  */,elementSize(fe.elementSize()/6) // element size for C3D8
        /* init  */,fe_ptr(&fe)

        {
            model::cout<<"Initializing IDmatch"<<std::endl;
            readInp();
            initNodes();
            if (readElementMap()==0 || readNodeMap()==0 )
            {
                // readInp();
                makeNodeMap();
                makeElementMap();
            }
            // intPointsAba = Eigen::Matrix<double,elementSize,dim>::Zero();
        }

        bool readNodeMap()
        {
            std::string nodeMapDir("./AbaqusInput/nodeMap.txt");
            std::ifstream nodeMapFile(nodeMapDir.c_str(),std::ifstream::in);
            if (nodeMapFile.is_open())
            {
                model::cout<<"reading node map from file ..."<<std::flush;
                double t0(clock());
                std::string line;
                while (std::getline(nodeMapFile, line))
                {
                    std::stringstream ss(line);
                    int abaID;
                    int ddID;
                    ss>>abaID;
                    ss>>ddID;
                    const bool success=NodeMap().emplace(abaID,ddID).second;
                    if(!success)
                    {
                        std::cout<<"Unable to insert node map "<<line<<std::endl;
                    }
                }
                                
                nodeMapFile.close();                
                std::cout<<" ("<<NodeMap().size()<<" node map) "<<std::flush;
                model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;
                return 1;
            }
            else
            {
                return 0;
            }
        }

        bool readElementMap()
        {
            std::string elementMapDir("./AbaqusInput/elementMap.txt");
            std::ifstream elementMapFile(elementMapDir.c_str(),std::ifstream::in);
            if (elementMapFile.is_open())
            {
                model::cout<<"reading element map from file ..."<<std::flush;
                double t0(clock());
                std::string line;
                while (std::getline(elementMapFile, line))
                {
                    std::stringstream ss(line);
                    int abaID;
                    int ddID;
                    ss>>ddID;
                    ss>>abaID;
                    const bool success=ElementMap().emplace(ddID,abaID).second;
                    if(!success)
                    {
                        std::cout<<"Unable to insert element map "<<line<<std::endl;
                    }
                }
                                
                elementMapFile.close();                
                std::cout<<" ("<<ElementMap().size()<<" element map) "<<std::flush;
                model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<defaultColor<<std::endl;

                return 1;
            }
            else
            {
                return 0;
            }
        }

// read C3D4 element
        // void readInp()
        // {
        //     int index=0; // indicate weather read success or not
        //     std::string line;
        //     std::ifstream fin;
        //     nodePosition.clear();
        //     fin.open("./AbaqusInput/mesh.inp");
        //     while (std::getline(fin,line))
        //     {
        //         if (line == "*NODE")
        //         {                    
        //             int nodeID; // node ID in ABAQUS
        //             VectorDim tempPos;
        //             model::cout<<"read ABAQUS nodes "<<std::endl;
        //             model::cout<<nodeSize<<" = nodeSize"<<std::endl;
        //             for (size_t k=0;k<nodeSize;k++)
        //             {
                        
        //                 VectorDim nodePos;
        //                 std::getline(fin,line,',');
        //                 std::stringstream ss(line);
        //                 ss >> nodeID;
        //                 //model::cout<<"nodeID = "<<nodeID<<std::endl;
        //                 for (int d=0;d<dim-1;d++)
        //                 {
        //                     std::getline(fin,line,',');
        //                    std::stringstream ss1(line);
        //                     ss1>> nodePos(d);
        //                     // tempPos(d)=nodePos(d);
        //                     //model::cout<<"nodePosition = "<<nodePos[d]<<std::endl;
        //                 }
        //                 std::getline(fin,line);
        //                 std::stringstream ss2(line);
        //                 ss2>>nodePos(dim-1);
        //                 tempPos=nodePos;
        //                 abaqusNodes.emplace(nodeID,nodePos);
        //                 nodePosition.push_back(tempPos);
        //                 //model::cout<<nodePosition[k]<<std::endl;
        //                 //model::cout<<abaqusNodes.size()<<std::endl;
        //                 //model::cout<<"ID = "<<nodeID<<"  position = "<<tempPos.transpose()<<std::endl;
                        
        //             }
        //             // for (const auto& iter : abaqusNodes )
        //             // {
        //             //     for (int i=0;i<dim;i++)
        //             //     {
        //             //         model::cout<<iter.second<<std::ends<<iter.first[i]<<std::ends;
        //             //     }
        //             //     model::cout<<std::endl;
        //             // }
        //             index+=1;
        //             // model::cout<<"index = "<<index<<std::endl;
        //         }

        //         if (line == "*ELEMENT, type=C3D4, ELSET=Volume1")
        //         {
        //              // barycenter position
        //             int elementID;
        //             model::cout<<"read ABAQUS elements, "<<std::flush;
        //             model::cout<<elementSize<<" = C3D4  element size"<<std::endl;
        //             // intPointsAba = Eigen::Matrix<double,Eigen::Dynamic,dim>::Zero(elementSize,dim);
        //             for (int k=0;k<elementSize;k++)
        //             {
        //                 int elementNodeID;
        //                 std::set<int> tempEleNodeList;
        //                 std::getline(fin,line,',');
        //                 VectorDim tempPos(VectorDim::Zero());
        //                 std::stringstream ss(line);
        //                 ss >> elementID; // element ID in ABAQUS, ID is not start from 1     
        //                 //model::cout<<"elementID"<<elementID<<std::endl;        
        //                 for (int d=0;d<dim;d++)
        //                 {
        //                     std::getline(fin,line,',');
        //                     std::stringstream ss1(line);
        //                     ss1 >> elementNodeID;
        //                     //model::cout<<"elementNodeID = "<<elementNodeID<<std::endl;
        //                     tempPos+=nodePosition[elementNodeID-1]*0.25;
        //                     tempEleNodeList.push_back(elementNodeID);
        //                 }
        //                 std::getline(fin,line);
        //                 std::stringstream ss2(line);
        //                 ss2 >> elementNodeID;
        //                 tempPos+=nodePosition[elementNodeID-1]*0.25;

        //                 tempEleNodeList.push_back(elementNodeID);

        //                 // intPointsAba.row(k)+=nodePosition[elementNodeID-1].transpose()*0.125;

        //                 VectorDim baryPos(tempPos);
        //                 // for (int d=0;d<dim;d++)
        //                 // {
        //                 //     baryPos(d)=tempPos(d);
        //                 // }

        //                 abaqusElements.emplace(k+1,baryPos);
        //                 // eleNodeList.emplace(elementID,tempEleNodeList);
        //                 eleNodeList.emplace(k+1,tempEleNodeList);
        //             }
        //             model::cout<<"read inp success "<<std::endl;
        //             index+=1;
        //             break;
                    
                    
        //         }
        //     }
        //     // model::cout<<"index = "<<index<<std::endl;
        //     assert(index==2 && "read inp file failed ");


        // }


        /**********************************************************************/
        // lizt, 2022-11-14
        // read C3D8 element from mesh.inp
        void readInp()
        {
            const auto t0= std::chrono::system_clock::now();
            int index=0; // indicate weather read success or not
            std::string line;
            std::ifstream fin;
            nodePosition.clear();
            fin.open("./AbaqusInput/mesh-hex.inp");
            while (std::getline(fin,line))
            {
                if (line == "*NODE")
                {                    
                    int nodeID; // node ID in ABAQUS
                    VectorDim tempPos;
                    model::cout<<"read ABAQUS nodes, "<<std::flush;
                    model::cout<<nodeSize<<" = nodeSize"<<std::endl;
                    for (size_t k=0;k<nodeSize;k++)
                    {
                        
                        VectorDim nodePos;
                        std::getline(fin,line,',');
                        std::stringstream ss(line);
                        ss >> nodeID;
                        //model::cout<<"nodeID = "<<nodeID<<std::endl;
                        for (int d=0;d<dim;d++)
                        {
                            if(d==dim-1)
                            {
                                std::getline(fin,line);
                            }
                            else
                            {
                                std::getline(fin,line,',');
                            }
                            std::stringstream ss1(line);
                            ss1>> nodePos(d);
                            // tempPos(d)=nodePos(d);
                            //model::cout<<"nodePosition = "<<nodePos[d]<<std::endl;
                        }
                        tempPos=nodePos;
                        abaqusNodes.emplace(nodeID,nodePos);
                        nodePosition.push_back(tempPos);
                        //model::cout<<nodePosition[k]<<std::endl;
                        //model::cout<<abaqusNodes.size()<<std::endl;
                        //model::cout<<"ID = "<<nodeID<<"  position = "<<tempPos.transpose()<<std::endl;
                        
                    }
                    // for (const auto& iter : abaqusNodes )
                    // {
                    //     for (int i=0;i<dim;i++)
                    //     {
                    //         model::cout<<iter.second<<std::ends<<iter.first[i]<<std::ends;
                    //     }
                    //     model::cout<<std::endl;
                    // }
                    index+=1;
                    // model::cout<<"index = "<<index<<std::endl;
                }

                if (line == "*ELEMENT, type=C3D8, ELSET=Volume1")
                {
                     // barycenter position
                    int elementID;
                    // int hexEleSize(elementSize);
                    model::cout<<"read ABAQUS elements, "<<std::flush;
                    model::cout<<elementSize<<" = hexahedron element size. "<<std::endl;
                    intPointsAba = Eigen::Matrix<double,Eigen::Dynamic,dim>::Zero(elementSize,dim);
                    for (int k=0;k<elementSize;k++)
                    {
                        int elementNodeID;
                        std::vector<int> tempEleNodeList;
                        std::getline(fin,line,',');
                        VectorDim tempPos(VectorDim::Zero());
                        std::stringstream ss(line);
                        ss >> elementID; // element ID in ABAQUS, ID is not start from 1     
                        //model::cout<<"elementID"<<elementID<<std::endl;        
                        for (int d=0;d<8;d++)
                        {
                            if(d==7)
                            {
                                std::getline(fin,line);
                            }
                            else
                            {
                                std::getline(fin,line,',');
                            }
                            std::stringstream ss1(line);
                            ss1 >> elementNodeID;
                            //model::cout<<"elementNodeID = "<<elementNodeID<<std::endl;
                            tempPos+=nodePosition[elementNodeID-1]*0.125;
                            tempEleNodeList.push_back(elementNodeID);
                        }


                        intPointsAba.row(k)=tempPos.transpose();

                        // eleNodeList.emplace(elementID,tempEleNodeList);
                        eleNodeList.emplace(k+1,tempEleNodeList);
                    }
                    model::cout<<"read inp success "<<std::endl;
                    index+=1;
                    break;
                    
                    
                }
            }
            // model::cout<<"eleNodeList size = "<<eleNodeList.size()<<std::endl;
            assert(index==2 && "read inp file failed ");
            // make node-element list from element-node list, ID in ABAQUS
            for(const auto& iter : eleNodeList)
            {
                // auto eleID = iter.first;
                std::vector<int> eleID(1,iter.first);
                for(const auto& eleNode : iter.second)  // a vector store nodeID of a element
                {
                    auto insert = nodeEleList.emplace(eleNode,eleID);
                    if(!insert.second)  // this node is already existed
                    {
                        insert.first->second.push_back(iter.first);
                    }
                }
            }
            model::cout<<"          Making node-element list "<<std::flush;
            std::ofstream outfile("./AbaqusInput/0nodeEleList.txt",std::ios::out);
            assert(nodeEleList.size()==nodeSize && "make node-element list failed! ");
            for(const auto& iter : nodeEleList)
            {
                outfile<<iter.first<<" ";   // nodeID
                for (const auto& nodalEleID : iter.second)
                {
                    outfile<<" "<<nodalEleID<<" ";      // elementID assosiate with a node
                }
                outfile<<std::endl;
                
            }
            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
        }

        /**********************************************************************/
       void  initNodes()
        {// get node info from DD, position and global ID
            const auto t0= std::chrono::system_clock::now();
            model::cout<<"          initializing Nodes"<<std::flush;
            Nodes.clear();
            for (const auto& node : fe_ptr->nodes() )
            {
                // Position temp;
                VectorDim temp(node.P0);
                // for (int d=0;d<dim;d++)
                // {
                //     temp[d]=node.P0(d);
                // }
                // Nodes.emplace(temp,node.gID);
                Nodes.emplace(node.gID,temp);
            }
            // test output
            // for(const auto& iter : Nodes)
            // {// output ABAQUS nodeID & DD nodeID
            //     std::cout<<iter.first<<"  "<<iter.second<<std::endl;
            // }    

            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
        }

        /**********************************************************************/
        mapContainer& getNodes()
        {
            return Nodes;
        }


        /**********************************************************************/
        void initElements()  
        {// get element from DD, barycenter position and element ID
            const auto t0= std::chrono::system_clock::now();
            model::cout<<"          initializing Elements"<<std::flush;
            int elementID(0);
            for (const auto& element : fe_ptr->elements() )
            {
                std::vector<int> temp; 
                for (int d=0;d<dim+1;d++)
                {
                    temp.push_back(element.second.node(d).gID);
                }
                Elements.emplace(elementID,temp);
                elementID++;
            }
            // test output
            // for(const auto& iter : Elements)
            // {// output ABAQUS nodeID & DD nodeID
            //     std::cout<<iter.first<<"  "<<iter.second[0]<<"  "<<iter.second[1]<<"  "<<iter.second[2]<<"  "<<iter.second[3]<<std::endl;
            // }  
            
            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
        }

        /**********************************************************************/
        listContainer& getElements()
        {
            return Elements;
        }

        /**********************************************************************/
        void makeNodeMap()
        {
            initNodes();
            const auto t0= std::chrono::system_clock::now();
            model::cout<<"          Making node map  "<<std::flush;

// #ifdef _OPENMP  // parallelize loop over Node
// #pragma omp parallel for
//             for (unsigned int k=0;k<AbaqusNodes().size();++k)
//             {
//                 auto snIter(AbaqusNodes().begin());
//                 std::advance(snIter,k);
//                 auto abaNode(*snIter);

//                 int abaID(-1);
//                 int ddID(-1);
//                 abaID=abaNode.first;
//                 for (const auto& ddNode : getNodes() )
//                 {
//                     bool equal(false);
//                     equal=( fabs((ddNode.second-abaNode.second).norm() < 0.001 ));
//                     if (equal)
//                     {
//                         ddID = ddNode.first;
//                         break;
//                     }
//                 }
//                 assert(ddID!=-1 && "DD nodeID wrong"); 
//                 NodeMap().emplace(abaID,ddID);
//             }
// #else
            for (const auto& abaNode : AbaqusNodes() )
            {
                int abaID(-1);
                int ddID(-1);
                abaID=abaNode.first;

                // 2021-5-14, lizt
                for (const auto& ddNode : getNodes() )
                {
                    bool equal(false);
                    equal=( fabs((ddNode.second-abaNode.second).norm() < 0.001 ));
                    if (equal)
                    {
                        ddID = ddNode.first;
                        break;
                    }
                }
                // model::cout<<abaID;
                // model::cout<<"ddID = "<<ddID<<std::endl;
                // model::cout<<"abaID = "<<abaID<<std::endl;
                assert(ddID!=-1 && "DD nodeID wrong"); 
                NodeMap().emplace(abaID,ddID);
            }
// #endif
            std::ofstream outfile("./AbaqusInput/nodeMap.txt",std::ios::out);
            for(const auto& iter : NodeMap())
            {// output ABAQUS nodeID & DD nodeID
                outfile<<iter.first<<"  "<<iter.second<<std::endl;
            }            
            
            //model::cout<<"Make node map success "<<std::endl;
            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;

            // std::vector<int> abaID;
            // std::vector<Position> aba;
            // for (const auto& abaNode : AbaqusNodes() )
            // {
            //     abaID.push_back(abaNode.second);
            //     aba.push_back(abaNode.first);
            // }
            // std::vector<int> ddID;
            // std::vector<Position> dd;
            // for (const auto& ddNode : getNodes() )
            // {
            //     ddID.push_back(ddNode.second);
            //     dd.push_back(ddNode.first);
            // }
            // for (int n=0;n<nodeSize;n++)
            // {
            //     nodeMap.emplace(abaID[n],ddID[n]);
            //     model::cout<<"abaID = "<<abaID[n]<<std::endl;
            //     model::cout<<"ddID = "<<ddID[n]<<std::endl;
            //     for(int d=0;d<dim;d++)
            //     {
            //         double ep(aba[n][d]-dd[n][d]);
            //         model::cout<<ep<<std::endl;
                    
            //         assert( (aba[n][d]-dd[n][d])<0.00001 &&"node position is different");
            //         //model::cout<<"aba = "<<aba[n][d]<<std::ends;
            //     }
            //     // for(int d=0;d<dim; d++)
            //     // {
            //     //     model::cout<<"dd = "<<dd[n][d]<<std::ends;
            //     // }
            //     //assert(aba[n]==dd[n] && "node position is different");
            //     model::cout<<std::endl<<n<<std::endl;
        }




        /**********************************************************************/
        void makeElementMap()
        {
            initElements();
            const auto t0= std::chrono::system_clock::now();
            model::cout<<"          Making element map"<<std::flush;
            auto& searchElement=getElements();

// #ifdef _OPENMP  // parallelize loop over Element
// #pragma omp parallel for
//             for (unsigned int k=0;k<EleNodeList().size();++k)
//             {
//                 auto snIter(EleNodeList().begin());
//                 std::advance(snIter,k);
//                 auto abaEle(*snIter);

//                 int abaID(-1);
//                 int ddID(-1);
//                 abaID=abaEle.first;

//                 std::set<int> abaEleNode;
//                 for (const auto& abaEleNodeID : abaEle.second)
//                 {
//                     abaEleNode.insert(NodeMap()[abaEleNodeID]); // node ID in DD of a C3D8 element
//                 }
//                 assert(abaEleNode.size()==8 && "ABAQUS element node list wrong"); 

//                 for (const auto& ddEle : searchElement )
//                 {                    
//                     bool matchSuccess(false);
//                     for ( auto ddEleNodeID : ddEle.second)
//                     {
//                         matchSuccess = !( abaEleNode.find(ddEleNodeID)==abaEleNode.end() );
//                         if (!matchSuccess)
//                         {
//                             break;
//                         }
//                     }
//                     if (matchSuccess)
//                     {
//                         ddID = ddEle.first;
//                         auto insert = ElementMap().emplace(ddID,abaID); 
//                         assert(insert.second && "element map insert failed");
//                         // searchElement.erase(ddID);
//                     }

//                 }
//                 assert(ddID!=-1 && "ddID wrong"); 
//             }
                
// #else
            for ( const auto& abaEle : EleNodeList() )
            {
                int abaID(-1);
                int ddID(-1);
                abaID=abaEle.first;

                std::set<int> abaEleNode;
                for (const auto& abaEleNodeID : abaEle.second)
                {
                    abaEleNode.insert(NodeMap()[abaEleNodeID]); // node ID in DD of a C3D8 element
                }
                assert(abaEleNode.size()==8 && "ABAQUS element node list wrong"); 

                for (const auto& ddEle : searchElement )
                {                    
                    bool matchSuccess(false);
                    for ( auto ddEleNodeID : ddEle.second)
                    {
                        matchSuccess = !( abaEleNode.find(ddEleNodeID)==abaEleNode.end() );
                        if (!matchSuccess)
                        {
                            break;
                        }
                    }
                    if (matchSuccess)
                    {
                        ddID = ddEle.first;
                        auto insert = ElementMap().emplace(ddID,abaID); 
                        assert(insert.second && "element map insert failed");
                        // searchElement.erase(ddID);
                    }

                }
                assert(ddID!=-1 && "ddID wrong"); 
            }
// #endif
            std::ofstream outfile("./AbaqusInput/elementMap.txt",std::ios::out);
            for(const auto& iter : ElementMap())
            {// output DD elementID & ABAQUS elementID
                outfile<<iter.first<<"  "<<iter.second<<std::endl;
            }     
                // outfile<<ddID<<"  "<<abaID<<std::endl;          
            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
            //model::cout<<"make element map success "<<std::endl;
            
        }

    

        /**********************************************************************/
        mapContainer& AbaqusNodes()
        {
            return abaqusNodes;
        }

        /**********************************************************************/
        // mapContainer& AbaqusElements()
        // {
        //     return abaqusElements;
        // }

        /**********************************************************************/
        std::map<int,int>& NodeMap()
        {
            return nodeMap;
        }

        /**********************************************************************/
        std::map<int,int>& ElementMap()
        {
            return elementMap;
        }

        /**********************************************************************/
        listContainer& EleNodeList()
        {
            return eleNodeList;
        }    

        /**********************************************************************/
        listContainer& NodeEleList()
        {
            return nodeEleList;
        }    



    };
    

    

    
}
#endif

