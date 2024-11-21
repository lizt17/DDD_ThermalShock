/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MeshRegion_H_
#define model_MeshRegion_H_

#include <set>
#include <deque>
#include <memory>
#include <assert.h>
#include <MPIcout.h>
#include <MeshRegionObserver.h>
#include <PlanarMeshFace.h>
#include <Simplex.h>
#include <TerminalColors.h>

namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    template<typename _SimplexType>
    struct MeshRegion : public std::set<const _SimplexType*>
    /*               */,public std::map<int,std::shared_ptr<PlanarMeshFace<_SimplexType::dim>>> // MeshRegionBoundary container
    {
        typedef _SimplexType SimplexType;
        static constexpr int dim=SimplexType::dim;
        typedef MeshRegion<SimplexType> MeshRegionType;
        typedef MeshRegionObserver<MeshRegionType> MeshRegionObserverType;
        typedef std::map<int,std::shared_ptr<PlanarMeshFace<dim>>> MeshFacesContainerType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
    private:
        
        /**********************************************************************/
        void buildSingleFace(const Simplex<dim,dim-1>* newStart,
                        std::shared_ptr<PlanarMeshFace<_SimplexType::dim>> newFace,
                        std::map<typename Simplex<dim,dim-1>::SimplexIDType,const Simplex<dim,dim-1>*>& allSimplices)
        {
            // model::cout<<" test = 9-7 "<<std::endl;
            allSimplices.erase(newStart->xID);
            // model::cout<<" allsimplex size =  "<<allSimplices.size()<<std::endl;
            const auto boundaryNeighbors(newStart->boundaryNeighbors());
            // model::cout<<" test = 9 "<<std::endl;
            for(const auto& neighbor : boundaryNeighbors)
            {
                const auto neighborRgnIDs(neighbor->regionIDs());
                // model::cout<<" test = 9-1 "<<std::endl;
                if(   neighborRgnIDs.find(regionID)!=neighborRgnIDs.end()  // neighbor in region
                   && allSimplices.find(neighbor->xID)!=allSimplices.end())   // neighbor found in allSimplices
                {
                    // model::cout<<" test = 9-2 "<<std::endl;
                    if((newStart->outNormal()-neighbor->outNormal()).norm()<FLT_EPSILON)
                    {// same plane
                        // model::cout<<" test = 9-3 "<<std::endl;
                        allSimplices.erase(neighbor->xID);
                        // model::cout<<" test = 9-4 "<<std::endl;
                        newFace->insert(neighbor);
                        // model::cout<<" test = 9-5 "<<std::endl;
                        buildSingleFace(neighbor,newFace,allSimplices);
                        // model::cout<<" test = 9-6 "<<std::endl;
                    }
                }
            }
            // model::cout<<" test = 10 "<<std::endl;
            faces().emplace(newFace->sID,newFace);
            // model::cout<<" test = 11 "<<std::endl;
        }
        
        /**********************************************************************/
        void buildFaces()
        {/*!Constructs and stores the external PlanarMeshFace(s) of this MeshRegion
          * This function supports non-convex regions.
          */
            const auto t0= std::chrono::system_clock::now();
            model::cout<<"MeshRegion "<<regionID<<" buiding faces "<<std::flush;


            faces().clear();
            // model::cout<<" test = 1 "<<std::endl;
            std::map<typename Simplex<dim,dim-1>::SimplexIDType,const Simplex<dim,dim-1>*> allSimplices;
            for(const auto& simplex : simplices())
            {// collect each possible boundary/region_boundary simplices
            // model::cout<<" test = 2 "<<std::endl;
                for(const auto& child : simplex->children())
                {
                    // model::cout<<" test = 3 "<<std::endl;
                    if(child->isBoundarySimplex())
                    {
                        allSimplices.emplace(child->xID,child.get());
                    }
                }
            }
            // model::cout<<"allSimplex size = "<<allSimplices.size()<<std::endl;

            while(allSimplices.size())
            {
                // model::cout<<" test = 5 "<<std::endl;
                std::shared_ptr<PlanarMeshFace<_SimplexType::dim>> newFace(new PlanarMeshFace<_SimplexType::dim>(allSimplices.begin()->second));
                // model::cout<<" test = 6 "<<std::endl;
                buildSingleFace(allSimplices.begin()->second,newFace,allSimplices);
                // model::cout<<" test = 7 "<<std::endl;
            }

            for(auto& face : faces())
            {
                
                face.second->finalize();
                
            }
            model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
        }
        
        MeshRegionObserverType& regionObserver;
        
        std::map<size_t,size_t> _parallelFaces;
        
    public:
        const int regionID;
        
        /**********************************************************************/
        MeshRegion(MeshRegionObserverType& ro,
                   const int& rID) :
        /* init */ regionObserver(ro),
        /* init */ regionID(rID)
        {
            const bool success=regionObserver.emplace(regionID,this).second;
            assert(success && "COULD NOT INSERT MeshRegion in MeshRegionObserver.");
        }
        
        /**********************************************************************/
        ~MeshRegion()
        {
            const size_t n=regionObserver.erase(regionID);
            assert(n==1 && "COULD NOT ERASE MeshRegion in MeshRegionObserver.");
        }
        
        /**********************************************************************/
        void update()
        {
            buildFaces();
        }
        
        /**********************************************************************/
        void identifyParallelFaces()
        {
            for(const auto& face1 : faces())
            {
                for(const auto& face2 : faces())
                {
                    if(face1.second.get()!=face2.second.get())
                    {
                        if(abs(face1.second->outNormal().dot(face2.second->outNormal())+1.0)<FLT_EPSILON)
                        {
                            _parallelFaces.emplace(face1.first,face2.first);
                        }
                    }
                }
            }
        }
        
        /**********************************************************************/
        const std::map<size_t,size_t>& parallelFaces() const
        {
            return _parallelFaces;
        }
        
        /**********************************************************************/
        const std::set<const _SimplexType*>& simplices() const
        {
            return *this;
        }
        
        /**********************************************************************/
        std::set<const _SimplexType*>& simplices()
        {
            return *this;
        }

        /**********************************************************************/
        const MeshFacesContainerType& faces() const
        {
            return *this;
        }
        
        /**********************************************************************/
        MeshFacesContainerType& faces()
        {
            return *this;
        }
        
        VectorDim outNormal(const std::set<size_t>& faceIDs) const
        {
            VectorDim temp(VectorDim::Zero());
            for(const auto& val : faceIDs)
            {
                temp+=faces().at(val)->outNormal();
            }
            const double tempNorm(temp.norm());
            return tempNorm>FLT_EPSILON? (temp/tempNorm).eval() : VectorDim::Zero();
        }
    };
    
}	// close namespace
#endif
