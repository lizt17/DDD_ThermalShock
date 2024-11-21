/* 2021-4-26, lizt
 *
 * store Abaqus solution
 * 
 * 
 * 
 */

#ifndef model_CoupleParameters_H_
#define model_CoupleParameters_H_

#include <TextFileParser.h>
#include <SimplicialMesh.h>

namespace model
{
    template <int dim>
    class CoupleParameters
    {
    private:
        /* data */
    public:

        const int layers; // search elements to localized plastic strain
        const double findShift;// distance in a step for finding simplex, relative with element size
        const SimplicialMesh<dim>& mesh;
        const int stepsBetweenBVPupdates;
        const int quadratureNumber; // point number in a link, -1 = use quadratureNumberPerlength
        const int searchOrder;

        CoupleParameters(const SimplicialMesh<dim>& mesh_in) :
        /* init  */ layers(TextFileParser("./inputFiles/DD.txt").readScalar<int>("layers",true))
        /* init  */,findShift(TextFileParser("./inputFiles/DD.txt").readScalar<double>("findShift",true))
        /* init  */,mesh(mesh_in)
        /* init  */,stepsBetweenBVPupdates(TextFileParser("./inputFiles/DD.txt").readScalar<int>("stepsBetweenBVPupdates",true))
        /* init  */,quadratureNumber(TextFileParser("./inputFiles/DD.txt").readScalar<int>("quadratureNumber",true))
        /* init  */,searchOrder(TextFileParser("./inputFiles/DD.txt").readScalar<int>("searchOrder",true))
        {
            model::cout<<"reading parameters for Abaqus ..."<<std::endl;
            model::cout<<"plastic strain zone thickness = "<<layers*findShift*2<<std::endl;
        }


    };
    

    
}


#endif