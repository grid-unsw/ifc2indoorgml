#ifndef INDOORGML_H
#define INDOORGML_H

#include "LCC_SpecialOps.h"

extern std::map<str, str> ori3CellIDs2IndoorGMLIds, ori3CellIDs2IndoorGMLIds;
extern std::map<str, Dart_handle> cellspace_dart, cellboundary_dart;

/////////////////////////////////////////////////////////////////////////
///                           IndoorGML1.1                            ///
/////////////////////////////////////////////////////////////////////////
namespace IndoorGML
{

    struct externalReferenceType
    {
        std::string informationSystem,
        name,
        uri;
    };

    class CellSpace;
    class CellSpaceBoundary;

    class State
    {
    public:
        std::string description,
                    name;

        State() : Id("St_" + LCCtools::generate_unique_ID())
        {}
        State(std::string id) : Id(id)
        {}

        std::string *getId(){return &Id;}
        std::string *getDuality(){return &Duality;}
        std::map<std::string, bool, LCCtools::cmp_string> *getConnects(){return &Connects;}
        Point_3 *getGeometry(){return &Geometry;}

        void setId(std::string id)
        { Id = id; }
        void setDuality(std::string dlt)
        { Duality = dlt; }
        void setConnects(std::map<std::string, bool, LCCtools::cmp_string> &cns)
        { Connects = cns; }
        void setGeometry(Point_3 pt)
        { Geometry = pt; }

        void addConnects(std::string cn)
        { Connects[cn] = true; }

    private:
        std::string Id;
        std::string Duality;
        std::map<std::string, bool, LCCtools::cmp_string> Connects; // The keys correspond to "Transition" Ids
        Point_3 Geometry;
    };



    class Transition
    {
    public:
        std::string description,
                    name;

        Transition() : Id("Ts_" + LCCtools::generate_unique_ID()),
                       Weight("1.0")
        {}
        Transition(std::string id, std::string wgt) : Id(id), Weight(wgt)
        {}

        std::string *getId(){return &Id;}
        std::string *getWeight(){return &Weight;}
        std::map<std::string, bool, LCCtools::cmp_string> *getConnects(){return &Connects;}
        std::string *getDuality(){return &Duality;}
        std::vector<Point_3> *getGeometry(){return &Geometry;}

        void setId(std::string id)
        { Id = id; }
        void setWeight( std::string w )
        { Weight = w; }
        void setDuality(std::string dlt)
        { Duality = dlt; }
        void setGeometry(std::vector<Point_3>& geom)
        { Geometry = geom; }

        void addConnects(std::string cn)
        { Connects[cn] = true; }
        void addGeometry( Point_3 pt )
        { Geometry.push_back( pt ); }


    private:
        std::string Id,
                    Weight;
        std::map<std::string, bool, LCCtools::cmp_string> Connects; // Pair of "State" Ids
        std::string Duality;
        std::vector<Point_3> Geometry; //vector because of the potential intermediate vertices of a path
    };



    class CellSpace
    {
    public:

        std::string className, // Tells whether the CellSpace is a General, Transfer, Transition, Connection or AnchorSpace
                    description, // Changes in the 1.0.3 to include information about the storey
                    name, // And name of the storey?
                    naviclass,
                    function,
                    usage;
        bool isNavigable; // Tells if the CellSpace is navigable or not

        CellSpace() : Id("CS_" + LCCtools::generate_unique_ID()),
                      isNavigable(true),
                      className("core:CellSpace"),
                      CellSpaceGeometry(LCC::null_handle)
        {}
        std::string *getId(){return &Id;}
        Dart_handle *getCellSpaceGeometry(){return &CellSpaceGeometry;}
        std::string *getDuality(){return &Duality;}
        std::map<std::string, bool, LCCtools::cmp_string> *getPartialboundedBy(){return &PartialboundedBy;}
        std::string *getLevel(){return &level;}

        void setId(std::string id)
        { Id = id; }
        void setCellSpaceGeometry(Dart_handle d)
        {
            CellSpaceGeometry = d;
            if( cellspace_dart.find(*getId()) == cellspace_dart.end() )
                cellspace_dart[ *getId() ] = d;
        }
        void setDuality(std::string dlt)
        { Duality = dlt; }
        void setPartialboundedBy(std::map<std::string, bool, LCCtools::cmp_string> res)
        { PartialboundedBy = res; }
        void setLevel(std::string l)
        { level = l; }

    private:
        std::string Id;
        // Dart of a 3-cell
        Dart_handle CellSpaceGeometry;
        std::string Duality;
        externalReferenceType ExternalReference;
        std::map<std::string, bool, LCCtools::cmp_string> PartialboundedBy; // The keys correspond to CellSpaceBoundary Ids
        std::string level; // Since 1.1
    };


    class CellSpaceBoundary
    {
    public:

        std::string className, // Tells whether the CellSpace is a General, Transfer, Transition, Connection or AnchorSpace
                    description, // Changes in the 1.0.3 to include information about the storey
                    name, // And name of the storey?
                    naviclass,
                    function,
                    usage;

        CellSpaceBoundary() :   Id("CSB_" + LCCtools::generate_unique_ID()),
                                className("core:CellSpaceBoundary")
        {}
        std::string *getId(){return &Id;}
        Dart_handle *getCellSpaceBoundaryGeometry(){return &CellSpaceBoundaryGeometry;}
        std::string *getDuality(){return &Duality;}
//        std::map<std::string, bool, LCCtools::cmp_string> *getBounds(){return &Bounds;}

        void setId(std::string id)
        { Id = id; }
        void setCellSpaceBoundaryGeometry(Dart_handle d)
        {
            CellSpaceBoundaryGeometry = d;
            cellboundary_dart[ *getId() ] = d;
        }
        void setDuality(std::string dlt)
        { Duality = dlt; }
//        void setBounds(std::map<std::string, bool, LCCtools::cmp_string> res)
//        { Bounds = res; }

    private:
        std::string Id;
        // Dart of a 2-cell
        Dart_handle CellSpaceBoundaryGeometry;
        std::string Duality;
        externalReferenceType ExternalReference;
//        std::map<std::string, bool, LCCtools::cmp_string> Bounds; // The keys correspond to CellSpace Ids
    };


    class SpaceLayer
    {

    public:
        SpaceLayer() : Id("SL_" + LCCtools::generate_unique_ID())
        {}
        SpaceLayer(std::string id) : Id(id)
        {}

        std::string *getId(){return &Id;}
        std::map<std::string, State, LCCtools::cmp_string> *getNodes(){return &nodes;}
        std::map<std::string, Transition, LCCtools::cmp_string> *getEdges(){return &edges;}

        void setId(std::string id)
        { Id = id; }
        void setEdges( std::map<std::string, Transition, LCCtools::cmp_string>& edg )
        { edges = edg; }

        void addNode( State st )
        { nodes[*(st.getId())] = st; }
        void addTransition ( Transition trs )
        { edges[*(trs.getId())] = trs; }

    private:
        std::string Id,
        Usage,
        TerminationDate,
        Function,
        CreationDate,
        Class;

        std::map<std::string, State, LCCtools::cmp_string> nodes;
        std::map<std::string, Transition, LCCtools::cmp_string> edges;
    };


    class SpaceLayers
    {
    public:
        SpaceLayers() : Id("SLs_" + LCCtools::generate_unique_ID())
        {}
        std::string *getId(){return &Id;}
        std::map<std::string, SpaceLayer, LCCtools::cmp_string> *getSpaceLayer(){return &spaceLayerMember;}

        void setId(std::string id)
        { Id = id; }
        void setSpaceLayer(std::map<std::string, SpaceLayer, LCCtools::cmp_string>& sl)
        { spaceLayerMember = sl; }

        void addSpaceLayer(SpaceLayer& sl)
        { spaceLayerMember[*(sl.getId())] = sl; }

    private:
        std::string Id;
        std::map<std::string, SpaceLayer, LCCtools::cmp_string> spaceLayerMember;
    };



    class InterLayerConnection
    {
    public:
        InterLayerConnection() : Id("ILC_" + LCCtools::generate_unique_ID())
        {}
        std::string *getId(){return &Id;}
        std::pair<State*, State*> *getInterConnects(){return &InterConnects;}
        std::pair<SpaceLayer*, SpaceLayer*> *getConnectedLayers(){return &ConnectedLayers;}

        void setId(std::string id)
        { Id = id; }

    private:
        std::string Id,
        typeOfTopoExpression,
        comment;

        std::pair<State*, State*> InterConnects;
        std::pair<SpaceLayer*, SpaceLayer*> ConnectedLayers;
    };


    class MultiLayeredGraph
    {
    public:
        MultiLayeredGraph() : Id("MLG_" + LCCtools::generate_unique_ID())
        {}
        std::string *getId(){return &Id;}
        std::map<std::string, SpaceLayers, LCCtools::cmp_string> *getSpaceLayers(){return &spaceLayers;}
        std::map<std::string, InterLayerConnection, LCCtools::cmp_string> *getInterEdges(){return &interEdges;}

        void setId(std::string id)
        { Id = id; }
        void setSpaceLayers(std::map<std::string, SpaceLayers, LCCtools::cmp_string>& sl)
        { spaceLayers = sl; }
        void setInterEdges(std::map<std::string, InterLayerConnection, LCCtools::cmp_string>& ie)
        { interEdges = ie; }


        void addSpaceLayers(SpaceLayers& sl)
        { spaceLayers[*(sl.getId())] = sl; }
        void addInterEdges(InterLayerConnection& ie)
        { interEdges[*(ie.getId())] = ie; }


    private:
        std::string Id;
        std::map<std::string, SpaceLayers, LCCtools::cmp_string> spaceLayers;
        std::map<std::string, InterLayerConnection, LCCtools::cmp_string> interEdges;
    };


    class PrimalSpaceFeatures
    {
    public:
        PrimalSpaceFeatures() : Id("PSF_" + LCCtools::generate_unique_ID()){}
        PrimalSpaceFeatures(std::string id) : Id(id){}

        std::string *getId(){return &Id;}
        std::map<std::string, CellSpace, LCCtools::cmp_string> *getCellSpaceMember(){return &cellSpaceMember;}
        std::map<std::string, CellSpaceBoundary, LCCtools::cmp_string> *getCellSpaceBoundaryMember(){return &cellSpaceBoundaryMember;}

        void setId(std::string id)
        { Id = id; }
        void setCellSpaceMember(std::map<std::string, CellSpace, LCCtools::cmp_string>& csm)
        { cellSpaceMember = csm; }
        void setCellSpaceBoundaryMember(std::map<std::string, CellSpaceBoundary, LCCtools::cmp_string>& csbm)
        { cellSpaceBoundaryMember = csbm; }

        void addCellSpaceMember(CellSpace& csm)
        { cellSpaceMember[*(csm.getId())] = csm; }
        void addCellSpaceBoundaryMember(CellSpaceBoundary& csbm)
        { cellSpaceBoundaryMember[*(csbm.getId())] = csbm; }

    private:
        std::string Id;
        std::map<std::string, CellSpace, LCCtools::cmp_string> cellSpaceMember;
        std::map<std::string, CellSpaceBoundary, LCCtools::cmp_string> cellSpaceBoundaryMember;

    };

    class IndoorFeatures
    {

    public:
        IndoorFeatures() : Id("IFt_" + LCCtools::generate_unique_ID())
        {
            cellspace_dart.clear();
            cellboundary_dart.clear();
        }

        std::string *getId(){return &Id;}
        PrimalSpaceFeatures *getprimalSpaceFeatures(){return &primalSpaceFeatures;}
        MultiLayeredGraph *getmultiLayeredGraph(){return &multiLayeredGraph;}
        std::map<std::string, std::string, LCCtools::cmp_string> *getheader(){return &header;}

        void setId(std::string id)
        { Id = id; }
        void setprimalSpaceFeatures(PrimalSpaceFeatures& psf)
        { primalSpaceFeatures = psf; }
        void setmultiLayeredGraph(MultiLayeredGraph& mlg)
        { multiLayeredGraph = mlg; }

        void addheader(std::pair<std::string, std::string> attrib)
        { header[attrib.first] = attrib.second; }


    private:
        std::string Id;
        PrimalSpaceFeatures primalSpaceFeatures;
        MultiLayeredGraph multiLayeredGraph;

        std::map<std::string, std::string, LCCtools::cmp_string> header;

    };



    /// Navigation classes
//    class NavigableSpace : CellSpace
//    {
//    public:
//        std::string function,
//                    usage;
//    };

//    class GeneralSpace : NavigableSpace
//    {};

//    class TransferSpace : NavigableSpace
//    {};

//    class TransitionSpace : TransferSpace
//    {};

//    class ConnectionSpace : TransferSpace
//    {};

//    class AnchorSpace : TransferSpace
//    {};


//    class NavigableBoundary : CellSpaceBoundary
//    {
//    public:
//        std::string function,
//                    usage;
//    };

//    class TransferBoundary : NavigableSpace
//    {};

//    class ConnectionBoundary : TransferBoundary
//    {};

//    class AnchorBoundary : TransferBoundary
//    {};


}








/////////////////////////////////////////////////////////////////////////
///                           IndoorGML2.0                            ///
/////////////////////////////////////////////////////////////////////////


namespace IndoorGML2
{

    struct externalReferenceType{
        std::string informationSystem,
        name,
        uri;
    };

    class CellSpace;
    class CellBoundary;

    class Node : public IndoorGML::State{

    public:
        Node() : IndoorGML::State("Nd_" + LCCtools::generate_unique_ID())
        {}
    };



    class Edge : public IndoorGML::Transition{

    public:
        Edge() : IndoorGML::Transition(("Ed_" + LCCtools::generate_unique_ID()), "1.0")
        {}
    };



    class CellSpace : public IndoorGML::CellSpace{

    public:
        CellSpace() : IndoorGML::CellSpace(),
                      PoI(false)
        {}

        std::map<std::string, bool, LCCtools::cmp_string> *getBoundedBy(){return getPartialboundedBy();}
//        std::string *getLevel(){return &level;}
        std::string *getName(){return &name;}
        bool getPoI(){return PoI;}

        void setBoundedBy(std::map<std::string, bool, LCCtools::cmp_string> res){ setPartialboundedBy(res); }
//        void setLevel(std::string l)
//        { level = l; }
        void setName(std::string n)
        { name = n; }
        void setPoI(bool p)
        { PoI = p; }

    private:
//        std::string level;
        std::string name;
        bool PoI;

    };


    class CellBoundary : IndoorGML::CellSpaceBoundary{

    public:
        CellBoundary() :   IndoorGML::CellSpaceBoundary(),
                           isVirtual(false){}

        bool getIsVirtual(){ return isVirtual; }

        void setIsVirtual(bool v)
        { isVirtual = v; }

    private:
        bool isVirtual;
    };



    class ThematicLayer;

    class InterLayerConnection
    {

    public:
        InterLayerConnection() : Id("ILC_" + LCCtools::generate_unique_ID())
        {}
        std::string *getId(){return &Id;}
        std::pair<Node*, Node*> *getInterConnects(){return &InterConnects;}
        std::pair<ThematicLayer*, ThematicLayer*> *getConnectedLayers(){return &ConnectedLayers;}

        void setId(std::string id)
        { Id = id; }

        std::string typeOfTopoExpression,
                    comment;
        std::pair<Node*, Node*> InterConnects;
        std::pair<ThematicLayer*, ThematicLayer*> ConnectedLayers;

    private:
        std::string Id;
    };


    class DualSpaceLayer : public IndoorGML::SpaceLayer
    {
    public:
        DualSpaceLayer() : IndoorGML::SpaceLayer(),
                           isLogical(false)
        {}

        void setCreationDate(std::string cd){creationDate = cd;}
        void setTerminationDate(std::string td){terminationDate = td;}
        void setIsLogical(bool l){isLogical = l;}

        std::string *getCreationDate(){return &creationDate;}
        std::string *getTerminationDate(){return &terminationDate;}
        bool getIsLogical(){return isLogical;}

    private:
        std::string creationDate,
                    terminationDate;
        bool isLogical;
    };


    class PrimalSpaceLayer : public IndoorGML::PrimalSpaceFeatures
    {
    public:
        PrimalSpaceLayer() : IndoorGML::PrimalSpaceFeatures(("PSL_" + LCCtools::generate_unique_ID()))
        {}

        void setCreationDate(std::string cd){creationDate = cd;}
        void setTerminationDate(std::string td){terminationDate = td;}
        void setFunction(std::string f){function = f;}

        std::string *getCreationDate(){return &creationDate;}
        std::string *getTerminationDate(){return &terminationDate;}
        std::string *getFunction(){return &function;}

    private:
        std::string creationDate,
                    terminationDate;
        std::string function;
    };



    class ThematicLayer{

    public:
        ThematicLayer() : Id("TL_" + LCCtools::generate_unique_ID()),
                          semanticExtension(false),
                          Theme("Physical")
        {}
        ThematicLayer(std::string id) : Id(id)
        {}

        bool semanticExtension;
        std::string Theme;
        // Map of all the other ThematicLayers connected to this one
        // Key = id of the connected TL, value = pointer to their ILC
        std::map<std::string, InterLayerConnection*> IL_connections;

        std::string *getId(){return &Id;}
        PrimalSpaceLayer *getprimalSpaceLayer(){return &primalSpaceLayer;}
        DualSpaceLayer *getdualSpaceLayer(){return &dualSpaceLayer;}

        void setId(std::string id)
        { Id = id; }
        void setprimalSpaceLayer(PrimalSpaceLayer& psl)
        { primalSpaceLayer = psl; }
        void setdualSpaceLayer(DualSpaceLayer& dsl)
        { dualSpaceLayer = dsl; }

    private:
        std::string Id;
        PrimalSpaceLayer primalSpaceLayer;
        DualSpaceLayer dualSpaceLayer;
    };


    class IndoorFeatures
    {

    public:
        IndoorFeatures() : Id("InFt_" + LCCtools::generate_unique_ID())
        {
            cellspace_dart.clear();
            cellboundary_dart.clear();
        }

        std::string *getId(){return &Id;}
        std::map<std::string, ThematicLayer> *getThematicLayerMap(){return &thematicLayerMap;}
        std::map<std::string, InterLayerConnection> *getInterLayerConnectionMap(){return &interLayerConnectionMap;}
        std::map<std::string, std::string, LCCtools::cmp_string> *getheader(){return &header;}

        void setId(std::string id)
        { Id = id; }

        void setThematicLayerMap( std::map<std::string, ThematicLayer>& tlm)
        { thematicLayerMap = tlm; }
        void setInterLayerConnectionMap(std::map<std::string, InterLayerConnection>& ilcm)
        { interLayerConnectionMap = ilcm; }

        void addThematicLayer (ThematicLayer& tl)
        { thematicLayerMap[ *tl.getId() ] = tl; }
        void addInterLayerConnection ( InterLayerConnection& ilc )
        { interLayerConnectionMap[ *ilc.getId() ] = ilc; }

        void addheader(std::pair<std::string, std::string> attrib)
        { header[attrib.first] = attrib.second; }


    private:
        std::string Id;
        std::map<std::string, ThematicLayer> thematicLayerMap;
        std::map<std::string, InterLayerConnection> interLayerConnectionMap;

        std::map<std::string, std::string, LCCtools::cmp_string> header;

    };



    /// Navigation classes
//    class NavigableSpace : CellSpace
//    {};

//    class GeneralSpace : NavigableSpace
//    {};

//    class TransferSpace : NavigableSpace
//    {};

//    class EdgeSpace : TransferSpace
//    {};

//    class ConnectionSpace : TransferSpace
//    {};

//    class AnchorSpace : TransferSpace
//    {};


//    class NavigableBoundary : CellBoundary
//    {
//    public:
//        std::string function,
//                    usage;
//    };

//    class TransferBoundary : NavigableSpace
//    {};

//    class ConnectionBoundary : TransferBoundary
//    {};

//    class AnchorBoundary : TransferBoundary
//    {};


}


#endif
