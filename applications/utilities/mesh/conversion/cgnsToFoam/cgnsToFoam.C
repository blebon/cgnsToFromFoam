/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright Hydro-Quebec - IREQ, 2008
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

  Description
  Conversion of CGNS files into Foam's mesh and fields

Authors
    Martin Beaudoin, Hydro-Quebec - IREQ, 2005
    Robert Magnan,   Hydro-Quebec - IREQ, 2007
    Bruno Santos, FSD blueCAPE Lda, 2018

\*---------------------------------------------------------------------------*/

// OpenFOAM includes
#include "argList.H"
#include "error.H"
#include "Time.H"
#include "IFstream.H"
#include "ISstream.H"
#include "volFields.H"
#include "fvMesh.H"
#include "pointFields.H"
#include "pointVolInterpolation.H"
#include "IOdictionary.H"

// cgns interface library
#include "libcgnsoo3/cgnsoo.H"
#include "libcgnsoo3/file.H"

// includes specific to this application
#include "CGNSBoundaryConditions.H"
#include "CGNSQuantity.H"
#include "CGNSElementType.H"
#include "CGNSQuantityConverter.H"
#include "ConnectivityMapper.H"
#include "SolutionConverter.H"
//#include "FoamReplaceChar.H"
#include "FoamQuantities.H"

// standard C++ library
#include <set>
#include <list>
#include <vector>
#include <map>
#include <string>
#include <cstring>

using namespace Foam;


// Taking care of rho
// Check if rho if it was specified as a command line argument
scalar getRhoValue( const argList& args, IOobject& ioobj )
{
    scalar rho = 1.0;  // Default value;

    // get rho from the command line
    if( args.optionReadIfPresent("rho", rho) )
    {
        Info << "User specified value for rho = " << rho << endl;
    }
    //Info << "Value for rho = " << rho << endl;
    return rho;
}

//------------------------------------------------------------------------------
// Subroutine to extract a list of word pairs from a list of region names
// This is used to parse the command line argument defining matching cyclic BC
// Returns true on error
//------------------------------------------------------------------------------
static bool parseMatchingCyclicBcNames( const std::string& args, 
list< pair<std::string,std::string> >& namepairs )
{
    // works on a modifiable c-string copy
    char* str = new char[args.length()+1];
    std::strcpy( str, args.c_str() );
    std::string tokens[2];
    int itok = 0;
    const char separators[] = { ':', ',', '\0' };
    char* psep = std::strtok( str, separators );
    while (psep != NULL)
    {
        tokens[itok++] = psep;
        if ( itok==2 )
        {
            itok = 0; // reset
            //Info << "parseMatchingCyclicBcNames: new BC pair: " << tokens[0] << " - " << tokens[1] << endl;
            namepairs.push_back( pair<std::string,std::string>(tokens[0],tokens[1]) );
        }
        psep = std::strtok( NULL, separators );
    }
    delete str;
    return (itok != 0); // error if odd number of tokens
}

//------------------------------------------------------------------------------
// Main program
//------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("CGNS file");

    argList::addBoolOption
    (
        "dryrun",
        "Just to test if the CGNS mesh makes sense"
    );
    argList::addBoolOption
    (
        "noMeshValidation",
        "Don't check if the mesh is sane"
    );
    argList::addBoolOption
    (
        "saveSolutions",
        "Also import the fields stored inside the CGNS file"
    );
    argList::addBoolOption
    (
        "separatePatches",
        "Use this option if the patches are to be imported"
    );
    argList::addBoolOption
    (
        "cfxCompatibility",
        "Deploy CFX compatibility"
    );
    argList::addBoolOption
    (
        "use2DElementsAsPatches",
        "Use 2D elements as patches"
    );
    argList::addOption
    (
        "matchCyclicBC",
        "bc1_name:bc2_name",
        "Define path of patches that should be matched as a pair of cyclics"
    );
    argList::addOption
    (
        "cyclicRotX",
        "value",
        "Cyclic rotation around X axis, in case of cyclic patches"
    );
    argList::addOption
    (
        "cyclicRotY",
        "value",
        "Cyclic rotation around Y axis, in case of cyclic patches"
    );
    argList::addOption
    (
        "cyclicRotZ",
        "value",
        "Cyclic rotation around Z axis, in case of cyclic patches"
    );
    argList::addOption
    (
        "mergeTolerance",
        "value",
        "Point merge tolerance for cyclic patches"
    );
    argList::addOption
    (
        "cyclicMatchFaceTolFactor",
        "value",
        "Face matching tolerance factor for cyclic patches"
    );
    argList::addOption
    (
        "rho",
        "value",
        "Constant density field value to be used when converting pressure field"
    );
    argList::addBoolOption
    (
        "debug",
        "For providing additional debugging information"
    );
#if DEFINED_MAP_UNKNOWN  // Work in progress
    argList::addBoolOption
    (
        "mapUnknown",
        ""
    );
    // Disabled for now. Work in progress @ Hydro-Quebec.
#endif

    argList::addNote
    (
        "To specify a rotation angle like 1/13 of 360 degree, you can use the"
        " Unix command bc to compute the value with full precision, eg:\n"
        "            -cyclicRotZ `echo \"360/13\" | bc -l`"
    );

    #include "addRegionOption.H"

    argList args(argc, argv);

    #include "createTime.H"

    word regionName = polyMesh::defaultRegion;
    fileName meshDir;

    if (args.optionReadIfPresent("region", regionName))
    {
        meshDir = regionName/polyMesh::meshSubDir;
    }
    else
    {
        meshDir = polyMesh::meshSubDir;
    }

    bool separatePatches = args.optionFound("separatePatches");
    bool cfxCompatibilityMode = args.optionFound("cfxCompatibility");
    bool use2DElementsAsPatches = args.optionFound("use2DElementsAsPatches");
    bool debug = args.optionFound("debug");

    bool mapUnknown = false;
#if DEFINED_MAP_UNKNOWN  // Work in progress
    mapUnknown = args.optionFound("mapUnknown");
#endif

    bool dryRun = args.optionFound("dryrun");
    if (dryRun)
    {
        Warning
            << "Option -dryrun was used: data will not be saved"
            << endl;
    }

    scalar mergeTolerance = 0.1;
    if ( args.optionReadIfPresent("mergeTolerance", mergeTolerance) )
    {
        Info << "MergeTolerance set to " << mergeTolerance << endl;
    }

    if( args.optionFound("rho") && !args.optionFound("saveSolutions") )
    {
        Warning
            << "Option -rho is useless since no solution will be saved (no "
            << "-saveSolutions option)"
            << endl;
    }

#if 0
    CgnsToFoamDictionary optionsDict(
        IOobject( "CgnsToFoamDict", 
        mesh.time().constant(), 
        mesh,
        IOobject::READ_IF_PRESENT, 
        IOobject::NO_WRITE )
    );
    if ( optionsDict.headerOk() )
    {
        Info << "CgnsToFoamOptions: " << endl;
        Info << "\tSplitting multiple cell into different zones : " << optionsDict.splitMixed() << endl;
        Info << "\tConversion Path is " << optionsDict.conversionDirectory() << endl;
    }
#endif
    //--------------------------------------------------

    std::string cgnsFilename(args.args()[1]);
    if ( debug )
    {
        Info << "Conversion of file " << cgnsFilename << endl;
    }
    CGNSOO::file cgnsFile( cgnsFilename, CGNSOO::file::READONLY );

    int nbases = cgnsFile.getNbBase();
    if ( nbases == 0 )
    {
        FatalErrorInFunction
            << "No Base node were found in the CGNS file" 
            << exit( FatalError );
    }
    else if ( nbases > 1 )
    {
        Warning
            << "More than one CGNS Base were found in the file." << endl
            << "Only the first one will be treated." << endl;
    }

    std::string baseName;
    int physicalDim, cellDim;
    CGNSOO::Base_t base = cgnsFile.readBase( 0, baseName, physicalDim, cellDim );
    int nbZones = base.getNbZone();
    if ( debug )
    {
        Info << "Mesh: total number of zones: " << nbZones << endl;
    }

    //--------------------------------------------------------------------------
    // Load family definitions into a map indexed by family name
    // We are only interested in family with BCs to resolve
    // boundary conditions of type FamilySpecified.
    //--------------------------------------------------------------------------
    int nfamilies = base.getNbFamily();
    map<std::string, CGNSOO::BCType_t> fammap;
    for ( int ifam=0 ; ifam<nfamilies ; ifam++ )
    {
        std::string famname;
        int ngeoref; // unused
        bool hasfbc;
        CGNSOO::Family_t fam = base.readFamily( ifam, famname, hasfbc, ngeoref );
        if ( hasfbc )
        {
            std::string fambcname; // not useful
            CGNSOO::BCType_t bctype;
            fam.readFamilyBC( fambcname, bctype );
            fammap[famname] = bctype;
        }
    }

    //--------------------------------------------------------------------------
    // Build a local to global node map table
    // Extract node coordinates
    //--------------------------------------------------------------------------
    ConnectivityMapper mapper;

    for( int indexZone=0 ; indexZone<nbZones ; indexZone++ )
    {
        // Zone
        std::string zonename;
        std::vector<int> nodesize, cellsize, bndrysize;
        CGNSOO::ZoneType_t zonetype;
        CGNSOO::Zone_t zone = base.readZone
        (
            indexZone,
            zonename,
            nodesize,
            cellsize,
            bndrysize,
            zonetype
        );

        // Read the mesh nodes
        std::vector<double> xCoordinates;
        std::vector<double> yCoordinates;
        std::vector<double> zCoordinates;

        int nbgrids = zone.getNbGridCoordinates();
        if ( nbgrids != 1 )
        {
            Warning << "CGNS file with more than one GridCoordinates"
                << " - only the first grid will be read"
                << endl;
        }

        // What type of coordinate system do we have?
        std::string coordName;
        CGNSOO::GridCoordinates_t gridcoo = zone.readGridCoordinates( 0, coordName );
        int nbcoords = gridcoo.getNbCoordinatesData();

        std::map<std::string,int> coomap;
        for ( int icoo=0 ; icoo<nbcoords ; icoo++ )
        {
            std::string     cooname;
            CGNSOO::DataType_t cootype;
            gridcoo.getCoordinatesDataInfo( icoo, cooname, cootype );
            coomap[cooname] = icoo;
        }

        if
        (
            nbcoords==3 &&
            coomap.find("CoordinateX") != coomap.end() &&
            coomap.find("CoordinateY") != coomap.end() &&
            coomap.find("CoordinateZ") != coomap.end()
        )
        {
            // Cartesian coordinate system
            CGNSOO::DataArray_t datax = gridcoo.readCoordinatesData( "CoordinateX", xCoordinates );
            CGNSOO::DataArray_t datay = gridcoo.readCoordinatesData( "CoordinateY", yCoordinates );
            CGNSOO::DataArray_t dataz = gridcoo.readCoordinatesData( "CoordinateZ", zCoordinates );
        }
        else if
        (
            nbcoords==3 &&
            coomap.find("CoordinateR") != coomap.end() &&
            coomap.find("CoordinateZ") != coomap.end() &&
            coomap.find("CoordinateTheta") != coomap.end()
        )
        {
            // lire en r,t,z et convertir en x,y,z
            std::vector<double> rCoo, tCoo;
            CGNSOO::DataArray_t datar = gridcoo.readCoordinatesData( "CoordinateR", rCoo );
            CGNSOO::DataArray_t datat = gridcoo.readCoordinatesData( "CoordinateTheta", tCoo );
            CGNSOO::DataArray_t dataz = gridcoo.readCoordinatesData( "CoordinateZ", zCoordinates );
            int ndata = rCoo.size();
            for ( int i=0 ; i<ndata ; i++ )
            {
                xCoordinates.push_back( rCoo[i]*std::cos(tCoo[i]) );
                yCoordinates.push_back( rCoo[i]*std::sin(tCoo[i]) );
            }
        }
        else
        {
            FatalError << "Unknown coordinate system ";
            for ( std::map<std::string,int>::const_iterator i  = coomap.begin() ; 
                  i != coomap.end() ; 
                  i++ )
                FatalError << " : " << (*i).first;
            FatalError << endl << exit(FatalError);
        }

        // Convert the 3 vectors of doubles into a list of points
        list<point> plist;
        int nnodes = xCoordinates.size();
        for( int i=0 ; i<nnodes ; i++)
        {
            point pt( xCoordinates[i],
            yCoordinates[i],
            zCoordinates[i] );
            plist.push_back( pt );
        }

        // Read the connectivity information
        // and let the merger handle all this mess!!
        switch ( zonetype )
        {
            case CGNSOO::Structured:
                // Connectivity is implicit ... it will be constructured internally
                mapper.addStructuredZone( nodesize[0],nodesize[1],nodesize[2],plist );
                break;
            case CGNSOO::Unstructured:
                {
                    // get all the element sections
                    int nbesections = zone.getNbElements();
                    std::vector< std::string        > sectionnames(nbesections);
                    std::vector< std::vector<int>   > connectivities(nbesections);
                    std::vector< CGNSOO::ElementType_t > etypes(nbesections);
                    for ( int ies=0 ; ies<nbesections ; ies++ )
                    {
                        //string sectionname;
                        int    start, end, nbndry;
                        bool   hasparent;
                        CGNSOO::Elements_t esection = zone.readElements( 
                            ies, sectionnames[ies], etypes[ies], start, end, nbndry, hasparent );
                        CGNSOO::DataArray_t connecda = esection.readConnectivity( connectivities[ies] );
                    }
                    mapper.addUnstructuredZone( nodesize[0], cellsize[0], plist,
                    sectionnames, etypes, connectivities );
                }
                break;
            default:
                FatalError << "Unrecognized zone type"
                    << " - only Structured and Unstructured are accepted"
                    << exit(FatalError);
                break;
        }
    }

    // Merge identical points at the interfaces
    // Normally, for an unstructured mesh, there should not be any change
    // For a structured one, the points at the boundary must be merged.
    if ( debug )
    {
        Info << "Mesh: total number of nodes before merge: " << mapper.getTotalNodes() << endl;
    }
    mapper.merge( mergeTolerance );

    // Validation
    int nCreatedCells = mapper.getTotalCells();
    if ( debug )
    {
        Info << "Mesh: total number of nodes after merge: " << mapper.getTotalNodes() << endl;
        Info << "Mesh: total number of cells after merge: " << nCreatedCells << endl;
    }

    //--------------------------------------------------------------------------
    // Read-in the boundary conditions
    //--------------------------------------------------------------------------

    CGNSBoundaryConditions bc_merger( mapper, fammap, debug );
    if ( use2DElementsAsPatches )
    {
        bc_merger.addPatchesFromElements( CGNSOO::BCTypeNull  );
    }
    else
    {
        for( int indexZone=0 ; indexZone<nbZones ; indexZone++ )
        {
            std::string zonename;
            std::vector<int> nodesize, cellsize, bndrysize;
            CGNSOO::ZoneType_t zonetype;
            CGNSOO::Zone_t    zone = base.readZone( indexZone, 
            zonename, 
            nodesize, cellsize, bndrysize, 
            zonetype );
            CGNSOO::ZoneBC_t zonebc = zone.readZoneBC();
            int nbBoco = zonebc.getNbBoundaryConditions();

            for( int indexBC=0 ; indexBC < nbBoco ; indexBC++ )
            {
                std::string         bcname;
                CGNSOO::BCType_t       bctype;
                CGNSOO::PointSetType_t psettype;
                CGNSOO::BC_t bc = zonebc.readBC( indexBC, bcname, bctype, psettype );

                bc_merger.addBoundaryPatch( base, indexZone, zonetype, bc, bcname, bctype, psettype );
            }
        }
    }
    bc_merger.buildPatches( separatePatches );

    // Lets deal with cyclic BC ....
    if (args.optionFound("matchCyclicBC"))
    {
        scalar cyclicRotX = 0;
        scalar cyclicRotY = 0;
        scalar cyclicRotZ = 0;
        scalar cyclicMatchFaceTolFactor = 0.1;  // Default value... one tenth of the smallest edge found in the mesh

        if ( debug )
            Info << "Processing cyclic boundary conditions" << endl;

        // Divide the argument to matchCyclicBC into pairs of names to
        // be matched one to another
        std::string argCyclicBC = args.options()["matchCyclicBC"];
        //std::string s_argCyclicBC = FoamReplaceChar( argCyclicBC, ' ', '_' );
        std::string s_argCyclicBC = argCyclicBC;
        list< pair<std::string,std::string> > l_matchingCyclicBcNames;
        if ( parseMatchingCyclicBcNames(s_argCyclicBC,l_matchingCyclicBcNames) )
        {
            FatalError  << "Unmatched bcname for cyclic boundaries"
                << " - number of names must be even."
                << exit(FatalError);
        }

        if( l_matchingCyclicBcNames.size() > 0 )
        {
            Foam::vector rotAxis(0,0,1);
            scalar rotAngle = 0;
            if( args.optionFound("cyclicRotX") )
            {
                cyclicRotX = readScalar(IStringStream(args.options()["cyclicRotX"])());
                rotAxis = Foam::vector(1,0,0);
                rotAngle = cyclicRotX;
            }
            if( args.optionFound("cyclicRotY") )
            {
                cyclicRotY = readScalar(IStringStream(args.options()["cyclicRotY"])());
                rotAxis = Foam::vector(0,1,0);
                rotAngle = cyclicRotY;
            }
            if( args.optionFound("cyclicRotZ") )
            {
                cyclicRotZ = readScalar(IStringStream(args.options()["cyclicRotZ"])());
                rotAxis = Foam::vector(0,0,1);
                rotAngle = cyclicRotZ;
            }
            if( args.optionFound("cyclicMatchFaceTolFactor") )
                cyclicMatchFaceTolFactor = readScalar(IStringStream(args.options()["cyclicMatchFaceTolFactor"])());

            // Validation
            if( cyclicRotX == 0 && cyclicRotY == 0 && cyclicRotZ == 0 )
            {
                FatalError << "Option -cyclicRot[XYZ] is missing. "
                    << "Please specify a rotation angle for the periodic cyclic BC." 
                    << exit(FatalError);
            }


            for ( list< pair<std::string,std::string> >::iterator p  = l_matchingCyclicBcNames.begin() ;
                  p != l_matchingCyclicBcNames.end() ;
                  p++ )
            {
                const std::string& firstPatchName  = p->first;
                const std::string& secondPatchName = p->second;

                int retcode = bc_merger.computeCyclicBC( firstPatchName, secondPatchName, rotAxis, rotAngle, cyclicMatchFaceTolFactor );
                if ( retcode < 0 )
                {
                    switch( retcode )
                    {
                        case -1:
                            FatalError << "Option -matchCyclicBC: Unknown BC name: "
                                << firstPatchName << endl;
                            break;
                        case -2:
                            FatalError << "Option -matchCyclicBC: Unknown BC name: "  
                                << secondPatchName << endl;
                            break;
                    }
                    FatalError.exit(1);
                }
            }
        }
    }

    //--------------------------------------------------------------------------
    // Output in the OpenFOAM format
    //--------------------------------------------------------------------------
    if ( debug )
    {
        Info << "Creation of Foam poly mesh" << endl;
    }

    polyMesh* pShapeMesh = bc_merger.buildFoamMesh( runTime );

    if (args.optionFound("noMeshValidation"))
    {
        if ( debug )
            Info << "Option -noMeshValidation activated: disabling the mesh validation with checkMesh" << endl;
    }
    else
    {
        if ( debug )
            Info << "Validation of the mesh" << endl;
        pShapeMesh->checkMesh();
    }

    if (!dryRun)
    {
            Info << "Output of mesh and boundary conditions" << endl;
        IOstream::defaultPrecision(10); // Set the precision of the points data to 10
        pShapeMesh->write();
    }

    // Save CGNS solutions as initial solution? 
    if ( args.optionFound("saveSolutions") )
    {
        fvMesh cell_mesh(*pShapeMesh);

        scalar rho = getRhoValue( args, runTime );

        CGNSQuantityConverter* qConverter = (cfxCompatibilityMode) 
            ? new CGNSQuantityConverter_CFXcompatibility()
            : new CGNSQuantityConverter();
        SolutionConverter solConverter( *pShapeMesh, base, *qConverter );
        solConverter.buildAndWriteFoamFields( mapper, base, rho, runTime, dryRun, mapUnknown);
        delete qConverter;
    }

    if ( debug )
        Info << "Please validate Boundary Conditions by running FoamX." << endl;

    return 0;
}

#ifdef MISSING_POINTVOLINTERPOLATION
# include "pointVolInterpolation.C"
#endif
