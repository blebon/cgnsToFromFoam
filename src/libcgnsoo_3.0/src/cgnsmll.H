/*------------------------------------------------------------------------------
Copyright (C) 2004-2007 Hydro-Quebec

This file is part of CGNSOO

CGNSOO is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2, or (at your option) any
later version.

CGNSOO is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with CGNSOO  If not, write to the Free Software Foundation,
Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
------------------------------------------------------------------------------*/
#ifndef CGNSMLL_H
#define CGNSMLL_H

namespace CGNSOO
{

#include "cgnslib.h"

//Have to create the non-CG_ prefixed names, because the current library might
//be using e.g. CG_ZoneType_t
#if CG_BUILD_SCOPE
#define MassUnits_t CG_MassUnits_t
#define LengthUnits_t CG_LengthUnits_t
#define TimeUnits_t CG_TimeUnits_t
#define TemperatureUnits_t CG_TemperatureUnits_t
#define AngleUnits_t CG_AngleUnits_t
#define ElectricCurrentUnits_t CG_ElectricCurrentUnits_t
#define SubstanceAmountUnits_t CG_SubstanceAmountUnits_t
#define LuminousIntensityUnits_t CG_LuminousIntensityUnits_t
#define DataClass_t CG_DataClass_t
#define GridLocation_t CG_GridLocation_t
#define BCDataType_t CG_BCDataType_t
#define GridConnectivityType_t CG_GridConnectivityType_t
#define PointSetType_t CG_PointSetType_t
#define GoverningEquationsType_t CG_GoverningEquationsType_t
#define ModelType_t CG_ModelType_t
#define BCType_t CG_BCType_t
#define DataType_t CG_DataType_t
#define ElementType_t CG_ElementType_t
#define ZoneType_t CG_ZoneType_t
#define RigidGridMotionType_t CG_RigidGridMotionType_t
#define ArbitraryGridMotionType_t CG_ArbitraryGridMotionType_t
#define SimulationType_t CG_SimulationType_t
#define WallFunctionType_t CG_WallFunctionType_t
#define AreaType_t CG_AreaType_t
#define AverageInterfaceType_t CG_AverageInterfaceType_t

#define Unstructured CG_Unstructured
#define Structured CG_Structured

#define PointRange CG_PointRange
#define PointList CG_PointList
#define ElementRange CG_ElementRange
#define ElementList CG_ElementList
#define RealSingle CG_RealSingle
#define RealDouble CG_RealDouble
#define Integer CG_Integer
#endif

typedef char cgnsstring[33]; //!< Representation of a text label in CGNS

inline void CGNSfree( char* p ) { cg_free(p); }
}

#endif
