/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, development version     *
*                (c) 2006-2017 INRIA, USTL, UJF, CNRS, MGH                    *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/



// File automatically generated by "generateTypedef"


#ifndef SOFA_TYPEDEF_Mass_float_H
#define SOFA_TYPEDEF_Mass_float_H

//Default files containing the declaration of the vector type
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/Mat.h>




#include <SofaBaseMechanics/DiagonalMass.h>
#include <SofaGeneralSimpleFem/HexahedralFEMForceFieldAndMass.h>
#include <SofaNonUniformFem/HexahedronCompositeFEMForceFieldAndMass.h>
#include <SofaGeneralSimpleFem/HexahedronFEMForceFieldAndMass.h>
#include <SofaMiscForceField/MatrixMass.h>
#include <SofaMiscForceField/MeshMatrixMass.h>
#include <SofaNonUniformFem/NonUniformHexahedralFEMForceFieldAndMass.h>
#include <SofaNonUniformFem/NonUniformHexahedronFEMForceFieldAndMass.h>
#include <SofaBaseMechanics/UniformMass.h>



//---------------------------------------------------------------------------------------------
//Typedef for DiagonalMass
typedef sofa::component::mass::DiagonalMass<sofa::defaulttype::StdRigidTypes<2, float>, sofa::defaulttype::RigidMass<2, float> > DiagonalMassRigid2f;
typedef sofa::component::mass::DiagonalMass<sofa::defaulttype::StdRigidTypes<3, float>, sofa::defaulttype::RigidMass<3, float> > DiagonalMassRigid3f;
typedef sofa::component::mass::DiagonalMass<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<1, float>, sofa::defaulttype::Vec<1, float>, float>, float> DiagonalMass1f;
typedef sofa::component::mass::DiagonalMass<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<2, float>, sofa::defaulttype::Vec<2, float>, float>, float> DiagonalMass2f;
typedef sofa::component::mass::DiagonalMass<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float>, float> DiagonalMass3f;



//---------------------------------------------------------------------------------------------
//Typedef for HexahedralFEMForceFieldAndMass
typedef sofa::component::forcefield::HexahedralFEMForceFieldAndMass<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > HexahedralFEMForceFieldAndMass3f;



//---------------------------------------------------------------------------------------------
//Typedef for HexahedronCompositeFEMForceFieldAndMass
typedef sofa::component::forcefield::HexahedronCompositeFEMForceFieldAndMass<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > HexahedronCompositeFEMForceFieldAndMass3f;



//---------------------------------------------------------------------------------------------
//Typedef for HexahedronFEMForceFieldAndMass
typedef sofa::component::forcefield::HexahedronFEMForceFieldAndMass<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > HexahedronFEMForceFieldAndMass3f;



//---------------------------------------------------------------------------------------------
//Typedef for MatrixMass
typedef sofa::component::mass::MatrixMass<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<1, float>, sofa::defaulttype::Vec<1, float>, float>, sofa::defaulttype::Mat<1, 1, float> > MatrixMass1f;
typedef sofa::component::mass::MatrixMass<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<2, float>, sofa::defaulttype::Vec<2, float>, float>, sofa::defaulttype::Mat<2, 2, float> > MatrixMass2f;
typedef sofa::component::mass::MatrixMass<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float>, sofa::defaulttype::Mat<3, 3, float> > MatrixMass3f;



//---------------------------------------------------------------------------------------------
//Typedef for MeshMatrixMass
typedef sofa::component::mass::MeshMatrixMass<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<1, float>, sofa::defaulttype::Vec<1, float>, float>, float> MeshMatrixMass1f;
typedef sofa::component::mass::MeshMatrixMass<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<2, float>, sofa::defaulttype::Vec<2, float>, float>, float> MeshMatrixMass2f;
typedef sofa::component::mass::MeshMatrixMass<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float>, float> MeshMatrixMass3f;



//---------------------------------------------------------------------------------------------
//Typedef for NonUniformHexahedralFEMForceFieldAndMass
typedef sofa::component::forcefield::NonUniformHexahedralFEMForceFieldAndMass<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > NonUniformHexahedralFEMForceFieldAndMass3f;



//---------------------------------------------------------------------------------------------
//Typedef for NonUniformHexahedronFEMForceFieldAndMass
typedef sofa::component::forcefield::NonUniformHexahedronFEMForceFieldAndMass<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > NonUniformHexahedronFEMForceFieldAndMass3f;



//---------------------------------------------------------------------------------------------
//Typedef for UniformMass
typedef sofa::component::mass::UniformMass<sofa::defaulttype::StdRigidTypes<2, float>, sofa::defaulttype::RigidMass<2, float> > UniformMassRigid2f;
typedef sofa::component::mass::UniformMass<sofa::defaulttype::StdRigidTypes<3, float>, sofa::defaulttype::RigidMass<3, float> > UniformMassRigid3f;
typedef sofa::component::mass::UniformMass<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<1, float>, sofa::defaulttype::Vec<1, float>, float>, float> UniformMass1f;
typedef sofa::component::mass::UniformMass<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<2, float>, sofa::defaulttype::Vec<2, float>, float>, float> UniformMass2f;
typedef sofa::component::mass::UniformMass<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float>, float> UniformMass3f;
typedef sofa::component::mass::UniformMass<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<6, float>, sofa::defaulttype::Vec<6, float>, float>, float> UniformMass6f;





#ifdef SOFA_FLOAT
typedef DiagonalMassRigid2f DiagonalMassRigid2;
typedef DiagonalMassRigid3f DiagonalMassRigid3;
typedef DiagonalMass1f DiagonalMass1;
typedef DiagonalMass2f DiagonalMass2;
typedef DiagonalMass3f DiagonalMass3;
typedef HexahedralFEMForceFieldAndMass3f HexahedralFEMForceFieldAndMass3;
typedef HexahedronCompositeFEMForceFieldAndMass3f HexahedronCompositeFEMForceFieldAndMass3;
typedef HexahedronFEMForceFieldAndMass3f HexahedronFEMForceFieldAndMass3;
typedef MatrixMass1f MatrixMass1;
typedef MatrixMass2f MatrixMass2;
typedef MatrixMass3f MatrixMass3;
typedef MeshMatrixMass1f MeshMatrixMass1;
typedef MeshMatrixMass2f MeshMatrixMass2;
typedef MeshMatrixMass3f MeshMatrixMass3;
typedef NonUniformHexahedralFEMForceFieldAndMass3f NonUniformHexahedralFEMForceFieldAndMass3;
typedef NonUniformHexahedronFEMForceFieldAndMass3f NonUniformHexahedronFEMForceFieldAndMass3;
typedef UniformMassRigid2f UniformMassRigid2;
typedef UniformMassRigid3f UniformMassRigid3;
typedef UniformMass1f UniformMass1;
typedef UniformMass2f UniformMass2;
typedef UniformMass3f UniformMass3;
typedef UniformMass6f UniformMass6;
#endif

#endif