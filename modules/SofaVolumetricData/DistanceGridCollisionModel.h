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
#ifndef SOFA_COMPONENT_COLLISION_DISTANCEGRIDCOLLISIONMODEL_H
#define SOFA_COMPONENT_COLLISION_DISTANCEGRIDCOLLISIONMODEL_H
#include "config.h"

#include <sofa/core/CollisionModel.h>
#include <SofaVolumetricData/DistanceGrid.h>
#include <SofaBaseMechanics/MechanicalObject.h>
#include <SofaMeshCollision/RigidContactMapper.h>
#include <sofa/defaulttype/Vec3Types.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <SofaBaseTopology/RegularGridTopology.h>
#include <sofa/core/objectmodel/DataFileName.h>
#include <SofaBaseTopology/SparseGridTopology.h>
#include <SofaMeshCollision/BarycentricContactMapper.h>


namespace sofa
{

namespace component
{

namespace collision
{

namespace _distancegridcollisionmodel_
{

using sofa::component::mapping::IdentityMapping;
using sofa::component::topology::RegularGridTopology;
using sofa::component::topology::SparseGridTopology;
using sofa::component::container::DistanceGrid;
using sofa::core::topology::BaseMeshTopology;
using sofa::core::TCollisionElementIterator;
using sofa::core::CollisionElementIterator;
using sofa::core::behavior::MechanicalState;
using sofa::core::CollisionModel;
using sofa::core::objectmodel::DataFileName;
using sofa::defaulttype::RigidTypes;
using sofa::defaulttype::Rigid3Types;
using sofa::defaulttype::Vec3Types;
using sofa::defaulttype::Matrix3;
using sofa::defaulttype::Matrix4;
using sofa::defaulttype::Mat3x3d;
using sofa::defaulttype::Vector3;
using sofa::defaulttype::Vector4;
using sofa::core::visual::VisualParams;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

class SOFA_VOLUMETRIC_DATA_API RigidDistanceGridCollisionModel;

class SOFA_VOLUMETRIC_DATA_API RigidDistanceGridCollisionElement : public TCollisionElementIterator<RigidDistanceGridCollisionModel>
{
public:
    RigidDistanceGridCollisionElement(RigidDistanceGridCollisionModel* model, int index);
    explicit RigidDistanceGridCollisionElement(const CollisionElementIterator& i);


    bool           isTransformed();
    const Matrix3& getRotation();
    const Vector3& getTranslation();
    bool           isFlipped();

    void           setGrid(DistanceGrid* surf);
    DistanceGrid*  getGrid();

    //// @name Previous state data
    //// Used to estimate velocity in case the distance grid itself is dynamic
    //// @{
    DistanceGrid*  getPrevGrid();
    const Matrix3& getPrevRotation();
    const Vector3& getPrevTranslation();
    double         getPrevDt();
    //// @}

    //// Set new grid and transform, keeping the old state to estimate velocity
    void setNewState(double dt,
                     DistanceGrid* grid,
                     const Matrix3& rotation,
                     const Vector3& translation);
};

class SOFA_VOLUMETRIC_DATA_API RigidDistanceGridCollisionModel : public CollisionModel
{
public:
    SOFA_CLASS(RigidDistanceGridCollisionModel, CollisionModel);

    typedef Rigid3Types InDataTypes;
    typedef Vec3Types DataTypes;
    typedef RigidDistanceGridCollisionElement Element;

public: //// Members
    /// Input data parameters
    DataFileName     fileRigidDistanceGrid;
    Data< double >   scale;
    Data< Vector3 >  translation;
    Data< Vector3 >  rotation;
    Data< double >   sampling;
    Data< helper::fixed_array<DistanceGrid::Coord,2> > box;
    Data< int >      nx;
    Data< int >      ny;
    Data< int >      nz;
    DataFileName     dumpfilename;

    Data< bool >     usePoints;
    Data< bool >     flipNormals;
    Data< bool >     showMeshPoints;
    Data< bool >     showGridPoints;
    Data< double >   showMinDist;
    Data< double >   showMaxDist;

public:
    MechanicalState<InDataTypes>* getRigidModel() { return rigid; }
    MechanicalState<InDataTypes>* getMechanicalState() { return rigid; }

    //////// virtual inheritance from Base
    void init() ;
    void draw(const core::visual::VisualParams*, int index);
    void draw(const core::visual::VisualParams* vparams);
    //////////////////////////////////////

    DistanceGrid* getGrid(int index=0)
    {
        return elems[index].grid;
    }

    bool isTransformed(int index=0) const
    {
        return elems[index].isTransformed;
    }

    const Matrix3& getRotation(int index=0) const
    {
        return elems[index].rotation;
    }

    const Vector3& getTranslation(int index=0) const
    {
        return elems[index].translation;
    }

    const Vector3& getInitTranslation() const
    {
        return translation.getValue();
    }

    const Matrix3 getInitRotation() const
    {
        SReal x = rotation.getValue()[0] * M_PI / 180;
        SReal y = rotation.getValue()[1] * M_PI / 180;
        SReal z = rotation.getValue()[2] * M_PI / 180;

        Matrix3 X(Vector3(1,0,0), Vector3(0, cos(x), -sin(x)), Vector3(0, sin(x), cos(x)));
        Matrix3 Y(Vector3(cos(y), 0, sin(y)), Vector3(0, 1, 0), Vector3(-sin(y), 0, cos(y)));
        Matrix3 Z(Vector3(cos(z), -sin(z), 0), Vector3(sin(z), cos(z), 0), Vector3(0, 0, 1));
        return X * Y * Z;
    }

    bool isFlipped() const
    {
        return flipNormals.getValue();
    }

    void setGrid(DistanceGrid* surf, int index=0);

    DistanceGrid* getPrevGrid(int index=0)
    {
        return elems[index].prevGrid;
    }
    const Matrix3& getPrevRotation(int index=0) const
    {
        return elems[index].prevRotation;
    }
    const Vector3& getPrevTranslation(int index=0) const
    {
        return elems[index].prevTranslation;
    }
    double getPrevDt(int index=0) const
    {
        return elems[index].prevDt;
    }

    //// Set new grid and transform, keeping the old state to estimate velocity
    void setNewState(int index, double dt, DistanceGrid* grid, const Matrix3& rotation, const Vector3& translation);

    //// @}

    //// Update transformation matrices from current rigid state
    void updateState();

    //// -- CollisionModel interface
    void resize(int size);

    //// Create or update the bounding volume hierarchy.
    void computeBoundingTree(int maxDepth=0);

protected:
    class ElementData
    {
    public:
        Matrix3       rotation;
        Vector3       translation;
        DistanceGrid* grid;

        //// @name Previous state data
        //// Used to estimate velocity in case the distance grid itself is dynamic
        //// @{
        DistanceGrid* prevGrid; ///< Previous grid
        Matrix3       prevRotation; ///< Previous rotation
        Vector3       prevTranslation; ///< Previous translation
        double        prevDt; ///< Time difference between previous and current state
        //// @}

        bool          isTransformed; ///< True if translation/rotation was set
        ElementData() : grid(NULL), prevGrid(NULL), prevDt(0.0), isTransformed(false) { rotation.identity(); prevRotation.identity(); }
    };

    helper::vector<ElementData>     elems;
    bool                            modified;
    MechanicalState<RigidTypes>*    rigid;

    RigidDistanceGridCollisionModel();
    ~RigidDistanceGridCollisionModel();

    void updateGrid();
};

inline RigidDistanceGridCollisionElement::RigidDistanceGridCollisionElement(RigidDistanceGridCollisionModel* model, int index)
    : core::TCollisionElementIterator<RigidDistanceGridCollisionModel>(model, index)
{}

inline RigidDistanceGridCollisionElement::RigidDistanceGridCollisionElement(const core::CollisionElementIterator& i)
    : core::TCollisionElementIterator<RigidDistanceGridCollisionModel>(static_cast<RigidDistanceGridCollisionModel*>(i.getCollisionModel()), i.getIndex())
{
}

inline DistanceGrid* RigidDistanceGridCollisionElement::getGrid() { return model->getGrid(index); }
inline void RigidDistanceGridCollisionElement::setGrid(DistanceGrid* surf) { return model->setGrid(surf, index); }

inline bool RigidDistanceGridCollisionElement::isTransformed() { return model->isTransformed(index); }
inline const Matrix3& RigidDistanceGridCollisionElement::getRotation() { return model->getRotation(index); }
inline const Vector3& RigidDistanceGridCollisionElement::getTranslation() { return model->getTranslation(index); }
inline bool RigidDistanceGridCollisionElement::isFlipped() { return model->isFlipped(); }

inline DistanceGrid* RigidDistanceGridCollisionElement::getPrevGrid() { return model->getPrevGrid(index); }
inline const Matrix3& RigidDistanceGridCollisionElement::getPrevRotation() { return model->getPrevRotation(index); }
inline const Vector3& RigidDistanceGridCollisionElement::getPrevTranslation() { return model->getPrevTranslation(index); }
inline double RigidDistanceGridCollisionElement::getPrevDt() { return model->getPrevDt(index); }

inline void RigidDistanceGridCollisionElement::setNewState(double dt, DistanceGrid* grid,
                                                           const Matrix3& rotation,
                                                           const Vector3& translation)
{
    return model->setNewState(this->getIndex(), dt, grid, rotation, translation);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

class FFDDistanceGridCollisionModel;

class FFDDistanceGridCollisionElement : public TCollisionElementIterator<FFDDistanceGridCollisionModel>
{
public:

    FFDDistanceGridCollisionElement(FFDDistanceGridCollisionModel* model, int index);
    explicit FFDDistanceGridCollisionElement(const core::CollisionElementIterator& i);

    DistanceGrid* getGrid();

    void setGrid(DistanceGrid* surf);
};

class SOFA_VOLUMETRIC_DATA_API FFDDistanceGridCollisionModel : public CollisionModel
{

public:
    SOFA_CLASS(FFDDistanceGridCollisionModel, CollisionModel);

    typedef SReal GSReal;
    typedef DistanceGrid::Coord GCoord;

    class SOFA_VOLUMETRIC_DATA_API DeformedCube
    {
    public:
        enum {C000 = 0+0+0,
                C100 = 1+0+0,
                C010 = 0+2+0,
                C110 = 1+2+0,
                C001 = 0+0+4,
                C101 = 1+0+4,
                C011 = 0+2+4,
                C111 = 1+2+4
             };


        enum {FX0 = 0+0,
                FX1 = 0+1,
                FY0 = 2+0,
                FY1 = 2+1,
                FZ0 = 4+0,
                FZ1 = 4+1
             };

        struct Point
        {
            GCoord bary; ///< Barycentric coordinates
            int index; ///< Index of corresponding point in DistanceGrid
        };

        typedef defaulttype::Vec<4,GSReal> Plane; ///< plane equation as defined by Plane.(x y z 1) = 0

        DistanceGrid* grid;
        DeformedCube() : grid(NULL) {}

        int                    elem; ///< Index of the corresponding element in the topology
        std::set<int>          neighbors; ///< Index of the neighbors (used for self-collisions)

        helper::vector<Point>  points; ///< barycentric coordinates of included points
        helper::vector<GCoord> normals; ///< normals in barycentric coordinates of included points
        GCoord                 initP0,initDP,invDP; ///< Initial corners position
        GCoord                 corners[8]; ///< Current corners position
        Plane                  faces[6]; ///< planes corresponding to the six faces (FX0,FX1,FY0,FY1,FZ0,FZ1)

        //// @name Precomputed deformation factors
        //// We have :
        ////   deform(b) = C000(1-b[0])(1-b[1])(1-b[2]) + C100(b[0])(1-b[1])(1-b[2]) + C010(1-b[0])(b[1])(1-b[2]) + C110(b[0])(b[1])(1-b[2])
        ////             + C001(1-b[0])(1-b[1])(  b[2]) + C101(b[0])(1-b[1])(  b[2]) + C011(1-b[0])(b[1])(  b[2]) + C111(b[0])(b[1])(  b[2])
        ////             = C000 + Dx b[0] + Dy b[1] + Dz b[2] + Dxy b[0]b[1] + Dxz b[0]b[2] + dyz b[1]b[2] + dxyz b[0]b[1]b[2]
        //// @{
        GCoord Dx;   ///< Dx = -C000+C100
        GCoord Dy;   ///< Dy = -C000+C010
        GCoord Dz;   ///< Dx = -C000+C001
        GCoord Dxy;  ///< Dxy = C000-C100-C010+C110 = C110-C010-Dx
        GCoord Dxz;  ///< Dxz = C000-C100-C001+C101 = C101-C001-Dx
        GCoord Dyz;  ///< Dyz = C000-C010-C001+C011 = C011-C001-Dy
        GCoord Dxyz; ///< Dxyz = - C000 + C100 + C010 - C110 + C001 - C101 - C011 + C111 = C001 - C101 - C011 + C111 - Dxy
        //// @}

        GCoord center; ///< current center;
        GSReal radius; ///< radius of enclosing sphere
        helper::vector<GCoord> deformedPoints; ///< deformed points
        helper::vector<GCoord> deformedNormals; ///< deformed normals
        bool pointsUpdated; ///< true the deformedPoints vector has been updated with the latest positions
        void updatePoints(); ///< Update the deformedPoints position if not done yet (i.e. if pointsUpdated==false)
        bool facesUpdated; ///< true the faces plane vector has been updated with the latest positions
        void updateFaces(); ///< Update the face planes if not done yet (i.e. if facesUpdated==false)

        //// Update the deformation precomputed values
        void updateDeform();

        //// Compute a plane equation given 4 corners
        Plane computePlane(int c00, int c10, int c01, int c11);

        //// Compute the barycentric coordinates of a point from its initial position
        DistanceGrid::Coord baryCoords(const GCoord& c) const
        {
            return GCoord( (c[0]-initP0[0])*invDP[0],
                    (c[1]-initP0[1])*invDP[1],
                    (c[2]-initP0[2])*invDP[2]);
        }

        //// Compute the initial position of a point from its barycentric coordinates
        GCoord initpos(const GCoord& b) const
        {
            return GCoord( initP0[0]+initDP[0]*b[0],
                    initP0[1]+initDP[1]*b[1],
                    initP0[2]+initDP[2]*b[2]);
        }

        //// Compute the deformed position of a point from its barycentric coordinates
        GCoord deform(const GCoord& b) const
        {
            return corners[C000] + Dx*b[0] + (Dy + Dxy*b[0])*b[1] + (Dz + Dxz*b[0] + (Dyz + Dxyz*b[0])*b[1])*b[2];
        }

        //// deform a direction relative to a point in barycentric coordinates
        GCoord deformDir(const GCoord& b, const GCoord& dir) const
        {
            GCoord r;
            /// dp/dx = Dx + Dxy*y + Dxz*z + Dxyz*y*z
            r  = (Dx + Dxy*b[1] + (Dxz + Dxyz*b[1])*b[2])*dir[0];
            /// dp/dy = Dy + Dxy*x + Dyz*z + Dxyz*x*z
            r += (Dy + Dxy*b[0] + (Dyz + Dxyz*b[0])*b[2])*dir[1];
            /// dp/dz = Dz + Dxz*x + Dyz*y + Dxyz*x*y
            r += (Dz + Dxz*b[0] + (Dyz + Dxyz*b[0])*b[1])*dir[2];
            return r;
        }

        //// Get the local jacobian matrix of the deformation
        Mat3x3d Jdeform(const GCoord& b) const
        {
            Mat3x3d J;
            for (int i=0; i<3; i++)
            {
                /// dp/dx = Dx + Dxy*y + Dxz*z + Dxyz*y*z
                J[i][0] = (Dx[i] + Dxy[i]*b[1] + (Dxz[i] + Dxyz[i]*b[1])*b[2]);
                /// dp/dy = Dy + Dxy*x + Dyz*z + Dxyz*x*z
                J[i][1] = (Dy[i] + Dxy[i]*b[0] + (Dyz[i] + Dxyz[i]*b[0])*b[2]);
                /// dp/dz = Dz + Dxz*x + Dyz*y + Dxyz*x*y
                J[i][2] = (Dz[i] + Dxz[i]*b[0] + (Dyz[i] + Dxyz[i]*b[0])*b[1]);
            }
            return J;
        }

        //// Compute an initial estimate to the barycentric coordinate of a point given its deformed position
        GCoord undeform0(const GCoord& p) const
        {
            GCoord b;
            for (int i=0; i<3; i++)
            {
                GSReal b0 = faces[2*i+0]*Plane(p,1);
                GSReal b1 = faces[2*i+1]*Plane(p,1);
                b[i] = b0 / (b0 + b1);
            }
            return b;
        }

        //// Undeform a direction relative to a point in barycentric coordinates
        GCoord undeformDir(const GCoord& b, const GCoord& dir) const
        {
            /// we want to find b2 so that deform(b2)-deform(b) = dir
            /// we can use Newton's method using the jacobian of the deformation.
            Mat3x3d m = Jdeform(b);
            Mat3x3d minv;
            minv.invert(m);
            return minv*dir;
        }

        static GSReal interp(GSReal coef, GSReal a, GSReal b)
        {
            return a+coef*(b-a);
        }

    };

public:
    typedef Vec3Types InDataTypes;
    typedef Vec3Types DataTypes;
    typedef RegularGridTopology Topology;
    typedef FFDDistanceGridCollisionElement Element;

    Data< bool > usePoints;
    Data< bool > singleContact;

    //// virtual inherited from Base.
    void init();
    void draw(const VisualParams*,int index);
    void draw(const VisualParams* vparams);

    MechanicalState<DataTypes>* getDeformModel() { return ffd; }
    BaseMeshTopology* getDeformGrid() { return ffdMesh; }

    /// alias used by ContactMapper
    MechanicalState<DataTypes>* getMechanicalState() { return ffd; }
    BaseMeshTopology* getMeshTopology() { return ffdMesh; }

    DistanceGrid* getGrid(int index=0)
    {
        return elems[index].grid;
    }

    DeformedCube& getDeformCube(int index=0)
    {
        return elems[index];
    }

    void setGrid(DistanceGrid* surf, int index=0);

    /// -- CollisionModel interface
    void resize(int size);

    //// Create or update the bounding volume hierarchy.
    void computeBoundingTree(int maxDepth=0);
    bool canCollideWithElement(int index, CollisionModel* model2, int index2);

protected:
    helper::vector<DeformedCube> elems;

    /// Input data parameters
    Data< double >              scale;
    Data< double >              sampling;
    Data< helper::fixed_array<DistanceGrid::Coord,2> > box;
    Data< int >                 nx;
    Data< int >                 ny;
    Data< int >                 nz;
    DataFileName                dumpfilename;
    DataFileName                fileFFDDistanceGrid;

    MechanicalState<Vec3Types>* ffd;
    BaseMeshTopology*           ffdMesh;
    RegularGridTopology*        ffdRGrid;
    SparseGridTopology*         ffdSGrid;

    void updateGrid();

    FFDDistanceGridCollisionModel();
    ~FFDDistanceGridCollisionModel();
};

inline FFDDistanceGridCollisionElement::FFDDistanceGridCollisionElement(FFDDistanceGridCollisionModel* model, int index)
    : TCollisionElementIterator<FFDDistanceGridCollisionModel>(model, index)
{}

inline FFDDistanceGridCollisionElement::FFDDistanceGridCollisionElement(const core::CollisionElementIterator& i)
    : TCollisionElementIterator<FFDDistanceGridCollisionModel>(static_cast<FFDDistanceGridCollisionModel*>(i.getCollisionModel()), i.getIndex())
{
}

inline DistanceGrid* FFDDistanceGridCollisionElement::getGrid() { return model->getGrid(index); }
inline void FFDDistanceGridCollisionElement::setGrid(DistanceGrid* surf) { return model->setGrid(surf, index); }

} /// _distancegridcollisionmodel_

using _distancegridcollisionmodel_::RigidDistanceGridCollisionModel ;
using _distancegridcollisionmodel_::RigidDistanceGridCollisionElement;
using _distancegridcollisionmodel_::FFDDistanceGridCollisionModel ;
using _distancegridcollisionmodel_::FFDDistanceGridCollisionElement;
using sofa::component::container::DistanceGrid ;

//// Mapper for FFDDistanceGridCollisionModel
template <class DataTypes>
class ContactMapper<FFDDistanceGridCollisionModel,DataTypes> : public BarycentricContactMapper<FFDDistanceGridCollisionModel,DataTypes>
{
public:
    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Coord;
    int addPoint(const Coord& P, int index, Real&)
    {
        defaulttype::Vector3 bary;
        int elem = this->model->getDeformCube(index).elem;
        bary = this->model->getDeformCube(index).baryCoords(P);
        return this->mapper->addPointInCube(elem,bary.ptr());
    }
};


//// Mapper for RigidDistanceGridCollisionModel
template <class DataTypes>
class ContactMapper<RigidDistanceGridCollisionModel,DataTypes> : public RigidContactMapper<RigidDistanceGridCollisionModel,DataTypes>
{
public:
    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef RigidContactMapper<RigidDistanceGridCollisionModel,DataTypes> Inherit;
    typedef typename Inherit::MMechanicalState MMechanicalState;
    typedef typename Inherit::MCollisionModel MCollisionModel;

    MMechanicalState* createMapping(const char* name="contactPoints")
    {
        MMechanicalState* outmodel = Inherit::createMapping(name);
        return outmodel;
    }

    int addPoint(const Coord& P, int index, Real& r)
    {
        Coord trans = this->model->getInitRotation() * this->model->getInitTranslation();
        int i = Inherit::addPoint(P+trans, index, r);
        if (!this->mapping)
        {
            MCollisionModel* model = this->model;
            MMechanicalState* outmodel = this->outmodel.get();
            {
                helper::WriteAccessor<Data<VecCoord> > xData = *outmodel->write(core::VecCoordId::position());
                Coord& x = xData.wref()[i];

                if (model->isTransformed(index))
                    x = model->getTranslation(index) + model->getRotation(index) * P;
                else
                    x = P;
            }
            helper::ReadAccessor<Data<VecCoord> >  xData = *outmodel->read(core::ConstVecCoordId::position());
            helper::WriteAccessor<Data<VecDeriv> > vData = *outmodel->write(core::VecDerivId::velocity());
            const Coord& x = xData.ref()[i];
            Deriv& v       = vData.wref()[i];
            v.clear();

            /// estimating velocity
            double gdt = model->getPrevDt(index);
            if (gdt > 0.000001)
            {
                if (model->isTransformed(index))
                {
                    v = (x - (model->getPrevTranslation(index) + model->    getPrevRotation(index) * P)) * (1.0/gdt);
                }
                DistanceGrid* prevGrid = model->getPrevGrid(index);
                //DistanceGrid* grid = model->getGrid(index);
                //if (prevGrid != NULL && prevGrid != grid && prevGrid->inGrid(P))
                {
                    DistanceGrid::Coord coefs;
                    int i = prevGrid->index(P, coefs);
                    SReal d = prevGrid->interp(i,coefs);
                    if (sofa::helper::rabs(d) < 0.3) /// todo : control threshold
                    {
                        DistanceGrid::Coord n = prevGrid->grad(i,coefs);
                        v += n * (d  / ( n.norm() * gdt));
                        //std::cout << "Estimated v at "<<P<<" = "<<v<<" using distance from previous model "<<d<<std::endl;
                    }
                }
            }
        }
        return i;
    }
};


#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_COLLISION_DISTANCEGRIDCOLLISIONMODEL_CPP)
extern template class SOFA_VOLUMETRIC_DATA_API ContactMapper<FFDDistanceGridCollisionModel, defaulttype::Vec3Types>;
extern template class SOFA_VOLUMETRIC_DATA_API ContactMapper<RigidDistanceGridCollisionModel, defaulttype::Vec3Types>;

#  ifdef _MSC_VER
/// Manual declaration of non-specialized members, to avoid warnings from MSVC.
extern template SOFA_VOLUMETRIC_DATA_API void BarycentricContactMapper<FFDDistanceGridCollisionModel, defaulttype::Vec3Types>::cleanup();
extern template SOFA_VOLUMETRIC_DATA_API core::behavior::MechanicalState<defaulttype::Vec3Types>* BarycentricContactMapper<FFDDistanceGridCollisionModel, defaulttype::Vec3Types>::createMapping(const char*);
extern template SOFA_VOLUMETRIC_DATA_API void RigidContactMapper<RigidDistanceGridCollisionModel, defaulttype::Vec3Types>::cleanup();
extern template SOFA_VOLUMETRIC_DATA_API core::behavior::MechanicalState<defaulttype::Vec3Types>* RigidContactMapper<RigidDistanceGridCollisionModel, defaulttype::Vec3Types>::createMapping(const char*);
#  endif
#endif


} /// namespace collision

} /// namespace component

} /// namespace sofa

#endif
