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
#ifndef SCRIPTFUNCTION_H
#define SCRIPTFUNCTION_H

#include <SofaPython27/config.h>
#include <string>

namespace sofa
{

namespace core
{

namespace objectmodel
{

class SOFA_SOFAPYTHON27_API ScriptFunctionParameter
{
protected:
	ScriptFunctionParameter();

public:
	virtual ~ScriptFunctionParameter();

};

class SOFA_SOFAPYTHON27_API ScriptFunctionResult
{
protected:
	ScriptFunctionResult();

public:
	virtual ~ScriptFunctionResult();

};

class SOFA_SOFAPYTHON27_API ScriptFunction
{
protected:
    ScriptFunction();

public:
	virtual ~ScriptFunction();

	void operator()(const ScriptFunctionParameter*, ScriptFunctionResult*) const;

protected:
	virtual void onCall(const ScriptFunctionParameter*, ScriptFunctionResult*) const = 0;

};

} // namespace objectmodel

} // namespace core

} // namespace sofa

#endif // SCRIPTFUNCTION_H
