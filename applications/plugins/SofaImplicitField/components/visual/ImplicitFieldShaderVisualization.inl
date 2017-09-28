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
#include "ImplicitFieldShaderVisualization.h"

namespace sofa
{

namespace component
{

namespace visual
{

namespace _pointcloudimplicitfieldvisualization_
{

ImplicitFieldShaderVisualization::ImplicitFieldShaderVisualization() :
    l_field(initLink("field", "The field to render.")),
    d_vertFilename(initData(&d_vertFilename, std::string("shaders/implicitShape.vert"), "fileVertexShaders", "Set the vertex shader filename to load")),
    d_fragFilename(initData(&d_fragFilename, std::string("shaders/implicitShape.frag"), "fileFragmentShaders", "Set the fragment shader filename to load"))
{
    f_listening.setValue(true);
    wheelDelta = 0.5;
}

ImplicitFieldShaderVisualization::~ImplicitFieldShaderVisualization()
{
    shader->TurnOff();
    shader->Release();
    delete shader;
}

void ImplicitFieldShaderVisualization::init()
{
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    float pixelSize[2];
    pixelSize[0] = viewport[2];
    pixelSize[1] = viewport[3];
    mouseX = pixelSize[0] /2.0;
    mouseY = pixelSize[1] /2.0;

    shader = new sofa::helper::gl::GLSLShader();

    shader->SetVertexShaderFileName(d_vertFilename.getFullPath());
    shader->SetFragmentShaderFileName(d_fragFilename.getFullPath());
    m_datatracker.trackData(*l_field.get()->findData("state"));
    shaderGenerationCodeHasChanged();

}

void ImplicitFieldShaderVisualization::reinit()
{
}

void ImplicitFieldShaderVisualization::initVisual()
{
    if (!sofa::helper::gl::GLSLShader::InitGLSL())
    {
        serr << "InitGLSL failed" << sendl;
        return;
    }
    shader->InitShaders();
}

void ImplicitFieldShaderVisualization::stop()
{
    if(shader->IsReady())
    {
        glDisable(GL_VERTEX_PROGRAM_TWO_SIDE);
        glClampColorARB(GL_CLAMP_VERTEX_COLOR, GL_TRUE);
        shader->TurnOff();
    }
}

void ImplicitFieldShaderVisualization::shaderGenerationCodeHasChanged()
{
    std::ofstream myfile;
    myfile.open (d_fragFilename.getFullPath());
    myfile << generateFragmentShader();
    myfile.close();

    myfile;
    myfile.open (d_vertFilename.getFullPath());
    myfile << generateVertexShader();
    myfile.close();
}

void ImplicitFieldShaderVisualization::start()
{
    sofa::helper::system::TemporaryLocale locale(LC_CTYPE, "");

    if(m_datatracker.isDirty())
    {
        std::ofstream myfile;
        myfile.open (d_fragFilename.getFullPath());
        myfile << generateFragmentShader();
        myfile.close();
        m_datatracker.clean();
    }

    if(shader->IsReady())
    {
        shader->TurnOn();

        GLint viewport[4];
        glGetIntegerv(GL_VIEWPORT, viewport);
        float pixelSize[2];
        pixelSize[0] = viewport[2];
        pixelSize[1] = viewport[3];

        shader->SetFloat2(shader->GetVariable("resolution"), pixelSize[0], pixelSize[1]);
        shader->SetFloat2(shader->GetVariable("mouse"), mouseX, mouseY);
        shader->SetFloat(shader->GetVariable("wheelDelta"), wheelDelta);


        //        Simulation* simu = sofa::simulation::getSimulation();
        //        simu->findData("");

        //        /// TODO se mettre d'accord sur le contenu des maps :/


        std::map<std::string, std::vector<GLSLCodeFragment>> glslMap = l_field->getGLSLCode();
        std::map<std::string, std::vector<GLSLCodeFragment>>::iterator itFind = glslMap.find("variable");
        if(itFind != glslMap.end())
        {
            std::vector<GLSLCodeFragment> uniforms = itFind->second;
            for( std::vector<GLSLCodeFragment>::iterator it = uniforms.begin(); it != uniforms.end(); it++)
            {
                GLSLCodeFragment tmpGLSLCode = *it;
                std::string uniformType = tmpGLSLCode.m_type;
                std::string uniformName = tmpGLSLCode.m_name;
                std::string uniformValue = tmpGLSLCode.m_value;

                std::regex rgx("\\s+");
                std::sregex_token_iterator iter(uniformValue.begin(), uniformValue.end(), rgx, -1);
                std::sregex_token_iterator end;

                if(uniformType.compare("float") == 0)
                    shader->SetFloat(shader->GetVariable(uniformName), std::atof(std::string((*iter)).c_str()));
                else if(uniformType.compare("vec2") == 0)
                {
                    float data[2];
                    data[0] = std::atof(std::string((*iter)).c_str());
                    data[1] = std::atof(std::string((*++iter)).c_str());
                    shader->SetFloat2(shader->GetVariable(uniformName), data[0], data[1]);
                }
                else if(uniformType.compare("vec3") == 0)
                {
                    float data[3];
                    data[0] = std::atof(std::string((*iter)).c_str());
                    data[1] = std::atof(std::string((*++iter)).c_str());
                    data[2] = std::atof(std::string((*++iter)).c_str());
                    shader->SetFloat3(shader->GetVariable(uniformName), data[0], data[1], data[2]);
                }
            }
        }

        glClampColorARB(GL_CLAMP_VERTEX_COLOR, GL_TRUE);
        glEnable(GL_VERTEX_PROGRAM_TWO_SIDE);
    }
}

bool ImplicitFieldShaderVisualization::isActive()
{
    return true;
}

void ImplicitFieldShaderVisualization::handleEvent(core::objectmodel::Event * event)
{
    if(MouseEvent* ev = dynamic_cast< MouseEvent *>(event))
    {
        if (ev->getState() == MouseEvent::Move)
        {
            mouseX = ev->getPosX();
            mouseY = ev->getPosY();
        }
        if (ev->getState() == MouseEvent::Wheel)
            wheelDelta += ev->getWheelDelta()/120 * 0.1;
    }
}

bool ImplicitFieldShaderVisualization::drawScene(VisualParams* vp)
{
    return true;
}

std::string ImplicitFieldShaderVisualization::generateVertexShader()
{
    std::string vertexShaderText;
    vertexShaderText = std::string() +
            "attribute vec2 g;\n"
            "void main() \n"
            "{\n"
            "    gl_Position = vec4(g.xy, 1.0, 1.0);\n"
            "}\n";
    return vertexShaderText;
}

std::string ImplicitFieldShaderVisualization::implicitFunction()
{
    std::string tmp;
    std::string implicitFunction;
    implicitFunction.append(
                "    res = minVec4(\n"
                "        vec4(sdPlane(pos), vec3(0.45, 0.45, 0.45)),\n"
                "        vec4(sdSphere(pos-vec3( -.0, 0.75, 0.0), .5), vec3(1.0, 1.0, 0.0))\n"
                "   );    \n"
                ); /// Default value if eval is empty

    std::map<std::string, std::vector<GLSLCodeFragment>> glslMap = l_field->getGLSLCode();
    std::map<std::string, std::vector<GLSLCodeFragment>>::iterator itFind = glslMap.find("eval");

    if(itFind != glslMap.end())
    {
        implicitFunction.clear();
        std::vector<GLSLCodeFragment> evals = itFind->second;
        for( std::vector<GLSLCodeFragment>::iterator it = evals.begin(); it != evals.end(); it++)
        {
            std::vector<GLSLCodeFragment> uniforms = itFind->second;
            GLSLCodeFragment data = *it;
            implicitFunction.append(
                        "    x = pos.x - evalPosition" + data.m_dataname + ".x;\n"
                                                                           "    y = pos.y - evalPosition" + data.m_dataname + ".y;\n"
                                                                                                                              "    z = pos.z - evalPosition" + data.m_dataname + ".z;\n"
                        );
            implicitFunction.append(
                        "    res = minVec4(\n"
                        "        res,\n"
                        );

            implicitFunction.append("\t\tvec4(" + data.m_value + ", evalColor" + data.m_dataname + ")\n");
            implicitFunction.append("   );    \n");
        }
    }

    tmp.append(
                "float sdPlane( vec3 p )\n"
                "{\n"
                "   return p.y;\n"
                "}\n"
                "\n"
                "float sdSphere( vec3 p, float s )\n"
                "{\n"
                "   return length(p)-s;\n"
                "}\n"
                "\n"
                "vec4 minVec4( vec4 d1, vec4 d2 )\n"
                "{\n"
                "    return (d1.x<d2.x) ? d1 : d2;\n"
                "}\n"
                "\n"
                );

    tmp.append(
                "vec4 map( in vec3 pos )\n"
                "{\n"
                "   float x = pos.x;\n"
                "   float y = pos.y;\n"
                "   float z = pos.z;\n"
                "   vec4 res = vec4(sdPlane(pos), vec3(0.45, 0.45, 0.45));\n"
                );

    tmp.append(implicitFunction);

    tmp.append(
                "   return res;\n"
                "}\n"
                );
    tmp.append("\n");

    return tmp;
}

std::string ImplicitFieldShaderVisualization::viewFunction()
{
    std::string tmp;
    tmp.append(
                "mat3 setCamera( in vec3 ro, in vec3 ta, float cr )\n"
                "{\n"
                "    vec3 cw = normalize(ta-ro);\n"
                "    vec3 cp = vec3(sin(cr), cos(cr),0.0);\n"
                "    vec3 cu = normalize( cross(cw,cp) );\n"
                "    vec3 cv = normalize( cross(cu,cw) );\n"
                "    return mat3( cu, cv, cw );\n"
                "}\n"
                );
    tmp.append("\n");
    return tmp;
}

std::string ImplicitFieldShaderVisualization::mainFragmenShaderFunction()
{
    std::string tmp;
    tmp.append(
                "void main()\n"
                "{\n"
                "    vec2 mouseR = mouse.xy/resolution.xy;\n"
                "    vec2 p = (-resolution.xy + 2.0*gl_FragCoord.xy)/resolution.y;\n"
                "    vec3 ro = vec3( 4.0 * cos(7.0*mouseR.x) * (wheelDelta), 4.0*mouseR.y, 4.0*sin(7.0*mouseR.x) * (wheelDelta) );\n"
                "    vec3 ta = vec3( 0.0, 0.0, 0.0 );\n"
                "    mat3 ca = setCamera( ro, ta, 0.0 );\n"
                "    vec3 rayDir = ca * normalize( vec3(p.xy, 2.0) );\n"
                "    gl_FragColor = vec4( render( ro, rayDir ), 1.0 );\n"
                "}\n"
                );
    tmp.append("\n");
    return tmp;
}

std::string ImplicitFieldShaderVisualization::renderFunction()
{
    std::string tmp;
    tmp.append(
                "float softshadow( in vec3 ro, in vec3 rd, in float mint, in float tmax )\n"
                "{\n"
                "   float res = 1.0;\n"
                "   float t = mint;\n"
                "   for( int i=0; i<16; i++ )\n"
                "   {\n"
                "       float h = map( ro + rd*t ).x;\n"
                "       res = min( res, 8.0*h/t );\n"
                "       t += clamp( h, 0.02, 0.10 );\n"
                "       if( h<0.001 || t>tmax ) break;\n"
                "   }\n"
                "   return clamp( res, 0.0, 1.0 );\n"
                "}\n"
                "\n"
                "vec3 estimateNormal(vec3 p) {\n"
                "    return normalize(vec3(\n"
                "        map(vec3(p.x + EPSILON, p.y, p.z)).x - map(vec3(p.x - EPSILON, p.y, p.z)).x,\n"
                "        map(vec3(p.x, p.y + EPSILON, p.z)).x - map(vec3(p.x, p.y - EPSILON, p.z)).x,\n"
                "        map(vec3(p.x, p.y, p.z  + EPSILON)).x - map(vec3(p.x, p.y, p.z - EPSILON)).x\n"
                "    ));\n"
                "}\n"
                "\n"
                "vec3 render( in vec3 ro, in vec3 rd )\n"
                "{\n"
                "    vec3 col = vec3(1.0, 1.0, 1.0) + rd.y;\n"
                "    vec4 res = castRay(ro,rd);\n"
                "    float t = res.x;\n"
                "    vec3 m = res.yzw;\n"

                "    vec3 pos = ro + t * rd;\n"
                "    vec3 nor = estimateNormal( pos );\n"
                "    vec3 ref = reflect( rd, nor );\n"

                "    col = m;\n"

                "    vec3  lig = normalize( vec3(-0.4, 0.7, -0.6) );\n"
                "    float amb = clamp( 0.5+0.5*nor.y, 0.0, 1.0 );\n"
                "    float dif = clamp( dot( nor, lig ), 0.0, 1.0 );\n"
                "    float bac = clamp( dot( nor, normalize(vec3(-lig.x,0.0,-lig.z))), 0.0, 1.0 )*clamp( 1.0-pos.y,0.0,1.0);\n"
                "    float dom = smoothstep( -0.1, 0.1, ref.y );\n"
                "    float fre = pow( clamp(1.0+dot(nor,rd),0.0,1.0), 2.0 );\n"
                "    float spe = pow(clamp( dot( ref, lig ), 0.0, 1.0 ),16.0);\n"
                "    dif *= softshadow( pos, lig, 0.02, 2.5 );\n"
                "    dom *= softshadow( pos, ref, 0.02, 2.5 );\n"
                "    vec3 lin = vec3(0.0);\n"
                "    lin += 1.30*dif*vec3(1.00,0.80,0.55);\n"
                "    lin += 2.00*spe*vec3(1.00,0.90,0.70)*dif;\n"
                "    lin += 0.40*amb*vec3(0.40,0.60,1.00);\n"
                "    lin += 0.50*dom*vec3(0.40,0.60,1.00);\n"
                "    col = col*lin;\n"
                "    col = mix( col, vec3(0.9,0.9,1.0), 1.0-exp( -0.0002*t*t*t ) );\n"
                "    return vec3( clamp(col,0.0,1.0) );\n"
                "}\n"
                );
    tmp.append("\n");
    return tmp;
}

std::string ImplicitFieldShaderVisualization::rayMarchingFunction()
{
    std::string tmp;
    tmp.append(
                "vec4 castRay( in vec3 ro, in vec3 rd )\n"
                "{\n"
                "    float depth = MIN_DIST;\n"
                "    vec3 m = vec3(-1.0, -1.0, -1.0);\n"
                "    for( int i = 0; i < MAX_MARCHING_STEPS; i++ )\n"
                "    {\n"
                "        vec4 res = map( ro+rd*depth );\n"
                "        if( res.x < EPSILON || depth > MAX_DIST ) break;\n"
                "        depth += res.x;\n"
                "            m = res.yzw;\n"
                "    }\n"
                "    return vec4( depth, m );\n"
                "}\n"
                );
    tmp.append("\n");
    return tmp;
}

std::string ImplicitFieldShaderVisualization::uniformsAndConst()
{
    std::string tmp;
    tmp.append(
                "uniform vec2 resolution;\n"
                "uniform vec2 mouse;\n"
                "uniform float wheelDelta;\n"

                "const int MAX_MARCHING_STEPS = 255;\n"
                "const float MIN_DIST = 0.0;\n"
                "const float MAX_DIST = 100.0;\n"
                "const float EPSILON = 0.00001;\n"
                );

    std::map<std::string, std::vector<GLSLCodeFragment>> glslMap = l_field->getGLSLCode();
    std::map<std::string, std::vector<GLSLCodeFragment>>::iterator itFind = glslMap.find("variable");
    if(itFind != glslMap.end())
    {
        std::vector<GLSLCodeFragment> uniforms = itFind->second;
        for( std::vector<GLSLCodeFragment>::iterator it = uniforms.begin(); it != uniforms.end(); it++)
        {
            GLSLCodeFragment tmpGLSLCode = *it;
            tmp.append("uniform " + tmpGLSLCode.m_type + " " + tmpGLSLCode.m_dataname + ";\n");
        }
    }
    tmp.append("\n");

    return tmp;
}

std::string ImplicitFieldShaderVisualization::generateFragmentShader()
{
    std::string fragmentShaderText;
    fragmentShaderText = std::string() + "";

    fragmentShaderText += uniformsAndConst();
    fragmentShaderText += implicitFunction();
    fragmentShaderText += rayMarchingFunction();
    fragmentShaderText += renderFunction();
    fragmentShaderText += viewFunction();
    fragmentShaderText += mainFragmenShaderFunction();

    return fragmentShaderText;
}


} /// namespace _pointcloudimplicitfieldvisualization_
} /// namespace visual
} /// namespace component
} /// namespace sofa
