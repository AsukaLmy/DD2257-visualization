/*********************************************************************
 *  Author  : Anke Friederici
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 **********************************************************************/

#include <inviwo/core/datastructures/geometry/basicmesh.h>
#include <inviwo/core/datastructures/volume/volumeram.h>
#include <labstreamlines/integrator.h>
#include <labutils/scalarvectorfield.h>
#include <labtopo/topology.h>
#include <labtopo/utils/gradients.h>

namespace inviwo
{

const vec4 Topology::ColorsCP[6] = {
    vec4(1, 0.95, 0.5, 1),  // Saddle - Yellow
    vec4(0.5, 0.5, 0.9, 1),  // AttractingNode - Sink - Blue
    vec4(0.9, 0.5, 0.5, 1),  // RepellingNode - Source - Red
    vec4(0.5, 0, 0.9, 1),// AttractingFocus - Purple
    vec4(0.9, 0.5, 0.0, 1),// RepellingFocus - Orange
    vec4(0.3, 0.6, 0.3, 1)   // Center - Green
};

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo Topology::processorInfo_{
    "org.inviwo.Topology",    // Class identifier
    "Vector Field Topology",// Display name
    "KTH Lab",                // Category
    CodeState::Experimental,// Code state
    Tags::None,               // Tags
};

const ProcessorInfo Topology::getProcessorInfo() const
{
    return processorInfo_;
}

Topology::Topology()
    : Processor()
    , inData("inData")
    , outMesh("meshOut")
    , meshBBoxOut("meshBBoxOut")
    , Maxdepth("Maxdepth", "Max depth", 9,1,10)
    , lengththreshold("lengththreshold","box length threshold",0,0,30)
    , modefythreshold("modefythreshold","threshold for point cluster", 0.05, 0, 0.5)
    , eigenstep("eigenstep","eigendirection step",0.01,0,0.5)
// TODO: Initialize additional properties
// propertyName("propertyIdentifier", "Display Name of the Propery",
// default value (optional), minimum value (optional), maximum value (optional), increment
// (optional)); propertyIdentifier cannot have spaces


{
    // Register Ports
    addPort(outMesh);
    addPort(inData);
    addPort(meshBBoxOut);

    // TODO: Register additional properties
    // addProperty(propertyName);
    addProperty(Maxdepth);
    addProperty(lengththreshold);
    addProperty(modefythreshold);
    addProperty(eigenstep);

}





void Topology::process()
{
    // Get input
    if (!inData.hasData())
    {
        return;
    }
    auto vol = inData.getData();

    // Retrieve data in a form that we can access it
    const VectorField2 vectorField = VectorField2::createFieldFromVolume(vol);

    // Add a bounding box to the mesh
    const dvec2& BBoxMin = vectorField.getBBoxMin();
    const dvec2& BBoxMax = vectorField.getBBoxMax();
    auto bboxMesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> bboxVertices;
    auto indexBufferBBox = bboxMesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    // Bounding Box vertex 0
    vec4 black = vec4(0, 0, 0, 1);
    Integrator::drawNextPointInPolyline(BBoxMin, black, indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(vec2(BBoxMin[0], BBoxMax[1]), black, indexBufferBBox.get(),
                                        bboxVertices);
    Integrator::drawNextPointInPolyline(BBoxMax, black, indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(vec2(BBoxMax[0], BBoxMin[1]), black, indexBufferBBox.get(),
                                        bboxVertices);
    // Connect back to the first point, to make a full rectangle
    indexBufferBBox->add(static_cast<std::uint32_t>(0));
    bboxMesh->addVertices(bboxVertices);
    meshBBoxOut.setData(bboxMesh);

    // Initialize mesh, vertices and index buffers for separatrices
    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;
    // Either add all line segments to this index buffer (one large buffer, two consecutive points
    // make up one line), or use several index buffers with connectivity type strip.
    auto indexBufferSeparatrices = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
    // auto indexBufferSeparatrices = mesh->addIndexBuffer(DrawType::Lines,
    // ConnectivityType::Strip);

    auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);

    // TODO: Compute the topological skeleton of the input vector field.
    // Find the critical points and color them according to their type.
    // Integrate all separatrices.

    size2_t dims = vectorField.getNumVerticesPerDim();
    std::vector<dvec2> criticalpointlist;
    std::vector<dvec2> boundarypointlist;
    int count = 0;
    // Looping through all values in the vector field.

    // extract the critical point
    for (size_t j = 0; j < dims[1]; ++j)
    {
        for (size_t i = 0; i < dims[0]; ++i)
        {
            //dvec2 vectorValue = vectorField.getValueAtVertex(size2_t(i, j));
            dvec2 pos_lb = vectorField.getPositionAtVertex(size2_t(i, j));
            dvec2 pos_rt = vectorField.getPositionAtVertex(size2_t(i+1, j+1));
            double boxlength_x = pos_rt[0] - pos_lb[0];
            double boxlength_y = pos_rt[1] - pos_lb[1];

            extractCriticalPoints(vectorField, criticalpointlist, 1, Maxdepth.get(), pos_lb[0], pos_lb[1], boxlength_x, boxlength_y, lengththreshold.get());
            modefyList(criticalpointlist,modefythreshold.get());
        }
    }

    //classification of points and draw streamline from saddle
    for (const auto& point : criticalpointlist)
    {
        dmat2 jacobian = vectorField.derive(point);
        auto eigenResult = util::eigenAnalysis(jacobian);
        double eigenindex = eigenstep.get();
        vec4 colorCenter = vec4();
        std::vector<dvec2> listForward;
        std::vector<dvec2> listBackward;
        // I1=I2=0
        if (eigenResult.eigenvaluesIm[0] == 0)
        {
            if (eigenResult.eigenvaluesRe[0] * eigenResult.eigenvaluesRe[1] < 0)
            {
                colorCenter = ColorsCP[static_cast<int>(TypeCP::Saddle)];

                // streamline
                dvec2 startpoint = dvec2((point[0] + eigenindex * eigenResult.eigenvectors[0][0]), (point[1] + eigenindex * eigenResult.eigenvectors[0][1]));
                listForward = Integrator::StreamLines(vectorField, startpoint, 0.01, 0, true, -1, 1000);
                listBackward = Integrator::StreamLines(vectorField, startpoint, 0.01, 1, true, -1, 1000);
                vec4 streamcolor = vec4(255, 255, 255, 1);
                dvec2 lastpoint = startpoint;
                for (const auto& streampoint : listForward)
                {
                    drawLineSegment(lastpoint, streampoint, streamcolor, indexBufferSeparatrices.get(), vertices);
                    lastpoint = streampoint;
                    //Integrator::drawPoint(streampoint, color, indexBufferPoints.get(), vertices);
                }
                lastpoint = startpoint;
                for (const auto& streampoint : listBackward)
                {
                    drawLineSegment(lastpoint, streampoint, streamcolor, indexBufferSeparatrices.get(), vertices);
                    lastpoint = streampoint;
                    //Integrator::drawPoint(streampoint, color, indexBufferPoints.get(), vertices);
                }

                //another eigenvector
                startpoint = dvec2((point[0] + eigenindex * eigenResult.eigenvectors[1][0]), (point[1] + eigenindex * eigenResult.eigenvectors[1][1]));
                listForward = Integrator::StreamLines(vectorField, startpoint, 0.01, 0, true, -1, 1000);
                listBackward = Integrator::StreamLines(vectorField, startpoint, 0.01, 1, true, -1, 1000);
                lastpoint = startpoint;
                for (const auto& streampoint : listForward)
                {
                    drawLineSegment(lastpoint, streampoint, streamcolor, indexBufferSeparatrices.get(), vertices);
                    lastpoint = streampoint;
                    //Integrator::drawPoint(streampoint, color, indexBufferPoints.get(), vertices);
                }
                lastpoint = startpoint;
                for (const auto& streampoint : listBackward)
                {
                    drawLineSegment(lastpoint, streampoint, streamcolor, indexBufferSeparatrices.get(), vertices);
                    lastpoint = streampoint;
                    //Integrator::drawPoint(streampoint, color, indexBufferPoints.get(), vertices);
                }
            }
            else if (eigenResult.eigenvaluesRe[0] > 0 && eigenResult.eigenvaluesRe[1] > 0)
            {
                colorCenter = ColorsCP[static_cast<int>(TypeCP::RepellingNode)];
            }
            else if (eigenResult.eigenvaluesRe[0] < 0 && eigenResult.eigenvaluesRe[1] < 0)
            {
                colorCenter = ColorsCP[static_cast<int>(TypeCP::AttractingNode)];
            }
        }
        //I1=-I2!=0
        else
        {
            if (eigenResult.eigenvaluesRe[0] == eigenResult.eigenvaluesRe[1] && eigenResult.eigenvaluesRe[0] == 0)
            {
                colorCenter = ColorsCP[static_cast<int>(TypeCP::Center)];
            }
            else if (eigenResult.eigenvaluesRe[0] == eigenResult.eigenvaluesRe[1] && eigenResult.eigenvaluesRe[0] > 0)
            {
                colorCenter = ColorsCP[static_cast<int>(TypeCP::RepellingFocus)];
            }
            else if (eigenResult.eigenvaluesRe[0] == eigenResult.eigenvaluesRe[1] && eigenResult.eigenvaluesRe[0] < 0)
            {
                colorCenter = ColorsCP[static_cast<int>(TypeCP::AttractingFocus)];
            }
        }
        
        Integrator::drawPoint(point, colorCenter, indexBufferPoints.get(), vertices);
        count++;
    }
    LogProcessorInfo(count);
    // Other helpful functions
    // dvec2 pos = vectorField.getPositionAtVertex(size2_t(i, j));
    // Computing the jacobian at a position
    // dmat2 jacobian = vectorField.derive(pos);
    // Doing the eigen analysis
    // auto eigenResult = util::eigenAnalysis(jacobian);
    // The result of the eigen analysis has attributed eigenvaluesRe eigenvaluesIm and
    // eigenvectors

    //boundary switch points
    boundaryswitchpoint(vectorField, boundarypointlist, BBoxMin[0], BBoxMin[1], BBoxMax[0], BBoxMax[1],0.01);
    modefyList(boundarypointlist, modefythreshold.get());
    for (const auto& point : boundarypointlist)
    {
        dmat2 jacobian = vectorField.derive(point);
        auto eigenResult = util::eigenAnalysis(jacobian);
        double eigenindex = eigenstep.get();
        
        std::vector<dvec2> listForward;
        std::vector<dvec2> listBackward;
        // streamline
        dvec2 startpoint = dvec2((point[0] + eigenindex * eigenResult.eigenvectors[0][0]), (point[1] + eigenindex * eigenResult.eigenvectors[0][1]));
        listForward = Integrator::StreamLines(vectorField, startpoint, 0.01, 0, true, -1, 1000);
        listBackward = Integrator::StreamLines(vectorField, startpoint, 0.01, 1, true, -1, 1000);
        vec4 streamcolor = vec4(128, 128, 128, 1);
        dvec2 lastpoint = startpoint;
        for (const auto& streampoint : listForward)
        {
            drawLineSegment(lastpoint, streampoint, streamcolor, indexBufferSeparatrices.get(), vertices);
            lastpoint = streampoint;
            //Integrator::drawPoint(streampoint, color, indexBufferPoints.get(), vertices);
        }
        lastpoint = startpoint;
        for (const auto& streampoint : listBackward)
        {
            drawLineSegment(lastpoint, streampoint, streamcolor, indexBufferSeparatrices.get(), vertices);
            lastpoint = streampoint;
            //Integrator::drawPoint(streampoint, color, indexBufferPoints.get(), vertices);
        }

        //another eigenvector
        startpoint = dvec2((point[0] + eigenindex * eigenResult.eigenvectors[1][0]), (point[1] + eigenindex * eigenResult.eigenvectors[1][1]));
        listForward = Integrator::StreamLines(vectorField, startpoint, 0.01, 0, true, -1, 1000);
        listBackward = Integrator::StreamLines(vectorField, startpoint, 0.01, 1, true, -1, 1000);
        lastpoint = startpoint;
        for (const auto& streampoint : listForward)
        {
            drawLineSegment(lastpoint, streampoint, streamcolor, indexBufferSeparatrices.get(), vertices);
            lastpoint = streampoint;
            //Integrator::drawPoint(streampoint, color, indexBufferPoints.get(), vertices);
        }
        lastpoint = startpoint;
        for (const auto& streampoint : listBackward)
        {
            drawLineSegment(lastpoint, streampoint, streamcolor, indexBufferSeparatrices.get(), vertices);
            lastpoint = streampoint;
            //Integrator::drawPoint(streampoint, color, indexBufferPoints.get(), vertices);
        }
        

        dvec4 colorboundary = dvec4(0, 0, 0, 1);
        Integrator::drawPoint(point, colorboundary, indexBufferPoints.get(), vertices);
    }
    //boundary switch points end




    // Accessing the colors
    vec4 colorCenter = ColorsCP[static_cast<int>(TypeCP::Center)];

    mesh->addVertices(vertices);
    outMesh.setData(mesh);
}


void Topology::boundaryswitchpoint(const VectorField2& vectorField, std::vector<dvec2>& boundarypointlist, double lb_x, double lb_y, double rt_x, double rt_y, double eps)
{
    double step_x = (rt_x - lb_x) / 100;
    double step_y = (rt_y - lb_y) / 100;
    for (double x = lb_x; x <= rt_x; x = x + step_x)
    {
        dvec2 vl = vectorField.interpolate(dvec2(x - eps, lb_y));
        dvec2 vr = vectorField.interpolate(dvec2(x + eps, lb_y));
        if (vl[1] * vr[1] < 0)
        {
            boundarypointlist.push_back(dvec2(x, lb_y));
        }
    }
    for (double x = lb_x; x <= rt_x; x = x + step_x)
    {
        dvec2 vl = vectorField.interpolate(dvec2(x - eps, rt_y));
        dvec2 vr = vectorField.interpolate(dvec2(x + eps, rt_y));
        if (vl[1] * vr[1] < 0)
        {
            boundarypointlist.push_back(dvec2(x, rt_y));
        }
    }

    for (double y = lb_y; y <= rt_y;  y = y + step_y)
    {
        dvec2 vl = vectorField.interpolate(dvec2(lb_x, y-eps));
        dvec2 vr = vectorField.interpolate(dvec2(lb_x, y+eps));
        if (vl[0] * vr[0] < 0)
        {
            boundarypointlist.push_back(dvec2(lb_x, y));
        }
    }
    for (double y = lb_y; y <= rt_y; y = y + step_y)
    {
        dvec2 vl = vectorField.interpolate(dvec2(rt_x, y - eps));
        dvec2 vr = vectorField.interpolate(dvec2(rt_x, y + eps));
        if (vl[0] * vr[0] < 0)
        {
            boundarypointlist.push_back(dvec2(rt_x, y));
        }
    }




}


bool Topology::zerodetection(const VectorField2& vectorField, double x, double y, double boxlength_x, double boxlength_y) {
    dvec2 v00 = vectorField.interpolate(dvec2(x, y));                            //left bottom
    dvec2 v01 = vectorField.interpolate(dvec2(x, y + boxlength_y));              //left top
    dvec2 v10 = vectorField.interpolate(dvec2(x + boxlength_x, y));              //right bottom
    dvec2 v11 = vectorField.interpolate(dvec2(x + boxlength_x, y + boxlength_y));//right top

    //check if any parameter changes
    bool changex = (v00[0] * v01[0] <= 0 || v00[0] * v10[0] <= 0 || v00[0] * v11[0] <= 0);
    bool changey = (v00[1] * v01[1] <= 0 || v00[1] * v10[1] <= 0 || v00[1] * v11[1] <= 0);

    if (changex && changey)
    {
        return true;
    }
    else
    {
        return false;
    }

}


void Topology::modefyList(std::vector<dvec2>& criticalpointlist, double threshold) {
    for (size_t i = 0; i < criticalpointlist.size(); ++i)
    {
        for (size_t j = i + 1; j < criticalpointlist.size(); ++j)
        {           
            double distance = glm::distance(criticalpointlist[i], criticalpointlist[j]);

            if (distance < threshold)
            {
                dvec2 mergedPoint = (criticalpointlist[i] + criticalpointlist[j]) / 2.0;
                criticalpointlist[i] = mergedPoint;
                criticalpointlist.erase(criticalpointlist.begin() + j);
                --j;
            }
        }
    }
}

void Topology::extractCriticalPoints(const VectorField2& vectorField, std::vector<dvec2>& criticalpointlist, int depth, int maxdepth, double x, double y, double boxlength_x, double boxlength_y,
                           double lengththreshold)
{

    // check if reaching maxdepth
    if (depth > maxdepth || boxlength_x < lengththreshold || boxlength_y < lengththreshold)
    {
        dvec2 criticalpoint = dvec2(x, y);
        criticalpointlist.push_back(criticalpoint);
        return;
        //return criticalpoint;
    }
    // find four point around box
    dvec2 v00 = vectorField.interpolate(dvec2(x, y));                            //left bottom
    dvec2 v01 = vectorField.interpolate(dvec2(x, y + boxlength_y));              //left top
    dvec2 v10 = vectorField.interpolate(dvec2(x + boxlength_x, y));              //right bottom
    dvec2 v11 = vectorField.interpolate(dvec2(x + boxlength_x, y + boxlength_y));//right top

    

    bool changex = (v00[0] * v01[0] < 0 || v00[0] * v10[0] < 0 || v00[0] * v11[0] < 0);
    bool changey = (v00[1] * v01[1] < 0 || v00[1] * v10[1] < 0 || v00[1] * v11[1] < 0);

    if (changex && changey)
    {
        extractCriticalPoints(vectorField, criticalpointlist, depth+1, maxdepth, x, y, boxlength_x / 2, boxlength_y / 2, lengththreshold);
        extractCriticalPoints(vectorField, criticalpointlist, depth+1, maxdepth, x + boxlength_x / 2, y, boxlength_x / 2, boxlength_y / 2, lengththreshold);
        extractCriticalPoints(vectorField, criticalpointlist, depth+1, maxdepth, x, y + boxlength_y / 2, boxlength_x / 2, boxlength_y / 2, lengththreshold);
        extractCriticalPoints(vectorField, criticalpointlist, depth+1, maxdepth, x + boxlength_x / 2, y + boxlength_y / 2, boxlength_x / 2, boxlength_y / 2, lengththreshold);
    }

}




void Topology::drawLineSegment(const dvec2& v1, const dvec2& v2, const vec4& color,
                               IndexBufferRAM* indexBuffer,
                               std::vector<BasicMesh::Vertex>& vertices)
{
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v1[0], v1[1], 0), vec3(0, 0, 1), vec3(v1[0], v1[1], 0), color});
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v2[0], v2[1], 0), vec3(0, 0, 1), vec3(v2[0], v2[1], 0), color});
}


}// namespace inviwo
