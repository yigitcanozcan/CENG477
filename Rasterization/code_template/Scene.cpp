#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>

#include "tinyxml2.h"
#include "Triangle.h"
#include "Helpers.h"
#include "Scene.h"

using namespace tinyxml2;
using namespace std;


/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *xmlElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *rootNode = xmlDoc.FirstChild();

	// read background color
	xmlElement = rootNode->FirstChildElement("BackgroundColor");
	str = xmlElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	xmlElement = rootNode->FirstChildElement("Culling");
	if (xmlElement != NULL)
	{
		str = xmlElement->GetText();

		if (strcmp(str, "enabled") == 0)
		{
			this->cullingEnabled = true;
		}
		else
		{
			this->cullingEnabled = false;
		}
	}

	// read cameras
	xmlElement = rootNode->FirstChildElement("Cameras");
	XMLElement *camElement = xmlElement->FirstChildElement("Camera");
	XMLElement *camFieldElement;
	while (camElement != NULL)
	{
		Camera *camera = new Camera();

		camElement->QueryIntAttribute("id", &camera->cameraId);

		// read projection type
		str = camElement->Attribute("type");

		if (strcmp(str, "orthographic") == 0)
		{
			camera->projectionType = ORTOGRAPHIC_PROJECTION;
		}
		else
		{
			camera->projectionType = PERSPECTIVE_PROJECTION;
		}

		camFieldElement = camElement->FirstChildElement("Position");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->position.x, &camera->position.y, &camera->position.z);

		camFieldElement = camElement->FirstChildElement("Gaze");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->gaze.x, &camera->gaze.y, &camera->gaze.z);

		camFieldElement = camElement->FirstChildElement("Up");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->v.x, &camera->v.y, &camera->v.z);

		camera->gaze = normalizeVec3(camera->gaze);
		camera->u = crossProductVec3(camera->gaze, camera->v);
		camera->u = normalizeVec3(camera->u);

		camera->w = inverseVec3(camera->gaze);
		camera->v = crossProductVec3(camera->u, camera->gaze);
		camera->v = normalizeVec3(camera->v);

		camFieldElement = camElement->FirstChildElement("ImagePlane");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &camera->left, &camera->right, &camera->bottom, &camera->top,
			   &camera->near, &camera->far, &camera->horRes, &camera->verRes);

		camFieldElement = camElement->FirstChildElement("OutputName");
		str = camFieldElement->GetText();
		camera->outputFilename = string(str);

		this->cameras.push_back(camera);

		camElement = camElement->NextSiblingElement("Camera");
	}

	// read vertices
	xmlElement = rootNode->FirstChildElement("Vertices");
	XMLElement *vertexElement = xmlElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (vertexElement != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = vertexElement->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = vertexElement->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		this->vertices.push_back(vertex);
		this->colorsOfVertices.push_back(color);

		vertexElement = vertexElement->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	xmlElement = rootNode->FirstChildElement("Translations");
	XMLElement *translationElement = xmlElement->FirstChildElement("Translation");
	while (translationElement != NULL)
	{
		Translation *translation = new Translation();

		translationElement->QueryIntAttribute("id", &translation->translationId);

		str = translationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		this->translations.push_back(translation);

		translationElement = translationElement->NextSiblingElement("Translation");
	}

	// read scalings
	xmlElement = rootNode->FirstChildElement("Scalings");
	XMLElement *scalingElement = xmlElement->FirstChildElement("Scaling");
	while (scalingElement != NULL)
	{
		Scaling *scaling = new Scaling();

		scalingElement->QueryIntAttribute("id", &scaling->scalingId);
		str = scalingElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		this->scalings.push_back(scaling);

		scalingElement = scalingElement->NextSiblingElement("Scaling");
	}

	// read rotations
	xmlElement = rootNode->FirstChildElement("Rotations");
	XMLElement *rotationElement = xmlElement->FirstChildElement("Rotation");
	while (rotationElement != NULL)
	{
		Rotation *rotation = new Rotation();

		rotationElement->QueryIntAttribute("id", &rotation->rotationId);
		str = rotationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		this->rotations.push_back(rotation);

		rotationElement = rotationElement->NextSiblingElement("Rotation");
	}

	// read meshes
	xmlElement = rootNode->FirstChildElement("Meshes");

	XMLElement *meshElement = xmlElement->FirstChildElement("Mesh");
	while (meshElement != NULL)
	{
		Mesh *mesh = new Mesh();

		meshElement->QueryIntAttribute("id", &mesh->meshId);

		// read projection type
		str = meshElement->Attribute("type");

		if (strcmp(str, "wireframe") == 0)
		{
			mesh->type = WIREFRAME_MESH;
		}
		else
		{
			mesh->type = SOLID_MESH;
		}

		// read mesh transformations
		XMLElement *meshTransformationsElement = meshElement->FirstChildElement("Transformations");
		XMLElement *meshTransformationElement = meshTransformationsElement->FirstChildElement("Transformation");

		while (meshTransformationElement != NULL)
		{
			char transformationType;
			int transformationId;

			str = meshTransformationElement->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			mesh->transformationTypes.push_back(transformationType);
			mesh->transformationIds.push_back(transformationId);

			meshTransformationElement = meshTransformationElement->NextSiblingElement("Transformation");
		}

		mesh->numberOfTransformations = mesh->transformationIds.size();

		// read mesh faces
		char *row;
		char *cloneStr;
		int v1, v2, v3;
		XMLElement *meshFacesElement = meshElement->FirstChildElement("Faces");
		str = meshFacesElement->GetText();
		cloneStr = strdup(str);

		row = strtok(cloneStr, "\n");
		while (row != NULL)
		{
			int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

			if (result != EOF)
			{
				mesh->triangles.push_back(Triangle(v1, v2, v3));
			}
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		this->meshes.push_back(mesh);

		meshElement = meshElement->NextSiblingElement("Mesh");
	}
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
			}

			this->image.push_back(rowOfColors);
		}
	}
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				this->image[i][j].r = this->backgroundColor.r;
				this->image[i][j].g = this->backgroundColor.g;
				this->image[i][j].b = this->backgroundColor.b;
			}
		}
	}
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}


/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open(camera->outputFilename.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFilename << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
	os_type == 1 		-> Ubuntu
	os_type == 2 		-> Windows
	os_type == other	-> No conversion
*/
void Scene::convertPPMToPNG(string ppmFileName, int osType)
{
	string command;

	// call command on Ubuntu
	if (osType == 1)
	{
		command = "./magick " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// call command on Windows
	else if (osType == 2)
	{
		command = "magick " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// default action - don't do conversion
	else
	{
	}
}


Vec4 homogenize(Vec4 vector){
	double normalizer = vector.t;
	
	double v0 = vector.x / normalizer;
	double v1 = vector.y / normalizer;
	double v2 = vector.z / normalizer;
	double v3 = 1;
	
	Vec4 homogenized = Vec4(v0, v1, v2, v3, vector.colorId);

	return homogenized;
}

Matrix4 Identity(){
	Matrix4 result;
	for(int i = 0 ; i < 4 ; i++ ){
		for(int j = 0; j < 4 ; j++){
			if(i==j){
				result.values[i][j]=1.0 ;
			}else{
				result.values[i][j]=0.0 ;
			}
		}
	}
	return result;
}

Matrix4 translationMatrix(double tx, double ty, double tz){
	Matrix4 result = Identity();
	result.values[0][3] = tx;
	result.values[1][3] = ty;
	result.values[2][3] = tz;

	return result;
		
}

Matrix4 scalingMatrix(double sx, double sy, double sz){
	Matrix4 result = Identity();
	result.values[0][0] = sx;
	result.values[1][1] = sy;
	result.values[2][2] = sz;

	return result;
		
}

Matrix4 rotationMatrix(double degrees_as_angles, double x, double y, double z){
	double degree = degrees_as_angles*M_PI/180 ;

	Vec3 u = Vec3(x,y,z);
	Vec3 v = Vec3();
	if(x==0 && y == 0){
		// if z is the only nonzero variable it should not be zeroed.
		v.x = -z;
		v.y = 0 ;
		v.z = x ;
	}else{
		// either y or x is nonzero, accomplishing a nonzero v.
		v.x = -y;
		v.y = x ;
		v.z = 0 ;
	}
	u = normalizeVec3(u);
	v = normalizeVec3(v);
	Vec3 w = crossProductVec3(u,v);

	// take to x axis, rotate and transform back.

	Matrix4 allignerMatrix = Identity();
	
	allignerMatrix.values[0][0] = u.x;
	allignerMatrix.values[0][1] = u.y;
	allignerMatrix.values[0][2] = u.z;

	allignerMatrix.values[1][0] = v.x;
	allignerMatrix.values[1][1] = v.y;
	allignerMatrix.values[1][2] = v.z;

	allignerMatrix.values[2][0] = w.x;
	allignerMatrix.values[2][1] = w.y;
	allignerMatrix.values[2][2] = w.z;



	Matrix4 pureRotationMatrix = Identity();
	pureRotationMatrix.values[1][1] = cos(degree);
	pureRotationMatrix.values[2][1] = sin(degree);
	pureRotationMatrix.values[1][2] = -sin(degree);
	pureRotationMatrix.values[2][2] = cos(degree);
	




	Matrix4 allignerMatrixInverse = Identity();
	
	allignerMatrixInverse.values[0][0] = u.x;
	allignerMatrixInverse.values[1][0] = u.y;
	allignerMatrixInverse.values[2][0] = u.z;

	allignerMatrixInverse.values[0][1] = v.x;
	allignerMatrixInverse.values[1][1] = v.y;
	allignerMatrixInverse.values[2][1] = v.z;

	allignerMatrixInverse.values[0][2] = w.x;
	allignerMatrixInverse.values[1][2] = w.y;
	allignerMatrixInverse.values[2][2] = w.z;


	return multiplyMatrixWithMatrix(multiplyMatrixWithMatrix(
		allignerMatrixInverse,
		pureRotationMatrix
		),
		allignerMatrix
	);
}


struct OurVertex{
	Vec3 position;
	Vec3 color;
};

struct OurTriangle
{
	struct OurVertex vertices [3];
	Vec3 normal;
};


//Rasterization
double **depthMap;

void fillTriangle(OurTriangle triangle,std::vector<std::vector<Color>> &image) {
	
	double xa = triangle.vertices[0].position.x;
	double ya = triangle.vertices[0].position.y;
	double za = triangle.vertices[0].position.z;
	double c0_red = triangle.vertices[0].color.x;
	double c0_green = triangle.vertices[0].color.y;
	double c0_blue = triangle.vertices[0].color.z;


	double xb = triangle.vertices[1].position.x;
	double yb = triangle.vertices[1].position.y;
	double zb = triangle.vertices[1].position.z;
	double c1_red = triangle.vertices[1].color.x;
	double c1_green = triangle.vertices[1].color.y;
	double c1_blue = triangle.vertices[1].color.z;


	double xc = triangle.vertices[2].position.x;
	double yc = triangle.vertices[2].position.y;
	double zc = triangle.vertices[2].position.z;
	double c2_red = triangle.vertices[2].color.x;
	double c2_green = triangle.vertices[2].color.y;
	double c2_blue = triangle.vertices[2].color.z;

	double ymin = min(ya,min(yb,yc));
	double ymax = max(ya,max(yb,yc));

	double xmin = min(xa,min(xb,xc));
	double xmax = max(xa,max(xb,xc));


    for (int y = ymin; y <= ymax; ++y) {
        for (int x = xmin; x <= xmax; ++x) {
			
			if(x < 0 || y < 0){
				continue;
			}

            double numerator_alpha = -(x-xb)*(yc-yb)+(y-yb)*(xc-xb);
			double denumerator_alpha = -(xa-xb)*(yc-yb)+(ya-yb)*(xc-xb);
			double numerator_beta = -(x-xc)*(ya-yc)+(y-yc)*(xa-xc);
			double denumerator_beta = -(xb-xc)*(ya-yc) + (yb-yc)*(xa-xc);
			
			double alpha = numerator_alpha/denumerator_alpha;
            double beta = numerator_beta / denumerator_beta;
            double gamma = 1- alpha-beta;

            if (alpha >= 0 && beta >= 0 && gamma >= 0) {
            	Color color = Color(
					(alpha * c0_red   + beta * c1_red   + gamma * c2_red),
					(alpha * c0_green + beta * c1_green + gamma * c2_green),
					(alpha * c0_blue  + beta * c1_blue  + gamma * c2_blue)
				);
				
				double depth = za*alpha+zb*beta+zc*gamma;
				if(depth < depthMap[x][y]){
					depthMap[x][y] = depth;
                	image[x][y] = color;	 
				}
            }
        }
    }
}

void drawLineBizzarev1(OurVertex v0, OurVertex v1, std::vector<std::vector<Color>> &image){
	
	double x0 = v0.position.x;
    double x1 = v1.position.x;
    double y0 = v0.position.y;
    double y1 = v1.position.y;
	double z0 = v0.position.z;
    double z1 = v1.position.z;

    double dx = x1 - x0;
    double dy = y1 - y0;

    double c0_red = v0.color.x;
    double c0_green = v0.color.y;
    double c0_blue = v0.color.z;

    double c1_red = v1.color.x;
    double c1_green = v1.color.y;
    double c1_blue = v1.color.z;

    double dc_red = (c1_red - c0_red) / dy;
    double dc_green = (c1_green - c0_green) / dy;
    double dc_blue = (c1_blue - c0_blue) / dy;

	
    double d = -dx + dy/2;
    int x = x0;

    double red = c0_red;
    double green = c0_green;
    double blue = c0_blue;
	
	for (int y = y0; y <= y1; y++) {
		//TODO:remove
		if(x < 0){
			throw "ERROR";
		}
		

		image[x][y] = Color(red, green, blue);

        if (d > 0) { 
            d -= dx;
        } else {     
            d += (dy - dx);
            x++;
        }

        red += dc_red;
        green += dc_green;
        blue += dc_blue;
    }

}

void drawLineBizzarev2(OurVertex v0, OurVertex v1, std::vector<std::vector<Color>> &image){
	
	double x0 = v0.position.x;
    double x1 = v1.position.x;
    double y0 = v0.position.y;
    double y1 = v1.position.y;
	double z0 = v0.position.z;
    double z1 = v1.position.z;

    double dx = x1 - x0;
    double dy = y1 - y0;

    double c0_red = v0.color.x;
    double c0_green = v0.color.y;
    double c0_blue = v0.color.z;

    double c1_red = v1.color.x;
    double c1_green = v1.color.y;
    double c1_blue = v1.color.z;

    double dc_red = (c1_red - c0_red) / dy;
    double dc_green = (c1_green - c0_green) / dy;
    double dc_blue = (c1_blue - c0_blue) / dy;

	
    double d = -dx - dy/2;
    int x = x0;

    double red = c0_red;
    double green = c0_green;
    double blue = c0_blue;
	
	for (int y = y0; y <= y1; y++) {
		if(x < 0){
			continue;
		}
		

		image[x][y] = Color(red, green, blue);

        if (d < 0) { 
            d -= dx;
        } else {     
            d += (-dy - dx);
            x--;
        }

        red += dc_red;
        green += dc_green;
        blue += dc_blue;
    }
}



void drawLineRegularv2(OurVertex v0, OurVertex v1, std::vector<std::vector<Color>> &image){

	double x0 = v0.position.x;
    double x1 = v1.position.x;
    double y0 = v0.position.y;
    double y1 = v1.position.y;
	double z0 = v0.position.z;
    double z1 = v1.position.z;

    double dx = x1 - x0;
    double dy = y1 - y0;

    double c0_red = v0.color.x;
    double c0_green = v0.color.y;
    double c0_blue = v0.color.z;

    double c1_red = v1.color.x;
    double c1_green = v1.color.y;
    double c1_blue = v1.color.z;

    double dc_red = (c1_red - c0_red) / dx;
    double dc_green = (c1_green - c0_green) / dx;
    double dc_blue = (c1_blue - c0_blue) / dx;


    double d = dy + dx/2;
    int y = y0;

    double red = c0_red;
    double green = c0_green;
    double blue = c0_blue;
	
	for (int x = x0; x <= x1; x++) {
		if(y < 0){
			continue;
		}
		

		image[x][y] = Color(red, green, blue);

        if (d < 0) { 
            d -= dy;
        } else {     
            d -= (dy + dx);
            y--;
        }

        red += dc_red;
        green += dc_green;
        blue += dc_blue;
    }

}



void drawLineRegularv1(OurVertex v0, OurVertex v1, std::vector<std::vector<Color>> &image){

	double x0 = v0.position.x;
    double x1 = v1.position.x;
    double y0 = v0.position.y;
    double y1 = v1.position.y;
	double z0 = v0.position.z;
    double z1 = v1.position.z;

    double dx = x1 - x0;
    double dy = y1 - y0;

    double c0_red = v0.color.x;
    double c0_green = v0.color.y;
    double c0_blue = v0.color.z;

    double c1_red = v1.color.x;
    double c1_green = v1.color.y;
    double c1_blue = v1.color.z;

    double dc_red = (c1_red - c0_red) / dx;
    double dc_green = (c1_green - c0_green) / dx;
    double dc_blue = (c1_blue - c0_blue) / dx;


    double d = dy - dx/2;
    int y = y0;

    double red = c0_red;
    double green = c0_green;
    double blue = c0_blue;
	
	for (int x = x0; x <= x1; x++) {
		
		

		image[x][y] = Color(red, green, blue);

        if (d < 0) { 
            d += dy;
        } else {     
            d += (dy - dx);
            y++;
        }

        red += dc_red;
        green += dc_green;
        blue += dc_blue;
    }

}


void drawLine(OurVertex v0, OurVertex v1, std::vector<std::vector<Color>> &image) {
    // Normalize vertices
    
	if (v0.position.x > v1.position.x || (v0.position.x == v1.position.x && v0.position.y > v1.position.y)) {
        std::swap(v0, v1);
    }

    double x0 = v0.position.x;
    double x1 = v1.position.x;
    double y0 = v0.position.y;
    double y1 = v1.position.y;
    double dx = x1 - x0;
    double dy = y1 - y0;

    if (dx == 0) {
        // Vertical line
        if (y0 > y1) std::swap(y0, y1); // Ensure y0 < y1
        for (int y = y0; y <= y1; y++) {
            image[x0][y] = Color(v0.color.x, v0.color.y, v0.color.z); // Flat color
        }
    } else if (dy == 0) {
        // Horizontal line
        for (int x = x0; x <= x1; x++) {
            image[x][y0] = Color(v0.color.x, v0.color.y, v0.color.z); // Flat color
        }
    } else {
        // General case
        double slope = dy / dx;
        if (-1 <= slope && slope <= 1) {
            if (slope >= 0) {
                drawLineRegularv1(v0, v1, image);
            } else {
                drawLineRegularv2(v0, v1, image);
            }
        } else {
            if (slope >= 0) {
                drawLineBizzarev1(v0, v1, image);
            } else {
                drawLineBizzarev2(v1, v0, image);
            }
        }
    }
}


std::vector<std::pair<OurVertex,OurVertex>> clippingPipeline(OurTriangle original){
	std::vector<OurVertex> i_l ;
	i_l.push_back(original.vertices[0]);
	i_l.push_back(original.vertices[1]);
	i_l.push_back(original.vertices[2]);

	//left
	std::vector<OurVertex> o_l_i_r ;
	for(int i = 0; i < i_l.size(); i++){
		OurVertex current_vertex = i_l[i];
		OurVertex next_vertex;
		if(i == i_l.size()-1){
			next_vertex = i_l[0];	
		}else{
			next_vertex = i_l[i+1];
		}

		if(current_vertex.position.x > -1 && next_vertex.position.x > -1){
			o_l_i_r.push_back(next_vertex);
		}else if(current_vertex.position.x < -1 && next_vertex.position.x < -1){
			//no addition
		}else{
			//add intermediate vertex
			double alpha = ((-1)-next_vertex.position.x)/(current_vertex.position.x-next_vertex.position.x);
			OurVertex intermediateVertex;
			intermediateVertex.position.x = alpha*current_vertex.position.x + (1-alpha)*next_vertex.position.x;
			intermediateVertex.position.y = alpha*current_vertex.position.y + (1-alpha)*next_vertex.position.y;
			intermediateVertex.position.z = alpha*current_vertex.position.z + (1-alpha)*next_vertex.position.z;
			
			intermediateVertex.color.x = alpha*current_vertex.color.x + (1-alpha)*next_vertex.color.x;
			intermediateVertex.color.y = alpha*current_vertex.color.y + (1-alpha)*next_vertex.color.y;
			intermediateVertex.color.z = alpha*current_vertex.color.z + (1-alpha)*next_vertex.color.z;
			
			o_l_i_r.push_back(intermediateVertex);		
			if(current_vertex.position.x > -1 && next_vertex.position.x < -1){

			}else if(current_vertex.position.x < -1 && next_vertex.position.x > -1){
				o_l_i_r.push_back(next_vertex);
			}else{
				throw "ERROR";
			}
		}
		
	}
	//right
	std::vector<OurVertex> o_r_i_t ;

	for(int i = 0; i < o_l_i_r.size(); i++){
		OurVertex current_vertex = o_l_i_r[i];
		OurVertex next_vertex;
		if(i == o_l_i_r.size()-1){
			next_vertex = o_l_i_r[0];	
		}else{
			next_vertex = o_l_i_r[i+1];
		}

		if(current_vertex.position.x < 1 && next_vertex.position.x < 1){
			o_r_i_t.push_back(next_vertex);
		}else if(current_vertex.position.x > 1 && next_vertex.position.x > 1){
			//no addition
		}else{
			//add intermediate vertex
			double alpha = ((1)-next_vertex.position.x)/(current_vertex.position.x-next_vertex.position.x);
			OurVertex intermediateVertex;
			intermediateVertex.position.x = alpha*current_vertex.position.x + (1-alpha)*next_vertex.position.x;
			intermediateVertex.position.y = alpha*current_vertex.position.y + (1-alpha)*next_vertex.position.y;
			intermediateVertex.position.z = alpha*current_vertex.position.z + (1-alpha)*next_vertex.position.z;
			
			intermediateVertex.color.x = alpha*current_vertex.color.x + (1-alpha)*next_vertex.color.x;
			intermediateVertex.color.y = alpha*current_vertex.color.y + (1-alpha)*next_vertex.color.y;
			intermediateVertex.color.z = alpha*current_vertex.color.z + (1-alpha)*next_vertex.color.z;
			
			o_r_i_t.push_back(intermediateVertex);		
			if(current_vertex.position.x < 1 && next_vertex.position.x > 1){

			}else if(current_vertex.position.x > 1 && next_vertex.position.x < 1){
				o_r_i_t.push_back(next_vertex);
			}else{
				throw "ERROR";
			}
		}
		
	}
	
	
	
	//top
	std::vector<OurVertex> o_t_i_b ;

	for(int i = 0; i < o_r_i_t.size(); i++){
		OurVertex current_vertex = o_r_i_t[i];
		OurVertex next_vertex;
		if(i == o_r_i_t.size()-1){
			next_vertex = o_r_i_t[0];	
		}else{
			next_vertex = o_r_i_t[i+1];
		}

		if(current_vertex.position.y < 1 && next_vertex.position.y < 1){
			o_t_i_b.push_back(next_vertex);
		}else if(current_vertex.position.y > 1 && next_vertex.position.y > 1){
			//no addition
		}else{
			//add intermediate vertex
			double alpha = ((1)-next_vertex.position.y)/(current_vertex.position.y-next_vertex.position.y);
			OurVertex intermediateVertex;
			intermediateVertex.position.x = alpha*current_vertex.position.x + (1-alpha)*next_vertex.position.x;
			intermediateVertex.position.y = alpha*current_vertex.position.y + (1-alpha)*next_vertex.position.y;
			intermediateVertex.position.z = alpha*current_vertex.position.z + (1-alpha)*next_vertex.position.z;
			
			intermediateVertex.color.x = alpha*current_vertex.color.x + (1-alpha)*next_vertex.color.x;
			intermediateVertex.color.y = alpha*current_vertex.color.y + (1-alpha)*next_vertex.color.y;
			intermediateVertex.color.z = alpha*current_vertex.color.z + (1-alpha)*next_vertex.color.z;
			
			o_t_i_b.push_back(intermediateVertex);		
			if(current_vertex.position.y < 1 && next_vertex.position.y > 1){

			}else if(current_vertex.position.y > 1 && next_vertex.position.y < 1){
				o_t_i_b.push_back(next_vertex);
			}else{
				throw "ERROR";
			}
		}
		
	}

	//bottom
	std::vector<OurVertex> o_b_i_f ;
	
	for(int i = 0; i < o_t_i_b.size(); i++){
		OurVertex current_vertex = o_t_i_b[i];
		OurVertex next_vertex;
		if(i == o_t_i_b.size()-1){
			next_vertex = o_t_i_b[0];	
		}else{
			next_vertex = o_t_i_b[i+1];
		}

		if(current_vertex.position.y > -1 && next_vertex.position.y > -1){
			o_b_i_f.push_back(next_vertex);
		}else if(current_vertex.position.y < -1 && next_vertex.position.y < -1){
			//no addition
		}else{
			//add intermediate vertex
			double alpha = ((-1)-next_vertex.position.y)/(current_vertex.position.y-next_vertex.position.y);
			OurVertex intermediateVertex;
			intermediateVertex.position.x = alpha*current_vertex.position.x + (1-alpha)*next_vertex.position.x;
			intermediateVertex.position.y = alpha*current_vertex.position.y + (1-alpha)*next_vertex.position.y;
			intermediateVertex.position.z = alpha*current_vertex.position.z + (1-alpha)*next_vertex.position.z;
			
			intermediateVertex.color.x = alpha*current_vertex.color.x + (1-alpha)*next_vertex.color.x;
			intermediateVertex.color.y = alpha*current_vertex.color.y + (1-alpha)*next_vertex.color.y;
			intermediateVertex.color.z = alpha*current_vertex.color.z + (1-alpha)*next_vertex.color.z;
			
			o_b_i_f.push_back(intermediateVertex);		
			if(current_vertex.position.y > -1 && next_vertex.position.y < -1){

			}else if(current_vertex.position.y < -1 && next_vertex.position.y > -1){
				o_b_i_f.push_back(next_vertex);
			}else{
				throw "ERROR";
			}
		}
		
	}
	//far
	std::vector<OurVertex> o_f_i_n ;

	for(int i = 0; i < o_b_i_f.size(); i++){
		OurVertex current_vertex = o_b_i_f[i];
		OurVertex next_vertex;
		if(i == o_b_i_f.size()-1){
			next_vertex = o_b_i_f[0];	
		}else{
			next_vertex = o_b_i_f[i+1];
		}

		if(current_vertex.position.z < 1 && next_vertex.position.z < 1){
			o_f_i_n.push_back(next_vertex);
		}else if(current_vertex.position.z > 1 && next_vertex.position.z > 1){
			//no addition
		}else{
			//add intermediate vertex
			double alpha = ((1)-next_vertex.position.z)/(current_vertex.position.z-next_vertex.position.z);
			OurVertex intermediateVertex;
			intermediateVertex.position.x = alpha*current_vertex.position.x + (1-alpha)*next_vertex.position.x;
			intermediateVertex.position.y = alpha*current_vertex.position.y + (1-alpha)*next_vertex.position.y;
			intermediateVertex.position.z = alpha*current_vertex.position.z + (1-alpha)*next_vertex.position.z;
			
			intermediateVertex.color.x = alpha*current_vertex.color.x + (1-alpha)*next_vertex.color.x;
			intermediateVertex.color.y = alpha*current_vertex.color.y + (1-alpha)*next_vertex.color.y;
			intermediateVertex.color.z = alpha*current_vertex.color.z + (1-alpha)*next_vertex.color.z;
			
			o_f_i_n.push_back(intermediateVertex);		
			if(current_vertex.position.z < 1 && next_vertex.position.z > 1){

			}else if(current_vertex.position.z > 1 && next_vertex.position.z < 1){
				o_f_i_n.push_back(next_vertex);
			}else{
				throw "ERROR";
			}
		}
		
	}

	//near
	std::vector<OurVertex> o_n ;
	for(int i = 0; i < o_f_i_n.size(); i++){
		OurVertex current_vertex = o_f_i_n[i];
		OurVertex next_vertex;
		if(i == o_f_i_n.size()-1){
			next_vertex = o_f_i_n[0];	
		}else{
			next_vertex = o_f_i_n[i+1];
		}

		if(current_vertex.position.z > -1 && next_vertex.position.z > -1){
			o_n.push_back(next_vertex);
		}else if(current_vertex.position.z < -1 && next_vertex.position.z < -1){
			//no addition
		}else{
			//add intermediate vertex
			double alpha = ((-1)-next_vertex.position.z)/(current_vertex.position.z-next_vertex.position.z);
			OurVertex intermediateVertex;
			intermediateVertex.position.x = alpha*current_vertex.position.x + (1-alpha)*next_vertex.position.x;
			intermediateVertex.position.y = alpha*current_vertex.position.y + (1-alpha)*next_vertex.position.y;
			intermediateVertex.position.z = alpha*current_vertex.position.z + (1-alpha)*next_vertex.position.z;
			
			intermediateVertex.color.x = alpha*current_vertex.color.x + (1-alpha)*next_vertex.color.x;
			intermediateVertex.color.y = alpha*current_vertex.color.y + (1-alpha)*next_vertex.color.y;
			intermediateVertex.color.z = alpha*current_vertex.color.z + (1-alpha)*next_vertex.color.z;
			
			o_n.push_back(intermediateVertex);		
			if(current_vertex.position.z > -1 && next_vertex.position.z < -1){

			}else if(current_vertex.position.z < -1 && next_vertex.position.z > -1){
				o_n.push_back(next_vertex);
			}else{
				throw "ERROR";
			}
		}
		
	}
	
	std::vector<std::pair<OurVertex,OurVertex>> result;
	for(int i = 0 ; i < o_n.size(); i++){
		if(i != o_n.size()-1){
			result.push_back(std::pair(o_n[i],o_n[i+1]));
		}else{
			result.push_back(std::pair(o_n[i],o_n[0]));
		}
	}
	return result;
}



/*
	Transformations, clipping, culling, rasterization are done here.
*/
bool insideCvv(Vec3 position){
	if(
		position.x <= -1 || position.x >= 1 ||
		position.y <= -1 || position.y >= 1 ||
		position.z <= -1 || position.z >= 1 
	){
		return false;
	}else{
		return true;
	}
}

void Scene::forwardRenderingPipeline(Camera *camera)
{
	if(image.size() != 0 && image[0].size() !=0 ){
		depthMap = new double*[image.size()];
		for(int i = 0 ; i < image.size(); i++){
			depthMap[i] = new double[image[0].size()];
		}
	}else{
		throw "ERROR";
	}
	for(int i = 0; i < image.size(); i++){
		for(int j = 0; j < image[0].size(); j++){
			depthMap[i][j] = 100;
		}
	}

	Matrix4 cameraTranslationInverter = translationMatrix(
		- camera->position.x,
		- camera->position.y,
		- camera->position.z
	);
	Matrix4 cameraRotationInverter = Identity();
	cameraRotationInverter.values[0][0] = camera->u.x;
	cameraRotationInverter.values[0][1] = camera->u.y;
	cameraRotationInverter.values[0][2] = camera->u.z;
	
	cameraRotationInverter.values[1][0] = camera->v.x;
	cameraRotationInverter.values[1][1] = camera->v.y;
	cameraRotationInverter.values[1][2] = camera->v.z;
	
	cameraRotationInverter.values[2][0] = camera->w.x;
	cameraRotationInverter.values[2][1] = camera->w.y;
	cameraRotationInverter.values[2][2] = camera->w.z;
	
	Matrix4 cameraTransformation = multiplyMatrixWithMatrix(cameraRotationInverter,cameraTranslationInverter);
	
	Matrix4 projectionTranformation = Identity();
	if(camera->projectionType == 0){
		//orthographic
		projectionTranformation.values[0][0] = 2 / (camera->right - camera->left);
		projectionTranformation.values[1][1] = 2 / (camera->top - camera->bottom);
		projectionTranformation.values[2][2] = -2 / (camera->far - camera->near);

		projectionTranformation.values[0][3] = -(camera->right + camera->left) / (camera->right - camera->left);
		projectionTranformation.values[1][3] = -(camera->top + camera->bottom) / (camera->top - camera->bottom);
		projectionTranformation.values[2][3] = -(camera->far + camera->near) / (camera->far - camera->near);

	}else if(camera->projectionType == 1){
		//perspective

		Matrix4 p20 = Identity();

		p20.values[0][0] = camera->near;
		p20.values[1][1] = camera->near;
		p20.values[2][2] = camera->far + camera->near;
		p20.values[2][3] = camera->far * camera->near;
		p20.values[3][3] = 0;
		p20.values[3][2] = -1;

		projectionTranformation.values[0][0] = 2 / (camera->right - camera->left);
		projectionTranformation.values[1][1] = 2 / (camera->top - camera->bottom);
		projectionTranformation.values[2][2] = -2 / (camera->far - camera->near);

		projectionTranformation.values[0][3] = -(camera->right + camera->left) / (camera->right - camera->left);
		projectionTranformation.values[1][3] = -(camera->top + camera->bottom) / (camera->top - camera->bottom);
		projectionTranformation.values[2][3] = -(camera->far + camera->near) / (camera->far - camera->near);


		projectionTranformation = multiplyMatrixWithMatrix(projectionTranformation, p20);
		
	}else{
		throw "ERROR" ;
	}
	
	Matrix4 ProjectionCameraUnified = multiplyMatrixWithMatrix(projectionTranformation,cameraTransformation);
	
	Matrix4 viewportTransformationMatrix = Identity();
	viewportTransformationMatrix.values[0][0] = ((double)camera->horRes)/2;
	viewportTransformationMatrix.values[1][1] = ((double)camera->verRes)/2;
	viewportTransformationMatrix.values[2][2] = 0.5;
	
	viewportTransformationMatrix.values[0][3] = ((double)camera->horRes-1)/2;
	viewportTransformationMatrix.values[1][3] = ((double)camera->verRes-1)/2;
	viewportTransformationMatrix.values[2][3] = 0.5;
	
	

	
	for(int meshId = 0 ; meshId< this->meshes.size();meshId++){
		Mesh* currentMesh = meshes[meshId];
		
		Matrix4 modelingTransformation = Identity();
		for(int transformationIndexer = 0; transformationIndexer < currentMesh->transformationIds.size(); transformationIndexer++){
			int  transformationId = currentMesh->transformationIds[transformationIndexer];
			char transformationType = currentMesh->transformationTypes[transformationIndexer];
			
			if(transformationType == 'r'){
				modelingTransformation = multiplyMatrixWithMatrix(	
					rotationMatrix(
						rotations[transformationId-1]->angle,
						rotations[transformationId-1]->ux,
						rotations[transformationId-1]->uy,
						rotations[transformationId-1]->uz
					),
					modelingTransformation
				);				
				
			}else if(transformationType == 't'){
				modelingTransformation = multiplyMatrixWithMatrix(	
					translationMatrix(
						translations[transformationId-1]->tx,
						translations[transformationId-1]->ty,
						translations[transformationId-1]->tz
					),
					modelingTransformation
				);
			}else if(transformationType == 's'){
				modelingTransformation = multiplyMatrixWithMatrix(	
					scalingMatrix(
						scalings[transformationId-1]->sx,
						scalings[transformationId-1]->sy,
						scalings[transformationId-1]->sz
					),
					modelingTransformation
				);
			}else{
				throw "ERROR";
			}



		}
		Matrix4 modeling_to_projection = multiplyMatrixWithMatrix(ProjectionCameraUnified,modelingTransformation);
		
		
		std::vector<std::pair<OurVertex,OurVertex>> ourLines;

		for(int i = 0; i < currentMesh->triangles.size(); i++){
			OurTriangle currentTriangle;

			currentTriangle.vertices[0].position = 
			*(
				this->vertices[
					currentMesh->triangles[i].vertexIds[0]-1
				]
			);
			
			
			currentTriangle.vertices[1].position = 
			*(
				this->vertices[
					currentMesh->triangles[i].vertexIds[1]-1
				]
			);

			currentTriangle.vertices[2].position = 
			*(
				this->vertices[
					currentMesh->triangles[i].vertexIds[2]-1
				]
			);




			currentTriangle.vertices[0].color.x = 
			(
				this->colorsOfVertices[
					currentMesh->triangles[i].vertexIds[0]-1
				]
			)->r;
			currentTriangle.vertices[0].color.y = 
			(
				this->colorsOfVertices[
					currentMesh->triangles[i].vertexIds[0]-1
				]
			)->g;
			currentTriangle.vertices[0].color.z = 
			(
				this->colorsOfVertices[
					currentMesh->triangles[i].vertexIds[0]-1
				]
			)->b;




			currentTriangle.vertices[1].color.x = 
			(
				this->colorsOfVertices[
					currentMesh->triangles[i].vertexIds[1]-1
				]
			)->r;
			currentTriangle.vertices[1].color.y = 
			(
				this->colorsOfVertices[
					currentMesh->triangles[i].vertexIds[1]-1
				]
			)->g;
			currentTriangle.vertices[1].color.z = 
			(
				this->colorsOfVertices[
					currentMesh->triangles[i].vertexIds[1]-1
				]
			)->b;



			currentTriangle.vertices[2].color.x = 
			(
				this->colorsOfVertices[
					currentMesh->triangles[i].vertexIds[2]-1
				]
			)->r;
			currentTriangle.vertices[2].color.y = 
			(
				this->colorsOfVertices[
					currentMesh->triangles[i].vertexIds[2]-1
				]
			)->g;
			currentTriangle.vertices[2].color.z = 
			(
				this->colorsOfVertices[
					currentMesh->triangles[i].vertexIds[2]-1
				]
			)->b;


			Vec4 position0 = Vec4(
				currentTriangle.vertices[0].position.x,
				currentTriangle.vertices[0].position.y,
				currentTriangle.vertices[0].position.z,
				1
			);

			Vec4 position1 = Vec4(
				currentTriangle.vertices[1].position.x,
				currentTriangle.vertices[1].position.y,
				currentTriangle.vertices[1].position.z,
				1
			);

			Vec4 position2 = Vec4(
				currentTriangle.vertices[2].position.x,
				currentTriangle.vertices[2].position.y,
				currentTriangle.vertices[2].position.z,
				1
			);

			position0 = multiplyMatrixWithVec4(modeling_to_projection,position0);
			position1 = multiplyMatrixWithVec4(modeling_to_projection,position1);
			position2 = multiplyMatrixWithVec4(modeling_to_projection,position2);

			position0 = homogenize(position0);
			position1 = homogenize(position1);
			position2 = homogenize(position2);
		
			if(! insideCvv(Vec3(position0.x,position0.y, position0.z)) &&
			   ! insideCvv(Vec3(position1.x,position1.y, position1.z)) && 
			   ! insideCvv(Vec3(position2.x,position2.y, position2.z))){
				continue;
			}

			
			currentTriangle.vertices[0].position.x = position0.x;
			currentTriangle.vertices[0].position.y = position0.y;
			currentTriangle.vertices[0].position.z = position0.z;

			currentTriangle.vertices[1].position.x = position1.x;
			currentTriangle.vertices[1].position.y = position1.y;
			currentTriangle.vertices[1].position.z = position1.z;

			currentTriangle.vertices[2].position.x = position2.x;
			currentTriangle.vertices[2].position.y = position2.y;
			currentTriangle.vertices[2].position.z = position2.z;

			if(this->cullingEnabled){
				currentTriangle.normal = crossProductVec3(
					subtractVec3(currentTriangle.vertices[1].position,currentTriangle.vertices[0].position),
					subtractVec3(currentTriangle.vertices[2].position,currentTriangle.vertices[0].position)
				);
						
						

				Vec3 center = Vec3 ((currentTriangle.vertices[0].position.x + currentTriangle.vertices[1].position.x + currentTriangle.vertices[2].position.x)/3,
									(currentTriangle.vertices[0].position.y + currentTriangle.vertices[1].position.y + currentTriangle.vertices[2].position.y)/3,
									(currentTriangle.vertices[0].position.z + currentTriangle.vertices[1].position.z + currentTriangle.vertices[2].position.z)/3
							);				
				Vec3 viewPosition = Vec3(0,0,-1);
				Vec3 viewVector = Vec3(
					center.x-viewPosition.x,
					center.y-viewPosition.y,
					center.z-viewPosition.z
				);

				//DANGER TODO
				if(dotProductVec3(viewVector,currentTriangle.normal) <= 0){
					continue;
				}

			}


			
			if(currentMesh->type==0){
				ourLines = clippingPipeline(currentTriangle);	
				for(int lineAdditionIndex = 0 ; lineAdditionIndex < ourLines.size(); lineAdditionIndex++){
					
					std::pair<OurVertex, OurVertex> assesedLine = ourLines[lineAdditionIndex];

					position0 = Vec4( 
						assesedLine.first.position.x,
						assesedLine.first.position.y,
						assesedLine.first.position.z,
						1
					);
					position1 = Vec4( 
						assesedLine.second.position.x,
						assesedLine.second.position.y,
						assesedLine.second.position.z,
						1
					);

					position0 = multiplyMatrixWithVec4(viewportTransformationMatrix,position0);
					position1 = multiplyMatrixWithVec4(viewportTransformationMatrix,position1);
					
					assesedLine.first.position = Vec3(
						position0.x,
						position0.y,
						position0.z
					);
					assesedLine.second.position = Vec3(
						position1.x,
						position1.y,
						position1.z
					);

					drawLine(assesedLine.first,assesedLine.second,image);

				}
			}else if(currentMesh->type==1){
					position0 = Vec4( 
						currentTriangle.vertices[0].position.x,
						currentTriangle.vertices[0].position.y,
						currentTriangle.vertices[0].position.z,
						1
					);
					position1 = Vec4( 
						currentTriangle.vertices[1].position.x,
						currentTriangle.vertices[1].position.y,
						currentTriangle.vertices[1].position.z,
						1
					);
					position2 = Vec4( 
						currentTriangle.vertices[2].position.x,
						currentTriangle.vertices[2].position.y,
						currentTriangle.vertices[2].position.z,
						1
					);

					position0 = multiplyMatrixWithVec4(viewportTransformationMatrix,position0);
					position1 = multiplyMatrixWithVec4(viewportTransformationMatrix,position1);
					position2 = multiplyMatrixWithVec4(viewportTransformationMatrix,position2);
					
					currentTriangle.vertices[0].position.x = position0.x;
					currentTriangle.vertices[0].position.y = position0.y;
					currentTriangle.vertices[0].position.z = position0.z;

					currentTriangle.vertices[1].position.x = position1.x;
					currentTriangle.vertices[1].position.y = position1.y;
					currentTriangle.vertices[1].position.z = position1.z;

					currentTriangle.vertices[2].position.x = position2.x;
					currentTriangle.vertices[2].position.y = position2.y;
					currentTriangle.vertices[2].position.z = position2.z;
					
					fillTriangle(currentTriangle,image);
			}else{

			}
		}

		
		
	}
}



//TODO: handle deallocations (do ctrl f for all news and confirm they are properly deallocated).
